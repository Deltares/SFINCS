! Non-hydrostatic pressure projection -- VERSION 2 (assembled-stencil, GPU-ready).
!
! Same physics as sfincs_nonhydrostatic.f90 (v1): single-layer Keller-box pressure
! projection, SPD operator A = (dt/rho) G^T diag(hu) G + diag(kbfac*dt/(rho*H)),
! solved with conjugate gradients, warm-started, with the same breaking (HFA /
! Froude), fade-in, filtering, capping and w_s/w_b bookkeeping. What changes is
! the SOLVER MACHINERY, redesigned for memory bandwidth and GPU offload:
!
!  1) ASSEMBLED 5-POINT STENCIL. v1 is matrix-free: every CG iteration re-reads
!     the per-face data (nh_hu, nh_cR, nh_cL twice -- once per side) and re-does
!     the products dtrho*hu*cR*cL, in TWO passes (face pass + row gather) with a
!     barrier between them. Here the symmetric 5-point operator is assembled ONCE
!     PER SOLVE into four off-diagonal arrays (nh2_aw/ae/as/an) with STATIC
!     neighbour indices (nh2_nbw/nbe/nbs/nbn). The matvec is then a single loop:
!         w(i) = r(i) + aw*r(nbw) + ae*r(nbe) + as*r(nbs) + an*r(nbn)
!     -> roughly half the memory traffic per iteration, one kernel instead of two.
!
!  2) SYMMETRIC JACOBI SCALING BAKED IN. The system is solved as
!         (S A S) y = S b,   pnh = S y,   S = diag(1/sqrt(diag(A))),
!     so the scaled operator has UNIT diagonal: the Jacobi preconditioner
!     disappears from the iteration (no cg_z array, no divisions) and the CG
!     identity <r,z> = <r,r> removes one reduction term.
!
!  3) BREAKING DIRICHLET BAKED IN. v1 zeroes the residual at breaking cells every
!     iteration (a branch in the hot loop). Here the row AND column of a breaking
!     cell are zeroed at assembly time (sclo = scl * activity, symmetric), the
!     scaled diagonal stays 1, and b=0 there -> y=0 is held exactly with NO
!     per-iteration branch.
!
!  4) GHOST ELEMENT 0. CG vectors run (0:nrows) with element 0 pinned at zero;
!     a missing neighbour has index 0 and a zero coefficient (sclo(0)=0), so the
!     matvec has NO conditionals at all -- fully vectorizable on CPU, coalesced
!     on GPU.
!
!  5) TWO FUSED KERNELS PER ITERATION (Chronopoulos-Gear single-reduction CG):
!         kernel B:  w = A r  fused with the dots  gam=<r,r>, del=<r,w>
!         kernel A:  p=r+beta*p ; s=w+beta*s ; y+=alpha*p ; r-=alpha*s
!     One host synchronisation point per iteration (the reduction). Dot products
!     accumulate in real*8 (the real*4 vectors are unchanged).
!
!  6) DUAL DIRECTIVES. Every kernel carries both !$omp (CPU, schedule(static) --
!     uniform work, keeps thread-to-data affinity across loops) and !$acc (GPU)
!     directives, following the sfincs_momentum.f90 convention. All module arrays
!     are made device-resident at initialization (enter data); per-step work then
!     runs entirely on the device. The host copy of pnh is refreshed only at map
!     output times ('update host(pnh) if_present' in sfincs_output, alongside
!     the other output transfers). NOTE: the convergence test uses the SCALED residual norm
!     ||S r|| / ||S b|| (the natural Jacobi-CG norm); v1 tests the unscaled one.
!     Same nh_tol applies.
!
! Grid/face mappings, bed slopes and the open-boundary fade-in are IDENTICAL to
! v1, so initialization simply calls initialize_nonhydrostatic() and reuses the
! v1 module state (row maps, face maps, pnh/ws/wb/wb0, which sfincs_ncoutput
! also reads). Selected with nh_version = 2 in sfincs.inp.
!
module sfincs_nonhydrostatic_v2
   !
   use sfincs_nonhydrostatic
   !
   implicit none
   !
   ! Static 5-point neighbour map: row index of the West/East/South/North
   ! neighbour of each row (0 = no nonh neighbour on that side -> ghost row 0).
   !
   integer, dimension(:), allocatable :: nh2_nbw, nh2_nbe, nh2_nbs, nh2_nbn   ! (nrows)
   !
   ! Assembled, symmetrically-scaled off-diagonals of the operator (unit diagonal),
   ! rebuilt every solve (wet/dry, depths and breaking change each step).
   !
   real*4, dimension(:), allocatable :: nh2_aw, nh2_ae, nh2_as, nh2_an        ! (nrows)
   real*4, dimension(:), allocatable :: nh2_scl    ! 1/sqrt(diag(A))            (nrows)
   real*4, dimension(:), allocatable :: nh2_sclo   ! scl * breaking-activity  (0:nrows), ghost 0
   !
   ! CG vectors, ghost element 0 pinned at zero (never written by any kernel).
   !
   real*4, dimension(:), allocatable :: nh2_y      ! scaled pressure iterate  (0:nrows)
   real*4, dimension(:), allocatable :: nh2_r      ! residual                 (0:nrows)
   real*4, dimension(:), allocatable :: nh2_p      ! search direction         (0:nrows)
   real*4, dimension(:), allocatable :: nh2_s      ! s = A p (recurrence)     (0:nrows)
   real*4, dimension(:), allocatable :: nh2_w      ! w = A r                  (0:nrows)
   real*4, dimension(:), allocatable :: nh2_b      ! scaled right-hand side   (0:nrows)
   !
contains
   !
   subroutine initialize_nonhydrostatic_v2()
   !
   ! v1 initialization builds everything grid-related (row<->nm maps, nonh face
   ! maps nh_faceL/R/uv + nh_cellface, per-row spacing nh_dxr/nh_dyr, per-face
   ! gradient weight nh_cf, bed slopes, open-boundary fade-in). Here we only add
   ! the static neighbour map and the v2 work arrays, and push everything the
   ! per-step kernels touch onto the GPU.
   !
   use sfincs_data
   use sfincs_log
   !
   implicit none
   !
   integer :: i, f
   !
   call initialize_nonhydrostatic()
   !
   if (nrows == 0) return
   !
   allocate(nh2_nbw(nrows))
   allocate(nh2_nbe(nrows))
   allocate(nh2_nbs(nrows))
   allocate(nh2_nbn(nrows))
   !
   ! Neighbour row of each row through its (up to 4) incident nonh faces:
   ! slot 1 = left face  -> neighbour is that face's L cell ; slot 2 = right -> R cell;
   ! slot 3 = below face -> L cell ; slot 4 = above -> R cell (see v1 nh_cellface).
   ! nh_faceL/R are 0 on a domain/open-boundary side -> neighbour 0 -> ghost.
   !
   do i = 1, nrows
      nh2_nbw(i) = 0 ; f = nh_cellface(1, i) ; if (f > 0) nh2_nbw(i) = nh_faceL(f)
      nh2_nbe(i) = 0 ; f = nh_cellface(2, i) ; if (f > 0) nh2_nbe(i) = nh_faceR(f)
      nh2_nbs(i) = 0 ; f = nh_cellface(3, i) ; if (f > 0) nh2_nbs(i) = nh_faceL(f)
      nh2_nbn(i) = 0 ; f = nh_cellface(4, i) ; if (f > 0) nh2_nbn(i) = nh_faceR(f)
   enddo
   !
   allocate(nh2_aw(nrows))
   allocate(nh2_ae(nrows))
   allocate(nh2_as(nrows))
   allocate(nh2_an(nrows))
   allocate(nh2_scl(nrows))
   allocate(nh2_sclo(0:nrows))
   allocate(nh2_y(0:nrows))
   allocate(nh2_r(0:nrows))
   allocate(nh2_p(0:nrows))
   allocate(nh2_s(0:nrows))
   allocate(nh2_w(0:nrows))
   allocate(nh2_b(0:nrows))
   !
   ! First-touch initialization (NUMA page placement matches the static loop
   ! partition the per-step kernels use). Ghost element 0 stays 0 forever.
   !
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      nh2_aw(i)   = 0.0
      nh2_ae(i)   = 0.0
      nh2_as(i)   = 0.0
      nh2_an(i)   = 0.0
      nh2_scl(i)  = 1.0
      nh2_sclo(i) = 1.0
      nh2_y(i)    = 0.0
      nh2_r(i)    = 0.0
      nh2_p(i)    = 0.0
      nh2_s(i)    = 0.0
      nh2_w(i)    = 0.0
      nh2_b(i)    = 0.0
   enddo
   !$omp end parallel do
   nh2_sclo(0) = 0.0
   nh2_y(0)    = 0.0
   nh2_r(0)    = 0.0
   nh2_p(0)    = 0.0
   nh2_s(0)    = 0.0
   nh2_w(0)    = 0.0
   nh2_b(0)    = 0.0
   !
   ! Device residency: v2 arrays plus the v1/module arrays the per-step kernels
   ! read or write. The global sfincs_data arrays (zs, zb, uv, q, kfuv, kcs,
   ! z_index_uv_*, uv_index_z_*) are already on the device via initialize_openacc.
   !
   !$acc enter data copyin( nh2_nbw, nh2_nbe, nh2_nbs, nh2_nbn, &
   !$acc                    nh2_aw, nh2_ae, nh2_as, nh2_an, nh2_scl, nh2_sclo, &
   !$acc                    nh2_y, nh2_r, nh2_p, nh2_s, nh2_w, nh2_b, &
   !$acc                    pnh, ws, wb, wb0, Dnm, nm_index_of_row, &
   !$acc                    nh_dvert, nh_fade, nh_brfac, nh_dxr, nh_dyr, nh_slbed, &
   !$acc                    nh_faceuv, nh_faceL, nh_faceR, nh_cellface, &
   !$acc                    nh_cf, nh_cR, nh_cL, nh_hu )
   !
   write(logstr,'(a,i0,a,i0,a)')'Non-hydrostatic solver v2 (assembled-stencil CG): ', nrows, ' cells, ', nhuv, ' faces'
   call write_log(logstr, 0)
   !
   end subroutine


   subroutine compute_nonhydrostatic_v2(dt, tloop)
   !
   ! Per-step non-hydrostatic pressure projection (v2). Physics identical to
   ! v1 compute_nonhydrostatic2; see module header for what changed in the
   ! solver machinery. Every numbered kernel below is one !$omp / !$acc loop.
   !
   use sfincs_data
   !
   implicit none
   !
   integer   :: count0, count1, count_rate, count_max
   real      :: tloop
   real*4    :: dt
   !
   integer   :: i, ip, ipuv, nm, nmn, nmu, iuv, iuvl, iuvr, j, jl, jr, f, iter, ipass, cnt, ineighbour
   real*4    :: dtrho, kbfac, hu, Dnm1, Dnmu, abf, fdep, tf, dval, sq, braw, accv, sumn
   real*4    :: qr, ql, dzdt, wmax, breform, pcap, gf, pL, pR, unh
   real*4    :: alpha, beta
   real*8    :: gam8, del8, bn8, gamold8, alpha8, beta8, pap8, bnorm8
   logical   :: bndok
   !
   call system_clock(count0, count_rate, count_max)
   !
   if (nrows == 0) return
   !
   kbfac = nh_disp
   dtrho = dt / rhow
   !
   ! 1) Layer depth per row and the vertical-acceleration diagonal.
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i, nm )
   do i = 1, nrows
      nm = nm_index_of_row(i)
      Dnm(i)      = max(zs(nm) - zb(nm), huthresh_nh)
      nh_dvert(i) = kbfac * dt / (rhow * Dnm(i))
   enddo
   !$omp end parallel do
   !
   ! 1b) Breaking criterion -> nh_brfac (1 = full nonh, <1 = hydrostatic pnh=0).
   !     Froude (NEOWAVE) when nh_brfr > 0, else steepness/HFA when nh_brsteep > 0.
   !     Same logic as v1; neighbour lookups via the static nh2_nb* map.
   !
   if (nh_brfr > 0.0) then
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, nm, iuvl, iuvr, ql, qr, dzdt )
      do i = 1, nrows
         nm = nm_index_of_row(i)
         ql = 0.0 ; qr = 0.0
         iuvl = z_index_uv_md(nm) ; iuvr = z_index_uv_mu(nm)
         if (iuvl > 0) ql = uv(iuvl)
         if (iuvr > 0) qr = uv(iuvr)
         dzdt = (0.5 * (ql + qr))**2                       ! u_cell^2
         if (nmax > 1) then
            ql = 0.0 ; qr = 0.0
            iuvl = z_index_uv_nd(nm) ; iuvr = z_index_uv_nu(nm)
            if (iuvl > 0) ql = uv(iuvl)
            if (iuvr > 0) qr = uv(iuvr)
            dzdt = dzdt + (0.5 * (ql + qr))**2             ! + v_cell^2
         endif
         dzdt = sqrt(dzdt) / sqrt(g * Dnm(i))              ! Fr = |U| / sqrt(g D)
         if (nh_brfac(i) >= 1.0) then                      ! was not breaking
            if (dzdt > nh_brfr) nh_brfac(i) = 0.0          ! onset
         else                                              ! was breaking -> release only when slow
            if (dzdt < 0.3 * nh_brfr) then
               nh_brfac(i) = 1.0
            else
               nh_brfac(i) = 0.0
            endif
         endif
         if (nh_brfac(i) < 1.0) pnh(i) = 0.0
      enddo
      !$omp end parallel do
      !
   elseif (nh_brsteep > 0.0) then
      !
      breform = 0.25 * nh_brsteep
      !
      ! (a) raw surface rate-of-rise dzdt = -div(q) into nh2_w
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, nm, iuvl, iuvr, ql, qr, dzdt )
      do i = 1, nrows
         nm   = nm_index_of_row(i)
         iuvr = z_index_uv_mu(nm) ; iuvl = z_index_uv_md(nm)
         qr = 0.0 ; ql = 0.0
         if (iuvr > 0) qr = q(iuvr)
         if (iuvl > 0) ql = q(iuvl)
         dzdt = - (qr - ql) * nh_dxr(i)
         if (nmax > 1) then
            iuvr = z_index_uv_nu(nm) ; iuvl = z_index_uv_nd(nm)
            qr = 0.0 ; ql = 0.0
            if (iuvr > 0) qr = q(iuvr)
            if (iuvl > 0) ql = q(iuvl)
            dzdt = dzdt - (qr - ql) * nh_dyr(i)
         endif
         nh2_w(i) = dzdt
      enddo
      !$omp end parallel do
      !
      ! (b) two [1 2 1]/4 smoothing passes (2dx noise -> connected front)
      !
      do ipass = 1, 2
         !$acc parallel loop default(present)
         !$omp parallel do schedule ( static ) private ( i, accv, j )
         do i = 1, nrows
            accv = nh2_w(i)
            j = nh2_nbw(i) ; if (j > 0) accv = accv + 0.25 * (nh2_w(j) - nh2_w(i))
            j = nh2_nbe(i) ; if (j > 0) accv = accv + 0.25 * (nh2_w(j) - nh2_w(i))
            j = nh2_nbs(i) ; if (j > 0) accv = accv + 0.25 * (nh2_w(j) - nh2_w(i))
            j = nh2_nbn(i) ; if (j > 0) accv = accv + 0.25 * (nh2_w(j) - nh2_w(i))
            nh2_p(i) = accv
         enddo
         !$omp end parallel do
         !$acc parallel loop default(present)
         !$omp parallel do schedule ( static ) private ( i )
         do i = 1, nrows
            nh2_w(i) = nh2_p(i)
         enddo
         !$omp end parallel do
      enddo
      !
      ! (c) snapshot last step's breaking state (no within-sweep contamination)
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh2_r(i) = nh_brfac(i)
      enddo
      !$omp end parallel do
      !
      ! (d) hysteretic state machine on the smoothed dzdt (onset / neighbour-spread /
      !     release only when the surface falls)
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, dzdt, ineighbour, j, wmax )
      do i = 1, nrows
         dzdt = nh2_w(i)
         ineighbour = 0
         j = nh2_nbw(i) ; if (j > 0) then ; if (nh2_r(j) < 1.0) ineighbour = 1 ; endif
         j = nh2_nbe(i) ; if (j > 0) then ; if (nh2_r(j) < 1.0) ineighbour = 1 ; endif
         j = nh2_nbs(i) ; if (j > 0) then ; if (nh2_r(j) < 1.0) ineighbour = 1 ; endif
         j = nh2_nbn(i) ; if (j > 0) then ; if (nh2_r(j) < 1.0) ineighbour = 1 ; endif
         wmax = sqrt(g * Dnm(i))
         if (nh2_r(i) >= 1.0) then                                  ! was not breaking
            if (dzdt > nh_brsteep * wmax) then
               nh_brfac(i) = 0.0
            elseif (dzdt > breform * wmax .and. ineighbour == 1) then
               nh_brfac(i) = 0.0
            else
               nh_brfac(i) = 1.0
            endif
         else                                                       ! was breaking
            if (dzdt < 0.0) then
               nh_brfac(i) = 1.0                                    ! release only when surface falls
            else
               nh_brfac(i) = 0.0
            endif
         endif
         if (nh_brfac(i) < 1.0) pnh(i) = 0.0
      enddo
      !$omp end parallel do
      !
      ! (e) close interior holes (a non-breaking cell flanked on both x-sides, or
      !     both y-sides, by breaking cells is filled in; 3 passes)
      !
      do ipass = 1, 3
         !$acc parallel loop default(present)
         !$omp parallel do schedule ( static ) private ( i )
         do i = 1, nrows
            nh2_r(i) = nh_brfac(i)
         enddo
         !$omp end parallel do
         !$acc parallel loop default(present)
         !$omp parallel do schedule ( static ) private ( i, ineighbour, jl, jr )
         do i = 1, nrows
            if (nh2_r(i) >= 1.0) then
               ineighbour = 0
               jl = nh2_nbw(i) ; if (jl > 0) then ; if (nh2_r(jl) < 1.0) ineighbour = ineighbour + 1 ; endif
               jr = nh2_nbe(i) ; if (jr > 0) then ; if (nh2_r(jr) < 1.0) ineighbour = ineighbour + 2 ; endif
               if (ineighbour == 3) then
                  nh_brfac(i) = 0.0 ; pnh(i) = 0.0
               endif
               if (nmax > 1 .and. nh_brfac(i) >= 1.0) then
                  ineighbour = 0
                  jl = nh2_nbs(i) ; if (jl > 0) then ; if (nh2_r(jl) < 1.0) ineighbour = ineighbour + 1 ; endif
                  jr = nh2_nbn(i) ; if (jr > 0) then ; if (nh2_r(jr) < 1.0) ineighbour = ineighbour + 2 ; endif
                  if (ineighbour == 3) then ; nh_brfac(i) = 0.0 ; pnh(i) = 0.0 ; endif
               endif
            endif
         enddo
         !$omp end parallel do
      enddo
      !
   else
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh_brfac(i) = 1.0
      enddo
      !$omp end parallel do
      !
   endif
   !
   ! 2) Per-face gradient coefficients (wet/dry, layer term abf, conveyance depth,
   !    open-boundary fade-in). Identical to v1 step 2.
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( ip, ipuv, nm, nmu, hu, Dnm1, Dnmu, abf, fdep, bndok )
   do ip = 1, nhuv
      ipuv = nh_faceuv(ip)
      nm   = uv_index_z_nm(ipuv)
      nmu  = uv_index_z_nmu(ipuv)
      nh_cR(ip) = 0.0
      nh_cL(ip) = 0.0
      nh_hu(ip) = 0.0
      ! Open (kcs 2/3) boundary faces are skipped (handled hydrostatically);
      ! closed walls / Neumann keep kfuv = 0 -> correct solid-wall BC.
      bndok = (kcs(nm) == 1 .and. kcs(nmu) == 1)
      if (kfuv(ipuv) == 1 .and. bndok) then
         if (zs(nm) > zb(nm) + huthresh_nh .and. zs(nmu) > zb(nmu) + huthresh_nh) then
            hu = max(zs(nm), zs(nmu)) - 0.5 * (zb(nm) + zb(nmu))
            if (hu > huthresh_nh) then
               Dnm1 = max(zs(nm)  - zb(nm),  huthresh_nh)
               Dnmu = max(zs(nmu) - zb(nmu), huthresh_nh)
               abf  = ( (zs(nmu) + zb(nmu)) - (zs(nm) + zb(nm)) ) / (Dnm1 + Dnmu)
               nh_cR(ip) = nh_cf(ip) * (1.0 + abf)
               nh_cL(ip) = nh_cf(ip) * (abf - 1.0)
               ! Conveyance depth q/uv keeps the corrected q and uv consistent
               if (uv(ipuv) /= 0.0) then
                  nh_hu(ip) = q(ipuv) / uv(ipuv)
               else
                  nh_hu(ip) = hu
               endif
               if (nh_fadein > 0) then
                  fdep = 1.0
                  if (nh_faceL(ip) > 0) fdep = min(fdep, nh_fade(nh_faceL(ip)))
                  if (nh_faceR(ip) > 0) fdep = min(fdep, nh_fade(nh_faceR(ip)))
                  nh_cR(ip) = nh_cR(ip) * fdep
                  nh_cL(ip) = nh_cL(ip) * fdep
               endif
            endif
         endif
      endif
   enddo
   !$omp end parallel do
   !
   ! 3) ASSEMBLY (once per solve). Per row: raw off-diagonals a_ij = dtrho*hu*cR*cL
   !    (one product per incident face, same both ways -> symmetric), the diagonal
   !    dval = dvert + sum dtrho*hu*coef^2, the scaling scl = 1/sqrt(dval), the
   !    breaking-masked scaling sclo, the SCALED right-hand side
   !    b = scl * [ G^T (hu u*) - (ws + wb0 - 2 wb) ], and the warm start
   !    y0 = pnh * sqrt(dval) (scaled previous pressure; 0 at breaking cells).
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i, dval, braw, f, tf, sq )
   do i = 1, nrows
      dval = nh_dvert(i)
      braw = - (ws(i) + wb0(i) - 2.0 * wb(i))
      f = nh_cellface(1, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         nh2_aw(i) = tf * nh_cR(f) * nh_cL(f)
         dval = dval + tf * nh_cR(f)**2
         braw = braw + nh_cR(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh2_aw(i) = 0.0
      endif
      f = nh_cellface(2, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         nh2_ae(i) = tf * nh_cL(f) * nh_cR(f)
         dval = dval + tf * nh_cL(f)**2
         braw = braw + nh_cL(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh2_ae(i) = 0.0
      endif
      f = nh_cellface(3, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         nh2_as(i) = tf * nh_cR(f) * nh_cL(f)
         dval = dval + tf * nh_cR(f)**2
         braw = braw + nh_cR(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh2_as(i) = 0.0
      endif
      f = nh_cellface(4, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         nh2_an(i) = tf * nh_cL(f) * nh_cR(f)
         dval = dval + tf * nh_cL(f)**2
         braw = braw + nh_cL(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh2_an(i) = 0.0
      endif
      sq = sqrt(dval)              ! dval >= dvert > 0 always (Dnm capped)
      nh2_scl(i) = 1.0 / sq
      if (nh_brfac(i) < 1.0) then
         ! breaking cell: zero row/column (via sclo), b=0, y0=0 -> pnh held at 0
         nh2_sclo(i) = 0.0
         nh2_b(i)    = 0.0
         nh2_y(i)    = 0.0
      else
         nh2_sclo(i) = nh2_scl(i)
         nh2_b(i)    = braw * nh2_scl(i)
         nh2_y(i)    = pnh(i) * sq
      endif
   enddo
   !$omp end parallel do
   !
   ! 3b) Scale the off-diagonals: a_ij <- a_ij * sclo(i) * sclo(j). Needs the
   !     completed sclo of the NEIGHBOUR -> separate kernel. sclo(0) = 0 zeroes
   !     coefficients to missing neighbours (the scaled diagonal is exactly 1).
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      nh2_aw(i) = nh2_aw(i) * nh2_sclo(i) * nh2_sclo(nh2_nbw(i))
      nh2_ae(i) = nh2_ae(i) * nh2_sclo(i) * nh2_sclo(nh2_nbe(i))
      nh2_as(i) = nh2_as(i) * nh2_sclo(i) * nh2_sclo(nh2_nbs(i))
      nh2_an(i) = nh2_an(i) * nh2_sclo(i) * nh2_sclo(nh2_nbn(i))
   enddo
   !$omp end parallel do
   !
   ! 4) CG solve of the scaled system (unit diagonal -> no preconditioner work).
   !    Chronopoulos-Gear single-reduction CG: per iteration ONE fused
   !    matvec+dots kernel and ONE fused vector-update kernel.
   !
   ! 4a) initial residual: w = A y0 (branch-free ghost matvec), r = b - w
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      nh2_w(i) = nh2_y(i) + nh2_aw(i) * nh2_y(nh2_nbw(i)) &
                          + nh2_ae(i) * nh2_y(nh2_nbe(i)) &
                          + nh2_as(i) * nh2_y(nh2_nbs(i)) &
                          + nh2_an(i) * nh2_y(nh2_nbn(i))
   enddo
   !$omp end parallel do
   !
   bn8 = 0.0d0
   !$acc parallel loop default(present) reduction(+:bn8)
   !$omp parallel do schedule ( static ) private ( i ) reduction(+:bn8)
   do i = 1, nrows
      nh2_r(i) = nh2_b(i) - nh2_w(i)
      nh2_p(i) = 0.0
      nh2_s(i) = 0.0
      bn8 = bn8 + nh2_b(i) * nh2_b(i)
   enddo
   !$omp end parallel do
   !
   bnorm8 = sqrt(bn8)
   if (bnorm8 <= 0.0d0) bnorm8 = 1.0d0
   !
   ! 4b) iterate. Convergence on the scaled residual: ||S r|| / ||S b|| < nh_tol.
   !
   iter    = 0
   gamold8 = 1.0d0
   alpha8  = 1.0d0
   !
   do
      !
      ! kernel B: w = A r fused with gam = <r,r>, del = <r,w> (real*8 accumulators)
      !
      gam8 = 0.0d0
      del8 = 0.0d0
      !$acc parallel loop default(present) reduction(+:gam8,del8)
      !$omp parallel do schedule ( static ) private ( i, accv ) reduction(+:gam8,del8)
      do i = 1, nrows
         accv = nh2_r(i) + nh2_aw(i) * nh2_r(nh2_nbw(i)) &
                         + nh2_ae(i) * nh2_r(nh2_nbe(i)) &
                         + nh2_as(i) * nh2_r(nh2_nbs(i)) &
                         + nh2_an(i) * nh2_r(nh2_nbn(i))
         nh2_w(i) = accv
         gam8 = gam8 + nh2_r(i) * nh2_r(i)
         del8 = del8 + nh2_r(i) * accv
      enddo
      !$omp end parallel do
      !
      if (sqrt(gam8) < nh_tol * bnorm8) exit
      if (iter >= nh_itermax) exit
      !
      ! scalar recurrences (Chronopoulos-Gear): beta, pAp, alpha
      !
      if (iter == 0) then
         beta8 = 0.0d0
         pap8  = del8
      else
         beta8 = gam8 / gamold8
         pap8  = del8 - beta8 * gam8 / alpha8     ! = <p_new, A p_new>
      endif
      if (pap8 <= 0.0d0) exit                     ! breakdown guard (SPD -> pAp > 0)
      alpha8  = gam8 / pap8
      gamold8 = gam8
      alpha   = real(alpha8, 4)
      beta    = real(beta8, 4)
      !
      ! kernel A: p = r + beta p ; s = w + beta s ; y += alpha p ; r -= alpha s
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh2_p(i) = nh2_r(i) + beta * nh2_p(i)
         nh2_s(i) = nh2_w(i) + beta * nh2_s(i)
         nh2_y(i) = nh2_y(i) + alpha * nh2_p(i)
         nh2_r(i) = nh2_r(i) - alpha * nh2_s(i)
      enddo
      !$omp end parallel do
      !
      iter = iter + 1
      !
   enddo
   !
   ! diagnostics (shared with v1; reported in the end-of-run timing summary)
   !
   nh_iter_total  = nh_iter_total + iter
   nh_solve_count = nh_solve_count + 1
   if (iter > nh_iter_max) nh_iter_max = iter
   !
   ! 4c) unscale: pnh = S y (breaking rows have y = 0 -> pnh = 0 held exactly)
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      pnh(i) = nh2_scl(i) * nh2_y(i)
   enddo
   !$omp end parallel do
   !
   ! 5b) Optional global 2dx filter on pnh (nh_filter > 0): one Jacobi
   !     neighbour-mean pass; snapshot into nh2_w first (race-free).
   !
   if (nh_filter > 0.0) then
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh2_w(i) = pnh(i)
      enddo
      !$omp end parallel do
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, sumn, cnt, j )
      do i = 1, nrows
         sumn = 0.0 ; cnt = 0
         j = nh2_nbw(i) ; if (j > 0) then ; sumn = sumn + nh2_w(j) ; cnt = cnt + 1 ; endif
         j = nh2_nbe(i) ; if (j > 0) then ; sumn = sumn + nh2_w(j) ; cnt = cnt + 1 ; endif
         j = nh2_nbs(i) ; if (j > 0) then ; sumn = sumn + nh2_w(j) ; cnt = cnt + 1 ; endif
         j = nh2_nbn(i) ; if (j > 0) then ; sumn = sumn + nh2_w(j) ; cnt = cnt + 1 ; endif
         if (cnt > 0) pnh(i) = (1.0 - nh_filter) * nh2_w(i) + nh_filter * sumn / real(cnt)
      enddo
      !$omp end parallel do
   endif
   !
   ! 5c) Localized 2dx smoothing in the marginal-nonh zones (fade-in zone and/or
   !     shallow run-up water); zero in the resolved interior. Same as v1.
   !
   if (nh_smoothbnd > 0.0 .and. (nh_fadein > 0 .or. nh_smoothdep > 0.0)) then
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh2_w(i) = pnh(i)
      enddo
      !$omp end parallel do
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, sumn, cnt, j, fdep, gf )
      do i = 1, nrows
         sumn = 0.0 ; cnt = 0
         j = nh2_nbw(i) ; if (j > 0) then ; sumn = sumn + nh2_w(j) ; cnt = cnt + 1 ; endif
         j = nh2_nbe(i) ; if (j > 0) then ; sumn = sumn + nh2_w(j) ; cnt = cnt + 1 ; endif
         j = nh2_nbs(i) ; if (j > 0) then ; sumn = sumn + nh2_w(j) ; cnt = cnt + 1 ; endif
         j = nh2_nbn(i) ; if (j > 0) then ; sumn = sumn + nh2_w(j) ; cnt = cnt + 1 ; endif
         if (cnt > 0) then
            fdep = 0.0
            if (nh_fadein > 0)      fdep = 1.0 - nh_fade(i)
            if (nh_smoothdep > 0.0) fdep = max(fdep, (nh_smoothdep - Dnm(i)) / nh_smoothdep)
            gf = nh_smoothbnd * max(0.0, min(1.0, fdep))
            pnh(i) = (1.0 - gf) * nh2_w(i) + gf * sumn / real(cnt)
         endif
      enddo
      !$omp end parallel do
   endif
   !
   ! 5e) Depth limiter: |pnh| <= nh_pmax * rho g H (clips wall/breaking/2dx spikes).
   !
   if (nh_pmax > 0.0) then
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, pcap )
      do i = 1, nrows
         pcap = nh_pmax * rhow * g * Dnm(i)
         pnh(i) = max(-pcap, min(pcap, pnh(i)))
      enddo
      !$omp end parallel do
   endif
   !
   ! 6) Correct fluxes / velocities with the pressure gradient -(dt/rho) G pnh.
   !    Each face writes only its own uv/q point -> race-free.
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( ip, ipuv, pL, pR, gf, unh )
   do ip = 1, nhuv
      if (nh_cR(ip) == 0.0 .and. nh_cL(ip) == 0.0) cycle
      ipuv = nh_faceuv(ip)
      pR = 0.0
      pL = 0.0
      if (nh_faceR(ip) > 0) pR = pnh(nh_faceR(ip))
      if (nh_faceL(ip) > 0) pL = pnh(nh_faceL(ip))
      gf  = nh_cR(ip) * pR + nh_cL(ip) * pL
      unh = - dtrho * gf
      uv(ipuv) = (1.0 - nh_fnudge) * uv(ipuv) + nh_fnudge * (uv(ipuv) + unh)
      q(ipuv)  = (1.0 - nh_fnudge) * q(ipuv) + nh_fnudge * (q(ipuv) + nh_hu(ip) * unh)
   enddo
   !$omp end parallel do
   !
   ! 7) Update surface/bottom vertical velocities for the next step's forcing
   !    (bottom kinematic condition w_b = u . d(zb)/dx with the slopes FROZEN at
   !    initialization in nh_slbed -- see v1; the wet checks use the current,
   !    possibly effective, zb). Identical to v1 step 7.
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i, nm, iuv, nmn )
   do i = 1, nrows
      nm = nm_index_of_row(i)
      !
      ! Non-hydrostatically dry cell: carry no nonh state (clean re-wetting).
      !
      if (zs(nm) - zb(nm) <= huthresh_nh) then
         wb0(i) = 0.0
         wb(i)  = 0.0
         ws(i)  = 0.0
         pnh(i) = 0.0
         nh_brfac(i) = 1.0
         cycle
      endif
      !
      wb0(i) = wb(i)
      wb(i)  = 0.0
      !
      ! w_b = u . d(zb)/dx over faces to non-hydrostatically wet neighbours only.
      !
      iuv = z_index_uv_md(nm)
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nm(iuv)
         if (zs(nmn) - zb(nmn) > huthresh_nh) then
            wb(i) = wb(i) - 0.5 * uv(iuv) * nh_slbed(1, i)
         endif
      endif
      iuv = z_index_uv_mu(nm)
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nmu(iuv)
         if (zs(nmn) - zb(nmn) > huthresh_nh) then
            wb(i) = wb(i) - 0.5 * uv(iuv) * nh_slbed(2, i)
         endif
      endif
      if (nmax > 1) then
         iuv = z_index_uv_nd(nm)
         if (kfuv(iuv) > 0) then
            nmn = uv_index_z_nm(iuv)
            if (zs(nmn) - zb(nmn) > huthresh_nh) then
               wb(i) = wb(i) - 0.5 * uv(iuv) * nh_slbed(3, i)
            endif
         endif
         iuv = z_index_uv_nu(nm)
         if (kfuv(iuv) > 0) then
            nmn = uv_index_z_nmu(iuv)
            if (zs(nmn) - zb(nmn) > huthresh_nh) then
               wb(i) = wb(i) - 0.5 * uv(iuv) * nh_slbed(4, i)
            endif
         endif
      endif
      !
      ws(i) = ws(i) - (wb(i) - wb0(i)) + (kbfac * dt / (rhow * Dnm(i))) * pnh(i)
   enddo
   !$omp end parallel do
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine

end module
