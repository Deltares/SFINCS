! Non-hydrostatic code now only works with regular grids (can still use quadtree netcdf file as long as there are no refinement levels).
! Now uses bicgstab_ilu to solve matrix. Both should ideally utilize CPU and GPU parallelization. Currently, this solver cannot be fully parallelized. 
! Should try to find a solver that can be.   

module sfincs_nonhydrostatic
   !   
   integer, dimension(:,:), allocatable :: index_sparse_matrix
   integer, dimension(:,:), allocatable :: nh_uv_index
   integer, dimension(:), allocatable   :: nm_index_of_row
   integer, dimension(:), allocatable   :: row_index_of_nm
   integer, dimension(:), allocatable   :: uv_index_of_nhuv
   integer, dimension(:), allocatable   :: col_idx
   integer, dimension(:), allocatable   :: row_ptr
   !
   real*4, dimension(:),   allocatable  :: pnh
   real*4, dimension(:),   allocatable  :: ws
   real*4, dimension(:),   allocatable  :: wb
   real*4, dimension(:),   allocatable  :: wb0
   !
   real*4,  dimension(:), allocatable   :: Dnm
   real*4,  dimension(:), allocatable   :: dzbdx
   real*4,  dimension(:), allocatable   :: dzbdy
   !
   ! Work arrays for the matrix-free SPD conjugate-gradient solver (nonh_solver = 2)
   !
   real*4,  dimension(:), allocatable   :: nh_dvert   ! vertical-accel diagonal term   (nrows)
   real*4,  dimension(:), allocatable   :: nh_diag    ! Jacobi diagonal of A           (nrows)
   real*4,  dimension(:), allocatable   :: nh_sponge_ramp ! open-boundary sponge ramp 0..1 (nrows)
   real*4,  dimension(:), allocatable   :: nh_spongediag  ! sponge diagonal contribution   (nrows)
   real*4,  dimension(:), allocatable   :: nh_fade        ! open-boundary nonh fade-in 0..1 (nrows)
   real*4,  dimension(:), allocatable   :: nh_brfac       ! HFA breaking factor 1=full nonh .. 0=hydrostatic (nrows)
   integer, dimension(:), allocatable   :: nh_isfront     ! 1 if cell has a dry neighbour  (nrows)
   real*4,  dimension(:), allocatable   :: cg_r       ! CG residual                    (nrows)
   real*4,  dimension(:), allocatable   :: cg_z       ! CG preconditioned residual     (nrows)
   real*4,  dimension(:), allocatable   :: cg_d       ! CG search direction            (nrows)
   real*4,  dimension(:), allocatable   :: cg_Ap      ! CG operator-times-direction    (nrows)
   !
   integer, dimension(:), allocatable   :: nh_faceuv  ! full uv index of each nh face  (nhuv)
   integer, dimension(:), allocatable   :: nh_faceL   ! row index of left/bottom cell  (nhuv, 0 = boundary)
   integer, dimension(:), allocatable   :: nh_faceR   ! row index of right/top cell    (nhuv, 0 = boundary)
   real*4,  dimension(:), allocatable   :: nh_cf      ! 0.5/dx static gradient weight  (nhuv)
   real*4,  dimension(:), allocatable   :: nh_cR      ! gradient coeff of p(R) per face(nhuv)
   real*4,  dimension(:), allocatable   :: nh_cL      ! gradient coeff of p(L) per face(nhuv)
   real*4,  dimension(:), allocatable   :: nh_hu      ! water depth at face this step  (nhuv)
   !
   real*4    :: huthresh_nh
   !
   integer   :: nrows
   integer   :: nr_vals_in_matrix
   integer   :: nhuv
   !
contains
   !
   subroutine initialize_nonhydrostatic()
   !
   ! Initialization of pardiso solver
   !
   use sfincs_data
   !
   implicit none
   !
   integer nm, k, nmu, nmd, num, ndm, irow, icol, inb, ip
   integer inhuv, iii
   !
   ! Temporary arrays
   !
   integer, dimension(:), allocatable   :: col_idx0
   integer, dimension(:,:), allocatable :: nh_nm_index
   !
   allocate(row_index_of_nm(np))
   !
   row_index_of_nm = 0
   !
   ! Count number of rows (= nr cols) in matrix
   !
   nrows = 0
   !
   do nm = 1, np
      !
      if (mask_nonh(nm) == 1) then
         !
         nrows = nrows + 1
         !
      endif
      !
   enddo
   !
   allocate(pnh(nrows))
   allocate(ws(nrows))
   allocate(wb(nrows))
   allocate(wb0(nrows))
   allocate(index_sparse_matrix(5, nrows))
   allocate(row_ptr(nrows + 1))
   allocate(col_idx0(5 * nrows))
   allocate(Dnm(nrows))
   allocate(dzbdx(nrows))
   allocate(dzbdy(nrows))
   allocate(nm_index_of_row(nrows))
   allocate(nh_uv_index(4, nrows))
   allocate(nh_nm_index(4, nrows))
   !
   ! Work arrays for the matrix-free CG solver (nonh_solver = 2)
   !
   allocate(nh_dvert(nrows))
   allocate(nh_diag(nrows))
   allocate(nh_sponge_ramp(nrows))
   allocate(nh_spongediag(nrows))
   allocate(nh_fade(nrows))
   nh_fade = 1.0
   allocate(nh_brfac(nrows))
   nh_brfac = 1.0
   allocate(nh_isfront(nrows))
   nh_isfront = 0
   allocate(cg_r(nrows))
   allocate(cg_z(nrows))
   allocate(cg_d(nrows))
   allocate(cg_Ap(nrows))
   nh_dvert = 0.0
   nh_diag = 0.0
   nh_sponge_ramp = 0.0
   nh_spongediag  = 0.0
   cg_r = 0.0
   cg_z = 0.0
   cg_d = 0.0
   cg_Ap = 0.0
   !
   huthresh_nh = max(huthresh, 0.01)
   !
   pnh = 0.0
   ws  = 0.0       
   wb  = 0.0       
   wb0 = 0.0       
   !
   index_sparse_matrix = 0
   row_ptr = 0
   col_idx = 0
   col_idx0 = 0
   irow = 0
   k = 0
   Dnm = 0.0
   dzbdx = 0.0
   dzbdy = 0.0
   nm_index_of_row = 0
   nh_uv_index = 0
   nh_nm_index = 0
   ! bnd = 0
   !
   ! Map row indices to nm indices and vice versa
   !
   irow = 0
   !
   do nm = 1, np
      if (mask_nonh(nm) == 1) then
         irow = irow + 1
         row_index_of_nm(nm) = irow
         nm_index_of_row(irow) = nm
      endif
   enddo      
   !
   ! Find velocity points needed for nh computations
   !
   !               4
   !        +------|------+
   !        |             | 
   !        |       irow  | 
   !     1  -      + 5    -  2
   !        |             | 
   !        |             | 
   !        +------|------+
   !               3
   !
   ! Loop through all uv points to get nh uv points (any uv point that touches nh cell)
   !
   ! First just count them
   !
   nhuv = 0
   !
   do ip = 1, npuv
      !
      if (mask_nonh(uv_index_z_nm(ip)) == 1 .or. mask_nonh(uv_index_z_nmu(ip)) == 1) then
         !
         nhuv = nhuv + 1
         !
      endif
      !
   enddo
   !
   allocate(uv_index_of_nhuv(nhuv))
   uv_index_of_nhuv = 0
   !
   ! Face arrays for the matrix-free CG solver (nonh_solver = 2)
   !
   allocate(nh_faceuv(nhuv))
   allocate(nh_faceL(nhuv))
   allocate(nh_faceR(nhuv))
   allocate(nh_cf(nhuv))
   allocate(nh_cR(nhuv))
   allocate(nh_cL(nhuv))
   allocate(nh_hu(nhuv))
   nh_faceuv = 0
   nh_faceL  = 0
   nh_faceR  = 0
   nh_cf     = 0.0
   nh_cR     = 0.0
   nh_cL     = 0.0
   nh_hu     = 0.0
   !
   ! Now find neigboring indices of nh uv points
   !
   inhuv = 0
   !
   do ip = 1, npuv
      !
      nm  = uv_index_z_nm(ip)
      nmu = uv_index_z_nmu(ip)
      !
      if (mask_nonh(nm) == 1 .or. mask_nonh(nmu) == 1) then
         !
         inhuv = inhuv + 1
         !
         uv_index_of_nhuv(inhuv) = ip ! uv index of this nhuv point
         !
         ! Get direction of point
         !
         if (uv_flags_dir(ip) == 0) then
            !
            ! x
            !
            if (mask_nonh(nm) == 1) then
               !
               ! nm point
               !
               irow = row_index_of_nm(nm)
               inb  = row_index_of_nm(nmu)
               !
               nh_uv_index(2, irow) = inhuv
               nh_nm_index(2, irow) = inb
               !
            endif   
            !
            if (mask_nonh(nmu) == 1) then
               !
               ! nmu point
               !
               inb  = row_index_of_nm(nm)
               irow = row_index_of_nm(nmu)
               !
               nh_uv_index(1, irow) = inhuv
               nh_nm_index(1, irow) = inb
               !
            endif   
            !
         else
            !
            ! y
            !
            if (mask_nonh(nm) == 1) then
               !
               ! nm point
               !
               irow = row_index_of_nm(nm)
               inb  = row_index_of_nm(nmu)
               !
               nh_uv_index(4, irow) = inhuv
               nh_nm_index(4, irow) = inb
               !
            endif   
            !
            if (mask_nonh(nmu) == 1) then
               !
               ! nmu point
               !
               inb  = row_index_of_nm(nm)
               irow = row_index_of_nm(nmu)
               !
               nh_uv_index(3, irow) = inhuv
               nh_nm_index(3, irow) = inb
               !
            endif   
            !
         endif            
         !
      endif
      !
   enddo   
   !
   ! Determine indices of sparse matrix AA
   !
   iii = 0
   do irow = 1, nrows
      !
      !      if (irow <= iii) then
      !         bnd(icol) = 1
      !      endif
      !
      nm = nm_index_of_row(irow)
      !
      ! Left
      !
      icol = nh_nm_index(1, irow)
      !
      if (icol > iii) then
         !
         ! Cell has a neighbor to the left
         !
         k = k + 1
         col_idx0(k) = icol
         index_sparse_matrix(1, irow) = k
         !
         if (row_ptr(irow) == 0) then ! first data in row
            !
            row_ptr(irow) = k
            !
         endif   
         !
      endif               
      !
      ! Bottom
      !
      icol = nh_nm_index(3, irow)
      !
      if (icol > iii) then
         !
         ! Cell has a neighbor below
         !
         k = k + 1
         col_idx0(k) = icol
         index_sparse_matrix(3, irow) = k
         !
         if (row_ptr(irow) == 0) then ! first data in row
            !
            row_ptr(irow) = k
            !
         endif   
         !
      endif
      !
      ! Centre
      !
      icol = irow
      !
      ! if (icol > iii) then
      !
      k = k + 1
      col_idx0(k) = icol
      index_sparse_matrix(5, irow) = k
      !
      if (row_ptr(irow) == 0) then ! first data in row
         !
         row_ptr(irow) = k
         !
      endif   
      !   
      !  endif
      !
      ! Top
      !
      icol = nh_nm_index(4, irow)
      !
      if (icol > iii) then
         !
         ! Cell has a neighbor above
         !
         k = k + 1
         col_idx0(k) = icol
         index_sparse_matrix(4, irow) = k
         !
         if (row_ptr(irow) == 0) then ! first data in row
            !
            row_ptr(irow) = k
            !
         endif   
         !
      endif
      !
      ! Right
      !
      icol = nh_nm_index(2, irow)
      !
      if (icol > iii) then
         !
         ! Cell has a neighbor to the right
         !
         k = k + 1
         col_idx0(k) = icol
         index_sparse_matrix(2, irow) = k
         !
         if (row_ptr(irow) == 0) then ! first data in row
            !
            row_ptr(irow) = k
            !
         endif   
         !
      endif
      !
   enddo
   !
   nr_vals_in_matrix = k
   allocate(col_idx(nr_vals_in_matrix))
   col_idx(1:k) = col_idx0(1:k)
   row_ptr(nrows + 1) = k + 1
   !
   ! Compute bed level slopes
   !
   do irow = 1, nrows
      !
      nm = nm_index_of_row(irow)
      !
      nmd = 0
      nmu = 0
      ndm = 0
      num = 0
      !
      ! Left
      !
      if (nh_nm_index(1, irow) > 0) then
         !
         nmd = nm_index_of_row(nh_nm_index(1, irow))
         dzbdx(irow) = dzbdx(irow) + 0.5 * (zb(nm) - zb(nmd)) * dxrinv(1)
         !
      endif   
      !
      ! Right
      !
      if (nh_nm_index(2, irow) > 0) then
         !
         nmu = nm_index_of_row(nh_nm_index(2, irow))
         dzbdx(irow) = dzbdx(irow) + 0.5 * (zb(nmu) - zb(nm)) * dxrinv(1)
         !
      endif
      !
      ! Below
      !
      if (nh_nm_index(3, irow) > 0) then
         !
         ndm = nm_index_of_row(nh_nm_index(3, irow))
         dzbdy(irow) = dzbdy(irow) + 0.5 * (zb(nm) - zb(ndm)) * dxrinv(1)
         !
      endif
      !
      ! Above
      !
      if (nh_nm_index(4, irow) > 0) then
         !
         num = nm_index_of_row(nh_nm_index(4, irow))
         dzbdy(irow) = dzbdy(irow) + 0.5 * (zb(num) - zb(nm)) * dxrinv(1)
         !
      endif
      !
   enddo
   !
   ! Static per-face data for the matrix-free CG solver (nonh_solver = 2).
   ! nh_cf = 0.5 / dx is the (depth-averaged) gradient weight; the factor 0.5
   ! reflects that the linearly-varying non-hydrostatic pressure has a depth
   ! mean of p_bed/2. Single refinement level here, so dxrinv(1) / dyrinv(1).
   !
   do ip = 1, nhuv
      !
      inhuv = uv_index_of_nhuv(ip)
      nh_faceuv(ip) = inhuv
      nh_faceL(ip)  = row_index_of_nm(uv_index_z_nm(inhuv))
      nh_faceR(ip)  = row_index_of_nm(uv_index_z_nmu(inhuv))
      !
      if (uv_flags_dir(inhuv) == 0) then
         nh_cf(ip) = 0.5 * dxrinv(1)
      else
         nh_cf(ip) = 0.5 * dyrinv(1)
      endif
      !
   enddo
   !
   ! Open-boundary non-hydrostatic sponge ramp (nonh_solver = 2). A relaxation
   ! layer nh_sponge_width cells deep inside each water-level / outflow (kcs 2/3)
   ! boundary, where a ramped term is added to the pressure-operator diagonal so
   ! that pnh is smoothly driven to zero towards the boundary. This absorbs the
   ! outgoing non-hydrostatic (dispersive) energy that would otherwise resonate.
   ! The cell-distance to the nearest open boundary is found by Bellman-Ford
   ! relaxation over the non-hydrostatic face graph.
   !
   nh_sponge_ramp = 0.0
   !
   if (nh_sponge_width > 0) then
      !
      block
         integer, dimension(:), allocatable :: idist
         integer :: pass, nmL, nmR, rL, rR, d
         real*4  :: rampval
         !
         allocate(idist(nrows))
         idist = nh_sponge_width + 1
         !
         ! Seed: nonh cells that border an open boundary cell (distance 0)
         !
         do ip = 1, nhuv
            inhuv = uv_index_of_nhuv(ip)
            nmL = uv_index_z_nm(inhuv)
            nmR = uv_index_z_nmu(inhuv)
            rL  = nh_faceL(ip)
            rR  = nh_faceR(ip)
            if (rL > 0 .and. rR == 0) then
               if (kcs(nmR) == 2 .or. kcs(nmR) == 3) idist(rL) = 0
            endif
            if (rR > 0 .and. rL == 0) then
               if (kcs(nmL) == 2 .or. kcs(nmL) == 3) idist(rR) = 0
            endif
         enddo
         !
         ! Propagate distance inward (nh_sponge_width passes suffice)
         !
         do pass = 1, nh_sponge_width
            do ip = 1, nhuv
               rL = nh_faceL(ip)
               rR = nh_faceR(ip)
               if (rL > 0 .and. rR > 0) then
                  if (idist(rL) + 1 < idist(rR)) idist(rR) = idist(rL) + 1
                  if (idist(rR) + 1 < idist(rL)) idist(rL) = idist(rR) + 1
               endif
            enddo
         enddo
         !
         ! Smooth (quadratic) ramp: 1 at the boundary, 0 at the inner edge
         !
         do irow = 1, nrows
            d = idist(irow)
            if (d < nh_sponge_width) then
               rampval = real(nh_sponge_width - d) / real(nh_sponge_width)
               nh_sponge_ramp(irow) = rampval * rampval
            endif
         enddo
         !
         deallocate(idist)
      end block
      !
   endif
   !
   ! Open-boundary non-hydrostatic FADE-IN (nonh_solver = 2). Over nh_fadein cells
   ! inside each open (water level / outflow) boundary, ramp the non-hydrostatic
   ! coupling from 0 (at the boundary) to full (nh_fadein cells in). The incident
   ! wave then enters hydrostatically at the boundary and the nonh turns on
   ! gradually, avoiding the sharp pnh=0 -> full-nonh transition that reflects an
   ! incoming dispersive wave. nh_fade is 1 everywhere by default (no fade).
   !
   if (nh_fadein > 0) then
      !
      block
         integer, dimension(:), allocatable :: idist
         integer :: pass, nmL, nmR, rL, rR
         !
         allocate(idist(nrows))
         idist = nh_fadein + 1
         do ip = 1, nhuv
            inhuv = uv_index_of_nhuv(ip)
            nmL = uv_index_z_nm(inhuv)
            nmR = uv_index_z_nmu(inhuv)
            rL  = nh_faceL(ip)
            rR  = nh_faceR(ip)
            if (rL > 0 .and. rR == 0) then
               if (kcs(nmR) == 2 .or. kcs(nmR) == 3) idist(rL) = 0
            endif
            if (rR > 0 .and. rL == 0) then
               if (kcs(nmL) == 2 .or. kcs(nmL) == 3) idist(rR) = 0
            endif
         enddo
         do pass = 1, nh_fadein
            do ip = 1, nhuv
               rL = nh_faceL(ip)
               rR = nh_faceR(ip)
               if (rL > 0 .and. rR > 0) then
                  if (idist(rL) + 1 < idist(rR)) idist(rR) = idist(rL) + 1
                  if (idist(rR) + 1 < idist(rL)) idist(rL) = idist(rR) + 1
               endif
            enddo
         enddo
         do irow = 1, nrows
            nh_fade(irow) = min(1.0, real(idist(irow)) / real(nh_fadein))
         enddo
         deallocate(idist)
      end block
      !
   endif
   !
   end subroutine


   subroutine compute_nonhydrostatic(dt, tloop)
   !
   ! Non-hydrostatic pressure correction on fluxes and velocities 
   !
   use sfincs_data
   use bicgstab_solver_ilu
   !
   implicit none
   !
   integer   :: count0
   integer   :: count1
   integer   :: count_rate
   integer   :: count_max
   real      :: tloop
   !
   real*4    :: dt
   !
   integer   :: ip
   integer   :: ipuv
   integer   :: nm
   integer   :: nmu
   integer   :: nmd
   integer   :: num
   integer   :: ndm
   integer   :: nmn
   integer   :: j
   integer   :: irow
   integer   :: nhnm
   integer   :: nhnmu
   integer   :: iuv
   !
   integer   :: iter
   !
   real*4    :: hu
   !
   real*4    :: Dnm1
   real*4    :: hnm
   real*4    :: Dnmu
   real*4    :: hnmu
   real*4    :: unh
   real*4    :: pnhnm
   real*4    :: pnhnmu
   !
   real*4    :: hnb
   !
   real*4    :: dtover2rhodx2
   real*4    :: dtover2rhodx
   !
   real*4, dimension(npuv)            :: AB
   real*4, dimension(:), allocatable  :: QQ
   real*4, dimension(:), allocatable  :: AA
   real*4                             :: relres
   !
   call system_clock(count0, count_rate, count_max)
   !
   allocate(QQ(nrows))
   allocate(AA(nr_vals_in_matrix))
   !
   ! Compute AB (A and B)
   !
   AB = 0.0
   !
   !$omp parallel &
   !$omp private ( nm, irow )
   !$omp do schedule ( dynamic, 256 )
   do irow = 1, nrows
      !
      nm = nm_index_of_row(irow)
      !
      Dnm(irow)  = max(zs(nm) - zb(nm), huthresh_nh)
      !
   enddo
   !$omp end do
   !$omp end parallel
   !   
   !$omp parallel &
   !$omp private ( ip, ipuv, nm, nmu, num, hnm, hnmu, Dnm1, Dnmu )
   !$omp do schedule ( dynamic, 256 )
   do ip = 1, nhuv
      !
      ! Get water levels of neighboring cells
      !
      ! First get index of 'complete' uv array
      !
      ipuv = uv_index_of_nhuv(ip)
      !
      nm   = uv_index_z_nm(ipuv)
      nmu  = uv_index_z_nmu(ipuv)
      !
      if (kfuv(ipuv) == 1) then
         if (zs(nm) > zb(nm) + huthresh_nh .and. zs(nmu) > zb(nmu) + huthresh_nh) then
            !
            hnm  = - zb(nm)
            Dnm1 = max(zs(nm) - zb(nm), huthresh_nh)
            !
            hnmu = - zb(nmu)
            Dnmu = max(zs(nmu) - zb(nmu), huthresh_nh)
            !          
            AB(ip) = ( (zs(nmu) - hnmu) - (zs(nm) - hnm) ) / (Dnm1 + Dnmu)
            !
         endif
      endif   
      !
   enddo
   !$omp end do
   !$omp end parallel
   !
   ! Compute non-hydrostatic pressure by solving matrix AA * PP = QQ, where AA is a sparse matrix, PP is the nonh pressure, and QQ is the forcing
   !
   AA = 0.0
   QQ = 0.0
   !
   dtover2rhodx2 = (dt * dxr2inv(1) / (2 * rhow))
   !
   ! Fill sparse matrix
   !
   !$omp parallel &
   !$omp private ( ip, nm, j, nmd, nmu, ndm, num )
   !$omp do schedule ( dynamic, 256 )
   do irow = 1, nrows
      !
      ! Indices in nh uv array of neighboring uv points
      !
      nmd = nh_uv_index(1, irow)
      nmu = nh_uv_index(2, irow)
      ndm = nh_uv_index(3, irow)
      num = nh_uv_index(4, irow)
      !
      ! Left
      !
      if (nmd > 0) then
         !
         j = index_sparse_matrix(1, irow)
         !
         if (j>0) then
            !
            AA(j) = dtover2rhodx2 * (-1.0 + AB(nmd))
            !
         endif
         !
      endif
      !
      ! Right
      !
      if (nmu > 0) then
         !
         j = index_sparse_matrix(2, irow)
         !
         if (j>0) then
            !
            AA(j) = dtover2rhodx2 * (-1.0 - AB(nmu))            
            !
         endif
         !
      endif
      !
      ! Bottom
      !
      if (ndm > 0) then
         !
         j = index_sparse_matrix(3, irow)
         !
         if (j>0) then
            !
            AA(j) = dtover2rhodx2 * (-1.0 + AB(ndm))
            !
         endif
         !
      endif
      !
      ! Top
      !
      if (num > 0) then
         !
         j = index_sparse_matrix(4, irow)
         !
         if (j>0) then
            !
            AA(j) = dtover2rhodx2 * (-1.0 - AB(num))            
            !
         endif
         !
      endif
      !
      ! Centre
      !
      j = index_sparse_matrix(5, irow)
      !
      ! if (j>0) then
      !
      AA(j) = 2 * dt / ( rhow * Dnm(irow)**2 )
      !
      ! endif
      !
      ! Forcing
      !
      QQ(irow) = 0.0
      !
      ! if (bnd(irow) == 0) then   
      !         
      QQ(irow) = - (ws(irow) + wb0(irow) - 2 * wb(irow)) / Dnm(irow)
      !
      if (nmd > 0) then
         !
         AA(j) = AA(j) + dtover2rhodx2 * (1.0 + AB(nmd))
         !
         QQ(irow) = QQ(irow) - (- uv(uv_index_of_nhuv(nmd))) * dxrinv(1)
         !
      endif
      !
      if (nmu > 0) then
         !
         AA(j) = AA(j) + dtover2rhodx2 * (1.0 - AB(nmu))
         !
         QQ(irow) = QQ(irow) - (uv(uv_index_of_nhuv(nmu))) * dxrinv(1)
         !
      endif
      !
      if (ndm > 0) then
         !
         AA(j) = AA(j) + dtover2rhodx2 * (1.0 + AB(ndm))
         !
         QQ(irow) = QQ(irow) - (- uv(uv_index_of_nhuv(ndm))) * dxrinv(1)
         !
      endif
      !
      if (num > 0) then
         !
         AA(j) = AA(j) + dtover2rhodx2 * (1.0 - AB(num))
         !
         QQ(irow) = QQ(irow) - (uv(uv_index_of_nhuv(num))) * dxrinv(1)
         !
      endif
      !
      ! endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !
   ! Solve matrix
   !
   call bicgstab_solve(nrows, AA, col_idx, row_ptr, QQ, pnh, nh_tol, nh_itermax, iter, relres, .true.)
   !
   ! Adjust fluxes   
   !
   dtover2rhodx = (dt * dxrinv(1) / (2 * rhow))
   !
   !$omp parallel &
   !$omp private ( ip, ipuv, nm, nmu, nhnm, nhnmu, pnhnm, pnhnmu, hu, unh )
   !$omp do schedule ( dynamic, 256 )
   do ip = 1, nhuv
      !
      ipuv = uv_index_of_nhuv(ip)
      !
      if (kfuv(ipuv) == 1) then
         !
         ! Indices of neighbors in full zs array
         !
         nm    = uv_index_z_nm(ipuv)
         nmu   = uv_index_z_nmu(ipuv)
         !
         ! Indices of neighbors in nonh zs array
         !
         nhnm  = row_index_of_nm(nm)
         nhnmu = row_index_of_nm(nmu)
         !
         ! Non-hydrostatic pressure at the two neighbours. A neighbour that is
         ! not itself a non-hydrostatic cell (row index 0) is an open boundary
         ! (water level / outflow) where the flow is hydrostatic: impose the
         ! Dirichlet condition pnh = 0 there. This is consistent with the matrix
         ! assembly, which already adds the centre-diagonal contribution for such
         ! a face without an off-diagonal term, and it removes the out-of-bounds
         ! pnh(0) read that previously occurred at the edge of the nonh region.
         !
         pnhnm  = 0.0
         pnhnmu = 0.0
         if (nhnm  > 0) pnhnm  = pnh(nhnm)
         if (nhnmu > 0) pnhnmu = pnh(nhnmu)
         !
         !hu = max(zs(nm), zs(nmu)) - 0.5 * (zb(nm) + zb(nmu))
         hu = max(zs(nm), zs(nmu)) - max(zb(nm), zb(nmu))
         !
         if (hu > huthresh_nh) then
            !
            unh = - 1.0 * dtover2rhodx * ( AB(ip) * (pnhnmu + pnhnm) + pnhnmu - pnhnm )
            !
            ! Do some nudging to avoid 2dx waves
            !
            q(ipuv)  = (1.0 - nh_fnudge) * q(ipuv) + nh_fnudge * (q(ipuv) + hu * unh)
            uv(ipuv) = (1.0 - nh_fnudge) * uv(ipuv) + nh_fnudge * (uv(ipuv) + unh)
            !
         endif
         !
      endif
      !
   enddo   
   !$omp end do
   !$omp end parallel   
   !
   ! Update vertical velocity ws and wb
   !
   !$omp parallel &
   !$omp private ( nm, iuv, irow, nmn, hnm, hnb)
   !$omp do schedule ( dynamic, 256 )
   do irow = 1, nrows
      !
      ! Copy wb0 from previous time step
      !
      wb0(irow) = wb(irow) 
      !
      wb(irow) = 0.0
      !
      nm = nm_index_of_row(irow)
      hnm  = zs(nm) - zb(nm)
      !
      ! This will not yet work for quadtree !
      !
      ! Indices of neighboring cells
      !
      ! Left
      !
      iuv = z_index_uv_md(nm)      ! uv index
      !
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nm(iuv)  ! nm index of neighbor
         !hnb = zs(nmn) - zb(nmn)   ! water depth at neighbor
         wb(irow) = wb(irow) - 0.5 * uv(iuv) * (zb(nmn) - zb(nm)) * dxrinv(1)
      endif
      !
      ! Right
      !
      iuv = z_index_uv_mu(nm)      ! uv index
      !
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nmu(iuv) ! nm index of neighbor
         !hnb = zs(nmn) - zb(nmn)   ! water depth at neighbor
         wb(irow) = wb(irow) - 0.5 * uv(iuv) * (zb(nm) - zb(nmn)) * dxrinv(1)
      endif
      !
      ! Bottom and Top contribute only when there is more than one row. For a
      ! 1-D model (nmax == 1) the cell has no y-neighbours and z_index_uv_nd/nu
      ! point at the dummy slot (npuv+ncuv+1), which would index kfuv (size npuv)
      ! out of bounds. The y-direction bed kinematic term is identically zero
      ! there anyway, so skip it.
      !
      if (nmax > 1) then
         !
         ! Bottom
         !
         iuv = z_index_uv_nd(nm)     ! uv index
         !
         if (kfuv(iuv) > 0) then
            nmn = uv_index_z_nm(iuv) ! nm index of neighbor
            !hnb = zs(nmn) - zb(nmn)  ! water depth at neighbor
            wb(irow) = wb(irow) - 0.5 * uv(iuv) * (zb(nmn) - zb(nm)) * dyrinv(1)
         endif
         !
         ! Top
         !
         iuv = z_index_uv_nu(nm)      ! uv index
         !
         if (kfuv(iuv) > 0) then
            nmn = uv_index_z_nmu(iuv) ! nm index of neighbor
            !hnb = zs(nmn) - zb(nmn)   ! water depth at neighbor
            wb(irow) = wb(irow) - 0.5 * uv(iuv) * (zb(nm) - zb(nmn)) * dyrinv(1)
         endif
         !
      endif
      !
      ws(irow) = ws(irow) - (wb(irow) - wb0(irow)) + (2 * dt / (rhow * Dnm(irow))) * pnh(irow) ! this is ws m+1 in the next time step
      !
   enddo   
   !$omp end do
   !$omp end parallel
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine


   subroutine compute_nonhydrostatic2(dt, tloop)
   !
   ! Matrix-free, symmetric-positive-definite non-hydrostatic pressure
   ! projection, solved with Jacobi-preconditioned conjugate gradients.
   !
   ! The pressure-gradient (momentum force) operator G and the divergence
   ! operator are exact negative transposes (D = -G^T), so the pressure operator
   !
   !     A = (dt/rho) G^T G  +  diag( 2 dt /(rho H^2) )
   !
   ! is SPD and CG applies. This replaces the legacy non-symmetric assembly and
   ! the sequential BiCGSTAB-ILU solve, and needs no anti-2dx nudging because the
   ! pressure gradient shares the momentum stencil. Single refinement level
   ! (quadtree transitions handled in a later phase).
   !
   use sfincs_data
   !
   implicit none
   !
   integer   :: count0, count1, count_rate, count_max
   real      :: tloop
   real*4    :: dt
   !
   integer   :: ip, ipuv, nm, nmu, nmn, irow, iuv, iter
   real*4    :: dtrho, abf, hu, Dnm1, Dnmu, unh, gf, pL, pR, fdep, sl, slmax
   real*4    :: rz, rznew, dAd, alpha, beta, bnorm, rnorm
   real*4    :: dzdt, wmax, breform, qr, ql, kbfac, bf, ratio, pcap
   integer   :: iuvr, iuvl, rl, rr, ineighbour
   logical   :: bndok
   real*4, dimension(:), allocatable :: bb
   !
   call system_clock(count0, count_rate, count_max)
   !
   if (nrows == 0) return
   !
   ! Keller-box vertical factor: 2.0 = linear-pressure single layer (solver 2).
   ! Solver 3 makes it tunable (nh_disp); >2 strengthens the short-wave pressure
   ! feedback, adding dispersion to flatten c(k) toward Airy.
   !
   kbfac = 2.0
   if (nonh_solver == 3) kbfac = nh_disp
   !
   allocate(bb(nrows))
   !
   dtrho = dt / rhow
   !
   ! 1) Layer depth per row
   !
   do irow = 1, nrows
      nm = nm_index_of_row(irow)
      Dnm(irow) = max(zs(nm) - zb(nm), huthresh_nh)
   enddo
   !
   ! 1b) Hydrostatic-front (wave-breaking) reduction -- a GRADUAL form of the XBeach
   !     HFA (Smit, Zijlema & Stelling 2013). A breaking wave becomes a bore, well
   !     described hydrostatically; the non-hydrostatic pressure there is wrong and
   !     destabilizing. We detect steepening from the surface rate of rise, built
   !     from the flux divergence (= d(zs)/dt by continuity = Miche steepness):
   !         d(zs)/dx = (d(zs)/dt)/c = dzdt / sqrt(g h),  threshold = nh_brsteep.
   !     Instead of XBeach's hard on/off (which abruptly decouples the front and makes
   !     the wave fall apart), we ramp a per-cell factor nh_brfac smoothly from 1 (full
   !     non-hydrostatic, at/below the threshold) to 0 (fully hydrostatic, at a steepness
   !     ratio of 1+nh_brwidth). nh_brfac then scales BOTH the horizontal face coupling
   !     (step 2) and the vertical source (step 4), so the cell transitions smoothly to
   !     a hydrostatic bore -- SPD-consistent, no abrupt collapse. nh_brwidth = 0 gives
   !     the old hard binary switch. Only the rising front (dzdt > 0) is reduced.
   !
   if (nh_brsteep > 0.0) then
      breform = nh_reformsteep
      if (breform <= 0.0) breform = 0.25 * nh_brsteep
      ! --- (a) raw surface rate-of-rise  dzdt = -div(q)  (= d(zs)/dt by continuity) into cg_z ---
      do irow = 1, nrows
         nm   = nm_index_of_row(irow)
         iuvr = z_index_uv_mu(nm) ; iuvl = z_index_uv_md(nm)
         qr = 0.0 ; ql = 0.0
         if (iuvr > 0) qr = q(iuvr)
         if (iuvl > 0) ql = q(iuvl)
         dzdt = - (qr - ql) * dxrinv(1)
         if (nmax > 1) then
            iuvr = z_index_uv_nu(nm) ; iuvl = z_index_uv_nd(nm)
            qr = 0.0 ; ql = 0.0
            if (iuvr > 0) qr = q(iuvr)
            if (iuvl > 0) ql = q(iuvl)
            dzdt = dzdt - (qr - ql) * dyrinv(1)
         endif
         cg_z(irow) = dzdt
      enddo
      ! --- (b) SMOOTH dzdt before thresholding. Raw dzdt carries a 2dx component, so the
      !     per-cell criterion trips on noise -> scattered, isolated pnh=0 cells (no spatial
      !     coherence; the non-breaking gaps in between keep their source-driven pnh -> a
      !     2dx sawtooth). Two [1 2 1]/4 diffusion passes over the nh faces make the criterion
      !     fire over a CONNECTED front -> a contiguous breaking block. cg_z = raw -> smoothed. ---
      do iter = 1, 2
         do irow = 1, nrows
            cg_d(irow) = cg_z(irow)
         enddo
         do ip = 1, nhuv
            ipuv = nh_faceuv(ip)
            rl = row_index_of_nm(uv_index_z_nm(ipuv))
            rr = row_index_of_nm(uv_index_z_nmu(ipuv))
            if (rl > 0 .and. rr > 0) then
               cg_d(rl) = cg_d(rl) + 0.25 * (cg_z(rr) - cg_z(rl))
               cg_d(rr) = cg_d(rr) + 0.25 * (cg_z(rl) - cg_z(rr))
            endif
         enddo
         do irow = 1, nrows
            cg_z(irow) = cg_d(irow)
         enddo
      enddo
      ! --- (c) snapshot last step's breaking state (cg_r is free until the solve) ---
      do irow = 1, nrows
         cg_r(irow) = nh_brfac(irow)
      enddo
      ! --- (d) hysteretic state machine on the SMOOTHED dzdt ---
      do irow = 1, nrows
         nm   = nm_index_of_row(irow)
         dzdt = cg_z(irow)
         ! is a neighbour already breaking? (read the snapshot -> no within-sweep contamination)
         ineighbour = 0
         iuvl = z_index_uv_md(nm)
         if (iuvl > 0) then ; rl = row_index_of_nm(uv_index_z_nm(iuvl))
            if (rl > 0) then ; if (cg_r(rl) < 1.0) ineighbour = 1 ; endif ; endif
         iuvr = z_index_uv_mu(nm)
         if (iuvr > 0) then ; rr = row_index_of_nm(uv_index_z_nmu(iuvr))
            if (rr > 0) then ; if (cg_r(rr) < 1.0) ineighbour = 1 ; endif ; endif
         if (nmax > 1) then
            iuvl = z_index_uv_nd(nm)
            if (iuvl > 0) then ; rl = row_index_of_nm(uv_index_z_nm(iuvl))
               if (rl > 0) then ; if (cg_r(rl) < 1.0) ineighbour = 1 ; endif ; endif
            iuvr = z_index_uv_nu(nm)
            if (iuvr > 0) then ; rr = row_index_of_nm(uv_index_z_nmu(iuvr))
               if (rr > 0) then ; if (cg_r(rr) < 1.0) ineighbour = 1 ; endif ; endif
         endif
         wmax = sqrt(g * Dnm(irow))
         ! Hysteretic state machine (XBeach nonh_break): onset at nh_brsteep, spread to a
         ! breaking neighbour at the lower nh_reformsteep, and -- the key for limiting the
         ! leading wave -- a breaking cell STAYS breaking (pnh=0) until the surface falls
         ! (dzdt < 0), so it remains hydrostatic through the whole crest passage, not just
         ! the steep front -> a sustained, connected hydrostatic roller that dissipates.
         if (cg_r(irow) >= 1.0) then                                   ! was not breaking
            if (dzdt > nh_brsteep * wmax) then
               nh_brfac(irow) = 0.0
            elseif (dzdt > breform * wmax .and. ineighbour == 1) then
               nh_brfac(irow) = 0.0
            else
               nh_brfac(irow) = 1.0
            endif
         else                                                          ! was breaking
            if (dzdt < 0.0) then
               nh_brfac(irow) = 1.0                                    ! release only when surface falls
            else
               nh_brfac(irow) = 0.0
            endif
         endif
         if (nh_brfac(irow) < 1.0) pnh(irow) = 0.0          ! breaking cell -> pnh = 0 Dirichlet (held through the solve)
      enddo
      !
      ! Close interior holes so the breaking region is CONTIGUOUS. dzdt has a 2dx
      ! component in the steep front, so the per-cell criterion flips breaking on/off
      ! cell-to-cell -> alternating pnh=0 / pnh!=0 -> a 2dx sawtooth in pnh and the water
      ! level. A non-breaking cell flanked by breaking cells (both x-sides, or both y-sides
      ! in 2D) is filled in, repeated a few times to close wider gaps. (Conservative: only
      ! fills internal holes, never grows the region outward.)
      do iter = 1, 3
         do irow = 1, nrows
            cg_r(irow) = nh_brfac(irow)
         enddo
         do irow = 1, nrows
            if (cg_r(irow) >= 1.0) then
               nm = nm_index_of_row(irow)
               ineighbour = 0
               iuvl = z_index_uv_md(nm)
               if (iuvl > 0) then ; rl = row_index_of_nm(uv_index_z_nm(iuvl))
                  if (rl > 0) then ; if (cg_r(rl) < 1.0) ineighbour = ineighbour + 1 ; endif ; endif
               iuvr = z_index_uv_mu(nm)
               if (iuvr > 0) then ; rr = row_index_of_nm(uv_index_z_nmu(iuvr))
                  if (rr > 0) then ; if (cg_r(rr) < 1.0) ineighbour = ineighbour + 2 ; endif ; endif
               if (ineighbour == 3) then        ! breaking on both x-sides -> fill the hole
                  nh_brfac(irow) = 0.0 ; pnh(irow) = 0.0
               endif
               if (nmax > 1 .and. nh_brfac(irow) >= 1.0) then
                  ineighbour = 0
                  iuvl = z_index_uv_nd(nm)
                  if (iuvl > 0) then ; rl = row_index_of_nm(uv_index_z_nm(iuvl))
                     if (rl > 0) then ; if (cg_r(rl) < 1.0) ineighbour = ineighbour + 1 ; endif ; endif
                  iuvr = z_index_uv_nu(nm)
                  if (iuvr > 0) then ; rr = row_index_of_nm(uv_index_z_nmu(iuvr))
                     if (rr > 0) then ; if (cg_r(rr) < 1.0) ineighbour = ineighbour + 2 ; endif ; endif
                  if (ineighbour == 3) then ; nh_brfac(irow) = 0.0 ; pnh(irow) = 0.0 ; endif
               endif
            endif
         enddo
      enddo
      !
      ! Dilate the breaking block outward by nh_brdilate cells. The block edge otherwise
      ! ends at the crest, so the first NON-breaking cell on the back face carries a steep
      ! nh source pinned against a pnh=0 neighbour -> a ~2x-hydrostatic one-cell pnh spike
      ! that radiates 2dx noise down the trailing wave. Growing the pnh=0 region one cell
      ! puts its edge in the gentle (dzdt<0) back, where the adjacent source is small.
      do iter = 1, nh_brdilate
         do irow = 1, nrows
            cg_r(irow) = nh_brfac(irow)
         enddo
         do irow = 1, nrows
            if (cg_r(irow) >= 1.0) then        ! non-breaking: break if ANY neighbour breaks
               nm = nm_index_of_row(irow)
               iuvl = z_index_uv_md(nm)
               if (iuvl > 0) then ; rl = row_index_of_nm(uv_index_z_nm(iuvl))
                  if (rl > 0) then ; if (cg_r(rl) < 1.0) then ; nh_brfac(irow) = 0.0 ; pnh(irow) = 0.0 ; endif ; endif ; endif
               iuvr = z_index_uv_mu(nm)
               if (iuvr > 0) then ; rr = row_index_of_nm(uv_index_z_nmu(iuvr))
                  if (rr > 0) then ; if (cg_r(rr) < 1.0) then ; nh_brfac(irow) = 0.0 ; pnh(irow) = 0.0 ; endif ; endif ; endif
               if (nmax > 1) then
                  iuvl = z_index_uv_nd(nm)
                  if (iuvl > 0) then ; rl = row_index_of_nm(uv_index_z_nm(iuvl))
                     if (rl > 0) then ; if (cg_r(rl) < 1.0) then ; nh_brfac(irow) = 0.0 ; pnh(irow) = 0.0 ; endif ; endif ; endif
                  iuvr = z_index_uv_nu(nm)
                  if (iuvr > 0) then ; rr = row_index_of_nm(uv_index_z_nmu(iuvr))
                     if (rr > 0) then ; if (cg_r(rr) < 1.0) then ; nh_brfac(irow) = 0.0 ; pnh(irow) = 0.0 ; endif ; endif ; endif
               endif
            endif
         enddo
      enddo
   else
      nh_brfac = 1.0
   endif
   !
   ! 2) Per-face gradient coefficients. These depend on the water level (through
   !    the bed-slope / layer term AB) and on which faces are wet and active this
   !    step; inactive faces get zero coefficients (no pressure coupling).
   !
   do ip = 1, nhuv
      ipuv = nh_faceuv(ip)
      nm   = uv_index_z_nm(ipuv)
      nmu  = uv_index_z_nmu(ipuv)
      nh_cR(ip) = 0.0
      nh_cL(ip) = 0.0
      nh_hu(ip) = 0.0
      ! Open (water-level/outflow, kcs 2/3) boundary faces. Two regimes:
      !  - sponge OFF (nh_sponge_width = 0): SKIP the face (homogeneous Neumann).
      !    The boundary flux is handled hydrostatically; this is reflective for the
      !    non-hydrostatic part (the incident-wave collapse), the legacy behaviour.
      !  - sponge ON  (nh_sponge_width > 0): keep the face ACTIVE. The boundary cell
      !    is row 0, so it carries pnh = 0 (true Dirichlet) while the interior
      !    smoothing stays intact, and the ramped sponge diagonal (step 3) absorbs
      !    the outgoing non-hydrostatic energy -> a radiation/absorbing boundary.
      ! Closed walls (kcs 0) and Neumann boundaries (kcs 6) keep kfuv = 0 / are
      ! excluded here -> homogeneous Neumann, which is the correct solid-wall BC.
      if (nh_sponge_width > 0) then
         bndok = (kcs(nm)  == 1 .or. kcs(nm)  == 2 .or. kcs(nm)  == 3) .and. &
                 (kcs(nmu) == 1 .or. kcs(nmu) == 2 .or. kcs(nmu) == 3)
      else
         bndok = (kcs(nm) == 1 .and. kcs(nmu) == 1)
      endif
      if (kfuv(ipuv) == 1 .and. bndok) then
         if (zs(nm) > zb(nm) + huthresh_nh .and. zs(nmu) > zb(nmu) + huthresh_nh) then
            hu = max(zs(nm), zs(nmu)) - 0.5 * (zb(nm) + zb(nmu))
            if (hu > huthresh_nh) then
               Dnm1 = max(zs(nm)  - zb(nm),  huthresh_nh)
               Dnmu = max(zs(nmu) - zb(nmu), huthresh_nh)
               abf  = ( (zs(nmu) + zb(nmu)) - (zs(nm) + zb(nm)) ) / (Dnm1 + Dnmu)
               nh_cR(ip) = nh_cf(ip) * (1.0 + abf)
               nh_cL(ip) = nh_cf(ip) * (abf - 1.0)
               !
               ! Conveyance depth for the flux correction and the depth-weighted
               ! operator. SFINCS sets uv = q/max(hu,huvmin) in momentum, so q/uv
               ! IS that depth. Using it keeps the corrected q and uv consistent
               ! (same sign) and matches the depth continuity uses for the flux.
               ! Falls back to the max-surface depth for a still (uv=0) face.
               !
               if (uv(ipuv) /= 0.0) then
                  nh_hu(ip) = q(ipuv) / uv(ipuv)
               else
                  nh_hu(ip) = hu
               endif
               !
               ! Run-up depth taper: ramp the non-hydrostatic response to zero in
               ! thin water (full hydrostatic at/below nh_runuph0, full nonh
               ! at/above nh_runuph1). Caps the single-layer nonh overshoot in the
               ! steep wall run-up without touching the deeper shelf/shoaling.
               !
               if (nh_runuph1 > nh_runuph0) then
                  fdep = (hu - nh_runuph0) / (nh_runuph1 - nh_runuph0)
                  fdep = max(0.0, min(1.0, fdep))
                  nh_cR(ip) = nh_cR(ip) * fdep
                  nh_cL(ip) = nh_cL(ip) * fdep
               endif
               !
               ! Open-boundary fade-in: ramp the coupling 0 -> full over nh_fadein
               ! cells from the boundary, so the incident wave enters hydrostatically.
               !
               if (nh_fadein > 0) then
                  fdep = 1.0
                  if (nh_faceL(ip) > 0) fdep = min(fdep, nh_fade(nh_faceL(ip)))
                  if (nh_faceR(ip) > 0) fdep = min(fdep, nh_fade(nh_faceR(ip)))
                  nh_cR(ip) = nh_cR(ip) * fdep
                  nh_cL(ip) = nh_cL(ip) * fdep
               endif
               !
               ! Gradual breaking (HFA): the face coupling is intentionally LEFT ACTIVE
               ! here (unlike a hard decoupling). The non-hydrostatic effect is removed
               ! gradually via the vertical SOURCE only (step 4, scaled by nh_brfac), so a
               ! breaking cell keeps its pressure-gradient coupling to its neighbours
               ! (XBeach-style pnh~0 Dirichlet) and the steep front still gets the one-
               ! sided shoreward push from the soliton behind -> it hands over to a
               ! propagating hydrostatic bore instead of decoupling and stalling.
               !
            endif
         endif
      endif
   enddo
   !
   ! 2b) Wet/dry front: a non-hydrostatic cell with a dry neighbour is only
   !     coupled on its wet side, so its pnh is under-constrained and spikes
   !     (the wall run-up overshoot + shed 2dx). Switch such cells to hydrostatic
   !     by zeroing every face coefficient that touches them. (nh_dryfront = 1)
   !
   if (nh_dryfront == 1) then
      do irow = 1, nrows
         nh_isfront(irow) = 0
      enddo
      do ip = 1, nhuv
         ipuv = nh_faceuv(ip)
         nm   = uv_index_z_nm(ipuv)
         nmu  = uv_index_z_nmu(ipuv)
         if (nh_faceL(ip) > 0 .and. zs(nmu) <= zb(nmu) + huthresh_nh) nh_isfront(nh_faceL(ip)) = 1
         if (nh_faceR(ip) > 0 .and. zs(nm)  <= zb(nm)  + huthresh_nh) nh_isfront(nh_faceR(ip)) = 1
      enddo
      do ip = 1, nhuv
         if (nh_faceL(ip) > 0) then
            if (nh_isfront(nh_faceL(ip)) == 1) then
               nh_cR(ip) = 0.0 ; nh_cL(ip) = 0.0
            endif
         endif
         if (nh_faceR(ip) > 0) then
            if (nh_isfront(nh_faceR(ip)) == 1) then
               nh_cR(ip) = 0.0 ; nh_cL(ip) = 0.0
            endif
         endif
      enddo
   endif
   !
   ! 3) Vertical-acceleration diagonal, Laplacian diagonal, open-boundary sponge,
   !    and the resulting Jacobi diagonal of A. nh_diag accumulates the Laplacian
   !    diagonal here; the sponge then adds a ramped multiple of that Laplacian
   !    diagonal near open boundaries (stored in nh_spongediag so apply_A_nh2
   !    reproduces the same operator).
   !
   !    DEPTH-WEIGHTED (flux-divergence) form, consistent with the continuity
   !    equation and with XBeach: each face Laplacian term carries the face depth
   !    nh_hu, and the whole per-cell equation is in flux units (so the vertical
   !    diagonal is 2 dt/(rho H), not 2 dt/(rho H^2)). This makes the projection
   !    null the same flux divergence that continuity uses - essential where the
   !    depth varies (shoaling / run-up).
   !
   do irow = 1, nrows
      nh_dvert(irow)     = kbfac * dt / (rhow * Dnm(irow))
      nh_diag(irow)      = 0.0
      nh_spongediag(irow) = 0.0
   enddo
   do ip = 1, nhuv
      if (nh_faceR(ip) > 0) nh_diag(nh_faceR(ip)) = nh_diag(nh_faceR(ip)) + dtrho * nh_hu(ip) * nh_cR(ip)**2
      if (nh_faceL(ip) > 0) nh_diag(nh_faceL(ip)) = nh_diag(nh_faceL(ip)) + dtrho * nh_hu(ip) * nh_cL(ip)**2
   enddo
   if (nh_sponge_width > 0) then
      do irow = 1, nrows
         nh_spongediag(irow) = nh_sponge_ramp(irow) * nh_sponge_coef * nh_diag(irow)
      enddo
   endif
   do irow = 1, nrows
      nh_diag(irow) = nh_diag(irow) + nh_dvert(irow) + nh_spongediag(irow)
   enddo
   !
   ! 4) Right-hand side (flux-divergence form):  b = G^T (hu u*) - (ws + wb0 - 2 wb)
   !    The divergence uses the flux hu*u* (depth-weighted), and the vertical
   !    source is in flux units (no division by H), matching the operator above.
   !
   do irow = 1, nrows
      bb(irow) = - (ws(irow) + wb0(irow) - 2.0 * wb(irow))
   enddo
   ! (breaking is applied post-solve by scaling pnh -> 0 at breaking cells, keeping the
   !  faces active so the front keeps its one-sided pressure-gradient correction)
   do ip = 1, nhuv
      if (nh_cR(ip) == 0.0 .and. nh_cL(ip) == 0.0) cycle
      unh = nh_hu(ip) * uv(nh_faceuv(ip))   ! provisional flux  hu*u*  at this face
      if (nh_faceR(ip) > 0) bb(nh_faceR(ip)) = bb(nh_faceR(ip)) + nh_cR(ip) * unh
      if (nh_faceL(ip) > 0) bb(nh_faceL(ip)) = bb(nh_faceL(ip)) + nh_cL(ip) * unh
   enddo
   !
   ! 5) Pressure solve. Solvers 2/3: Jacobi-preconditioned CG (warm start).
   !     Solver 4: reduction-free weighted-Jacobi relaxation -- a fixed number of
   !     LOCAL sweeps (no global dot-products, no convergence check), warm-started
   !     from the previous step. This is the explicit artificial-compressibility /
   !     hyperbolized pressure update: embarrassingly parallel, GPU / quadtree-native
   !     (the global CG reductions are the parallel bottleneck this removes).
   !
   if (nonh_solver == 4) then
      do iter = 1, nh_subiter
         call apply_A_nh2(pnh, cg_Ap, dtrho)
         do irow = 1, nrows
            if (nh_brfac(irow) >= 1.0) &      ! breaking cells held at pnh = 0 (Dirichlet)
               pnh(irow) = pnh(irow) + nh_omega * (bb(irow) - cg_Ap(irow)) / nh_diag(irow)
         enddo
      enddo
   else
   !
   call apply_A_nh2(pnh, cg_Ap, dtrho)
   do irow = 1, nrows
      cg_r(irow) = bb(irow) - cg_Ap(irow)
      if (nh_brfac(irow) < 1.0) cg_r(irow) = 0.0   ! pnh=0 Dirichlet at breaking cells (faces stay active)
      cg_z(irow) = cg_r(irow) / nh_diag(irow)
      cg_d(irow) = cg_z(irow)
   enddo
   rz    = 0.0
   bnorm = 0.0
   do irow = 1, nrows
      rz    = rz    + cg_r(irow) * cg_z(irow)
      bnorm = bnorm + bb(irow) * bb(irow)
   enddo
   bnorm = sqrt(bnorm)
   if (bnorm <= 0.0) bnorm = 1.0
   !
   iter = 0
   do
      rnorm = 0.0
      do irow = 1, nrows
         rnorm = rnorm + cg_r(irow) * cg_r(irow)
      enddo
      rnorm = sqrt(rnorm)
      if (rnorm / bnorm < nh_tol) exit
      if (iter >= nh_itermax) exit
      !
      call apply_A_nh2(cg_d, cg_Ap, dtrho)
      dAd = 0.0
      do irow = 1, nrows
         dAd = dAd + cg_d(irow) * cg_Ap(irow)
      enddo
      if (dAd <= 0.0) exit            ! safeguard against round-off breakdown
      alpha = rz / dAd
      do irow = 1, nrows
         pnh(irow)  = pnh(irow)  + alpha * cg_d(irow)
         cg_r(irow) = cg_r(irow) - alpha * cg_Ap(irow)
         if (nh_brfac(irow) < 1.0) cg_r(irow) = 0.0   ! hold breaking cells at pnh=0
         cg_z(irow) = cg_r(irow) / nh_diag(irow)
      enddo
      rznew = 0.0
      do irow = 1, nrows
         rznew = rznew + cg_r(irow) * cg_z(irow)
      enddo
      beta = rznew / rz
      do irow = 1, nrows
         cg_d(irow) = cg_z(irow) + beta * cg_d(irow)
      enddo
      rz = rznew
      iter = iter + 1
   enddo
   endif   ! end pressure solve (CG for solvers 2/3, Jacobi relaxation for solver 4)
   !
   ! 5b) Spatial 2dx filter on pnh (nh_filter > 0). Damps the persistent grid-
   !     scale (2dx) mode without touching the resolved smooth pressure: the
   !     neighbour-mean of a 2dx pattern is -pnh (strongly damped), of a smooth
   !     field +pnh (passes through). One Jacobi smoothing pass over nonh
   !     neighbours; cg_d/cg_r/cg_z are free scratch after the CG solve.
   !
   if (nh_filter > 0.0) then
      do irow = 1, nrows
         cg_d(irow) = pnh(irow)
         cg_r(irow) = 0.0
         cg_z(irow) = 0.0
      enddo
      do ip = 1, nhuv
         if (nh_faceL(ip) > 0 .and. nh_faceR(ip) > 0) then
            cg_r(nh_faceL(ip)) = cg_r(nh_faceL(ip)) + cg_d(nh_faceR(ip))
            cg_z(nh_faceL(ip)) = cg_z(nh_faceL(ip)) + 1.0
            cg_r(nh_faceR(ip)) = cg_r(nh_faceR(ip)) + cg_d(nh_faceL(ip))
            cg_z(nh_faceR(ip)) = cg_z(nh_faceR(ip)) + 1.0
         endif
      enddo
      do irow = 1, nrows
         if (cg_z(irow) > 0.0) then
            pnh(irow) = (1.0 - nh_filter) * cg_d(irow) + nh_filter * cg_r(irow) / cg_z(irow)
         endif
      enddo
   endif
   !
   ! 5c) Localized 2dx smoothing of pnh in the marginal-nonh zones. The grid-scale
   !     (2dx) pressure mode the vertical source excites is normally suppressed by
   !     the horizontal Laplacian; it re-appears wherever that coupling is weak:
   !       - the open-boundary FADE-IN zone (nh_fade < 1), and
   !       - SHALLOW run-up water near a wall (D < nh_smoothdep), where the layer
   !         term dominates and the wave is near-hydrostatic anyway.
   !     One Jacobi neighbour-mean pass, weighted per cell by nh_smoothbnd times the
   !     larger of the two triggers, so it is identically zero in the resolved
   !     interior (nh_fade = 1 and D >= nh_smoothdep) and never touches the physical
   !     trailing waves in deeper water. A 2dx checkerboard has neighbour-mean =
   !     -pnh, so weight w damps it by (1-2w); a smooth field passes through.
   !     cg_d/cg_r/cg_z are free scratch after the CG solve.
   !
   if (nh_smoothbnd > 0.0 .and. (nh_fadein > 0 .or. nh_smoothdep > 0.0)) then
      do irow = 1, nrows
         cg_d(irow) = pnh(irow)
         cg_r(irow) = 0.0
         cg_z(irow) = 0.0
      enddo
      do ip = 1, nhuv
         if (nh_faceL(ip) > 0 .and. nh_faceR(ip) > 0) then
            cg_r(nh_faceL(ip)) = cg_r(nh_faceL(ip)) + cg_d(nh_faceR(ip))
            cg_z(nh_faceL(ip)) = cg_z(nh_faceL(ip)) + 1.0
            cg_r(nh_faceR(ip)) = cg_r(nh_faceR(ip)) + cg_d(nh_faceL(ip))
            cg_z(nh_faceR(ip)) = cg_z(nh_faceR(ip)) + 1.0
         endif
      enddo
      do irow = 1, nrows
         if (cg_z(irow) > 0.0) then
            fdep = 0.0
            if (nh_fadein > 0)        fdep = 1.0 - nh_fade(irow)               ! fade-zone trigger
            if (nh_smoothdep > 0.0)   fdep = max(fdep, (nh_smoothdep - Dnm(irow)) / nh_smoothdep)  ! shallow trigger
            gf = nh_smoothbnd * max(0.0, min(1.0, fdep))   ! local weight; 0 in resolved interior
            pnh(irow) = (1.0 - gf) * cg_d(irow) + gf * cg_r(irow) / cg_z(irow)
         endif
      enddo
   endif
   !
   ! 5d) (Breaking is now applied AS A DIRICHLET pnh=0 CONSTRAINT INSIDE the solve --
   !      breaking cells are held at pnh=0 with their faces left active, so the pressure
   !      field goes smoothly to zero at the front and the gradient is the consistent
   !      projection. This avoids the artificial pnh cliff a post-solve zeroing created.)
   !
   ! 5e) Depth limiter: cap the non-hydrostatic (bed) pressure at a fraction nh_pmax of
   !     the hydrostatic bed pressure rho*g*H. The total bed pressure rho*g*H + pnh then
   !     stays >= 0 for nh_pmax <= 1 (no "suction"). Resolved waves have |pnh| << rho*g*H
   !     so are untouched; this only clips the unphysical spikes at steep wall run-up,
   !     breaking fronts and grid-scale (2dx) modes -- a physically-scaled safety cap.
   !
   if (nh_pmax > 0.0) then
      do irow = 1, nrows
         pcap = nh_pmax * rhow * g * Dnm(irow)
         pnh(irow) = max(-pcap, min(pcap, pnh(irow)))
      enddo
   endif
   !
   ! 6) Correct fluxes / velocities with the pressure gradient. This is the
   !    momentum force -(dt/rho) G pnh integrated over the step; no nudging.
   !
   do ip = 1, nhuv
      if (nh_cR(ip) == 0.0 .and. nh_cL(ip) == 0.0) cycle
      ipuv = nh_faceuv(ip)
      pR = 0.0
      pL = 0.0
      if (nh_faceR(ip) > 0) pR = pnh(nh_faceR(ip))
      if (nh_faceL(ip) > 0) pL = pnh(nh_faceL(ip))
      gf  = nh_cR(ip) * pR + nh_cL(ip) * pL          ! (G pnh)_face
      unh = - dtrho * gf                              ! velocity increment
      
!      uv(ipuv) = uv(ipuv) + nh_relax * unh
!      q(ipuv)  = q(ipuv)  + nh_hu(ip) * nh_relax * unh
 
      uv(ipuv) = (1.0 - nh_fnudge) * uv(ipuv) + nh_fnudge * (uv(ipuv) + unh)
      q(ipuv)  = (1.0 - nh_fnudge) * q(ipuv) + nh_fnudge * (q(ipuv) + nh_hu(ip) * unh)

   
   
   enddo
   !
   ! 7) Update surface/bottom vertical velocities for the next step's forcing.
   !    The bottom kinematic condition is w_b = u . d(zb)/dx : the flow follows
   !    the (static) BED slope, built from the bed levels zb directly. The bed
   !    slope is optionally capped at +/- nh_dzbmax: a near-vertical wall step
   !    (|d(zb)/dx| ~ 20-80 here) would otherwise produce an enormous spurious
   !    w_b when a thin wall cell wets, triggering a 2dt explicit instability.
   !    Real bed slopes (<= 1:13 here) are well below the cap and untouched.
   !
   slmax = nh_dzbmax
   if (slmax <= 0.0) slmax = huge(1.0)   ! 0 = no cap
   !
   do irow = 1, nrows
      nm  = nm_index_of_row(irow)
      !
      ! Non-hydrostatically dry cell (h <= huthresh_nh): carry NO non-hydrostatic
      ! state. Reset pnh/ws/wb to zero so the cell is fully hydrostatic and starts
      ! clean when it next wets. (A just-wetted thin cell otherwise re-enters the
      ! solve with stale ws/wb/pnh, which is the residual wet/dry artifact at the
      ! repeatedly wetting/drying wall cell.)
      !
      if (zs(nm) - zb(nm) <= huthresh_nh) then
         wb0(irow) = 0.0
         wb(irow)  = 0.0
         ws(irow)  = 0.0
         pnh(irow) = 0.0
         nh_brfac(irow) = 1.0   ! clear breaking so a re-wetting cell starts clean
         cycle
      endif
      !
      wb0(irow) = wb(irow)
      wb(irow)  = 0.0
      !
      ! Only sum faces to NON-HYDROSTATICALLY-wet neighbours (h > huthresh_nh).
      ! Gating on kfuv alone (the SFINCS wet flag, h > huthresh ~ 1e-4) lets a
      ! repeatedly wetting/drying thin neighbour toggle a term in wb every step.
      !
      iuv = z_index_uv_md(nm)
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nm(iuv)
         if (zs(nmn) - zb(nmn) > huthresh_nh) then
            sl  = (zb(nmn) - zb(nm)) * dxrinv(1)
            sl  = max(-slmax, min(slmax, sl))
            wb(irow) = wb(irow) - 0.5 * uv(iuv) * sl
         endif
      endif
      iuv = z_index_uv_mu(nm)
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nmu(iuv)
         if (zs(nmn) - zb(nmn) > huthresh_nh) then
            sl  = (zb(nm) - zb(nmn)) * dxrinv(1)
            sl  = max(-slmax, min(slmax, sl))
            wb(irow) = wb(irow) - 0.5 * uv(iuv) * sl
         endif
      endif
      if (nmax > 1) then
         iuv = z_index_uv_nd(nm)
         if (kfuv(iuv) > 0) then
            nmn = uv_index_z_nm(iuv)
            if (zs(nmn) - zb(nmn) > huthresh_nh) then
               sl  = (zb(nmn) - zb(nm)) * dyrinv(1)
               sl  = max(-slmax, min(slmax, sl))
               wb(irow) = wb(irow) - 0.5 * uv(iuv) * sl
            endif
         endif
         iuv = z_index_uv_nu(nm)
         if (kfuv(iuv) > 0) then
            nmn = uv_index_z_nmu(iuv)
            if (zs(nmn) - zb(nmn) > huthresh_nh) then
               sl  = (zb(nm) - zb(nmn)) * dyrinv(1)
               sl  = max(-slmax, min(slmax, sl))
               wb(irow) = wb(irow) - 0.5 * uv(iuv) * sl
            endif
         endif
      endif
      !
      ws(irow) = ws(irow) - (wb(irow) - wb0(irow)) + (kbfac * dt / (rhow * Dnm(irow))) * pnh(irow)
   enddo
   !
   deallocate(bb)
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine


   subroutine apply_A_nh2(p, Ap, dtrho)
   !
   ! Matrix-free SPD operator (depth-weighted / flux-divergence form, XBeach-
   ! consistent):  Ap = (dt/rho) G^T diag(hu) G p  +  diag( 2 dt/(rho H) ) p .
   ! Each face Laplacian term carries the face depth nh_hu, so the operator nulls
   ! the same flux divergence the continuity equation uses. Still SPD (hu > 0,
   ! Gram form). Boundary cells (row index 0) carry p = 0 (Dirichlet), realised
   ! by not gathering/scattering them.
   !
   implicit none
   !
   real*4, intent(in)  :: p(:)
   real*4, intent(out) :: Ap(:)
   real*4, intent(in)  :: dtrho
   !
   integer :: ip, rL, rR, irow
   real*4  :: pL, pR, hf
   !
   do irow = 1, nrows
      Ap(irow) = (nh_dvert(irow) + nh_spongediag(irow)) * p(irow)
   enddo
   !
   do ip = 1, nhuv
      if (nh_cR(ip) == 0.0 .and. nh_cL(ip) == 0.0) cycle
      rL = nh_faceL(ip)
      rR = nh_faceR(ip)
      pL = 0.0
      pR = 0.0
      if (rR > 0) pR = p(rR)
      if (rL > 0) pL = p(rL)
      hf = dtrho * nh_hu(ip) * (nh_cR(ip) * pR + nh_cL(ip) * pL)
      if (rR > 0) Ap(rR) = Ap(rR) + nh_cR(ip) * hf
      if (rL > 0) Ap(rL) = Ap(rL) + nh_cL(ip) * hf
   enddo
   !
   end subroutine

end module
