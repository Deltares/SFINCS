module sfincs_groundwater
   !
   ! Unconfined groundwater for SFINCS: a depth-averaged Boussinesq aquifer solved in the SAME
   ! semi-implicit system as the free surface, so surface and subsurface heads come out of one
   ! matrix rather than being exchanged between two solvers.
   !
   ! The aquifer equation is
   !
   !    Sy dh/dt = div( K (h - zb) grad h ) + R + q_exchange
   !
   ! Transmissivity K*(h - zb) depends on the solution, so the system is mildly nonlinear. That
   ! nonlinearity sits in the OFF-diagonals and is lagged (Picard), which keeps the matrix
   ! symmetric and lets CG stay. A Python spike measured Picard converging in under 3 outer
   ! iterations at the SFINCS operating point, with roughly 150x margin to first failure, because
   ! at a 30 s timestep the head moves ~1e-4 m and the lagged transmissivity is nearly exact.
   ! MODFLOW's Picard failures come from daily stress periods where the head moves metres.
   !
   ! The outer iteration is NOT built here. sfincs_semi_implicit already runs one for the subgrid
   ! storage term, with a stagnation exit; groundwater reassembles inside it.
   !
   use sfincs_data
   use sfincs_log
   !
   implicit none
   !
   private
   public :: initialize_groundwater, gw_face_transmissivity, gw_exchange_terms, gw_seepage_terms
   public :: gw_cell_storage, gw_diffusion_number, gw_budget_add, gw_budget_report
   public :: gw_subgrid_level, gw_explicit_step, gw_drain_terms
   !
   ! Cumulative volumes since the start of the run, m3. Signed so that a positive value is water
   ! ENTERING the aquifer. real*8 throughout: these are running totals over ~1e5 timesteps and
   ! real*4 would lose the small per-step increments entirely.
   !
   real*8 :: gw_vol_initial   = 0.0d0   ! storage at t = 0
   real*8 :: gw_vol_recharge  = 0.0d0   ! from gw_recharge, negative where it is drainage
   real*8 :: gw_vol_exchange  = 0.0d0   ! from the surface via leakance, negative = exfiltration
   real*8 :: gw_vol_boundary  = 0.0d0   ! lateral flux across faces to cells that are not unknowns
   real*8 :: gw_vol_ceiling   = 0.0d0   ! forced out by the topographic ceiling (seepage)
   real*8 :: gw_vol_drain     = 0.0d0   ! removed by the drain boundary, negative
   !
   ! Water that actually moved, m3: the sum of the MAGNITUDES of every elementary contribution --
   ! each face, each cell, each timestep. This is the scale the closure error has to be judged
   ! against, and it cannot be reconstructed from the signed totals above.
   !
   ! The signed totals hide two different cancellations, and both are the normal case rather than
   ! the exception. In TIME, a tidal aquifer takes water in on the flood and gives it back on the
   ! ebb, so the boundary total over a run is near zero while the water crossing the boundary is
   ! enormous -- Ferris read 12% closure that way against a residual of a few parts in ten
   ! thousand. In SPACE, a steady seepage problem has water entering one boundary and leaving the
   ! other in equal measure, so even the per-timestep net is near zero -- Dupuit read 0.86%
   ! against a true 0.08%. Only accumulating magnitudes at the point each contribution is computed
   ! sees through both.
   !
   real*8 :: gw_vol_gross     = 0.0d0
   integer, parameter :: gw_maxsub = 10000
   integer :: gw_nsub_max = 0
   !
contains
   !
   subroutine initialize_groundwater()
   !
   implicit none
   !
   integer :: nm
   real*4, dimension(:), allocatable :: rtmp4   ! real*4 buffer for the flat head file
   !
   allocate(gw_head(np))
   allocate(gw_head_n(np))
   allocate(gw_kh(np))
   allocate(gw_sy(np))
   allocate(gw_zbase(np))
   allocate(gw_recharge(np))
   allocate(gw_zdrain(np))
   !
   gw_kh       = gw_kh_uniform
   gw_sy       = gw_sy_uniform
   gw_zbase    = gw_zbase_uniform
   gw_recharge = gw_recharge_uniform
   gw_zdrain   = gw_zdrain_uniform
   !
   ! Optional spatial aquifer properties, each a flat binary over active cells in internal order
   ! -- the same convention manningfile and the head file use (sfincs_domain.f90:2005).
   !
   ! These are read BEFORE the initial head below, because the head is clamped to gw_zbase there
   ! and to the zbase FIELD if one was given. Reading them afterwards would clamp against the
   ! uniform value and then silently leave the head below the aquifer base wherever the field is
   ! higher, which shows up much later as a cell that stores nothing.
   !
   ! A distinct unit number per read: a stale handle here would read the wrong file into a
   ! parameter that is never printed, and nothing downstream would look wrong until the answer
   ! was.
   !
   if (gwkhfile(1:4) /= 'none') then
      write(*,'(a,a)') ' Groundwater: reading conductivity file ', trim(gwkhfile)
      open(unit = 504, file = trim(gwkhfile), form = 'unformatted', access = 'stream')
      read(504) gw_kh
      close(504)
   endif
   !
   if (gwsyfile(1:4) /= 'none') then
      write(*,'(a,a)') ' Groundwater: reading specific yield file ', trim(gwsyfile)
      open(unit = 505, file = trim(gwsyfile), form = 'unformatted', access = 'stream')
      read(505) gw_sy
      close(505)
   endif
   !
   if (gwzbasefile(1:4) /= 'none') then
      write(*,'(a,a)') ' Groundwater: reading aquifer base file ', trim(gwzbasefile)
      open(unit = 506, file = trim(gwzbasefile), form = 'unformatted', access = 'stream')
      read(506) gw_zbase
      close(506)
   endif
   !
   ! Initial head: gw_zsini if given, otherwise the initial surface level.
   !
   do nm = 1, np
      if (gw_zsini > -999.0) then
         gw_head(nm) = gw_zsini
      else
         gw_head(nm) = real(zs(nm))
      endif
      gw_head(nm) = max(gw_head(nm), gw_zbase(nm))
   enddo
   !
   ! Optional spatial initial head, flat binary over active cells in internal order, the same
   ! convention manningfile uses (sfincs_domain.f90:2005). Cells masked as boundary (kcs == 2)
   ! are not unknowns, so whatever they are given here they keep -- that is how a fixed-head
   ! aquifer boundary is imposed.
   !
   if (gwheadfile(1:4) /= 'none') then
      write(*,'(a,a)') ' Groundwater: reading head file ', trim(gwheadfile)
      open(unit = 501, file = trim(gwheadfile), form = 'unformatted', access = 'stream')
      !
      ! Through a real*4 buffer. The file is a flat real*4 stream and gw_head is real*8, so
      ! reading straight into it would consume two file records per cell and produce a plausible
      ! but wrong field -- the same class of failure as reading a netCDF initial condition as
      ! binary.
      !
      allocate(rtmp4(np))
      read(501) rtmp4
      gw_head = dble(rtmp4)
      deallocate(rtmp4)
      close(501)
      do nm = 1, np
         gw_head(nm) = max(gw_head(nm), gw_zbase(nm))
      enddo
   endif
   !
   !
   ! Optional spatial recharge, same flat-binary-over-active-cells convention as the head file.
   ! A NEGATIVE value is a sink: that is how a polder's ditch network is represented, as drainage
   ! removed from the aquifer over the drained area rather than from the surface.
   !
   if (gwrechargefile(1:4) /= 'none') then
      write(*,'(a,a)') ' Groundwater: reading recharge file ', trim(gwrechargefile)
      open(unit = 503, file = trim(gwrechargefile), form = 'unformatted', access = 'stream')
      read(503) gw_recharge
      close(503)
   endif
   !
   ! Optional spatial drain level, same convention.
   !
   if (gwzdrainfile(1:4) /= 'none') then
      write(*,'(a,a)') ' Groundwater: reading drain level file ', trim(gwzdrainfile)
      open(unit = 503, file = trim(gwzdrainfile), form = 'unformatted', access = 'stream')
      read(503) gw_zdrain
      close(503)
   endif
   !
   gw_drain_active = (gw_cdrain > 0.0 .and. (gw_zdrain_uniform > -998.0 .or. gwzdrainfile(1:4) /= 'none'))
   if (gw_drain_active) then
      write(*,'(a,e12.4,a)') ' Groundwater: drain boundary on, gw_cdrain = ', gw_cdrain, ' 1/s'
   endif
   !
   !
   ! The seepage face only exists on the semi-implicit path. The explicit path already handles the
   ! ceiling by moving the excess volume to gw_qsurf, and this flag is what keeps that path bit
   ! identical: it gates the storage-derivative floor in gw_cell_storage below.
   !
   gw_seepage_active = (semi_implicit .and. gw_seepage_fac > 0.0)
   !
   ! The ceiling the explicit path used last step. It has to be remembered, because the ceiling
   ! moves with the surface: when a pond drains the ground under it is exposed again and the
   ! ceiling FALLS. Storage above the new ceiling is then real water that has to seep out, and
   ! measuring the old head against the new ceiling drops it instead. See gw_explicit_step.
   !
   allocate(gw_zceil_n(np))
   do nm = 1, np
      if (subgrid) then
         gw_zceil_n(nm) = max(subgrid_z_zmax(nm), real(zs(nm)))
      else
         gw_zceil_n(nm) = max(zb(nm), real(zs(nm)))
      endif
   enddo
   !
   if (gw_seepage_active) then
      write(*,'(a,f8.3)') ' Groundwater: seepage face active, gw_seepage_fac = ', gw_seepage_fac
   else
      write(*,'(a)')      ' Groundwater: seepage face OFF - water is lost where the aquifer saturates'
   endif
   !
   gw_head_n = gw_head
   !
   ! Check what the spatial fields actually contain.
   !
   ! A zero or negative conductivity silently removes a cell from the lateral system without
   ! removing it from the matrix: the row keeps its storage term, so the solve still converges
   ! and the cell simply never exchanges water with its neighbours. That is far harder to
   ! diagnose from the answer than a stop is from the log. A negative specific yield is worse --
   ! it puts a negative number on the diagonal and breaks the SPD property CG depends on.
   !
   do nm = 1, np
      if (kcs(nm) == 0) cycle
      if (gw_kh(nm) <= 0.0) then
         write(*,*) 'Error: gw_kh <= 0 at cell ', nm, ' value ', gw_kh(nm)
         stop
      endif
      if (gw_sy(nm) < 0.0) then
         write(*,*) 'Error: gw_sy < 0 at cell ', nm, ' value ', gw_sy(nm)
         stop
      endif
   enddo
   !
   ! Guard the eigenvalue floor.
   !
   ! The smallest eigenvalue of the coupled system is set by the storage term plus the exchange
   ! conductance. A Python spike measured it: with Sy = 0 alone the exchange still catches it and
   ! the condition number roughly doubles, but with BOTH zero the floor collapses to 1.79e-05 and
   ! the condition number rises 700x to 6.4e6. Dry or confined cells are the realistic way in.
   !
   if (gw_leakance <= 0.0) then
      do nm = 1, np
         if (kcs(nm) > 0 .and. gw_sy(nm) <= 0.0) then
            write(*,*) 'Error: groundwater cell ', nm, ' has gw_sy = 0 and gw_leakance = 0.'
            write(*,*) '       The coupled matrix has no storage there and is badly conditioned.'
            stop
         endif
      enddo
   endif
   !
   ! Record the initial stored volume so the budget can be closed later. Through the same helper
   ! the report uses, so the two ends of the balance are measured the same way -- including the
   ! subgrid storage area and the topographic ceiling, which the old inline Sy*(h-zbase)*A ignored.
   !
   call gw_total_storage(gw_vol_initial)
   !
   write(*,'(a,i10,a,e12.4,a,e12.4)') ' Groundwater: cells ', np, &
      '  initial storage ', gw_vol_initial, ' m3   leakance ', gw_leakance
   !
   end subroutine initialize_groundwater
   !
   !
   subroutine gw_cell_area(nm, acell)
   !
   implicit none
   integer, intent(in)  :: nm
   real*4,  intent(out) :: acell
   !
   if (crsgeo) then
      acell = cell_area_m2(nm)
   else
      acell = cell_area(z_flags_iref(nm))
   endif
   !
   end subroutine gw_cell_area
   !
   !
   subroutine gw_face_transmissivity(ip, tface)
   !
   ! Face transmissivity, UPSTREAM weighted: saturated thickness taken from whichever side has
   ! the higher head, times that side's conductivity.
   !
   ! Not a harmonic mean. A spike found that with a harmonic mean a dry cell can never rewet --
   ! T = 0 on every one of its faces freezes it permanently, so the nonlinearity switches itself
   ! off and the cell is stuck for the rest of the run.
   !
   implicit none
   !
   integer, intent(in)  :: ip
   real*4,  intent(out) :: tface
   !
   integer :: nm, nmu
   real*4  :: b
   !
   nm  = uv_index_z_nm(ip)
   nmu = uv_index_z_nmu(ip)
   !
   if (gw_head(nm) >= gw_head(nmu)) then
      b     = gw_head(nm) - gw_zbase(nm)
      tface = gw_kh(nm) * max(b, 0.0)
   else
      b     = gw_head(nmu) - gw_zbase(nmu)
      tface = gw_kh(nmu) * max(b, 0.0)
   endif
   !
   end subroutine gw_face_transmissivity
   !
   !
   subroutine gw_exchange_terms(nm, zs_k, h_k, csym, qexpl)
   !
   ! Surface/aquifer exchange, split into a symmetric implicit part and a lagged remainder.
   !
   ! The exchange is NOT simply C*(zs - h). A cell whose surface is dry has no water to give: with
   ! the unconditional form, a dry cell sitting above its water table pumps water that does not
   ! exist into the aquifer and drives zs below the bed. Measured on a flat dry test the water
   ! table climbed from -5 m to -0.833 m, exactly the level it would reach by equilibrating with a
   ! phantom reservoir at the bed.
   !
   ! The physical form is the MODFLOW river/drain switch, written against the level zref at which
   ! the cell starts to hold water:
   !
   !    Q = C * ( max(zs, zref) - max(h, zref) )        positive = surface into aquifer
   !
   ! which gives, in the four cases that matter:
   !    wet surface, connected water table   Q = C (zs - h)     ordinary two-way exchange
   !    wet surface, deep water table        Q = C (zs - zref)  infiltration set by ponded depth
   !    dry surface, deep water table        Q = 0              nothing to infiltrate
   !    dry surface, water table above ground  Q = C (zref - h) < 0   seepage out of the ground
   !
   ! Those max() branches make the two cross-derivatives unequal whenever only one side is
   ! connected, and an asymmetric matrix would cost us CG. So only the part that IS symmetric goes
   ! into the matrix, and the rest is lagged onto the right-hand side:
   !
   !    Q = csym * (zs - h)  +  qexpl
   !
   ! with csym = C when BOTH sides are connected and zero otherwise. qexpl is whatever the exact
   ! Q is minus what the implicit part already accounts for, so at outer convergence the pair
   ! reproduces Q exactly. Because the same qexpl is subtracted from the surface and added to the
   ! aquifer, the split stays conservative whichever branch it is in.
   !
   implicit none
   !
   integer, intent(in)  :: nm
   real*8,  intent(in)  :: zs_k, h_k
   real*4,  intent(out) :: csym, qexpl
   !
   real*4 :: acell, awet, cexch
   real*8 :: vsurf, zref, qk
   !
   call gw_cell_area(nm, acell)
   !
   if (subgrid) then
      !
      ! Partial seepage area. Surface water and the aquifer are only in contact over the part of
      ! the cell that is actually wet, so a cell holding water in a single channel exchanges over
      ! that channel's footprint, not over the whole cell. A floor keeps the coupled matrix from
      ! losing its exchange term entirely as a cell dries, which is what the eigenvalue guard in
      ! initialize_groundwater relies on.
      !
      call gw_subgrid_level(nm, zs_k, vsurf, awet)
      cexch = gw_leakance * min(max(awet, gw_awet_floor * acell), acell)
      zref  = dble(subgrid_z_zmin(nm))
      !
   else
      !
      cexch = gw_leakance * acell
      zref  = dble(zb(nm))
      !
   endif
   !
   qk = dble(cexch) * (max(zs_k, zref) - max(h_k, zref))
   !
   if (zs_k > zref .and. h_k > zref) then
      csym = cexch
   else
      csym = 0.0
   endif
   !
   ! The two terms are nearly equal and their difference is the whole signal, so the
   ! subtraction happens in real*8 and only the small remainder is rounded back down.
   !
   qexpl = real(qk - dble(csym) * (zs_k - h_k))
   !
   end subroutine gw_exchange_terms
   !
   !
   subroutine gw_seepage_terms(nm, zs_k, h_k, dt, cseep, zceil)
   !
   ! Seepage face: the outflow a saturated cell needs so that water arriving at a full aquifer
   ! becomes surface water instead of disappearing.
   !
   !    Q = cseep * (h - zceil)     for h > zceil, positive OUT of the aquifer
   !
   ! zceil is the level above which there is no pore space left, and it is the same expression
   ! gw_explicit_step uses for its hcap -- the two paths have to agree on where the ceiling is or
   ! they cannot be compared. It is max(ground, surface level) rather than just the ground,
   ! because the ground beneath standing water is saturated and the table can stand as high as the
   ! free surface there.
   !
   ! cseep * dt = Sy * A means one timestep removes exactly the volume that would have been stored
   ! above the ceiling, which is the implicit statement of what the explicit path does when it
   ! moves the excess to gw_qsurf. That leaves no free parameter at the default.
   !
   implicit none
   !
   integer, intent(in)  :: nm
   real*8,  intent(in)  :: zs_k, h_k
   real*4,  intent(in)  :: dt
   real*4,  intent(out) :: cseep
   real*8,  intent(out) :: zceil
   !
   real*4 :: acell
   !
   call gw_cell_area(nm, acell)
   !
   if (subgrid) then
      zceil = max(dble(subgrid_z_zmax(nm)), zs_k)
   else
      zceil = max(dble(zb(nm)), zs_k)
   endif
   !
   if (h_k > zceil .and. gw_seepage_active) then
      cseep = gw_seepage_fac * gw_sy(nm) * acell / dt
   else
      cseep = 0.0
   endif
   !
   end subroutine gw_seepage_terms
   !
   !
   subroutine gw_drain_terms(nm, h_k, cdrn, zdrn)
   !
   ! Drain boundary: Q = cdrn * (h - zdrn) for h > zdrn, positive OUT of the aquifer, gone from
   ! the model. The switch is evaluated at the outer iterate, the conductance is implicit on
   ! the aquifer diagonal, and there is no surface partner, so symmetry is untouched.
   !
   implicit none
   !
   integer, intent(in)  :: nm
   real*8,  intent(in)  :: h_k
   real*4,  intent(out) :: cdrn
   real*8,  intent(out) :: zdrn
   !
   real*4 :: acell
   !
   call gw_cell_area(nm, acell)
   zdrn = dble(gw_zdrain(nm))
   !
   if (gw_drain_active .and. zdrn > dble(gw_zbase(nm)) .and. h_k > zdrn) then
      cdrn = gw_cdrain * acell
   else
      cdrn = 0.0
   endif
   !
   end subroutine gw_drain_terms
   !
   !
   subroutine gw_subgrid_level(nm, z, vsurf, awet)
   !
   ! Surface water volume and wet area of a cell at level z, read from the subgrid table.
   !
   ! The table stores the level as a function of UNIFORMLY BINNED volume, so within a bin the
   ! volume is linear in z and the wet area dV/dz is piecewise constant. Both are returned from
   ! the same bin, which keeps awet the exact derivative of vsurf -- the residual form below
   ! depends on that consistency.
   !
   implicit none
   !
   integer, intent(in)  :: nm
   real*8,  intent(in)  :: z
   real*8,  intent(out) :: vsurf
   real*4,  intent(out) :: awet
   !
   integer :: ilevel, ivol
   real*4  :: acell, dz
   real*8  :: dzvol, zmn, zmx
   !
   call gw_cell_area(nm, acell)
   !
   zmn   = subgrid_z_zmin(nm)
   zmx   = subgrid_z_zmax(nm)
   dzvol = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
   !
   if (z <= zmn) then
      !
      ! Below the lowest point in the cell: nothing wet, the aquifer has the whole footprint.
      !
      vsurf = 0.0d0
      awet  = 0.0
      !
   elseif (z >= zmx) then
      !
      ! Above the highest point: the cell is fully flooded and there is no unsaturated ground
      ! left to store water in. This is the topographic ceiling.
      !
      ! The table stops at zmax, so the volume has to be extended by hand above it, and the
      ! extension is the same one sfincs_continuity uses when it inverts a full cell:
      ! volmax + acell*(z - zmax). Returning a FLAT volmax here -- which is what this did -- is
      ! inconsistent with returning awet = acell, because awet is supposed to be d(vsurf)/dz. The
      ! inconsistency lands in gw_cell_storage, whose subgrid branch computes the aquifer volume
      ! as Sy*(acell*b - vsurf): with vsurf flat, that keeps GROWING at Sy*acell per metre above
      ! the ceiling while dvol correctly reports zero. So a submerged subgrid cell could store
      ! unbounded groundwater above its own ground, which is the mirror image of the ceiling
      ! defect on the non-subgrid side. With the extension the aquifer volume comes out constant
      ! above zmax, which is what "no unsaturated ground left" means.
      !
      vsurf = dble(subgrid_z_volmax(nm)) + dble(acell) * (z - zmx)
      awet  = acell
      !
   else
      !
      ivol = 1
      do ilevel = 2, subgrid_nlevels
         if (subgrid_z_dep(ilevel, nm) > z) then
            ivol = ilevel - 1
            exit
         endif
      enddo
      !
      dz    = max(subgrid_z_dep(ivol + 1, nm) - subgrid_z_dep(ivol, nm), 1.0e-6)
      awet  = real(dzvol / dz)
      vsurf = (ivol - 1) * dzvol + dble(awet) * (z - dble(subgrid_z_dep(ivol, nm)))
      !
   endif
   !
   awet = min(max(awet, 0.0), acell)
   !
   end subroutine gw_subgrid_level
   !
   !
   subroutine gw_cell_storage(nm, head, vol, dvol, zs_in)
   !
   ! Stored groundwater volume and its derivative at a given head.
   !
   ! Without subgrid the aquifer fills the whole cell footprint, so storage is Sy * area per
   ! metre of head.
   !
   ! With subgrid it does not. At level z the cell is partly flooded: over the wet area the water
   ! is SURFACE water, and only over the remaining dry area is there unsaturated ground able to
   ! store groundwater. So the storage area is acell - awet(z), not awet(z), and the stored
   ! volume is its integral from the aquifer base up to the head:
   !
   !    vol(h) = Sy * [ acell * (h - zbase) - ( vsurf(h) - vsurf(zbase) ) ]
   !
   ! which is exact given the table, because vsurf is exactly the integral of awet. Written this
   ! way the topographic ceiling falls out on its own: as the head approaches the highest point
   ! in the cell the storage area goes to zero, so any further recharge has nowhere to go and
   ! must leave as exfiltration rather than pushing the water table above the ground.
   !
   ! Only the DERIVATIVE is floored, never the volume. The residual form carries vol exactly on
   ! the right-hand side, so flooring dvol changes how fast the outer loop converges, not the
   ! answer it converges to -- while leaving it unfloored puts a zero on the diagonal exactly
   ! when a cell saturates.
   !
   implicit none
   !
   integer, intent(in)  :: nm
   real*8,  intent(in)  :: head
   real*8,  intent(out) :: vol
   real*4,  intent(out) :: dvol
   real*8,  intent(in), optional :: zs_in
   !
   ! zs_in overrides the surface level the ceiling is taken from. The semi-implicit assembly needs
   ! it: zs(nm) still holds the level at time n during assembly, so capping the NEW iterate's
   ! storage at it freezes the ceiling for the whole step. All the recharge arriving at a
   ! saturated cell is then ejected as seepage, and the storage the aquifer gains under the
   ! deepening pond -- Sy*A*dzs, real water -- appears at the next step with no flux having
   ! supplied it. On the ceiling case that is 77.3 m3 over the run and a 5.8 % closure error.
   ! gw_total_storage measures the new-time cap, so the row has to as well.
   !
   ! Absent, the level is zs(nm), which is what the explicit path wants: it sub-steps within one
   ! surface level, so time n is the only level it has.
   !
   real*4  :: acell, a0, a1, adry
   real*8  :: b, v0, v1, zcap, zsurf
   !
   call gw_cell_area(nm, acell)
   b = max(head - dble(gw_zbase(nm)), 0.0d0)
   !
   if (present(zs_in)) then
      zsurf = zs_in
   else
      zsurf = zs(nm)
   endif
   !
   if (.not. subgrid) then
      !
      ! Topographic ceiling. Without subgrid the hypsometry is a step at the bed, so the ceiling
      ! is a single level -- but it is NOT simply zb. Where the ground is exposed the water table
      ! cannot rise above it, because water above the ground is surface water and has to seep out
      ! instead. Where the cell is flooded there is no such limit: the ground beneath a pond is
      ! saturated, and the table can stand as high as the free surface.
      !
      ! So the cap is max(zb, zs). Capping at zb alone would be wrong for a submerged cell and
      ! would freeze the aquifer under standing water; leaving it uncapped altogether let the
      ! table climb 1.70 m above dry ground in the compound case, held back only by how fast
      ! leakance could drain it.
      !
      zcap = max(dble(zb(nm)), zsurf)
      b    = max(min(head, zcap) - dble(gw_zbase(nm)), 0.0d0)
      vol  = dble(gw_sy(nm)) * b * dble(acell)
      !
      ! Above the cap the stored volume genuinely stops changing, so the true derivative is zero.
      ! The floor below is not physics -- it exists so the matrix diagonal does not vanish. Once
      ! the seepage face is carrying cseep*dt on that diagonal the floor is unnecessary, and it is
      ! actively harmful: it is a store the budget cannot see, and it is what let the head reach
      ! +232 m on the ceiling case.
      !
      if (head < zcap) then
         dvol = gw_sy(nm) * acell
      elseif (gw_seepage_active) then
         dvol = 0.0
      else
         dvol = gw_sy(nm) * gw_awet_floor * acell
      endif
      !
   else
      !
      call gw_subgrid_level(nm, dble(gw_zbase(nm)), v0, a0)
      call gw_subgrid_level(nm, head,         v1, a1)
      !
      vol  = dble(gw_sy(nm)) * max(dble(acell) * b - max(v1 - v0, 0.0d0), 0.0d0)
      !
      if (gw_seepage_active) then
         adry = min(max(acell - a1, 0.0), acell)
      else
         adry = min(max(acell - a1, gw_awet_floor * acell), acell)
      endif
      dvol = gw_sy(nm) * adry
      !
   endif
   !
   end subroutine gw_cell_storage
   !
   !
   subroutine gw_diffusion_number(dt, numax_seen)
   !
   ! Largest diffusion number Nu = K b dt / (Sy dx^2) anywhere in the domain.
   !
   ! This is what limits the timestep once groundwater is on. A spike measured Picard converging
   ! in under 3 outer iterations at Nu ~ 0.03, first missing between Nu 3.9 and 5.6, and diverging
   ! outright by 55.6. Semi-implicit exists to escape the gravity-wave CFL, so as soon as dt grows
   ! the groundwater diffusion limit becomes the binding constraint -- the opposite of the usual
   ! assumption that groundwater is the slow, cheap part.
   !
   implicit none
   !
   real*4, intent(in)  :: dt
   real*4, intent(out) :: numax_seen
   !
   integer :: nm
   real*4  :: b, dxr, nu
   !
   numax_seen = 0.0
   !
   do nm = 1, np
      if (kcs(nm) /= 1) cycle
      if (gw_sy(nm) <= 0.0) cycle
      b   = max(gw_head(nm) - gw_zbase(nm), 0.0)
      dxr = 1.0 / dxrinv(z_flags_iref(nm))
      nu  = gw_kh(nm) * b * dt / (gw_sy(nm) * dxr * dxr)
      numax_seen = max(numax_seen, nu)
   enddo
   !
   end subroutine gw_diffusion_number
   !
   !
   subroutine gw_budget_add(v_recharge, v_exchange, v_boundary, v_ceiling, v_drain, v_gross)
   !
   ! Accumulate one timestep's worth of aquifer volume terms. Called by whichever solver path is
   ! active, with volumes it has already computed -- recomputing them here would risk the budget
   ! measuring something subtly different from what the solver did, which is the one thing a
   ! budget must not do.
   !
   implicit none
   !
   real*8, intent(in) :: v_recharge, v_exchange, v_boundary, v_ceiling, v_drain
   real*8, intent(in) :: v_gross
   !
   gw_vol_recharge = gw_vol_recharge + v_recharge
   gw_vol_exchange = gw_vol_exchange + v_exchange
   gw_vol_boundary = gw_vol_boundary + v_boundary
   gw_vol_ceiling  = gw_vol_ceiling  + v_ceiling
   gw_vol_drain    = gw_vol_drain    + v_drain
   !
   gw_vol_gross = gw_vol_gross + v_gross
   !
   end subroutine gw_budget_add
   !
   !
   subroutine gw_total_storage(vtot)
   !
   ! Total stored groundwater volume over the control volume, m3.
   !
   ! The control volume is the interior cells only. A cell with kcs == 2 carries a prescribed head
   ! and is not an unknown, so it sits OUTSIDE the balance and the water crossing into it is
   ! counted as lateral boundary flux instead.
   !
   ! Through gw_cell_storage rather than Sy*(h - zbase)*A, so that the subgrid storage area and
   ! the topographic ceiling are included -- the budget has to measure the same volume the solver
   ! conserves, not an idealisation of it.
   !
   implicit none
   !
   real*8, intent(out) :: vtot
   !
   integer :: nm
   real*8  :: vol
   real*4  :: dvol
   !
   vtot = 0.0d0
   !
   do nm = 1, np
      if (kcs(nm) /= 1) cycle
      call gw_cell_storage(nm, gw_head(nm), vol, dvol)
      vtot = vtot + vol
   enddo
   !
   end subroutine gw_total_storage
   !
   !
   subroutine gw_explicit_step(dt)
   !
   ! Advance the aquifer explicitly, for use with the explicit surface solver.
   !
   ! The aquifer is a diffusion problem, so an explicit update is limited by
   !
   !     dt <= Sy dx^2 / (4 K b)
   !
   ! which sounds restrictive but is not, because the explicit surface solver is already crawling
   ! at the gravity-wave CFL and water diffuses through an aquifer far more slowly than a shallow
   ! water wave propagates. Measured against dt = alfa dx / sqrt(g h) the aquifer limit is 93x
   ! looser on the tidal island case, 3000x on the polder, and 6800x on a sandy coast. Only very
   ! high conductivity on a fine grid brings the two together, and that is what the sub-cycling
   ! below is for: the aquifer block is small, so taking several sub-steps of it per surface step
   ! is cheap.
   !
   ! Geometry is deliberately identical to what the semi-implicit assembly uses -- same face
   ! widths, same centre-to-centre distances, same 1.5x factor across a refinement transition. The
   ! whole point of having two paths is that they can be compared, and that only works if the
   ! discretisation is the same and only the time integration differs.
   !
   ! Conservative by construction: each face is visited once and its flux is added to one cell and
   ! subtracted from the other, and any water that will not fit under the topographic ceiling is
   ! handed to the surface as seepage rather than discarded.
   !
   use sfincs_data
   !
   implicit none
   !
   real*4, intent(in) :: dt
   !
   integer :: ip, nm, nmu, nsub, it
   real*4  :: tface, wface, dinv, qface, dtsub, nurate, tsub
   real*4  :: acell, dvol
   real*8  :: vol, volcap, hcap, excess, hnew
   real*4  :: csym, qexpl, qex
   real*8  :: volh
   real*4  :: dvolh
   real*8  :: bv_rech, bv_exch, bv_bnd, bv_ceil, bv_drain, bv_gross
   real*4  :: cdrn
   real*8  :: zdrn, qdrn
   !
   if (.not. allocated(gw_dvol)) allocate(gw_dvol(np))
   if (.not. allocated(gw_qsurf)) allocate(gw_qsurf(np))
   !
   ! How many sub-steps does stability demand?
   !
   ! Where the aquifer meets a prescribed water level, its head is that water level. The
   ! semi-implicit path does this when it seeds the outer iterate; the explicit path has no outer
   ! iterate, so it has to be done here. Without it a tidal boundary drives nothing: the Ferris
   ! case came back with a completely flat aquifer, no decay length at all, because the boundary
   ! head never moved.
   !
   if (gw_bnd_from_zs) then
      do nm = 1, np
         if (kcs(nm) == 2) gw_head(nm) = zs(nm)
      enddo
   endif
   !
   ! Sub-step to whatever explicit stability demands, recomputing the limit each pass.
   !
   ! Forward-in-time centred-in-space diffusion is stable only for
   !
   !     Nu = K b dt / (Sy dx^2) <= 1/4   in two dimensions
   !
   ! Recomputed every sub-step rather than once, because the saturated thickness b moves with the
   ! head and a step that was stable at the start of the surface timestep need not stay so.
   ! Wflow.jl does the same in its groundwater module, calling the coefficient alpha and citing
   ! Chu & Willis (1984); its default is 0.25, which is where gw_numax's default now sits too.
   !
   ! It did NOT sit there originally. gw_numax defaulted to 4.0, chosen for the Picard convergence
   ! limit of the semi-implicit path -- which is unconditionally stable and never calls this
   ! routine. Sub-stepping an EXPLICIT scheme down to Nu = 4 is 16x past the limit. No test caught
   ! it because every case already ran at Nu below 0.5, so the sub-stepping never engaged; raising
   ! Dupuit's conductivity to 1e-1 m/s puts Nu at 5.0 and the head reaches 1.6e6 m.
   !
   gw_qsurf = 0.0
   tsub = 0.0
   nsub = 0
   !
   ! Budget terms for this surface timestep, signed as water entering the aquifer. Accumulated
   ! over the sub-steps and handed over once, so the budget sees the same volumes the scheme
   ! actually moved rather than a re-derivation of them.
   !
   ! Not accumulated: the vol = max(vol, 0.0) floor further down, and the three-step Newton
   ! inversion of the storage relation. Those are numerics, not fluxes, and leaving them out is
   ! deliberate -- they are exactly what the reported closure error is there to expose.
   !
   bv_rech  = 0.0d0
   bv_exch  = 0.0d0
   bv_bnd   = 0.0d0
   bv_ceil  = 0.0d0
   bv_drain = 0.0d0
   bv_gross = 0.0d0
   !
   do while (tsub < dt)
      !
      call gw_diffusion_number(1.0, nurate)     ! diffusion number per second
      if (nurate > 0.0) then
         dtsub = gw_numax / nurate
      else
         dtsub = dt
      endif
      dtsub = min(dtsub, dt - tsub)
      !
      nsub = nsub + 1
      if (nsub > gw_maxsub) then
         write(*,*) 'Error: groundwater explicit sub-stepping exceeded ', gw_maxsub, ' steps.'
         write(*,*) '       Diffusion number per second is ', nurate
         write(*,*) '       The aquifer is much stiffer than the surface here; use semi_implicit = 1.'
         stop
      endif
      !
      gw_dvol = 0.0d0
      !
      ! Lateral flux, one pass over the faces.
      !
      do ip = 1, npuv
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         if (nm == 0 .or. nmu == 0) cycle
         if (kcs(nm) == 0 .or. kcs(nmu) == 0) cycle
         !
         call gw_face_transmissivity(ip, tface)
         if (tface <= 0.0) cycle
         !
         if (uv_flags_dir(ip) == 0) then
            wface = dyrm(uv_flags_iref(ip))
            dinv  = dxrinv(uv_flags_iref(ip))
         else
            wface = dxrm(uv_flags_iref(ip))
            dinv  = dyrinv(uv_flags_iref(ip))
         endif
         if (uv_flags_type(ip) /= 0) dinv = dinv / 1.5
         !
         qface = tface * wface * dinv * (gw_head(nm) - gw_head(nmu))
         !
         gw_dvol(nm)  = gw_dvol(nm)  - dble(qface) * dble(dtsub)
         gw_dvol(nmu) = gw_dvol(nmu) + dble(qface) * dble(dtsub)
         !
         ! Lateral boundary flux. A cell with kcs == 2 holds a prescribed head and is not part of
         ! the control volume, so a face touching one carries water across the boundary. qface is
         ! signed from nm towards nmu, so it ENTERS the control volume when nm is the outside
         ! cell and LEAVES when nmu is. A face with kcs == 2 on both sides nets to zero, which is
         ! right: neither cell is inside.
         !
         if (kcs(nm)  == 2) then
            bv_bnd   = bv_bnd + dble(qface) * dble(dtsub)
            bv_gross = bv_gross + abs(dble(qface) * dble(dtsub))
         endif
         if (kcs(nmu) == 2) then
            bv_bnd   = bv_bnd - dble(qface) * dble(dtsub)
            bv_gross = bv_gross + abs(dble(qface) * dble(dtsub))
         endif
         !
      enddo
      !
      ! Recharge, exchange with the surface, and the new head.
      !
      do nm = 1, np
         !
         if (kcs(nm) /= 1) cycle
         !
         call gw_cell_area(nm, acell)
         gw_dvol(nm) = gw_dvol(nm) + dble(acell) * dble(gw_recharge(nm)) * dble(dtsub)
         bv_rech  = bv_rech + dble(acell) * dble(gw_recharge(nm)) * dble(dtsub)
         bv_gross = bv_gross + abs(dble(acell) * dble(gw_recharge(nm)) * dble(dtsub))
         !
         ! Drain boundary, explicit in the sub-step's head. Leaves the model.
         !
         call gw_drain_terms(nm, gw_head(nm), cdrn, zdrn)
         qdrn = dble(cdrn) * (gw_head(nm) - zdrn) * dble(dtsub)
         gw_dvol(nm) = gw_dvol(nm) - qdrn
         bv_drain = bv_drain - qdrn
         bv_gross = bv_gross + abs(qdrn)
         !
         call gw_exchange_terms(nm, dble(zs(nm)), gw_head(nm), csym, qexpl)
         qex = real(dble(csym) * (zs(nm) - gw_head(nm))) + qexpl   ! positive: surface into aquifer
         gw_dvol(nm) = gw_dvol(nm) + dble(qex) * dble(dtsub)
         gw_qsurf(nm) = gw_qsurf(nm) - qex * dtsub             ! and the surface loses it
         bv_exch  = bv_exch + dble(qex) * dble(dtsub)
         bv_gross = bv_gross + abs(dble(qex) * dble(dtsub))
         !
         ! Convert the volume change into a head, honouring the topographic ceiling. Anything
         ! that will not fit below the ceiling has nowhere to go underground and becomes surface
         ! water, which is what a seepage face is.
         !
         ! Against the ceiling this cell had at the END of the last step, not the one it has
         ! now. The two differ whenever the surface moved, and taking the new one here loses the
         ! difference silently -- it never reaches the excess test below and never becomes
         ! seepage. With the ceiling FALLING, which is what a draining pond does, that is water
         ! destroyed: -5.6 % on the seepslope case before this line was made explicit.
         !
         call gw_cell_storage(nm, gw_head(nm), vol, dvol, gw_zceil_n(nm))
         vol = vol + gw_dvol(nm)
         !
         if (subgrid) then
            hcap = max(dble(subgrid_z_zmax(nm)), zs(nm))
         else
            hcap = max(dble(zb(nm)), zs(nm))
         endif
         call gw_cell_storage(nm, hcap, volcap, dvol)
         !
         if (vol > volcap) then
            excess = vol - volcap
            vol    = volcap
            gw_qsurf(nm) = gw_qsurf(nm) + real(excess)
            bv_ceil  = bv_ceil - excess               ! leaves the aquifer, so negative
            bv_gross = bv_gross + abs(excess)
         endif
         !
         vol = max(vol, 0.0d0)
         !
         ! Invert the storage relation. It is piecewise linear, so a few Newton steps are exact
         ! to round-off; the guard on the derivative only matters at a saturated cell.
         !
         hnew = gw_head(nm)
         do it = 1, 3
            call gw_cell_storage(nm, hnew, volh, dvolh)
            hnew = hnew + (vol - volh) / dble(max(dvolh, 1.0e-12))
            hnew = min(max(hnew, dble(gw_zbase(nm))), hcap)
         enddo
         gw_head(nm) = hnew
         gw_zceil_n(nm) = hcap
         !
      enddo
      !
      tsub = tsub + dtsub
      !
   enddo
   !
   call gw_budget_add(bv_rech, bv_exch, bv_bnd, bv_ceil, bv_drain, bv_gross)
   !
   gw_nsub_max = max(gw_nsub_max, nsub)
   !
   ! Hand the surface its share. Volume, not level, so that the subgrid path stays consistent
   ! with how continuity converts one to the other.
   !
   do nm = 1, np
      if (kcs(nm) /= 1 .or. gw_qsurf(nm) == 0.0) cycle
      call gw_cell_area(nm, acell)
      if (subgrid) then
         z_volume(nm) = max(z_volume(nm) + dble(gw_qsurf(nm)), 0.0d0)
         call gw_level_from_volume(nm, acell)
      else
         zs(nm) = max(zs(nm) + dble(gw_qsurf(nm) / acell), dble(zb(nm)))
      endif
   enddo
   !
   end subroutine gw_explicit_step
   !
   !
   subroutine gw_level_from_volume(nm, acell)
   !
   ! Water level from cell volume, mirroring sfincs_continuity.f90:562-583. Duplicated rather
   ! than shared because continuity does it inline inside its own loop; if that logic changes,
   ! this has to follow.
   !
   implicit none
   !
   integer, intent(in) :: nm
   real*4,  intent(in) :: acell
   !
   integer :: iuv
   real*4  :: dzvol, facint
   !
   if (z_volume(nm) >= subgrid_z_volmax(nm) * 0.999) then
      zs(nm) = max(subgrid_z_zmax(nm), -20.0) + (z_volume(nm) - subgrid_z_volmax(nm)) / acell
   elseif (z_volume(nm) <= 1.0e-6) then
      zs(nm) = max(subgrid_z_zmin(nm), -20.0)
   else
      dzvol  = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
      iuv    = int(z_volume(nm) / dzvol) + 1
      facint = (z_volume(nm) - (iuv - 1) * dzvol) / dzvol
      zs(nm) = subgrid_z_dep(iuv, nm) &
             + (subgrid_z_dep(iuv + 1, nm) - subgrid_z_dep(iuv, nm)) * facint
   endif
   !
   end subroutine gw_level_from_volume
   !
   !
   subroutine gw_budget_report()
   !
   ! Closure is the only number here that matters. Everything else is context for it.
   !
   ! Every term is signed as water ENTERING the aquifer, so the storage change must equal their
   ! sum. The residual is judged against gw_vol_gross, the water that actually moved, and not
   ! against the net inflow -- see the comment on gw_vol_gross for why the two differ so much on
   ! anything tidal.
   !
   ! Through write_log rather than write(*,*), so the balance lands in sfincs.log next to the
   ! rest of the run summary. On the cluster stdout is not always kept, and a diagnostic that only
   ! exists in a terminal that has since closed is not a diagnostic.
   !
   implicit none
   !
   real*8 :: vnow, dstore, vin, resid
   !
   call gw_total_storage(vnow)
   !
   dstore = vnow - gw_vol_initial
   vin    = gw_vol_recharge + gw_vol_exchange + gw_vol_boundary + gw_vol_ceiling + gw_vol_drain
   resid  = dstore - vin
   !
   call write_log('', 1)
   call write_log(' ---------- Groundwater water balance ----------', 1)
   write(logstr,'(a,e14.6,a)') ' Initial storage      : ', gw_vol_initial,  ' m3'
   call write_log(logstr, 1)
   write(logstr,'(a,e14.6,a)') ' Final storage        : ', vnow,            ' m3'
   call write_log(logstr, 1)
   write(logstr,'(a,e14.6,a)') ' Storage change       : ', dstore,          ' m3'
   call write_log(logstr, 1)
   write(logstr,'(a,e14.6,a)') '   recharge           : ', gw_vol_recharge, ' m3'
   call write_log(logstr, 1)
   write(logstr,'(a,e14.6,a)') '   exchange w/ surface: ', gw_vol_exchange, ' m3'
   call write_log(logstr, 1)
   write(logstr,'(a,e14.6,a)') '   lateral boundary   : ', gw_vol_boundary, ' m3'
   call write_log(logstr, 1)
   write(logstr,'(a,e14.6,a)') '   ceiling seepage    : ', gw_vol_ceiling,  ' m3'
   call write_log(logstr, 1)
   write(logstr,'(a,e14.6,a)') '   drainage           : ', gw_vol_drain,    ' m3'
   call write_log(logstr, 1)
   write(logstr,'(a,e14.6,a)') ' Throughput           : ', gw_vol_gross,    ' m3'
   call write_log(logstr, 1)
   write(logstr,'(a,e14.6,a)') ' Closure error        : ', resid,           ' m3'
   call write_log(logstr, 1)
   if (gw_vol_gross > 0.0d0) then
      write(logstr,'(a,f12.6,a)') ' Closure error        : ', &
         100.0d0 * resid / gw_vol_gross, ' % of throughput'
      call write_log(logstr, 1)
   endif
   call write_log(' -----------------------------------------------', 1)
   call write_log('', 1)
   !
   end subroutine gw_budget_report
   !
end module sfincs_groundwater
