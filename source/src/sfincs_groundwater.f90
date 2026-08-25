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
   !
   implicit none
   !
   private
   public :: initialize_groundwater, gw_face_transmissivity, gw_exchange_conductance
   public :: gw_cell_storage, gw_diffusion_number, gw_budget_update, gw_budget_report
   public :: gw_write_output
   !
   ! Cumulative water budget, m3. Recharge and exfiltration are integrated as they are applied.
   !
   real*8 :: gw_vol_recharge = 0.0d0
   real*8 :: gw_vol_exfiltration = 0.0d0
   real*8 :: gw_vol_initial = 0.0d0
   !
contains
   !
   subroutine initialize_groundwater()
   !
   implicit none
   !
   integer :: nm
   real*4  :: acell
   !
   allocate(gw_head(np))
   allocate(gw_head_n(np))
   allocate(gw_kh(np))
   allocate(gw_sy(np))
   allocate(gw_zbase(np))
   allocate(gw_recharge(np))
   allocate(gw_qexch(np))
   !
   gw_kh       = gw_kh_uniform
   gw_sy       = gw_sy_uniform
   gw_zbase    = gw_zbase_uniform
   gw_recharge = gw_recharge_uniform
   gw_qexch    = 0.0
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
      read(501) gw_head
      close(501)
      do nm = 1, np
         gw_head(nm) = max(gw_head(nm), gw_zbase(nm))
      enddo
   endif
   !
   gw_head_n = gw_head
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
   ! Record the initial stored volume so the budget can be closed later.
   !
   gw_vol_initial = 0.0d0
   do nm = 1, np
      if (kcs(nm) > 0) then
         call gw_cell_area(nm, acell)
         gw_vol_initial = gw_vol_initial + &
            dble(gw_sy(nm)) * dble(max(gw_head(nm) - gw_zbase(nm), 0.0)) * dble(acell)
      endif
   enddo
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
   subroutine gw_exchange_conductance(nm, cexch)
   !
   ! Surface/aquifer exchange, as a conductance in m3/s per m of head difference.
   !
   ! The SAME value is written to both off-diagonal blocks of the coupled matrix and added to
   ! both diagonals. That symmetry is what keeps the system SPD: a spike measured the coupled
   ! matrix at exactly 0.000e+00 symmetry residual with this construction, and anything one-sided
   ! leaves CG converging confidently to the wrong answer.
   !
   implicit none
   !
   integer, intent(in)  :: nm
   real*4,  intent(out) :: cexch
   !
   real*4 :: acell
   !
   call gw_cell_area(nm, acell)
   cexch = gw_leakance * acell
   !
   end subroutine gw_exchange_conductance
   !
   !
   subroutine gw_cell_storage(nm, head, vol, dvol)
   !
   ! Stored volume and its derivative at a given head.
   !
   ! Without subgrid the aquifer fills the whole cell, so storage is Sy * area per metre. With
   ! subgrid the water table only occupies the part of the cell below it, so the derivative is
   ! Sy times the WET AREA at that level, taken from the same table the surface uses.
   !
   ! The derivative is bounded to [gw_awet_floor, 1] x cell area. Bounding is safe because it
   ! only linearises: the residual form carries the volume exactly on the right-hand side, so a
   ! bound changes how fast the outer loop converges, not the answer it converges to. Leaving it
   ! unbounded put a millionfold spread on the diagonal in an earlier version and stalled CG
   ! completely.
   !
   implicit none
   !
   integer, intent(in)  :: nm
   real*4,  intent(in)  :: head
   real*4,  intent(out) :: vol, dvol
   !
   integer :: ilevel, ivol
   real*4  :: acell, b, awet, dzvol, dz, zmn, zmx
   !
   call gw_cell_area(nm, acell)
   b = max(head - gw_zbase(nm), 0.0)
   !
   if (.not. subgrid) then
      !
      vol  = gw_sy(nm) * b * acell
      dvol = gw_sy(nm) * acell
      !
   else
      !
      zmn = max(subgrid_z_zmin(nm), -20.0)
      zmx = max(subgrid_z_zmax(nm), -20.0)
      !
      if (head >= zmx) then
         awet = acell
      elseif (head <= zmn) then
         dzvol = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
         dz    = max(subgrid_z_dep(2, nm) - subgrid_z_dep(1, nm), 0.001)
         awet  = dzvol / dz
      else
         ivol = 1
         do ilevel = 2, subgrid_nlevels
            if (subgrid_z_dep(ilevel, nm) > head) then
               ivol = ilevel - 1
               exit
            endif
         enddo
         dzvol = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
         dz    = max(subgrid_z_dep(ivol + 1, nm) - subgrid_z_dep(ivol, nm), 0.001)
         awet  = dzvol / dz
      endif
      !
      awet = min(max(awet, gw_awet_floor * acell), acell)
      !
      vol  = gw_sy(nm) * b * awet
      dvol = gw_sy(nm) * awet
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
   subroutine gw_budget_update(dt)
   !
   implicit none
   real*4, intent(in) :: dt
   integer :: nm
   real*4  :: acell
   !
   do nm = 1, np
      if (kcs(nm) /= 1) cycle
      call gw_cell_area(nm, acell)
      gw_vol_recharge = gw_vol_recharge + dble(gw_recharge(nm)) * dble(dt) * dble(acell)
      gw_vol_exfiltration = gw_vol_exfiltration + dble(gw_qexch(nm)) * dble(dt) * dble(acell)
   enddo
   !
   end subroutine gw_budget_update
   !
   !
   subroutine gw_write_output(tnow)
   !
   ! Append the aquifer head to a flat stream file: one record of (time, head over all active
   ! cells in internal order) per map output time.
   !
   ! Deliberately not routed through sfincs_ncoutput. The aquifer is still being verified against
   ! analytical solutions, and this keeps that verification independent of the netCDF module.
   ! Promote it once the solver is settled.
   !
   implicit none
   !
   real*8, intent(in) :: tnow
   !
   logical, save :: opened = .false.
   !
   if (.not. opened) then
      open(unit = 502, file = 'gw_head.dat', form = 'unformatted', access = 'stream', &
           status = 'replace')
      opened = .true.
   endif
   !
   write(502) tnow
   write(502) gw_head
   flush(502)
   !
   end subroutine gw_write_output
   !
   !
   subroutine gw_budget_report()
   !
   ! Close the budget: recharge in must equal storage change plus exfiltration out.
   !
   implicit none
   integer :: nm
   real*4  :: acell
   real*8  :: vnow, resid, scale
   !
   vnow = 0.0d0
   do nm = 1, np
      if (kcs(nm) > 0) then
         call gw_cell_area(nm, acell)
         vnow = vnow + dble(gw_sy(nm)) * dble(max(gw_head(nm) - gw_zbase(nm), 0.0)) * dble(acell)
      endif
   enddo
   !
   resid = gw_vol_recharge - (vnow - gw_vol_initial) - gw_vol_exfiltration
   scale = max(abs(gw_vol_recharge), abs(vnow - gw_vol_initial), 1.0d0)
   !
   write(*,'(a)')        ' Groundwater budget (m3)'
   write(*,'(a,e14.6)')  '   recharge in       : ', gw_vol_recharge
   write(*,'(a,e14.6)')  '   storage change    : ', vnow - gw_vol_initial
   write(*,'(a,e14.6)')  '   exfiltration out  : ', gw_vol_exfiltration
   write(*,'(a,e14.6,a,f9.5,a)') '   residual          : ', resid, &
      '   (', 100.0d0 * resid / scale, ' %)'
   !
   end subroutine gw_budget_report
   !
end module sfincs_groundwater
