module sfincs_src_structures
   !
   ! Point structures that move water between two grid cells by user-specified
   ! rules rather than by momentum conservation:
   !    type 1 - pump          (fixed discharge)
   !    type 2 - culvert       (bidirectional, weir-like)
   !    type 3 - check valve   (unidirectional culvert)
   !    type 4 - controlled gate, water-level triggered
   !    type 5 - controlled gate, schedule triggered
   !
   ! These used to live in sfincs_discharges.f90 alongside the river point
   ! discharges read from src/dis/netsrcdis. They have been split out so that
   ! each module has a single responsibility.
   !
   ! Runtime handoff to the continuity module is via the cell-wise qsrc(np)
   ! array (in sfincs_data): this module accumulates qq on intake (nmindrn_in)
   ! and outfall (nmindrn_out) cells. Per-structure signed discharge is also
   ! stored in qdrain(ndrn) for his output.
   !
   ! Concurrency: qsrc updates use atomic because two structures (or a river
   ! source and a structure) can land in the same cell.
   !
   use sfincs_log
   use sfincs_error

contains
   !
   subroutine initialize_src_structures()
   !
   ! Parse drnfile and populate drainage_type/_params/_status/_distance/
   ! _fraction_open, nmindrn_in(ndrn), nmindrn_out(ndrn), and the output
   ! buffer qdrain(ndrn).
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   real*4, dimension(:), allocatable :: xsrc_drn, ysrc_drn
   real*4, dimension(:), allocatable :: xsnk,     ysnk
   real*4    :: dummy, xsnk_tmp, ysnk_tmp, xsrc_tmp, ysrc_tmp
   integer   :: idrn, nmq, stat, npars
   logical   :: ok
   character(len=256) :: drainage_line
   !
   ndrn = 0
   !
   if (drnfile(1:4) == 'none') return
   !
   ok = check_file_exists(drnfile, 'Drainage points drn file', .true.)
   !
   ! Count lines
   !
   open(501, file=trim(drnfile))
   do while (.true.)
      read(501, *, iostat=stat) dummy
      if (stat < 0) exit
      ndrn = ndrn + 1
   enddo
   rewind(501)
   !
   if (ndrn <= 0) then
      close(501)
      return
   endif
   !
   write(logstr,'(a,a,a,i0,a)')' Reading ', trim(drnfile), ' (', ndrn, ' drainage points found) ...'
   call write_log(logstr, 0)
   !
   allocate(xsrc_drn(ndrn))
   allocate(ysrc_drn(ndrn))
   allocate(xsnk(ndrn))
   allocate(ysnk(ndrn))
   !
   allocate(nmindrn_in(ndrn))
   allocate(nmindrn_out(ndrn))
   allocate(qdrain(ndrn))
   allocate(drainage_type(ndrn))
   allocate(drainage_params(ndrn, 6))
   allocate(drainage_status(ndrn))
   allocate(drainage_distance(ndrn))
   allocate(drainage_fraction_open(ndrn))
   !
   nmindrn_in             = 0
   nmindrn_out            = 0
   qdrain                 = 0.0
   drainage_params        = 0.0
   drainage_distance      = 0.0
   drainage_fraction_open = 1.0   ! initially fully open (could be refined from zmin/zmax)
   drainage_status        = 1     ! 0=closed, 1=open, 2=closing, 3=opening
   !
   do idrn = 1, ndrn
      !
      read(501, '(a)') drainage_line
      !
      ! Determine drainage type first (5th integer in the line)
      !
      read(drainage_line, *, iostat=stat) xsnk_tmp, ysnk_tmp, xsrc_tmp, ysrc_tmp, drainage_type(idrn)
      !
      npars = 0
      !
      if (drainage_type(idrn) == 1 .or. drainage_type(idrn) == 2 .or. drainage_type(idrn) == 3) then
         npars = 1      ! pump, culvert, or check valve
      elseif (drainage_type(idrn) == 4 .or. drainage_type(idrn) == 5) then
         npars = 6      ! controlled gate (width, sill, manning, zmin/tclose, zmax/topen, closing time)
      endif
      !
      if (npars == 0) then
         write(logstr,'(a,i0,a)') 'Drainage type ', drainage_type(idrn), ' not recognized !'
         call stop_sfincs(logstr, -1)
      endif
      !
      if (npars == 1) then
         read(drainage_line, *, iostat=stat) xsnk(idrn), ysnk(idrn), xsrc_drn(idrn), ysrc_drn(idrn), &
              drainage_type(idrn), drainage_params(idrn, 1)
      elseif (npars == 6) then
         read(drainage_line, *, iostat=stat) xsnk(idrn), ysnk(idrn), xsrc_drn(idrn), ysrc_drn(idrn), &
              drainage_type(idrn), drainage_params(idrn, 1), drainage_params(idrn, 2), &
              drainage_params(idrn, 3), drainage_params(idrn, 4), drainage_params(idrn, 5), &
              drainage_params(idrn, 6)
      endif
      !
      if (stat /= 0) then
         write(logstr,'(a,i0,a,i0,a)') 'Drainage type ', drainage_type(idrn), ' requires ', npars, ' parameters !'
         call stop_sfincs(logstr, -1)
      endif
      !
   enddo
   !
   close(501)
   !
   ! Map intake / outfall points to cell indices and compute centre-to-centre
   ! distance (needed by controlled-gate types 4 and 5).
   !
   do idrn = 1, ndrn
      !
      nmq = find_quadtree_cell(xsnk(idrn), ysnk(idrn))
      if (nmq > 0) nmindrn_in(idrn)  = index_sfincs_in_quadtree(nmq)
      !
      nmq = find_quadtree_cell(xsrc_drn(idrn), ysrc_drn(idrn))
      if (nmq > 0) nmindrn_out(idrn) = index_sfincs_in_quadtree(nmq)
      !
      if (nmindrn_in(idrn) > 0 .and. nmindrn_out(idrn) > 0) then
         xsnk_tmp = z_xz(nmindrn_in(idrn))
         ysnk_tmp = z_yz(nmindrn_in(idrn))
         xsrc_tmp = z_xz(nmindrn_out(idrn))
         ysrc_tmp = z_yz(nmindrn_out(idrn))
         drainage_distance(idrn) = sqrt( (xsrc_tmp - xsnk_tmp)**2 + (ysrc_tmp - ysnk_tmp)**2 )
      endif
      !
   enddo
   !
   deallocate(xsrc_drn)
   deallocate(ysrc_drn)
   deallocate(xsnk)
   deallocate(ysnk)
   !
   if (any(nmindrn_in == 0) .or. any(nmindrn_out == 0)) then
      write(logstr,'(a)') 'Warning ! For some sink/source drainage points no matching active grid cell was found!'
      call write_log(logstr, 0)
      write(logstr,'(a)') 'Warning ! These points will be skipped, please check your input!'
      call write_log(logstr, 0)
   endif
   !
   end subroutine
   !
   !
   subroutine update_src_structures(t, dt, tloop)
   !
   ! Compute discharges through each drainage structure, accumulate them
   ! into qsrc(np) (intake: -qq, outfall: +qq), and store per-structure
   ! signed discharge in qdrain(ndrn) for his output.
   !
   ! Called AFTER update_discharges, which zeros qsrc first.
   !
   ! Atomic updates on qsrc(nm) guard against two structures (or a river
   ! and a structure) writing to the same cell under parallel execution.
   !
   use sfincs_data
   !
   implicit none
   !
   real*8  :: t
   real*4  :: dt
   real    :: tloop
   !
   integer :: count0, count1, count_rate, count_max
   integer :: idrn, nmin, nmout
   real*4  :: qq, qq0
   real*4  :: dzds, frac, wdt, zsill, zmin, zmax, mng, hgate, dfrac, tcls, topen, tclose
   !
   if (ndrn <= 0) return
   !
   call system_clock(count0, count_rate, count_max)
   !
   !$acc parallel loop present( z_volume, zs, zb, qsrc, qdrain, &
   !$acc                        nmindrn_in, nmindrn_out, &
   !$acc                        drainage_type, drainage_params, &
   !$acc                        drainage_distance, drainage_status, drainage_fraction_open ) &
   !$acc              private( nmin, nmout, qq, qq0, dzds, frac, wdt, zsill, &
   !$acc                       zmin, zmax, mng, hgate, dfrac, tcls, topen, tclose )
   !$omp parallel do &
   !$omp   private( nmin, nmout, qq, qq0, dzds, frac, wdt, zsill, &
   !$omp            zmin, zmax, mng, hgate, dfrac, tcls, topen, tclose ) &
   !$omp   schedule ( static )
   do idrn = 1, ndrn
      !
      nmin  = nmindrn_in(idrn)
      nmout = nmindrn_out(idrn)
      !
      if (nmin > 0 .and. nmout > 0) then
         !
         select case(drainage_type(idrn))
            !
            case(1)
               !
               ! Pump
               !
               qq = drainage_params(idrn, 1)
               !
            case(2)
               !
               ! Culvert (bidirectional)
               !
               if (zs(nmin) > zs(nmout)) then
                  qq =  drainage_params(idrn, 1) * sqrt(zs(nmin)  - zs(nmout))
               else
                  qq = -drainage_params(idrn, 1) * sqrt(zs(nmout) - zs(nmin))
               endif
               !
            case(3)
               !
               ! Check valve (culvert, but flow only from intake to outfall)
               !
               if (zs(nmin) > zs(nmout)) then
                  qq =  drainage_params(idrn, 1) * sqrt(zs(nmin)  - zs(nmout))
               else
                  qq = -drainage_params(idrn, 1) * sqrt(zs(nmout) - zs(nmin))
               endif
               qq = max(qq, 0.0)
               !
            case(4)
               !
               ! Controlled gate - opens when intake water level is between zmin and zmax.
               !
               wdt   = drainage_params(idrn, 1)                        ! width
               zsill = drainage_params(idrn, 2)                        ! sill elevation
               mng   = drainage_params(idrn, 3)                        ! Manning's n
               zmin  = drainage_params(idrn, 4)                        ! min water level for open
               zmax  = drainage_params(idrn, 5)                        ! max water level for open
               tcls  = drainage_params(idrn, 6)                        ! closing time (s)
               !
               dzds  = (zs(nmout) - zs(nmin)) / drainage_distance(idrn)
               frac  = drainage_fraction_open(idrn)
               hgate = max(max(zs(nmin), zs(nmout)) - zsill, 0.0)
               dfrac = dt / tcls
               !
               qq0 = qdrain(idrn) / (wdt * max(frac, 0.001))           ! previous discharge per unit width, ignoring fraction
               !
               if (drainage_status(idrn) == 0) then
                  if (zs(nmin) > zmin .and. zs(nmin) < zmax) drainage_status(idrn) = 3
               elseif (drainage_status(idrn) == 1) then
                  if (zs(nmin) <= zmin .or. zs(nmin) >= zmax) drainage_status(idrn) = 2
               endif
               !
               if (drainage_status(idrn) == 2) then
                  frac = frac - dfrac
                  if (frac < 0.0) then
                     frac = 0.0
                     drainage_status(idrn) = 0
                  endif
               elseif (drainage_status(idrn) == 3) then
                  frac = frac + dfrac
                  if (frac > 1.0) then
                     frac = 1.0
                     drainage_status(idrn) = 1
                  endif
               endif
               !
               drainage_fraction_open(idrn) = frac
               !
               qq = (qq0 - g * hgate * dzds * dt) / (1.0 + g * mng**2 * dt * abs(qq0) / hgate**(7.0 / 3.0))
               qq = qq * wdt * frac
               !
            case(5)
               !
               ! Controlled gate - schedule triggered (one open/close window).
               !
               wdt    = drainage_params(idrn, 1)                       ! width
               zsill  = drainage_params(idrn, 2)                       ! sill elevation
               mng    = drainage_params(idrn, 3)                       ! Manning's n
               tclose = drainage_params(idrn, 4)                       ! time wrt tref to close
               topen  = drainage_params(idrn, 5)                       ! time wrt tref to open
               tcls   = drainage_params(idrn, 6)                       ! closing time (s)
               !
               dzds  = (zs(nmout) - zs(nmin)) / drainage_distance(idrn)
               frac  = drainage_fraction_open(idrn)
               hgate = max(max(zs(nmin), zs(nmout)) - zsill, 0.0)
               dfrac = dt / tcls
               !
               qq0 = qdrain(idrn) / (wdt * max(frac, 0.001))
               !
               if (drainage_status(idrn) == 0) then
                  if (t >= topen) drainage_status(idrn) = 3
               elseif (drainage_status(idrn) == 1) then
                  if (t >= tclose .and. t < topen) drainage_status(idrn) = 2
               endif
               !
               if (drainage_status(idrn) == 2) then
                  frac = frac - dfrac
                  if (frac < 0.0) then
                     frac = 0.0
                     drainage_status(idrn) = 0
                  endif
               elseif (drainage_status(idrn) == 3) then
                  frac = frac + dfrac
                  if (frac > 1.0) then
                     frac = 1.0
                     drainage_status(idrn) = 1
                  endif
               endif
               !
               drainage_fraction_open(idrn) = frac
               !
               qq = (qq0 - g * hgate * dzds * dt) / (1.0 + g * mng**2 * dt * abs(qq0) / hgate**(7.0 / 3.0))
               qq = qq * wdt * frac
               !
         end select
         !
         ! Relaxation: blend new and previous discharge to damp oscillations.
         !
         qq = 1.0 / (structure_relax / dt) * qq + (1.0 - (1.0 / (structure_relax / dt))) * qdrain(idrn)
         !
         ! Limit discharge by available volume in the intake / outfall cell.
         !
         if (subgrid) then
            if (qq > 0.0) then
               qq = min(qq,  max(z_volume(nmin),  0.0) / dt)
            else
               qq = max(qq, -max(z_volume(nmout), 0.0) / dt)
            endif
         else
            if (qq > 0.0) then
               qq = min(qq,  max((zs(nmin)  - zb(nmin))  * cell_area(z_flags_iref(nmin)),  0.0) / dt)
            else
               qq = max(qq, -max((zs(nmout) - zb(nmout)) * cell_area(z_flags_iref(nmout)), 0.0) / dt)
            endif
         endif
         !
         qdrain(idrn) = qq
         !
         ! Accumulate into cell-wise qsrc. Atomic guards against multiple
         ! structures (or a river and a structure) in the same cell.
         !
         !$acc atomic update
         !$omp atomic
         qsrc(nmin)  = qsrc(nmin)  - qq
         !$acc atomic update
         !$omp atomic
         qsrc(nmout) = qsrc(nmout) + qq
         !
      endif
      !
   enddo
   !$omp end parallel do
   !$acc end parallel loop
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0 * (count1 - count0) / count_rate
   !
   end subroutine

end module
