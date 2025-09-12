module sfincs_boundaries
   
   use sfincs_log
   use sfincs_error

contains

   subroutine read_boundary_data()
   !
   ! Reads bnd, bzs etc. files
   !
   use sfincs_ncinput
   use sfincs_data
   !
   implicit none
   !
   integer n, itb, ib, stat, ifreq, iok
   logical :: ok
   !
   real*4 dummy,r
   !
   ! Read water level boundaries
   !
   nbnd      = 0
   ntbnd     = 0
   itbndlast = 1
   !
   if (netbndbzsbzifile(1:4) /= 'none') then    ! FEWS compatible Netcdf water level time-series input
      !
      call read_netcdf_boundary_data()
      !
      if ((t_bnd(1) > (t0 + 1.0)) .or. (t_bnd(ntbnd) < (t1 - 1.0))) then
         !
         write(logstr,'(a)')' WARNING! Times in boundary conditions file do not cover entire simulation period!'
         call write_log(logstr, 1)
         !
      endif   
      !
   elseif (bndfile(1:4) /= 'none') then    ! Normal ascii input files
      !
      call write_log('Info    : reading water level boundaries', 0)
      !
      open(500, file=trim(bndfile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         nbnd = nbnd + 1
      enddo
      rewind(500)
      allocate(x_bnd(nbnd))
      allocate(y_bnd(nbnd))
      do n = 1, nbnd
         read(500,*)x_bnd(n),y_bnd(n)
      enddo
      close(500)
      !
      ! Read water level boundary conditions file
      !
      open(500, file=trim(bzsfile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         ntbnd = ntbnd + 1
      enddo
      rewind(500)
      !
      allocate(t_bnd(ntbnd))
      allocate(zs_bnd(nbnd,ntbnd))
      allocate(zst_bnd(nbnd))
      !
      do itb = 1, ntbnd
         read(500,*)t_bnd(itb),(zs_bnd(ib, itb), ib = 1, nbnd)
      enddo
      !
      close(500)
      !
      if (bzifile(1:4) /= 'none') then
         !
         ! Incoming infragravity waves
         !
         allocate(zsi_bnd(nbnd,ntbnd))
         allocate(zsit_bnd(nbnd))
         !
         open(500, file=trim(bzifile))
         do itb = 1, ntbnd
            read(500,*)dummy,(zsi_bnd(ib, itb), ib = 1, nbnd)
         enddo
         close(500)
         !
      endif
      !      
      if ((t_bnd(1) > (t0 + 1.0)) .or. (t_bnd(ntbnd) < (t1 - 1.0))) then
         ! 
         write(logstr,'(a)')'Warning! Times in boundary conditions file do not cover entire simulation period !'
         call write_log(logstr, 1)
         !
         if (t_bnd(1) > (t0 + 1.0)) then
            ! 
            write(logstr,'(a)')'Warning! Adjusting first time in boundary conditions time series !'
            call write_log(logstr, 1)
            !
            t_bnd(1) = t0 - 1.0
            !
         else
            ! 
            write(logstr,'(a)')'Warning! Adjusting last time in boundary conditions time series !'
            call write_log(logstr, 1)
            !
            t_bnd(ntbnd) = t1 + 1.0
            !
         endif
         !
      endif   
      !
   elseif (water_level_boundaries_in_mask) then   
      !
      write(logstr,'(a)')'Warning! Boundary cells found in mask, without boundary conditions. Using water level of 0.0 m at these points.'
      call write_log(logstr, 1)
      !
      nbnd = 1
      allocate(x_bnd(nbnd))
      allocate(y_bnd(nbnd))
      x_bnd(1) = 0.0
      y_bnd(1) = 0.0
      !
      ntbnd = 2
      !
      allocate(t_bnd(ntbnd))
      allocate(zs_bnd(nbnd,ntbnd))
      allocate(zst_bnd(nbnd))
      !
      t_bnd(1) = t0
      t_bnd(2) = t1
      zs_bnd(1,1) = 0.0
      zs_bnd(2,1) = 0.0
      zs_bnd(1,2) = 0.0
      zs_bnd(2,2) = 0.0
      zst_bnd(1)  = 0.0      
      !
   endif
   !
   ! Check for 'weird' values
   !
   iok = 1
   ! 
   do ib = 1, nbnd
      do itb = 1, ntbnd
         !
         if (zs_bnd(ib, itb)<-99.0 .or. zs_bnd(ib, itb)>990.0) then
            !
            iok = 0
            !
            zs_bnd(ib, itb) = 0.0
            !
         endif
         !
      enddo 
   enddo    
   !
   if (iok == 0) then
      ! 
      write(logstr,'(a)')'Warning! Very low, very high or NaN values found in boundary conditions file ! These have now been replaced with zeros. Please check !'
      call write_log(logstr, 1)
      ! 
   endif
   !
   ! Now the downstream river boundaries. No time series, just points, slope and direction
   !
   nbdr = 0
   !
   if (downstream_river_boundaries_in_mask) then
      !      
      if (bdrfile(1:4) /= 'none') then
         !
         write(logstr,'(a)')'Info    : reading downstream river boundaries'
         call write_log(logstr, 0)
         !
         open(500, file=trim(bdrfile))
         do while(.true.)
            read(500,*,iostat = stat)dummy
            if (stat<0) exit
            nbdr = nbdr + 1
         enddo
         rewind(500)
         allocate(x_bdr(nbdr))
         allocate(y_bdr(nbdr))
         allocate(slope_bdr(nbdr))
         allocate(azimuth_bdr(nbdr))      
         do n = 1, nbdr
            read(500,*)x_bdr(n), y_bdr(n), slope_bdr(n), azimuth_bdr(n)
         enddo
         close(500)
         !
      else
         !
         ! This really should not happen
         !
         call stop_sfincs('Error! Downstream river points found in mask, without boundary conditions (missing bdrfile in sfincs.inp) !', 1)
         !
      endif
      !
   else   
      !
      if (bdrfile(1:4) /= 'none') then
         !
         write(logstr,'(a)')'Warning : Found bdr file in sfincs.inp, but no downstream river points in were found in the mask!'
         call write_log(logstr, 0)
         !  
      endif
      !
   endif   
   !
   end subroutine


   subroutine find_boundary_indices()
   !
   use sfincs_data
   !
   implicit none
   !
   ! For each grid boundary point (kcs=2, ) :
   !
   ! Determine indices and weights of boundary points
   ! For tide and surge, these are the indices and weights of the points in the bnd file
   !
   integer nm, m, n, nb, ib1, ib2, ib, ic, ibnd, ibdr, iref
   integer nm_1, nm_2, nm_3, nmi
   !
   real x, y, dst1, dst2, dst, phi_riv, phi_r
   real*4 :: w_1, w_2, w_3, d_1, d_2, d_3 
   !
   ! Check there are points from bnd file and/or bdr file and that kcs mask contains 2/3/5/6
   !
   if (ngbnd == 0) then
      !
      if (nbnd > 0 .or. nbdr > 0) then
         !
         write(logstr,'(a)')'Warning : no open boundary points found in mask!'
         call write_log(logstr, 1)
         !
      endif
      !
      return
      !
   endif   
   !
   if (nbnd == 0 .and. nbdr == 0) then
      !
      !write(logstr,'(a)')'Warning : no open boundary points found in mask!'
      !call write_log(logstr, 1)
      !
      return
      !
   endif
   !   
   ! Allocate boundary arrays
   !
   allocate(ind1_bnd_gbp(ngbnd))
   allocate(ind2_bnd_gbp(ngbnd))
   allocate(fac_bnd_gbp(ngbnd))
   !
   ind1_bnd_gbp = 0
   ind2_bnd_gbp = 0
   fac_bnd_gbp  = 0.0
   !
   if (downstream_river_boundaries_in_mask) then
      !
      ! There are downstream river boundaries, so we need indices, weights and distances of upstream points
      !
      allocate(slope_gbp(ngbnd))
      allocate(nm_nbr_gbp(3, ngbnd))
      allocate(w_nbr_gbp(3, ngbnd))
      allocate(d_nbr_gbp(3, ngbnd))
      !
      nm_nbr_gbp = 0
      w_nbr_gbp  = 0.0
      d_nbr_gbp  = 0.0
      !
   endif
   !
   if (neumann_boundaries_in_mask) then
      !
      ! There are Neumann boundaries, so we need nm indices of internal points
      !
      allocate(nmi_gbp(ngbnd))
      nmi_gbp = 0
      !
   endif
   !
   ! Find two closest boundary condition points for each boundary point
   !
   nb = 0
   !
   ! Loop through all grid boundary points
   !
   do ib = 1, ngbnd
      !
      nm = nmindbnd(ib)
      !
      x = z_xz(nm)
      y = z_yz(nm)
      !
      if (kcs(nm) == 2) then ! This cell is a water level boundary point
         !
         if (nbnd > 1) then
            !
            ! Multiple points in bnd file, so use distance-based weighting
            !
            dst1 = 1.0e10
            dst2 = 1.0e10
            ib1 = 0
            ib2 = 0
            !
            ! Loop through all water level boundary points in bnd file
            !
            do ibnd = 1, nbnd
               !
               ! Compute distance of this point to grid boundary point
               !
               dst = sqrt((x_bnd(ibnd) - x)**2 + (y_bnd(ibnd) - y)**2)
               !
               if (dst < dst1) then
                  !
                  ! Nearest point found
                  !
                  dst2 = dst1
                  ib2  = ib1
                  dst1 = dst
                  ib1  = ibnd
                  !
               elseif (dst < dst2) then
                  !
                  ! Second nearest point found
                  !
                  dst2 = dst
                  ib2  = ibnd
                  !
               endif
            enddo
            !
            ind1_bnd_gbp(ib) = ib1
            ind2_bnd_gbp(ib) = ib2
            fac_bnd_gbp(ib)  = dst2 / max(dst1 + dst2, 1.0e-9)
            !
         else
            !
            ind1_bnd_gbp(ib)  = 1
            ind2_bnd_gbp(ib)  = 1
            fac_bnd_gbp(ib)   = 1.0
            !
         endif
         !
      elseif (kcs(nm) == 5) then  ! This cell is a downstream river boundary point
         !
         ! We just look up nearest point in bdr file
         !
         dst1 = 1.0e10
         ib1 = 0
         !
         ! Loop through all points in bdr file
         !
         do ibdr = 1, nbdr
            !
            ! Compute distance of this point to grid boundary point
            !
            dst = sqrt((x_bdr(ibdr) - x)**2 + (y_bdr(ibdr) - y)**2)
            !
            if (dst < dst1) then
               !
               ! Nearest point found
               !
               dst1 = dst
               ib1  = ibdr
               !
            endif
            !
         enddo
         !
         ! Now we need to determine info on how to obtain zs for this cell at each time step based on slope and direction in bdr file.
         ! The problem is that a grid boundary cell can have up to three internal neighbors inside the domain.
         ! Of each of these points, we need: nm index, projected distance to line perpendicular to river, and weight.
         ! In case of three neighbors, zsb will then be computed in update_boundary_conditions as:
         !
         ! zsb = w1 * (zs(i1) - slope*dx1) + w2 * (zs(i2) - slope*dx2) + w3 * (zs(i3) - slope*dx3)
         !
         ! Quadtree refinement is not allowed near river outflow boundaries! This would get too complicated for now. Also, grid cells must be square.
         !
         slope_gbp(ib) = slope_bdr(ib1)
         !
         ! Indices of boundary cell
         !
         n = z_index_z_n(nm) 
         m = z_index_z_m(nm)
         iref = z_flags_iref(nm)
         !
         ! Virtually rotate grid to 0.0 and river direction (where the river is flowing from) accordingly
         !
         phi_riv = modulo(2 * pi * (90.0 - azimuth_bdr(ib1)) - rotation, 2 * pi)
         !
         if (phi_riv < 0.5 * pi) then
            !
            ! Neighbors in NE
            !
            nm_1 = find_sfincs_cell(n    , m + 1, iref)
            nm_2 = find_sfincs_cell(n + 1, m + 1, iref)
            nm_3 = find_sfincs_cell(n + 1, m    , iref)
            !
            phi_r = phi_riv
            !
         elseif (phi_riv < pi) then
            !
            ! Neighbors in NW
            !
            nm_1 = find_sfincs_cell(n + 1, m    , iref)
            nm_2 = find_sfincs_cell(n + 1, m - 1, iref)
            nm_3 = find_sfincs_cell(n    , m - 1, iref)
            !
            phi_r = phi_riv - 0.5 * pi
            !
         elseif (phi_riv < 1.5 * pi) then
            !
            ! Neighbors in SE
            !
            nm_1 = find_sfincs_cell(n    , m - 1, iref)
            nm_2 = find_sfincs_cell(n - 1, m - 1, iref)
            nm_3 = find_sfincs_cell(n - 1, m    , iref)
            !
            phi_r = phi_riv - 1.0 * pi
            !
         else
            !
            ! Neighbors in SW
            !
            nm_1 = find_sfincs_cell(n - 1, m    , iref)
            nm_2 = find_sfincs_cell(n - 1, m + 1, iref)
            nm_3 = find_sfincs_cell(n    , m + 1, iref)
            !
            phi_r = phi_riv - 1.5 * pi
            !
         endif
         !
         ! Now compute weights
         !
         w_1 = cos(phi_r) / (cos(phi_r) + sin(2 * phi_r) + sin(phi_r))
         w_2 = sin(2 * phi_r) / (cos(phi_r) + sin(2 * phi_r) + sin(phi_r))
         w_3 = sin(phi_r) / (cos(phi_r) + sin(2 * phi_r) + sin(phi_r))
         !
         ! And the distances
         !
         d_1 = cos(phi_r) * dxr(iref)
         d_2 = (1.0 + (sqrt(2.0) - 1.0) * sin(2 * phi_r)) * dxr(iref)
         d_3 = sin(phi_r) * dxr(iref)
         !
         ! Fill nbr arrays
         !
         nm_nbr_gbp(1, ib) = nm_1
         nm_nbr_gbp(2, ib) = nm_2
         nm_nbr_gbp(3, ib) = nm_3
         !
         w_nbr_gbp(1, ib) = w_1
         w_nbr_gbp(2, ib) = w_2
         w_nbr_gbp(3, ib) = w_3
         !
         d_nbr_gbp(1, ib) = d_1
         d_nbr_gbp(2, ib) = d_2
         d_nbr_gbp(3, ib) = d_3
         !
      elseif (kcs(nm) == 6) then  ! This cell is a Neumann boundary point
         !
         ! Indices of boundary cell
         !
         n = z_index_z_n(nm) 
         m = z_index_z_m(nm)
         iref = z_flags_iref(nm)
         !
         ! Get index of internal point
         !
         nmi = find_sfincs_cell(n, m + 1, iref)
         if (nmi > 0) then
            if (kcs(nmi) == 1) nmi_gbp(ib) = nmi
         endif
         !
         nmi = find_sfincs_cell(n + 1, m, iref)
         if (nmi > 0) then 
            if (kcs(nmi) == 1) nmi_gbp(ib) = nmi
         endif
         !
         nmi = find_sfincs_cell(n, m - 1, iref)
         if (nmi > 0) then 
            if (kcs(nmi) == 1) nmi_gbp(ib) = nmi
         endif
         !
         nmi = find_sfincs_cell(n - 1, m, iref)
         if (nmi > 0) then 
            if (kcs(nmi) == 1) nmi_gbp(ib) = nmi
         endif
         !
      endif
   enddo
   !
   end subroutine


   subroutine update_boundary_points(t)
   !
   ! Update values at boundary points
   !
   use sfincs_data
   !
   implicit none
   !
   integer ib, itb, itb0, itb1
   !
   real*8 t
   !
   real*4 zstb, tbfac, hs, tp, wd, tb
   !
   if (nbnd == 0) return
   !
   ! Interpolate boundary conditions in time
   !
   if (t_bnd(1) > (t - 1.0e-3)) then ! use first time in boundary conditions
      !
      itb0 = 1
      itb1 = 1
      tb   = t_bnd(itb0)
      !
   elseif (t_bnd(ntbnd) < (t + 1.0e-3)) then  ! use last time in boundary conditions       
      !
      itb0 = ntbnd
      itb1 = ntbnd
      tb   = t_bnd(itb0)
      !
   else
      !
      do itb = itbndlast, ntbnd ! Loop in time
         if (t_bnd(itb) > (t + 1.0e-6)) then
            itb0 = itb - 1
            itb1 = itb
            tb   = t
            itbndlast = itb - 1
            exit
         endif
      enddo 
      !
   endif            
   !
   tbfac  = (tb - t_bnd(itb0))/max(t_bnd(itb1) - t_bnd(itb0), 1.0e-6)
   !
   do ib = 1, nbnd ! Loop along boundary points
      !
      ! Tide and surge
      !
      zstb = zs_bnd(ib, itb0) + (zs_bnd(ib, itb1) - zs_bnd(ib, itb0))*tbfac
      !
      zst_bnd(ib) = zstb
      !
      if (bzifile(1:4) /= 'none') then
         !
         ! Incoming infragravity waves
         !
         zsit_bnd(ib) = zsi_bnd(ib, itb0) + (zsi_bnd(ib, itb1) - zsi_bnd(ib, itb0))*tbfac
         !
      endif
      !
   enddo
   !
   end subroutine
   
   subroutine update_boundary_conditions(t, dt)
   !
   ! Update water level at boundary grid points
   !
   use sfincs_data
   !
   implicit none
   !
   integer ib, ifreq, ic, nm, ibuv, nmi, nmb, indb, iw
   !
   real*8                            :: t
   real                              :: dt
   !
   real*8 zst
   real*4 hst, l0t, wdt, ksi, zsetup, zig, a, angdiff, angfac, hm0ig, tmmin01ig, fp
   real*4 zs0act
   real*4 smfac
   real*4 zs0smooth
   real*4 sumw
   !
   ! Set water level in all boundary points on grid
   !
   do ib = 1, ngbnd
      !
      nmb = nmindbnd(ib)
      !
      ! kcs = 1 : regular point
      ! kcs = 2 : water level boundary point
      ! kcs = 3 : outflow boundary point
      ! kcs = 4 : wave maker point
      ! kcs = 5 : river outflow point (dzs/dx = i)
      ! kcs = 6 : lateral (coastal) boundary point (Neumann dzs/dx = 0.0)
      !
      if (kcs(nmb) == 2) then
         !
         ! Regular water level boundary point
         !
         ! Get water levels (surge + tide) from time series boundary conditions
         !
         if (nbnd > 1) then
            !
            ! Interpolation of nearby points
            !
            zst   = zst_bnd(ind1_bnd_gbp(ib)) * fac_bnd_gbp(ib) + zst_bnd(ind2_bnd_gbp(ib)) * (1.0 - fac_bnd_gbp(ib))
            !
         else
            !
            ! Just use the value of the one boundary point
            !
            zst   = zst_bnd(1)
            !
         endif
         !
         if (patmos .and. pavbnd>1.0) then
            !
            ! Barometric pressure correction
            !
            zst = zst + (pavbnd - patmb(ib)) / (rhow * 9.81)
            !
         endif
         !
         zsetup = 0.0
         zig    = 0.0
         !
         ! Incoming IG waves from file (this will overrule IG signal computed before)
         !
         if (bzifile(1:4) /= 'none') then
            !
            if (nbnd>1) then
               !
               ! Interpolation of nearby points
               !
               zig   = zsit_bnd(ind1_bnd_gbp(ib)) * fac_bnd_gbp(ib)  + zsit_bnd(ind2_bnd_gbp(ib)) * (1.0 - fac_bnd_gbp(ib))
               !
            else
               !
               ! Just use the value of the one boundary point
               !
               zig   = zsit_bnd(1)
               !
            endif
            !
         endif
         !
         if (t < (tspinup - 1.0e-3)) then
            !
            smfac = 1.0 - (t - t0) / (tspinup - t0)
            !
            zs0act = zst + zsetup
            call weighted_average(zini, zs0act, smfac, 1, zs0smooth)
            zsb0(ib) = zs0smooth           ! Still water level at kcs=2 point
            !
            zs0act = zst + zsetup + zig
            call weighted_average(zini, zs0act, smfac, 1, zs0smooth)
            zsb(ib) = zs0smooth            ! Total water level at kcs=2 point
            !
         else
            !
            zsb0(ib) = zst + zsetup        ! Still water level at kcs=2 point
            zsb(ib)  = zst + zsetup + zig  ! Total water level at kcs=2 point
            !
         endif
         !
         if (subgrid) then                  ! Check on waterlevels minimally equal to z_zmin
            zsb0(ib) = max(zsb0(ib), subgrid_z_zmin(nmb))
            zsb(ib)  = max(zsb(ib),  subgrid_z_zmin(nmb))
         else                               ! Check on waterlevels minimally equal to zb           
            zsb0(ib) = max(zsb0(ib), zb(nmb))
            zsb(ib)  = max(zsb(ib),  zb(nmb))
         endif
         !
      elseif (kcs(nmb) == 5) then
         !
         ! Downstream river point
         !
         ! Get water levels from inside model, and adjust for slope.
         !
         zst = 0.0
         sumw = 0.0
         !         
         do iw = 1, 3
            !
            nm = nm_nbr_gbp(iw, ib) ! nm index of internal neighbor
            !
            if (nm > 0) then
               !
               ! Check that this point is not dry. If so, skip.
               !
               if (subgrid) then
                  if (zs(nm) < subgrid_z_zmin(nmb) + 0.01) cycle
               else
                  if (zs(nm) < zb(nmb) + huthresh) cycle
               endif
               !
               zst = zst + w_nbr_gbp(iw, ib) * (zs(nm) - slope_gbp(ib) * d_nbr_gbp(iw, ib))
               sumw = sumw + w_nbr_gbp(iw, ib)
               !
            endif
            !
         enddo
         !
         if (sumw > 1.0e-6) then
            !
            zst = zst / sumw
            !
         else
            !
            ! No wet upstream point found, so set water level equal to bed level
            !
            zst = -99999.0
            !
         endif   
         !
         ! Make sure water level is not below bed level
         !
         if (subgrid) then
            zst = max(zst, subgrid_z_zmin(nmb))
         else
            zst = max(zst, zb(nmb))
         endif
         !
         zsb(ib) = zst
         zsb0(ib) = zst
         !
      elseif (kcs(nmb) == 6) then
         !
         ! Lateral coastal (Neumann) boundary
         !
         ! Set water level at boundary point equal to water level inside model.
         ! No need to set zsb and zsb0, as flux for this type of boundary is solved in sfincs_momentum.f90.
         ! Lateral boundary u/v points have kcuv=6. They are skipped in update_boundary_fluxes.
         !
         zs(nmb) = zs(nmi_gbp(ib)) ! nm index of internal point. Technically there can be more than one internal point. This always uses the last point that was found.
         !         
      endif
      !
   enddo
   !
   end subroutine


   subroutine update_boundary_fluxes(dt)
   !
   ! Update fluxes qx and qy at boundary points
   !
   use sfincs_data
   !
   implicit none
   !
   integer ib, nm, nmi, nmb, iuv, indb, ip
   real*4  hnmb, dt, zsnmi, zsnmb, zs0nmb, facrel
   real*4  factime, one_minus_factime
   !
   real*4 ui, ub, dzuv, facint, zsuv, depthuv
   !
   !$acc update device( zsb0, zsb ), async(1)
   !
   factime = min(dt/btfilter, 1.0)
   one_minus_factime = 1.0 - factime
   facrel  = 1.0 - min(dt/btrelax, 1.0)
   !
   ! UV fluxes at boundaries
   !
   !$acc kernels present(index_kcuv2, nmikcuv2, nmbkcuv2, ibkcuv2, kcuv, zs, z_volume, q, uvmean, uv, zb, zbuv, zsb, zsb0, &
   !$acc                 subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_havg_zmax, subgrid_z_zmin, ibuvdir, zsmax ), async(1)
   !$acc loop independent, private(ib)
   do ib = 1, nkcuv2
      !
      indb   = ibkcuv2(ib)
      !
      nmb    = nmbkcuv2(ib)     ! nm index of kcs=2/3/5/6 boundary point
      !
      if (kcs(nmb) == 6) cycle  ! Lateral boundary point. Fluxes computed in sfincs_momentum.f90
      !
      nmi    = nmikcuv2(ib)     ! Index of kcs=1 point
      !
      ip     = index_kcuv2(ib)  ! Index in uv array of kcuv=2 velocity point
      !
      zsnmi  = zs(nmi)          ! total water level inside model
      zsnmb  = zsb(indb)        ! total water level at boundary
      zs0nmb = zsb0(indb)       ! average water level inside model (without waves)
      !
      zsnmb  = zsb(indb)     ! total water level at boundary
      zs0nmb = zsb0(indb)    ! average water level inside model (without waves)
      !
      if (bndtype == 1) then
         !
         ! Weakly reflective boundary (default)
         !          
         if (subgrid) then
            !
            zsuv = max(zsnmb, zsnmi)
            !
            if (zsuv>=subgrid_uv_zmax(ip) - 1.0e-4) then
               !
               ! Entire cell is wet, no interpolation from table needed
               !
               depthuv  = subgrid_uv_havg_zmax(ip) + zsuv
               !
            elseif (zsuv>subgrid_uv_zmin(ip)) then
               !
               ! Interpolation required
               !
               dzuv    = (subgrid_uv_zmax(ip) - subgrid_uv_zmin(ip)) / (subgrid_nlevels - 1)
               iuv     = int((zsuv - subgrid_uv_zmin(ip))/dzuv) + 1
               facint  = (zsuv - (subgrid_uv_zmin(ip) + (iuv - 1)*dzuv) ) / dzuv
               depthuv = subgrid_uv_havg(iuv, ip) + (subgrid_uv_havg(iuv + 1, ip) - subgrid_uv_havg(iuv, ip))*facint
               !
            else
               !
               depthuv = 0.0
               !
            endif
            !
            hnmb   = max(depthuv, huthresh)
            zsnmb  = max(zsnmb,  subgrid_z_zmin(nmb))
            zs0nmb = max(zs0nmb, subgrid_z_zmin(nmb))
            !
         else
            !
            hnmb   = max(0.5 * (zsnmb + zsnmi) - zbuv(ip), huthresh)
            zsnmb  = max(zsnmb,  zb(nmb))
            zs0nmb = max(zs0nmb, zb(nmb))
            !
         endif      
         !
         if (hnmb < huthresh + 1.0e-6 .or. kcuv(ip) == 3) then
            !
            ! Very shallow or also a structure point.
            !
            q(ip)      = 0.0
            uv(ip)     = 0.0
            uvmean(ib) = 0.0
            !
         else
            !
            ui = sqrt(g / hnmb) * (zsnmb - zs0nmb)
            ub = ibuvdir(ib) * (2 * ui - sqrt(g / hnmb) * (zsnmi - zs0nmb))
            !
            q(ip) = ub*hnmb + uvmean(ib)            
            !
            if (subgrid) then
               !
               ! Sub-grid
               !
               if (z_volume(nmi) <= 0.0) then
                  if (ibuvdir(ib) == 1) then
                     q(ip) = max(q(ip), 0.0) ! Nothing can flow out
                  else
                     q(ip) = min(q(ip), 0.0) ! Nothing can flow out
                  endif
               endif
               !
               if (zsnmb - subgrid_z_zmin(nmb) < huthresh) then
                  if (ibuvdir(ib) == 1) then
                     q(ip) = min(q(ip), 0.0) ! Nothing can flow in
                  else
                     q(ip) = max(q(ip), 0.0) ! Nothing can flow in
                  endif
               endif
               !
            else
               !
               ! Regular
               !
               if (zsnmi - zb(nmi)<=huthresh) then                  
                  if (ibuvdir(ib)==1) then
                     q(ip) = max(q(ip), 0.0) ! Nothing can flow out
                  else
                     q(ip) = min(q(ip), 0.0) ! Nothing can flow out
                  endif
               endif
               !
               if (zsnmb - zb(nmb)<huthresh) then
                  if (ibuvdir(ib)==1) then
                     q(ip) = min(q(ip), 0.0) ! Nothing can flow in
                  else
                     q(ip) = max(q(ip), 0.0) ! Nothing can flow in
                  endif
               endif
            endif
            !
            ! Limit velocities (this does not change fluxes, but may prevent advection term from exploding in the next time step)
            !
            uv(ip)  = max(min(q(ip)/hnmb, 4.0), -4.0)
            !
         endif
         !
         if (btfilter>=-1.0e-6) then
            !
            ! Added a little bit of relaxation in uvmean to avoid persistent jets shooting into the model
            ! Using: facrel = 1.0 - min(dt/btrelax, 1.0)
            ! Default: btrelax = 3600 s.
            !
            uvmean(ib) = factime * q(ip) + facrel * one_minus_factime * uvmean(ib)
            !
         else
            !
            uvmean(ib) = 0.0
            !
         endif
         !
         ! Set value on kcs=2 point to boundary condition
         !
         zs(nmb) = zsb(indb)
         !
         ! Store maximum water levels also on the boundary
         !
         if (store_maximum_waterlevel) then
            zsmax(nmb) = max(zsmax(nmb), zs(nmb))
         endif
         !
      endif
      !
   enddo
   !
   !$acc end kernels
   !
   end subroutine

   
   
   subroutine update_boundaries(t, dt, tloop)
   !
   ! Update all boundary conditions
   !
   use sfincs_data
   !
   implicit none
   !
   integer  :: count0
   integer  :: count1
   integer  :: count_rate
   integer  :: count_max
   real     :: tloop
   !
   real*8           :: t
   real*4           :: dt
   !
   call system_clock(count0, count_rate, count_max)
   !
   if (boundaries_in_mask) then
      !
      if (nbnd > 0) then
         !
         ! Update boundary conditions at boundary points from time series
         !
         call update_boundary_points(t)
         !
      endif
      !
      ! Update boundary conditions at grid points (water levels)
      !
      call update_boundary_conditions(t, dt)
      !
      ! Update boundary fluxes()
      !
      call update_boundary_fluxes(dt)
      !
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine
   !
   !
   !
   subroutine weighted_average(val1,val2,fac,iopt,val3)
   !
   implicit none
   !
   integer, intent(in)    :: iopt
   real*4,  intent(in)    :: val1
   real*4,  intent(in)    :: val2
   real*4,  intent(in)    :: fac
   real*4,  intent(out)   :: val3
   !
   real*4                 :: u1
   real*4                 :: v1
   real*4                 :: u2
   real*4                 :: v2
   real*4                 :: u
   real*4                 :: v
   !
   if (iopt==1) then
      !
      ! Regular
      !
      val3 = val1*fac  + val2*(1.0 - fac)
      !
   else
      !
      ! Angles (input must be in radians!)
      !
      u1 = cos(val1)
      v1 = sin(val1)
      u2 = cos(val2)
      v2 = sin(val2)
      !
      u = u1*fac  + u2*(1.0 - fac)
      v = v1*fac  + v2*(1.0 - fac)
      !
      val3 = atan2(v, u)
      !
   endif
   !
   end subroutine

   
   function find_sfincs_cell(n, m, iref) result (nm)
   !
   ! Find nm index for cell n, m, iref
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   integer, intent(in)  :: n
   integer, intent(in)  :: m
   integer, intent(in)  :: iref
   integer              :: nm
   !
   integer :: nmq
   !
   nmq = find_quadtree_cell_by_index(n, m, iref)
   nm = index_sfincs_in_quadtree(nmq) 
   !
   end function
   
   
end module
