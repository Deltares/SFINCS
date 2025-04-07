! Non-hydrostatic code now only works with regular grids (can still use quadtree netcdf file as long as there are no refinement levels).
! Still to do: add nonh_mask to netcdf file. Now the non-hydrostatic corrections are applied to the entire SFINCS grid.
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
   !$omp private ( ip, ipuv, nm, nmu, nhnm, nhnmu, hu, unh )
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
         hu = max(zs(nm), zs(nmu)) - 0.5 * (zb(nm) + zb(nmu))
         !
         if (hu > huthresh_nh) then
            !
            unh = - 1.0 * dtover2rhodx * ( AB(ip) * (pnh(nhnmu) + pnh(nhnm)) + pnh(nhnmu) - pnh(nhnm) )
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
         hnb = zs(nmn) - zb(nmn)   ! water depth at neighbor
         wb(irow) = wb(irow) - 0.5 * uv(iuv) * (hnm - hnb) * dxrinv(1)
      endif
      !
      ! Right
      !
      iuv = z_index_uv_mu(nm)      ! uv index
      !
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nmu(iuv) ! nm index of neighbor
         hnb = zs(nmn) - zb(nmn)   ! water depth at neighbor
         wb(irow) = wb(irow) - 0.5 * uv(iuv) * (hnb - hnm) * dxrinv(1)
      endif
      !
      ! Bottom
      !      
      iuv = z_index_uv_nd(nm)     ! uv index
      !
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nm(iuv) ! nm index of neighbor
         hnb = zs(nmn) - zb(nmn)  ! water depth at neighbor
         wb(irow) = wb(irow) - 0.5 * uv(iuv) * (hnm - hnb) * dyrinv(1)
      endif
      !
      ! Top
      !
      iuv = z_index_uv_nu(nm)      ! uv index
      !
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nmu(iuv) ! nm index of neighbor
         hnb = zs(nmn) - zb(nmn)   ! water depth at neighbor
         wb(irow) = wb(irow) - 0.5 * uv(iuv) * (hnb - hnm) * dyrinv(1)
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

end module
