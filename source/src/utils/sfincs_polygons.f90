module sfincs_polygons
   !
   ! Minimal polygon helper module for SFINCS.
   !
   ! Provides:
   !   t_polygon               Derived type: one named polygon with vertex arrays.
   !   read_tek_polygons       Read a Delft3D "tek" polyline/polygon file
   !                           (name line + "nrows ncols" + rows of x y) and
   !                           return an array of t_polygon, each closed so
   !                           that the last vertex equals the first. Called
   !                           from initialize_urban_drainage (sfincs_urban_drainage).
   !   point_in_polygon        Scalar ray-casting test for a single point.
   !                           Called from points_in_polygon_omp.
   !   points_in_polygon_omp   OMP-parallel sweep of a point array against one
   !                           polygon, writing a logical mask. Called from
   !                           initialize_urban_drainage (sfincs_urban_drainage).
   !
   use sfincs_log
   use sfincs_error
   !
   implicit none
   !
   private
   public :: t_polygon, read_tek_polygons, point_in_polygon, points_in_polygon_omp
   !
   type :: t_polygon
      character(len=64)                  :: name = ''
      integer                            :: n = 0
      real*4, dimension(:), allocatable  :: x
      real*4, dimension(:), allocatable  :: y
   end type t_polygon
   !
contains
   !
   subroutine read_tek_polygons(filename, polygons, ierr)
      !
      ! Read a Delft3D "tek" polyline file into an array of t_polygon.
      !
      ! File format (one or more blocks):
      !    <name line>
      !    <nrows> <ncols>
      !    <x_1> <y_1>
      !    ...
      !    <x_nrows> <y_nrows>
      !
      ! The file is swept twice: first to count blocks, then to read them.
      ! Each polygon is auto-closed: if the last vertex is not equal to the
      ! first, a copy of the first vertex is appended so downstream
      ! point-in-polygon tests can treat the last-to-first edge uniformly.
      !
      ! Called from: initialize_urban_drainage (sfincs_urban_drainage).
      !
      implicit none
      !
      character(len=*),               intent(in)  :: filename
      type(t_polygon), allocatable,   intent(out) :: polygons(:)
      integer,                        intent(out) :: ierr
      !
      integer            :: unit, stat, npoly, nrows, ncols, irow, ipoly
      character(len=256) :: name_line
      real*4             :: dummy
      logical            :: ok
      !
      ierr = 0
      !
      ok = check_file_exists(filename, 'Urban drainage polygon file', .true.)
      !
      ! First pass: count polygons.
      !
      unit = 501
      open(unit, file=trim(filename), status='old', action='read', iostat=stat)
      if (stat /= 0) then
         write(logstr,'(a,a)')' Error ! Could not open polygon file ', trim(filename)
         call write_log(logstr, 1)
         ierr = 1
         return
      endif
      !
      npoly = 0
      do
         read(unit, '(a)', iostat=stat) name_line
         if (stat /= 0) exit
         if (len_trim(name_line) == 0) cycle
         read(unit, *, iostat=stat) nrows, ncols
         if (stat /= 0) exit
         npoly = npoly + 1
         do irow = 1, nrows
            read(unit, *, iostat=stat) dummy
            if (stat /= 0) exit
         enddo
      enddo
      rewind(unit)
      !
      if (npoly == 0) then
         close(unit)
         allocate(polygons(0))
         return
      endif
      !
      allocate(polygons(npoly))
      !
      ! Second pass: read each polygon.
      !
      do ipoly = 1, npoly
         !
         read(unit, '(a)', iostat=stat) name_line
         if (stat /= 0) exit
         do while (len_trim(name_line) == 0)
            read(unit, '(a)', iostat=stat) name_line
            if (stat /= 0) exit
         enddo
         !
         read(unit, *, iostat=stat) nrows, ncols
         if (stat /= 0) then
            write(logstr,'(a,i0,a,a)')' Error ! Missing shape line for polygon ', ipoly, &
                 ' in ', trim(filename)
            call write_log(logstr, 1)
            ierr = 1
            close(unit)
            return
         endif
         !
         ! Reserve one extra slot so we can auto-close the ring if needed.
         !
         allocate(polygons(ipoly)%x(nrows + 1))
         allocate(polygons(ipoly)%y(nrows + 1))
         polygons(ipoly)%name = trim(adjustl(name_line))
         !
         do irow = 1, nrows
            read(unit, *, iostat=stat) polygons(ipoly)%x(irow), polygons(ipoly)%y(irow)
            if (stat /= 0) then
               write(logstr,'(a,i0,a,i0,a,a)')' Error ! Failed reading vertex ', irow, &
                    ' of polygon ', ipoly, ' in ', trim(filename)
               call write_log(logstr, 1)
               ierr = 1
               close(unit)
               return
            endif
         enddo
         !
         ! Close the ring if the user omitted it.
         !
         if (polygons(ipoly)%x(nrows) /= polygons(ipoly)%x(1) .or. &
             polygons(ipoly)%y(nrows) /= polygons(ipoly)%y(1)) then
            polygons(ipoly)%x(nrows + 1) = polygons(ipoly)%x(1)
            polygons(ipoly)%y(nrows + 1) = polygons(ipoly)%y(1)
            polygons(ipoly)%n = nrows + 1
         else
            polygons(ipoly)%n = nrows
         endif
         !
      enddo
      !
      close(unit)
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   pure function point_in_polygon(xp, yp, xv, yv, nv) result(inside)
      !
      ! Classic even-odd ray-casting point-in-polygon test.
      !
      ! Returns .true. if (xp, yp) lies inside the closed polygon defined by
      ! the first nv vertices of (xv, yv). The polygon is assumed closed by
      ! the caller (last vertex == first vertex).
      !
      ! Called from: points_in_polygon_omp (this module).
      !
      implicit none
      !
      real*4, intent(in)              :: xp, yp
      integer, intent(in)             :: nv
      real*4, intent(in)              :: xv(nv), yv(nv)
      logical                         :: inside
      !
      integer :: i, j
      !
      inside = .false.
      j = nv - 1
      if (j < 1) j = nv
      do i = 1, nv
         if (((yv(i) > yp) .neqv. (yv(j) > yp)) .and. &
             (xp < (xv(j) - xv(i)) * (yp - yv(i)) / (yv(j) - yv(i) + tiny(1.0)) + xv(i))) then
            inside = .not. inside
         endif
         j = i
      enddo
      !
   end function
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine points_in_polygon_omp(xp, yp, np_pts, poly, inside)
      !
      ! OMP-parallel sweep of an (xp, yp) point array against a single polygon.
      ! inside(:) must be allocated by the caller with size np_pts. Points
      ! already flagged .true. on entry are preserved (so the caller can
      ! accumulate hits across multiple polygons if desired); this routine
      ! only flips .false. to .true..
      !
      ! Called from: initialize_urban_drainage (sfincs_urban_drainage).
      !
      implicit none
      !
      integer,        intent(in)    :: np_pts
      real*4,         intent(in)    :: xp(np_pts), yp(np_pts)
      type(t_polygon), intent(in)   :: poly
      logical,        intent(inout) :: inside(np_pts)
      !
      integer :: i
      real*4  :: xmin, xmax, ymin, ymax
      !
      if (poly%n < 3) return
      !
      xmin = minval(poly%x(1:poly%n))
      xmax = maxval(poly%x(1:poly%n))
      ymin = minval(poly%y(1:poly%n))
      ymax = maxval(poly%y(1:poly%n))
      !
      !$omp parallel do default(shared) private(i) schedule(static)
      do i = 1, np_pts
         if (inside(i)) cycle
         if (xp(i) < xmin .or. xp(i) > xmax .or. yp(i) < ymin .or. yp(i) > ymax) cycle
         if (point_in_polygon(xp(i), yp(i), poly%x, poly%y, poly%n)) then
            inside(i) = .true.
         endif
      enddo
      !$omp end parallel do
      !
   end subroutine
   !
end module sfincs_polygons
