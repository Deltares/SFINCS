module geometry
   
contains

   function check_line_through_square(xx0, yy0, xx1, yy1, x11, y11, x12, y12) result (iok)
   !
   implicit none
   !
   real*4, intent(in)   :: xx0, yy0, xx1, yy1, x11, y11, x12, y12
   real*4               :: x21,x22,y21,y22
   logical              :: iok
   !
   if (x11>xx0 .and. x11<xx1 .and. y11>yy0 .and. y11<yy1) then
      !
      ! First point
      !
      iok = .true.
      !
   elseif (x12>xx0 .and. x12<xx1 .and. y12>yy0 .and. y12<yy1) then
      !
      ! Second point
      !
      iok = .true.
      !
   else
      !
      ! lower
      !
      x21 = xx0
      x22 = xx1
      y21 = yy0
      y22 = yy0
      !    
      iok = check_intersect(x11,y11,x12,y12,x21,y21,x22,y22)
      !
      if (.not. iok) then
         !
         ! left
         !
         x21 = xx0
         x22 = xx0
         y21 = yy0
         y22 = yy1
         !        
         iok = check_intersect(x11,y11,x12,y12,x21,y21,x22,y22)
         !        
      endif
      !    
      if (.not. iok) then
         !
         ! top
         !
         x21 = xx0
         x22 = xx1
         y21 = yy1
         y22 = yy1
         !        
         iok = check_intersect(x11,y11,x12,y12,x21,y21,x22,y22)
         !        
      endif
      !
      if (.not. iok) then
         !
         ! right
         !
         x21 = xx1
         x22 = xx1
         y21 = yy0
         y22 = yy1
         !        
         iok = check_intersect(x11,y11,x12,y12,x21,y21,x22,y22)
         !        
      endif
   endif
   !
   end function

   function check_intersect(x11,y11,x12,y12,x21,y21,x22,y22) result(iok)
   !
   implicit none
   !
   real*4, intent(in)   :: x11,y11,x12,y12,x21,y21,x22,y22
   real*4               :: a1,b1,c1,a2,b2,c2,xi,yi,det
   real*4               :: dst,eps
   logical              :: iok
   !
   iok = .false.
   eps = 1.0e-1
   !
   a1 = y12 - y11
   b1 = x11 - x12
   c1 = a1*x11 + b1*y11
   !
   a2 = y22 - y21
   b2 = x21 - x22
   c2 = a2*x21 + b2*y21
   !
   det = a1*b2 - a2*b1
   !
   dst = sqrt((x12 - x11)**2 + (y12 - y11)**2)
   dst = min(dst, sqrt((x22 - x21)**2 + (y22 - y21)**2))
   !
   eps = dst/100
   !
   if (det==0.0) then
      ! parallel
   else
      !
      xi = (b2*c1 - b1*c2)/det
      yi = (a1*c2 - a2*c1)/det
      !
      if (xi>=min(x11, x12) - eps .and. xi<=max(x11, x12) + eps .and. yi>=min(y11, y12) - eps .and. yi<=max(y11, y12) + eps) then
         if (xi>=min(x21, x22) - eps .and. xi<=max(x21, x22) + eps .and. yi>=min(y21, y22) - eps .and. yi<=max(y21, y22) + eps) then
            !
            iok = .true.
            !
         endif
      endif
   endif
   !
   end function


   function direction(ax, ay, bx, by, cx, cy) result(dr)
   !
   real    :: val
   integer :: dr
   !
   val = (by - ay)*(cx - bx) - (bx - ax)*(cy - by)
   !   
   if (val<-1.0e-6) then
      dr = -1
   elseif (val>1.0e-6) then
      dr = 1
   else
      dr = 0
   endif   
   !
   end function

   function on_line(x1,y1,x2,y2,xp,yp) result(iok)
   !
   implicit none
   !
   real*4, intent(in)   :: x1,y1,x2,y2,xp,yp
   logical              :: iok
   !
   iok = .false.
   !
   if (xp <= max(x1, x2) .and. xp >= min(x1, x2) .and. yp <= max(y1, y2) .and. yp >= min(y1, y2)) then
      iok = .true.
   endif  
   !
   end function

   
   function cross(x11,y11,x12,y12,x21,y21,x22,y22) result(iok)
   !
   implicit none
   !
   real*4, intent(in)   :: x11,y11,x12,y12,x21,y21,x22,y22
   integer              :: dir1, dir2, dir3, dir4
   logical              :: iok
   !
   dir1 = direction(x11, y11, x12, y12, x21, y21)
   dir2 = direction(x11, y11, x12, y12, x22, y22)
   dir3 = direction(x21, y21, x22, y22, x11, y11)
   dir4 = direction(x21, y21, x22, y22, x12, y12)
   !
   iok = .false.
   !
   if (dir1/=dir2 .and. dir3/=dir4) then
      iok = .true.
      return
   endif   
   !
!   if (dir1 == 0 .and. on_line(x11,y11,x12,y12,x21,y21)) then
!      iok = .true.
!      return
!   endif   
!   if (dir2 == 0 .and. on_line(x11,y11,x12,y12,x22,y22)) then
!      iok = .true.
!      return
!   endif   
!   if (dir3 == 0 .and. on_line(x21,y21,x22,y22,x11,y11)) then
!      iok = .true.
!      return
!   endif   
!   if (dir4 == 0 .and. on_line(x21,y21,x22,y22,x12,y12)) then
!      iok = .true.
!      return
!   endif   
   !
   end function
   
   function prj(x1,y1,x2,y2,x,y) result (on_line)
   !
   real    :: x1, y1, x2, y2
   real    :: x, y, inner_product, dx, dy
   logical :: on_line
   !
   ! Checks if a point (x,y) can be projected onto a line segment (x1,y1)-(x2,y2)
   !
   dx = x2 - x1
   dy = y2 - y1
   inner_product =  (x - x1)*dx + (y - y1)*dy
   on_line = .false.
   if (inner_product > 0.0 .and. inner_product < dx**2 + dy**2) then
      on_line = .true.
   endif
   !
   end function prj

   
   function distance_between_points_projected_on_line_segment(x0a,y0a,x0b,y0b,x1,y1,x2,y2,mxdst) result (dst)
   !
   ! Computes distance between two points (x0a, y0a) and (x0b, y0b) projected on a line segment
   !
   real :: x0a,y0a,x0b,y0b,x1,y1,x2,y2,mxdst
   real :: b1,b2,c1,c2,aa,a1,a2,p1,p2,dst
   !
   ! First point
   !
   ! Distance of point to line
   !
   b1 = abs((x2-x1)*(y1-y0a) - (x1-x0a)*(y2-y1))/(sqrt((x2-x1)**2+(y2-y1)**2))
   c1 = sqrt((x0a-x1)**2+(y0a-y1)**2)
   c2 = sqrt((x0a-x2)**2+(y0a-y2)**2)
   aa = sqrt((x2-x1)**2+(y2-y1)**2)
   a1 = sqrt(c1**2 - b1**2)
   a2 = sqrt(c2**2 - b1**2)
   !
   if (a1>aa .or. a2>aa) then
      ! Point not on line
      if (a1<a2) then
          p1 = 0.0 
      else
          p1 = 1.0
      endif
   else
      p1 = a1/aa
   endif
   !
   ! Second point
   !
   ! Distance of point to line
   !
   b2 = abs((x2-x1)*(y1-y0b) - (x1-x0b)*(y2-y1))/(sqrt((x2-x1)**2+(y2-y1)**2))
   c1 = sqrt((x0b-x1)**2+(y0b-y1)**2)
   c2 = sqrt((x0b-x2)**2+(y0b-y2)**2)
   aa = sqrt((x2-x1)**2+(y2-y1)**2)
   a1 = sqrt(c1**2 - b2**2)
   a2 = sqrt(c2**2 - b2**2)
   !
   if (a1>aa .or. a2>aa) then
      ! Point not on line
      if (a1<a2) then
          p2 = 0.0 
      else
          p2 = 1.0
      endif
   else
      p2 = a1/aa
   endif
   !
   if (b1<mxdst .and. b2<mxdst) then
      dst = aa*abs(p2 - p1)
   else
      dst = 0.0
   endif
   !
   end function
   
   
end module
