   subroutine deg2utm(la,lo,xx,yy,utmzone)
   !
   real         :: la
   real         :: lo
   real         :: xx
   real         :: yy
   real         :: pi
   real         :: sa
   real         :: sb
   real         :: e2
   real         :: e2cuadrada
   real         :: c
   real         :: lat
   real         :: lon
   real         :: deltaS
   real         :: S
   real         :: a
   real         :: epsilon
   real         :: nu
   real         :: v
   real         :: ta
   real         :: a1
   real         :: a2
   real         :: j2
   real         :: j4
   real         :: j6
   real         :: alfa
   real         :: beta
   real         :: gama
   real         :: Bm
   ! 
   integer      :: Huso
   !
   character*3  :: utmzone
   
   pi = 3.141592653589793
   sa = 6378137.000000
   sb = 6356752.314245
   !      
   read(utmzone(1:2), *)Huso   
   !
   e2 = ( ( ( sa ** 2 ) - ( sb ** 2 ) ) ** 0.5 ) / sb
   e2cuadrada = e2 ** 2
   c = ( sa ** 2 ) / sb
   !
   lat = la * ( pi / 180 )
   lon = lo * ( pi / 180 )
   !
!   Huso = int( ( lo / 6 ) + 31.0)
   S = ( ( Huso * 6.0 ) - 183.0 )
   deltaS = lon - ( S * ( pi / 180 ) )
   !
   a = cos(lat) * sin(deltaS)
   epsilon = 0.5 * log( ( 1.0 +  a) / ( 1.0 - a ) )
   nu = atan( tan(lat) / cos(deltaS) ) - lat
   v = ( c / ( ( 1.0 + ( e2cuadrada * ( cos(lat) ) ** 2 ) ) ) ** 0.5 ) * 0.9996
   ta = ( e2cuadrada / 2 ) * epsilon ** 2 * ( cos(lat) ) ** 2
   a1 = sin( 2 * lat )
   a2 = a1 * ( cos(lat) ) ** 2
   j2 = lat + ( a1 / 2 )
   j4 = ( ( 3 * j2 ) + a2 ) / 4
   j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat) ) ** 2) ) / 3
   alfa = ( 3.0 / 4.0 ) * e2cuadrada
   beta = ( 5.0 / 3.0 ) * alfa ** 2
   gama = ( 35.0 / 27.0 ) * alfa ** 3
   Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 )
   xx = epsilon * v * ( 1 + ( ta / 3 ) ) + 500000
   yy = nu * v * ( 1.0 + ta ) + Bm
   !
   if (yy<0) then
       yy = 9999999.0 + yy
   endif
   !
!   utm_zone.number(i,:) = num2str(Huso);
!   utm_zone.letter(i,:) = Letra;
!   utm_zone.NS(i,:)     = northsouth;
   !
   end subroutine
    