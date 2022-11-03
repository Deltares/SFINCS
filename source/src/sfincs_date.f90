module sfincs_date

   
  IMPLICIT NONE

  ! Routines to handle Julian Dates and Times
  !
  ! References:
  !
  !   Montenbruck, Oliver, "Practical Ephmeris Calculations", Ch. 2, pp 33-40.
  !   The 1992 Astronomical Almanac, page B4.
  !
  ! The Julian Date is defined as the the number of days which have elapsed
  ! since the 1st of January of the year 4713 BC 12:00 Universal Time.
  !
  ! Up to 4th October 1582 AD the Julian Calendar was in force with leap
  ! years every four years, but from then on the Gregorian Calendar carried 
  ! on from 15th October 1582 with the leap years defined by the rule:
  !
  !  "Leap year is every year whose yearly number is divisible by four, but
  !   not by a hundred, or is divisible by four hundred."
  !
  ! At midday on 4th October 1582, 2,299,160.5 Julian days had elapsed.
  ! The Gregorian Calendar then carried on at this point from 15th October
  ! 1582, so its beginning occured on the Julian date 2,299,160.5.
  ! 
  ! Note: the astronomical year -4712 corresponds to the year 4713 BC, the
  !       year 0 to the year 1 BC; thereafter the astronomical year match
  !       the year AD, e.g. year 1 = year 1 AD.
  !
  !       This routines work for the years -5877402 BC until 5868098 AD. This 
  !       dates are neveetheless coverable by the current GRIB edition 1. Which
  !       can cover dates between 1 AD and 25599 AD.  
  !
  !       The "Modified Julian Date" is the Julian Date - 2400000.5 and so 
  !       has a zero point of 17th November 1858 AD 00:00 Universal Time.
  ! 
  !       Due to the small area coverable by the GRIB output there is no need
  !       to use a "Modified Julian Date".
  !
  ! The Julian day number is stored as two doubles to guarantee sufficient
  ! precision on all machines. The complete value is (day + fraction)
  ! although this addition will sometimes lose precision. Note that day
  ! might not be an integer value and fraction might be greater than one.

  ! This is the clean version !!!!!!!

  !INTEGER, PRIVATE :: days_in_month(0:12) = (/ &
!&      0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  ! This is the version compatible with the old code ...

  INTEGER, PRIVATE :: days_in_month(0:12) = (/ &
&      365, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)



  TYPE j_date
     REAL :: day
     REAL :: fraction
  END TYPE j_date

CONTAINS

  FUNCTION SetYMD (ky, km, kd, zfraction) RESULT (julian_day)

    TYPE (j_date) :: julian_day
    
    REAL, OPTIONAL :: zfraction

    INTEGER, INTENT(IN) :: ky, km, kd

    INTEGER :: ib, iy, im

    REAL :: zf

    IF (.NOT. PRESENT(zfraction)) THEN
       julian_day%fraction = 0.0
    ELSE 
       julian_day%fraction = zfraction
    ENDIF

    IF (km <= 2) THEN
        iy = ky-1
        im = km+12
     ELSE 
        iy = ky
        im = km
    ENDIF

    IF (ky > 1582 .OR. (ky == 1582 .AND. km > 10 &
&        .OR. (km == 10 .AND. kd >= 15))) THEN

       ! 15th October 1582 AD or later

       ib = INT(iy/400)-INT(iy/100)
    ELSE 

       ! 4th October 1582 AD or earlier

       ib = -2
    ENDIF

    julian_day%day = FLOOR(365.25*iy)+INT(30.6001*(im+1))+ib+1720996.5+kd

    zf = julian_day%day-AINT(julian_day%day)+julian_day%fraction
    julian_day%day = AINT(julian_day%day)
    if (zf >= 1.0) THEN
       julian_day%fraction = zf-AINT(zf)
       julian_day%day = julian_day%day+AINT(zf)
    ELSE
       julian_day%fraction = zf-AINT(zf)
    ENDIF 

  END FUNCTION SetYMD

  SUBROUTINE YMD (julian_day, ky, km, kd, zfraction)
    
    TYPE (j_date), INTENT(IN) :: julian_day
    INTEGER, INTENT(OUT) :: km, kd, ky
    REAL, INTENT(OUT) :: zfraction
    
    REAL :: za, zb, zc, zd, ze, zf

    za = FLOOR(julian_day%day+julian_day%fraction+0.5)

    IF (za < 2299161.0) THEN
       zc = za+1524.0
    ELSE 
       zb = FLOOR((za-1867216.25)/36524.25)
       zc = za+zb-FLOOR(zb/4)+1525
    ENDIF

    zd = FLOOR((zc-122.1)/365.25)
    ze = FLOOR(365.25*zd)
    zf = FLOOR((zc-ze)/30.6001)

    kd = INT(zc-ze- FLOOR(30.6001*zf))
    km = INT(zf-1-12*FLOOR((zf+0.0001)/14))
    ky = INT(zd-4715-((7+km)/10))
    zfraction = (julian_day%day+0.5-za)+julian_day%fraction;

  END SUBROUTINE YMD

  SUBROUTINE YD (julian_day, ky, kd, zfraction)
    
    TYPE (j_date), INTENT(IN) :: julian_day
    INTEGER, INTENT(OUT) :: kd, ky
    REAL, INTENT(OUT) :: zfraction
    
    INTEGER :: i, im
    REAL :: za, zb, zc, zd, ze, zf

    za = FLOOR(julian_day%day+julian_day%fraction+0.5)

    IF (za < 2299161.0) THEN
       zc = za+1524.0
    ELSE 
       zb = FLOOR((za-1867216.25)/36524.25)
       zc = za+zb-FLOOR(zb/4)+1525
    ENDIF

    zd = FLOOR((zc-122.1)/365.25)
    ze = FLOOR(365.25*zd)
    zf = FLOOR((zc-ze)/30.6001)

    kd = INT(zc-ze- FLOOR(30.6001*zf))
    im = INT(zf-1-12*FLOOR((zf+0.0001)/14))
    ky = INT(zd-4715-((7+im)/10))
    zfraction = (julian_day%day+0.5-za)+julian_day%fraction

    DO i = 1, im
       IF (im == 2 .AND. (MOD(ky,4) == 0 .AND. MOD(ky,100) /= 0) &
&           .OR. MOD(ky,400) == 0) THEN
            kd = kd+29
            CYCLE
        ENDIF
        kd = kd + days_in_month(i-1)
    ENDDO

  END SUBROUTINE YD
   ! 
   subroutine time_difference(datespw,datesim,dtsec)
       !
       integer ijul1, ijul2     
       integer yyyy1,mm1,dd1,hh1,mn1,ss1,yyyy2,mm2,dd2,hh2,mn2,ss2
       integer dtsec,sec1,sec2
       !
       character*15  :: datespw
       character*15  :: datesim
       !
       read(datespw,'(I4,2I2,1X,3I2)')yyyy1,mm1,dd1,hh1,mn1,ss1
       read(datesim,'(I4,2I2,1X,3I2)')yyyy2,mm2,dd2,hh2,mn2,ss2
       !
       ! Compute Julian days
       !
       ijul1 = julian_date (yyyy1,mm1,dd1)
       ijul2 = julian_date (yyyy2,mm2,dd2)
       !
       sec1  = hh1*3600 + mn1*60 + ss1
       sec2  = hh2*3600 + mn2*60 + ss2
       dtsec = (ijul2 - ijul1)*86400 + sec2 - sec1
       !
   end subroutine
   !
   FUNCTION julian_date (yyyy, mm, dd) RESULT (julian)
       ! converts calendar date to Julian date
       ! cf Fliegel & Van Flandern, CACM 11(10):657, 1968
       ! example: julian_date(1970,1,1)=2440588
       INTEGER,INTENT(IN) :: yyyy,mm,dd
       INTEGER :: julian
       julian = dd-32075+1461*(yyyy+4800+(mm-14)/12)/4 + &
       367*(mm-2-((mm-14)/12)*12)/12- &
       3*((yyyy + 4900 + (mm - 14)/12)/100)/4
   END FUNCTION julian_date
   !
   function date_to_iso8601(date_string) result (date_iso8601)
       character(len=*), intent(in) :: date_string
       character*256 :: date_iso8601
       integer yyyy,mm,dd,hh,mn,ss
       ! assume an input format of yyyymmdd HHMMSS
       read(date_string,'(I4,2I2,1X,3I2)')yyyy,mm,dd,hh,mn,ss 
       ! return an iso8601 formatted date
       !                       yyyy -  mm   -  dd      HH   :  MM   :  SS
       write(date_iso8601, '(I4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)')  yyyy, '-', mm, '-', dd, ' ', hh, ':', mn, ':', ss  !in DFM/FEWS without the T
       !write(date_iso8601, '(I4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)')  yyyy, '-', mm, '-', dd, 'T', hh, ':', mn, ':', ss   
       !
   end function   
   !   
   function convert_fewsdate(timein, timent, treftimefews, trefstr) result (timeout)
       !
       implicit none
       !
       character*41 :: treftimefews
       character*15 :: trefstr
       !
       integer ijul1, ijul2, itb, timent     
       integer yyyy1,mm1,dd1,hh1,mn1,ss1,yyyy2,mm2,dd2,hh2,mn2,ss2, tmp
       integer*8 dtsec,sec1,sec2
       !
       real*4, dimension(:),     allocatable :: timein             
       real*4, dimension(:),     allocatable :: timeout       
       !       
       read(trefstr,'(I4,2I2,1X,3I2)')yyyy1,mm1,dd1,hh1,mn1,ss1
       !
       ! Fews time input: minutes since 1970-01-01 00:00:00.0 +0000       
       ! yyyy -  mm   -  dd      HH   :  MM   :  SS       
       read(treftimefews, '(A14,I4,A1,I2,A1,I2,A1,I2,A1,I2,A1,I2)')tmp,yyyy2,tmp,mm2,tmp,dd2,tmp,hh2,tmp,mn2,tmp,ss2 !assume time in UTC, whole seconds 
       !
       ijul1 = julian_date (yyyy1,mm1,dd1)
       ijul2 = julian_date (yyyy2,mm2,dd2)
       !
       sec1  = hh1*3600 + mn1*60 + ss1
       sec2  = hh2*3600 + mn2*60 + ss2
       dtsec = (ijul2 - ijul1)*86400 + sec2 - sec1       
       !       
       allocate(timeout(timent))
       !
       timeout = real((int(timein) * 60) + dtsec)  ! time from fews is in minutes, then correct for use wrt sfincs treftime
       !       
   end function      
   !   
   function time_to_string(t_sec, tref_string) result (date_string)
       !   
       character*15                 :: date_string
       character*15                 :: tref_string
       real*8                       :: t_sec
       integer yyyy,mm,dd,hh,mn,ss,ndays_add
       real zfraction,secs_time,sec_add,ref_fraction
       type(j_date) :: jd
       !
       ! assume an input format of yyyymmdd HHMMSS
       !
       read(tref_string,'(I4,2I2,1X,3I2)')yyyy,mm,dd,hh,mn,ss 
       ref_fraction = real(hh)/24 + real(mn)/1440 + real(ss)/86400
       jd = SetYMD (yyyy, mm, dd, 0.0)
       !
       sec_add = ref_fraction*86400 + t_sec
       ndays_add = int(sec_add/86400)
       secs_time = sec_add - ndays_add*86400.0       
       jd%day = jd%day + 1.0*ndays_add
       !
       call YMD(jd, yyyy, mm, dd, zfraction)
       hh = int(secs_time/3600)
       mn = int((secs_time - hh*3600)/60)
       ss = int(secs_time - hh*3600 - mn*60)
       !
       ! return a formatted date
       !
       write(date_string, '(I4,I0.2,I0.2,A1,I0.2,I0.2,I0.2)')  yyyy, mm, dd, '.', hh, mn, ss
       !
   end function   
   !
end module