MODULE DE_utils

! Defines obliquity and utc_to_tt functions
  USE rot

contains

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! PLEPH is using lun = 12 to access the ephemerides file.
!
! DO NOT USE LUN + 12 IN FILE ACCESS IN YOUR PROGRAM
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! All computations are done in Equatorial Barycentric coordinates.
! Convert to and from Ecliptic, and to and from Heliocentric as needed.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

  subroutine init_DE()

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine opens and reads in the survey description file.
! Angles are returned in radian.
! Time is returned in TDB. Initially read as UTC in survey description
! file. SO need to apply an offset.
! We don't use TCB, so no linear shift.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : January 2010
! Version 2 : February 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! OUTPUT
!     common block /de_const/
!     masses: Mass of the Sun and the planets (10*R8)
!     ss    : initial, final time and time step of ephemerides (3*R8)
!     obli  : Obliquity of Earth on JD = 2451544.5 (degree) (R8)
!     idtarg: index of Sun and planets for call to Pleph (10*I4)
!     nctr  : Index for origin of coordinates in Pleph (I4)
!     first : Logical. If .true. init_DE has not been called yet.
!                1 = SUN                     7 = JUPITER
!                2 = MERCURY                 8 = SATURN
!                3 = VENUS                   9 = URANUS 
!                4 = EARTH                  10 = NEPTUNE
!                5 = MOON                   11 = PLUTO
!                6 = MARS                   12 = EMB
!
! See PLEPH for meaning of idtarg and nctr.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    USE DE_const
    USE DE_ephhdr
    USE DE_chrhdr
    USE DE_size

    implicit none
    integer (kind=4) :: nvs, i
!    real (kind=8) :: obliquity
!    external obliquity

    if (first) then
       call const (CNAM, vals, ss, nvs)
       do i = 1, nvs
          if (CNAM(i)(1:4) .eq. 'GMS ') masses(1) = vals(i)
          if (CNAM(i)(1:4) .eq. 'GM1 ') masses(2) = vals(i)
          if (CNAM(i)(1:4) .eq. 'GM2 ') masses(3) = vals(i)
          if (CNAM(i)(1:4) .eq. 'GMB ') masses(12) = vals(i)
          if (CNAM(i)(1:5) .eq. 'EMRAT') masses(15) = vals(i)
          if (CNAM(i)(1:4) .eq. 'GM4 ') masses(6) = vals(i)
          if (CNAM(i)(1:4) .eq. 'GM5 ') masses(7) = vals(i)
          if (CNAM(i)(1:4) .eq. 'GM6 ') masses(8) = vals(i)
          if (CNAM(i)(1:4) .eq. 'GM7 ') masses(9) = vals(i)
          if (CNAM(i)(1:4) .eq. 'GM8 ') masses(10) = vals(i)
          if (CNAM(i)(1:4) .eq. 'GM9 ') masses(11) = vals(i)
          if (CNAM(i)(1:3) .eq. 'AU ') au = vals(i)
       end do
! Since masses(12) = masses(4) + masses(5)
! and   masses(15) = masses(4)/masses(5),
! one finds masses(5) = masses(12)/(1.+masses(15))
! and       masses(4) = masses(5)*masses(15)
       masses(5) = masses(12)/(1.d0+masses(15))
       masses(4) = masses(5)*masses(15)
! If we want to have J2000 ecliptic reference frame
       obli = obliquity(2451544.5d0)
       first = .false.
    end if

    return

  end subroutine init_DE
!
!++++++++++++++++++++++++
!
  subroutine FSIZER2(NRECL,KSIZE,NRFILE,NAMFIL)
!
!++++++++++++++++++++++++

!     FSIZER2 opens the file, 'NAMFIL', with an arbitrarily large record length,
!     reads the first record, and uses the information read to compute KSIZE,
!     the number of single-precision words in a record.
!
    USE DE_const
    USE DE_ephhdr
    USE DE_chrhdr

    implicit none

    integer (kind=4) :: KSIZE,NRECL,NRFILE
    integer (kind=4) :: I,J,K,ND,KHI,KMX,MRECL

    character (len=80) :: NAMFIL

!  *****************************************************************
!  *****************************************************************

!     The parameters NRECL, NRFILE, and NAMFIL are to be set by the user

!  *****************************************************************

!     NRECL=1 if "RECL" in the OPEN statement is the record length in 
!                       in single precision words
!     (for VAX/VMS, NRECL is probably 1)

!     NRECL=4 if "RECL" in the OPEN statement is the record length in bytes
!     (for UNIX, NRECL is probably 4)

    NRECL=4

!  *****************************************************************

!     NRFILE is the internal unit number used for the ephemeris file

    NRFILE=12

!  *****************************************************************

!     NAMFIL is the external name of the binary ephemeris file

    NAMFIL='binEphem.DE' 

!  *****************************************************************
!  *****************************************************************

!     Open the direct-access files and get the pointers in order to
!     Determine the size of the ephemeris record.

!     Starting with DE430, the number of ephemeris constants used
!     in the integration exceeded the old maximum number of 400
!     so the number of constants NCON is read first, to find
!     the location of all the pointers.

    MRECL = NRECL*10000

    open(NRFILE, &
         FILE=NAMFIL, &
         ACCESS='DIRECT', &
         FORM='UNFORMATTED', &
         RECL=MRECL, &
         STATUS='OLD')

    read(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON

    if (NCON .le. OLDMAX) then
       read(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT, &
            ((IPT(I,J),I=1,3),J=1,12),DENUM,(IPT(I,13),I=1,3), &
            (IPT(I,14),I=1,3), &
            (IPT(I,15),I=1,3)
    else
       if(NCON .gt. NMAX)then
          write(*,*)'Number of ephemeris constants too big for FSIZER2'
          stop
       endif
       read(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT, &
            ((IPT(I,J),I=1,3),J=1,12),DENUM,(IPT(I,13),I=1,3), &
            (CNAM(J),J=OLDMAX+1,NCON), &
            (IPT(I,14),I=1,3), &
            (IPT(I,15),I=1,3)
    endif

    close(NRFILE)

!     Find the number of data coefficients from the pointers

    KMX = 0
    KHI = 0

    do I = 1,15
       if (IPT(1,I) .ge. KMX) then
          KMX = IPT(1,I)
          KHI = I
       endif
    enddo

    ND = 3
    if (KHI .EQ. 12) ND=2
    if (KHI .EQ. 15) ND=1

    KSIZE = IPT(1,KHI) - 1 + ND*IPT(2,KHI)*IPT(3,KHI)
    KSIZE = 2*KSIZE

    return

  end subroutine FSIZER2
!
!++++++++++++++++++++++++
!
  subroutine FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)
!
!++++++++++++++++++++++++

!     FSIZER3 requires the user to set the values of KSIZE
!     as well as the values of NRECL, NRFILE, and NAMFIL.

!     KSIZE is listed in the first line of the file header.xxx
    implicit none

    integer (kind=4) :: KSIZE,NRECL,NRFILE
    character (len=80) :: NAMFIL

!  *****************************************************************
!  *****************************************************************

!     The parameters NRECL, NRFILE, NAMFIL and KSIZE are to be set by the user

!  *****************************************************************

!     NRECL=1 if "RECL" in the OPEN statement is the record length in 
!             single-precisions words
!     (for VAX/VMS, NRECL is probably 1)
!
!     NRECL=4 if "RECL" in the OPEN statement is the record length in bytes
!     (for Unix/Linux, NRECL is probably 4)

    NRECL = 4

!  *****************************************************************

!     NRFILE is the internal unit number used for the ephemeris file

    NRFILE = 12

!  *****************************************************************

!     NAMFIL is the external name of the binary ephemeris file

    NAMFIL = 'binEphem.DE'

!  *****************************************************************

!     KSIZE must be set by the user according to the ephemeris to be read

!     For  de200, set KSIZE to 1652
!     For  de405, set KSIZE to 2036
!     For  de406, set KSIZE to 1456
!     For  de414 through de429,  set KSIZE to 2036
!     For  de430 & de431, versions without TT-TDB have KSIZE = 2036
!                         versions with    TT-TDB have KSIZE = 1964

    KSIZE = 2244

!  *****************************************************************
!  *******************************************************************

    return

  end subroutine FSIZER3
!
!++++++++++++++++++++++++
!
  subroutine READHD( &
       AU_KM,DAY_SEC,IAU_AU, &
       NMAX1,NCON,CNAM,VALS,SS,IPT, &
       FACTE,FACTM,XSCALE,VSCALE, &
       NRFILE,NCOEFF)
!
!     Read the first record of a JPL binary planetary ephemeris file
!     to determine the record length, return the ephemeris constants,
!     and calculate scale factors to be used in PLEPH.
!
!$ Disclaimer
!
!     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
!     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
!     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
!     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
!     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
!     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
!     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
!     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
!     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
!     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
!
!     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
!     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
!     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
!     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
!     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
!     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!
!     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
!     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
!     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
!     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
!
    USE DE_const
    USE DE_size

    implicit none

!     Input parameters:

    logical :: AU_KM, DAY_SEC, IAU_AU
    integer (kind=4) :: NMAX1

!     Output parameters: 

    integer (kind=4) :: NCOEFF
    integer (kind=4) :: NRFILE
    integer (kind=4) :: IPT(3,15), DENUM, NCON
    real (kind=8) :: SS(3), VALS(*), AU, EMRAT
    CHARACTER (len=6) :: TTL(14,3), CNAM(*)
    real (kind=8) :: FACTE,FACTM
    real (kind=8) :: XSCALE,VSCALE

!     Local variables

    integer (kind=4) :: NRECL
    data NRECL /0/
    integer (kind=4) :: KSIZE
    character (len=80) :: NAMFIL

    real (kind=8) :: iau
    data iau /149597870.700d0/

    integer (kind=4) :: I,J,K,IRECSZ

! ************************************************************************
! ************************************************************************

!     The user must select one of the following by deleting the 'C' in column 1

! ************************************************************************

    call FSIZER2(NRECL,KSIZE,NRFILE,NAMFIL)
!    call FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)

    if(NRECL .EQ. 0) WRITE(*,*)'  ***** FSIZER IS NOT WORKING *****'

! JMP    write(*,*) 'KSIZE =', KSIZE

! ************************************************************************
! ************************************************************************

    IRECSZ = NRECL*KSIZE
    NCOEFF = KSIZE/2

    open(NRFILE, &
         FILE   = NAMFIL, &
         ACCESS ='DIRECT', &
         FORM   ='UNFORMATTED', &
         RECL   = IRECSZ, &
         STATUS ='OLD')

    read(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT

    if (NCON .le. OLDMAX) then

       read(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT, &
            ((IPT(I,J),I=1,3),J=1,12),DENUM,(IPT(I,13),I=1,3), &
            (IPT(I,14),I=1,3), &
            (IPT(I,15),I=1,3)

       read(NRFILE,REC=2)(VALS(I),I=1,OLDMAX)

    else

       if(NCON .gt. NMAX1) then
          write(*,*)'Number of ephemeris constants too big in READHD'
          stop
       endif

       read(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT, &
            ((IPT(I,J),I=1,3),J=1,12),DENUM,(IPT(I,13),I=1,3), &
            (CNAM(J),J=OLDMAX+1,NCON), &
            (IPT(I,14),I=1,3), &
            (IPT(I,15),I=1,3)

       read(NRFILE,REC=2)(VALS(I),I=1,NCON)

    endif

    FACTE = 1.D0/(1.D0+EMRAT)
    FACTM = FACTE - 1.D0

    if(FACTE .eq. 0.D0) then
       write(*,*)'Invalid value of EMRAT from file in READHD'
       stop
    endif

    if(AU_KM)then
       if(IAU_AU)then
          XSCALE = 1.d0/IAU
       else
          XSCALE = 1.d0/AU
       endif
    else
       XSCALE = 1.d0
    endif

    if(DAY_SEC) then
       VSCALE = XSCALE*86400.d0
    else
       VSCALE = XSCALE
    endif

    if(DENUM .eq. 0) stop 'DENUM not found by READHD in constants'
! JMP    write(6,*)
! JMP    write(6,'(a27,i3.3)')' JPL planetary ephemeris DE',DENUM
! JMP    write(6,*)'Requested output units are :'
! JMP    if(AU_KM)then
! JMP       if(IAU_AU)then
! JMP          write(6,*)'IAU au for distance'
! JMP          if(DAY_SEC)then
! JMP             write(6,*)'IAU au/day for velocity'
! JMP          else
! JMP             write(6,*)'IAU au/sec for velocity'
! JMP          endif
! JMP       else
! JMP          write(6,'(a2,i3.3,a16)')'DE',DENUM,' au for distance'
! JMP          if(DAY_SEC)then
! JMP             write(6,'(a2,i3.3,a20)')'DE',DENUM,' au/day for velocity'
! JMP          else
! JMP             write(6,'(a2,i3.3,a20)')'DE',DENUM,' au/sec for velocity'
! JMP          endif
! JMP       endif
! JMP    else
! JMP       write(6,*)'km for distance'
! JMP       if(DAY_SEC)then
! JMP          write(6,*)'km/day for velocity'
! JMP       else
! JMP          write(6,*)'km/sec for velocity'
! JMP       endif
! JMP    endif

    return
  end subroutine READHD
!
!++++++++++++++++++++++++
!
  subroutine INTCHB( BUF, T, LINT, NCF, NCM, NSC, REQ, PV)
!
!-----------------------------------------------------------------------
!  INTCHB (INTerpolate CHeByshev polynomial) computes the components of
!  position by interpolating Chebyshev polynomials, and if requested, it
!  also computes the components of velocity by differentiating the
!  Chebyshev polynomials and then interpolating.
!
!  Inputs:
!
!   BUF(1..NCF, 1..NCM, 1..NSC)  Chebyshev coefficients
!   T     Fractional time within the interval covered by the
!         coefficients at which the interpolation is required.
!         T must satisfy 0 .LE. T .LE. 1
!   LINT  Length of the interval in input time units
!   NCF   Number of coefficients per component
!   NCM   Number of components per set of coefficients
!   NSC   Number of sets of coefficients within interval
!   REQ   Request indicator:
!           REQ= 1 for position components only,
!           REQ= 2 for both position and velocity components.
!
!  Output:
!
!   PV(i)    Computed position and velocity components:
!            1 to NCM are position components.
!            NCM+1 to 2*NCM  are velocity components.
!-----------------------------------------------------------------------
!
!$ Disclaimer
!
!     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
!     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
!     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
!     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
!     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
!     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
!     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
!     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
!     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
!     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
!
!     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
!     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
!     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
!     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
!     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
!     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!
!     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
!     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
!     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
!     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
!
    implicit none

!     Input variables:

    integer (kind=4) :: NCF, NCM, NSC, REQ
    real (kind=8) :: BUF(NCF,NCM,NSC), T, LINT

!     Output variables

    real (kind=8) :: PV(*)

!     Local variable

    real (kind=8) :: PC(18), VC(2:18), TEMP, TC, TTC, BMA
    integer (kind=4) :: NS, L, NP, NV, NC, I, J

    data       PC(1) /1.D0/, PC(2) /2.D0/, VC(2) /1.D0/

    save

!     Compute set number within interval (L), and scaled Chebyshev time
!     within that set (TC).

    NS   = NSC
    TEMP = T * dble(NS)
    TC   = 2.d0 * ( TEMP-DINT(TEMP) ) - 1.d0
    L    = TEMP + 1
    if (L .gt. NS) then
       L  = L - 1
       TC = TC + 2.d0
    endif

!-- If the Chebyshev time has changed, compute new polynomial values.

    if (TC .ne. PC(2)) then
       NP = 3
       NV = 4
       PC(2) = TC
       TTC   = TC + TC
       VC(3) = TTC + TTC
    endif

!-- Compute position polynomial values.

    NC = NCF
    do  NP = NP,NC
       PC(NP) = TTC*PC(NP-1) - PC(NP-2)
    enddo

!-- Compute position components.

    do I = 1,NCM
       TEMP = 0.D0
       do  J = NC,1,-1
          TEMP= TEMP + PC(J)*BUF(J,I,L)
       enddo
       PV(I) = TEMP
    enddo

!-- If only positions are requested, exit

    IF (REQ .LE. 1) RETURN

!-- Compute velocity polynomial values.

    do  NV = NV,NC
       VC(NV)= TTC*VC(NV-1) + PC(NV-1) + PC(NV-1) - VC(NV-2)
    enddo

!-- Compute velocity components.

    BMA = DBLE( 2*NS ) / LINT
    do  I = 1,NCM
       TEMP= 0.D0
       do  J = NC,2,-1
          TEMP= TEMP + VC(J)*BUF(J,I,L)
       enddo
       PV(I+NCM) = TEMP*BMA
    enddo

    return

  end subroutine INTCHB
!
!++++++++++++++++++++++++
!
  subroutine PLEPH( TDB, NTARG, NCENT, PV)
!
!++++++++++++++++++++++++++
!
!     PLEPH reads the JPL planetary ephemeris file.
!     For 'bodies' (Sun, Moon, planets, Earth-Moon barycenter, 
!     solar-system-barycenter) the position and velocity of the 'body' NTARG 
!     are given with respect to the 'body' NCENT. 
!     
!     Auxiliary values are 1980 IAU nutation angles, lunar libration angles,
!     lunar Euler angle rates, and TT-TDB may also be available.
!
!     TDB is Ephemeris coordinate time, expressed in Julian days
!
!     The numbering for NTARG and NCENT is:
!
!                 1 = Mercury
!                 2 = Venus
!                 3 = Earth
!                 4 = Mars system barycenter
!                 5 = Jupiter system barycenter
!                 6 = Saturn system barycenter
!                 7 = Uranus system barycenter
!                 8 = Neptune system barycenter
!                 9 = Pluto system barycenter
!                10 = Moon (of Earth)
!                11 = Sun
!                12 = Solar-System Barycenter
!                13 = Earth-Moon barycenter
!                14 = Nutations (Longitude and Obliquity)
!                15 = Lunar Euler angles: phi, theta, psi
!                16 = Lunar angular velocity: omegax, omegay, omegaz
!                17 = TT - TDB
!
!     Note that not all ephemerides include all of the above quantities.
!     When a quantity is requested that is not on the file,
!     a warning is printed and the components of PV are set to -99.d99,
!     which is not a valid value for any quantity.
!
!      For nutations, librations, and TT-TDB,  NCENT is ignored     
!
!      PV(6)    Returned values; 
! 
!            For 'bodies', x,y,z, and vx,vy,vz are returned.
!            The values stored on the files are in units of km and km/s.
!            By default, the positions are scaled the value of the astronomical
!            unit used in the ephemeris integration to return units of 
!            au and au/day.
!            The user can override the default units by calling PLUNITS 
!            as described below.
!
!            For nutation angles, the returned values in PV are 
!                 [1-longitude, 2-obliquity, 3-long. rate, 4-oblq. rate]
!                 in units of radians and radians/day.
!            
!            For libration angles, the returned values in PV are 
!                 [1-omegax, 2-omega-y, 3=omega-z, 4-6 rates]
!                 in units of radians/day and radians/day**2.
!            
!            For TDB-TT, PV(1) is TT-TDB in seconds, 
!                        PV(2) rate of change of TT-TDB in sec/day
!
!     In addition to the main entry PLEPH, there are optional entry points.
!
!
!     entry DPLEPH( TDB2, NTARG, NCENT, PV)
!
!     Inputs TDB2, NTARG, NCENT and output PV are as described above for PLEPH,
!     except that TDB2 is a two-element array allowing the possibility
!     of specifying the TDB coordinate time with greater numerical precision.
!     Typically TDB2(1) is set to the integer number of the Julian day,
!     and TDB(2) is set to the fraction of the Julian day.
!
!     Note that with ephemeris data in 32-day long blocks, the 
!     precision of the ephemeris lookup time is limited to approximately
!     ~1.e-15*(32 days) = ~3 nanoseconds
!     since the Chebyshev coefficients are interpolated using a
!     double precision representation of a fraction of the block time span.
!
!
!     entry CONST(CNAM,VALS,SS,NCON)
!
!     No input arguments, output quantities are:
!
!        CNAM       6-character names of ephemeris constants
!        VALS       double precision values of ephemeris constants
!        SS         earliest TDB time on ephemeris (Julian days),
!                   latest TDB time on ephemeris   (Julian days), and
!                   number of days covered by each ephemeris record
!        NCON       number of ephemeris constants
!
!
!     entry PL_UNITS(AU_KM,DAY_SEC,IAU_AU)
!
!     When called before first call to PLEPH or DPLEPH or CONST, this entry
!     allows overriding the default units returned for positions and velocities.
!
!     AU_KM        True  means output length unit is astronomical units
!                  False means output length unit is km
!
!     DAY_SEC      True means velocities are returned in length/day
!                  False means velocities are returned in length/second
!
!     IAU_AU       True means use value of astronomical unit adopted in 2012
!                       (149597870.700 km) if AU_KM is True
!                  False means use value of astronomical units used at time
!                        of ephemeris integration.
!                  (Note that the value of AU returned in the ephemeris
!                   constants will not be changed by this option).
!
!     Defaults correspond to PL_UNITS(.TRUE. , .TRUE. , .FALSE.)
!
!     Last updated 14 March 2015
!
!$ Disclaimer
!
!     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
!     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
!     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
!     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
!     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
!     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
!     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
!     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
!     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
!     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
!
!     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
!     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
!     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
!     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
!     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
!     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!
!     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
!     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
!     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
!     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
!     
!++++++++++++++++++++++++++
!
    USE DE_const
    USE DE_size
    USE DE_ephhdr
    USE DE_chrhdr

    implicit none

!++++++++++++++++++++++++++

!     arguments for main call PLEPH

!++++++++++++++++++++++++++

    real (kind=8) :: TDB
    integer (kind=4) :: NTARG,NCENT

    real (kind=8) :: PV(6)

!     argument for DPLEPH

    real (kind=8) :: TDB2(2)

!++++++++++++++++++++++++++

!     arguments for return CONST

!++++++++++++++++++++++++++

    character (len=6) :: NAMS1(*)
    real (kind=8) :: VALS1(*), SS1(*)
    integer (kind=4) :: NCON1

!++++++++++++++++++++++++++

!     arguments for PL_UNITS (overwrite default values in PLEPH)
      
!++++++++++++++++++++++++++

    logical :: AU_KM1, DAY_SEC1, IAU_AU1

!++++++++++++++++++++++++++

!     internal variables

!++++++++++++++++++++++++++

    real (kind=8) :: data(NMAX)
    real (kind=8) :: PV1(6),PV2(6),PVM(6),PVB(6)

    real (kind=8) :: T1,T2,TF,TF1,TF2
    real (kind=8) :: FACTE,FACTM,XSCALE,VSCALE,SUMT
    integer (kind=4) :: NRFILE,NCOEFF

    real (kind=8) :: SECSPAN
    real (kind=8) :: SECDAY
    data SECDAY /86400.D0/

    integer (kind=4) :: IPVA
    data IPVA /2/

    integer (kind=4) :: KODE(15)
    data KODE        /1,2,13,4,5,6,7,8,9,13,11,0,3,3,15/

    logical :: AU_KM
    data AU_KM /.true./
    logical :: DAY_SEC
    data DAY_SEC /.true./
    logical :: IAU_AU
    data IAU_AU/.false./

!    logical :: FIRST    /.true./
    logical :: RDCONST
    data RDCONST /.false./

    integer (kind=4) :: I, j
    integer (kind=4) :: NRREC
    data NRREC / 0/

    save

!++++++++++++++++++++++++++

!     Main entry point PLEPH

    T1 = TDB
    T2 = 0.D0

    go to 1

!++++++++++++++++++++++++++

    entry DPLEPH( TDB2, NTARG, NCENT, PV)

!     Typically TDB2(1) is set to the integer number of the Julian day,
!     and TDB(2) is set to the fraction of the Julian day.

    T1 = TDB2(1)
    T2 = TDB2(2)

    go to 1

!++++++++++++++++++++++++++

    entry CONST(NAMS1,VALS1,SS1,NCON1)
!
!     This entry allows one to retrieve all of the constants associated with
!     the ephemeris file. The values returned are copies of the values
!     read off the ephemeris file and saved within PLEPH.
!
!     There is no INPUT
!
!     OUTPUT:
!
!     NCON           [integ.] : number of ephemeris constants
!     CNAM           [char*6] : names of NCON ephemeris constants
!     VALS           [d.p.]   : values of NCON ephemeris constants
!     SS(3)          [d.p.]   : SS(1) is starting JED of the ephemeris file
!                               SS(2) is  ending JED of the ephemeris file
!                               SS(3) is the length(in days) of a file record
!
    RDCONST = .true.

!++++++++++++++++++++++++++

1   continue

    if (FIRST) then
       call READHD( &
            AU_KM,DAY_SEC,IAU_AU, &
            NMAX,NCON,CNAM,VALS,SS,IPT, &
            FACTE,FACTM,XSCALE,VSCALE, &
            NRFILE,NCOEFF)
       if(NCOEFF .gt. NMAX)stop 'Coefficient array too small in PLEPH' 
       FIRST = .false.
       SECSPAN = SECDAY * SS(3)
       do I = 1,NMAX
          DATA(I) = -999999999999999.d0
       enddo
    endif

    if(RDCONST) then

       NCON1 = NCON
       do I = 1,NCON
          NAMS1(I) = CNAM(I)
          VALS1(I) = VALS(I)
       enddo

       do I = 1,3
          SS1(I) = SS(I)
       enddo
       RDCONST = .false.
       return
    endif

!-----------------------------------------------------------------------
!     start main ephemeris lookup operation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     first check to see if NTARG and NCENT in allowed ranges
!-----------------------------------------------------------------------

    if (NTARG .le. 0)  stop 'invalid NTARG < 0 in PLEPH'
    if (NTARG .gt. 17) stop 'invalid NTARG >17 in PLEPH'
    if (NCENT .le. 0)then
       if(NTARG .lt. 14) stop 'invalid NCENT < 0 in PLEPH'
    endif
    if (NCENT .gt. 13)then
       if(NTARG .lt. 14) stop 'invalid NCENT > 13 in PLEPH'
    endif

    if(NTARG .le. 13 .and. NCENT .ge. 14)then
       stop 'invalid NCENT >13 in for body PLEPH'
    endif

    if(NTARG .le. 13 .and. NTARG .eq. NCENT)then
       do i=1,6
          PV(I) = 0.d0
       enddo
       return
    endif

!-----------------------------------------------------------------------
!     determine which data record needed; read in if not already in DATA array
!-----------------------------------------------------------------------

    SUMT = T1+T2

    if (SUMT .lt. SS(1)) stop 'input time before file start in PLEPH'
    if (SUMT .gt. SS(2)) stop 'input time after  file start in PLEPH'

    I = INT( (SUMT-SS(1))/SS(3) ) + 3
    if(SUMT .eq. SS(2)) I = I - 1

    if (I .ne. NRREC)then
       NRREC = I
       read (NRFILE, rec = NRREC, ERR = 99)(DATA(I),I=1,NCOEFF)
    endif

    TF1 = dble(int(T1))
    TF2 = T1-TF1
    TF  = (TF1-DATA(1))/SS(3)
    TF  = TF + ((TF2+T2)/SS(3))

!-----------------------------------------------------------------------
!     do lookups for nutations, Euler angles, omega, or TT-TDB
!     which do not require a center, and might not be on file
!-----------------------------------------------------------------------

!     nutation

    if(NTARG .eq. 14) then
       if(IPT(1,12) .gt. 0  .and. IPT(2,12)*IPT(3,12) .ne. 0)then
          call INTCHB( DATA(IPT(1,12)), TF, SECSPAN, &
               IPT(2,12), 2, IPT(3,12), IPVA, PV)
          PV(3) = SECDAY*PV(3)   ! always convert radian/second to radian/day
          PV(4) = SECDAY*PV(4)
          PV(5) = -99.D99
          PV(6) = -99.D99
          return
       else
          write(6,*)'Nutations requested but not found on file'
          do I = 1,6
             PV(I) = -99.d99
          enddo
       endif
       return
    endif

!     Lunar mantle Euler angles

    if(NTARG .eq. 15) then
       if(IPT(1,13) .gt. 0  .and. IPT(2,13)*IPT(3,13) .ne. 0)then
          call INTCHB( DATA(IPT(1,13)), TF, SECSPAN, &
               IPT(2,13), 3, IPT(3,13), IPVA, PV)
          PV(4) = SECDAY*PV(4)   ! always convert radian/second to radian/day
          PV(5) = SECDAY*PV(5)
          PV(6) = SECDAY*PV(6)
          return
       else
          write(6,*)'Mantle Euler angles requested but not on file'
          do I = 1,6
             PV(I) = -99.d99
          enddo
       endif
       return
    endif

!     Lunar mantle angular velocity

    if(NTARG .eq. 16) then
       if(IPT(1,14) .gt. 0  .and. IPT(2,14)*IPT(3,14) .ne. 0)then
          call INTCHB( DATA(IPT(1,14)), TF, SECSPAN, &
               IPT(2,14), 3, IPT(3,14), IPVA, PV)
          PV(4) = SECDAY*PV(4)   ! always convert to radian/day**2
          PV(5) = SECDAY*PV(5)
          PV(6) = SECDAY*PV(6)
          return
       else
          write(6,*)'Mantle angular velocity requested but not on file'
          do I = 1,6
             PV(I) = -99.d99
          enddo
       endif
       return
    endif

!     TT-TDB

    if(NTARG .eq. 17) then
       if(IPT(1,15) .gt. 0  .and. IPT(2,15)*IPT(3,15) .ne. 0)then
          call INTCHB( DATA(IPT(1,15)), TF, SECSPAN, &
               IPT(2,15), 1, IPT(3,15), IPVA, PV)
          PV(2) = SECDAY*PV(2)   ! always convert second/second to second/day
          do I = 3,6
             PV(I) = -99.d99
          enddo
          return
       else
          write(6,*)'TT-TDB requested but not found on file'
          do I = 1,6
             PV(I) = -99.d99
          enddo
       endif
       return
    endif

!-----------------------------------------------------------------------
!     do lookups for bodies
!-----------------------------------------------------------------------

    if(NTARG .eq. 10 .and. NCENT .eq. 3)then
       call INTCHB( DATA(IPT(1,10)), TF, SECSPAN, &
            IPT(2,10), 3, IPT(3,10), IPVA, PVM)
       do I = 1,3
          PV(i)   = PVM(i)*XSCALE
          PV(i+3) = PVM(i+3)*VSCALE
       enddo
       return
    endif

    if(NTARG .eq. 3 .and. NCENT .eq. 10)then
       call INTCHB( DATA(IPT(1,10)), TF, SECSPAN, &
            IPT(2,10), 3, IPT(3,10), IPVA, PVM)
       do I = 1,3
          PV(i)   = -PVM(i)*XSCALE
          PV(i+3) = -PVM(i+3)*VSCALE
       enddo
       return
    endif

    do I = 1,6
       PV1(I) = 0.d0
       PV2(I) = 0.d0
    enddo

    if(NTARG .eq. 3 .or. NTARG .eq. 10 .or. &
         NCENT .eq. 3 .or. NCENT .eq. 10)then

       call INTCHB( DATA(IPT(1,10)), &
            TF, SECSPAN, IPT(2,10), 3, IPT(3,10), IPVA, PVM)
       call INTCHB( DATA(IPT(1,3)), &
            TF, SECSPAN, IPT(2,3),  3, IPT(3,3),  IPVA, PVB)

       if(NTARG .eq. 3 .or. NTARG .eq. 10)then

          if(NCENT .lt. 12)then
             call INTCHB( DATA(IPT(1,NCENT)),TF, SECSPAN, &
                  IPT(2,NCENT), 3, IPT(3,NCENT), IPVA, PV2)
          else if(NCENT .eq. 13)then
             call INTCHB( DATA(IPT(1,3)),TF, SECSPAN, &
                  IPT(2,3), 3, IPT(3,3), IPVA, PV2)
          endif

          if(NTARG .eq. 3)then
             do I = 1,6
                PV(I) = PVB(I) - FACTE*PVM(I) - PV2(I)
             enddo
          else
             do I=1,6
                PV(I) = PVB(I) - FACTM*PVM(I) - PV2(I)
             enddo
          endif

       else 

          if(NTARG .lt. 12) then
             call INTCHB( DATA(IPT(1,NTARG)),TF, SECSPAN, &
                  IPT(2,NTARG), 3, IPT(3,NTARG), IPVA, PV1)
          else if (NTARG .eq. 13)then
             call INTCHB( DATA(IPT(1,3)),TF, SECSPAN, &
                  IPT(2,3), 3, IPT(3,3), IPVA, PV1)
          endif

          if(NCENT .eq. 3)then
             do I = 1,6
                PV(I) = PV1(I) - (PVB(I) - FACTE*PVM(I))
             enddo
          else
             do I=1,6
                PV(I) = PV1(I) - (PVB(I) - FACTM*PVM(I))
             enddo
          endif
       endif

       do I = 1,3
          PV(i)   = PV(i)*XSCALE
          PV(i+3) = PV(i+3)*VSCALE
       enddo
       return
    endif

! -- Only cases left involve neither Earth nor Moon

    if(NTARG .lt. 12) then
       call INTCHB( DATA(IPT(1,NTARG)),TF, SECSPAN, &
            IPT(2,NTARG), 3, IPT(3,NTARG), IPVA, PV1)
    else if (NTARG .eq.13)then
       call INTCHB( DATA(IPT(1,3)),TF, SECSPAN, &
            IPT(2,3), 3, IPT(3,3), IPVA, PV1)
    endif
        
    if(NCENT .lt. 12)then
       call INTCHB( DATA(IPT(1,NCENT)),TF, SECSPAN, &
            IPT(2,NCENT), 3, IPT(3,NCENT), IPVA, PV2)
    else if (NCENT .eq. 13)then
       call INTCHB( DATA(IPT(1,3)),TF, SECSPAN, &
            IPT(2,3), 3, IPT(3,3), IPVA, PV2)
    endif

    do I = 1,3
       PV(i)   = (PV1(I) - PV2(I) )     * XSCALE
       PV(i+3) = (PV1(I+3) - PV2(I+3) ) * VSCALE
    enddo
    return

!-----------------------------------------------------------------------

    entry PL_UNITS(AU_KM1, DAY_SEC1, IAU_AU1)

    AU_KM   = AU_KM1
    DAY_SEC = DAY_SEC1
    IAU_AU  = IAU_AU1

    return

!-----------------------------------------------------------------------

99  continue

    write(6,*)'Error reading ephemeris data record in PLEPH'
    write(6,*)'T1,T2 = ',T1,T2
    write(6,*)'Record not found = ',NRREC
    stop

  end subroutine PLEPH

end MODULE DE_utils
