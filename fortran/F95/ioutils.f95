module ioutils

  use datadec

contains

  subroutine  hms(str,val)
!
!...Crack String And Create Value
!
!f2py intent(in) str
!f2py intent(out) val
    implicit none
    character(*), intent(in) :: str
    real (kind=8), intent(out) :: val
    character(1) :: c
    real (kind=8) :: piece(3), dp, sgn, z
    integer :: nstr, i, j, dpfind
!
!...Initialization
!
100 val = 0.0D00
    piece = 0.0d0
    j = 1
    dpfind = 0
    sgn = 1.0D00
    nstr = LEN(str)
    IF (nstr.le.0) RETURN
!
!...Loop Over The String
!
    DO i=1,nstr
       c = str(i:i)
!
!...Parse
!
       IF ((c.eq.'-').or.(c.eq.'e').or.(c.eq.'E') &
            .or.(c.eq.'s').or.(c.eq.'S')) THEN
          sgn = -1.0D00
       ELSEIF ((c.eq.'+').or.(c.eq.'w').or.(c.eq.'W') &
            .or.(c.eq.'n').or.(c.eq.'N')) THEN
          sgn = 1.0D00
       ELSEIF ((c.eq.':').or.(c.eq.',').or.(c.eq.' ')) THEN
          j = j+1
          dpfind = 0
          IF (j.gt.3) GO TO 110
       ELSEIF (c.eq.'.') THEN
          dpfind = 1
          dp = 1.0D00
       ELSEIF ((c.ge.'0').and.(c.le.'9')) THEN
          z = ICHAR(c)-ICHAR('0')
          IF (dpfind.eq.0) THEN
             piece(j) = 10.0D00*piece(j) + z
          ELSE
             dp = 0.1D00*dp
             piece(j) = piece(j) + dp*z
          ENDIF
       ENDIF
    ENDDO
!
!...Return
!
110 val = piece(1) + piece(2)/60.0D00 + piece(3)/3600.0D00
    val = val*sgn
    RETURN
  END subroutine hms

  subroutine trim (base_name, i1, i2, finished, len)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine trims the input character string 'base_name' and returns
! the 
!
! Author: Jean-Marc Petit
! version 1: December 2019
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     base_name: Input string with file name in it (CH)
!     len   : Length of input string (I4)
!
! OUPUT
!     i1    : Index of first character of file name (I4)
!     i2    : Index of last character of file name (I4)
!     finished: Whether it reached the end of the string when searching
!               for first character
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) base_name
!f2py intent(out) i1
!f2py intent(out) i2
!f2py intent(out) finished
!f2py intent(in) len
    implicit none

    integer, intent(in) :: len
    integer, intent(out) :: i1, i2
    character(*), intent(in) :: base_name
    logical, intent(out) :: finished

    finished = .false.
    i1 = 1
100 continue
    if ((base_name(i1:i1) .eq. char(0)) &
         .or. (base_name(i1:i1) .eq. char(9)) &
         .or. (base_name(i1:i1) .eq. ' ')) then
       i1 = i1 + 1
       if (i1 .eq. len) then
          finished = .true.
          return
       end if
       goto 100
    end if
101 continue

    i2 = len
110 continue
    if ((base_name(i2:i2) .eq. char(0)) &
         .or. (base_name(i2:i2) .eq. char(9)) &
         .or. (base_name(i2:i2) .eq. ' ')) then
       if (i2 .eq. i1) goto 111
       i2 = i2 - 1
       goto 110
    end if
111 continue

    return
  end subroutine trim

  subroutine read_file_name (base_name, i1, i2, finished, len)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine trims the input character string 'base_name' and returns
! the 
!
! Author: Jean-Marc Petit
! version 1: September 2017
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     base_name: Input string with file name in it (CH)
!     len   : Length of input string (I4)
!
! OUPUT
!     i1    : Index of first character of file name (I4)
!     i2    : Index of last character of file name (I4)
!     finished: Whether it reached the end of the string when searching
!               for first character
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) base_name
!f2py intent(out) i1
!f2py intent(out) i2
!f2py intent(out) finished
!f2py intent(in) len
    implicit none

    integer, intent(in) :: len
    integer, intent(out) :: i1, i2
    character(*), intent(in) :: base_name
    logical,intent(out) :: finished

    finished = .false.
    i1 = 1
100 continue
    if ((base_name(i1:i1) .eq. char(0)) &
         .or. (base_name(i1:i1) .eq. char(9)) &
         .or. (base_name(i1:i1) .eq. ' ')) then
       i1 = i1 + 1
       if (i1 .eq. len) then
          finished = .true.
          return
       end if
       goto 100
    end if
101 continue

    i2 = i1 + 1
110 continue
    if ((base_name(i2:i2) .ne. char(0)) &
         .and. (base_name(i2:i2) .ne. char(9)) &
         .and. (base_name(i2:i2) .ne. ' ')) then
       if (i2 .eq. len) goto 111
       i2 = i2 + 1
       goto 110
    end if
    i2 = i2 - 1
111 continue

    return
  end subroutine read_file_name

! \subroutine{parse}

! Parses a line returns a list of words by getting rid of space characters

  subroutine parse (comd, nwmax, nw, word, lw)
! \subsection{Arguments}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|comd| = command line to parse \\
! \verb|lw()| = word lengthes \\
! \verb|nw| = number of words in command line \\
! \verb|nwmax| = maximum allowed number of words \\
! \verb|word()| = words in command line
! \end{verse}
!
!f2py intent(in) comd
!f2py intent(in) nwmax
!f2py intent(out) nw
!f2py intent(out) word
!f2py intent(out) lw

! \subsubsection{Declarations}
    implicit none

    integer, intent(in) :: nwmax
    integer, intent(out) :: lw(:), nw
    character(*), intent(in) :: comd
    character(*), intent(out) :: word(:)

! \subsection{Variables}
! \subsubsection{Internal variables}
! \begin{verse}
! \verb|k| = dummy index \\
! \verb|lc| = length of command line \\
! \end{verse}

! \subsubsection{Intrinsic Fortran functions used}
! \begin{verse}
! \verb|index|
! \end{verse}

! \subsubsection{Declarations}

    integer :: k, lc, lw0
    character(1024) :: command

! \subsection{Parsing}

    do nw = 1, nwmax
       lw(nw) = 0
    end do
    lc = len(comd)
    command(1:lc) = comd(1:lc)
1000 continue
    if ((command(lc:lc) .eq. char(0)) &
         .or. (command(lc:lc) .eq. char(9)) &
         .or. (command(lc:lc) .eq. ' ')) then
       lc = lc - 1
       if (lc .eq. 0) goto 1001
       goto 1000
    end if
1001 continue
    nw = 0
    do k = 1, nwmax
       word(k) = ' '
    end do

1100 continue
    if (lc .gt. 0) then

! Gets rid of leading space characters
       if (nw .ge. nwmax) then
          write (6, *) command(1:lc)
          write (6, *) 'parse: too many words in command line.'
          stop
       end if
1050   continue
       if (command(1:1) .eq. ' ') then
          command = command (2:lc)
          lc = lc - 1
          goto 1050
       end if

! Finds a word

       nw = nw + 1
       lw0 = index(command(1:lc), ' ') - 1
       if (lw0 .le. 0) then
          lw(nw) = lc
          word(nw) = command(1:lc)
          lc = -1
       else
          word (nw) = command(1:lw0)
          lw(nw) = lw0
          command = command(lw0+2:lc)
          lc = lc - lw0 - 1
       end if
       goto 1100
    end if

    return
  end subroutine parse

  subroutine Format (angle, incode, outcod, string, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine formats an angle (in rd) into deg, min, sec or hour, min,
! sec. Output is a string.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : September 2003
! Version 2 : March 2004
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     angle : angle to format, radians (r8)
!     incode: 1 input in radian; 0 input in decimal degrees (i)
!     outcod: 1 converts to hours,min,sec; 0 converts to deg.,min,sec (i)
!
! OUTPUT
!     string: Output string (CH)
!     ierr  : Error code
!                0 : nominal run
!               10 : input data code
!               20 : wrong conversion code
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directive to create a Python module
!
!f2py intent(in) angle
!f2py intent(in) incode
!f2py intent(in) outcod
!f2py intent(out) string
!f2py intent(out) ierr
    implicit none

    integer, intent(in) :: incode, outcod
    integer, intent(out) :: ierr
    real (kind=8), intent(in) :: angle
    character(13), intent(out) :: string
    real (kind=8), parameter :: raddeg = 180.0d0/Pi
    integer :: deg, mn, si
    real (kind=8) :: sec, rm, w

    ierr = 0
    if (incode .eq. 1) then
       w = angle*raddeg
    else if (incode .eq. 0) then
       w = angle
    else
       ierr = 10
       return
    end if
    if (outcod .eq. 1) then
       w = w/15.d0
    else if (outcod .ne. 0) then
       ierr = 20
       return
    end if
    if (w .lt. 0.d0) then
       si = -1
       w = abs(w)
    else
       si = 1
    end if
    deg = int(w)
    rm = (w - deg)*60.d0
    mn = int(rm)
    sec = (rm - mn)*60.d0
    write (string, '(i4.2, 1x, i2.2, 1x, f5.2)') deg, mn, sec
    if (string(9:9) .eq. ' ') string(9:9) = '0'
    if (string(10:10) .eq. ' ') string(10:10) = '0'
    if (si .eq. -1) then
       if (deg .ge. 100) then
          string(1:1) = '-'
       else
          string(2:2) = '-'
       end if
    end if
    deg = si*deg

    return
  end subroutine Format

  character(100) function strip_comment(str)
    character (*), intent(in) :: str
    integer c_idx, i1, i2
    logical finished
    strip_comment = str
    c_idx = index(strip_comment, '!')
    if (c_idx /= 0) then
       strip_comment = strip_comment(1:c_idx-1)
    end if
    call trim(strip_comment, i1, i2, finished, len(strip_comment))
    strip_comment = strip_comment(i1:i2)

  end function strip_comment

  subroutine read_observatories(filename)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine reads in a file describing the observatories. Format is
! the one derived from MPC, as used by BK2000:
! <code (I)> <longitude (deg) (CH)> <latitude (deg) (CH)> ...
!  ... <altitude (m) (R)> <name (CH)>
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     filename: Name of file with observatory list (CH)
!
! OUTPUT
!     sitelist: List of topographique informations for sites (n*t_site)
!                 %lon: longitude (hour)
!                 %lat: latitude (rd)
!                 %altitude: altitude (m)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) filename
    implicit none

    character (len=*), intent(in) :: filename
    character (len=80) :: word(nw_max)
    character (len=256) :: line
    integer (kind=4) :: lun_in, nw, lw(nw_max), i, j
    data lun_in /15/

    open(unit=lun_in, file=filename, status='old', err=1000)
    nsites = 0

150 continue
    read (lun_in, '(a)', err=150, end=3000) line
    if (line(1:1) .eq. '#') goto 150
    call parse(line, nw_max, nw, word, lw)
    nsites = nsites + 1
    read(word(1), *, err=2000) sitelist(nsites)%code
    call hms(word(2), sitelist(nsites)%lon)
    sitelist(nsites)%lon = sitelist(nsites)%lon/15.0d0
    call hms(word(3), sitelist(nsites)%lat)
    sitelist(nsites)%lat = sitelist(nsites)%lat*drad
    read(word(4), *, err=2000) sitelist(nsites)%altitude
    j = 1
    do i = 5, nw
       sitelist(nsites)%name(j:j-1+lw(i)) = word(i)(:lw(i))
       sitelist(nsites)%name(j+lw(i):j+lw(i)) = ' '
       j = j + lw(i) + 1
    end do
    goto 150

1000 continue
    return

2000 continue
    nsites = 0

3000 continue
    close(lun_in)
    return

  end subroutine read_observatories

  subroutine read_jpl_csv(iunit, jd, pos, vel, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! lookup line in ephemeris file (pointed to by iunit) with data close to
! jd and then use velocity to adjust the pos values to that requested.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! JJ Kavelaars National Research Council of Canada
! Version 1 : November 2022
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     iunit : Logical unit number of CSV file (I4)
!     jd    : Epoch of observation (Julian day) (R8)
!
! OUTPUT
!     pos   : Position of spacecraft at epoch JD (AU) (3*R8)
!     vel   : Velocity of spacecraft at epoch JD (AU/day) (3*R8)
!     ierr  : Return code
!                0: OK
!               10: Error
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    implicit none

    integer, intent(in) :: iunit
    real(kind=8), intent(in) :: jd
    type(t_v3d), intent(out) :: pos, vel
    integer, intent(out) :: ierr

    character(len = 512) :: line
    character(len = 30) :: date
    real(kind=8) ejd
    integer :: iend, header_offset, ferr, offset
    logical :: read_header

! only read the header the first time we are called
    data read_header /.true./
    data header_offset /0/
    save read_header, header_offset

    ierr = 0
    if (read_header) then
       do
          read(iunit, '(A512)', end=999) line
          iend = len_trim(line)
          if ( line == '$$SOE' ) then
             read_header = .false.
             header_offset = FTELL(iunit)
             exit
          end if
       end do
    end if

! Loop through the ephemeris lines to get to the desired JD
! starting from line after the header
    offset = header_offset - FTELL(iunit) 
    CALL FSEEK(iunit, offset, 1, ferr)
    do
       read(iunit, '(A512)', end=999) line
       iend = len_trim(line)
       if ( line == '$$EOE' ) then
          write (0,*) "Date ", jd," outside range of JPL Ephemeris provided."
          ierr = 10
          return 
       end if
       read(line(1:iend), '(F17.8,2X,A30,6(2X,F22.16))', err=998) ejd,  date, &
            pos%x, pos%y, pos%z, vel%x, vel%y, vel%z
       if ( ejd > jd) then
          pos%x = pos%x + vel%x*(jd-ejd)
          pos%y = pos%y + vel%y*(jd-ejd)
          pos%z = pos%z + vel%z*(jd-ejd)
          return
       end if
    end do

998 continue
999 continue
    ierr = 10

    return

  end subroutine read_jpl_csv

end module ioutils

