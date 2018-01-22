module ioutils

  use datadec

contains

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

    integer :: i1, i2, len
    character :: base_name*(*)
    logical :: finished

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

  subroutine parse (command, nwmax, nw, word, lw)

! \subsection{Arguments}
! \subsubsection{Definitions}
! \begin{verse}
! \verb|command| = command line to parse \\
! \verb|lw()| = word lengthes \\
! \verb|nw| = number of words in command line \\
! \verb|nwmax| = maximum allowed number of words \\
! \verb|word()| = words in command line
! \end{verse}
!
!f2py intent(in) command
!f2py intent(in) nwmax
!f2py intent(out) nw
!f2py intent(out) word
!f2py intent(out) lw

! \subsubsection{Declarations}
    implicit none

    integer :: lw(*), nw, nwmax
    character :: command*(*), word(1)*(*)

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

! \subsection{Parsing}

    do nw = 1, nwmax
       lw(nw) = 0
    end do
    lc = len(command)
1000 continue
    if ((command(lc:lc) .eq. char(0)) .or. (command(lc:lc) .eq. ' ')) then
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
          write (6, *) command
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
       lw0 = index(command, ' ') - 1
       if (lw0 .le. 0) then
          lw(nw) = lc
          word(nw) = command(1:lc)
          lc = -1
       else
          word (nw) = command (1:lw0)
          lw(nw) = lw0
          command = command (lw0+2:lc)
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

    real (kind=8), parameter :: Pi = 3.141592653589793d0, raddeg = 180.0d0/Pi
    integer :: deg, mn, ierr, incode, outcod, si
    real (kind=8) :: angle, sec, rm, w
    character :: string*13

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

end module ioutils
