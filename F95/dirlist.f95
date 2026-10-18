module dirlist
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Directory helpers backed by list_dirs.c (POSIX).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
  use, intrinsic :: iso_c_binding
  use parameters
  implicit none

  private
  public :: list_subdirs, file_exists

  interface
     function c_list_subdirs(root, names_buf, names_cap, name_wid) &
          bind(C, name='ossssim_list_subdirs')
       import :: c_char, c_int
       character(kind=c_char), intent(in) :: root(*)
       character(kind=c_char), intent(out) :: names_buf(*)
       integer(c_int), value, intent(in) :: names_cap, name_wid
       integer(c_int) :: c_list_subdirs
     end function c_list_subdirs

     function c_file_exists(path) bind(C, name='ossssim_file_exists')
       import :: c_char, c_int
       character(kind=c_char), intent(in) :: path(*)
       integer(c_int) :: c_file_exists
     end function c_file_exists
  end interface

contains

  subroutine f_to_c_string(fstr, cstr)
    character(*), intent(in) :: fstr
    character(kind=c_char), intent(out) :: cstr(*)
    integer :: i, n
    n = len_trim(fstr)
    do i = 1, n
       cstr(i) = fstr(i:i)
    end do
    cstr(n+1) = c_null_char
  end subroutine f_to_c_string

  logical function file_exists(path)
    character(*), intent(in) :: path
    character(kind=c_char) :: cpath(path_len+1)
    call f_to_c_string(path(1:len_trim(path)), cpath)
    file_exists = (c_file_exists(cpath) /= 0)
  end function file_exists

  subroutine list_subdirs(root, names, n_out, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! List immediate subdirectory basenames of root (exclude . and ..).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    character(*), intent(in) :: root
    character(*), intent(out) :: names(:)
    integer, intent(out) :: n_out, ierr
    character(kind=c_char) :: cpath(path_len+1)
    character(kind=c_char), allocatable :: buf(:)
    integer :: i, j, wid, cap, rc

    n_out = 0
    ierr = 0
    do i = 1, size(names)
       names(i) = ' '
    end do

    wid = len(names(1))
    cap = size(names)
    allocate (buf(cap * wid))
    buf = c_char_' '

    call f_to_c_string(root(1:len_trim(root)), cpath)
    rc = c_list_subdirs(cpath, buf, int(cap, c_int), int(wid, c_int))
    if (rc < 0) then
       write (6, *) 'list_subdirs: failed for ', root(1:len_trim(root))
       ierr = 20
       deallocate (buf)
       return
    end if

    n_out = rc
    do i = 1, n_out
       do j = 1, wid
          names(i)(j:j) = buf((i-1)*wid + j)
       end do
    end do
    deallocate (buf)
    return
  end subroutine list_subdirs

end module dirlist
