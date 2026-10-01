! Created by  on 2025-11-11.

module debug
  IMPLICIT NONE
  INTEGER, PARAMETER :: STRLEN = 256

  ! Runtime state stays private so f90wrap wrappers cannot clash with
  ! Detos1 / Driver dummy argument names.
  LOGICAL, PRIVATE :: dbg_active = .FALSE.
  INTEGER, PRIVATE :: debug_lvl  = 0
  INTEGER, PRIVATE :: log_unit   = 6       ! default to stdout (*)
  LOGICAL, PRIVATE :: log_opened = .FALSE.

  PUBLIC :: debug_set, debug_init_from_args, debug_set_logfile
  PUBLIC :: dbg_enabled, dbg_print, assert_ok, debug_is_active

CONTAINS

  SUBROUTINE debug_set(on, level)
    LOGICAL, INTENT(IN) :: on
    INTEGER, INTENT(IN), OPTIONAL :: level
    dbg_active = on
    IF (PRESENT(level)) debug_lvl = level
  END SUBROUTINE debug_set

  LOGICAL FUNCTION debug_is_active()
    debug_is_active = dbg_active
  END FUNCTION debug_is_active

  LOGICAL FUNCTION dbg_enabled(level)
    INTEGER, INTENT(IN) :: level
    dbg_enabled = (dbg_active .AND. level <= debug_lvl)
  END FUNCTION dbg_enabled

  SUBROUTINE dbg_print(level, msg)
    INTEGER,     INTENT(IN) :: level
    CHARACTER(*),INTENT(IN) :: msg
    CHARACTER(len=8) :: date
    CHARACTER(len=10) :: time
    CHARACTER(len=5) :: zone
    INTEGER :: tvalues(8)
    CALL DATE_AND_TIME(date, time, zone, tvalues)
    IF (dbg_enabled(level)) THEN
      WRITE (log_unit,'(A)') '[DBG-'//date//'T'//time//':'//trim(level_tag(level))//'] '//TRIM(msg)
      CALL flush_log()
    END IF
  END SUBROUTINE dbg_print

  SUBROUTINE assert_ok(cond, msg)
    LOGICAL,      INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN), OPTIONAL :: msg
    IF (.NOT. cond) THEN
      IF (PRESENT(msg)) THEN
        WRITE (log_unit,'(A)') 'ASSERTION FAILED: '//TRIM(msg)
      ELSE
        WRITE (log_unit,'(A)') 'ASSERTION FAILED.'
      END IF
      CALL flush_log()
      STOP 1
    END IF
  END SUBROUTINE assert_ok

  SUBROUTINE debug_set_logfile(filename)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER :: ios
    IF (log_opened) CLOSE(log_unit)
    OPEN(NEWUNIT=log_unit, FILE=TRIM(filename), STATUS='UNKNOWN', IOSTAT=ios)
    IF (ios /= 0) THEN
      ! Fall back to stdout on failure
      log_unit = 6
      log_opened = .FALSE.
    ELSE
      log_opened = .TRUE.
    END IF
  END SUBROUTINE debug_set_logfile

  SUBROUTINE debug_init_from_args()
    ! Enable via:
    !   --debug              (level defaults to 1)
    !   --debug=2            (set level)
    !   --log=filename.txt   (send logs to file)
    INTEGER :: n, i, eqpos, ival, ios
    CHARACTER(STRLEN) :: arg, val

    n = IARGC()
    DO i = 1, n
      CALL GETARG(i, arg)
      CALL trim_inplace(arg)

      IF (starts_with(arg, '--debug')) THEN
        eqpos = INDEX(arg, '=')
        IF (eqpos > 0) THEN
          val = arg(eqpos+1:)
          CALL trim_inplace(val)
          READ(val, *, IOSTAT=ios) ival
          IF (ios == 0) THEN
            CALL debug_set(.TRUE., ival)
          ELSE
            CALL debug_set(.TRUE., 1)
          END IF
        ELSE
          CALL debug_set(.TRUE., 1)
        END IF

      ELSEIF (starts_with(arg, '--log=')) THEN
        val = arg(7:)
        CALL trim_inplace(val)
        IF (LEN_TRIM(val) > 0) CALL debug_set_logfile(val)
      END IF
    END DO
  END SUBROUTINE debug_init_from_args

  ! -------- helpers (F95-friendly) --------

  SUBROUTINE trim_inplace(s)
    CHARACTER(*), INTENT(INOUT) :: s
    INTEGER :: lt
    lt = LEN_TRIM(s)
    IF (lt < LEN(s)) s(lt+1:) = ' '
  END SUBROUTINE trim_inplace

  LOGICAL FUNCTION starts_with(s, prefix)
    CHARACTER(*), INTENT(IN) :: s, prefix
    INTEGER :: lp

    lp = LEN_TRIM(prefix)
    IF (LEN_TRIM(s) < lp) THEN
      starts_with = .FALSE.
    ELSE
      starts_with = (s(1:lp) == prefix(1:lp))
    END IF
  END FUNCTION starts_with

  CHARACTER(6) FUNCTION level_tag(level)
    INTEGER, INTENT(IN) :: level
    WRITE(level_tag,'(I0)') level
  END FUNCTION level_tag

  SUBROUTINE flush_log()
    ! FLUSH is F2003, but many F95 compilers supported it; if not available,
    ! comment the next line out or keep it inside a vendor conditional.
    FLUSH(log_unit)
  END SUBROUTINE flush_log

end module debug