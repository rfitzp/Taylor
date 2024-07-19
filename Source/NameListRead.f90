! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read TAYLOR namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (Q_LARGE_E, Q_SMALL_E, Q_I, D,&
     P_PHI, P_PERP, IOTA_E, LARGE_SIGMA,&
     TMAX, NT,&
     SMALL_SIGMA, OMAX, NO,&
     PSTART, PEND, NP)&
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  real    (kind = c_double), intent (inout) :: Q_LARGE_E
  real    (kind = c_double), intent (inout) :: Q_SMALL_E
  real    (kind = c_double), intent (inout) :: Q_I
  real    (kind = c_double), intent (inout) :: D
  real    (kind = c_double), intent (inout) :: P_PHI
  real    (kind = c_double), intent (inout) :: P_PERP
  real    (kind = c_double), intent (inout) :: IOTA_E
  real    (kind = c_double), intent (inout) :: LARGE_SIGMA

  real    (kind = c_double), intent (inout) :: TMAX
  integer (kind = c_int),    intent (inout) :: NT

  real    (kind = c_double), intent (inout) :: SMALL_SIGMA
  real    (kind = c_double), intent (inout) :: OMAX
  integer (kind = c_int),    intent (inout) :: NO

  real    (kind = c_double), intent (inout) :: PSTART
  real    (kind = c_double), intent (inout) :: PEND
  integer (kind = c_int),    intent (inout) :: NP

  namelist /TAYLOR_CONTROL/ Q_LARGE_E, Q_SMALL_E, Q_I, D,&
       P_PHI, P_PERP, IOTA_E, LARGE_SIGMA,&
       TMAX, NT,&
       SMALL_SIGMA, OMAX, NO,&
       PSTART, PEND, NP

  open  (unit = 100, file = 'Inputs/Namelist.nml', status = 'old')
  read  (unit = 100, nml  = TAYLOR_CONTROL)
  close (unit = 100)

endsubroutine NameListRead
