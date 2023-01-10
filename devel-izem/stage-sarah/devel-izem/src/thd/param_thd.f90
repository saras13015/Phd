module paramthd


  implicit none


  integer, parameter :: IA = kind ( 0                 )     !4          !kind(1)
  integer, parameter :: RA = kind ( 0.0d0             )     !8          !kind(1.0d0)
  integer, parameter :: CA = kind ( ( 0.0d0 , 0.0d0 ) )     !(8,8)      !kind((1.0d0,1.0d0))


  real(RA), parameter :: epsm2_asp   = 1.e-2_RA
  real(RA), parameter :: epsm3_asp   = 1.e-3_RA
  real(RA), parameter :: epsm4_asp   = 1.e-4_RA
  real(RA), parameter :: epsm5_asp   = 1.e-5_RA
  real(RA), parameter :: epsm6_asp   = 1.e-6_RA
  real(RA), parameter :: epsm8_asp   = 1.e-8_RA
  real(RA), parameter :: epsm10_asp  = 1.e-10_RA
  real(RA), parameter :: epsm12_asp  = 1.e-12_RA


  integer (IA) , parameter :: Ndim_thd          = 2
  integer (IA) , parameter :: Nspcmax_thd       = 29
  integer (IA) , parameter :: Nelemax_thd       = 10
  integer (IA) , parameter :: Npolymax_thd      = 7
  integer (IA) , parameter :: NTranspolymax_thd = 4
  integer (IA) , parameter :: Length_thd        = 50


end module paramthd
