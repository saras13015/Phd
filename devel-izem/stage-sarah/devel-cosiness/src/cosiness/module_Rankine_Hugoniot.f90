!------------------------------------------------------------------------------
!MODULE: Rankine_Hugoniot
!------------------------------------------------------------------------------
!> \brief Provide Rankine-Hugoniot relations for a multi-especies gas mixture.
!!
!! Resolution of pre- and post-shock conditions for an oblique shock
!! wave, with a angle of incidence \f$ \beta \f$ from the ground, and
!! spreads left (index 1) to right (index 2) at an average speed \f$
!! v_1 \f$. In the reference frame associated with shock, in STEADY
!! state, we calculate the flow thermodynamic state after the shock
!! passage \f$ (p_2, T_2, v_2) \f$ and the angle of deflection of the
!! velocity after impact \f$ \theta \f$.
!!
!! To resolve Rankine-Hugoniot's equations, we use the
!! Newton-Raphson's iterative method with the initial value \f$
!! \phi=T_2/T_1 \f$, calculated with Rankine-Hugoniot's relation for
!! an ideal gas and thermodynamics properties independent to
!! temperature. If this method fails, we then use the dichotomy
!! method.
!!
!!                          shock
!!                            ||
!!          <----             ||                <----
!!           (2)              ||                 (1)
!!      u2, T2, P2, rho2      ||           u1, T1, P1, rho1
!!
!! References:
!!
!! -# R. E. Mitchell, R. J. Kee. _A general-purpose computer code for
!!    predicting chemical kinetic behavior behind incident and
!!    reflected shocks_. Sandia Report, p. 13, (1982).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module Rankine_Hugoniot

  use parameters
  use type_thd
  use thermodynamics

  implicit none

  integer (ip) , parameter , private        :: itmax   =  50        ! maximum number of iterations for the method
  real (dp)    , parameter , private        :: eps_phi =  1.0e-6_dp ! solution error tolerance
  real (dp)    , parameter , private        :: dphi    =  1.0e-5_dp ! small epsilon to add to the current value of phi
  real (dp)    , parameter , private        :: neg_tol = -1.0e-2_dp ! maximum negative tolerance to avoid square of negative numbers


contains


!> \brief Selector to update the boundary.
!!
!! This subroutine reads the keyword to select the corresponding boundary condition.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine RH_variables ( thd , Y1 , p1 , T1 , rho1 , v1 , h1 , beta_ , &
                            p2 , T2 , rho2 , v2 , h2 , theta_ )


    type (thd_type) , intent (in)                            :: thd    !< thermodynamic derived type
    real (dp) , dimension (:) , intent (in)                  :: Y1     !< species mass fractions before (and after) shock
    real (dp) , intent (in)                                  :: p1     !< pressure before shock
    real (dp) , intent (in)                                  :: T1     !< temperature before shock
    real (dp) , intent (in)                                  :: v1     !< velocity before shock
    real (dp) , intent (in)                                  :: beta_  !< angle between wall and the shock
    real (dp) , intent (inout)                               :: rho1   !< density before shock
    real (dp) , intent (inout)                               :: h1     !< enthalpy before shock
    real (dp) , intent (inout)                               :: p2     !< pressure after shock
    real (dp) , intent (inout)                               :: T2     !< temperature after shock
    real (dp) , intent (inout)                               :: v2     !< velocity after shock
    real (dp) , intent (inout)                               :: rho2   !< density after shock
    real (dp) , intent (inout)                               :: h2     !< enthalpy after shock
    real (dp) , intent (inout)                               :: theta_ !< angle between wall and velocity after shock


    logical                                 :: dichotomy
    integer (ip)                            :: it , l
    real (dp)                               :: W_i_1 , cp_1 , u1 , gamma1 , cs1 , M1n
    real (dp)                               :: u2
    real (dp)                               :: phi   , phi1   , phi2   , & ! = T2/T1
                                               psi   , psi1   , psi2   , & ! = p2/p1
                                               F_phi , F_phi1 , F_phi2
    real (dp)                               :: dF_phi !, dpsi
    real (dp) , dimension (nrv)             :: ha1 , ha2


    dichotomy = .false. ! use only when Newton's method fails


    ! Initilize the state 1


    call Wmix_i_scalar ( thd , Y1 , W_i_1 )
    call cp_scalar     ( thd , T1 , Y1 , cp_1 )

    rho1   = P1 / ( T1 * W_i_1 )
    u1     = v1 * sin (beta_)             ! normal component of velocity relative to the shock
    gamma1 = thd % gam2 * W_i_1 / cp_1
    gamma1 = 1.0_dp / ( 1.0_dp - gamma1 ) ! = (1-(gamma-1)/gamma * W_i_1 / cp_1)^-1
    cs1    = sqrt ( gamma1 * p1 / rho1 )  ! speed of sound
    M1n    = u1 / cs1                     ! normal Mach number
    phi    = RH_phi0 ( gamma1 , M1n )     ! temperature ratio initialisation

    call ha_scalar ( thd , T1 , ha1 )
    h1 = 0.0_dp
    do l = 1 , nrv
        h1 = h1 + ha1 (l) * Y1 (l)
    end do


    ! Newton-Raphson's iterative method to calculate the value of phi


    it = 0
    do


       it = it + 1
       if ( it > itmax ) then
          if ( rank == rank_default ) then
             write (*,*) 'WARNING: failure in phi convergence in Newton method'
             write (*,*) 'Max number of iterations' , itmax
             write (*,*) 'error' , dF_phi
          end if
          dichotomy = .true. ! if Newton's method fails, try using dichiotomy
          exit
       end if


       ! METHOD 1: Resolution using the analytic derivative.
       ! Sometimes, in some cases, using the analytical derivative,
       ! the Newton's method failed because at one moment, it
       ! calculates a root of one negative term. The numerical
       ! derivative (method 2) is more stable (A.Techer)


       ! ! mixing enthalpies (thd, temp, nrv, Ya) and its derivative
       ! call ha_scalar ( thd , T1 * phi , ha2 )
       ! h2 = 0.0_dp
       ! do l = 1 , nrv
       !    h2 = h2 + ha2 (l) * Y1 (l)
       ! end do

       ! ! psi and its derivative
       ! call RH_psi ( rho1 , u1 , p1 , phi , psi , dichotomy )
       ! dpsi = RH_dpsi ( rho1 , u1 , p1 , phi )

       ! ! F and its derivative
       ! F_phi  = RH_F_phi  ( h1 , h2 , u1 , phi , psi )
       ! dF_phi = RH_dF_phi ( u1 , phi , psi , dpsi )

       ! ! error function estimation
       ! if ( abs ( F_phi / dF_phi ) < eps_phi ) then
       !    exit
       ! end if

       ! ! update temperature ratio
       ! phi = phi - F_phi / dF_phi


       ! METHOD 2:  Resolution using the numerical derivative


       call ha_scalar ( thd , T1 * phi , ha2 )
       h2 = 0.0_dp
       do l = 1 , nrv
          h2 = h2 + ha2 (l) * Y1 (l)
       end do

       call RH_psi ( rho1 , u1 , p1 , phi , psi , dichotomy )
       if (dichotomy) exit

       F_phi = RH_F_phi ( h1 , h2 , u1 , phi , psi )


       call ha_scalar ( thd , T1 * (phi+dphi) , ha2 )
       h2 = 0.0_dp
       do l = 1 , nrv
           h2 = h2 + ha2 (l) * Y1 (l)
       end do

       call RH_psi ( rho1 , u1 , p1 , phi+dphi , psi2 , dichotomy )
       if (dichotomy) exit

       F_phi2 = RH_F_phi ( h1 , h2 , u1 , phi+dphi , psi2 )


       dF_phi = (F_phi2 - F_phi) / dphi


       ! error function estimation
       if ( abs ( F_phi / dF_phi ) < eps_phi ) exit


       ! update temperature ratio
       phi = phi - F_phi / dF_phi


    end do


    ! Dichotomy iterative method to calculate the value of phi (only when Newton's fails)


    if (dichotomy) then


       dichotomy = .false.                  ! if the two methods fail, then dichotomy will become true again to stop the program
       phi1 = 1.0_dp                        ! minimum physical value of phi (T2/T1) == unity
       phi2 = RH_max_phi ( rho1 , u1 , p1 ) ! maximum ?physical? value of phi (T2/T1)


       call ha_scalar ( thd , T1 * phi1 , ha2 )
       h2 = 0.0_dp
       do l = 1 , nrv
          h2 = h2 + ha2 (l) * Y1 (l)
       end do

       call RH_psi ( rho1 , u1 , p1 , phi1 , psi1 , dichotomy )
       if (dichotomy) call abort_mpi ('ERROR: Dichitomy has failed at psi1')

       F_phi1 = RH_F_phi ( h1 , h2 , u1 , phi1 , psi1 )


       call ha_scalar ( thd , T1 * phi2 , ha2 )
       h2 = 0.0_dp
       do l = 1 , nrv
          h2 = h2 + ha2 (l) * Y1 (l)
       end do

       call RH_psi ( rho1 , u1 , p1 , phi2 , psi2 , dichotomy )
       if (dichotomy) call abort_mpi ('ERROR: Dichitomy has failed at psi2')

       F_phi2 = RH_F_phi ( h1 , h2 , u1 , phi2 , psi2 )


       ! checking interval for dichotomy method
       if ( F_phi1 * F_phi2 > 0.0_dp ) call abort_mpi ('ERROR: No zero found between phi_min and phi_max')


       it = 0
       do


          it = it + 1
          if ( it > itmax ) then
             if ( rank == rank_default ) then
                write (*,*) 'ERROR: failure in phi convergence in Dichotomy method'
                write (*,*) 'Max number of iterations' , itmax
                write (*,*) 'error' , dF_phi
             end if
             call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
          end if


          phi = 0.5_dp * ( phi1 + phi2 )


          call ha_scalar ( thd , T1 * phi , ha2 )
          h2 = 0.0_dp
          do l = 1 , nrv
             h2 = h2 + ha2 (l) * Y1 (l)
          end do

          call RH_psi ( rho1 , u1 , p1 , phi , psi , dichotomy )
          if (dichotomy) call abort_mpi ('ERROR: Dichitomy has failed at psi_i')

          F_phi = RH_F_phi ( h1 , h2 , u1 , phi , psi )


          if ( abs(F_phi) < eps_phi ) exit


          if ( F_phi1 * F_phi < 0.0_dp ) then
             phi2 = phi ; F_phi2 = F_phi ! zero is left
          else
             phi1 = phi ; F_phi1 = F_phi ! zero is right
          endif


       end do


    end if


    ! calculate conditions at state 2


    T2     = T1 * phi
    p2     = p1 * psi
    rho2   = P2 / ( T2 * W_i_1 )
    u2     = rho1 * u1 / rho2                       ! normal component of velocity relative to the shock
    theta_ = beta_ - atan ( u2 * tan (beta_) / u1 )
    v2     = u2 / sin ( beta_ - theta_ )            ! velocity magnitude


  end subroutine RH_variables


!> \brief Calculate the initial value of the temperature ratio.
!!
!! This uses the Rankine-Hugoniot's relations for a pure air flow
!! assuming that thermodynamics properties are independant to
!! temperature.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function RH_phi0 ( gamma1 , M1n )


    real (dp) , intent (in)                                        :: gamma1  !< ratio of specific heats
    real (dp) , intent (in)                                        :: M1n     !< normal component of the Mach relative to the shock
    real (dp)                                                      :: RH_phi0 !< temperature ratio


    real (dp)                                                      :: wrk


    wrk = 0.5_dp * ( gamma1 + 1.0_dp ) * M1n

    RH_phi0 = gamma1 * M1n * M1n - 0.5_dp * ( gamma1 - 1.0_dp )
    RH_phi0 = RH_phi0 * ( 1.0_dp + 0.5_dp * ( gamma1 - 1.0_dp ) * M1n * M1n )
    RH_phi0 = RH_phi0 / ( wrk * wrk )


  end function RH_phi0


!> \brief Calculate the maximum ?physical? value of the temperature ratio.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function RH_max_phi ( rho1 , u1 , p1 )


    real (dp) , intent (in)                                        :: rho1       !< density
    real (dp) , intent (in)                                        :: u1         !< velocity
    real (dp) , intent (in)                                        :: p1         !< pressure
    real (dp)                                                      :: RH_max_phi !< temperature ratio


    RH_max_phi = 1.0_dp + rho1 * u1 * u1 / p1
    RH_max_phi = RH_max_phi * RH_max_phi
    RH_max_phi = RH_max_phi * p1 / (4.0_dp * rho1 * u1 * u1 )


  end function RH_max_phi


!> \brief Calculate the pressure ratio.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine RH_psi ( rho1 , u1 , p1 , phi , psi , dichotomy )


    real (dp) , intent (in)                                        :: rho1      !< density
    real (dp) , intent (in)                                        :: u1        !< velocity
    real (dp) , intent (in)                                        :: p1        !< pressure
    real (dp) , intent (in)                                        :: phi       !< temperature ratio
    real (dp) , intent (inout)                                     :: psi       !< pressure ratio
    logical , intent (inout)                                       :: dichotomy !< pressure ratio


    psi = 1.0_dp + rho1 * u1 * u1 / p1
    psi = psi * psi
    psi = psi - 4.0_dp * rho1 * u1 * u1 * phi / p1


    ! filter negative pressure ratio
    if ( psi < neg_tol ) then
       if ( rank == rank_default ) write (*,*) 'WARNING: ratio of pressure is too negative ' , psi
       dichotomy = .true.
    else if ( psi > neg_tol .and. psi < 0.0_dp ) then
       if ( rank == rank_default ) write (*,*) 'WARNING: ratio of pressure is a bit negative ' , psi
       psi = 0.0_dp
    end if


    psi = 0.5_dp * ( sqrt (psi) + 1.0_dp + rho1 * u1 * u1 / p1 )


  end subroutine RH_psi


!> \brief Calculate the derivative of the pressure ratio.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function RH_dpsi ( rho1 , u1 , p1 , phi )


    real (dp) , intent (in)                                        :: rho1    !< density
    real (dp) , intent (in)                                        :: u1      !< velocity
    real (dp) , intent (in)                                        :: p1      !< pressure
    real (dp) , intent (in)                                        :: phi     !< temperature ratio
    real (dp)                                                      :: RH_dpsi !< derivative of the pressure ratio


    RH_dpsi = 1.0_dp + rho1 * u1 * u1 / p1
    RH_dpsi = RH_dpsi * RH_dpsi
    RH_dpsi = RH_dpsi - 4.0_dp * rho1 * u1 * u1 * phi / p1


    ! when using the analytical derivative, the ratio can not be zero nor negative at all!
    if ( RH_dpsi <= 0.0_dp ) call abort_mpi ('ERROR: ratio of pressure derivative is negative') 


    RH_dpsi = 1.0_dp / ( sqrt (RH_dpsi) * p1 )
    RH_dpsi = RH_dpsi * ( - rho1 * u1 * u1 )


  end function RH_dpsi


!> \brief Calculate the error function.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function RH_F_phi ( h1 , h2 , u1 , phi , psi )


    real (dp) , intent (in)                                        :: h1       !< enthalpy at state 1
    real (dp) , intent (in)                                        :: h2       !< enthalpy at state 2
    real (dp) , intent (in)                                        :: u1       !< velocity at state 1
    real (dp) , intent (in)                                        :: phi      !< temperature ratio
    real (dp) , intent (in)                                        :: psi      !< pressure ratio
    real (dp)                                                      :: RH_F_phi !< error function


    RH_F_phi = h1 - h2 + 0.5_dp * u1 * u1 * ( 1.0_dp - phi*phi / (psi*psi) )


  end function RH_F_phi


!> \brief Calculate the derivative of the error function.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function RH_dF_phi ( u1 , phi , psi , dpsi )


    real (dp) , intent (in)                                        :: u1        !< density
    real (dp) , intent (in)                                        :: phi       !< temperature ratio
    real (dp) , intent (in)                                        :: psi       !< pressure ratio
    real (dp) , intent (in)                                        :: dpsi      !< derivative of the pressure ratio
    real (dp)                                                      :: RH_dF_phi !< derivative of the error function


    RH_dF_phi = phi * phi * dpsi / ( psi * psi * psi ) - &
                phi / ( psi * psi )

    RH_dF_phi = RH_dF_phi * u1 * u1


  end function RH_dF_phi


end module Rankine_Hugoniot
