!------------------------------------------------------------------------------
! MODULE: solver_reaction
!------------------------------------------------------------------------------
!> \brief Solver for the chemical production rate in species equations.
!!
!! This module contains the subroutines and utilities needed to
!! solve the chemical production rate $\omega_\alpha$ by using:\n
!!   * the PSR (Perfectly Stired Reactor) with vode solver,\n
!!   * the MIL (Model Intermittent Lagrangien).\n
!!   * the hybrid PSR-MIL.\n
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module solver_reaction

  use parameters
  use parallel
  use input
  use adim
  use type_thd
  use BCs
  use thermodynamics , only : upd_prim_var_woT
  use tools

  implicit none

  ! tolerances for mass fraction : production rate, fuel and oxidizer. For equivalent ratio
  real (dp) , parameter , private :: eps_w = 1.0e-15_dp , eps_f = 1.0e-2_dp , eps_o = 1.0e-6_dp , &
                                     phi_min = 1.0e-20_dp , phi_max = 1.0e6_dp

  real (dp) , parameter , private :: eps8 = 1.0e-8_dp , eps4 = 1.0e-4_dp , eps3 = 1.0e-3_dp , eps2 = 1.0e-2_dp , eps = 1.0e-15_dp


contains


!> \brief Selector for the combustion model/solver
!!
!! Reads the keyword to select the corresponding combustion solver.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine reaction_selec ( inp , thd , adi , grid , time , dt , mu_SGS , T , W_i , cp , ha , v , DTmax )


    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid   !< grid derived type
    real (dp) , intent (in)                                        :: time   !< time
    real (dp) , intent (in)                                        :: dt     !< time step
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: mu_SGS !< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T      !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp     !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha     !< species enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v      !< conserved variables array
    real (dp) , intent (inout)                                     :: DTmax  !< maximum temperature increase


    if      ( inp % reac_selec == PSR ) then
       call reaction_PSR ( inp , thd , adi , dt , W_i , T , v , DTmax )

    else if ( inp % reac_selec == MIL ) then
       call reaction_MIL ( inp , thd , adi , grid , dt , W_i , mu_SGS , T , v )
       call upd_prim_var_domain ( thd , T , W_i , cp , ha , v )

    else if ( inp % reac_selec == hybrid_PSR_MIL ) then
       call reaction_MIL ( inp , thd , adi , grid , dt , W_i , mu_SGS , T , v )
!       call upd_prim_var_domain ( thd , T , W_i , cp , ha , v )
       call reaction_PSR ( inp , thd , adi , dt , W_i , T , v , DTmax )

    else
       call abort_mpi ('Reaction selector ' // trim (inp % reac_selec) // ' not defined')

    end if

    call upd_prim_var_woT ( thd , 0 , v , T , W_i , cp , ha )
    call comm_cons (v)
    call upd_boundaries ( inp , adi , thd , time , dt , grid , T , W_i , cp , ha , v )
    call upd_prim_var_ghost ( thd , T , W_i , cp , ha , v )


  end subroutine reaction_selec


!> \brief Reaction domain
!!
!! Defines the reaction domain boundaries.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine domain_reaction ( ix , fx , iy , fy , iz , fz )


    integer (ip) , intent (inout) :: ix , fx !< start/end points in x-direction
    integer (ip) , intent (inout) :: iy , fy !< start/end points in y-direction
    integer (ip) , intent (inout) :: iz , fz !< start/end points in z-direction


    ! domain to apply combustion (not in every boundary)
    if ( ndim >= 1 ) then
      ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 ) ! added ng
!      ix = sx                ; fx = ex
      iy = sy                ; fy = ey
      iz = sz                ; fz = ez
    end if 

    if ( ndim >= 2 ) then
!      iy = max ( ng+1 , sy ) ; fy = min ( ey , nty-ng-1 ) ! added ng
      iy = sy                ; fy = min ( ey , nty-ng-1 ) ! added ng
!      iy = sy                ; fy = ey
      iz = sz                ; fz = ez
    end if

    if ( ndim == 3 ) then
      iz = max ( ng+1 , sz ) ; fz = min ( ez , ntz-ng-1 ) ! added ng
    end if


  end subroutine domain_reaction


!> \brief Reaction detector.
!!
!! New criteria that seems to work well with the presence or not of
!! radicals:
!!
!! #- The fuel stream is composed by hydrogen and nitrogen. If the
!!    inlet quantities are changed we are not in the fuel stream but
!!    in a mixture region.
!! #- The oxidizer stream is always composed by a certain amount of
!!    nitrogen and zero hydrogen.
!!
!! If these two quantities change, then we are not in the oxidizer
!! stream.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine reaction_det ( T , rho , P , Ya , dt , call_dvode )

    real (dp)                       , intent (in)     :: T          !< temperature in CGS system [K]
    real (dp)                       , intent (in)     :: rho        !< density in CGS system [g/cm^3]
    real (dp)                       , intent (in)     :: P          !< pressure in CGS system [dynes/cm^2=g/cm/s^2]
    real (dp) , dimension (1:nrv)   , intent (in)     :: Ya         !< mass fraction [g/g]
    real (dp)                       , intent (in)     :: dt         !< time step in CGS system [s]
    logical                         , intent (inout)  :: call_dvode !< true if the PSR model should be called

    real (dp) , dimension (1:nrv)  :: DeltaY      ! mass fractions variation array: Yfinal-Yinitial [g/g]
    integer (ip)  :: l
    logical       :: call_dvode_w , call_dvode_f , call_dvode_o



    ! Initialize some variables
    call_dvode   = .false.
    call_dvode_w = .false.
    call_dvode_f = .false.
    call_dvode_o = .false.


!    ! Verify the production rate of each species
!    DeltaY (:) = Ya (:)
!    call psr_simp_chemkin ( T , rho , DeltaY , dt )
!    do l = 1 , nrv
!       if ( abs ( DeltaY(l) ) > eps_w ) call_dvode_w = .true.
!    end do

!    ! Verify the pure fuel
!    if ( abs ( Ya(inp % index_H2) - inp % Y1f(inp % index_H2) ) > eps_f .and. &
!         abs ( Ya(inp % index_N2) - inp % Y1f(inp % index_N2) ) > eps_f ) call_dvode_f = .true.

!    ! Verify the pure oxidizer
!    if ( abs ( Ya(inp % index_H2) - inp % Y0o(inp % index_H2) ) > eps_o .and. &
!         abs ( Ya(inp % index_N2) - inp % Y0o(inp % index_N2) ) > eps_o ) call_dvode_o = .true.

    ! Verify the equivalent ratio
!    phi = 8.0_dp * Ya(inp % index_H2) / ( Ya(inp % index_O2) + eps_o )
!    if ( phi_min < phi .and. phi < phi_max ) call_dvode_f = .true.


    ! Verify if there is chemical reaction or not
!    if ( call_dvode_w )       call_dvode = .true.
!    if ( .not. call_dvode_f ) call_dvode = .false.
!    if ( .not. call_dvode_o ) call_dvode = .false.


    ! Verify the production rate of each species (A.Techer: This only criteria seems to work well also: checked !)
    DeltaY (:) = Ya (:)
    call psr_simp_chemkin ( T , rho , P , DeltaY , dt )
    do l = 1 , nrv
       if ( abs ( DeltaY(l) ) > eps_w ) call_dvode = .true. ; exit
    end do


  end subroutine reaction_det


!> \brief Calculate the reactive system of equations by using PSR model 
!! and DVODE solver.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine reaction_PSR ( inp , thd , adi , dt , W_i , T , v , DTmax )


    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    real (dp) , intent (in)                                        :: dt    !< time step
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i   !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v     !< conserved variables array
    real (dp) , intent (inout)                                     :: DTmax !< maximum temperature increase


    integer (ip)                  :: ix , fx , &
                                     iy , fy , &
                                     iz , fz
    integer (ip)                  :: i , j , k , l , istate
    real (dp)                     :: dt_phys , rho_phys , T_phys , P_phys , rho_i , rho0 , T0_phys , sumY , Tref_i , rho_cgs , P_cgs
    real (dp)                     :: segz
    logical                       :: call_dvode_subdomain , call_dvode
    real (dp) , dimension (1:nrv+npv+nvv) :: Ya

! ==== A.Techer test =======================
    real (dp)                     :: phi
    real (dp) , dimension (1:nrv) :: Xa

    ! New criteria that seems to work well with the presence or not of radicals:
    ! 1. The fuel stream is composed by hydrogen and nitrogen. If the inlet quantities are
    ! changed we are not in the fuel stream but in a mixture region.
    ! 2. The oxidizer stream is always composed by a certain amount of nitrogen and zero hydrogen.
    ! If these two quantities varies then we are not in the oxidizer stream.


    dt_phys = dt * adi % time_ref
    Tref_i  = 1.0_dp / adi % T_ref
    DTmax   = eps
    rho_cgs = adi % rho_ref * 1.0e-3_dp
    P_cgs   = adi % P_ref * 10.0_dp


    call domain_reaction ( ix , fx , iy , fy , iz , fz )


    if ( inp % filter_reaction ) then ! DVODE solver will be applied at certains points to save computational costs


       call_dvode_subdomain = .false. ! useless ?


       do k = iz , fz
          do j = iy , fy
             do i = ix , fx


                rho_i = 1.0_dp / v (i,j,k,1)
                do l = 1 , nrv+npv+nvv
                   Ya (l) = v (i,j,k, niv+l ) * rho_i
                end do


                rho0     = v (i,j,k,1)
                rho_phys = rho0 * rho_cgs
                T0_phys  = T (i,j,k) * adi % T_ref
                T_phys   = T0_phys
                P_phys   = rho0 * T (i,j,k) * W_i (i,j,k) * P_cgs


                call reaction_det ( T_phys , rho_phys , P_phys , Ya , dt_phys , call_dvode )


                if (call_dvode) then


                   call_dvode_subdomain = .true. ! useless ?


                   segz = 0.0_dp

                   if ( inp % reac_selec == hybrid_PSR_MIL .and. npv > 0 .and. nvv > 0 ) then
                      segz = ( Ya (nrv+1) + eps8 ) * ( 1.0_dp - Ya (nrv+1) - eps8 )
                      segz = Ya (nrv+npv+1) / segz
                      segz = max ( min ( segz , 1.0_dp ) , 0.0_dp )
                   end if


                   call psr_chemkin ( T_phys , rho_phys , P_phys , Ya , dt_phys , segz , istate )
                   if ( istate <= -2 ) then !.and. istate /= -5 ) then
                      write (*,'(A,3(1X,I5))') "error: failed integration of chemical source term in DVODE at (i,j,k) = " , &
                                   i , j , k

                      Xa (1:nrv) = Ya (1:nrv)
                      call Y_to_X ( thd , Xa )
                      phi = 8.0_dp * Ya(inp % index_H2) / Ya(inp % index_O2)
                      write (*,'(A,F15.2,F10.2,10(1PE15.5))') "P, T, phi, X = ", &
                      &     rho0*T (i,j,k)*W_i(i,j,k)*adi%P_ref , T0_phys , phi, Xa (1:nrv)
                      call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
                   end if


                   T (i,j,k) = T_phys * Tref_i


                   ! DTmax has physical dimensions: degrees Kelvin
                   DTmax = max ( DTmax , abs ( T_phys - T0_phys ) ) ! T_phys is updated


                   sumY = 0.0_dp
                   do l = 1 , nrv-1
                      Ya (l)          = min ( max ( Ya (l) , 0.0_dp ) , 1.0_dp )
                      v (i,j,k,niv+l) = rho0 * Ya (l)
                      sumY            = sumY + Ya (l)
                   end do
                   Ya (nrv)          = 1.0_dp - sumY
                   Ya (nrv)          = min ( max ( Ya(nrv) , 0.0_dp ) , 1.0_dp )
                   v (i,j,k,niv+nrv) = rho0 * Ya (nrv)


                end if


             end do
          end do
       end do


    else ! DVODE solver is applied in every point of the computational domain


       do k = iz , fz
          do j = iy , fy
             do i = ix , fx


                rho_i = 1.0_dp / v (i,j,k,1)
                do l = 1 , nrv+npv+nvv
                   Ya (l) = v (i,j,k, niv+l ) * rho_i
                end do

                rho0     = v (i,j,k,1)
                rho_phys = rho0 * rho_cgs ! in g/cm^3
                T0_phys  = T (i,j,k) * adi % T_ref
                T_phys   = T0_phys
                P_phys   = rho0 * T (i,j,k) * W_i (i,j,k) * p_cgs


                segz = 0.0_dp

                if ( inp % reac_selec == hybrid_PSR_MIL .and. npv > 0 .and. nvv > 0 ) then
                   segz = ( Ya (nrv+1) + eps8 ) * ( 1.0_dp - Ya (nrv+1) - eps8 )
                   segz = Ya (nrv+npv+1) / segz
                   segz = max ( min ( segz , 1.0_dp ) , 0.0_dp )
                end if


                call psr_chemkin ( T_phys , rho_phys , P_phys , Ya , dt_phys , segz , istate )
                if ( istate <= -2 ) then !.and. istate /= -5 ) then
                   write (*,'(A,3(1X,I5))') "error: failed integration of chemical source term in DVODE at (i,j,k) = " , &
                                i , j , k

                   Xa (1:nrv) = Ya (1:nrv)
                   call Y_to_X ( thd , Xa )
                   phi = 8.0_dp * Ya(inp % index_H2) / Ya(inp % index_O2)
                   write (*,'(A,F15.2,F10.2,10(1PE15.5))') "P, T, phi, X = ", &
                   &    rho0*T (i,j,k)*W_i(i,j,k)*adi%P_ref , T0_phys , phi, Xa (1:nrv)
                   call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
                end if


                T (i,j,k) = T_phys * Tref_i


                ! DTmax has physical dimensions: degrees Kelvin
                DTmax = max ( DTmax , abs ( T_phys - T0_phys ) ) ! T_phys is updated


                ! Updating the solution
                sumY = 0.0_dp
                do l = 1 , nrv-1
                   Ya (l)          = min ( max ( Ya (l) , 0.0_dp ) , 1.0_dp )
                   v (i,j,k,niv+l) = rho0 * Ya (l)
                   sumY            = sumY + Ya (l)
                end do
                Ya (nrv)          = 1.0_dp - sumY
                Ya (nrv)          = min ( max ( Ya(nrv) , 0.0_dp ) , 1.0_dp )
                v (i,j,k,niv+nrv) = rho0 * ( 1.0_dp - sumY )


             end do
          end do
       end do


    end if


  end subroutine reaction_PSR


!> \brief Calculate the reactive system of equations by using MIL (Model Intermittent Lagrangien)
!!
!! References:
!! -# A. Mura and J.-F. Izard., Numerical simulation of supersonic nonpremixed turbulent
!!    combustion in a scramjet combustor model, Journal of Propulsion and Power (2010).
!!
!! -# L. Gomet, V. Robin, and A. Mura., Influence of residence and scalar mixing time scales
!!    in non-premixed combustion in supersonic turbulent flows,
!!    Combustion Science and Technology (2012).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine reaction_MIL ( inp , thd , adi , grid , dt , W_i , mu_SGS , T , v )


    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid  !< grid derived type
    real (dp) , intent (in)                                        :: dt    !< time step
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i   !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: mu_SGS!< turbulent viscosity
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v     !< conserved variables array


    integer (ip)                  :: ix , fx , &
                                     iy , fy , &
                                     iz , fz
    integer (ip)                  :: i , j , k , l
    real (dp)                     :: rho , rho_i

    ! Pressure, Temperatures
    real (dp)                     :: T_phys , Tref_i
    real (dp)                     :: P_phys

    ! Mass fractions
    real (dp) , dimension (1:nrv+npv+nvv) :: Ya
    real (dp)                             :: sumY , yom

    ! Mixture fraction, variance, segregation
    real (dp)                     :: zm , zp2 , segz

    ! Chemical production rate
    real (dp) , dimension (1:nrv)   :: wreac

    ! Others
    real (dp)                     :: dm_SGS , Sc_sgs_i , taumix , rvd0 , wrk1
    real (dp) , parameter         :: theta = 0.0_dp , & ! implicitation factor (0 for explicite, 1 for implicite, 0.5 for Crank-Nicolson)
                                     thetam1 = theta - 1.0_dp

    ! For additionnal output variables
    integer (ip)  :: call_reac
    real (dp)     :: c , zjm , zjp , int_zjmp_pdfbeta , yoequm , yoburnt , &
                     prob_mil , wo2mil , wo2sauts , wo2bur

    if ( npv > 1 ) &
       call abort_mpi ('ERROR: subroutine reaction_MIL npv > 1')


    rvd0     = adi % sqgmr / adi % Sc ! = Ma sqrt(gamma) / ( Re Sc )
    Sc_sgs_i = 1.0_dp / inp % Sc_sgs
    Tref_i   = 1.0_dp / adi % T_ref


    call domain_reaction ( ix , fx , iy , fy , iz , fz )


    do k = iz , fz
       do j = iy , fy
          do i = ix , fx


             rho   = v (i,j,k,1)
             rho_i = 1.0_dp / rho
             do l = 1 , nrv+npv+nvv
                Ya (l) = v (i,j,k, niv+l ) * rho_i
             end do

             yom      = Ya (inp % index_O2)              ! mean oxygen mass fraction
             zm       = Ya (nrv+npv)                     ! mean mixture fractionn value
             zp2      = Ya (nrv+npv+1)                   ! mixture fraction variance value

             segz = 1.0_dp
             if ( inp % reac_selec == hybrid_PSR_MIL ) then
                segz = ( zm + eps8 ) * ( 1.0_dp - zm - eps8 )
                segz = zp2 / segz
                segz = max ( min ( segz , 1.0_dp ) , 0.0_dp )
             end if

             P_phys   = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k) * adi % P_ref ! pressure (Pa)
             T_phys   = T (i,j,k) * adi % T_ref                             ! temperature (K)

             ! Relaxation time ( + eps to avoid divisions by zero )
             dm_SGS   = adi % Sc * mu_SGS (i,j,k) * rho_i * Sc_SGS_i ! (-)
             taumix   = grid % delta (i,j,k) * grid % delta (i,j,k) / ( rvd0 * dm_SGS + eps ) ! (-)

             taumix   = taumix * adi % time_ref ! dimensionalize in second


             call wreac_MIL ( inp , adi , thd , taumix , yom , zm , zp2 , &
                              P_phys , T_phys , & ! for reactiviy rate
                              c , call_reac , zjm , zjp , int_zjmp_pdfbeta , yoequm , yoburnt , & ! useless (because for post-processing)
                              prob_mil , wo2mil , wo2sauts , wo2bur , & ! useless (because for post-processing)
                              wreac )


             if ( call_reac == 1 ) then

                wreac (:) = wreac (:) * adi % time_ref ! adimensionalize (s^-1 * s)
                taumix    = taumix / adi % time_ref

                ! Updating the solution
                sumY = 0.0_dp
                do l = 1 , nrv-1
                   ! TEST avec integration temporelle explicite ordre 1
!                   Ya (l)          = Ya (l) + dt * segz * wreac (l)

                   ! TEST avec integration temporelle linearisation + semi-implicite (Crank-Nicolson)
!                   wrk1            = dt * segz * wreac (l)
!                   wrk2            = 1.0_dp - theta * wrk1 / Ya (l)
!                   Ya (l)          = ( Ya (l) + thetam1 * wrk1 ) / wrk2

                   ! TEST avec integration temporelle linearisation + semi-analytique
                   wrk1            = dt * segz * wreac (l)
                   Ya (l)          = Ya (l) * exp ( wrk1 / Ya (l) )

                   ! TEST integration linearisation + semi-analytique avec introduction tau_mix
!                   wrk1            = - 1.0_dp / taumix
!                   wrk2            = ( segz * wreac (l) - wrk1 * Ya (l) ) / wrk1
!                   Ya (l)          = ( Ya (l) + wrk2 ) * exp ( dt * wrk1 ) - wrk2

                   Ya (l)          = min ( max ( Ya (l) , min_Y ) , max_Y )
                   v (i,j,k,niv+l) = rho * Ya (l)
                   sumY            = sumY + Ya (l)
                end do
                Ya (nrv)          = 1.0_dp - sumY
                Ya (nrv)          = min ( max ( Ya(nrv) , min_Y ) , max_Y )
                v (i,j,k,niv+nrv) = rho * ( 1.0_dp - sumY )

             end if

          end do
       end do
    end do


  end subroutine reaction_MIL


!> \brief Calculate the chemical production rate of each species equations by using MIL
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine wreac_MIL ( inp , adi , thd , taumix , yom , zm , zp2 , &
                         P , T , &
                         c , call_reac , zjm , zjp , int_zjmp_pdfbeta , yoequm , yoburnt , &
                         prob_mil , wo2mil , wo2sauts , wo2bur , &
                         wreac )


    type (inp_type) , intent (in)                     :: inp      !< input derived type
    type (adi_type) , intent (in)                     :: adi      !< non-dimensional derived type
    type (thd_type) , intent (in)                     :: thd      !< thermodynamic derived type
    real (dp) , intent (in)                           :: taumix   !< mixing time scale (s)
    real (dp) , intent (in)                           :: yom      !< oxydizer mass fraction (-)
    real (dp) , intent (in)                           :: zm       !< mixture fraction (-)
    real (dp) , intent (in)                           :: zp2      !< SGS variance (-)

    ! For reactivity rate
    real (dp) , intent (in)                           :: P        !< pressure (Pa)
    real (dp) , intent (in)                           :: T        !< temperature (K)

    ! For additionnal output variables
!    real (dp) , intent (inout)                        :: c         !< other param: progress variable (-)
!    integer (ip) , intent (inout)                     :: call_reac !< other param: binary calling reaction MIL
!    real (dp) , intent (inout)                        :: zjm , zjp !< other param: progress variable (-)
!    real (dp) , intent (inout)                        :: int_zjmp_pdfbeta !< other param: PDF between zjm-zjp
!    real (dp) , intent (inout)                        :: yoequm    !< other param: mean equilibrium mixture fraction (-)
!    real (dp) , intent (inout)                        :: yoburnt   !< other param: mean MIL mixture fraction (-)
!    real (dp) , intent (inout)                        :: prob_mil  !< other param: propotion MIL for conservativity property
!    real (dp) , intent (inout)                        :: wo2mil , wo2sauts , wo2bur !< other param: chemical production rate fraction (1/s)

    real (dp) , dimension (1:nrv) , intent (inout)    :: wreac    !< MAIN OUTPUT: chemical production rate (1/s)


    real (dp) , parameter         :: Da_c_diff = 1.0_dp , & ! critical diffusion     Damkohler number (control parameter)
                                     Da_c_auto = 1.0_dp     ! critical auto-ignition Damkohler number (control parameter)
    integer (ip)                  :: cas , call_reac

    ! Mass fractions
    real (dp)                     :: yoequ , yomix , yoiem , yoequm , yoburnt , yomil , yoinf
    ! Mixture fraction, its variance, jump positions
    real (dp)                     :: zjm , zjp , zp2max , zignim , zignip
    ! Chemical production rate
    real (dp)                     :: wo2mil , wo2bur , wo2sauts , wo2 , prob_mil
    ! Reversed mixing time scale
    real (dp)                     :: taumix_i

    logical                       :: twophase , call_mil , call_self , verbose
    ! Two-phase variables (change of variables)
    real (dp)                     :: max_Y_zi , yfmax_zi , ystmax_zi , & ! mass fraction
                                     zst_zi     ! mixture fraction

    ! Beta-PDF
    real (dp)                     :: xb , beta1 , beta2

    real (dp)                     :: c , & ! progress variable
                                     int_zjmp_pdfbeta ! ignition probability

    ! Epsilon values ! Be caution... they are very sensitive ! (A.Techer)
    real (dp) , parameter         :: eps_z = 1.0e-5_dp
    real (dp) , parameter         :: eps_c = eps3
    real (dp) , parameter         :: eps_s = eps3
    real (dp)                     :: eps_y


    !====================================
    ! MIL marker and initialization
    !====================================
    verbose   = .true.
    twophase  = .false.
    call_mil  = .false.
    call_self = .false.
    call_reac = 0
    wo2       = 0.0_dp
    wreac     = 0.0_dp
    yoinf     = inp % Y0o ( inp % index_O2 ) ! Oxygen mass fraction in oxydizer stream (= 1 if PURE oxygen)
    eps_y     = yoinf * eps_z / zst
!    eps_y     = eps2 * yoinf


    ! === Other output variables (e.g. for post-processing)
    c = 0.0_dp
    yoequm = 0.0_dp
    yoburnt = 0.0_dp
    int_zjmp_pdfbeta = 0.0_dp
    prob_mil = 0.0_dp


    taumix_i = Da_c_diff / taumix

    zp2max   = zm * ( max_Z - zm )
!    zp2      = min ( max ( zp2 , eps ) , zp2max )

    ! === Jump positions initialization
    zjm      = zm
    zjp      = zm

    ! === (Pure) Mixing line Yo at zm
    yomix = yoinf * ( 1.0_dp - zm ) 
    ! yomix = \tilde{Y}^{MIX} = Y^{MIX}(\tilde{\xi}) because Y^{MIX} is linear !!!
    !       = \int_0^1 Y^{MIX}(\xi) \tilde{P}(\xi)d\xi

    ! === Micro-mixing line = IEM (Interaction by Exchange with the Mean) line Yo at zm
    yoiem = yom                                                                                       ! NOT USED !?

    ! === Chemical equilibrium line Yo at zm
    yoequ = 0.0_dp ! These LG (A.3)
    if ( zm <= zst ) yoequ = ( 1.0_dp - zm / zst ) * yoinf ! These LG (A.3)
    yoequ = min ( max ( yoequ , min_Y ) , yoinf ) ! Normally, yoequ is already bounded...


    ! Too close to the boundaries (A.Techer) it happens !
    if ( yom < yoequ + eps_y .or. yom > yomix - eps_y ) return
    if ( yomix - yoequ < eps_y ) return ! used also for progress variable

!    !====================================
!    ! check if AUTOIGNITION
!    !====================================

!    call auto_ign ( zm , yom , yh2 , taumix , Ps , Ec , Htot , ac_crit_MIL , &                       ! ac_crit_MIL non défini c'est pas le Damkohler critique diffusif ? = 1 = parametre de controle ?
!                    zignip , zignim , Da_c_auto , call_self , taumix_i , tignimin , &
!                    tau_sej_ZM , tau_sej_AIR )


!    if ( .not. call_self ) then ! Exit only if NO AUTOIGNITION, (increase algorithm efficiency)

       if ( zm <= ( min_Z + eps_z ) .or. zm >= ( max_Z - eps_z ) ) then
!          if (verbose) write (*,*)'MIL: Too high / small mixture fraction'
          return
       end if

!       if ( ( yomix > eps8 ) .and. ( ( yomix - yoequ ) > eps8 ) ) then
!          ! First  condition: Outside PURE oxydizer stream
!          ! Second condition: Non zero departure from yomix and yoequ
!          c = abs ( ( yomix - yom ) / ( yomix - yoequ + eps ) )                                       ! ORIGINAL
!       end if


       c = ( yomix - yom ) / ( yomix - yoequ ) ! not need to add 'eps' because of the previous 'if' conditions
       c = max ( min ( c , 1.0_dp ) , 0.0_dp )
       if ( c < eps_c .or. c > 1.0_dp - eps_c ) then
!          if (verbose) write (*,*)'MIL: Too high / small progress variable'
          return
       end if


!       if ( taumix_i > inp % ext % rctmax ) then
!!          if (verbose) write (*,*)'MIL: arret pas de sauts'
!          call_mil = .false.
!          return
!       else
!          call_mil = .true.
!       end if

!    else ! AUTOIGNITION possible test if diffusion flame possible

!       call_mil = .true.

!       if ( taumix < ac_crit_MIL / inp % ext % rctmax &
!            .or. yom <= yoequ                         &
!            .or. zm  <= 0.0025_dp                     &
!            .or. zm  >= 0.99_dp                       &
!            .or. yom >= 0.9975_dp * yoinf             &
!            .or. yom <  0.001_dp * yoinf              ) call_mil = .false.

!       if ( yomix > eps8 .and. (yomix-yoequ) > eps8 ) ) then
!          if ( zm >= zst ) then
!             c = abs ( ( yomix - yom ) / ( yomix + eps ) )                                           ! NOT USED !?
!          else
!             c = abs ( ( yomix - yom ) / ( yomix - yoequ + eps ) )
!          end if
!       end if

!       if ( c < 0.005_dp ) call_mil = .false.

!    end if ! .not. call_self


!    !====================================
!    ! Jump positions assessment
!    !====================================

!    if ( call_self .and. call_mil ) then ! MIL possible: AUTOIGNITION jumps comparison

!       call jump_mp ( inp , taumix_i , call_mil , zjm , zjp )
!       call jump_mp_reactivity ( inp , adi , thd , P , T , taumix_i , call_mil , zjm , zjp )
!       if ( .not. call_mil ) return

!       zjm = min ( zignim , zjm , zm )
!       zjp = max ( zignip , zjm , zm )

!       cas = 1

!    else if ( .not. call_self .and. call_mil ) then

!       call jump_mp ( inp , taumix_i , call_mil , zjm , zjp )
       call jump_mp_reactivity ( inp , adi , thd , P , T , taumix_i , call_mil , zjm , zjp )
       if ( .not. call_mil ) return

       zjm = min ( zjm , zm )
       zjp = max ( zjp , zm )

       if ( zjp - zjm < eps_z ) return

!    else if ( call_self .and. .not. call_mil ) then

!       wo2 = wo_auto ( zignim , zignip , zm , max_Z , zp2 , yom , yoinf , zst , zone )               ! TO ADD...
!       wo2 = wo2 * taumix_i

!       cas = 3
!       call_reac = 0

!    else ! ( .not. call_self .and. .not. call_mil )

!       !wo2 = 0.0_dp ! already initialize to zero

!       cas = 0
!       call_reac = 0

!    end if ! ( call_self .and. call_mil )


!    !====================================
!    if ( twophase ) then
!    !====================================

!       ! Two-phase flow injection take into account the transient heating of the droplets
!       ! PAS PRISE EN COMPTE EVAPORATION COMBUSTIBLE (a modifier dans yo_mil et womil)
!       ! seulement cas evaporation lox traite

!       ! change of variables
!       max_Y_zi  = 1.0_dp - ( 1.0_dp - yom ) * min_Z / zm ! = 1 for gas flow because min_Z = 0
!       yfmax_zi  = yom / ( 1.0_dp - zm ) * ( 1.0_dp - max_Z ) ! = 0 for gas flow because max_Z = 1
!       ystmax_zi = 1.0_dp * ( 1.0_dp - min_Z / zm )

!       dz_max_min_i = 1.0_dp / ( max_Z - min_Z )

!       zm  = ( zm - min_Z ) * dz_max_min_i
!       zm  = max ( eps8 , zm )                          ! rescale zm in [0,1]

!       zp2 = zp2 * dz_max_min_i * dz_max_min_i
!       zp2 = min ( zp2 , zm * ( 1.0_dp - zm ) )         ! rescale zp2 in [0,zm(1-zm)]

!       zp2max = zp2max * dz_max_min_i * dz_max_min_i
!       zp2max = min ( zp2max , zm * ( 1.0_dp - zm ) )   ! rescale zp2max in [0,zm(1-zm)]

!       zjm = ( zjm - min_Z ) * dz_max_min_i
!       zjm = max ( esp8 , zjm )                         ! rescale zjm in [0,1]

!       zjp = ( zjp - min_Z ) * dz_max_min_i
!       zjp = min ( 1.0_dp - eps8 , zjp )                ! rescale zjp in [0,1]

!       zst_zi = ( zst - min_Z ) * dz_max_min_i
!       zst_zi = max ( esp8 , zst_zi  )

!    else

       max_Y_zi  = 1.0_dp ! = 1 for gas flow
       yfmax_zi  = 0.0_dp ! = 0 for gas flow
       ystmax_zi = 1.0_dp
       zst_zi    = zst

!    end if ! twophase


    if ( zp2 <= eps_s * zp2max .or. zp2 >= zp2max * ( 1.0_dp - eps_s ) ) then
!       if (verbose) write (*,*) "MIL: too small/high variance: zp2 , zp2max" , zp2 , zp2max
       return
    end if

    !====================================
    ! Beta-PDF shape parameters
    !====================================

    xb = zp2max / zp2 - 1.0_dp
    ! not need to add eps in division because of the previous "if" condition
    if ( xb <= 0.0_dp ) xb = eps8                                             ! ORIGINALY
!    if ( xb <= 0.0_dp ) return                                                  ! (A.T. addition)
    beta1 = zm * xb
    beta2 = ( 1.0_dp - zm ) * xb

    !==============================================================
    ! Burnt gases composition integration from min_Z to max_Z
    !==============================================================

    call yo_mil ( zm , zp2 , 0.0_dp , 1.0_dp , zst_zi , yom , yom , max_Y_zi , yfmax_zi , ystmax_zi , yoequm , yoburnt )
    ! yoburnt = \tilde{Y}^{EQU} = \int_0^1 Y^{EQU}(\xi) \tilde{P}(\xi)d\xi

    if ( yoequm /= yoequm .or. yoburnt /= yoburnt ) then ! (A.T. addition)
       ! Decomment "call yo_mil_debug" to debug the subroutine MIL (A.T.)
!       call yo_mil_debug ( zm , zp2 , 0.0_dp , 1.0_dp , zst_zi , yom , yom , max_Y_zi , yfmax_zi , ystmax_zi , yoequm , yoburnt )
       return
    end if

    if ( yoequm  < yoequ ) yoequm  = yoequ ! happens rarely
    if ( yoburnt < yoequ ) yoburnt = yoequ ! happens rarely


    !===========================================================
    ! PDF between zjm-zjp for possible premixed reaction rate
    !===========================================================

    int_zjmp_pdfbeta = intmg ( zjm , zjp , beta1 , beta2 )
!    yoequm  = 0.0_dp


    !===============================================
    ! MIL composition integration from zjm to zjp
    !===============================================

    call yo_mil ( zm , zp2 , zjm , zjp , zst_zi , yom , yom , max_Y_zi , yfmax_zi , ystmax_zi , yoequm , yomil ) ! yoequm sert à qqch à partir d'ici ? non ?
    ! yomil   = \tilde{Y}^{MIL} ? = \int_0^1 Y^{MIL}(\xi) \tilde{P}(\xi)d\xi ? entre 0 et 1 ?! car je ne suis pas sur que ça soit le cas...

    if ( yomil /= yomil ) return ! (A.T. addition)

    !======================
    ! Model conservation
    !======================

    if ( yomil > yom ) then

!       if (verbose) write (*,*) "MIL: yomil > yom"
       prob_mil = ( yom - yoburnt ) / ( yomil - yoburnt + eps8 )
       prob_mil = min ( max ( prob_mil , 0.0_dp ) , 1.0_dp )

       wo2bur = wo_mil ( zm , beta1 , beta2 , 0.0_dp , 1.0_dp , yom , zst_zi , ystmax_zi )
       wo2bur = wo2bur * taumix_i

       cas = 21

    else

!       if (verbose) write (*,*) "MIL: yomil < yom"
       prob_mil = ( yom - yomix ) / ( yomil - yomix + eps8 )
       prob_mil = min ( max ( prob_mil , 0.0_dp ) , 1.0_dp )

       wo2bur = 0.0_dp

       cas=22

    end if

    wo2mil = wo_mil ( zm , beta1 , beta2 , zjm , zjp , yom , zst_zi , ystmax_zi )
    wo2mil = wo2mil * taumix_i

    wo2sauts = wo_saut ( zjm , zm , zp2 , yom , max_Y_zi , zst_zi , yfmax_zi , ystmax_zi ) &
             + wo_saut ( zjp , zm , zp2 , yom , max_Y_zi , zst_zi , yfmax_zi , ystmax_zi )
    wo2sauts = wo2sauts * taumix_i

    wo2 = prob_mil * ( wo2mil + wo2sauts ) + ( 1.0_dp - prob_mil ) * wo2bur        ! OPTIM : Supprimer les "* taumix_i" pour ne le garder qu'à la fin.


    ! physical limitation of production rate

!    wo2=max(wo2 , -(yom-yoequm)/dt_phys)
!    print* , "yom" , yom
!    print* , "yoequm" , yoequm

    ! back to the physical space, no change of variables to take into account the lower bound of the mixture fraction

!    wo2 = min ( wo2 , 0.0_dp ) ! wo2 <= 0
    if ( wo2 >= - eps2 .or. wo2 /= wo2 ) return

    ! Save jump positions, if autoignition zignim = zignim (because return)
    zignim = zjm
    zignip = zjp

    wreac ( inp % index_O2 ) = wo2

    wreac ( inp % index_H2 ) = wo2 * thd % Wc (inp % index_H2) * thd % Wc_i (inp % index_O2)
    wreac ( inp % index_H2 ) = wreac ( inp % index_H2 ) + wreac ( inp % index_H2 )

    wreac ( inp % index_H2O ) = - wo2 * thd % Wc (inp % index_H2O) * thd % Wc_i (inp % index_O2)
    wreac ( inp % index_H2O ) = wreac ( inp % index_H2O ) + wreac ( inp % index_H2O )

    call_reac = 1

!    if (verbose) then
!       write (*,*) "progress c = " , c
!       write (*,*) "call_reac = " , call_reac"
!       write (*,*) "zjm , zjp = " , zjm , zjp
!       write (*,*) "yoequm = " , yoequm
!       write (*,*) "yoburnt = " , yoburnt
!       write (*,*) "prob_mil = " , prob_mil
!       write (*,*) "wo2mil , wo2sauts , wo2bur = " , wo2mil , wo2sauts , wo2bur
!       write (*,*) "wreac(:) = " , wreac (:)
!    end if



  end subroutine wreac_MIL


!> \brief Calculate the mixture fraction left and right jump positions
!!
!! Assess by interpolation the left and right jump positions from the inversed mixing
!! time scale 'taumix_i' and the table of the reversed chemical time scale.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine jump_mp ( inp , taumix_i , call_mil , zjm , zjp )


    type (inp_type) , intent (in)              :: inp      !< input derived type
    real (dp) , intent (in)                    :: taumix_i !< inversed mixing time scale (1/s)
    logical   , intent (inout)                 :: call_mil !< call MIL ?
    real (dp) , intent (inout)                 :: zjm      !< mixture fraction left  jump position
    real (dp) , intent (inout)                 :: zjp      !< mixture fraction right jump position

    integer (ip)                               :: zj_index , iminus , iplus , i
    real (dp)                                  :: xorigin , yorigin , pente_i


    if ( taumix_i > inp % ext % rctmax ) then
       call_mil = .false.
       return
    end if


    !=================================================
    ! Mixture fraction left jump position assessment
    !=================================================

    iminus = 1                         ! left (minus) index
    iplus  = inp % ext % index_rctmax  ! right (plus) index

    do

       ! sweep from left to right
       zj_index = int(iminus + (iplus-iminus) * 0.5_dp)

       if ( inp % ext % table (2,zj_index) < taumix_i ) then
          iminus = zj_index ! move forward
       else
          iplus  = zj_index ! move backward
       end if

       if ( iplus-iminus == 1 ) exit

    end do


    ! linear interpolation zjm
    i = 1

    if ( 1 < zj_index .and. zj_index <= inp % ext % index_rctmax ) then

       if ( taumix_i > inp % ext % table (2,zj_index) ) then
          i = zj_index
       else
          i = zj_index - 1
       end if

    end if

    xorigin = inp % ext % table (1,i)
    yorigin = inp % ext % table (2,i)

    pente_i = ( inp % ext % table (1,i+1) - xorigin ) / ( inp % ext % table (2,i+1) - yorigin )
    zjm     = ( taumix_i - yorigin ) * pente_i + xorigin


    !=================================================
    ! Mixture fraction right jump position assessment
    !=================================================

    iminus = inp % ext % index_rctmax  ! left (minus) index
    iplus  = inp % ext % nval          ! right (plus) index

    do

       ! sweep from right to left
       zj_index = int(iminus + (iplus-iminus) * 0.5_dp)

       if ( inp % ext % table (2,zj_index) < taumix_i ) then
          iplus  = zj_index ! move backward
       else
          iminus = zj_index ! move forward
       end if

       if ( iplus-iminus == 1 ) exit

    end do


    ! linear interpolation zjp
    i = inp % ext % nval

    if ( inp % ext % index_rctmax <= zj_index .and. zj_index < inp % ext % nval ) then

       if ( taumix_i > inp % ext % table (2,zj_index) ) then
          i = zj_index
       else
          i = zj_index + 1
       end if

    end if

    xorigin = inp % ext % table (1,i)
    yorigin = inp % ext % table (2,i)

    pente_i = ( xorigin - inp % ext % table (1,i-1) ) / ( yorigin - inp % ext % table (2,i-1) )
    zjp     = ( taumix_i - yorigin ) * pente_i + xorigin


   end subroutine jump_mp


!> \brief Calculate the mixture fraction left and right jump positions
!!
!! Assess the left and right jump positions from the inversed mixing time scale 'taumix_i' 
!! and the reactivity rate of Boivin, which is inversely proportional to an 
!! auto-ignition time.
!!
!! References:
!! -# Bridel-Bertomeu, Boivin, Explicit chemical timescale as a substitute 
!!    for tabultad chemistry in a H2-O2 turbulent flame simulation, 
!!    Combustion Science and Technology (2015).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine jump_mp_reactivity ( inp , adi , thd , P , T , taumix_i , call_mil , zjm , zjp )


    type (inp_type) , intent (in)      :: inp      !< input derived type
    type (adi_type) , intent (in)      :: adi      !< non-dimensional derived type
    type (thd_type) , intent (in)      :: thd      !< thermodynamic derived type
    real (dp)       , intent (in)      :: P        !< pressure (Pa)
    real (dp)       , intent (in)      :: T        !< temperature (K)
    real (dp)       , intent (in)      :: taumix_i !< inversed mixing time scale (1/s)
    logical         , intent (inout)   :: call_mil !< call MIL ?
    real (dp)       , intent (inout)   :: zjm      !< mixture fraction left  jump position
    real (dp)       , intent (inout)   :: zjp      !< mixture fraction right jump position

    integer (ip)                               :: i , l 
    integer (ip) , parameter                   :: npt = 200 ! nb points in the table
    real (dp)                                  :: rho , dim_rho , W_i , wrk_z

    real (dp) , dimension (1:nrv)              :: Ya
    real (dp) , parameter                      :: z_min = 1.0e-5_dp , &  ! min and max mixture_fraction
                                                  z_max = 0.99_dp
    real (dp)                                  :: zmr , lambda_max ! most reactive mixture fraction and max reactivity rate

    ! Dichotomy method (with Newton method, it is longest !!! cheched by CPU time - A.Techer)
    integer (ip)                               :: i_zmr , zj_index , iminus , iplus
    real (dp)                                  :: xorigin , yorigin , pente_i
    real (dp) , dimension (1:npt)              :: z , lambda


    wrk_z = log ( z_max / z_min ) / float ( npt - 1 )
    dim_rho = adi % W_ref / adi % R_ref


    ! Build table
    !=================================================
!    open ( 13 , file = 'lambda_xi.dat' )
    do i = 1 , npt

       z (i) = z_min * exp ( float ( i- 1 ) * wrk_z )
       do l = 1 , nrv
          Ya (l) = z (i) * inp % Y1f (l) + ( 1.0_dp - z (i) ) * inp % Y0o (l)
       end do
       call Wmix_i_scalar ( thd , Ya , W_i )
       rho = P / ( T * W_i ) * dim_rho

       call reactivity_0D ( inp , adi , thd , rho , T , Ya , lambda (i) )
!       write (13,*) z(i) , lambda (i)

    end do
!    close (13)


    ! Maximum reactivity
    !=================================================
    i_zmr = maxloc ( lambda (:) , dim=1 ) ! index most reactive mixture fraction
    zmr        = z      (i_zmr)
    lambda_max = lambda (i_zmr)
!    write (*,'(A,I5,5(1PE15.5))') ' i_zmr, zmr, lambda_max = ' , i_zmr , zmr , lambda_max


    ! Continued call MIL ?
    !=================================================
    call_mil = .true.
    if ( taumix_i > lambda_max ) then
       call_mil = .false.
!       write (*,*) 'pas de saut : return'
       return
    end if
!    write (*,*) call_mil


    ! Mixture fraction left jump position assessment
    !=================================================
    iminus = 1      ! left (minus) index
    iplus  = i_zmr  ! right (plus) index
!    nit = 0
    do
!       nit = nit + 1
       ! sweep from left to right
       zj_index = int(iminus + (iplus-iminus)*0.5_dp)
       if ( lambda (zj_index) < taumix_i ) then
!          write (*,*) "move forward"
          iminus = zj_index ! move forward
       else
!          write (*,*) "move backward"
          iplus  = zj_index ! move backward
       end if
       if ( iplus-iminus == 1 ) exit
    end do
!    write (*,*) "nit left jump = " , nit

    ! linear interpolation zjm
    i = 1
    if ( 1 < zj_index .and. zj_index <= i_zmr ) then
       if ( taumix_i > lambda (zj_index) ) then
          i = zj_index
       else
          i = zj_index - 1
       end if
    end if

    xorigin = z (i)
    yorigin = lambda (i)

    pente_i = ( z (i+1) - xorigin ) / ( lambda (i+1) - yorigin )
    zjm     = ( taumix_i - yorigin ) * pente_i + xorigin
    zjm     = max ( zjm , min_Z )


    ! Mixture fraction right jump position assessment
    !=================================================
    iminus = i_zmr  ! left (minus) index
    iplus  = npt    ! right (plus) index
!    nit = 0
    do
!       nit = nit + 1
       ! sweep from right to left
       zj_index = int(iminus + (iplus-iminus) * 0.5_dp)
       if ( lambda (zj_index) < taumix_i ) then
!          write (*,*) "move backward"
          iplus  = zj_index ! move backward
       else
!          write (*,*) "move forward"
          iminus = zj_index ! move forward
       end if
       if ( iplus-iminus == 1 ) exit
    end do
!    write (*,*) "nit right jump = " , nit

    ! linear interpolation zjp
    i = npt
    if ( i_zmr <= zj_index .and. zj_index < npt ) then
       if ( taumix_i > lambda (zj_index) ) then
          i = zj_index
       else
          i = zj_index + 1
       end if
    end if

    xorigin = z (i)
    yorigin = lambda (i)

    pente_i = ( xorigin - z (i-1) ) / ( yorigin - lambda (i-1) )
    zjp     = ( taumix_i - yorigin ) * pente_i + xorigin
    zjp     = min ( zjp , max_Z )


!    write (*,*) "zjm , zjp = " , zjm , zjp


   end subroutine jump_mp_reactivity


!> \brief Calculate the scalar mean values in composition space.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine yo_mil ( zm , zp2 , zjm , zjp , zstoic , yom , yo2sub , yomax , yfmax , ystmax , &
                      yoequm , yo2mil )


    real (dp) , intent (in)                    :: zm     !< mixture fraction
    real (dp) , intent (in)                    :: zp2    !< SGS variance
    real (dp) , intent (in)                    :: zjm    !< mixture fraction left  jump position
    real (dp) , intent (in)                    :: zjp    !< mixture fraction right jump position
    real (dp) , intent (in)                    :: zstoic !< stoichiometric mixture fraction (modified if twophase)
    real (dp) , intent (in)                    :: yom    !< mean oxidizer mass fraction
    real (dp) , intent (in)                    :: yo2sub !< ???
    real (dp) , intent (in)                    :: yomax  !< maximum oxidizer mass fraction (modified if twophase)
    real (dp) , intent (in)                    :: yfmax  !< ??? Yo max in fuel stream ?
    real (dp) , intent (in)                    :: ystmax !< ??? Yo stoechiometric max = 1 (onephase) ?
    real (dp) , intent (inout)                 :: yoequm !< mean value on equilibrium trajectory
    real (dp) , intent (inout)                 :: yo2mil !< mean value on MIL trajectory

    real (dp) :: aeq , beqm , beqp , yz1 , yz2 , zstoic_i
    real (dp) :: int1 , int2 , int3 , ycalc , ymeq
    real (dp) :: xb , beta1 , beta2


    ymeq = 0.0_dp

    xb = ( zm * ( 1.0_dp - zm ) ) / zp2 - 1.0_dp
    if ( xb <= 0.0_dp ) xb = eps8
    beta1 = zm * xb
    beta2 = ( 1.0_dp - zm ) * xb

    if ( ( max_Z - zm ) <= eps8 ) then

       ycalc  = yom
       yo2mil = ycalc
       int1   = 0.0_dp
       int2   = 0.0_dp
       int3   = 0.0_dp

    else

       zstoic_i = 1.0_dp / zstoic
       aeq  = ystmax - ystmax * zm * zstoic_i
       beqm = ystmax * zstoic_i
       beqp = - ystmax * zstoic_i
       yz1  = (yomax-yo2sub) / zm
       yz2  = ( - yo2sub + yfmax ) / (1.0_dp-zm)

       if ( zstoic < zjm ) then

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = 0.0_dp

       else if ( zstoic < zm ) then

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zstoic , beta1 , beta2 ) &
               + beqm * intzmg ( zjm , zstoic , zm , beta1 , beta2 )

       else if ( zstoic < zjp ) then

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zm , beta1 , beta2 ) &
               + beqm * intzmg ( zjm , zm , zm , beta1 , beta2 ) &
               + aeq * intpg ( zm , zstoic , beta1 , beta2 ) &
               + beqp * intzpg ( zm , zstoic , zm , beta1 , beta2 )

       else if ( zstoic < max_Z ) then

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zm , beta1 , beta2 ) &
               + beqm * intzmg ( zjm , zm , zm , beta1 , beta2 ) &
               + aeq * intpg ( zm , zjp , beta1 , beta2 ) &
               + beqp * intzpg ( zm , zjp , zm , beta1 , beta2 )

       else

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zm , beta1 , beta2 ) &
               + beqm * intzmg ( zjm , zm , zm , beta1 , beta2 ) &
               + aeq * intpg ( zm , zjp , beta1 , beta2 ) &
               + beqp * intzpg ( zm , zjp , zm , beta1 , beta2 )

       end if

       ycalc  = yo2sub * int1 + yz1 * int2 + yz2 * int3 + ymeq
       yo2mil = ycalc

    end if

    yoequm = ymeq


  end subroutine yo_mil



!> \brief Calculate the scalar mean values in composition space.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine yo_mil_debug ( zm , zp2 , zjm , zjp , zstoic , yom , yo2sub , yomax , yfmax , ystmax , &
                            yoequm , yo2mil )


    real (dp) , intent (in)                    :: zm     !< mixture fraction
    real (dp) , intent (in)                    :: zp2    !< SGS variance
    real (dp) , intent (in)                    :: zjm    !< mixture fraction left  jump position
    real (dp) , intent (in)                    :: zjp    !< mixture fraction right jump position
    real (dp) , intent (in)                    :: zstoic !< stoichiometric mixture fraction (modified if twophase)
    real (dp) , intent (in)                    :: yom    !< mean oxidizer mass fraction
    real (dp) , intent (in)                    :: yo2sub !< ???
    real (dp) , intent (in)                    :: yomax  !< maximum oxidizer mass fraction (modified if twophase)
    real (dp) , intent (in)                    :: yfmax  !< ??? Yo max in fuel stream ?
    real (dp) , intent (in)                    :: ystmax !< ??? Yo stoechiometric max = 1 (onephase) ?
    real (dp) , intent (inout)                 :: yoequm !< mean value on equilibrium trajectory
    real (dp) , intent (inout)                 :: yo2mil !< mean value on MIL trajectory

    real (dp) :: aeq , beqm , beqp , yz1 , yz2 , zstoic_i
    real (dp) :: int1 , int2 , int3 , ycalc , ymeq
    real (dp) :: xb , beta1 , beta2


    ymeq = 0.0_dp
write (*,*) "zm , zp2 = " , zm , zp2
    xb = ( zm * ( 1.0_dp - zm ) ) / zp2 - 1.0_dp
write (*,*) "xb = " , xb
    if ( xb <= 0.0_dp ) xb = eps8
    beta1 = zm * xb
    beta2 = ( 1.0_dp - zm ) * xb
write (*,*) "beta1 , beta2 = " , beta1 , beta2

    if ( ( max_Z - zm ) <= eps8 ) then

       ycalc  = yom
       yo2mil = ycalc
       int1   = 0.0_dp
       int2   = 0.0_dp
       int3   = 0.0_dp

    else

       zstoic_i = 1.0_dp / zstoic
       aeq  = ystmax - ystmax * zm * zstoic_i
       beqm = ystmax * zstoic_i
       beqp = - ystmax * zstoic_i
write (*,*) "aeq , beqm , beqp = " , aeq , beqm , beqp
       yz1  = (yomax-yo2sub) / zm
       yz2  = ( - yo2sub + yfmax ) / (1.0_dp-zm)
write (*,*) "yz1 , yz2 = " , yz1 , yz2

       if ( zstoic < zjm ) then

write (*,*) "case: zstoic < zjm"
          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = 0.0_dp

       else if ( zstoic < zm ) then

write (*,*) "case: zstoic < zm"
          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zstoic , beta1 , beta2 ) &
               + beqm * intzmg ( zjm , zstoic , zm , beta1 , beta2 )

       else if ( zstoic < zjp ) then

write (*,*) "case: zstoic < zjp"
          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zm , beta1 , beta2 ) &
               + beqm * intzmg ( zjm , zm , zm , beta1 , beta2 ) &
               + aeq * intpg ( zm , zstoic , beta1 , beta2 ) &
               + beqp * intzpg ( zm , zstoic , zm , beta1 , beta2 )

       else if ( zstoic < max_Z ) then

write (*,*) "case: zstoic < max_Z"
          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zm , beta1 , beta2 ) &
               + beqm * intzmg ( zjm , zm , zm , beta1 , beta2 ) &
               + aeq * intpg ( zm , zjp , beta1 , beta2 ) &
               + beqp * intzpg ( zm , zjp , zm , beta1 , beta2 )

       else

write (*,*) "case: else"
          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 ) &
               + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zm , beta1 , beta2 ) &
               + beqm * intzmg ( zjm , zm , zm , beta1 , beta2 ) &
               + aeq * intpg ( zm , zjp , beta1 , beta2 ) &
               + beqp * intzpg ( zm , zjp , zm , beta1 , beta2 )

       end if

write (*,*) "int1 , int2 , int3 , ymeq = " , int1 , int2 , int3 , ymeq
       ycalc  = yo2sub * int1 + yz1 * int2 + yz2 * int3 + ymeq
       yo2mil = ycalc

    end if

    yoequm = ymeq
write (*,*) "ycalc , yo2mil , yoequm = " , ycalc , yo2mil , yoequm
stop


  end subroutine yo_mil_debug

!> \brief Calculate ....
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function wo_mil ( zm , beta1 , beta2 , zjm , zjp , yostar , zst , ystmax )


    real (dp) , intent (in)   :: zm     !< 
    real (dp) , intent (in)   :: beta1  !< 
    real (dp) , intent (in)   :: beta2  !< 
    real (dp) , intent (in)   :: zjm    !< 
    real (dp) , intent (in)   :: zjp    !< 
    real (dp) , intent (in)   :: yostar !< 
    real (dp) , intent (in)   :: zst    !< 
    real (dp) , intent (in)   :: ystmax !< 
    real (dp)                 :: wo_mil !< 

    real(dp)                  :: aeq , wm , wp


    if ( ( 1.0_dp - zm ) <= eps8 ) then

       wo_mil = 0.0_dp

    else

       aeq = ystmax - ystmax * zm / zst
       wm  = aeq - yostar
       wp  = - yostar

       if ( zst < zjm ) then
          wo_mil = wp * intmg ( zjm , zm , beta1 , beta2 ) &
                 + wp * intpg ( zm , zjp , beta1 , beta2 )
       else if ( zst < zm ) then
          wo_mil = wm * intmg ( zjm , zst , beta1 , beta2 ) &
                 + wp * intmg ( zst , zm , beta1 , beta2 ) &
                 + wp * intpg ( zm , zjp , beta1 , beta2 )
       else if ( zst < zjp ) then
          wo_mil = wm * intmg ( zjm , zm , beta1 , beta2 ) &
                 + wm * intpg ( zm , zst , beta1 , beta2 ) &
                 + wp * intpg ( zst , zjp , beta1 , beta2 )
       else if ( zst < max_Z ) then  
          wo_mil = wm * intmg ( zjm , zm , beta1 , beta2 ) &
                 + wm * intpg ( zm , zjp , beta1 , beta2 )
       else
          wo_mil = wm * intmg ( zjm , zm , beta1 , beta2 ) &
                 + wm * intpg ( zm , zjp , beta1 , beta2 )
       end if

       wo_mil = min ( wo_mil , 0.0_dp )

    end if


  end function wo_mil


!> \brief Calculate ....
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function wo_saut ( saut , zm , zp2 , yo2sub , yomax , zst , yfmax , ystmax )


    real (dp) , intent (in)   :: saut    !< 
    real (dp) , intent (in)   :: zm      !< 
    real (dp) , intent (in)   :: zp2     !< 
    real (dp) , intent (in)   :: yo2sub  !< 
    real (dp) , intent (in)   :: yomax   !< 
    real (dp) , intent (in)   :: zst     !< 
    real (dp) , intent (in)   :: yfmax   !< 
    real (dp) , intent (in)   :: ystmax  !< 
    real (dp)                 :: wo_saut !< 

    real(dp)                  :: psaut , zp2max , xb , a , b
    real(dp)                  :: aeq , beqm , yz1 , yz2 , yo2iem , yo2equ
    real(dp)                  :: zst_i , zm_i , wrk


    psaut  = 0.0_dp
    zp2max = zm * ( 1.0_dp - zm )

    if ( zp2 > ( 0.01_dp * zp2max ) .and. ( zp2 < 0.99_dp * zp2max ) ) then
       xb = zp2max / zp2 - 1.0_dp
       if ( xb <= 0.0_dp ) xb = eps8
       a = zm * xb
       b = ( 1.0_dp - zm ) * xb
       call pdfbeta ( saut , a , b , psaut )
    end if


    if ( ( max_Z - zm ) <= eps8 ) then

       wo_saut = 0.0_dp

    else

       zm_i = 1.0_dp / zm
       zst_i = 1.0_dp / zst
       wrk = 1.0_dp / ( 1.0_dp - zm )
       aeq  = yomax - yomax * zm * zst_i
       beqm = yomax * zst_i
       yz1  = ( yomax - yo2sub ) * zm_i
       yz2  = - yo2sub * wrk

       if ( saut < zm ) then
         yo2iem = ( yo2sub - yomax ) * zm_i * saut + yomax
       else
         yo2iem = wrk * ( yfmax * ( saut - zm ) + yo2sub * ( 1.0_dp - saut ) )
       end if

       if ( saut < zst ) then
         yo2equ = ystmax * ( 1.0_dp - saut * zst_i )
       else
         yo2equ = 0.0_dp
       end if

       if ( zm >= saut ) then
         wo_saut = ( yo2equ - yo2iem ) * ( zm - saut ) * psaut
       else
         wo_saut = ( yo2equ - yo2iem ) * ( saut - zm ) * psaut
       end if

       wo_saut = min ( wo_saut , 0.0_dp )

    end if


  end function wo_saut


!> \brief Calculate ....
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function intmg ( a , b , beta1 , beta2 )

    use nr_subroutines , only : betai

    real (dp) , intent (in)   :: a     !< 
    real (dp) , intent (in)   :: b     !< 
    real (dp) , intent (in)   :: beta1 !< 
    real (dp) , intent (in)   :: beta2 !< 
    real (dp)                 :: intmg !< 

    intmg = betai ( beta1 , beta2 , b ) - betai ( beta1 , beta2 , a )

  end function intmg


!> \brief Calculate ....
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function intpg ( a , b , beta1 , beta2 )

    use nr_subroutines , only : betai

    real (dp) , intent (in)   :: a     !< 
    real (dp) , intent (in)   :: b     !< 
    real (dp) , intent (in)   :: beta1 !< 
    real (dp) , intent (in)   :: beta2 !< 
    real (dp)                 :: intpg !< 

    intpg = betai ( beta1 , beta2 , b ) - betai ( beta1 , beta2 , a )

  end function intpg


!> \brief Calculate ....
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function intzmg ( a , b , zm , beta1 , beta2 )

    use nr_subroutines , only : betai , beta_s

    real (dp) , intent (in)   :: a      !< 
    real (dp) , intent (in)   :: b      !< 
    real (dp) , intent (in)   :: zm     !< 
    real (dp) , intent (in)   :: beta1  !< 
    real (dp) , intent (in)   :: beta2  !< 
    real (dp)                 :: intzmg !< 

    intzmg = zm * ( betai ( beta1 , beta2 , b ) - betai ( beta1 , beta2 , a ) )
    intzmg = intzmg - beta_s ( beta1 + 1.0_dp , beta2 ) / beta_s ( beta1 , beta2 ) &
           * ( betai ( beta1 + 1.0_dp , beta2 , b ) - betai ( beta1 + 1.0_dp , beta2 , a ) )

  end function intzmg


!> \brief Calculate ....
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function intzpg ( a , b , zm , beta1 , beta2 )

    use nr_subroutines , only : betai , beta_s

    real (dp) , intent (in)   :: a      !< 
    real (dp) , intent (in)   :: b      !< 
    real (dp) , intent (in)   :: zm     !< 
    real (dp) , intent (in)   :: beta1  !< 
    real (dp) , intent (in)   :: beta2  !< 
    real (dp)                 :: intzpg !< 

    intzpg = - zm * ( betai ( beta1 , beta2 , b ) - betai ( beta1 , beta2 , a ) )
    intzpg = intzpg + beta_s ( beta1 + 1.0_dp , beta2 ) / beta_s ( beta1 , beta2 ) &
           * ( betai ( beta1 + 1.0_dp , beta2 , b ) - betai ( beta1 + 1.0_dp , beta2 , a ) )

  end function intzpg


!> \brief Calculate ....
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function intz ( a , b , beta1 , beta2 )

    use nr_subroutines , only : betai , beta_s

    real (dp) , intent (in)   :: a      !< 
    real (dp) , intent (in)   :: b      !< 
    real (dp) , intent (in)   :: beta1  !< 
    real (dp) , intent (in)   :: beta2  !< 
    real (dp)                 :: intz   !< 

    intz = beta_s ( beta1 + 1.0_dp , beta2 ) / beta_s ( beta1 , beta2 ) &
         * ( betai ( beta1 + 1.0_dp , beta2 , b ) - betai ( beta1 + 1.0_dp , beta2 , a ) )

  end function intz


!> \brief Calculate ....
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function intzsq ( a , b , beta1 , beta2 )

    use nr_subroutines , only : betai , beta_s

    real (dp) , intent (in)   :: a      !< 
    real (dp) , intent (in)   :: b      !< 
    real (dp) , intent (in)   :: beta1  !< 
    real (dp) , intent (in)   :: beta2  !< 
    real (dp)                 :: intzsq !< 

    intzsq = beta_s ( beta1 + 2.0_dp , beta2 ) / beta_s ( beta1 , beta2 ) &
           * ( betai ( beta1 + 2.0_dp , beta2 , b ) - betai ( beta1 + 2.0_dp , beta2 , a ) )

  end function intzsq


!> \brief Load table values of the reverse chemical time scale in function of the mixture fraction
!!
!! Load table values from the file file_tabchem = 'tabrtchem.dat' which looks like :
!!
!!   #   mixture_fraction   1/tau_chem
!!     0.000000000000000E+00     3.710000490256525E-03
!!     1.010101010101010E-02     9.540475481540168E+00
!!     .....................     .....................
!!     1.121115708232785E-02     1.871266080573139E+01
!!
!! so, it contains 2 columns (C1=mixture_fraction (-), C2=1/tau_chem (s^-1).
!! This table is used by the MIL non-premixed combustion model.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine load_table (inp)


    type (inp_type) , intent (inout)                               :: inp   !< input derived type

    real (dp)                     :: dummy
    integer (ip)                  :: i , ok


    open ( unit_tabchem , file = file_tabchem , status = 'old' , action = 'read' , iostat = ok )
    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (file_tabchem))

       read (unit_tabchem,*) ! file header


       ! count the number of values in the tab
       inp % ext % nval = 0
       do
          read ( unit_tabchem , * , iostat=ok ) dummy
          if ( ok /=0 ) exit
          inp % ext % nval = inp % ext % nval + 1
       end do

       rewind (unit_tabchem) ! go back up to the file

       read (unit_tabchem,*) ! file header

       allocate  ( inp % ext % table ( 2 , inp % ext % nval )           , &
                   stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate load_table')

       do i = 1 , inp % ext % nval
          read (unit_tabchem,*,iostat=ok) inp % ext % table (:,i)
          if ( ok /=0 ) call abort_mpi ('error: reading file ' // trim (file_tabchem))
       end do


    close (unit_tabchem)


    inp % ext % index_rctmax = maxloc ( inp % ext % table (2,:) , dim=1 )
    inp % ext % z_rctmax     = inp % ext % table (1,inp % ext % index_rctmax)
    inp % ext % rctmax       = inp % ext % table (2,inp % ext % index_rctmax)


  end subroutine load_table


!> \brief Calculate the reactivity [Boivin et al., 2011].
!!
!! Assess the reactivity rate \lambda (1/s) defined by Boivin et al. for this
!!  3-Step Reduced Mechanism for H2-O2 Systems
!!  5(+1) species: H2  O2  H2O  H  HO2  (N2 diluent)
!!  2(+1) SPECIES: O  OH  (H2O2 not used)
!!                 3H2 + O2   = 2H2O + 2H    (R1)
!!                 H + H + M  = H2 + M       (R2)
!!                 H2 + O2    = HO2 + H      (R3)
!!
!! References:
!! -# P. Boivin, C. Jiménez, A.L. Sánchez, F.A. Williams, 
!!    An explicit reduced mechanism for H2-air combustion,
!!    Proceedings of the Combustion Institute, 2011, 33 (1), 517-523
!! -# P. Boivin, A. Dauptain, C. Jiménez, B. Cuenot
!!    Simulation of a supersonic hydrogen–air autoignition-stabilized flame using reduced chemistry,
!!    Combustion and Flame, 2012, 159 (4), 1779-1790
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine reactivity_0D ( inp , adi , thd , rho , T , Ya , lambda )


    type (inp_type) , intent (in)                  :: inp    !< input derived type
    type (adi_type) , intent (in)                  :: adi    !< non-dimensional derived type
    type (thd_type) , intent (in)                  :: thd    !< thermodynamic derived type
    real (dp)                     , intent (in)    :: rho    !< density (kg/m3)
    real (dp)                     , intent (in)    :: T      !< temperature (K)
    real (dp) , dimension (1:nrv) , intent (in)    :: Ya     !< mass fractions array (g/g)
    real (dp)                     , intent (inout) :: lambda !< reactivity rate (1/s)

    real (dp) , parameter :: small = 1.0e-50_dp , &
                             RU  = 8.314510e7_dp , & ! ideal gas constant in (g*cm^2/s^-2)/mol/K
                             RUc = 1.987_dp          ! ideal gas constant in cal/mol/K

    real (dp) , dimension (1:4)   :: RF
    real (dp) , dimension (1:nrv) :: Xcon ! Mass fraction and Molar concentration

    integer (ip)  :: l
    real (dp)     :: W_cgs_i
    real (dp)     :: rho_cgs , rho_phys , T_phys
    real (dp)     :: CH2 , CO2 , CH2O , CMb , sumC ! Molar concentration
    real (dp)     :: PR , PRLOG , FC , FCLOG , X , F , PCOR , CPRLOG , K0 , Kinf ! TROE Falloff
    real (dp)     :: WF1_CH , WF2_CO , WF3_COH , WF4_CH
    real (dp)     :: lnT , RTR
    real (dp)     :: A0 , A1 , A2 , AA


    rho_cgs  = 1.0e-3_dp
    rho_phys = rho * rho_cgs ! in g/cm^3
    W_cgs_i  = 1.0e-3_dp / adi % W_ref
    T_phys   = T ! in K


    sumC = 0.0_dp
    do l = 1 , nrv
       Xcon (l) = rho_phys * Ya(l) * thd % Wc_i (l) * W_cgs_i ! = rho*Y_k/W_k = molar concentration (mol/cm^-3)
       sumC = sumC + Xcon (l)
       Xcon (l) = max ( Xcon (l) , small )
    end do


    ! SET LOCAL VALUES FOR THE MOLAR CONCENTRATIONS
    ! H2  O2  H2O  H  HO2  N2

    CH2  =  Xcon (inp % index_H2)
    CO2  =  Xcon (inp % index_O2)
    CH2O =  Xcon (inp % index_H2O)

    ! EFFECTIVE THIRD BODY FOR ALL REACTIONS
    !   Chaperon efficiencies H2/2.5 /H2O/16.0

    CMb = sumC + 1.5_dp*CH2 + 15.0_dp*CH2O

    lnT = log(T_phys)                 ! ln(T)
    RTR = 1.0e3_dp / ( RUc * T_phys ) ! 10^3 / (cal/mol) ==> 10^3 because Ea is given in kcal

    ! SET THE ELEMENTARY RATES
    !  RF = k_forward = A_k * exp ( n_k * ln(T) - Ea / RT )
    !  with A_k  = Frequency factor in cm, s, K et mol
    !       n_k  = Temperature exponent (-)
    !       Ea_k = Activation energy in kcal/mol --> ATTENTION in KILOcal !

    RF(1) = 3.52e16_dp * exp (  -0.7_dp*lnT - 17.07_dp*RTR )
    RF(2) = 5.06e04_dp * exp (  2.67_dp*lnT -  6.29_dp*RTR )
    RF(3) = 1.17e09_dp * exp (  1.30_dp*lnT -  3.64_dp*RTR )

    ! SET THE ELEMENTARY RATES INCLUDING THIRD BODIES: TROE FALLOFF CORRECTION

    K0    = 5.75e19_dp * exp (  -1.4_dp*lnT                )
    Kinf  = 4.65e12_dp * exp (  0.44_dp*lnT                )

    PR = K0*CMb / Kinf
    PRLOG = DLOG10(MAX(PR,SMALL))
    FC = 0.5_dp
    FCLOG = DLOG10(MAX(FC,SMALL))
    CPRLOG = PRLOG - (0.4_dp + 0.67_dp*FCLOG)
    X = CPRLOG/(0.75_dp-1.27_dp*FCLOG-0.14_dp*CPRLOG)
    F = 10.0_dp**( FCLOG/(1.0_dp+X*X) )
    PCOR = PR*F/(1.0_dp+PR)
    RF(4) = Kinf*PCOR

    WF1_CH = RF(1)*CO2
    WF2_CO = RF(2)*CH2
    WF3_COH = RF(3)*CH2
    WF4_CH = RF(4)*CO2

    ! Reactivity in CF(2012)eq.5 = PCI(2011)eq.8 + 4th reaction
!    WF4_CH = 0.0_dp                                            ! Uncomment this line to assess the reactivity defined in PCI(2011)eq.8

    call SLAMBDA_COEFF (WF1_CH, WF2_CO, WF3_COH, WF4_CH, &
                        A0, A1, A2, AA)
     ! A2 /= 0 normally, it's zero if you have CO2 .AND CH2 = 0
     ! AA CAN BE < 0, if T<Tc. Cf. CF_2013 for correction

    call CK_SLAMBDA (AA < 0.0_dp, WF1_CH, WF2_CO, WF3_COH, WF4_CH, &
                     A0, A1, A2, AA, lambda )
    call CK_SLAMBDA (lambda < 0.0_dp, WF1_CH, WF2_CO, WF3_COH, WF4_CH, &
                     A0, A1, A2, AA, lambda )

    if ( lambda < 0.0_dp ) &
       write(*,*) "ERROR: in reactivity calculating lambda negative ", lambda
    lambda = max ( 0.0_dp , lambda )


  end subroutine reactivity_0D


end module solver_reaction
