!------------------------------------------------------------------------------
! MODULE: solver_reaction
!------------------------------------------------------------------------------
!> \brief Solver for the chemical production rate in species equations.
!!
!! This module contains the subroutines and utilities needed to
!! solve the chemical production rate $\omega_\alpha$ by using:\n
!!   * the vode solver,\n
!!   * the MIL (Model Intermittent Lagrangien).\n
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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

  real (dp)    , parameter , private        :: eps8 = 1.0e-8_dp , eps4 = 1.0e-4_dp , eps2 = 1.0e-2_dp , eps = 1.0e-15_dp


contains


!> \brief Selector for the combustion model/solver
!!
!! Reads the keyword to select the corresponding combustion solver.
!!
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine reaction_selec ( inp , thd , adi , time , dt , x , y , z ,      &
                              dx_i , dy_i , dz_i , mu_SGS , T , W_i , cp , ha , v , DTmax )


    type (inp_type) , intent (inout)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    real (dp) , intent (in)                                        :: time  !< time
    real (dp) , intent (in)                                        :: dt    !< time step
    real (dp) , allocatable , dimension (:) , intent (in)          :: x     !< x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: y     !< y-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: z     !< z-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i  !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i  !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i  !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: mu_SGS!< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i   !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp    !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha    !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v     !< conserved variables array
    real (dp) , intent (inout)                                     :: DTmax !< maximum temperature increase
    integer(ip)                                                    :: ite


    if ( inp % reac_selec == dvode ) then
       call reaction_dvode ( inp , thd , adi , dt , T , v , DTmax )

    else if ( inp % reac_selec == MIL ) then
       call reaction_MIL ( inp , thd , adi , dt , dx_i , dy_i , dz_i , mu_SGS , T , v , DTmax )

    else
       call abort_mpi ('Reaction selector ' // trim (inp % reac_selec) // ' not defined')

    end if

    call upd_prim_var_woT ( thd , 0 , v , T , W_i , cp , ha )
    call comm_cons (v)
    call upd_boundaries ( inp , adi , thd , time , dt , x , y , z ,    &
                          dx_i , dy_i , dz_i , T , W_i , cp , ha , ite , v )
    call upd_prim_var_ghost ( thd , T , W_i , cp , ha , v )


  end subroutine reaction_selec



!> \brief Calculate the reactive system of equations by using VODE solver
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine reaction_dvode ( inp , thd , adi , dt , T , v , DTmax )


    type (inp_type) , intent (inout)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    real (dp) , intent (in)                                        :: dt    !< time step
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v     !< conserved variables array
    real (dp) , intent (inout)                                     :: DTmax !< maximum temperature increase


    real (dp) , parameter         :: eps_w = 1.0e-12_dp , eps_f = 1.0e-6_dp , eps_o = 1.0e-6_dp
    integer (ip)                  :: ix , fx , &
                                     iy , fy , &
                                     iz , fz
    integer (ip)                  :: i , j , k , l
    real (dp)                     :: dt_phys , rho_phys , T_phys , rho_i , rho0 , T0_phys , sumY , Tref_i , rho_cgs
    logical                       :: call_dvode_subdomain , call_dvode , call_dvode_w , call_dvode_f , call_dvode_o
    real (dp) , dimension (1:nrv) :: Ya , DeltaY , Yf , Yo


    ! New criteria that seems to work well with the presence or not of radicals:
    ! 1. The fuel stream is composed by hydrogen and nitrogen. If the inlet quantities are
    ! changed we are not in the fuel stream but in a mixture region.
    ! 2. The oxidizer stream is always composed by a certain amount of nitrogen and zero hydrogen.
    ! If these two quantities varies then we are not in the oxidizer stream.


    dt_phys = dt * adi % time_ref
    Tref_i  = 1.0_dp / adi % T_ref
    DTmax   = 0.0_dp
    rho_cgs = adi % rho_ref * 1.0e-3_dp


    ! domain to apply combustion (not in every boundary)
    if      ( ndim == 3 ) then
       ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 ) ! added ng
       iy = max ( ng+1 , sy ) ; fy = min ( ey , nty-ng-1 ) ! added ng
       iz = max ( ng+1 , sz ) ; fz = min ( ey , ntz-ng-1 ) ! added ng
    else if ( ndim == 2 ) then
       ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 ) ! added ng
       iy = max ( ng+1 , sy ) ; fy = min ( ey , nty-ng-1 ) ! added ng
       iz = sz                ; fz = ez
    else
       ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 ) ! added ng
       iy = sy                ; fy = ey
       iz = sz                ; fz = ez
    end if


    if ( dvodesimp ) then ! DVODE will be applied at certains points to save computational costs


       ! specify the pure fuel
       Yf (1:nrv) = inp % Y1f (1:nrv)
!       call X_to_Y ( thd , Yf )

       ! specify the pure oxidizer
       Yo (1:nrv) = inp % Y0o (1:nrv)
!       call X_to_Y ( thd , Yo )


       call_dvode_subdomain = .false.


       do k = iz , fz
          do j = iy , fy
             do i = ix , fx


                call_dvode   = .false.
                call_dvode_w = .false.
                call_dvode_f = .false.
                call_dvode_o = .false.


                rho_i = 1.0_dp / v (i,j,k,1)
                do l = 1 , nrv
                   Ya (l) = v (i,j,k, niv+l ) * rho_i
                end do


                rho0     = v (i,j,k,1)
                rho_phys = rho0 * rho_cgs
                T0_phys  = T (i,j,k) * adi % T_ref
                T_phys   = T0_phys


                ! verify the production rate of each species
                DeltaY (:) = Ya (:)
                call psr_simp_chemkin ( T_phys , rho_phys , DeltaY , dt_phys ) ! here dt_phys
                ! here DeltaY is the mass fractions variation array: Yf-Yi
                do l = 1 , nrv
                   if ( abs ( DeltaY(l) ) > eps_w ) call_dvode_w = .true.
                end do


                ! verify the pure fuel
                ! do l = 1 , nrv
                !    if ( abs ( Ya(l)-Yf(l) ) > eps2 ) call_dvode_f = .true.
                ! end do
                if ( abs ( Ya(inp % index_H2) - Yf(inp % index_H2) ) > eps_f .and. &
                     abs ( Ya(inp % index_N2) - Yf(inp % index_N2) ) > eps_f ) call_dvode_f = .true.


                ! verify the pure oxidizer
                ! do l = 1 , nrv
                !    if ( abs ( Ya(l)-Yo(l) ) > eps2 ) call_dvode_o = .true.
                ! end do
                if ( abs ( Ya(inp % index_H2) - Yo(inp % index_H2) ) > eps_o .and. &
                     abs ( Ya(inp % index_N2) - Yo(inp % index_N2) ) > eps_o ) call_dvode_o = .true.


                if (call_dvode_w)         call_dvode = .true.
                if ( .not. call_dvode_f ) call_dvode = .false.
                if ( .not. call_dvode_o ) call_dvode = .false.


                if (call_dvode) then


                   call_dvode_subdomain = .true.


                   call psr_chemkin ( T_phys , rho_phys , Ya , dt_phys )


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


    else ! DVODE is applied in every point of the computational domain


       do k = iz , fz
          do j = iy , fy
             do i = ix , fx


                rho_i = 1.0_dp / v (i,j,k,1)
                do l = 1 , nrv
                   Ya (l) = v (i,j,k, niv+l ) * rho_i
                end do

                rho0     = v (i,j,k,1)
                rho_phys = rho0 * rho_cgs ! in g/cm^3
                T0_phys  = T (i,j,k) * adi % T_ref
                T_phys   = T0_phys


                call psr_chemkin ( T_phys , rho_phys , Ya , dt_phys )


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
                v (i,j,k,niv+nrv) = rho0 * ( 1.0_dp - sumY )


             end do
          end do
       end do


    end if


  end subroutine reaction_dvode


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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine reaction_MIL ( inp , thd , adi , dt , dx_i , dy_i , dz_i , mu_SGS , T , v , DTmax )


    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    real (dp) , intent (in)                                        :: dt    !< time step
    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i  !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i  !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i  !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: mu_SGS!< turbulent viscosity
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v     !< conserved variables array
    real (dp) , intent (inout)                                     :: DTmax !< maximum temperature increase


    integer (ip)                  :: ix , fx , &
                                     iy , fy , &
                                     iz , fz
    integer (ip)                  :: i , j , k , l , ok
    real (dp)                     :: rho , rho_i

    ! Temperatures
    real (dp)                     :: T_phys , Tref_i , T0_phys

    ! Mass fractions
    real (dp) , dimension (1:nrv+npv+nvv) :: Ya
    real (dp)                             :: sumY , yo2

    ! Mixture fraction, its variance
    real (dp)                     :: zm , zp2

    ! Chemical production rate
    real (dp) , dimension (nrv)   :: wreac

    ! Others
    real (dp)                     :: dm_SGS , Sc_sgs_i , taumix , rvd0

    real (dp) , allocatable , dimension (:,:,:)   :: delta_2


    if ( npv > 1 ) then
       write (*,*) 'the code is not prepared to support more than one passive scalar'
       write (*,*) 'please reconsider this parameter'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if

    rvd0     = adi % sqgmr / adi % Sc ! = Ma sqrt(gamma) / ( Re Sc )
    Sc_sgs_i = 1.0_dp / inp % Sc_sgs
    Tref_i   = 1.0_dp / adi % T_ref
    DTmax    = 0.0_dp


    ! domain to apply combustion (not in every boundary)
    if      ( ndim == 3 ) then
       ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 ) ! added ng
       iy = max ( ng+1 , sy ) ; fy = min ( ey , nty-ng-1 ) ! added ng
       iz = max ( ng+1 , sz ) ; fz = min ( ey , ntz-ng-1 ) ! added ng
    else if ( ndim == 2 ) then
       ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 ) ! added ng
       iy = max ( ng+1 , sy ) ; fy = min ( ey , nty-ng-1 ) ! added ng
       iz = sz                ; fz = ez
    else
       ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 ) ! added ng
       iy = sy                ; fy = ey
       iz = sz                ; fz = ez
    end if


    ! Square filter-width
    allocate  ( delta_2 ( sx:ex , sy:ey , sz:ez )                   , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate reaction_MIL')

    call filter_width ( 0 , dx_i , dy_i , dz_i , delta_2 )
    delta_2 = delta_2 * delta_2


    do k = iz , fz
       do j = iy , fy
          do i = ix , fx

             rho   = v (i,j,k,1)
             rho_i = 1.0_dp / rho
             do l = 1 , nrv+npv+nvv
                Ya (l) = v (i,j,k, niv+l ) * rho_i
             end do

             yo2      = Ya (inp % index_O2)              ! mean oxygen mass fraction
             zm       = Ya (nrv+npv)                     ! mean mixture fractionn value
             zp2      = Ya (nrv+npv+1)                   ! mixture fraction variance value

if ( yo2 /= yo2 .or. zm /= zm .or. zp2 /= zp2 ) then
   write (*,*) 'subroutine reaction_MIL'
   write (*,'(A,3(1X,I5))') '(i,j,k) = ' , i,j,k
   write (*,'(A,15(1X,1PE10.3))') 'rho , Ya(1:nrv+npv) = ' , rho , Ya(1:nrv+npv)
   write (*,'(A,15(1X,1PE10.3))') 'yo2 , zm , zp2 = ' , yo2 , zm , zp2
end if

             T0_phys  = T (i,j,k) * adi % T_ref
             T_phys   = T0_phys

             ! Relaxation time ( + eps to avoid divisions by zero )
             dm_SGS   = adi % Sc * mu_SGS (i,j,k) * rho_i * Sc_SGS_i
             taumix   = delta_2 (i,j,k) / ( rvd0 * dm_SGS + eps )

             taumix   = taumix * adi % time_ref ! dimensionalize in second

             call wreac_MIL ( inp , thd , taumix , yo2 , zm , zp2 , wreac )

             wreac (:) = wreac (:) * adi % time_ref ! adimensionalize (s^-1 * s)


             T (i,j,k) = T_phys * Tref_i

             ! DTmax has physical dimensions: degrees Kelvin
             DTmax = max ( DTmax , abs ( T_phys - T0_phys ) ) ! T_phys is updated


             sumY = 0.0_dp
             do l = 1 , nrv-1
                Ya (l)          = Ya (l) + dt * wreac (l)
                Ya (l)          = min ( max ( Ya (l) , min_Y ) , max_Y )
                v (i,j,k,niv+l) = rho * Ya (l)
                sumY            = sumY + Ya (l)
             end do
             Ya (nrv)          = 1.0_dp - sumY
             Ya (nrv)          = min ( max ( Ya(nrv) , min_Y ) , max_Y )
             v (i,j,k,niv+nrv) = rho * ( 1.0_dp - sumY )

          end do
       end do
    end do


    deallocate  ( delta_2 )


  end subroutine reaction_MIL


!> \brief Calculate the chemical production rate of each species equations by using MIL
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine wreac_MIL ( inp , thd , taumix , yo2 , zm , zp2 , wreac )


    type (inp_type) , intent (in)                     :: inp   !< input derived type
    type (thd_type) , intent (in)                     :: thd   !< thermodynamic derived type
    real (dp) , intent (in)                           :: taumix!< mixing time scale (in s)
    real (dp) , intent (in)                           :: yo2   !< oxydizer mass fraction
    real (dp) , intent (in)                           :: zm    !< mixture fraction
    real (dp) , intent (in)                           :: zp2   !< SGS variance
    real (dp) , dimension (nrv) , intent (inout)      :: wreac !< chemical production rate (in s^-1)

    real (dp) , parameter         :: Da_c_diff = 1.0_dp ! critical diffusion Damkohler number

    integer (ip)                  :: cas

    ! Mass fractions
    real (dp)                     :: yoequ , yomix , yoiem , yoeqm , yoburnt , yomil

    ! Mixture fraction, its variance, jump positions
    real (dp)                     :: zjm , zjp , zp2max , zignim , zignip !, dzj

    ! Chemical production rate
    real (dp)                     :: wo2mil , wo2bur , wo2sauts , wo2 , prob_mil

    ! Reversed mixing time scale
    real (dp)                     :: taumix_i

    logical                       :: call_mil , call_self , call_reac

    ! Two-phase variables (change of variables)
    real (dp)                     :: max_Y_zi , yfmax_zi , ystmax_zi , & ! mass fraction
                                     zst_zi , dz_max_min_i               ! mixture fraction

    ! Beta-PDF
    real (dp)                     :: xb , beta1 , beta2

!Not used ?
    real (dp)                     :: c


    ! MIL marker and initialization
    call_mil  = .false.
    call_self = .false.
    call_reac = .true.
    wo2       = 0.0_dp
    wreac     = 0.0_dp

    taumix_i = 1.0_dp / taumix

    zp2max   = zm * ( max_Z - zm )

    ! Jump positions
    zjm      = zm
    zjp      = zm
!    dzj      = zjp - zjm ! NOT USED !?!?

    ! Chemical equilibrium line
    yoequ = 0.0_dp
    if ( zm <= zst ) yoequ = ( 1.0_dp - zm / zst ) * max_Y
    yoequ = min ( max ( yoequ , min_Y ) , max_Y )

    ! Mixing line
    yomix = max_Y * ( 1.0_dp - zm )

    ! IEM (Interaction by Exchange with the Mean) line
    yoiem = yo2


!    !====================================
!    ! check if AUTOIGNITION
!    !====================================

!    call auto_ign ( zm , yo2 , yh2 , taumix , Ps , Ec , Htot , ac_crit_MIL , &
!                    zignip , zignim , ac_crit_auto , call_self , taumix_i , tignimin , &
!                    tau_sej_ZM , tau_sej_AIR )                                                       ! TO ADD...


!    if ( .not. call_self ) then ! Exit only if NO AUTOIGNITION, (increase algorithm efficiency)

!       if ( yo2 < yoequ         &
!            .or. zm <= 0.001_dp &
!            .or. zm >= 0.999_dp ) return

       if ( ( yomix > eps8 ) .and. ( (yomix-yoequ) > eps8 ) ) then
          c = abs ( ( yomix - yo2 ) / ( yomix - yoequ + eps ) )                                       ! NOT USED !?
       end if

!       yomix = max_Y * ( 1.0_dp - zm )
!       ! ac inverse d'un Da
!       if ( inp % ext % rctmax < eps8 ) return

!       if ( taumix < ac_crit_MIL / inp % ext % rctmax ) then
!          if ( itest == 1 ) write (*,*)'Arret pas de sauts'
!          return
!       else
          call_mil = .true.
!       end if

!    else ! AUTOIGNITION possible test if diffusion flame possible

!       call_mil = .true.

!       if ( taumix < ac_crit_MIL / inp % ext % rctmax &
!            .or. yo2 <= yoequ                         &
!            .or. zm  <= 0.0025_dp                     &
!            .or. zm  >= 0.99_dp                       &
!            .or. yo2 >= 0.9975_dp * max_Y             &
!            .or. yo2 <  0.001_dp * max_Y              ) call_mil = .false.

!       if ( yomix > eps8 .and. (yomix-yoequ) > eps8 ) ) then
!          if ( zm >= zst ) then
!             c = abs ( ( yomix - yo2 ) / ( yomix + eps ) )                                           ! NOT USED !?
!          else
!             c = abs ( ( yomix - yo2 ) / ( yomix - yoequ + eps ) )
!          end if
!       end if

!       if ( c < 0.005_dp ) call_mil = .false.

!    end if ! .not. call_self


!    !====================================
!    ! Jump positions assessment
!    !====================================

!    if ( call_self .and. call_mil ) then ! MIL possible: AUTOIGNITION jumps comparison

!       call jump_mp ( inp , taumix , Da_c_diff , zjm , zjp )

!       zjm = min ( zignim , zjm , zm )
!       zjp = max ( zignip , zjm , zm )

!       cas = 1

!    else if ( .not. call_self .and. call_mil ) then

       call jump_mp ( inp , taumix , Da_c_diff , zjm , zjp )

       zjm = min ( zjm , zm )
       zjp = max ( zjp , zm )

!    else if ( call_self .and. .not. call_mil ) then

!       wo2 = wo_auto ( zignim , zignip , zm , max_Z , zp2 , yo2 , max_Y , zst , zone )               ! TO ADD...
!       wo2 = wo2 * taumix_i

!       cas = 3
!       call_reac = .false.

!    else ! ( .not. call_self .and. .not. call_mil )

!       !wo2 = 0.0_dp ! already initialize to zero

!       cas = 0
!       call_reac = .false.

!    end if ! ( call_self .and. call_mil )


!    !====================================
!    if ( twophase ) then
!    !====================================

!       ! Two-phase flow injection take into account the transient heating of the droplets
!       ! PAS PRISE EN COMPTE EVAPORATION COMBUSTIBLE (a modifier dans yo_mil et womil)
!       ! seulement cas evaporation lox traite

!       ! change of variables
!       max_Y_zi  = 1.0_dp - ( 1.0_dp - yo2 ) * min_Z / ( zm + eps ) ! = 1 for gas flow
!       yfmax_zi  = yo2 / ( 1.0_dp - zm + eps ) * ( 1.0_dp - max_Z ) ! = 0 for gas flow
!       ystmax_zi = 1.0_dp * ( 1.0_dp - min_Z / ( zm + eps ) )

!       dz_max_min_i = 1.0_dp / ( max_Z - min_Z )

!       zm  = ( zm - min_Z ) * dz_max_min_i
!       zm  = max ( 0.0_dp , zm )                        ! rescale zm in [0,1]

!       zp2 = zp2 * dz_max_min_i * dz_max_min_i
!       zp2 = min ( zp2 , zm * ( 1.0_dp - zm ) )         ! rescale zp2 in [0,zm(1-zm)]

!       zp2max = zp2max * dz_max_min_i * dz_max_min_i
!       zp2max = min ( zp2max , zm * ( 1.0_dp - zm ) )   ! rescale zp2max in [0,zm(1-zm)]

!       zjm = ( zjm - min_Z ) * dz_max_min_i
!       zjm = max ( 0.0_dp , zjm )                       ! rescale zjm in [0,1]

!       zjp = ( zjp - min_Z ) * dz_max_min_i
!       zjp = min ( 1.0_dp , zjp )                       ! rescale zjp in [0,1]

!       zst_zi = ( zst - min_Z ) * dz_max_min_i
!       zst_zi = max ( 0.0_dp , zst_zi  )

!    else

       max_Y_zi  = 1.0_dp ! = 1 for gas flow
       yfmax_zi  = 0.0_dp ! = 0 for gas flow
       ystmax_zi = 1.0_dp
       zst_zi    = zst

!    end if


    if ( zp2 <= eps4 * zp2max ) then ! around the average value

!       if ( zm < zst_zi ) then
!          wo2 = - ( yo2 * zst_zi - max_Y_zi * ( zst_zi - zm ) ) / ( taumix * zst_zi + eps )
!       else
!          wo2 = - yo2 * taumix_i
!       end if

       wo2 = 0.0_dp
       call_reac = .false.

    else if ( zp2 >= ( 1.0_dp - eps4 ) * zp2max ) then ! around the 2 extrem values

       wo2 = 0.0_dp
       call_reac = .false.

    end if


    if ( call_reac ) then

       !====================================
       ! Beta-PDF shape parameters
       !====================================

       xb = zp2max / zp2 - 1.0_dp
       ! not need to add eps in division because of the previous "if" condition
       if ( xb <= 0.0_dp ) xb = eps8
       beta1 = zm * xb
       beta2 = ( 1.0_dp - zm ) * xb
if ( xb /= xb .or. beta1 /= beta1 .or. beta2 /= beta2 ) then
   write (*,*) 'subroutine wreac_MIL'
   write (*,'(A,5(1X,1PE10.3))') 'zp2 , zp2max , 1/Seg = ' , zp2 , zp2max , zp2max / zp2
   write (*,'(A,5(1X,1PE10.3))') 'beta1 , beta2 = ' , beta1 , beta2
end if

       !==============================================================
       ! Burnt gases composition integration from min_Z to max_Z
       !==============================================================

       call yo_mil ( zm , zp2 , 0.0_dp , 1.0_dp , zst_zi , &
                     yo2 , yo2 , max_Y_zi , yfmax_zi , ystmax_zi , yoeqm , yoburnt )

       if ( yoeqm   < yoequ ) yoeqm   = yoequ
       if ( yoburnt < yoequ ) yoburnt = yoequ


       !===========================================================
       ! PDF between zjm-zjp for possible premixed reaction rate
       !===========================================================

!       intsmsp = intmg ( zjm , zjp , beta1 , beta2 )                                           ! intsmsp NOT USED and DEFINED !!!!
!       yoeqm   = 0.0_dp


       !===============================================
       ! MIL composition integration from zjm to zjp
       !===============================================

       call yo_mil ( zm , zp2 , zjm , zjp , zst_zi , &
                     yo2 , yo2 , max_Y_zi , yfmax_zi , ystmax_zi , yoeqm , yomil )


       if ( yomil > yo2 ) then

          prob_mil = ( yo2 - yoburnt ) / ( yomil - yoburnt + eps8 )
          prob_mil = min ( max ( yoequ , 0.0_dp ) , 1.0_dp )

          wo2bur = wo_mil ( zm , beta1 , beta2 , 0.0_dp , 1.0_dp , yo2 , zst_zi , ystmax_zi )
          wo2bur = wo2bur * taumix_i

          cas = 21

       else

          prob_mil = ( yo2 - yomix ) / ( yomil - yomix + eps8 )
          prob_mil = min ( max ( yoequ , 0.0_dp ) , 1.0_dp )

          wo2bur = 0.0_dp

          cas=22

       end if

       wo2mil = wo_mil ( zm , beta1 , beta2 , zjm , zjp , yo2, zst_zi , ystmax_zi )
       wo2mil = wo2mil * taumix_i

       wo2sauts = wo_saut ( zjm , zm , zp2 , yo2 , max_Y_zi , zst_zi , yfmax_zi , ystmax_zi )
       wo2sauts = wo2sauts + wo_saut ( zjp , zm , zp2 , yo2 , max_Y_zi , zst_zi , yfmax_zi , ystmax_zi )
       wo2sauts = wo2sauts * taumix_i

       wo2 = prob_mil * ( wo2mil + wo2sauts ) + ( 1.0_dp - prob_mil ) * wo2bur

    end if ! ( call_reac )


    ! physical limitation of production rate

!    wo2=max(wo2 , -(yo2-yoeqm)/dt_phys)
!    print* , "yo2" , yo2
!    print* , "yoeqm" , yoeqm

    ! back to the physical space, no change of variables to take into account the lower bound of the mixture fraction

    wo2 = min ( wo2 , 0.0_dp )

    ! Save jump positions, if autoignition zignim = zignim (because return)
    zignim = zjm
    zignip = zjp

    wreac ( inp % index_O2 ) = wo2

    wreac ( inp % index_H2 ) = wo2 * thd % Wc (inp % index_H2) * thd % Wc_i (inp % index_O2)
    wreac ( inp % index_H2 ) = wreac ( inp % index_H2 ) + wreac ( inp % index_H2 )

    wreac ( inp % index_H2O ) = - wo2 * thd % Wc (inp % index_H2O) * thd % Wc_i (inp % index_O2)
    wreac ( inp % index_H2O ) = wreac ( inp % index_H2O ) + wreac ( inp % index_H2O )


  end subroutine wreac_MIL


!> \brief Calculate the mixture fraction left and right jump positions
!!
!! Assess by extrapolation the left and right jump positions from the mixing time scale 'taumix' 
!! and the table of the reversed chemical time scale.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine jump_mp ( inp , taumix , ac , zjm , zjp )


    type (inp_type) , intent (in)              :: inp    !< input derived type
    real (dp) , intent (in)                    :: taumix !< mixing time scale
    real (dp) , intent (in)                    :: ac     !< control parameter
    real (dp) , intent (inout)                 :: zjm    !< mixture fraction left  jump position
    real (dp) , intent (inout)                 :: zjp    !< mixture fraction right jump position

    integer (ip)                               :: zj_index , iminus , iplus , i
    real (dp)                                  :: taumix_i , xorigin , yorigin , pente_i


    taumix_i = ac / taumix
    if ( taumix_i > inp % ext % rctmax ) return                                                         ! ce test doit etre fait BIEN AVANT de chercher a calculer les positions de sauts ! 
                                                                                                        ! Car si la condition n'est pas realise, il n'y a tout simplement pas de reaction !


    !=================================================
    ! Mixture fraction left jump position assessment
    !=================================================

    iminus = 1                         ! left (minus) index
    iplus  = inp % ext % index_rctmax  ! right (plus) index

    do

       ! sweep from left to right
       zj_index = iminus + (iplus-iminus) * 0.5

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
       zj_index = iminus + (iplus-iminus) * 0.5

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


!> \brief Calculate the scalar mean values in composition space.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine yo_mil ( zm , zp2 , zjm , zjp , zstoic , yo2m , yo2sub , yomax , yfmax , ystmax , &
                      yeqm , yo2mil )


    real (dp) , intent (in)                    :: zm     !< mixture fraction
    real (dp) , intent (in)                    :: zp2    !< SGS variance
    real (dp) , intent (in)                    :: zjm    !< mixture fraction left  jump position
    real (dp) , intent (in)                    :: zjp    !< mixture fraction right jump position
    real (dp) , intent (in)                    :: zstoic !< stoichiometric mixture fraction (modified if twophase)
    real (dp) , intent (in)                    :: yo2m   !< oxidizer mass fraction
    real (dp) , intent (in)                    :: yo2sub !< ???
    real (dp) , intent (in)                    :: yomax  !< maximum oxidizer mass fraction (modified if twophase)
    real (dp) , intent (in)                    :: yfmax  !< ???
    real (dp) , intent (in)                    :: ystmax !< ???
    real (dp) , intent (inout)                 :: yeqm   !< mean value on equilibrium trajectory
    real (dp) , intent (inout)                 :: yo2mil !< mean value on MIL trajectory

    real (dp) :: aeq , beqm , beqp , yz1 , yz2
    real (dp) :: int1 , int2 , int3 , ycalc , ymeq
    real (dp) :: xb , beta1 , beta2


    ymeq = 0.0_dp

    xb = ( zm * ( 1.0_dp - zm ) ) / zp2 - 1.0_dp
    if ( xb <= 0.0_dp ) xb = eps8
    beta1 = zm * xb
    beta2 = ( 1.0_dp - zm ) * xb

    if ( ( max_Z - zm ) <= eps8 ) then

       ycalc  = yo2m
       yo2mil = ycalc
       int1   = 0.0_dp
       int2   = 0.0_dp
       int3   = 0.0_dp

    else

       aeq  = ystmax - ystmax * zm / zstoic
       beqm = ystmax / zstoic
       beqp = - ystmax / zstoic
       yz1  = (yomax-yo2sub) / zm
       yz2  = - yo2sub / (1.0_dp-zm) + yfmax / (1.0_dp-zm)

       if ( zstoic < zjm ) then

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 )
          int1 = int1 + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = 0.0_dp

       else if ( zstoic < zm ) then

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 )
          int1 = int1 + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zstoic , beta1 , beta2 )
          ymeq = ymeq + beqm * intzmg ( zjm , zstoic , zm , beta1 , beta2 )

       else if ( zstoic < zjp ) then

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 )
          int1 = int1 + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zm , beta1 , beta2 )
          ymeq = ymeq + beqm * intzmg ( zjm , zm , zm , beta1 , beta2 )
          ymeq = ymeq + aeq * intpg ( zm , zstoic , beta1 , beta2 )
          ymeq = ymeq + beqp * intzpg ( zm , zstoic , zm , beta1 , beta2 )

       else if ( zstoic < max_Z ) then

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 )
          int1 = int1 + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zm , beta1 , beta2 )
          ymeq = ymeq + beqm * intzmg ( zjm , zm , zm , beta1 , beta2 )
          ymeq = ymeq + aeq * intpg ( zm , zjp , beta1 , beta2 )
          ymeq = ymeq + beqp * intzpg ( zm , zjp , zm , beta1 , beta2 )

       else

          int1 = intmg ( 0.0_dp , zjm , beta1 , beta2 )
          int1 = int1 + intpg ( zjp , max_Z , beta1 , beta2 )
          int2 = intzmg ( 0.0_dp , zjm , zm , beta1 , beta2 )
          int3 = intzpg ( zjp , max_Z , zm , beta1 , beta2 )

          ymeq = aeq * intmg ( zjm , zm , beta1 , beta2 )
          ymeq = ymeq + beqm * intzmg ( zjm , zm , zm , beta1 , beta2 )
          ymeq = ymeq + aeq * intpg ( zm , zjp , beta1 , beta2 )
          ymeq = ymeq + beqp * intzpg ( zm , zjp , zm , beta1 , beta2 )

       end if

       ycalc  = yo2sub * int1 + yz1 * int2 + yz2 * int3 + ymeq
       yo2mil = ycalc

    end if

    yeqm = ymeq


  end subroutine yo_mil


!> \brief Calculate ....
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
          wo_mil = wp * intmg ( zjm , zm , beta1 , beta2 )
          wo_mil = wo_mil + wp * intpg ( zm , zjp , beta1 , beta2 )
       else if ( zst < zm ) then
          wo_mil = wm * intmg ( zjm , zst , beta1 , beta2 )
          wo_mil = wo_mil + wp * intmg ( zst , zm , beta1 , beta2 )
          wo_mil = wo_mil + wp * intpg ( zm , zjp , beta1 , beta2 )
       else if ( zst < zjp ) then
          wo_mil = wm * intmg ( zjm , zm , beta1 , beta2 )
          wo_mil = wo_mil + wm * intpg ( zm , zst , beta1 , beta2 )
          wo_mil = wo_mil + wp * intpg ( zst , zjp , beta1 , beta2 )
       else if ( zst < max_Z ) then  
          wo_mil = wm * intmg ( zjm , zm , beta1 , beta2 )
          wo_mil = wo_mil + wm * intpg ( zm , zjp , beta1 , beta2 )
       else
          wo_mil = wm * intmg ( zjm , zm , beta1 , beta2 )
          wo_mil = wo_mil + wm * intpg ( zm , zjp , beta1 , beta2 )
       end if

       wo_mil = min ( wo_mil , 0.0_dp )

    end if


  end function wo_mil


!> \brief Calculate ....
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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

       aeq  = yomax - yomax * zm / zst
       beqm = yomax / zst
       yz1  = ( yomax - yo2sub ) / zm
       yz2  = - yo2sub / ( 1.0_dp - zm )

       if ( saut < zm ) then
         yo2iem = ( yo2sub - yomax ) / zm * saut + yomax
       else
         yo2iem = yfmax / ( 1.0_dp - zm ) * ( saut - zm ) + yo2sub / ( 1.0_dp - zm ) * ( 1.0_dp - saut )
       end if

       if ( saut < zst ) then
         yo2equ = ystmax * ( 1.0_dp - saut / zst )
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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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


end module solver_reaction
