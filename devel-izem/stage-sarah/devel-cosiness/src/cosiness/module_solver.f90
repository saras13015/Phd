!------------------------------------------------------------------------------
! MODULE: solver
!------------------------------------------------------------------------------
!> \brief Solver for the Navier-Stokes equations.
!!
!! This module contains all the subroutines and utilities needed to
!! solve the system of equations.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module solver

  use parameters
  use parallel
  use input
  use adim
  use type_thd
  use BCs
  use weno
  use viscflux
  use eg_lib
  use solver_reaction

  implicit none

contains


!> \brief Read an existing restart file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine read_restart ( inp , adi , thd , time , grid , T , W_i , cp , ha , v )


    type (inp_type) , intent (in)                                  :: inp  !< input derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , intent (in)                                        :: time !< time
    type (inp_grid) , intent (in)                                  :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip)                                                   :: ok , i , j , k , l
    integer (kind=8)                                               :: reclmax , rec_index
    character (len_default)                                        :: number_ascii , rank_ascii , adr_file
    real    (dp)                                                   :: dtmin
    real    (dp) , dimension (4)                                   :: dt_4
    real (dp) , allocatable , dimension (:,:,:)                    :: mu_SGS !< turbulent viscosity (for LES)

    dt_4 = 1.0_dp

    allocate ( v       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
               ha      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
               T       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               W_i     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               cp      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               mu_SGS  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
               stat = ok )


    if ( ok > 0 ) call abort_mpi ('error allocate read_restart')

    v   = 0.0_dp
    ha  = 0.0_dp
    T   = 0.0_dp
    W_i = 0.0_dp
    cp  = 0.0_dp
    mu_SGS = 0.0_dp


    write ( number_ascii , format_restart ) inp % number_restart
    adr_file = trim (dir_parent) // trim (inp % name_restart) // '_' // trim (number_ascii)
    if ( rank == rank_default ) write (*,*) 'reading previous restart file: ' , trim(adr_file)


    if (ind_files) then ! reading individual binary files per process

       reclmax = int ( (ex-sx+1) * (ey-sy+1) * (ez-sz+1) * nv * nbit , kind = 8 )

       write ( rank_ascii , format_restart ) rank
       adr_file = trim (adr_file) // '_' // trim (rank_ascii)

       rec_index = 1

    else ! reading one single binary file

       reclmax = int ( nxmax * nymax * nzmax * nv * nbit , kind = 8 )

       rec_index = rank + 1

    end if


    open ( unit_restart , file   = trim (adr_file) , &
                          access = 'direct'        , &
                          recl   = reclmax         , &
                          form   = 'unformatted'   , &
                          status = 'old'           , &
                          iostat =  ok             )

    if ( ok /= 0 .and. rank == rank_default ) call abort_mpi ('error opening ' // trim (adr_file))
    call mpi_barrier ( MPI_COMM_WORLD , mpicode )

    read ( unit_restart , rec = rec_index ) & ! it must start at 1
         (((( v (i,j,k,l) , i = sx , ex ) , &
                            j = sy , ey ) , &
                            k = sz , ez ) , &
                            l = 1  , nv )

    close (unit_restart)


!    if (LES) call mu_SGS_selector_post ( inp , adi , dx_i , dy_i , dz_i , v , mu_SGS )


    call comm_cons (v)
    call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
    call upd_extrapolation (v)
    call timestep ( inp , adi , thd , grid , cp , W_i , T , v , mu_SGS , 1.0_dp , dt_4 , dtmin )
    call upd_boundaries ( inp , adi , thd , time , dtmin , grid , T , W_i , cp , ha , v )
    do l = 1 , ndim+ndim
       call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
    end do

    deallocate ( mu_SGS )


  end subroutine read_restart


!> \brief Calculate the time step.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine timestep ( inp , adi , thd , grid , cp , W_i , T , v , mu_SGS , DTmax , dt_4 , dtmin )


    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (adi_type)                                                :: adi   !< non-dimensional derived type
    type (thd_type)                                                :: thd   !< thermodynamic derived type
    type (inp_grid) , intent (in)                                  :: grid  !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp    !< heat capacity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i   !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v     !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: mu_SGS!< turbulent viscosity (for LES)
    real (dp) , intent (in)                                        :: DTmax !< maximum temperature increase
    real (dp) , dimension (:) , intent (inout)                     :: dt_4  !< array of different time steps
    real (dp) , intent (inout)                                     :: dtmin !< minimum time step


    integer (ip)                                          :: ok , i , j , k , l
    real (dp)                                             :: rho_i   , ux       , vy     , wz       , &
                                                             dtcx_i  , dtcy_i   , dtcz_i , &
                                                             dtc     , diffmax  , gamma  , P        , &
                                                             csound
    real (dp) , dimension (4)                             :: dt_tmp
    real (dp) , allocatable , dimension (:,:,:,:)         :: Xa , dm
    real (dp) , allocatable , dimension (:,:,:)           :: ct , delta_2

    integer (ip)                                          :: ldmax
    real (dp)                                             :: rvd0_i , rlamcp0_i , Pr_sgs_i


    ! fixed time step
    if (inp % fixeddt) then
       dtmin = inp % dt / adi % time_ref
       dt_4 = dtmin
       return
    end if

    rvd0_i     = adi % Sc / adi % sqgmr
    rlamcp0_i  = adi % Pr / adi % sqgmr
    Pr_sgs_i   = 1.0_dp / inp % Pr_sgs

    dtc        = 1.0_dp
    dt_tmp (:) = 1.0_dp
    dt_4 (:)   = 1.0_dp


    ! Convective criteria
    do k = sz , ez
       do j = sy , ey
          do i = sx , ex

             rho_i  = 1.0_dp / v (i,j,k,1)
             ux     = v (i,j,k,2) * rho_i
             vy     = v (i,j,k,3) * rho_i
             wz     = v (i,j,k,4) * rho_i

             gamma  = thd % gam2 * W_i (i,j,k) / cp (i,j,k)
             gamma  = 1.0_dp / ( 1.0_dp - gamma )
             P      = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             csound = sqrt ( gamma * P * rho_i )

             dtcx_i = ( abs (ux) + csound ) * grid % dx_i (i)
             dtcy_i = ( abs (vy) + csound ) * grid % dy_i (j)
             dtcz_i = ( abs (wz) + csound ) * grid % dz_i (k)

             dtc = max ( dtc , dtcx_i , dtcy_i , dtcz_i )

          end do
       end do
    end do

!    if ( .not. inp % read_restart ) then
!       ite_cst_i = 500.0_dp ! ite from which the CFL is constant to inp % CFL
!       ite_cst_i = 5.0_dp / ite_cst_i
!       CFL0 = CFLmin
!
!       CFL = ( inp % CFL - CFL0 ) * ( 1.0_dp - exp ( - ite * ite_cst_i ) ) + CFL0
!    end if

    dt_4 (1) = CFL / dtc


    ! To estimate the time step, the viscous criteria will be
    ! the one from simple transport even if EGLIB is being used
    if ( vis .and. vis_dt ) then

       allocate ( ct       ( sx:ex , sy:ey , sz:ez )         , &
                  Xa       ( sx:ex , sy:ey , sz:ez , 1:nrv ) , &
                  dm       ( sx:ex , sy:ey , sz:ez , 1:nrv ) , &
                  delta_2  ( sx:ex , sy:ey , sz:ez )         , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate timestep')


       call molefrac            ( 0 , thd , W_i , v , Xa )
       call prim_vis_vars_wo_mu ( 0 , thd , v , W_i , T , Xa , dm , ct )
       delta_2 (:,:,:) = grid % delta (sx:ex,sy:ey,sz:ez)
       delta_2 = delta_2 * delta_2


       ! Mass diffusivity
       ldmax = 1
       diffmax = dm (sx,sy,sz,1)
       do l = 1 , nrv
          do k = sz , ez
             do j = sy , ey
                do i = sx , ex
                   if ( diffmax < dm (i,j,k,l) ) then
                      ldmax = l
                      diffmax = dm (i,j,k,l)
                   end if
                end do
             end do
          end do
       end do

       dt_4 (2) = rvd0_i * Fo                                           &
                * minval ( delta_2 (sx:ex,sy:ey,sz:ez)                  &
                         / ( dm (sx:ex,sy:ey,sz:ez,ldmax)               &
                           + adi % Sc * mu_SGS (sx:ex,sy:ey,sz:ez)      &
                           / ( v (sx:ex,sy:ey,sz:ez,1) * inp % Sc_sgs )))

       ! Thermal diffusivity

       dt_4 (3) = rlamcp0_i * Fo                                                 &
                * minval ( delta_2 (sx:ex,sy:ey,sz:ez) * v (sx:ex,sy:ey,sz:ez,1) &
                         / ( ct (sx:ex,sy:ey,sz:ez) / cp (sx:ex,sy:ey,sz:ez)     &
                           + adi % Pr * mu_SGS (sx:ex,sy:ey,sz:ez) * Pr_sgs_i )  )


       deallocate ( ct , Xa , dm , delta_2 )

    end if


    ! chemical criteria
    if (reaction) then

       ! if 2 convective and 1 reactive: dtmin and dTmax is somehow overstimated
       ! if 1 convective and 2 reactive: 0.5 * dtmin and dTmax is exact

       dt_4 (4) = dtmin * CPR / DTmax       ! maximum CPR degrees Kelvin of temperature
       dt_4 (4) = min ( dt_4 (4) , 1.0_dp )

       ! The minimum chemical time step is CFLmin
       dt_4 (4) = max ( dt_4 (4) , dt_4 (1) * CFLmin / CFL )

    end if

    ! MPI global time step evaluation
    call mpi_allreduce ( dt_4 , dt_tmp , 4 , MPI_DOUBLE_PRECISION , &
                         MPI_MIN , MPI_COMM_WORLD , mpicode )
!    call mpi_allreduce ( dt_4 , dt_tmploc , 4 , MPI_DOUBLE_PRECISION , &
!                         MPI_MINLOC , MPI_COMM_WORLD , mpicode )

!write (*,*) 'dtmin index' , dt_temploc


    dt_4  = dt_tmp

    dtmin = minval (dt_4)

    if ( inp % dtlimit .and. &
         dtmin * adi % time_ref < inp % dtmin_allowed ) dtmin = inp % dtmin_allowed / adi % time_ref


  end subroutine timestep


!> \brief Calculate the source/sink terms of variance equation (for LES)
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine source_sink ( inp , thd , adi , time , dt , grid , mu_SGS , T , W_i , cp , ha , v )


    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    real (dp) , intent (in)                                        :: time  !< time
    real (dp) , intent (in)                                        :: dt    !< time step
    type (inp_grid) , intent (in)                                  :: grid  !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: mu_SGS!< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i   !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp    !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha    !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v     !< conserved variables array


    integer (ip)                  :: ok , i , j , k , l
    integer (ip)                  :: domain_id , i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                     :: rho , rho_i , tau

    real (dp) , allocatable , dimension (:,:,:)   :: ps

    ! Constants
    real (dp) , parameter         :: twothirds = 2.0_dp / 3.0_dp , &
                                     C_xi      = 1.0_dp          , &
                                     eps       = 1.0e-15_dp
    ! Non-dimensionnal values
    real (dp)                     :: Pr_sgs_i , Sc_sgs_i , rvd0 ,  rlamcp0

    real (dp)                     :: grad2 , grad2max , wrkcst , wrktau , wrkexp 

    real (dp) , allocatable , dimension (:,:,:)   :: dx_ps , dy_ps , dz_ps , delta_2 , ct

    real (dp) , allocatable , dimension (:,:,:,:) :: Xa

    ! Dissipation and Production term
    real (dp)                     :: dissip
!    real (dp) , allocatable , dimension (:,:,:)   :: dissip , prod


    Pr_sgs_i = 1.0_dp / inp % Pr_sgs
    Sc_sgs_i = 1.0_dp / inp % Sc_sgs

    rvd0    = adi % sqgmr / adi % Sc ! = Ma sqrt(gamma) / ( Re Sc )
    rlamcp0 = adi % sqgmr / adi % Pr ! = Ma sqrt(gamma) / ( Re Pr )


    allocate  ( ps      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                dx_ps   ( sx:ex , sy:ey , sz:ez )                   , &
                dy_ps   ( sx:ex , sy:ey , sz:ez )                   , &
                dz_ps   ( sx:ex , sy:ey , sz:ez )                   , &
                delta_2 ( sx:ex , sy:ey , sz:ez )                   , &
!                prod    ( sx:ex , sy:ey , sz:ez )                   , &
!                dissip  ( sx:ex , sy:ey , sz:ez )                   , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate source_sink')

    allocate  ( Xa      ( sx:ex , sy:ey , sz:ez , nrv )             , &
                ct      ( sx:ex , sy:ey , sz:ez )                   , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate source_sink 2')

    call molefrac ( 0 , thd , W_i , v , Xa )
    call prim_vis_ct ( 0 , thd , T , Xa , ct )


    dx_ps = 0.0_dp ; dy_ps = 0.0_dp ; dz_ps = 0.0_dp

    if ( npv > 1 ) &
       call abort_mpi ('ERROR: subroutine source_sink npv > 1')


    ! Loop to define ps in the inner domain AND on the ghost points, 
    ! to calculate the first derivative.
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ps (i,j,k) = v (i,j,k,niv+nrv+npv) * rho_i ! Passive scalar (need to change here if npv >1)
             end do
          end do
       end do

    end do


    ! Passive scalar fisrt derivative
    if ( ndim >= 1 ) call dx ( grid % dx_i , ps , dx_ps )
    if ( ndim >= 2 ) call dy ( grid % dy_i , ps , dy_ps )
    if ( ndim == 3 ) call dz ( grid % dz_i , ps , dz_ps )


    ! Square filter-width
    delta_2 (:,:,:) = grid % delta (sx:ex,sy:ey,sz:ez)
    delta_2 = delta_2 * delta_2


!    do l = niv+nrv+npv+1 , niv+nrv+npv+nvv
       do k = sz , ez
          do j = sy , ey
             do i = sx , ex


                rho = v (i,j,k,1)


                ! Squared gradient
                ! if ( dim == 1 ) then ... else ... else ...
                grad2 =   dx_ps (i,j,k) * dx_ps (i,j,k) & 
                        + dy_ps (i,j,k) * dy_ps (i,j,k) &
                        + dz_ps (i,j,k) * dz_ps (i,j,k)

                ! if ( dim == 1 ) then ... else ... else ...
                grad2max = grid % dx_i (i) * grid % dx_i (i) &
                         + grid % dy_i (j) * grid % dy_i (j) &
                         + grid % dz_i (k) + grid % dz_i (k)

                if ( grad2 < 0.0_dp ) call abort_mpi ('error: source_sink, grad2 < 0')
                grad2 = min ( grad2 , grad2max )


                ! Relaxation time ( + eps to avoid divisions by zero )
                tau = delta_2 (i,j,k) * inp % Sc_sgs * rho / ( adi % Sc * mu_SGS (i,j,k) + eps )

                wrktau = - rvd0 * 2.0_dp * C_xi / tau ! independent to SGS variance => necessarily < 0
                if ( wrktau > 0.0_dp ) call abort_mpi ('error: source_sink, wrktau > 0')

                wrkexp = exp ( wrktau * dt )


               !=======================================================
               ! SGS Variance equation
               !=======================================================

                l = niv+6

!                prod   = 0.0_dp
!                dissip = 0.0_dp

!                ! Production term
!                prod = rvd0 * 2.0_dp * adi % Sc * Sc_sgs_i * mu_SGS (i,j,k) * grad2
!                if ( prod < 0.0_dp ) write (*,*) 'error: source_sink/SGSvariance, prod < 0'

!                ! Dissipation term
!                dissip = wrktau * v (i,j,k,l)

!                xiv1        = v (i,j,k,l) + dt * ( prod + dissip )

                   ! 2eme facon d'ecrire la meme chose...
                   !=======================================================

                wrkcst = rvd0 * 2.0_dp * adi % Sc * Sc_sgs_i * mu_SGS (i,j,k) * grad2

                v (i,j,k,l) = v (i,j,k,l) * wrkexp - wrkcst / wrktau * ( 1.0_dp - wrkexp )


                !=======================================================
                ! Squared filtered passive scalar equation \tilde{\xi}^2
                !=======================================================

                l = niv+8

                dissip = 0.0_dp

                ! Dissipation term
                dissip = - 2.0_dp * ( rlamcp0 * ct (i,j,k) / cp (i,j,k) &
                                    + rvd0 * adi % Sc * Sc_sgs_i * mu_SGS (i,j,k) ) * grad2
!                if ( dissip > 0.0_dp ) call abort_mpi ('error: source_sink/Squared filtered scalar, dissip > 0')
                if ( dissip > 0.0_dp ) write (*,*) 'error: source_sink/Squared filtered scalar, dissip > 0'

                v (i,j,k,l) = v (i,j,k,l) + dt * dissip


                !=======================================================
                ! Filtered squared passive scalar equation \tilde{\xi^2}
                !=======================================================

                l = niv+7

!                dissip = 0.0_dp

!                ! Dissipation term
!                dissip = - rlamcp0 * 2.0_dp * ct (i,j,k) / cp (i,j,k) * grad2
!!                if ( dissip > 0.0_dp ) call abort_mpi ('error: source_sink/Filtered squared scalar_1, dissip > 0')
!                if ( dissip > 0.0_dp ) write (*,*) 'error: source_sink/Filtered squared scalar_1, dissip > 0'

!                wrk = wrktau * ( v (i,j,k,niv+7) - rho * ( ps (i,j,k) * ps (i,j,k) ) )
!                dissip = dissip + wrk

!                xiv1        = v (i,j,k,l) + dt * dissip

                    ! 2eme methode d'integration...
                    !=======================================================

                wrkcst = - rlamcp0 * 2.0_dp * ct (i,j,k) / cp (i,j,k) * grad2 - wrktau * rho * ( ps (i,j,k) * ps (i,j,k) )

                v (i,j,k,l) = v (i,j,k,l) * wrkexp - wrkcst / wrktau * ( 1.0_dp - wrkexp )


                !=======================================================
                ! Departure from maximal variance equation \Delta_\xi = \xi*(1-\xi) - variance
                !=======================================================

                l = niv+9
!
!                prod = 0.0_dp
!
!                ! Production term
!                prod = rlamcp0 * 2.0_dp * ct (i,j,k) / cp (i,j,k) * grad2
!                if ( prod < 0.0_dp ) call abort_mpi ('error: source_sink/Dvar_1, prod < 0')
!
!                wrk = - wrktau * ( rho * ps (i,j,k) * ( 1.0_dp - ps (i,j,k) ) - v (i,j,k,niv+9) )
!                if ( wrk < 0.0_dp ) write (*,*) 'error: source_sink/Dvar_2, prod < 0'
!                prod = prod + wrk
!
!                xiv1        = v (i,j,k,l) + dt * prod
!
!                    ! 2eme methode d'integration...
!                    !=======================================================
!
                wrkcst = rlamcp0 * 2.0_dp * ct (i,j,k) / cp (i,j,k) * grad2 - &
                &        wrktau * ( v (i,j,k,niv+5) - rho * ( ps (i,j,k) * ps (i,j,k) ) )

                v (i,j,k,l) = v (i,j,k,l) * wrkexp - wrkcst / wrktau * ( 1.0_dp - wrkexp )


             end do
          end do
       end do
!    end do


    ! Updating the solution
!    do l = niv+nrv+npv+1 , niv+nrv+npv+nvv
!       v (sx:ex,sy:ey,sz:ez,l) = v (sx:ex,sy:ey,sz:ez,l) & 
!                                 + dt * ( prod (sx:ex,sy:ey,sz:ez) + dissip (sx:ex,sy:ey,sz:ez) )
!    end do


    deallocate ( Xa , ct )

    deallocate ( dx_ps , dy_ps , dz_ps )
    deallocate ( delta_2 , ps )
!    deallocate ( prod , dissip )


    call comm_cons (v)
    call upd_boundaries ( inp , adi , thd , time , dt , grid , T , W_i , cp , ha , v )
    call upd_prim_var_ghost ( thd , T , W_i , cp , ha , v )


  end subroutine source_sink


!> \brief 3rd-order Runge-Kutta loop.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine rk3 ( inp , thd , adi , time0 , dt , grid , mu_SGS , rweg , iweg , T , W_i , cp , ha , v )


    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    real (dp) , intent (in)                                        :: time0 !< time
    real (dp) , intent (in)                                        :: dt    !< time step
    type (inp_grid) , intent (in)                                  :: grid  !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mu_SGS!< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:) , intent (inout)       :: rweg  !< real EGLIB work array
    integer (ip) , allocatable , dimension (:) , intent (inout)    :: iweg  !< integer EGLIB work array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i   !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp    !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha    !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v     !< conserved variables array


    integer (ip) , parameter :: nk = 3

    real (dp) , dimension (nk) , parameter :: &

       bk3 = (/ 0.00_dp , 0.75_dp , 1.00_dp / 3.00_dp /) , &
       ck3 = (/ 1.00_dp , 0.25_dp , 2.00_dp / 3.00_dp /) , &
       tk3 = (/ 0.0_dp , 0.5_dp , 1.0_dp /)

    integer (ip) :: ok , ik , i , j , k , l

    real (dp) :: timek

    real (dp) , allocatable , dimension (:,:,:,:)                   :: v0 , f


    allocate ( v0 (sx:ex,sy:ey,sz:ez,1:nv) , &
               f  (sx:ex,sy:ey,sz:ez,1:nv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate rk3')


    do l = 1 , nv
       do k = sz , ez
          do j = sy , ey
             do i = sx , ex
                v0 (i,j,k,l) = v (i,j,k,l)
!                f  (i,j,k,l) = 0.0_dp      ! This is not necessary if euler_LLF_x starts
             end do
          end do
       end do
    end do


    do ik = 1 , nk


       timek = time0 + dt * tk3 (ik) ! time at irk-th RK stage


       ! Eulerian inviscid fluxes
       if ( ndim >= 1 ) call euler_LLF_x ( thd , adi , grid , T , W_i , cp , ha , v , f )
       if ( ndim >= 2 ) call euler_LLF_y ( thd , adi , grid , T , W_i , cp , ha , v , f )
       if ( ndim == 3 ) call euler_LLF_z ( thd , adi , grid , T , W_i , cp , ha , v , f )


       ! Viscous fluxes
       if (vis) then

          if ( .not. eglib ) then ! without EGLIB

             if (.not. LES) then ! without LES

                if ( ndim == 1 ) call viscflux_x   ( thd , adi , grid , T , W_i , cp , ha , v , f )
                if ( ndim == 2 ) call viscflux_xy  ( thd , adi , grid , T , W_i , cp , ha , v , f )
                if ( ndim == 3 ) call viscflux_xyz ( thd , adi , grid , T , W_i , cp , ha , v , f )

             else if (LES) then

                if ( ndim == 1 ) call viscflux_LES_x   ( inp , thd , adi , grid , T , W_i , cp , ha , v , mu_SGS , f )
                if ( ndim == 2 ) call viscflux_LES_xy  ( inp , thd , adi , grid , T , W_i , cp , ha , v , mu_SGS , f )
                if ( ndim == 3 ) call viscflux_LES_xyz ( inp , thd , adi , grid , T , W_i , cp , ha , v , mu_SGS , f )

             end if

          else if (eglib) then ! without EGLIB

             if (.not. LES) then ! without LES

                if ( ndim == 1 ) call viscflux_eg_x   ( thd , adi , grid , T , W_i , cp , ha , v , rweg , iweg , f )
                if ( ndim == 2 ) call viscflux_eg_xy  ( thd , adi , grid , T , W_i , cp , ha , v , rweg , iweg , f )
                if ( ndim == 3 ) call viscflux_eg_xyz ( thd , adi , grid , T , W_i , cp , ha , v , rweg , iweg , f )

             else if (LES) then

                if ( ndim == 1 ) call viscflux_eg_LES_x   ( inp , thd , adi , &
                  & grid , T , W_i , cp , ha , v , rweg , iweg , mu_SGS , f )
                if ( ndim == 2 ) call viscflux_eg_LES_xy  ( inp , thd , adi , &
                  & grid , T , W_i , cp , ha , v , rweg , iweg , mu_SGS , f )
                if ( ndim == 3 ) call viscflux_eg_LES_xyz ( inp , thd , adi , &
                  & grid , T , W_i , cp , ha , v , rweg , iweg , mu_SGS , f )

             end if

          end if

       end if


       ! Non reflecting boundary conditions
       call NR_boundaries ( thd , grid , T , W_i , cp , ha , v , f )


       ! Updating the solution
       do l = 1 , nv
          do k = sz , ez
             do j = sy , ey
                do i = sx , ex
                   v (i,j,k,l) = v0 (i,j,k,l) * bk3 (ik) + ck3 (ik) * &  ! v is a flux: do not update it
                               ( v  (i,j,k,l) - dt * f (i,j,k,l) )
                end do
             end do
          end do
       end do


       ! v is not a flux anymore but the array of conserved variables: update it
       call upd_prim_var_domain ( thd , T , W_i , cp , ha , v )
       if ( inp % perturb ) call perturbation ( adi , thd , grid % x , grid % y , grid % z , T , W_i , v )
       if ( inp % filter  ) call smooth_filter( adi , inp , grid % x , v )
       call comm_cons (v)
       call upd_boundaries ( inp , adi , thd , timek , dt , grid , T , W_i , cp , ha , v )
       call upd_prim_var_ghost ( thd , T , W_i , cp , ha , v )


    end do ! end of RK loop


    deallocate ( v0 , f )


  end subroutine rk3


!> \brief Add white-noise perturbation for the shear layer.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine perturbation ( adi , thd , x , y , z , T , W_i , v )


    type (adi_type) , intent (in)                                  :: adi !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: x   !< x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: y   !< x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: z   !< x-coordinate array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T   !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v   !< conserved variables array


    integer (ip)                                                   :: i , j , k , l , nperiod
    real (dp)                                                      :: u1_ , u2_ , cs1 , cs2 , uc , alpha , expo
    real (dp)                                                      :: epsy , epsphasey
    real (dp)                                                      :: profily , vpy
    real (dp)                                                      :: x0 , y0 , z0 , dw0 , dy02_i
    real (dp)                                                      :: P , rho_i , hm
    real (dp)                                                      :: has (nrv)


    if ( ndim == 1 ) then
       call abort_mpi ('no perturbation in 1D problems')
    else if ( ndim == 2 ) then
       call random_number (epsy)
       epsy = epsy + epsy - 1.0_dp
    else if ( ndim == 3 ) then
       call random_number (epsy)
       call random_number (epsphasey)
       epsy      = epsy + epsy - 1.0_dp
       epsphasey = epsphasey * pi
    end if


    ! nperiod is the number of periods in the z-direction
    nperiod = 3

    u1_ = 100.0_dp / adi % u_ref ; cs1 = 340.0_dp / adi % u_ref
    u2_ = 50.0_dp / adi % u_ref  ; cs2 = 340.0_dp / adi % u_ref
    uc  = ( cs1 * u1_ + cs2 * u2_ ) / ( cs1 + cs2 )

    dw0 = 1.6e-3_dp / adi % L_ref
    if (reaction) then
       alpha  = 1.0e-3_dp
       dy02_i = dw0
    else
       alpha  = 1.0e-2_dp
       dy02_i = 0.25_dp * dw0
    end if
    dy02_i = 1.0_dp / ( dy02_i * dy02_i )

!    Lz = 30.072_dp * dw0 ! (j'ai deja defini les differentes longueurs du domaine de calcul)

    x0 = 4.0_dp * dw0
    y0 = 0.0_dp
    z0 = 0.5_dp * Lz


    if ( ndim == 2 ) then


       do k = sz , ez
          do j = sy , ey
             do i = sx , ex

                expo = ( x(i) - x0 ) * ( x(i) - x0 ) + ( y(j) - y0 ) * ( y(j) - y0 )
                expo = exp ( - expo * dy02_i )
                vpy  = epsy * alpha * uc * expo

                rho_i = 1.0_dp / v (i,j,k,1)
                P     = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)

                call ha_scalar ( thd , T (i,j,k) , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * v (i,j,k,niv+l) * rho_i
                end do


                v (i,j,k,3) = v (i,j,k,3) + v (i,j,k,1) * vpy
                v (i,j,k,5) = v (i,j,k,1) * hm - P + 0.5_dp * rho_i * &
                            ( v (i,j,k,2) * v (i,j,k,2) +             &
                              v (i,j,k,3) * v (i,j,k,3) +             &
                              v (i,j,k,4) * v (i,j,k,4) )

             end do
          end do
       end do


    else if ( ndim == 3 ) then


       do k = sz , ez

          profily = cos ( (pi+pi) * ( z(k) - z0 ) * nperiod / Lz + epsphasey )

          do j = sy , ey
             do i = sx , ex

                expo = ( x(i) - x0 ) * ( x(i) - x0 ) + ( y(j) - y0 ) * ( y(j) - y0 )
                expo = exp ( - expo * dy02_i )

                vpy = epsy * alpha * uc * expo * profily

                rho_i = 1.0_dp / v (i,j,k,1)
                P     = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)

                call ha_scalar ( thd , T (i,j,k) , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * v (i,j,k,niv+l) * rho_i
                end do


                v (i,j,k,3) = v (i,j,k,3) + v (i,j,k,1) * vpy
                v (i,j,k,5) = v (i,j,k,1) * hm - P + 0.5_dp * rho_i * &
                            ( v (i,j,k,2) * v (i,j,k,2) +             &
                              v (i,j,k,3) * v (i,j,k,3) +             &
                              v (i,j,k,4) * v (i,j,k,4) )

             end do
          end do
       end do


    end if


  end subroutine perturbation


!> \brief Add smooth filter to the solution.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine smooth_filter ( adi , inp , x , v )


    type (inp_type)                                  , intent (in)    :: inp !< input derived type
    type (adi_type)                                  , intent (in)    :: adi !< non-dimensional derived type
    real    (dp) , allocatable , dimension (:)       , intent (in)    :: x   !< x-coordinate array
    real    (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v   !< conserved variables array


    integer (ip) :: csx , cex
    real (dp)    :: xini  , xend

    ! parameter for parabolic eps (can be introduce in input)
    xini  = inp % fil_xini / adi % L_ref
    xend  = inp % fil_xend / adi % L_ref

    call comm_cons (v)

    call loc_smooth_x ( x , xini , csx , cex  )
    if ( csx  /= 0 ) call filterx ( x , csx  , cex  , xini  , xend  , v )


  end subroutine smooth_filter


!> \brief x-direction filter for the sponge zone.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine filterx ( x , csx , cex , xini , xend , v )


    real    (dp) , allocatable , dimension (:)       , intent (in)    :: x
    integer (ip)                                     , intent (in)    :: csx
    integer (ip)                                     , intent (in)    :: cex
    real    (dp)                                     , intent (in)    :: xini
    real    (dp)                                     , intent (in)    :: xend
    real    (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v


    real    (dp) , allocatable , dimension (:,:,:,:)                  :: wrk
    integer (ip) :: i , j, k , l , ok
    real (dp)    :: a_W , b_E , c_S , d_N , e_B , f_F , mid , sum , epsx , vmoy

    allocate  ( wrk ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv ) , stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate filter')


    do l = 1 , nv
       do k = sz , ez
          do j = sy , ey
             do i = csx , cex

                a_W = 1.0_dp
                b_E = 1.0_dp
                c_S = 1.0_dp
                d_N = 1.0_dp
                e_B = 1.0_dp
                f_F = 1.0_dp
                mid = 1.0_dp

                if ( i==1 )   a_W = 0.0_dp
                if ( i==ntx ) b_E = 0.0_dp
                if ( j==1 )   c_S = 0.0_dp
                if ( j==nty ) d_N = 0.0_dp
                if ( k==1 )   e_B = 0.0_dp
                if ( k==ntz ) f_F = 0.0_dp

                sum = a_W + b_E + c_S + d_N + e_B + f_F + mid

                vmoy = a_W * v (i-1,j,k,l) + b_E * v (i+1,j,k,l) + &
                       c_S * v (i,j-1,k,l) + d_N * v (i,j+1,k,l) + &
                       e_B * v (i,j,k-1,l) + f_F * v (i,j,k+1,l) + &
                       mid * v (i,j,k,l)

                vmoy = vmoy / sum

                epsx = ( ( x(i) - xini ) / ( xend - xini ) ) * &  ! Puissance 1
                       ( ( x(i) - xini ) / ( xend - xini ) ) * &  ! Puissance 2
                       ( ( x(i) - xini ) / ( xend - xini ) ) * &  ! Puissance 3
                       ( ( x(i) - xini ) / ( xend - xini ) )      ! Puissance 4

                wrk (i,j,k,l) = ( 1.0_dp - epsx ) * v(i,j,k,l) + epsx * vmoy

             end do
          end do
       end do
    end do

    do l = 1 , nv
       do k = sz , ez
          do j = sy , ey
             do i = csx , cex
                v (i,j,k,l) = wrk (i,j,k,l)
             end do
          end do
       end do
    end do

    deallocate (wrk)


  end subroutine filterx


!> \brief locate x-coordinate of the filter for the sponge zone.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine loc_smooth_x ( x , s_x , csx , cex )


    real (dp) , dimension (:) , allocatable , intent (in)    :: x
    real (dp)                               , intent (in)    :: s_x
    integer (ip)                            , intent (inout) :: csx
    integer (ip)                            , intent (inout) :: cex


    integer (ip) :: i

    cex = 0
    csx = 0

    if ( s_x < x (sx) ) then
       csx = sx
       cex = ex
    else
       do i = ex , sx , -1
          if ( x (i) > s_x ) then
             csx = i
             cex = ex
          end if
       end do
    end if

    if ( neigh (E) == MPI_PROC_NULL ) cex = cex - 1


  end subroutine loc_smooth_x


!> \brief Filter for the outlet boundary condition. (Subroutine not used...?!)
!!
!! This filter is used whenever the solution accumulates errors at the
!! outlet BC. It provides a _zero gradient_ BC through a determined
!! number of grid points near that boundary.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
!  subroutine clean_exit ( inp , adi , thd , time , dt , x , y , z ,    &
!                          dx_i , dy_i , dz_i , T , W_i , cp , ha , v )
!
!
!    type (inp_type) , intent (in)                                  :: inp  !< input derived type
!    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
!    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
!    real (dp) , intent (in)                                        :: time !< time
!    real (dp) , intent (in)                                        :: dt   !< time step
!    real (dp) , allocatable , dimension (:) , intent (in)          :: x    !< x-coordinate array
!    real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
!    real (dp) , allocatable , dimension (:) , intent (in)          :: z    !< z-coordinate array
!    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i !< inverted dx array
!    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i !< inverted dy array
!    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i !< inverted dz array
!    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
!    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
!    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
!    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
!    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array
!
!
!    integer (ip) :: i , j , k , l
!
!
!    if ( ndim < ndimmax ) return ! only in 3D problems
!
!
!    if ( rank == rank_default ) write (*,*) 'applying clean_exit ...'
!
!
!    if ( neigh (E) == MPI_PROC_NULL ) then
!
!
!       do l = 1 , nv
!          do k = sz , ez
!             do j = sy , ey
!                do i = 1 , 11
!                   v (ex-11+i,j,k,l) = v (ex-11,j,k,l)
!                end do
!             end do
!          end do
!       end do
!
!
!    end if
!
!
!    call upd_prim_var_domain ( thd , T , W_i , cp , ha , v )
!    call comm_cons (v)
!    call upd_boundaries ( inp , adi , thd , time , dt , x , y , z ,    &
!                          dx_i , dy_i , dz_i , T , W_i , cp , ha , v )
!    call upd_prim_var_ghost ( thd , T , W_i , cp , ha , v )
!
!
!  end subroutine clean_exit


end module solver
