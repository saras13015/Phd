!------------------------------------------------------------------------------
! MODULE: eg_lib
!------------------------------------------------------------------------------
!> \brief EGLIB implementation.
!!
!! This module contains all the subroutines necessary to call EGLIB
!! functions and calculate the viscous fluxes.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module eg_lib


  use parameters
  use parallel , only : domain_select
  use type_thd
  use adim
  use input
  use thermodynamics
  use deriv
  use SGS_models
  use tools , only : shock_det_slv

  implicit none


contains


!> \brief Initialize EGLIB library.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine init_eglib ( inp , rweg , iweg )


    type (inp_type) , intent (in)                               :: inp  !< input derived type
    real (dp) , allocatable , dimension (:) , intent (inout)    :: rweg !< real EGLIB work array
    integer (ip) , allocatable , dimension (:) , intent (inout) :: iweg !< integer EGLIB work array


    integer (ip) , parameter :: lout = 6
    integer (ip)             :: ok , np , iflag , itls


    ! number of points (including ghost points) in each block
    np = 1


    ! evaluating iflag , itls and the lenghts of work arrays
    if ( inp % eglib_precision == 1 ) then
       itls  = 1
       iflag = 7
       liweg = nrv
       lrweg = 23 + 14*nrv + 32*nrv*nrv + &
               13 * np + 14 * np * nrv +  &
               np * nrv * nrv
    else if ( inp % eglib_precision == 2 ) then
       itls  = 2
       iflag = 7
       liweg = nrv
       lrweg = 23 + 14*nrv + 32*nrv*nrv + &
               13 * np + 21 * np * nrv +  &
               np * ( nrv*nrv + nrv*nrv + &
               ( nrv*(nrv+1)/2 ) )
    else if ( inp % eglib_precision == 3 ) then
       itls  = 2
       iflag = 7
       liweg = nrv
       lrweg = 23 + 14*nrv + 32*nrv*nrv + &
               13 * np + 21 * np * nrv +  &
               np * ( nrv*nrv + nrv*nrv + &
               ( nrv*(nrv+1)/2 ) )
    else
       write (*,*) 'ERROR: eglib_precision = ' , eglib_precision
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if

    allocate ( iweg (liweg) , rweg (lrweg) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate init_eglib')

    if ( rank == rank_default ) write (*,*) 'initializing EGLIB...'
    call egini ( np , lout , iflag , itls , rweg , lrweg , iweg , liweg )


  end subroutine init_eglib


!> \brief Calculate transport variables from EGLIB library.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine egvars ( thd , T , W_i , v , iweg , rweg , rd , mu , kpa , ct , tdr )


    type (thd_type) , intent (in)                                    :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)        :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)        :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)      :: v    !< conserved variables array
    integer (ip) , allocatable , dimension (:) , intent (inout)      :: iweg !< integer EGLIB work array
    real (dp) , allocatable , dimension (:) , intent (inout)         :: rweg !< real EGLIB work array
    real (dp) , allocatable , dimension (:,:,:,:,:) , intent (inout) :: rd   !< density times diffusion
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)     :: mu   !< dynamic viscosity
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)     :: kpa  !< bulk viscosity
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)     :: ct   !< thermal conductivity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout)   :: tdr  !< thermal diffusivity


    integer (ip) , parameter                      :: npoly = 4 ! order of the polynom
    integer (ip) , parameter                      :: io1=1 , io2=2 , io3=3 , io4=4 , io5=5 , io6=6
    integer (ip)                                  :: is , domain_id
    integer (ip)                                  :: i , j , k
    integer (ip)                                  :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                     :: eta_adim , lbda_adim , rd_adim
    real (dp)                                     :: rho_i , Wm , wrk , cp_cgs , ww_cgs
    real (dp) , dimension (npoly+1)               :: tempkelvin
    ! EGLIB variables
    real (dp)                                     :: T_eg , ww_eg , ct_eg , mu_eg , kpa_eg
    real (dp) , dimension (nrv)                   :: Y_eg , X_eg , cp_eg , tdr_eg
    real (dp) , dimension (nrv,nrv)               :: rd_eg

    real (dp)                                     :: cp0  ,& ! heat capacity at Tmin or Tmax
                                                     dcp0    ! heat capacity derivative at Tmin or Tmax


    cp_cgs   = 1.0e4_dp * thd % cpinf
    ww_cgs   = 1.0e3_dp

    eta_adim  = 1.0e-1_dp / thd % muinf
    lbda_adim = 1.0e-5_dp / thd % laminf
    rd_adim   = 1.0e-1_dp / ( thd % rhoinf * thd % Dinf (1) )


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1


                ! first part


                rho_i = 1.0_dp / v (i,j,k,1)
                Wm    = 1.0_dp / W_i (i,j,k)

                ww_eg = Wm * ww_cgs
                wrk   = rho_i * Wm

                do is = 1 , nrv

                   tempkelvin(io1) = T(i,j,k)        * thd % Tinf      ! T
                   tempkelvin(io2) = tempkelvin(io1) * tempkelvin(io1) ! T^2
                   tempkelvin(io3) = tempkelvin(io2) * tempkelvin(io1) ! T^3
                   tempkelvin(io4) = tempkelvin(io3) * tempkelvin(io1) ! T^4

                   T_eg  = tempkelvin (io1)

                   Y_eg (is) = v (i,j,k,niv+is) * rho_i
                   X_eg (is) = v (i,j,k,niv+is) * wrk * thd % Wc_i (is)

                   if ( tempkelvin (io1) < thd % Tpolymin (is) ) then ! T < Tmin linear fit from Tmin

                      tempkelvin (io1) = thd % Tpolymin (is)                 ! T
                      tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
                      tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
                      tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4

                      cp0 =   ( thd % HeatCap (io1,1,is)                                     + &
                                thd % HeatCap (io2,1,is) * tempkelvin (io1)                  + &
                                thd % HeatCap (io3,1,is) * tempkelvin (io2)                  + &
                                thd % HeatCap (io4,1,is) * tempkelvin (io3)                  + &
                                thd % HeatCap (io5,1,is) * tempkelvin (io4) )

                      dcp0 =  ( thd % HeatCap (io2,1,is)                                     + &
                                thd % HeatCap (io3,1,is) * tempkelvin (io1) * 2.0_dp         + &
                                thd % HeatCap (io4,1,is) * tempkelvin (io2) * 3.0_dp         + &
                                thd % HeatCap (io5,1,is) * tempkelvin (io3) * 4.0_dp )

                      cp_eg (is) = thd % gam2 * thd % Wc_i (is)                              * &
                              ( dcp0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymin (is) ) + cp0 )

                   else if ( tempkelvin (io1) < thd % Tpolylim (is) ) then ! T < Tlim

                      cp_eg (is) = thd % gam2 * thd % Wc_i (is)                              * &
                              ( thd % HeatCap (io1,1,is)                                     + &
                                thd % HeatCap (io2,1,is) * tempkelvin (io1)                  + &
                                thd % HeatCap (io3,1,is) * tempkelvin (io2)                  + &
                                thd % HeatCap (io4,1,is) * tempkelvin (io3)                  + &
                                thd % HeatCap (io5,1,is) * tempkelvin (io4) )

                   else if ( tempkelvin (io1) < thd % Tpolymax (is) ) then ! T < Tmax

                      cp_eg (is) = thd % gam2 * thd % Wc_i (is)                              * &
                              ( thd % HeatCap (io1,2,is)                                     + &
                                thd % HeatCap (io2,2,is) * tempkelvin (io1)                  + &
                                thd % HeatCap (io3,2,is) * tempkelvin (io2)                  + &
                                thd % HeatCap (io4,2,is) * tempkelvin (io3)                  + &
                                thd % HeatCap (io5,2,is) * tempkelvin (io4) )

                   else ! T > Tmax linear fit of cp from Tmax

                      tempkelvin (io1) = thd % Tpolymax (is)                 ! T
                      tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
                      tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
                      tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4

                      cp0 =   ( thd % HeatCap (io1,2,is)                                     + &
                                thd % HeatCap (io2,2,is) * tempkelvin (io1)                  + &
                                thd % HeatCap (io3,2,is) * tempkelvin (io2)                  + &
                                thd % HeatCap (io4,2,is) * tempkelvin (io3)                  + &
                                thd % HeatCap (io5,2,is) * tempkelvin (io4) )

                      dcp0 =  ( thd % HeatCap (io2,2,is)                                     + &
                                thd % HeatCap (io3,2,is) * tempkelvin (io1) * 2.0_dp         + &
                                thd % HeatCap (io4,2,is) * tempkelvin (io2) * 3.0_dp         + &
                                thd % HeatCap (io5,2,is) * tempkelvin (io3) * 4.0_dp )

                      cp_eg (is) = thd % gam2 * thd % Wc_i (is)                              * &
                              ( dcp0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymax (is) ) + cp0 )

                   end if

                end do

                cp_eg (:) = cp_eg (:) * cp_cgs


                ! second part


                call egspar ( T_eg , X_eg , Y_eg , cp_eg , rweg , iweg )
                if ( eglib_precision == 1 ) then
                   call egsl2   ( T_eg , Y_eg , rweg , iweg , ct_eg )
                   call egse2   ( T_eg , Y_eg , rweg , mu_eg )
                   call egsk2   ( T_eg , Y_eg , rweg , kpa_eg )
                   ! replace egslct1 by egslct2 because of problems in 2D and 3D
                   call egslct2 ( T_eg , Y_eg , rweg , iweg , ct_eg , tdr_eg )
                   !
                   call egsdr1  ( T_eg , Y_eg , ww_eg , rweg , rd_eg )
                else if ( eglib_precision == 2 ) then
                   call egsl4   ( T_eg , Y_eg , rweg , iweg , ct_eg )
                   call egse3   ( T_eg , Y_eg , rweg , mu_eg )
                   call egsk5   ( T_eg , Y_eg , rweg , kpa_eg )
                   call egslct3 ( T_eg , Y_eg , rweg , iweg , ct_eg , tdr_eg )
                   call egsdr1  ( T_eg , Y_eg , ww_eg , rweg , rd_eg )
                else if ( eglib_precision == 3 ) then
                   call egsl5   ( T_eg , Y_eg , rweg , iweg , ct_eg )
                   call egse4   ( T_eg , Y_eg , rweg , mu_eg )
                   call egsk6   ( T_eg , Y_eg , rweg , kpa_eg )
                   call egslct4 ( T_eg , Y_eg , rweg , iweg , ct_eg , tdr_eg )
                   call egsdr2  ( T_eg , Y_eg , ww_eg , rweg , rd_eg )
                end if

                ct (i,j,k)     = ct_eg * lbda_adim
                mu (i,j,k)     = mu_eg * eta_adim
                kpa (i,j,k)    = kpa_eg * eta_adim
                tdr (i,j,k,:)  = tdr_eg
                rd (i,j,k,:,:) = rd_eg (:,:) * rd_adim


             end do
          end do
       end do

    end do


  end subroutine egvars


!> \brief 1D viscous flux.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_eg_x ( thd , adi , grid , T , W_i , cp , ha , v , &
                             rweg , iweg , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:) , intent (inout)       :: rweg !< real EGLIB work array
    integer (ip) , allocatable , dimension (:) , intent (inout)    :: iweg !< integer EGLIB work array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< viscous flux


    real (dp) , parameter                           :: twothirds = 2.0_dp / 3.0_dp
    integer (ip)                                    :: nfluid
    integer (ip)                                    :: ok , domain_id , i , j , k , l , m
    integer (ip)                                    :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                       :: rmu , rmk , rvd0 , rlam0 , rlamcp0 , bulk
    real (dp)                                       :: rho_i , Ps , Ts , uxs
    real (dp) , dimension (nrv)                     :: d_x
    real (dp)                                       :: wrk1 , wrk2
    ! viscous variables (eglib)
    real (dp) , dimension (:,:,:) , allocatable     :: mu , kpa , ct
    real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
    real (dp) , dimension (:,:,:,:,:) , allocatable :: rd
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: ux , P
    real (dp) , dimension (:,:,:,:) , allocatable   :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: q_x  , &
                                                       sx_x
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: wrk_x
    real (dp) , dimension (:,:,:,:) , allocatable   :: wrka_x
    ! array containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable   :: fd


    nfluid  = ndim * ( ndim + 2 )
    rvd0    = adi % sqgmr / adi % Sc
    rlam0   = adi % sqgmr * adi % ggmopr
    rlamcp0 = adi % sqgmr / adi % Pr


    allocate ( rd  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv,nrv) , &
               tdr (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)     , &
               mu  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               kpa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               ct  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_x')


    call egvars ( thd , T , W_i , v , iweg , rweg , rd , mu , kpa , ct , tdr )


    ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_x 2')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                P (i,j,k)  = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             end do
          end do
       end do

    end do


    ! passive scalar equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do l = nrv+1,nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      Xa (i,j,k,l) = v (i,j,k,niv+l) / v (i,j,k,1)
                   end do
                end do
             end do
          end do

       end do

    end if


    call dx_fixed1 ( grid % dx_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,1) )

    call dx_fixed1 ( grid % dx_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) )

    call dx_fixed1 ( grid % dx_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,3) ) ! nfluid = 3 in 1D

    do l = 1 , nrv+npv
       call dx_fixed1 ( grid % dx_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l) )
    end do


    deallocate ( ux )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species
    ! and/or non-pressure, you must specity the temperature and/or the
    ! velocity and/or the species and/or the pressure normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       fd (sx,sy:ey,sz:ez,1) = 0.0_dp ! normal velocity gradient
       fd (sx,sy:ey,sz:ez,2) = 0.0_dp ! normal temperature gradient
       fd (sx,sy:ey,sz:ez,3) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       fd (ex,sy:ey,sz:ez,1) = 0.0_dp ! normal velocity gradient
       fd (ex,sy:ey,sz:ez,2) = 0.0_dp ! normal temperature gradient
       fd (ex,sy:ey,sz:ez,3) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l) = 0.0_dp ! normal species gradient
       end do
    end if


    ! communicate the first derivative
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( sx_x  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_x 3')


    ! Reynolds tensor
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rmu = - adi % sqgmr * mu (i,j,k)
                rmk = - adi % sqgmr * kpa (i,j,k)

                bulk = rmk * fd (i,j,k,1)

                sx_x (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) )

             end do
          end do
       end do

    end do

    deallocate ( mu , kpa )


    ! rYV_x
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                Ps    = 1.0_dp / P (i,j,k)
                Ts    = 1.0_dp / T (i,j,k)

                wrk1 = fd (i,j,k,2) * Ts

                wrk2 = fd (i,j,k,3) * Ps

                do l = 1 , nrv
                   ! d_x
                   d_x (l) = fd (i,j,k,nfluid+ndim*l) +                 &
                           ( Xa (i,j,k,l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk2
                end do

                do l = 1 , nrv

                   fd (i,j,k,nfluid+ndim*l) = 0.0_dp ! rYV_x

                   do m = 1 , nrv
                      fd (i,j,k,nfluid+ndim*l) = fd (i,j,k,nfluid+ndim*l) +                      &
                                                 rd (i,j,k,l,m) *                                &
                                               ( d_x (m) + Xa (i,j,k,m) * tdr (i,j,k,m) * wrk1 )
                   end do

                   fd (i,j,k,nfluid+ndim*l) = - rvd0 * fd (i,j,k,nfluid+ndim*l)

                end do

             end do
          end do
       end do

    end do


    ! passive scalar equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   P (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) ! P = -rlamcp
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l)
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( P , Xa , rd )


    ! components of the viscous stresses
    allocate ( q_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_x 4')


    ! q_x
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                uxs   = v (i,j,k,2) * rho_i

                wrk1 = -rlam0 * ct (i,j,k)

                q_x (i,j,k) = 0.0_dp

                do l = 1 , nrv
                   q_x (i,j,k) = q_x (i,j,k) + fd (i,j,k,nfluid+ndim*l) *                    &
                               ( ha (i,j,k,l) + T (i,j,k) * thd % Wc_i (l) * tdr (i,j,k,l) )
                end do

                q_x (i,j,k) = q_x (i,j,k) + wrk1 * fd (i,j,k,2) + &
                              uxs * sx_x (i,j,k)

             end do
          end do
       end do

    end do


    deallocate ( ct , tdr )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_x 4')


    call dx ( grid % dx_i , sx_x , wrk_x )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,2) = fl (i,j,k,2) + &
                          ( wrk_x (i,j,k) )
          end do
       end do
    end do

    deallocate (sx_x)


    call dx ( grid % dx_i , q_x , wrk_x )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,5) = fl (i,j,k,5) + &
                          ( wrk_x (i,j,k) )
          end do
       end do
    end do

    deallocate (q_x)
    deallocate (wrk_x)


    allocate ( wrka_x (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_x 5')

    do l = 1 , nrv+npv
       call dx_fixed2 ( grid % dx_i                                                , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                               )
    end do

    do l = 1 , nrv+npv
       do k = sz , ez
          do j = sy , ey
             do i = sx , ex
                fl (i,j,k,niv+l) = fl (i,j,k,niv+l) + &
                                 ( wrka_x (i,j,k,l) )
             end do
          end do
       end do
    end do

    deallocate (fd)
    deallocate (wrka_x)


  end subroutine viscflux_eg_x


!> \brief 1D viscous flux with LES models.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_eg_LES_x ( inp , thd , adi , grid , T , W_i , cp , ha , v , &
                                 rweg , iweg , mu_SGS , fl )


    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid   !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T      !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp     !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha     !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v      !< conserved variables array
    real (dp) , allocatable , dimension (:) , intent (inout)       :: rweg   !< real EGLIB work array
    integer (ip) , allocatable , dimension (:) , intent (inout)    :: iweg   !< integer EGLIB work array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mu_SGS !< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl     !< viscous flux


    real (dp) , parameter                           :: twothirds = 2.0_dp / 3.0_dp
    integer (ip)                                    :: nfluid , nSGS
    integer (ip)                                    :: ok , domain_id , i , j , k , l , m
    integer (ip)                                    :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                       :: rmu , rmu_SGS , rmk , rvd0 , rlam0 , rlam , rlamcp0 , bulk , &
                                                       Sc_SGS_i , Pr_SGS_i , gm1_g , ct_SGS , dm_SGS
    real (dp)                                       :: rho_i , Ps , Ts , uxs
    real (dp) , dimension (nrv)                     :: d_x
    real (dp)                                       :: wrk , wrk1 , wrk2
    ! viscous variables (eglib)
    real (dp) , dimension (:,:,:) , allocatable     :: mu , kpa , ct
    real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
    real (dp) , dimension (:,:,:,:,:) , allocatable :: rd
    ! SGS additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: tau_iso_SGS
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: ux , P
    real (dp) , dimension (:,:,:,:) , allocatable   :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: q_x  , &
                                                       sx_x
    real (dp) , dimension (:,:,:) , allocatable     :: sx_x_SGS
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: wrk_x
    real (dp) , dimension (:,:,:,:) , allocatable   :: wrka_x
    ! array containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable   :: fd


    nfluid  = ndim * ( ndim + 2 )
    nSGS    = nfluid + ndim * ( nrv + npv )

    rvd0    = adi % sqgmr / adi % Sc     ! M * sqrt(gamma) / ( Re * Sc )
    rlam0   = adi % sqgmr * adi % ggmopr ! M * sqrt(gamma) * gamma / ( Re * (gamma-1) * Pr )
    rlamcp0 = adi % sqgmr / adi % Pr     ! M * sqrt(gamma) / ( Re * Pr )

    Sc_SGS_i = 1.0_dp / inp % Sc_SGS
    Pr_SGS_i = 1.0_dp / inp % Pr_SGS
    gm1_g    = adi % gm1 / adi % gamma


    allocate ( rd  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv,nrv) , &
               tdr (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)     , &
               mu  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               kpa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               ct  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_x')


    call egvars ( thd , T , W_i , v , iweg , rweg , rd , mu , kpa , ct , tdr )


    ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_x 2')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                P (i,j,k)  = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             end do
          end do
       end do

    end do


    ! passive scalar equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do l = nrv+1,nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      Xa (i,j,k,l) = v (i,j,k,niv+l) / v (i,j,k,1)
                   end do
                end do
             end do
          end do

       end do

    end if


    call dx_fixed1 ( grid % dx_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,1) ) ! du/dx = fd (:,:,:,1)

    call dx_fixed1 ( grid % dx_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) ) ! dT/dx = fd (:,:,:,2)

    call dx_fixed1 ( grid % dx_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,3) ) ! dp/dx = fd (:,:,:,3) nfluid = 3 in 1D

    do l = 1 , nrv+npv
       call dx_fixed1 ( grid % dx_i                                            , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)             , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l) ) ! dXa/dx
    end do


    ! SGS: W_i derivatives need to be communicated
    call dx_fixed1 ( grid % dx_i , W_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nSGS+1) ) ! dW_i/dx

    ! SGS: shock detector criteria needs to be communicated
    call shock_det_slv ( v , T , W_i , fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) )


    deallocate ( ux )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species
    ! and/or non-pressure, you must specity the temperature and/or the
    ! velocity and/or the species and/or the pressure normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       if ( bc (W) == symmetryplane ) then
          fd (sx,sy:ey,sz:ez,1) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx,sy:ey,sz:ez,2) = 0.0_dp ! normal temperature gradient
       fd (sx,sy:ey,sz:ez,3) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l) = 0.0_dp ! normal species gradient
       end do
       fd (sx,sy:ey,sz:ez,nSGS+1) = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       if ( bc (E) == symmetryplane ) then
          fd (ex,sy:ey,sz:ez,1) = 0.0_dp ! normal velocity gradient
       end if
       fd (ex,sy:ey,sz:ez,2) = 0.0_dp ! normal temperature gradient
       fd (ex,sy:ey,sz:ez,3) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l) = 0.0_dp ! normal species gradient
       end do
       fd (ex,sy:ey,sz:ez,nSGS+1) = 0.0_dp ! normal mixture molecular weight gradient
    end if


    ! communicate the first derivative
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( sx_x         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_x_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               tau_iso_SGS  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_x 3')

    ! components of the viscous stresses
    allocate ( q_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_x 4')


    ! SGS: viscosity
    do domain_id = -ndim , ndim
       call mu_SGS_selector ( inp , adi , domain_id , grid , v , fd , mu_SGS )
    end do

    ! SGS: tau_iso_SGS derivatives need to be COMMUNICATED
    if ( inp % tau_iso_SGS_switch ) then
       do domain_id = -ndim , ndim
          call tau_iso_SGS_selector ( inp , domain_id , grid , v , fd , tau_iso_SGS )
       end do
    else
       tau_iso_SGS (:,:,:) = 0.0_dp
    end if


    ! Reynolds tensor
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rmu          = - adi % sqgmr * mu (i,j,k)
                rmu_SGS      = - adi % sqgmr * mu_SGS (i,j,k) / rmu
                rmk          = - adi % sqgmr * kpa (i,j,k)

                bulk = rmk * fd (i,j,k,1)

                sx_x (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) )
                sx_x_SGS (i,j,k) = bulk + ( sx_x (i,j,k) - bulk ) * ( 1.0_dp + rmu_SGS ) + tau_iso_SGS (i,j,k)

             end do
          end do
       end do

    end do

    deallocate ( mu , kpa )
    deallocate ( sx_x )

    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)

                uxs   = v (i,j,k,2) * rho_i

                Ps    = 1.0_dp / P (i,j,k)
                Ts    = 1.0_dp / T (i,j,k)

                wrk1 = fd (i,j,k,2) * Ts

                wrk2 = fd (i,j,k,3) * Ps

                dm_SGS = adi % Sc * mu_SGS (i,j,k) * rho_i * Sc_SGS_i
                ct_SGS = adi % Pr * mu_SGS (i,j,k) * cp (i,j,k) * Pr_SGS_i

                ! Diffusive mass flux
                ! part : d_{\beta,j} = dXb/dxj + (Xb-Yb)/p * dp/dxj
                do l = 1 , nrv ! \beta loop
                   wrk = Xa (i,j,k,l) - v (i,j,k,niv+l) * rho_i

                   d_x (l) = fd (i,j,k,nfluid+ndim*l) + wrk * wrk2
                end do

                ! Heat flux
                q_x (i,j,k) = 0.0_dp

                do l = 1 , nrv ! \alpha loop

                   ! Diffusive mass flux
                   ! rho * Ya * Vaj = - \sum_{\beta=1,N} rho * D_ab * ( d_bj + Xb * \chi_b / T * dT/dxj )
                   wrk2 = 0.0_dp ! rYV_x

                   do m = 1 , nrv ! \beta loop

                      wrk  = Xa (i,j,k,m) * tdr (i,j,k,m)
                      wrk2 = wrk2 + rd (i,j,k,l,m) * ( d_x (m) + wrk * wrk1 )

                   end do

                   ! + SGS diffusive mass flux
                   ! rho * D_sgs * dYa/dxj = rho * D_SGS * Wa * ( W_i * dXa/dxj + Xa * dW_i/dxj )
                   wrk1 = v(i,j,k,1) * dm_SGS * thd % Wc (l)

                   fd (i,j,k,nfluid+ndim*l) = - rvd0 * ( wrk2 + wrk1 * ( W_i (i,j,k) * fd (i,j,k,nfluid+ndim*l) + &
                                                                         Xa (i,j,k,l) * fd (i,j,k,nSGS+1)        ))

                   ! Heat flux part 1
                   ! \sum_{\alpha=1,N} rho * Ya * Vaj * ( ha + R*T*\chi_a / Wa )
                   wrk1 = ha (i,j,k,l) + T (i,j,k) * thd % Wc_i (l) * tdr (i,j,k,l)
                   q_x (i,j,k) = q_x (i,j,k) + wrk2 * wrk1

                end do

                ! Heat flux part 2 (with SGS terms)
                ! qj = - lambda * dT/dxj                                             &
                !      + \sum_{\alpha=1,N} rho * Ya * Vaj * ( ha + R*T*\chi_a / Wa ) &
                !      + tau_ij * u_i                                                &
                !      - lambda_SGS * dT/dxj                                         & (SGS term)
                !      - tau_ij_SGS * u_i                                              (SGS term)
                rlam = - rlam0 * ( ct (i,j,k) + ct_SGS )

                q_x (i,j,k) = q_x (i,j,k) +                  &
                              rlam * fd (i,j,k,10) +         &
                              uxs * sx_x_SGS (i,j,k)

             end do
          end do
       end do

    end do


    deallocate ( tdr )
    deallocate ( tau_iso_SGS )
    deallocate ( Xa , rd )


    ! passive scalar equations: viscous term
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   dm_SGS = adi % Sc * mu_SGS (i,j,k) / v (i,j,k,1) * Sc_SGS_i
                   P (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) & ! > laminar/molecular term ( P = -rlamcp )
                               - rvd0 * v (i,j,k,1) * dm_sgs         ! > SGS term for LES
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l)
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( ct , P )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_x 5')


    ! viscous flux in rho*u equation -> fl (:,:,:,2)
    call dx ( grid % dx_i , sx_x_SGS , wrk_x )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,2) = fl (i,j,k,2) + &
                          ( wrk_x (i,j,k) )
          end do
       end do
    end do

    deallocate (sx_x_SGS)


    ! viscous flux in rho*e_tot equation -> fl (:,:,:,5)
    call dx ( grid % dx_i , q_x , wrk_x )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,5) = fl (i,j,k,5) + &
                          ( wrk_x (i,j,k) )
          end do
       end do
    end do

    deallocate (q_x)
    deallocate (wrk_x)


    allocate ( wrka_x (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_x 6')


    ! viscous flux in rho*Y_alpha equation -> fl (:,:,:,niv+1:niv+npv+nvv)
    do l = 1 , nrv+npv
       call dx_fixed2 ( grid % dx_i                                                , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                               )
    end do

    do l = 1 , nrv+npv
       do k = sz , ez
          do j = sy , ey
             do i = sx , ex
                fl (i,j,k,niv+l) = fl (i,j,k,niv+l) + &
                                 ( wrka_x (i,j,k,l) )
             end do
          end do
       end do
    end do

    deallocate (fd)
    deallocate (wrka_x)


  end subroutine viscflux_eg_LES_x


!> \brief 2D viscous flux.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_eg_xy ( thd , adi , grid , T , W_i , cp , ha , v , &
                              rweg , iweg , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:) , intent (inout)       :: rweg !< real EGLIB work array
    integer (ip) , allocatable , dimension (:) , intent (inout)    :: iweg !< integer EGLIB work array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< viscous flux


    real (dp) , parameter                           :: twothirds = 2.0_dp / 3.0_dp
    integer (ip)                                    :: nfluid
    integer (ip)                                    :: ok , domain_id , i , j , k , l , m
    integer (ip)                                    :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                       :: rmu , rmk , rvd0 , rlam0 , rlamcp0 , bulk
    real (dp)                                       :: rho_i , Ps , Ts , uxs , vys
    real (dp) , dimension (nrv)                     :: d_x , d_y
    real (dp)                                       :: wrk1 , wrk2 , wrk3 , wrk4
    ! viscous variables (eglib)
    real (dp) , dimension (:,:,:) , allocatable     :: mu , kpa , ct
    real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
    real (dp) , dimension (:,:,:,:,:) , allocatable :: rd
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: ux , vy , P
    real (dp) , dimension (:,:,:,:) , allocatable   :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: q_x  , q_y  , &
                                                       sx_x , sx_y , &
                                                       sy_y
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: wrk_x , wrk_y
    real (dp) , dimension (:,:,:,:) , allocatable   :: wrka_x , wrka_y
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable   :: fd


    nfluid  = ndim * ( ndim + 2 )
    rvd0    = adi % sqgmr / adi % Sc
    rlam0   = adi % sqgmr * adi % ggmopr
    rlamcp0 = adi % sqgmr / adi % Pr


    allocate ( rd  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv,nrv) , &
               tdr (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)     , &
               mu  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               kpa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               ct  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xy')


    call egvars ( thd , T , W_i , v , iweg , rweg , rd , mu , kpa , ct , tdr )


    ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xy 2')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
                P (i,j,k)  = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             end do
          end do
       end do

    end do


    ! passive scalar equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do l = nrv+1,nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      Xa (i,j,k,l) = v (i,j,k,niv+l) / v (i,j,k,1)
                   end do
                end do
             end do
          end do

       end do

    end if


    call dx_fixed1 ( grid % dx_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,1) )
    call dy_fixed1 ( grid % dy_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) )

    call dx_fixed1 ( grid % dx_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,3) )
    call dy_fixed1 ( grid % dy_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,4) )

    call dx_fixed1 ( grid % dx_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,5) )
    call dy_fixed1 ( grid % dy_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,6) )

    call dx_fixed1 ( grid % dx_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,7) )
    call dy_fixed1 ( grid % dy_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,8) ) ! nfluid = 8 in 2D

    do l = 1 , nrv+npv
       call dx_fixed1 ( grid % dx_i                                            , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)             , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) )
       call dy_fixed1 ( grid % dy_i                                            , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)             , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) )
    end do

    deallocate ( ux , vy )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species
    ! and/or non-pressure, you must specity the temperature and/or the
    ! velocity and/or the species and/or the pressure normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       fd (sx,sy:ey,sz:ez,1) = 0.0_dp ; fd (sx,sy:ey,sz:ez,3) = 0.0_dp ! normal velocity gradient
       fd (sx,sy:ey,sz:ez,5) = 0.0_dp ! normal temperature gradient
       fd (sx,sy:ey,sz:ez,7) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       fd (ex,sy:ey,sz:ez,1) = 0.0_dp ; fd (ex,sy:ey,sz:ez,3) = 0.0_dp ! normal velocity gradient
       fd (ex,sy:ey,sz:ez,5) = 0.0_dp ! normal temperature gradient
       fd (ex,sy:ey,sz:ez,7) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (S) == symmetryplane .or. bc (S) == adiabaticwall ) then
       fd (sx:ex,sy,sz:ez,2) = 0.0_dp ; fd (sx:ex,sy,sz:ez,4) = 0.0_dp ! normal velocity gradient
       fd (sx:ex,sy,sz:ez,6) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,sy,sz:ez,8) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy,sz:ez,nfluid+ndim*l ) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (N) == symmetryplane .or. bc (N) == adiabaticwall ) then
       fd (sx:ex,ey,sz:ez,2) = 0.0_dp ; fd (sx:ex,ey,sz:ez,4) = 0.0_dp ! normal velocity gradient
       fd (sx:ex,ey,sz:ez,6) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,ey,sz:ez,8) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,ey,sz:ez,nfluid+ndim*l ) = 0.0_dp ! normal species gradient
       end do
    end if


    ! communicate the first derivative
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( sx_x  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               sx_y  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               sy_y  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xy 3')


    ! Reynolds tensor
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rmu = - adi % sqgmr * mu (i,j,k)
                rmk = - adi % sqgmr * kpa (i,j,k)

                bulk = rmk * ( fd (i,j,k,1) + fd (i,j,k,4) )

                sx_x (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) - fd (i,j,k,4) )

                sx_y (i,j,k) = rmu * ( fd (i,j,k,3) + fd (i,j,k,2) )

                sy_y (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,4) + fd (i,j,k,4) - fd (i,j,k,1) )

             end do
          end do
       end do

    end do

    deallocate ( mu , kpa )


    ! rYV_x , rYV_y
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                Ps    = 1.0_dp / P (i,j,k)
                Ts    = 1.0_dp / T (i,j,k)

                wrk1 = fd (i,j,k,5) * Ts
                wrk2 = fd (i,j,k,6) * Ts

                wrk3 = fd (i,j,k,7) * Ps
                wrk4 = fd (i,j,k,8) * Ps

                do l = 1 , nrv
                   ! d_x
                   d_x (l) = fd (i,j,k,nfluid+ndim*l-1) +               &
                           ( Xa (i,j,k,l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk3
                   ! d_y
                   d_y (l) = fd (i,j,k,nfluid+ndim*l  ) +               &
                           ( Xa (i,j,k,l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk4
                end do

                do l = 1 , nrv

                   fd (i,j,k,nfluid+ndim*l-1) = 0.0_dp ! rYV_x
                   fd (i,j,k,nfluid+ndim*l  ) = 0.0_dp ! rYV_y

                   do m = 1 , nrv
                      fd (i,j,k,nfluid+ndim*l-1) = fd (i,j,k,nfluid+ndim*l-1) +                    &
                                                   rd (i,j,k,l,m) *                                &
                                                 ( d_x (m) + Xa (i,j,k,m) * tdr (i,j,k,m) * wrk1 )
                      fd (i,j,k,nfluid+ndim*l  ) = fd (i,j,k,nfluid+ndim*l  ) +                    &
                                                   rd (i,j,k,l,m) *                                &
                                                 ( d_y (m) + Xa (i,j,k,m) * tdr (i,j,k,m) * wrk2 )
                   end do

                   fd (i,j,k,nfluid+ndim*l-1) = - rvd0 * fd (i,j,k,nfluid+ndim*l-1)
                   fd (i,j,k,nfluid+ndim*l  ) = - rvd0 * fd (i,j,k,nfluid+ndim*l  )

                end do

             end do
          end do
       end do

    end do


    ! passive scalar equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   P (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) ! P = -rlamcp
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l-1) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l-1)
                      fd (i,j,k,nfluid+ndim*l  ) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l  )
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( P , Xa , rd )


    ! components of the viscous stresses
    allocate ( q_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               q_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xy 4')


    ! q_x , q_y
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                uxs   = v (i,j,k,2) * rho_i
                vys   = v (i,j,k,3) * rho_i

                wrk1 = -rlam0 * ct (i,j,k)

                q_x (i,j,k) = 0.0_dp
                q_y (i,j,k) = 0.0_dp

                do l = 1 , nrv
                   q_x (i,j,k) = q_x (i,j,k) + fd (i,j,k,nfluid+ndim*l-1) *                  &
                               ( ha (i,j,k,l) + T (i,j,k) * thd % Wc_i (l) * tdr (i,j,k,l) )
                   q_y (i,j,k) = q_y (i,j,k) + fd (i,j,k,nfluid+ndim*l  ) *                  &
                               ( ha (i,j,k,l) + T (i,j,k) * thd % Wc_i (l) * tdr (i,j,k,l) )
                end do

                q_x (i,j,k) = q_x (i,j,k) + wrk1 * fd (i,j,k,5) +     &
                              uxs * sx_x (i,j,k) + vys * sx_y (i,j,k)
                q_y (i,j,k) = q_y (i,j,k) + wrk1 * fd (i,j,k,6) +     &
                              uxs * sx_y (i,j,k) + vys * sy_y (i,j,k)

             end do
          end do
       end do

    end do


    deallocate ( ct , tdr )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               wrk_y (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xy 4')


    call dx ( grid % dx_i , sx_x , wrk_x )
    call dy ( grid % dy_i , sx_y , wrk_y )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,2) = fl (i,j,k,2) +                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) )
          end do
       end do
    end do

    deallocate (sx_x)


    call dx ( grid % dx_i , sx_y , wrk_x )
    call dy ( grid % dy_i , sy_y , wrk_y )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,3) = fl (i,j,k,3) +                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) )
          end do
       end do
    end do

    deallocate ( sx_y , sy_y )


    call dx ( grid % dx_i , q_x , wrk_x )
    call dy ( grid % dy_i , q_y , wrk_y )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,5) = fl (i,j,k,5) +                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) )
          end do
       end do
    end do

    deallocate ( q_x , q_y  )
    deallocate ( wrk_x , wrk_y )


    allocate ( wrka_x (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               wrka_y (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xy 5')

    do l = 1 , nrv+npv
       call dx_fixed2 ( grid % dx_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                                 )
       call dy_fixed2 ( grid % dy_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) , &
                        wrka_y (sx:ex,sy:ey,sz:ez,l)                                 )
    end do

    do l = 1 , nrv+npv
       do k = sz , ez
          do j = sy , ey
             do i = sx , ex
                fl (i,j,k,niv+l) = fl (i,j,k,niv+l) +                    &
                                 ( wrka_x (i,j,k,l) + wrka_y (i,j,k,l) )
             end do
          end do
       end do
    end do

    deallocate (fd)
    deallocate ( wrka_x , wrka_y )


  end subroutine viscflux_eg_xy


!> \brief 2D viscous flux with LES models.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_eg_LES_xy ( inp , thd , adi , grid , T , W_i , cp , ha , v , &
                                  rweg , iweg , mu_SGS , fl )


    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid   !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T      !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp     !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha     !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v      !< conserved variables array
    real (dp) , allocatable , dimension (:) , intent (inout)       :: rweg   !< real EGLIB work array
    integer (ip) , allocatable , dimension (:) , intent (inout)    :: iweg   !< integer EGLIB work array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mu_SGS !< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl     !< viscous flux


    real (dp) , parameter                           :: twothirds = 2.0_dp / 3.0_dp
    integer (ip)                                    :: nfluid , nSGS
    integer (ip)                                    :: ok , domain_id , i , j , k , l , m
    integer (ip)                                    :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                       :: rmu , rmu_SGS , rmk , rvd0 , rlam0 , rlam , rlamcp0 , bulk , &
                                                       Sc_SGS_i , Pr_SGS_i , gm1_g , ct_SGS , dm_SGS
    real (dp)                                       :: rho_i , Ps , Ts , uxs , vys
    real (dp) , dimension (nrv)                     :: d_x , d_y
    real (dp)                                       :: wrk , wrk1 , wrk2 , wrk3 , wrk4
    ! viscous variables (eglib)
    real (dp) , dimension (:,:,:) , allocatable     :: mu , kpa , ct
    real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
    real (dp) , dimension (:,:,:,:,:) , allocatable :: rd
    ! SGS additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: tau_iso_SGS
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: ux , vy , P
    real (dp) , dimension (:,:,:,:) , allocatable   :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: q_x  , q_y  , &
                                                       sx_x , sx_y , &
                                                       sy_y
    real (dp) , dimension (:,:,:) , allocatable     :: sx_x_SGS , sx_y_SGS , &
                                                       sy_y_SGS
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: wrk_x , wrk_y
    real (dp) , dimension (:,:,:,:) , allocatable   :: wrka_x , wrka_y
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable   :: fd


    nfluid  = ndim * ( ndim + 2 )
    nSGS    = nfluid + ndim * ( nrv + npv )

    rvd0    = adi % sqgmr / adi % Sc     ! M * sqrt(gamma) / ( Re * Sc )
    rlam0   = adi % sqgmr * adi % ggmopr ! M * sqrt(gamma) * gamma / ( Re * (gamma-1) * Pr )
    rlamcp0 = adi % sqgmr / adi % Pr     ! M * sqrt(gamma) / ( Re * Pr )

    Sc_SGS_i = 1.0_dp / inp % Sc_SGS
    Pr_SGS_i = 1.0_dp / inp % Pr_SGS
    gm1_g    = adi % gm1 / adi % gamma


    allocate ( rd  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv,nrv) , &
               tdr (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)     , &
               mu  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               kpa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               ct  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xy')


    call egvars ( thd , T , W_i , v , iweg , rweg , rd , mu , kpa , ct , tdr )


    ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xy 2')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
                P (i,j,k)  = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             end do
          end do
       end do

    end do


    ! passive scalar equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do l = nrv+1,nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      Xa (i,j,k,l) = v (i,j,k,niv+l) / v (i,j,k,1)
                   end do
                end do
             end do
          end do

       end do

    end if


    call dx_fixed1 ( grid % dx_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,1) ) ! du/dx = fd (:,:,:,1)
    call dy_fixed1 ( grid % dy_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) ) ! du/dy = fd (:,:,:,2)

    call dx_fixed1 ( grid % dx_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,3) ) ! dv/dx = fd (:,:,:,3)
    call dy_fixed1 ( grid % dy_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,4) ) ! dv/dy = fd (:,:,:,4)

    call dx_fixed1 ( grid % dx_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,5) ) ! dT/dx = fd (:,:,:,5)
    call dy_fixed1 ( grid % dy_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,6) ) ! dT/dy = fd (:,:,:,6)

    call dx_fixed1 ( grid % dx_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,7) ) ! dp/dx = fd (:,:,:,7)
    call dy_fixed1 ( grid % dy_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)   , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,8) ) ! dp/dy = fd (:,:,:,8) nfluid = 8 in 2D

    do l = 1 , nrv+npv
       call dx_fixed1 ( grid % dx_i                                              , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)               , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) ) ! dXa/dx
       call dy_fixed1 ( grid % dy_i                                              , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)               , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) ) ! dXa/dy
    end do


    ! SGS: W_i derivatives need to be communicated
    call dx_fixed1 ( grid % dx_i , W_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nSGS+1) ) ! dW_i/dx
    call dy_fixed1 ( grid % dy_i , W_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nSGS+2) ) ! dW_i/dy

    ! SGS: shock detector criteria needs to be communicated
    call shock_det_slv ( v , T , W_i , fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) )


    deallocate ( ux , vy )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species
    ! and/or non-pressure, you must specity the temperature and/or the
    ! velocity and/or the species and/or the pressure normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       if ( bc (W) == symmetryplane ) then
          fd (sx,sy:ey,sz:ez,1) = 0.0_dp
          fd (sx,sy:ey,sz:ez,3) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx,sy:ey,sz:ez,5) = 0.0_dp ! normal temperature gradient
       fd (sx,sy:ey,sz:ez,7) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
       fd (sx,sy:ey,sz:ez,nSGS+1) = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       if ( bc (E) == symmetryplane ) then
          fd (ex,sy:ey,sz:ez,1) = 0.0_dp
          fd (ex,sy:ey,sz:ez,3) = 0.0_dp ! normal velocity gradient
       end if
       fd (ex,sy:ey,sz:ez,5) = 0.0_dp ! normal temperature gradient
       fd (ex,sy:ey,sz:ez,7) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
       fd (ex,sy:ey,sz:ez,nSGS+1) = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (S) == symmetryplane .or. bc (S) == adiabaticwall ) then
       if ( bc (S) == symmetryplane ) then
          fd (sx:ex,sy,sz:ez,2) = 0.0_dp
          fd (sx:ex,sy,sz:ez,4) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,sy,sz:ez,6) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,sy,sz:ez,8) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy,sz:ez,nfluid+ndim*l ) = 0.0_dp ! normal species gradient
       end do
       fd (sx:ex,sy,sz:ez,nSGS+2) = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (N) == symmetryplane .or. bc (N) == adiabaticwall ) then
       if ( bc (N) == symmetryplane ) then
          fd (sx:ex,ey,sz:ez,2) = 0.0_dp
          fd (sx:ex,ey,sz:ez,4) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,ey,sz:ez,6) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,ey,sz:ez,8) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,ey,sz:ez,nfluid+ndim*l ) = 0.0_dp ! normal species gradient
       end do
       fd (sx:ex,ey,sz:ez,nSGS+2) = 0.0_dp ! normal mixture molecular weight gradient
    end if


    ! communicate the first derivative
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( sx_x         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_y         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sy_y         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_x_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_y_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sy_y_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               tau_iso_SGS  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xy 3')

    ! components of the viscous stresses
    allocate ( q_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               q_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xy 4')


    ! SGS: viscosity
    do domain_id = -ndim , ndim
       call mu_SGS_selector ( inp , adi , domain_id , grid , v , fd , mu_SGS )
    end do

    ! SGS: tau_iso_SGS derivatives need to be COMMUNICATED
    if ( inp % tau_iso_SGS_switch ) then
       do domain_id = -ndim , ndim
          call tau_iso_SGS_selector ( inp , domain_id , grid , v , fd , tau_iso_SGS )
       end do
    else
       tau_iso_SGS (:,:,:) = 0.0_dp
    end if


    ! Reynolds tensor
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rmu          = - adi % sqgmr * mu (i,j,k)
                rmu_SGS      = - adi % sqgmr * mu_SGS (i,j,k) / rmu
                rmk          = - adi % sqgmr * kpa (i,j,k)

                bulk = rmk * ( fd (i,j,k,1) + fd (i,j,k,4) )

                sx_x (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) - fd (i,j,k,4) )
                sx_x_SGS (i,j,k) = bulk + ( sx_x (i,j,k) - bulk ) * ( 1.0_dp + rmu_SGS ) + tau_iso_SGS (i,j,k)

                sx_y (i,j,k) = rmu * ( fd (i,j,k,3) + fd (i,j,k,2) )
                sx_y_SGS (i,j,k) = sx_y (i,j,k) * ( 1.0_dp + rmu_SGS )

                sy_y (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,4) + fd (i,j,k,4) - fd (i,j,k,1) )
                sy_y_SGS (i,j,k) = bulk + ( sy_y (i,j,k) - bulk ) * ( 1.0_dp + rmu_SGS ) + tau_iso_SGS (i,j,k)

             end do
          end do
       end do

    end do

    deallocate ( mu , kpa )
    deallocate ( sx_x , sx_y , &
                        sy_y )

    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)

                uxs   = v (i,j,k,2) * rho_i
                vys   = v (i,j,k,3) * rho_i

                Ps    = 1.0_dp / P (i,j,k)
                Ts    = 1.0_dp / T (i,j,k)

                wrk1 = fd (i,j,k,5) * Ts
                wrk2 = fd (i,j,k,6) * Ts

                wrk3 = fd (i,j,k,7) * Ps
                wrk4 = fd (i,j,k,8) * Ps

                dm_SGS = adi % Sc * mu_SGS (i,j,k) * rho_i * Sc_SGS_i
                ct_SGS = adi % Pr * mu_SGS (i,j,k) * cp (i,j,k) * Pr_SGS_i

                ! Diffusive mass flux
                ! part : d_{\beta,j} = dXb/dxj + (Xb-Yb)/p * dp/dxj
                do l = 1 , nrv ! \beta loop
                   wrk = Xa (i,j,k,l) - v (i,j,k,niv+l) * rho_i

                   d_x (l) = fd (i,j,k,nfluid+ndim*l-1) + wrk * wrk3
                   d_y (l) = fd (i,j,k,nfluid+ndim*l  ) + wrk * wrk4
                end do

                ! Heat flux
                q_x (i,j,k) = 0.0_dp
                q_y (i,j,k) = 0.0_dp

                do l = 1 , nrv ! \alpha loop

                   ! Diffusive mass flux
                   ! rho * Ya * Vaj = - \sum_{\beta=1,N} rho * D_ab * ( d_bj + Xb * \chi_b / T * dT/dxj )
                   wrk3 = 0.0_dp ! rYV_x
                   wrk4 = 0.0_dp ! rYV_y

                   do m = 1 , nrv ! \beta loop

                      wrk  = Xa (i,j,k,m) * tdr (i,j,k,m)
                      wrk3 = wrk3 + rd (i,j,k,l,m) * ( d_x (m) + wrk * wrk1 )
                      wrk4 = wrk4 + rd (i,j,k,l,m) * ( d_y (m) + wrk * wrk2 )

                   end do

                   ! + SGS diffusive mass flux
                   ! rho * D_sgs * dYa/dxj = rho * D_SGS * Wa * ( W_i * dXa/dxj + Xa * dW_i/dxj )
                   wrk1 = v(i,j,k,1) * dm_SGS * thd % Wc (l)

                   fd (i,j,k,nfluid+ndim*l-1) = - rvd0 * ( wrk3 + wrk1 * ( W_i (i,j,k) * fd (i,j,k,nfluid+ndim*l-1) + &
                                                                           Xa (i,j,k,l) * fd (i,j,k,nSGS+1)          ))

                   fd (i,j,k,nfluid+ndim*l  ) = - rvd0 * ( wrk4 + wrk1 * ( W_i (i,j,k) * fd (i,j,k,nfluid+ndim*l  ) + &
                                                                           Xa (i,j,k,l) * fd (i,j,k,nSGS+2)          ))

                   ! Heat flux part 1
                   ! \sum_{\alpha=1,N} rho * Ya * Vaj * ( ha + R*T*\chi_a / Wa )
                   wrk1 = ha (i,j,k,l) + T (i,j,k) * thd % Wc_i (l) * tdr (i,j,k,l)
                   q_x (i,j,k) = q_x (i,j,k) + wrk3 * wrk1
                   q_y (i,j,k) = q_y (i,j,k) + wrk4 * wrk1

                end do

                ! Heat flux part 2 (with SGS terms)
                ! qj = - lambda * dT/dxj                                             &
                !      + \sum_{\alpha=1,N} rho * Ya * Vaj * ( ha + R*T*\chi_a / Wa ) &
                !      + tau_ij * u_i                                                &
                !      - lambda_SGS * dT/dxj                                         & (SGS term)
                !      - tau_ij_SGS * u_i                                              (SGS term)
                rlam = - rlam0 * ( ct (i,j,k) + ct_SGS )

                q_x (i,j,k) = q_x (i,j,k) +                  &
                              rlam * fd (i,j,k,10) +         &
                              uxs * sx_x_SGS (i,j,k) +       &
                              vys * sx_y_SGS (i,j,k)

                q_y (i,j,k) = q_y (i,j,k) +                  &
                              rlam * fd (i,j,k,11) +         &
                              uxs * sx_y_SGS (i,j,k) +       &
                              vys * sy_y_SGS (i,j,k)

             end do
          end do
       end do

    end do


    deallocate ( tdr )
    deallocate ( tau_iso_SGS )
    deallocate ( Xa , rd )


    ! passive scalar equations: viscous term
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   dm_SGS = adi % Sc * mu_SGS (i,j,k) / v (i,j,k,1) * Sc_SGS_i
                   P (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) & ! > laminar/molecular term ( P = -rlamcp )
                               - rvd0 * v (i,j,k,1) * dm_sgs         ! > SGS term for LES
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l-1) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l-1)
                      fd (i,j,k,nfluid+ndim*l  ) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l  )
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( ct , P )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               wrk_y (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xy 5')


    ! viscous flux in rho*u equation -> fl (:,:,:,2)
    call dx ( grid % dx_i , sx_x_SGS , wrk_x )
    call dy ( grid % dy_i , sx_y_SGS , wrk_y )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,2) = fl (i,j,k,2) +                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) )
          end do
       end do
    end do

    deallocate (sx_x_SGS)


    ! viscous flux in rho*v equation -> fl (:,:,:,3)
    call dx ( grid % dx_i , sx_y_SGS , wrk_x )
    call dy ( grid % dy_i , sy_y_SGS , wrk_y )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,3) = fl (i,j,k,3) +                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) )
          end do
       end do
    end do

    deallocate ( sx_y_SGS , sy_y_SGS )


    ! viscous flux in rho*e_tot equation -> fl (:,:,:,5)
    call dx ( grid % dx_i , q_x , wrk_x )
    call dy ( grid % dy_i , q_y , wrk_y )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,5) = fl (i,j,k,5) +                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) )
          end do
       end do
    end do

    deallocate ( q_x , q_y  )
    deallocate ( wrk_x , wrk_y )


    allocate ( wrka_x (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               wrka_y (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xy 6')


    ! viscous flux in rho*Y_alpha equation -> fl (:,:,:,niv+1:niv+npv+nvv)
    do l = 1 , nrv+npv
       call dx_fixed2 ( grid % dx_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                                 )
       call dy_fixed2 ( grid % dy_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) , &
                        wrka_y (sx:ex,sy:ey,sz:ez,l)                                 )
    end do

    do l = 1 , nrv+npv
       do k = sz , ez
          do j = sy , ey
             do i = sx , ex
                fl (i,j,k,niv+l) = fl (i,j,k,niv+l) +                    &
                                 ( wrka_x (i,j,k,l) + wrka_y (i,j,k,l) )
             end do
          end do
       end do
    end do

    deallocate (fd)
    deallocate ( wrka_x , wrka_y )


  end subroutine viscflux_eg_LES_xy


!> \brief 3D viscous flux.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_eg_xyz ( thd , adi , grid , T , W_i , cp , ha , v , &
                               rweg , iweg , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:) , intent (inout)       :: rweg !< real EGLIB work array
    integer (ip) , allocatable , dimension (:) , intent (inout)    :: iweg !< integer EGLIB work array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< viscous flux


    real (dp) , parameter                           :: twothirds = 2.0_dp / 3.0_dp
    integer (ip)                                    :: nfluid
    integer (ip)                                    :: ok , domain_id , i , j , k , l , m
    integer (ip)                                    :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                       :: rmu , rmk , rvd0 , rlam0 , rlamcp0 , bulk
    real (dp)                                       :: rho_i , Ps , Ts , uxs , vys , wzs
    real (dp) , dimension (nrv)                     :: d_x , d_y , d_z
    real (dp)                                       :: wrk1 , wrk2 , wrk3 , wrk4 , wrk5 , wrk6
    ! viscous variables (eglib)
    real (dp) , dimension (:,:,:) , allocatable     :: mu , kpa , ct
    real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
    real (dp) , dimension (:,:,:,:,:) , allocatable :: rd
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: ux , vy , wz , P
    real (dp) , dimension (:,:,:,:) , allocatable   :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: q_x  , q_y  , q_z  , &
                                                       sx_x , sx_y , sx_z , &
                                                       sy_y , sy_z , sz_z
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: wrk_x , wrk_y , wrk_z
    real (dp) , dimension (:,:,:,:) , allocatable   :: wrka_x , wrka_y , wrka_z
    ! array containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable   :: fd


    nfluid  = ndim * ( ndim + 2 )
    rvd0    = adi % sqgmr / adi % Sc
    rlam0   = adi % sqgmr * adi % ggmopr
    rlamcp0 = adi % sqgmr / adi % Pr


    allocate ( rd  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv,nrv) , &
               tdr (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)     , &
               mu  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               kpa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               ct  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xyz')


    call egvars ( thd , T , W_i , v , iweg , rweg , rd , mu , kpa , ct , tdr )


    ! Allocation MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xyz 2')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
                wz (i,j,k) = v (i,j,k,4) * rho_i
                P (i,j,k)  = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             end do
          end do
       end do

    end do


    ! passive scalar equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do l = nrv+1,nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      Xa (i,j,k,l) = v (i,j,k,niv+l) / v (i,j,k,1)
                   end do
                end do
             end do
          end do

       end do

    end if


    call dx_fixed1 ( grid % dx_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,1) )
    call dy_fixed1 ( grid % dy_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) )
    call dz_fixed1 ( grid % dz_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,3) )

    call dx_fixed1 ( grid % dx_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,4) )
    call dy_fixed1 ( grid % dy_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,5) )
    call dz_fixed1 ( grid % dz_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,6) )

    call dx_fixed1 ( grid % dx_i , wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,7) )
    call dy_fixed1 ( grid % dy_i , wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,8) )
    call dz_fixed1 ( grid % dz_i , wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                   fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,9) )

    call dx_fixed1 ( grid % dx_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,10) )
    call dy_fixed1 ( grid % dy_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,11) )
    call dz_fixed1 ( grid % dz_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,12) )

    call dx_fixed1 ( grid % dx_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,13) )
    call dy_fixed1 ( grid % dy_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,14) )
    call dz_fixed1 ( grid % dz_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,15) ) ! nfluid = 15 in 3D

    do l = 1 , nrv+npv
       call dx_fixed1 ( grid % dx_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-2) )
       call dy_fixed1 ( grid % dy_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) )
       call dz_fixed1 ( grid % dz_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) )
    end do

    deallocate ( ux , vy , wz )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species
    ! and/or non-pressure, you must specity the temperature and/or the
    ! velocity and/or the species and/or the pressure normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       fd (sx,sy:ey,sz:ez,1) = 0.0_dp ; fd (sx,sy:ey,sz:ez,4) = 0.0_dp ; fd (sx,sy:ey,sz:ez,7) = 0.0_dp ! normal velocity gradient
       fd (sx,sy:ey,sz:ez,10) = 0.0_dp ! normal temperature gradient
       fd (sx,sy:ey,sz:ez,13) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l-2) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       fd (ex,sy:ey,sz:ez,1) = 0.0_dp ; fd (ex,sy:ey,sz:ez,4) = 0.0_dp ; fd (ex,sy:ey,sz:ez,7) = 0.0_dp ! normal velocity gradient
       fd (ex,sy:ey,sz:ez,10) = 0.0_dp ! normal temperature gradient
       fd (ex,sy:ey,sz:ez,13) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l-2) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (S) == symmetryplane .or. bc (S) == adiabaticwall ) then
       fd (sx:ex,sy,sz:ez,2) = 0.0_dp ; fd (sx:ex,sy,sz:ez,5) = 0.0_dp ; fd (sx:ex,sy,sz:ez,8) = 0.0_dp ! normal velocity gradient
       fd (sx:ex,sy,sz:ez,11) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,sy,sz:ez,14) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (N) == symmetryplane .or. bc (N) == adiabaticwall ) then
       fd (sx:ex,ey,sz:ez,2) = 0.0_dp ; fd (sx:ex,ey,sz:ez,5) = 0.0_dp ; fd (sx:ex,ey,sz:ez,8) = 0.0_dp ! normal velocity gradient
       fd (sx:ex,ey,sz:ez,11) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,ey,sz:ez,14) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (B) == symmetryplane .or. bc (B) == adiabaticwall ) then
       fd (sx:ex,sy:ey,sz,3) = 0.0_dp ; fd (sx:ex,sy:ey,sz,6) = 0.0_dp ; fd (sx:ex,sy:ey,sz,9) = 0.0_dp ! normal velocity gradient
       fd (sx:ex,sy:ey,sz,12) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,sy:ey,sz,15) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy:ey,sz,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (F) == symmetryplane .or. bc (F) == adiabaticwall ) then
       fd (sx:ex,sy:ey,ez,3) = 0.0_dp ; fd (sx:ex,sy:ey,ez,6) = 0.0_dp ; fd (sx:ex,sy:ey,ez,9) = 0.0_dp ! normal velocity gradient
       fd (sx:ex,sy:ey,ez,12) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,sy:ey,ez,15) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy:ey,ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
    end if


    ! communicate the first derivative
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( sx_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               sx_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               sx_z (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               sy_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               sy_z (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               sz_z (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xyz 3')


    ! Reynolds tensor
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rmu = - adi % sqgmr * mu (i,j,k)
                rmk = - adi % sqgmr * kpa (i,j,k)

                bulk = rmk * ( fd (i,j,k,1) + fd (i,j,k,5) + fd (i,j,k,9) )

                sx_x (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) - &
                                                        ( fd (i,j,k,5) + fd (i,j,k,9) ) )

                sx_y (i,j,k) = rmu * ( fd (i,j,k,4) + fd (i,j,k,2) )

                sx_z (i,j,k) = rmu * ( fd (i,j,k,7) + fd (i,j,k,3) )

                sy_y (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,5) + fd (i,j,k,5) - &
                                                        ( fd (i,j,k,1) + fd (i,j,k,9) ) )

                sy_z (i,j,k) = rmu * ( fd (i,j,k,8) + fd (i,j,k,6) )

                sz_z (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,9) + fd (i,j,k,9) - &
                                                        ( fd (i,j,k,1) + fd (i,j,k,5) ) )

             end do
          end do
       end do

    end do

    deallocate ( mu , kpa )


    ! rYV_x , rYV_y , rYV_z
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                Ps    = 1.0_dp / P (i,j,k)
                Ts    = 1.0_dp / T (i,j,k)

                wrk1 = fd (i,j,k,10) * Ts
                wrk2 = fd (i,j,k,11) * Ts
                wrk3 = fd (i,j,k,12) * Ts

                wrk4 = fd (i,j,k,13) * Ps
                wrk5 = fd (i,j,k,14) * Ps
                wrk6 = fd (i,j,k,15) * Ps

                do l = 1 , nrv
                   ! d_x
                   d_x (l) = fd (i,j,k,nfluid+ndim*l-2) +               &
                           ( Xa (i,j,k,l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk4
                   ! d_y
                   d_y (l) = fd (i,j,k,nfluid+ndim*l-1) +               &
                           ( Xa (i,j,k,l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk5
                   ! d_z
                   d_z (l) = fd (i,j,k,nfluid+ndim*l  ) +               &
                           ( Xa (i,j,k,l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk6
                end do

                do l = 1 , nrv

                   fd (i,j,k,nfluid+ndim*l-2) = 0.0_dp ! rYV_x
                   fd (i,j,k,nfluid+ndim*l-1) = 0.0_dp ! rYV_y
                   fd (i,j,k,nfluid+ndim*l  ) = 0.0_dp ! rYV_z

                   do m = 1 , nrv

                      fd (i,j,k,nfluid+ndim*l-2) = fd (i,j,k,nfluid+ndim*l-2) +                    &
                                                   rd (i,j,k,l,m) *                                &
                                                 ( d_x (m) + Xa (i,j,k,m) * tdr (i,j,k,m) * wrk1 )
                      fd (i,j,k,nfluid+ndim*l-1) = fd (i,j,k,nfluid+ndim*l-1) +                    &
                                                   rd (i,j,k,l,m) *                                &
                                                 ( d_y (m) + Xa (i,j,k,m) * tdr (i,j,k,m) * wrk2 )
                      fd (i,j,k,nfluid+ndim*l  ) = fd (i,j,k,nfluid+ndim*l  ) +                    &
                                                   rd (i,j,k,l,m) *                                &
                                                 ( d_z (m) + Xa (i,j,k,m) * tdr (i,j,k,m) * wrk3 )

                   end do

                   fd (i,j,k,nfluid+ndim*l-2) = - rvd0 * fd (i,j,k,nfluid+ndim*l-2)
                   fd (i,j,k,nfluid+ndim*l-1) = - rvd0 * fd (i,j,k,nfluid+ndim*l-1)
                   fd (i,j,k,nfluid+ndim*l  ) = - rvd0 * fd (i,j,k,nfluid+ndim*l  )

                end do

             end do
          end do
       end do

    end do


    ! passive scalar equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   P (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) ! P = -rlamcp
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l-2) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l-2)
                      fd (i,j,k,nfluid+ndim*l-1) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l-1)
                      fd (i,j,k,nfluid+ndim*l  ) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l  )
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( P , Xa , rd )


    ! components of the viscous stresses
    allocate ( q_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               q_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               q_z (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xyz 4')


    ! q_x , q_y , q_z
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                uxs   = v (i,j,k,2) * rho_i
                vys   = v (i,j,k,3) * rho_i
                wzs   = v (i,j,k,4) * rho_i

                wrk1 = -rlam0 * ct (i,j,k)

                q_x (i,j,k) = 0.0_dp
                q_y (i,j,k) = 0.0_dp
                q_z (i,j,k) = 0.0_dp

                do l = 1 , nrv
                   q_x (i,j,k) = q_x (i,j,k) + fd (i,j,k,nfluid+ndim*l-2) *                  &
                               ( ha (i,j,k,l) + T (i,j,k) * thd % Wc_i (l) * tdr (i,j,k,l) )
                   q_y (i,j,k) = q_y (i,j,k) + fd (i,j,k,nfluid+ndim*l-1) *                  &
                               ( ha (i,j,k,l) + T (i,j,k) * thd % Wc_i (l) * tdr (i,j,k,l) )
                   q_z (i,j,k) = q_z (i,j,k) + fd (i,j,k,nfluid+ndim*l  ) *                  &
                               ( ha (i,j,k,l) + T (i,j,k) * thd % Wc_i (l) * tdr (i,j,k,l) )
                end do

                q_x (i,j,k) = q_x (i,j,k) + wrk1 * fd (i,j,k,10) +                         &
                              uxs * sx_x (i,j,k) + vys * sx_y (i,j,k) + wzs * sx_z (i,j,k)
                q_y (i,j,k) = q_y (i,j,k) + wrk1 * fd (i,j,k,11) +                         &
                              uxs * sx_y (i,j,k) + vys * sy_y (i,j,k) + wzs * sy_z (i,j,k)
                q_z (i,j,k) = q_z (i,j,k) + wrk1 * fd (i,j,k,12) +                         &
                              uxs * sx_z (i,j,k) + vys * sy_z (i,j,k) + wzs * sz_z (i,j,k)

             end do
          end do
       end do

    end do

    deallocate ( ct , tdr )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               wrk_y (sx:ex,sy:ey,sz:ez) , &
               wrk_z (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xyz 5')


    call dx ( grid % dx_i , sx_x , wrk_x )
    call dy ( grid % dy_i , sx_y , wrk_y )
    call dz ( grid % dz_i , sx_z , wrk_z )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,2) = fl (i,j,k,2) +                                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) + wrk_z (i,j,k) )
          end do
       end do
    end do

    deallocate (sx_x)


    call dx ( grid % dx_i , sx_y , wrk_x )
    call dy ( grid % dy_i , sy_y , wrk_y )
    call dz ( grid % dz_i , sy_z , wrk_z )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,3) = fl (i,j,k,3) +                                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) + wrk_z (i,j,k) )
          end do
       end do
    end do

    deallocate ( sx_y , sy_y )


    call dx ( grid % dx_i , sx_z , wrk_x )
    call dy ( grid % dy_i , sy_z , wrk_y )
    call dz ( grid % dz_i , sz_z , wrk_z )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,4) = fl (i,j,k,4) +                                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) + wrk_z (i,j,k) )
          end do
       end do
    end do

    deallocate ( sx_z , sy_z , sz_z )


    call dx ( grid % dx_i , q_x , wrk_x )
    call dy ( grid % dy_i , q_y , wrk_y )
    call dz ( grid % dz_i , q_z , wrk_z )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,5) = fl (i,j,k,5) +                                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) + wrk_z (i,j,k) )
          end do
       end do
    end do

    deallocate ( q_x , q_y , q_z )
    deallocate ( wrk_x , wrk_y , wrk_z )


    allocate ( wrka_x (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               wrka_y (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               wrka_z (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_xyz 6')

    do l = 1 , nrv+npv
       call dx_fixed2 ( grid % dx_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-2) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                                 )
       call dy_fixed2 ( grid % dy_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) , &
                        wrka_y (sx:ex,sy:ey,sz:ez,l)                                 )
       call dz_fixed2 ( grid % dz_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) , &
                        wrka_z (sx:ex,sy:ey,sz:ez,l)                                 )
    end do

    do l = 1 , nrv+npv
       do k = sz , ez
          do j = sy , ey
             do i = sx , ex
                fl (i,j,k,niv+l) = fl (i,j,k,niv+l) +                                       &
                                 ( wrka_x (i,j,k,l) + wrka_y (i,j,k,l) + wrka_z (i,j,k,l) )
             end do
          end do
       end do
    end do

    deallocate (fd)
    deallocate ( wrka_x , wrka_y , wrka_z )


  end subroutine viscflux_eg_xyz


!> \brief 3D viscous flux with LES models.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_eg_LES_xyz ( inp , thd , adi , grid , T , W_i , cp , ha , v , &
                                   rweg , iweg , mu_SGS , fl )


    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid   !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T      !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp     !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha     !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v      !< conserved variables array
    real (dp) , allocatable , dimension (:) , intent (inout)       :: rweg   !< real EGLIB work array
    integer (ip) , allocatable , dimension (:) , intent (inout)    :: iweg   !< integer EGLIB work array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mu_SGS !< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl     !< viscous flux


    real (dp) , parameter                           :: twothirds = 2.0_dp / 3.0_dp
    integer (ip)                                    :: nfluid , nSGS
    integer (ip)                                    :: ok , domain_id , i , j , k , l , m
    integer (ip)                                    :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                       :: rmu , rmu_SGS , rmk , rvd0 , rlam0 , rlam , rlamcp0 , bulk , &
                                                       Sc_SGS_i , Pr_SGS_i , gm1_g , ct_SGS , dm_SGS
    real (dp)                                       :: rho_i , Ps , Ts , uxs , vys , wzs
    real (dp) , dimension (nrv)                     :: d_x , d_y , d_z
    real (dp)                                       :: wrk , wrk1 , wrk2 , wrk3 , wrk4 , wrk5 , wrk6
    ! viscous variables (eglib)
    real (dp) , dimension (:,:,:) , allocatable     :: mu , kpa , ct
    real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
    real (dp) , dimension (:,:,:,:,:) , allocatable :: rd
    ! SGS additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: tau_iso_SGS
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: ux , vy , wz , P
    real (dp) , dimension (:,:,:,:) , allocatable   :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: q_x  , q_y  , q_z  , &
                                                       sx_x , sx_y , sx_z , &
                                                       sy_y , sy_z , sz_z
    real (dp) , dimension (:,:,:) , allocatable     :: sx_x_SGS , sx_y_SGS , sx_z_SGS , &
                                                       sy_y_SGS , sy_z_SGS , sz_z_SGS
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable     :: wrk_x , wrk_y , wrk_z
    real (dp) , dimension (:,:,:,:) , allocatable   :: wrka_x , wrka_y , wrka_z
    ! array containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable   :: fd


    nfluid  = ndim * ( ndim + 2 )
    nSGS    = nfluid + ndim * ( nrv + npv )

    rvd0    = adi % sqgmr / adi % Sc     ! M * sqrt(gamma) / ( Re * Sc )
    rlam0   = adi % sqgmr * adi % ggmopr ! M * sqrt(gamma) * gamma / ( Re * (gamma-1) * Pr )
    rlamcp0 = adi % sqgmr / adi % Pr     ! M * sqrt(gamma) / ( Re * Pr )

    Sc_SGS_i = 1.0_dp / inp % Sc_SGS
    Pr_SGS_i = 1.0_dp / inp % Pr_SGS
    gm1_g    = adi % gm1 / adi % gamma


    allocate ( rd          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv,nrv) , &
               tdr         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)     , &
               mu          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               kpa         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               ct          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xyz')


    call egvars ( thd , T , W_i , v , iweg , rweg , rd , mu , kpa , ct , tdr )


    ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xyz 2')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
                wz (i,j,k) = v (i,j,k,4) * rho_i
                P (i,j,k)  = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             end do
          end do
       end do

    end do


    ! passive scalar equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do l = nrv+1,nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      Xa (i,j,k,l) = v (i,j,k,niv+l) / v (i,j,k,1)
                   end do
                end do
             end do
          end do

       end do

    end if


    call dx_fixed1 ( grid % dx_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,1) ) ! du/dx = fd (:,:,:,1)
    call dy_fixed1 ( grid % dy_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) ) ! du/dy = fd (:,:,:,2)
    call dz_fixed1 ( grid % dz_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,3) ) ! du/dz = fd (:,:,:,3)

    call dx_fixed1 ( grid % dx_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,4) ) ! dv/dx = fd (:,:,:,4)
    call dy_fixed1 ( grid % dy_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,5) ) ! dv/dy = fd (:,:,:,5)
    call dz_fixed1 ( grid % dz_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,6) ) ! dv/dz = fd (:,:,:,6)

    call dx_fixed1 ( grid % dx_i , wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,7) ) ! dw/dx = fd (:,:,:,7)
    call dy_fixed1 ( grid % dy_i , wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,8) ) ! dw/dy = fd (:,:,:,8)
    call dz_fixed1 ( grid % dz_i , wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,9) ) ! dw/dz = fd (:,:,:,9)

    call dx_fixed1 ( grid % dx_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,10) ) ! dT/dx = fd (:,:,:,10)
    call dy_fixed1 ( grid % dy_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,11) ) ! dT/dy = fd (:,:,:,11)
    call dz_fixed1 ( grid % dz_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,12) ) ! dT/dz = fd (:,:,:,12)

    call dx_fixed1 ( grid % dx_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,13) ) ! dp/dx = fd (:,:,:,13)
    call dy_fixed1 ( grid % dy_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,14) ) ! dp/dy = fd (:,:,:,14)
    call dz_fixed1 ( grid % dz_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,15) ) ! dp/dz = fd (:,:,:,15) nfluid = 15 in 3D

    do l = 1 , nrv+npv
       call dx_fixed1 ( grid % dx_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-2) ) ! dXa/dx
       call dy_fixed1 ( grid % dy_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) ) ! dXa/dy
       call dz_fixed1 ( grid % dz_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) ) ! dXa/dz
    end do


    ! SGS: W_i derivatives need to be communicated
    call dx_fixed1 ( grid % dx_i , W_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nSGS+1) )    ! dW_i/dx
    call dy_fixed1 ( grid % dy_i , W_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nSGS+2) )    ! dW_i/dy
    call dz_fixed1 ( grid % dz_i , W_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nSGS+3) )    ! dW_i/dz


    ! SGS: shock detector criteria needs to be communicated
    call shock_det_slv ( v , T , W_i , fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) )


    deallocate ( ux , vy , wz )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species
    ! and/or non-pressure, you must specity the temperature and/or the
    ! velocity and/or the species and/or the pressure normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       if ( bc (W) == symmetryplane ) then
          fd (sx,sy:ey,sz:ez,1) = 0.0_dp
          fd (sx,sy:ey,sz:ez,4) = 0.0_dp
          fd (sx,sy:ey,sz:ez,7) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx,sy:ey,sz:ez,10) = 0.0_dp ! normal temperature gradient
       fd (sx,sy:ey,sz:ez,13) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l-2) = 0.0_dp ! normal species gradient
       end do
       fd (sx,sy:ey,sz:ez,nSGS+1) = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       if ( bc (E) == symmetryplane ) then
          fd (ex,sy:ey,sz:ez,1) = 0.0_dp
          fd (ex,sy:ey,sz:ez,4) = 0.0_dp
          fd (ex,sy:ey,sz:ez,7) = 0.0_dp ! normal velocity gradient
       end if
       fd (ex,sy:ey,sz:ez,10) = 0.0_dp ! normal temperature gradient
       fd (ex,sy:ey,sz:ez,13) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l-2) = 0.0_dp ! normal species gradient
       end do
       fd (ex,sy:ey,sz:ez,nSGS+1) = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (S) == symmetryplane .or. bc (S) == adiabaticwall ) then
       if ( bc (S) == symmetryplane ) then
          fd (sx:ex,sy,sz:ez,2) = 0.0_dp
          fd (sx:ex,sy,sz:ez,5) = 0.0_dp
          fd (sx:ex,sy,sz:ez,8) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,sy,sz:ez,11) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,sy,sz:ez,14) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
       fd (sx:ex,sy,sz:ez,nSGS+2) = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (N) == symmetryplane .or. bc (N) == adiabaticwall ) then
       if ( bc (N) == symmetryplane ) then
          fd (sx:ex,ey,sz:ez,2) = 0.0_dp
          fd (sx:ex,ey,sz:ez,5) = 0.0_dp
          fd (sx:ex,ey,sz:ez,8) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,ey,sz:ez,11) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,ey,sz:ez,14) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
       fd (sx:ex,ey,sz:ez,nSGS+2) = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (B) == symmetryplane .or. bc (B) == adiabaticwall ) then
       if ( bc (B) == symmetryplane ) then
          fd (sx:ex,sy:ey,sz,3) = 0.0_dp
          fd (sx:ex,sy:ey,sz,6) = 0.0_dp
          fd (sx:ex,sy:ey,sz,9) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,sy:ey,sz,12) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,sy:ey,sz,15) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy:ey,sz,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
       fd (sx:ex,sy:ey,sz,nSGS+3) = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (F) == symmetryplane .or. bc (F) == adiabaticwall ) then
       if ( bc (F) == symmetryplane ) then
          fd (sx:ex,sy:ey,ez,3) = 0.0_dp
          fd (sx:ex,sy:ey,ez,6) = 0.0_dp
          fd (sx:ex,sy:ey,ez,9) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,sy:ey,ez,12) = 0.0_dp ! normal temperature gradient
       fd (sx:ex,sy:ey,ez,15) = 0.0_dp ! normal pressure gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy:ey,ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
       fd (sx:ex,sy:ey,ez,nSGS+3) = 0.0_dp ! normal mixture molecular weight gradient
    end if


    ! communicate the first derivative
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( sx_x         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_y         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_z         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sy_y         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sy_z         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sz_z         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_x_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_y_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_z_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sy_y_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sy_z_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sz_z_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               tau_iso_SGS  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xyz 3')

    ! components of the viscous stresses
    allocate ( q_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               q_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               q_z (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xyz 4')


    ! SGS: viscosity
    do domain_id = -ndim , ndim
       call mu_SGS_selector ( inp , adi , domain_id , grid , v , fd , mu_SGS )
    end do

    ! SGS: tau_iso_SGS derivatives need to be COMMUNICATED
    if ( inp % tau_iso_SGS_switch ) then
       do domain_id = -ndim , ndim
          call tau_iso_SGS_selector ( inp , domain_id , grid , v , fd , tau_iso_SGS )
       end do
    else
       tau_iso_SGS (:,:,:) = 0.0_dp
    end if


    ! Reynolds tensor
    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rmu          = - adi % sqgmr * mu (i,j,k)
                rmu_SGS      = - adi % sqgmr * mu_SGS (i,j,k) / rmu
                rmk          = - adi % sqgmr * kpa (i,j,k)

                bulk = rmk * ( fd (i,j,k,1) + fd (i,j,k,5) + fd (i,j,k,9) )

                sx_x (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) - &
                                                        ( fd (i,j,k,5) + fd (i,j,k,9) ) )
                sx_x_SGS (i,j,k) = bulk + ( sx_x (i,j,k) - bulk ) * ( 1.0_dp + rmu_SGS ) + tau_iso_SGS (i,j,k)

                sx_y (i,j,k) = rmu * ( fd (i,j,k,4) + fd (i,j,k,2) )
                sx_y_SGS (i,j,k) = sx_y (i,j,k) * ( 1.0_dp + rmu_SGS )

                sx_z (i,j,k) = rmu * ( fd (i,j,k,7) + fd (i,j,k,3) )
                sx_z_SGS (i,j,k) = sx_z (i,j,k) * ( 1.0_dp + rmu_SGS )

                sy_y (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,5) + fd (i,j,k,5) - &
                                                        ( fd (i,j,k,1) + fd (i,j,k,9) ) )
                sy_y_SGS (i,j,k) = bulk + ( sy_y (i,j,k) - bulk ) * ( 1.0_dp + rmu_SGS ) + tau_iso_SGS (i,j,k)

                sy_z (i,j,k) = rmu * ( fd (i,j,k,8) + fd (i,j,k,6) )
                sy_z_SGS (i,j,k) = sy_z (i,j,k) * ( 1.0_dp + rmu_SGS )

                sz_z (i,j,k) = bulk + twothirds * rmu * ( fd (i,j,k,9) + fd (i,j,k,9) - &
                                                        ( fd (i,j,k,1) + fd (i,j,k,5) ) )
                sz_z_SGS (i,j,k) = bulk + ( sz_z (i,j,k) - bulk ) * ( 1.0_dp + rmu_SGS ) + tau_iso_SGS (i,j,k)

             end do
          end do
       end do

    end do

    deallocate ( mu , kpa )
    deallocate ( sx_x , sx_y , sx_z , &
                        sy_y , sy_z , &
                               sz_z )

    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)

                uxs   = v (i,j,k,2) * rho_i
                vys   = v (i,j,k,3) * rho_i
                wzs   = v (i,j,k,4) * rho_i

                Ps    = 1.0_dp / P (i,j,k)
                Ts    = 1.0_dp / T (i,j,k)

                wrk1 = fd (i,j,k,10) * Ts
                wrk2 = fd (i,j,k,11) * Ts
                wrk3 = fd (i,j,k,12) * Ts

                wrk4 = fd (i,j,k,13) * Ps
                wrk5 = fd (i,j,k,14) * Ps
                wrk6 = fd (i,j,k,15) * Ps

                dm_SGS = adi % Sc * mu_SGS (i,j,k) * rho_i * Sc_SGS_i
                ct_SGS = adi % Pr * mu_SGS (i,j,k) * cp (i,j,k) * Pr_SGS_i

                ! Diffusive mass flux
                ! part : d_{\beta,j} = dXb/dxj + (Xb-Yb)/p * dp/dxj
                do l = 1 , nrv ! \beta loop
                   wrk = Xa (i,j,k,l) - v (i,j,k,niv+l) * rho_i

                   d_x (l) = fd (i,j,k,nfluid+ndim*l-2) + wrk * wrk4
                   d_y (l) = fd (i,j,k,nfluid+ndim*l-1) + wrk * wrk5
                   d_z (l) = fd (i,j,k,nfluid+ndim*l  ) + wrk * wrk6
                end do

                ! Heat flux
                q_x (i,j,k) = 0.0_dp
                q_y (i,j,k) = 0.0_dp
                q_z (i,j,k) = 0.0_dp

                do l = 1 , nrv ! \alpha loop

                   ! Diffusive mass flux
                   ! rho * Ya * Vaj = - \sum_{\beta=1,N} rho * D_ab * ( d_bj + Xb * \chi_b / T * dT/dxj )
                   wrk4 = 0.0_dp ! rYV_x
                   wrk5 = 0.0_dp ! rYV_y
                   wrk6 = 0.0_dp ! rYV_z

                   do m = 1 , nrv ! \beta loop

                      wrk  = Xa (i,j,k,m) * tdr (i,j,k,m)
                      wrk4 = wrk4 + rd (i,j,k,l,m) * ( d_x (m) + wrk * wrk1 )
                      wrk5 = wrk5 + rd (i,j,k,l,m) * ( d_y (m) + wrk * wrk2 )
                      wrk6 = wrk6 + rd (i,j,k,l,m) * ( d_z (m) + wrk * wrk3 )

                   end do

                   ! + SGS diffusive mass flux
                   ! rho * D_sgs * dYa/dxj = rho * D_SGS * Wa * ( W_i * dXa/dxj + Xa * dW_i/dxj )
                   wrk1 = v(i,j,k,1) * dm_SGS * thd % Wc (l)

                   fd (i,j,k,nfluid+ndim*l-2) = - rvd0 * ( wrk4 + wrk1 * ( W_i (i,j,k) * fd (i,j,k,nfluid+ndim*l-2) + &
                                                                           Xa (i,j,k,l) * fd (i,j,k,nSGS+1)          ))

                   fd (i,j,k,nfluid+ndim*l-1) = - rvd0 * ( wrk5 + wrk1 * ( W_i (i,j,k) * fd (i,j,k,nfluid+ndim*l-1) + &
                                                                           Xa (i,j,k,l) * fd (i,j,k,nSGS+2)          ))

                   fd (i,j,k,nfluid+ndim*l  ) = - rvd0 * ( wrk6 + wrk1 * ( W_i (i,j,k) * fd (i,j,k,nfluid+ndim*l  ) + &
                                                                           Xa (i,j,k,l) * fd (i,j,k,nSGS+3)          ))

                   ! Heat flux part 1
                   ! \sum_{\alpha=1,N} rho * Ya * Vaj * ( ha + R*T*\chi_a / Wa )
                   wrk1 = ha (i,j,k,l) + T (i,j,k) * thd % Wc_i (l) * tdr (i,j,k,l)
                   q_x (i,j,k) = q_x (i,j,k) + wrk4 * wrk1
                   q_y (i,j,k) = q_y (i,j,k) + wrk5 * wrk1
                   q_z (i,j,k) = q_z (i,j,k) + wrk6 * wrk1

                end do

                ! Heat flux part 2 (with SGS terms)
                ! qj = - lambda * dT/dxj                                             &
                !      + \sum_{\alpha=1,N} rho * Ya * Vaj * ( ha + R*T*\chi_a / Wa ) &
                !      + tau_ij * u_i                                                &
                !      - lambda_SGS * dT/dxj                                         & (SGS term)
                !      - tau_ij_SGS * u_i                                              (SGS term)
                rlam = - rlam0 * ( ct (i,j,k) + ct_SGS )

                q_x (i,j,k) = q_x (i,j,k) +                  &
                              rlam * fd (i,j,k,10) +         &
                              uxs * sx_x_SGS (i,j,k) +       &
                              vys * sx_y_SGS (i,j,k) +       &
                              wzs * sx_z_SGS (i,j,k)

                q_y (i,j,k) = q_y (i,j,k) +                  &
                              rlam * fd (i,j,k,11) +         &
                              uxs * sx_y_SGS (i,j,k) +       &
                              vys * sy_y_SGS (i,j,k) +       &
                              wzs * sy_z_SGS (i,j,k)

                q_z (i,j,k) = q_z (i,j,k) +                  &
                              rlam * fd (i,j,k,12) +         &
                              uxs * sx_z_SGS (i,j,k) +       &
                              vys * sy_z_SGS (i,j,k) +       &
                              wzs * sz_z_SGS (i,j,k)

             end do
          end do
       end do

    end do


    deallocate ( tdr )
    deallocate ( tau_iso_SGS )
    deallocate ( Xa , rd )


    ! passive scalar equations: viscous term
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   dm_SGS = adi % Sc * mu_SGS (i,j,k) / v (i,j,k,1) * Sc_SGS_i
                   P (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) & ! > laminar/molecular term ( P = -rlamcp )
                               - rvd0 * v (i,j,k,1) * dm_sgs         ! > SGS term for LES
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l-2) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l-2)
                      fd (i,j,k,nfluid+ndim*l-1) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l-1)
                      fd (i,j,k,nfluid+ndim*l  ) = P (i,j,k) * fd (i,j,k,nfluid+ndim*l  )
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( ct , P )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               wrk_y (sx:ex,sy:ey,sz:ez) , &
               wrk_z (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xyz 5')


    ! viscous flux in rho*u equation -> fl (:,:,:,2)
    call dx ( grid % dx_i , sx_x_SGS , wrk_x )
    call dy ( grid % dy_i , sx_y_SGS , wrk_y )
    call dz ( grid % dz_i , sx_z_SGS , wrk_z )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,2) = fl (i,j,k,2) +                                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) + wrk_z (i,j,k) )
          end do
       end do
    end do

    deallocate (sx_x_SGS)


    ! viscous flux in rho*v equation -> fl (:,:,:,3)
    call dx ( grid % dx_i , sx_y_SGS , wrk_x )
    call dy ( grid % dy_i , sy_y_SGS , wrk_y )
    call dz ( grid % dz_i , sy_z_SGS , wrk_z )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,3) = fl (i,j,k,3) +                                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) + wrk_z (i,j,k) )
          end do
       end do
    end do

    deallocate ( sx_y_SGS , sy_y_SGS )


    ! viscous flux in rho*w equation -> fl (:,:,:,4)
    call dx ( grid % dx_i , sx_z_SGS , wrk_x )
    call dy ( grid % dy_i , sy_z_SGS , wrk_y )
    call dz ( grid % dz_i , sz_z_SGS , wrk_z )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,4) = fl (i,j,k,4) +                                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) + wrk_z (i,j,k) )
          end do
       end do
    end do

    deallocate ( sx_z_SGS , sy_z_SGS , sz_z_SGS )


    ! viscous flux in rho*e_tot equation -> fl (:,:,:,5)
    call dx ( grid % dx_i , q_x , wrk_x )
    call dy ( grid % dy_i , q_y , wrk_y )
    call dz ( grid % dz_i , q_z , wrk_z )

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             fl (i,j,k,5) = fl (i,j,k,5) +                                  &
                          ( wrk_x (i,j,k) + wrk_y (i,j,k) + wrk_z (i,j,k) )
          end do
       end do
    end do

    deallocate ( q_x , q_y , q_z )
    deallocate ( wrk_x , wrk_y , wrk_z )


    allocate ( wrka_x (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               wrka_y (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               wrka_z (sx:ex,sy:ey,sz:ez,nrv+npv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_eg_LES_xyz 6')


    ! viscous flux in rho*Y_alpha equation -> fl (:,:,:,niv+1:niv+npv+nvv)
    do l = 1 , nrv+npv
       call dx_fixed2 ( grid % dx_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-2) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                                 )
       call dy_fixed2 ( grid % dy_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) , &
                        wrka_y (sx:ex,sy:ey,sz:ez,l)                                 )
       call dz_fixed2 ( grid % dz_i                                                  , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) , &
                        wrka_z (sx:ex,sy:ey,sz:ez,l)                                 )
    end do

    do l = 1 , nrv+npv
       do k = sz , ez
          do j = sy , ey
             do i = sx , ex
                fl (i,j,k,niv+l) = fl (i,j,k,niv+l) +                                       &
                                 ( wrka_x (i,j,k,l) + wrka_y (i,j,k,l) + wrka_z (i,j,k,l) )
             end do
          end do
       end do
    end do

    deallocate (fd)
    deallocate ( wrka_x , wrka_y , wrka_z )


  end subroutine viscflux_eg_LES_xyz


end module eg_lib
