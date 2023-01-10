!------------------------------------------------------------------------------
! MODULE: weno
!------------------------------------------------------------------------------
!> \brief Inviscid fluxes through WENO.
!!
!! This module allows to calculate the inviscid fluxes using a
!! 7th-order WENO scheme.
!!
!! Reference article (pure air): C.-W. Shu. _Essentially
!!   non-oscillatory and weighted essentially non-oscillatory schemes
!!   for hyperbolic conservation laws_. Technical report, ICASE
!!   (1997).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module weno

  use parameters
  use parallel
  use type_thd
  use adim
  use thermodynamics , only : Wmix_i_scalar , cp_scalar , ha_scalar

  implicit none

  !> WENO internal parameters
  real (dp) , parameter , private :: eps_weno = 1.0e-10_dp
  real (dp) , private :: cweno  ( -1:ng-1 , 0:ng-1 ) , &
                         dweno  ( 0:ng-1 ) ,           &
                         dtweno ( 0:ng-1 )
  real (dp) , private :: za , zb , zc , zd , ze , zf , zg , zh , &
                         zi , zl , zm , zn , zo , zp , zq , zr , &
                         zs , zt , zu
  real (dp) , private :: cnr ( 2*ng-1 )

contains


!> \brief Definition of WENO parameters.
!!
!! Definition of matrix c(r,j) of Shu (page 8) end of the coefficients
!! for the smoothness indicators
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine wenopar ( opt , optord )


    logical      , intent (in) :: opt    !< optimized WENO coefficients to reduce dispersion
    integer (ip) , intent (in) :: optord !< order of the optimization


    if ( ng == 4 ) then ! coefficients for 7-th order


       if ( .not. opt ) then ! non-optimized scheme


          cweno (-1,0) =   25.0_dp / 12.0_dp
          cweno (-1,1) = - 23.0_dp / 12.0_dp
          cweno (-1,2) =   13.0_dp / 12.0_dp
          cweno (-1,3) = -  1.0_dp /  4.0_dp

          cweno (0,0) =    1.0_dp /  4.0_dp
          cweno (0,1) =   13.0_dp / 12.0_dp
          cweno (0,2) = -  5.0_dp / 12.0_dp
          cweno (0,3) =    1.0_dp / 12.0_dp

          cweno (1,0) = -  1.0_dp / 12.0_dp
          cweno (1,1) =    7.0_dp / 12.0_dp
          cweno (1,2) =    7.0_dp / 12.0_dp
          cweno (1,3) = -  1.0_dp / 12.0_dp

          cweno (2,0) =    1.0_dp / 12.0_dp
          cweno (2,1) = -  5.0_dp / 12.0_dp
          cweno (2,2) =   13.0_dp / 12.0_dp
          cweno (2,3) =    1.0_dp /  4.0_dp

          cweno (3,0) = -  1.0_dp /  4.0_dp
          cweno (3,1) =   13.0_dp / 12.0_dp
          cweno (3,2) = - 23.0_dp / 12.0_dp
          cweno (3,3) =   25.0_dp / 12.0_dp

          ! Ideal weights d dtilde of Shu (page 17)
          ! (r=0,1,2,3)

          dweno (0)   =    4.0_dp / 35.0_dp
          dweno (1)   =   18.0_dp / 35.0_dp
          dweno (2)   =   12.0_dp / 35.0_dp
          dweno (3)   =    1.0_dp / 35.0_dp

          dtweno (0)  = dweno (3)
          dtweno (1)  = dweno (2)
          dtweno (2)  = dweno (1)
          dtweno (3)  = dweno (0)


       else ! optimized scheme


          if ( optord == 1 ) then


             cweno (-1,0) =   1.8839434_dp
             cweno (-1,1) = - 1.5416604_dp
             cweno (-1,2) =   0.89823611_dp
             cweno (-1,3) = - 0.24051909_dp

             cweno (0,0)  =   0.289506_dp
             cweno (0,1)  =   1.02041_dp
             cweno (0,2)  = - 0.404792_dp
             cweno (0,3)  =   0.0948743_dp

             cweno (1,0)  = - 0.100769_dp
             cweno (1,1)  =   0.600769_dp
             cweno (1,2)  =   0.600769_dp
             cweno (1,3)  = - 0.100769_dp

             cweno (2,0)  =   0.0948743_dp
             cweno (2,1)  = - 0.404792_dp
             cweno (2,2)  =   1.02041_dp
             cweno (2,3)  =   0.289506_dp

             cweno (3,0)  = - 0.240519_dp
             cweno (3,1)  =   0.898236_dp
             cweno (3,2)  = - 1.54166_dp
             cweno (3,3)  =   1.88394_dp

             ! Ideal weights d dtilde of Shu (page 17)
             ! (r=0,1,2,3)

             dweno (0)    =   0.150244_dp
             dweno (1)    =   0.480178_dp
             dweno (2)    =   0.32988_dp
             dweno (3)    =   0.0396983_dp

             dtweno (0)   = dweno (3)
             dtweno (1)   = dweno (2)
             dtweno (2)   = dweno (1)
             dtweno (3)   = dweno (0)


          else if ( optord == 3 ) then


             cweno (-1,0) =   1.9228346_dp
             cweno (-1,1) = - 1.6250794_dp
             cweno (-1,2) =   0.98165507_dp
             cweno (-1,3) = - 0.27941025_dp

             cweno (0,0)  =   0.28418590_dp
             cweno (0,1)  =   1.0318226_dp
             cweno (0,2)  = - 0.41620299_dp
             cweno (0,3)  =   0.10019444_dp

             cweno (1,0)  = - 0.10076912_dp
             cweno (1,1)  =   0.60076912_dp
             cweno (1,2)  =   0.60076912_dp
             cweno (1,3)  = - 0.10076912_dp

             cweno (2,0)  =   0.10019444_dp
             cweno (2,1)  = - 0.41620299_dp
             cweno (2,2)  =   1.0318226_dp
             cweno (2,3)  =   0.28418590_dp

             cweno (3,0) = - 0.27941025_dp
             cweno (3,1) =   0.98165507_dp
             cweno (3,2) = - 1.6250794_dp
             cweno (3,3) =   1.9228346_dp

             ! Ideal weights d dtilde of Shu (page 17)
             ! (r=0,1,2,3)

             dweno (0)   =   0.14150117_dp
             dweno (1)   =   0.48616615_dp
             dweno (2)   =   0.33383476_dp
             dweno (3)   =   0.038497919_dp

             dtweno (0)  =   dweno (3)
             dtweno (1)  =   dweno (2)
             dtweno (2)  =   dweno (1)
             dtweno (3)  =   dweno (0)


          else if ( optord == 5 ) then


             cweno (-1,0) =   2.0177891_dp
             cweno (-1,1) = - 1.7200339_dp
             cweno (-1,2) =   0.88670057_dp
             cweno (-1,3) = - 0.18445575_dp

             cweno (0,0)  =   0.25866239_dp
             cweno (0,1)  =   1.0573462_dp
             cweno (0,2)  = - 0.39067948_dp
             cweno (0,3)  =   0.074670939_dp

             cweno (1,0)  = - 0.083333333_dp
             cweno (1,1)  =   0.58333333_dp
             cweno (1,2)  =   0.58333333_dp
             cweno (1,3)  = - 0.083333333_dp

             cweno (2,0)  =   0.074670939_dp
             cweno (2,1)  = - 0.39067948_dp
             cweno (2,2)  =   1.0573462_dp
             cweno (2,3)  =   0.25866239_dp

             cweno (3,0)  = - 0.18445575_dp
             cweno (3,1)  =   0.88670057_dp
             cweno (3,2)  = - 1.7200339_dp
             cweno (3,3)  =   2.0177891_dp

             ! Ideal weights d dtilde of Shu (page 17)
             ! (r=0,1,2,3)

             dweno (0)    =   0.14196688_dp
             dweno (1)    =   0.51976365_dp
             dweno (2)    =   0.31535440_dp
             dweno (3)    =   0.022915068_dp

             dtweno (0)   =   dweno (3)
             dtweno (1)   =   dweno (2)
             dtweno (2)   =   dweno (1)
             dtweno (3)   =   dweno (0)


          else


             write (*,*) 'not implemented in wenopar'
             call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )


          end if


       end if


    else


       write (*,*) 'not implemented in wenopar'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )


    end if


    ! Constants for evaluation of the beta's
    za =     547.0_dp
    zb = -  1854.0_dp
    zc =    4642.0_dp
    zd = -  3882.0_dp
    ze = - 17246.0_dp
    zf =    2107.0_dp
    zg =    7043.0_dp
    zh = -  9402.0_dp
    zi =   11003.0_dp
    zl =    7042.0_dp
    zm =    3443.0_dp
    zn = -  5966.0_dp
    zo = -  2522.0_dp
    zp =    1602.0_dp
    zq =    2843.0_dp
    zr =     267.0_dp
    zs =    1922.0_dp
    zt = -   494.0_dp
    zu = -  1642.0_dp


    ! Definition of parameters for upwind derivatives at the boundaries
    cnr (1) = - 49.0_dp / 20.0_dp
    cnr (2) =    6.0_dp
    cnr (3) = - 15.0_dp / 2.0_dp
    cnr (4) =   20.0_dp / 3.0_dp
    cnr (5) = - 15.0_dp / 4.0_dp
    cnr (6) =    6.0_dp / 5.0_dp
    cnr (7) = -  1.0_dp / 6.0_dp


  end subroutine wenopar


!> \brief WENO reconstruction using non-linear weights.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine wenorec_nw ( vp , vm , vminus , vplus )


    real (dp) , dimension (:) , intent (in)          :: vp     !< plus value (in)
    real (dp) , dimension (:) , intent (in)          :: vm     !< minus value (in)
    real (dp) , intent (inout)                       :: vminus !< minus value (out)
    real (dp) , intent (inout)                       :: vplus  !< plus value (out)

    real (dp) , parameter :: denom = 1.0_dp / 420.0_dp


    vminus = ( - 3.0_dp * vp(1) + 25.0_dp * vp(2) - 101.0_dp * vp(3) + &
               319.0_dp * vp(4) + 214.0_dp * vp(5) - 38.0_dp * vp(6) + 4.0_dp * vp(7) ) * denom

    vplus  = ( - 3.0_dp * vm(8) + 25.0_dp * vm(7) - 101.0_dp * vm(6) + &
               319.0_dp * vm(5) + 214.0_dp * vm(4) - 38.0_dp * vm(3) + 4.0_dp * vm(2) ) * denom


  end subroutine wenorec_nw


!> \brief WENO reconstruction using the centered stencil (max accuracy and min dispersion).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine wenorec ( vp , vm , vminus , vplus )


    real (dp) , dimension (:) , intent (in)          :: vp     !< plus value (in)
    real (dp) , dimension (:) , intent (in)          :: vm     !< minus value (in)
    real (dp) , intent (inout)                       :: vminus !< minus value (out)
    real (dp) , intent (inout)                       :: vplus  !< plus value (out)

    integer (ip) :: i , j , r

    real (dp) :: sum

    real (dp) , dimension (0:ng-1) :: vpr , vmr , alfa , alfat , beta , omega , omegat


    !     Evaluation of the '-' state at the interface i+1/2


    i = 4 ! (i-3 = 1, i+3 = 7)


    do r = 0 , 3
       vpr(r) = 0.0_dp
       do j = 0 , 3
          vpr(r) = vpr(r) + cweno(r,j) * vp(i-r+j)
       end do
    end do


    !     Evaluation of beta[r] for r = 0,1,2,3


    beta (0) = vp(i)   * ( zf*vp(i) + zh*vp(i+1) + zl*vp(i+2) + zb*vp(i+3) ) + &
               vp(i+1) * (            zi*vp(i+1) + ze*vp(i+2) + zc*vp(i+3) ) + &
               vp(i+2) * (                         zg*vp(i+2) + zd*vp(i+3) ) + &
               vp(i+3) * (                                      za*vp(i+3) )

    beta (1) = vp(i-1) * ( za*vp(i-1) + zo*vp(i) + zs*vp(i+1) + zt*vp(i+2) ) + &
               vp(i)   * (              zm*vp(i) + zn*vp(i+1) + zp*vp(i+2) ) + &
               vp(i+1) * (                         zq*vp(i+1) + zu*vp(i+2) ) + &
               vp(i+2) * (                                      zr*vp(i+2) )

    beta (2) = vp(i-2) * ( zr*vp(i-2) + zu*vp(i-1) + zp*vp(i) + zt*vp(i+1) ) + &
               vp(i-1) * (              zq*vp(i-1) + zn*vp(i) + zs*vp(i+1) ) + &
               vp(i)   * (                           zm*vp(i) + zo*vp(i+1) ) + &
               vp(i+1) * (                                      za*vp(i+1) )

    beta (3) = vp(i-3) * ( za*vp(i-3) + zd*vp(i-2) + zc*vp(i-1) + zb*vp(i) ) + &
               vp(i-2) * (              zg*vp(i-2) + ze*vp(i-1) + zl*vp(i) ) + &
               vp(i-1) * (                           zi*vp(i-1) + zh*vp(i) ) + &
               vp(i)   * (                                        zf*vp(i) )


    !     Evaluation of the alfa's for r = 0,1,2,3


    sum = 0.0_dp
    do r = 0 , 3
       alfa(r) = dweno(r) / ( ( eps_weno + beta(r) ) * ( eps_weno + beta(r) ) )
       sum     = sum + alfa(r)
    end do


    !     Evaluation of the omega's for r = 0,1,2,3


    sum = 1.0_dp / sum
    do r = 0 , 3
       omega(r) = alfa(r) * sum
    end do


    !     Formula (2.52) of Shu, i.e. the 7-th order reconstruction
    !     Evaluate v(-)_i+1/2 = vminus


    vminus = 0.0_dp
    do r = 0 , 3
       vminus = vminus + omega(r) * vpr(r)
    end do


    !     Evaluation of the '+' state at the interface i+1/2


    i = 5 ! (i-3 = 2, i+3 = 8)


    do r = 0 , 3
       vmr(r) = 0.0_dp
       do j = 0 , 3
          vmr(r) = vmr(r) + cweno(r-1,j) * vm(i-r+j)
       end do
    end do


    !     Evaluation of beta[r] for r = 0,1,2,3


    beta (0) = vm(i)   * ( zf*vm(i) + zh*vm(i+1) + zl*vm(i+2) + zb*vm(i+3) ) + &
               vm(i+1) * (            zi*vm(i+1) + ze*vm(i+2) + zc*vm(i+3) ) + &
               vm(i+2) * (                         zg*vm(i+2) + zd*vm(i+3) ) + &
               vm(i+3) * (                                      za*vm(i+3) )

    beta (1) = vm(i-1) * ( za*vm(i-1) + zo*vm(i) + zs*vm(i+1) + zt*vm(i+2) ) + &
               vm(i)   * (              zm*vm(i) + zn*vm(i+1) + zp*vm(i+2) ) + &
               vm(i+1) * (                         zq*vm(i+1) + zu*vm(i+2) ) + &
               vm(i+2) * (                                      zr*vm(i+2) )

    beta (2) = vm(i-2) * ( zr*vm(i-2) + zu*vm(i-1) + zp*vm(i) + zt*vm(i+1) ) + &
               vm(i-1) * (              zq*vm(i-1) + zn*vm(i) + zs*vm(i+1) ) + &
               vm(i)   * (                           zm*vm(i) + zo*vm(i+1) ) + &
               vm(i+1) * (                                      za*vm(i+1) )

    beta (3) = vm(i-3) * ( za*vm(i-3) + zd*vm(i-2) + zc*vm(i-1) + zb*vm(i) ) + &
               vm(i-2) * (              zg*vm(i-2) + ze*vm(i-1) + zl*vm(i) ) + &
               vm(i-1) * (                           zi*vm(i-1) + zh*vm(i) ) + &
               vm(i)   * (                                        zf*vm(i) )


    !     Evaluation of the alfatilde's for r = 0,1,2,3


    sum = 0.0_dp
    do r = 0 , 3
       alfat(r) = dtweno(r) / ( ( eps_weno + beta(r) ) * ( eps_weno + beta(r) ) )
       sum      = sum + alfat(r)
    end do


    !     Evaluation of the omegatilde's for r = 0,1,2,3


    sum = 1.0_dp / sum
    do r = 0 , 3
       omegat(r) = alfat(r) * sum
    end do


    !     Formula (2.52)i of Shu, i.e. the 7-th order reconstruction


    vplus = 0.0_dp
    do r = 0 , 3
       vplus = vplus + omegat(r) * vmr(r)
    end do


  end subroutine wenorec


!> \brief WENO local Lax-Friedrich inviscid x-component flux.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine euler_LLF_x ( thd , adi , dx_i , x , y , z , T , W_i , cp , ha , shock_sensor , v , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (in)          :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    logical , allocatable , dimension (:,:,:) , intent (inout)     :: shock_sensor !< shock criterion
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< x-component flux


    integer (ip) , parameter                  :: stencil_m1 = ng+ng

    integer (ip)                              :: ok

    integer (ip)                              :: st , i_s , i_c , i_l , i_r , j , k , l , m , mm , s1 , s2

    real (dp)                                 :: ux_c , vy_c , wz_c , T_c , ht_c , et_c , W_i_c , gam_c , cs_c , cp_c

    real (dp) , dimension (nrv)               :: psi_c , ha_c

    real (dp) , dimension (nrv+npv+nvv)       :: Y_c

    real (dp) , dimension (nv)                :: ghat

    real (dp)                                 :: drho , dpres , eigenmax , df , gl , gr

    real (dp)                                 :: r , r1 , cs_c_i , cs_c_2 , b1 , b2 , q , xi , fl_W , fl_E

    real (dp)                                 :: wrk

    real (dp) , dimension (nv,nv)             :: er , el

    logical                                   :: shk , neg_pres

    integer (ip) , dimension (ndimmax)        :: neg_pres_coords

    real (dp) , dimension (stencil_m1,nv)     :: f_s , L_s

    real (dp) , dimension (stencil_m1)        :: gplus , gminus , wc , gc

    real (dp) , dimension (:,:) , allocatable :: fhat , Ya

    real (dp) , dimension (:) , allocatable   :: rho_i , gam , ux , vy , wz , cs , P , ht , et


    ! real (dp) , dimension (nv,nv) :: matrix         , jacob     , &
    !                                  eigenvalmatrix , jacobtest


    allocate ( fhat  (sx-1:ex,nv)            , &
               vy    (sx-1:ex+1)             , &
               wz    (sx-1:ex+1)             , &
               ht    (sx-1:ex+1)             , &
               et    (sx-1:ex+1)             , &
               Ya    (sx-1:ex+1,nrv+npv+nvv) , &
               ux    (sx-ng:ex+ng)           , &
               cs    (sx-ng:ex+ng)           , &
               P     (sx-ng:ex+ng)           , &
               gam   (sx-ng:ex+ng)           , &
               rho_i (sx-ng:ex+ng)           , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate euler_LLF_x'


    fl_W = 1.0_dp
    if ( bc (W) == noreflection ) fl_W = 0.0_dp
    fl_E = 1.0_dp
    if ( bc (E) == noreflection ) fl_E = 0.0_dp


    neg_pres = .false.


    do k = sz , ez ! loop in the z-direction
       do j = sy , ey ! loop in the y-direction


          ! presolving some variables (1/2)
          do i_c = sx-ng , ex+ng

             rho_i (i_c) = 1.0_dp / v (i_c,j,k,1)
             ux (i_c)    = v (i_c,j,k,2) * rho_i (i_c)
             P (i_c)     = v (i_c,j,k,1) * T (i_c,j,k) * W_i (i_c,j,k)
             if ( P (i_c) <= 0.0_dp .and. .not. neg_pres ) then
                neg_pres = .true.
                neg_pres_coords (1) = i_c
                neg_pres_coords (2) = j
                neg_pres_coords (3) = k
             end if
             gam (i_c) = thd % gam2 * W_i (i_c,j,k) / cp (i_c,j,k)
             gam (i_c) = 1.0_dp / ( 1.0_dp - gam (i_c) )
             cs (i_c)  = sqrt ( gam (i_c) * P (i_c) * rho_i (i_c) )

          end do


          ! presolving some variables (2/2)
          do i_c = sx-1 , ex+1

             vy (i_c) = v (i_c,j,k,3) * rho_i (i_c)
             wz (i_c) = v (i_c,j,k,4) * rho_i (i_c)
             do l = 1 , nrv+npv+nvv
                Ya (i_c,l) = v (i_c,j,k,niv+l) * rho_i (i_c)
             end do
             ht (i_c) = 0.0_dp
             do l = 1 , nrv
                ht (i_c) = ht (i_c) + ha (i_c,j,k,l) * Ya (i_c,l)
             end do
             ht (i_c) = ht (i_c) + 0.5_dp * ( ux(i_c)*ux(i_c) + vy(i_c)*vy(i_c) + wz(i_c)*wz(i_c) )
             et (i_c) = v (i_c,j,k,5) * rho_i (i_c)

          end do


          do i_l = sx-1 , ex ! loop on the cell faces


             ! shk = .false.
             i_r = i_l + 1


             do st = 1 , stencil_m1

                i_s = i_l + st - ng

                f_s (st,1) = v (i_s,j,k,2)
                f_s (st,2) = v (i_s,j,k,2) * ux (i_s) + P (i_s)
                f_s (st,3) = v (i_s,j,k,3) * ux (i_s)
                f_s (st,4) = v (i_s,j,k,4) * ux (i_s)
                f_s (st,5) = ux (i_s) * ( v (i_s,j,k,5) + P (i_s) )
                do l = niv+1 , nv
                   f_s (st,l) = v (i_s,j,k,l) * ux (i_s)
                end do

                L_s (st,1) = abs ( ux (i_s) - cs (i_s) )
                L_s (st,2) = abs ( ux (i_s) )
                L_s (st,3) = abs ( ux (i_s) + cs (i_s) )
                L_s (st,4) = L_s (st,2)
                L_s (st,5) = L_s (st,2)
                do l = niv+1 , nv
                   L_s (st,l) = L_s (st,2)
                end do

                ! ! density criteria
                ! drho = abs ( v (min(i_s+1,ex),j,k,1) - v (i_s,j,k,1) ) * &
                !        min ( rho_i (i_s) , rho_i (min(i_s+1,ex)) )

                ! ! pressure criteria
                ! dpres = abs ( P (min(i_s+1,ex)) - P (i_s) ) / &
                !         min ( P (i_s) , P (min(i_s+1,ex)) )

                ! if ( drho > max_rel_weight .and. dpres > max_rel_weight ) shk = .true.

             end do


!             ! non-reflecting conditions
!              if ( i_l < sx+ng .and. bc (W) == noreflection ) shk = .true.
!              if ( i_l > ex-ng .and. bc (E) == noreflection ) shk = .true.


! !             ! activate 2D WENO at the boundaries (not 3D because of periodicity)
!              if ( i_l < 1+ng )   shk = .true.
!              if ( i_l > ntx-ng ) shk = .true.
!              if ( j   < 1+ng )   shk = .true.
!              if ( j   > nty-ng ) shk = .true.


             ! Roe's average state (by definition: velocities, mass fractios and total enthalpy)
             r     = sqrt ( v (i_r,j,k,1) * rho_i (i_l) )
             r1    = 1.0_dp / ( r + 1.0_dp )
             ux_c  = ( r * ux (i_r)  + ux (i_l)  ) * r1
             vy_c  = ( r * vy (i_r)  + vy (i_l)  ) * r1
             wz_c  = ( r * wz (i_r)  + wz (i_l)  ) * r1
             do l = 1 , nrv+npv+nvv
                Y_c (l) = ( r * ( Ya (i_r,l) ) + Ya (i_l,l) ) * r1
             end do
             ht_c  = ( r * ht (i_r)  + ht (i_l)  ) * r1
             et_c  = ( r * et (i_r)  + et (i_l)  ) * r1 ! but also total energy (Dieterding)


             ! bad average
             ! T_c   = ( r * T (i_r,j,k) + T (i_l,j,k) ) * r1
             ! gam_c = ( r * gam (i_r) + gam (i_l) ) * r1
             ! cs_c  = ( r * cs (i_r)  + cs (i_l)  ) * r1
             ! wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T_c
             ! do l = 1 , nrv
             !    psi_c (l) = - ( r * ha (i_r,j,k,l) + ha (i_l,j,k,l) ) * r1
             !    psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             ! end do


             ! average state (Roe approx. compatible with EOS)
             call Wmix_i_scalar ( thd , Y_c , W_i_c )
             T_c = ( ht_c - et_c ) / W_i_c ! averaged temperature compatible with EOS and previous Roe's averages
             call cp_scalar ( thd , T_c , Y_c , cp_c ) ! cp compatible with EOS and previous Roe's averages
             gam_c = thd % gam2 * W_i_c / cp_c
             gam_c = 1.0_dp / ( 1.0_dp - gam_c ) ! gamma compatible with EOS
             cs_c  = sqrt ( gam_c * T_c * W_i_c ) ! sound speed compatible with EOS
             call ha_scalar ( thd , T_c , ha_c )  ! species enthalpies compatible with EOS
             wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T_c
             do l = 1 , nrv
                psi_c (l) = - ( r * ha_c (l) + ha_c (l) ) * r1
                psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             end do


             ! auxiliary variables to fill the matrices
             cs_c_i = 1.0_dp / cs_c
             cs_c_2 = cs_c * cs_c
             q      = 0.5_dp * ( ux_c*ux_c + vy_c*vy_c + wz_c*wz_c )
             b2     = ( gam_c - 1.0_dp ) / cs_c_2
             b1     = b2 * q
             xi     = b1 - b2 * ht_c

             ! left eigenvectors matrix (at Roe's state)

             ! 1 : FIRST
             el(1,1)   =   0.5_dp * ( b1        + ux_c * cs_c_i )
             el(1,2)   = - 0.5_dp * ( b2 * ux_c +        cs_c_i )
             el(1,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(1,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(1,5)   =   0.5_dp * b2

             el(2,1)   =   1.0_dp - b1
             el(2,2)   =   b2 * ux_c
             el(2,3)   =   b2 * vy_c
             el(2,4)   =   b2 * wz_c
             el(2,5)   = - b2

             el(3,1)   =   0.5_dp * ( b1        - ux_c * cs_c_i )
             el(3,2)   = - 0.5_dp * ( b2 * ux_c -        cs_c_i )
             el(3,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(3,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(3,5)   =   0.5_dp * b2

             el(4,1)   = - vy_c
             el(4,2)   =   0.0_dp
             el(4,3)   =   1.0_dp
             el(4,4)   =   0.0_dp
             el(4,5)   =   0.0_dp

             el(5,1)   = - wz_c
             el(5,2)   =   0.0_dp
             el(5,3)   =   0.0_dp
             el(5,4)   =   1.0_dp
             el(5,5)   =   0.0_dp

             ! 2 : SECOND
             wrk = 0.5_dp * b2
             do l = niv+1 , nv-npv-nvv
                el ( 1 , l ) = wrk * psi_c ( l-niv )
             end do
             el ( 2   , niv+1:nv-npv-nvv ) = - el ( 1 , niv+1:nv-npv-nvv ) - el ( 1 , niv+1:nv-npv-nvv )
             el ( 3   , niv+1:nv-npv-nvv ) =   el ( 1 , niv+1:nv-npv-nvv )
             el ( 4:5 , niv+1:nv-npv-nvv ) =   0.0_dp
             if ( npv > 0 ) then
                el ( 1:5 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             ! 3 : THIRD
             el(6,1)   =   ( 1 + xi ) * q
             el(6,2)   = - ( 1 + xi ) * ux_c
             el(6,3)   = - ( 1 + xi ) * vy_c
             el(6,4)   = - ( 1 + xi ) * wz_c
             el(6,5)   =   ( 1 + xi )
             do l = niv+1 , nv-npv-nvv
                el ( niv+1 , l ) = xi * psi_c ( l-niv )
             end do
             if ( npv > 0 ) then
                el ( niv+1 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             if ( nv-npv-nvv >= 7 ) then
                ! 4 : FOURTH
                do l = 7 , nv-npv-nvv
                   el(l,1) = - b1        * Y_c (l-6)
                   el(l,2) =   b2 * ux_c * Y_c (l-6)
                   el(l,3) =   b2 * vy_c * Y_c (l-6)
                   el(l,4) =   b2 * wz_c * Y_c (l-6)
                   el(l,5) = - b2        * Y_c (l-6)
                end do
                if ( npv > 0 ) then
                   do l = nv-npv-nvv+1 , nv
                      el(l,1) = - b1        * Y_c (l-5)
                      el(l,2) =   b2 * ux_c * Y_c (l-5)
                      el(l,3) =   b2 * vy_c * Y_c (l-5)
                      el(l,4) =   b2 * wz_c * Y_c (l-5)
                      el(l,5) = - b2        * Y_c (l-5)
                   end do
                end if

                ! 5 : FIFTH
                do s1 = 1 , nv-npv-nvv-niv-1
                   do s2 = 1 , nv-npv-nvv-niv
                      el ( s1+6 , s2+5 ) = - b2 * Y_c (s1) * psi_c (s2)
                   end do
                end do
                do l = 1 , nv-npv-nvv-niv-1
                   el ( l+6 , l+5 ) = el ( l+6 , l+5 ) + 1.0_dp
                end do
                if ( npv > 0 ) then
                   el (7:nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , nv-npv-nvv-niv
                         el ( s1+nv-npv-nvv , s2+5 ) = - b2 * Y_c (s1+nv-npv-nvv-niv) * psi_c (s2)
                      end do
                   end do
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , npv+nvv
                         el ( s1+nv-npv-nvv , s2+nv-npv-nvv ) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      el ( l+nv-npv-nvv , l+nv-npv-nvv ) = 1.0_dp
                   end do
                end if
             end if


             ! right eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             er(1,1)   =  1.0_dp
             er(1,2)   =  1.0_dp
             er(1,3)   =  1.0_dp
             er(1,4)   =  0.0_dp
             er(1,5)   =  0.0_dp

             er(2,1)   =  ux_c - cs_c
             er(2,2)   =  ux_c
             er(2,3)   =  ux_c + cs_c
             er(2,4)   =  0.0_dp
             er(2,5)   =  0.0_dp

             er(3,1)   =  vy_c
             er(3,2)   =  vy_c
             er(3,3)   =  vy_c
             er(3,4)   =  1.0_dp
             er(3,5)   =  0.0_dp

             er(4,1)   =  wz_c
             er(4,2)   =  wz_c
             er(4,3)   =  wz_c
             er(4,4)   =  0.0_dp
             er(4,5)   =  1.0_dp

             er(5,1)   =  ht_c - ux_c * cs_c
             er(5,2)   =  q
             er(5,3)   =  ht_c + ux_c * cs_c
             er(5,4)   =  vy_c
             er(5,5)   =  wz_c

             ! 2 : SECOND
             er(1:4,niv+1:nv) = 0.0_dp

             ! 3 : THIRD
             do l = niv+1 , nv
                er(l,1) =  Y_c (l-5)
                er(l,2) =  0.0_dp
                er(l,3) =  Y_c (l-5)
                er(l,4) =  0.0_dp
                er(l,5) =  0.0_dp
             end do

             ! 4 : FOURTH
             do s1 = niv , nv-npv-nvv
                do s2 = niv+1 , nv
                   er(s1,s2) = 0.0_dp
                end do
             end do
             do l = niv+1 , nv-npv-nvv
                er(l-1,l) = 1.0_dp
             end do

             ! 5 : FIFTH
             er ( nv-npv-nvv , 6 ) = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             if ( nv-npv-nvv >= 7 ) then
                wrk = - 1.0_dp / psi_c (nv-npv-nvv-niv)
                do l = 1 , nv-npv-nvv-niv-1
                   er (nv-npv-nvv,l+niv+1) = wrk * psi_c (l)
                end do
                if ( npv > 0 ) then
                   er (nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = nv-npv-nvv+1 , nv
                      do s2 = niv+1 , nv
                         er(s1,s2) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      er(l+nv-npv-nvv,l+nv-npv-nvv) = 1.0_dp
                   end do
                end if
             end if

             do m = 1 , nv ! loop on the five char. fields

                eigenmax = -1.0_dp
                do st = 1 , stencil_m1 ! LLF
                   eigenmax = max ( L_s (st,m) , eigenmax )
                end do

                do st = 1 , stencil_m1 ! loop over the stencil centered at face i
                   wc(st) = 0.0_dp
                   gc(st) = 0.0_dp
                   i_s = i_l + st - ng
                   do mm = 1 , nv
                      wc(st) = wc(st) + el(m,mm) * v (i_s,j,k,mm)
                      gc(st) = gc(st) + el(m,mm) * f_s (st,mm)
                   end do
                   gplus (st) = 0.5_dp * ( gc(st) + eigenmax * wc(st) )
                   gminus(st) = gc(st) - gplus (st)
                end do

                ! Reconstruction of the '+' and '-' fluxes (page 32)
                if (any(shock_sensor(i_l-ng+1:i_l+ng,j,k))  .and. weno_avg ) then
                   call wenorec    ( gplus , gminus , gl , gr )
                else
                   call wenorec_nw ( gplus , gminus , gl , gr )
                end if

                ghat(m) = gl + gr ! char. flux

             end do

             ! Evaluation of fhat: the aim of this loop
             do m = 1 , nv
                fhat (i_l,m) = 0.0_dp
                do mm = 1 , nv
                   fhat (i_l,m) = fhat(i_l,m) + er(m,mm) * ghat(mm)
                end do
             end do


          end do ! end of loop on the cell faces


          ! evaluation of the flux


          do m = 1 , nv
             do i_c = sx , sx
                i_l = i_c - 1
                df = ( fhat(i_c,m) - fhat(i_l,m) ) * dx_i (i_c)
                fl (i_c,j,k,m) = df * fl_W ! instead of fl = fl + df
             end do
          end do

          do m = 1 , nv
             do i_c = sx+1 , ex-1 ! loop on the inner nodes
                i_l = i_c - 1
                df = ( fhat(i_c,m) - fhat(i_l,m) ) * dx_i (i_c)
                fl (i_c,j,k,m) = df ! instead of fl = fl + df
             end do
          end do

          do m = 1 , nv
             do i_c = ex , ex
                i_l = i_c - 1
                df = ( fhat(i_c,m) - fhat(i_l,m) ) * dx_i (i_c)
                fl (i_c,j,k,m) = df * fl_E ! instead of fl = fl + df
             end do
          end do


       end do ! end of j-loop
    end do ! end of k-loop


    if (neg_pres) then
       write (*,'(1X,A,I9,5(1X,I10))') 'WARNING: negative pressure at X-dir (rank,i,j,k)   =' , &
                    rank, neg_pres_coords(1) , neg_pres_coords(2) , neg_pres_coords(3)
       write (*,'(42X,A,9X,5(1X,1PE10.3))') '(x,y,z)   =' , &
                                      x (neg_pres_coords(1)) * adi % L_ref , &
                                      y (neg_pres_coords(2)) * adi % L_ref , &
                                      z (neg_pres_coords(3)) * adi % L_ref

       call abort_mpi (' ')
    end if


    deallocate ( fhat , vy , wz , ht , Ya , et , &
                 ux , cs , P , gam , rho_i )


  end subroutine euler_LLF_x


!> \brief WENO local Lax-Friedrich inviscid y-component flux.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine euler_LLF_y ( thd , adi , dy_i , x , y , z , T , W_i , cp , ha , shock_sensor , v , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (in)          :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    logical , allocatable , dimension (:,:,:) , intent (inout)     :: shock_sensor !< shock criterion
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< y-component flux


    integer (ip) , parameter                  :: stencil_m1 = ng+ng

    integer (ip)                              :: ok

    integer (ip)                              :: st , j_s , j_c , j_l , j_r , i , k , l , m , mm , s1 , s2

    real (dp)                                 :: ux_c , vy_c , wz_c , T_c , ht_c , et_c , W_i_c , gam_c , cs_c , cp_c

    real (dp) , dimension (nrv)               :: psi_c , ha_c

    real (dp) , dimension (nrv+npv+nvv)       :: Y_c

    real (dp) , dimension (nv)                :: ghat

    real (dp)                                 :: drho , dpres , eigenmax , df , gl , gr

    real (dp)                                 :: r , r1 , cs_c_i , cs_c_2 , b1 , b2 , q , xi , fl_S , fl_N

    real (dp)                                 :: wrk

    real (dp) , dimension (nv,nv)             :: er , el

    logical                                   :: shk , neg_pres

    integer (ip) , dimension (ndimmax)        :: neg_pres_coords

    real (dp) , dimension (stencil_m1,nv)     :: f_s , L_s

    real (dp) , dimension (stencil_m1)        :: gplus , gminus , wc , gc

    real (dp) , dimension (:,:) , allocatable :: fhat , Ya

    real (dp) , dimension (:) , allocatable   :: rho_i , gam , ux , vy , wz , cs , P , ht , et


    ! real (dp) , dimension (nv,nv) :: matrix         , jacob     , &
    !                                  eigenvalmatrix , jacobtest


    allocate ( fhat  (sy-1:ey,nv)            , &
               ux    (sy-1:ey+1)             , &
               wz    (sy-1:ey+1)             , &
               ht    (sy-1:ey+1)             , &
               et    (sy-1:ey+1)             , &
               Ya    (sy-1:ey+1,nrv+npv+nvv) , &
               vy    (sy-ng:ey+ng)           , &
               cs    (sy-ng:ey+ng)           , &
               P     (sy-ng:ey+ng)           , &
               gam   (sy-ng:ey+ng)           , &
               rho_i (sy-ng:ey+ng)           , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate euler_LLF_y'


    fl_S = 1.0_dp
    if ( bc (S) == noreflection ) fl_S = 0.0_dp
    fl_N = 1.0_dp
    if ( bc (N) == noreflection ) fl_N = 0.0_dp


    neg_pres = .false.


    do k = sz , ez ! loop in the z-direction
       do i = sx , ex ! loop in the x-direction


          ! presolving some variables (1/2)
          do j_c = sy-ng , ey+ng

             rho_i (j_c) = 1.0_dp / v (i,j_c,k,1)
             vy (j_c)    = v (i,j_c,k,3) * rho_i (j_c)
             P (j_c)     = v (i,j_c,k,1) * T (i,j_c,k) * W_i (i,j_c,k)
             if ( P (j_c) <= 0.0_dp .and. .not. neg_pres ) then
                neg_pres = .true.
                neg_pres_coords (1) = i
                neg_pres_coords (2) = j_c
                neg_pres_coords (3) = k
             end if
             gam (j_c) = thd % gam2 * W_i (i,j_c,k) / cp (i,j_c,k)
             gam (j_c) = 1.0_dp / ( 1.0_dp - gam (j_c) )
             cs (j_c)  = sqrt ( gam (j_c) * P (j_c) * rho_i (j_c) )

          end do


          ! presolving some variables (2/2)
          do j_c = sy-1 , ey+1

             ux (j_c) = v (i,j_c,k,2) * rho_i (j_c)
             wz (j_c) = v (i,j_c,k,4) * rho_i (j_c)
             do l = 1 , nrv+npv+nvv
                Ya (j_c,l) = v (i,j_c,k,niv+l) * rho_i (j_c)
             end do
             ht (j_c) = 0.0_dp
             do l = 1 , nrv
                ht (j_c) = ht (j_c) + ha (i,j_c,k,l) * Ya (j_c,l)
             end do
             ht (j_c) = ht (j_c) + 0.5_dp * ( ux(j_c)*ux(j_c) + vy(j_c)*vy(j_c) + wz(j_c)*wz(j_c) )
             et (j_c) = v (i,j_c,k,5) * rho_i (j_c)

          end do


          do j_l = sy-1 , ey ! loop on the cell faces


             ! shk = .false.
             j_r = j_l + 1


             do st = 1 , stencil_m1

                j_s = j_l + st - ng

                f_s (st,1) = v (i,j_s,k,3)
                f_s (st,2) = v (i,j_s,k,2) * vy (j_s)
                f_s (st,3) = v (i,j_s,k,3) * vy (j_s) + P (j_s)
                f_s (st,4) = v (i,j_s,k,4) * vy (j_s)
                f_s (st,5) = vy (j_s) * ( v (i,j_s,k,5) + P (j_s) )
                do l = niv+1 , nv
                   f_s (st,l) = v (i,j_s,k,l) * vy (j_s)
                end do

                L_s (st,1) = abs ( vy (j_s) - cs (j_s) )
                L_s (st,2) = abs ( vy (j_s) )
                L_s (st,3) = abs ( vy (j_s) + cs (j_s) )
                L_s (st,4) = L_s (st,2)
                L_s (st,5) = L_s (st,2)
                do l = niv+1 , nv
                   L_s (st,l) = L_s (st,2)
                end do

                ! ! density criteria
                ! drho = abs ( v (i,min(j_s+1,ey),k,1) - v (i,j_s,k,1) ) * &
                !        min ( rho_i (j_s) , rho_i (min(j_s+1,ey)) )

                ! ! pressure criteria
                ! dpres = abs ( P (min(j_s+1,ey)) - P (j_s) ) / &
                !         min ( P (j_s) , P (min(j_s+1,ey)) )

                ! if ( drho > max_rel_weight .and. dpres > max_rel_weight ) shk = .true.

             end do


!             ! non-reflecting conditions
!              if ( j_l <= sy+ng .and. bc (S) == noreflection ) shk = .true.
!              if ( j_l >= ey-ng .and. bc (N) == noreflection ) shk = .true.


! !             ! activate 2D WENO at the boundaries (not 3D because of periodicity)
!              if ( i   < 1+ng )   shk = .true.
!              if ( i   > ntx-ng ) shk = .true.
!              if ( j_l < 1+ng )   shk = .true.
!              if ( j_l > nty-ng ) shk = .true.


             ! Roe's average state (by definition: velocities, mass fractios and total enthalpy)
             r     = sqrt ( v (i,j_r,k,1) * rho_i (j_l) )
             r1    = 1.0_dp / ( r + 1.0_dp )
             ux_c  = ( r * ux (j_r)  + ux (j_l)  ) * r1
             vy_c  = ( r * vy (j_r)  + vy (j_l)  ) * r1
             wz_c  = ( r * wz (j_r)  + wz (j_l)  ) * r1
             do l = 1 , nrv+npv+nvv
                Y_c (l) = ( r * ( Ya (j_r,l) ) + Ya (j_l,l) ) * r1
             end do
             ht_c  = ( r * ht (j_r)  + ht (j_l)  ) * r1
             et_c  = ( r * et (j_r)  + et (j_l)  ) * r1 ! but also total energy (Dieterding)


             ! bad average
             ! T_c   = ( r * T (i,j_r,k) + T (i,j_l,k) ) * r1
             ! gam_c = ( r * gam (j_r) + gam (j_l) ) * r1
             ! cs_c  = ( r * cs (j_r)  + cs (j_l)  ) * r1
             ! wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T_c
             ! do l = 1 , nrv
             !    psi_c (l) = - ( r * ha (i,j_r,k,l) + ha (i,j_l,k,l) ) * r1
             !    psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             ! end do


             ! average state (Roe approx. compatible with EOS)
             call Wmix_i_scalar ( thd , Y_c , W_i_c )
             T_c = ( ht_c - et_c ) / W_i_c ! averaged temperature compatible with EOS and previous Roe's averages
             call cp_scalar ( thd , T_c , Y_c , cp_c ) ! cp compatible with EOS and previous Roe's averages
             gam_c = thd % gam2 * W_i_c / cp_c
             gam_c = 1.0_dp / ( 1.0_dp - gam_c ) ! gamma compatible with EOS
             cs_c  = sqrt ( gam_c * T_c * W_i_c ) ! sound speed compatible with EOS
             call ha_scalar ( thd , T_c , ha_c )  ! species enthalpies compatible with EOS
             wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T_c
             do l = 1 , nrv
                psi_c (l) = - ( r * ha_c (l) + ha_c (l) ) * r1
                psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             end do


             ! auxiliary variables to fill the matrices
             cs_c_i = 1.0_dp / cs_c
             cs_c_2 = cs_c * cs_c
             q      = 0.5_dp * ( ux_c*ux_c + vy_c*vy_c + wz_c*wz_c )
             b2     = ( gam_c - 1.0_dp ) / cs_c_2
             b1     = b2 * q
             xi     = b1 - b2 * ht_c


             ! write (*,*) ux_c , vy_c , wz_c
             ! write (*,*)
             ! write (*,*) T_c , gam_c , cs_c  , ht_c , cs_c_i , cs_c_2 , q , b1 , b2 , xi
             ! write (*,*)
             ! write (*,*) Y_c (:)
             ! write (*,*)
             ! write (*,*) ha (i,j_l,k,:)
             ! write (*,*)
             ! write (*,*) psi_c
             ! write (*,*)
             ! write (*,*) f_s (1,:)
             ! write (*,*)
             ! read(*,*)


             ! left eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             el(1,1)   =   0.5_dp * ( b1        + vy_c * cs_c_i )
             el(1,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(1,3)   = - 0.5_dp * ( b2 * vy_c +        cs_c_i )
             el(1,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(1,5)   =   0.5_dp * b2

             el(2,1)   =   1.0_dp - b1
             el(2,2)   =   b2 * ux_c
             el(2,3)   =   b2 * vy_c
             el(2,4)   =   b2 * wz_c
             el(2,5)   = - b2

             el(3,1)   =   0.5_dp * ( b1        - vy_c * cs_c_i )
             el(3,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(3,3)   = - 0.5_dp * ( b2 * vy_c -        cs_c_i )
             el(3,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(3,5)   =   0.5_dp * b2

             el(4,1)   = - ux_c
             el(4,2)   =   1.0_dp
             el(4,3)   =   0.0_dp
             el(4,4)   =   0.0_dp
             el(4,5)   =   0.0_dp

             el(5,1)   = - wz_c
             el(5,2)   =   0.0_dp
             el(5,3)   =   0.0_dp
             el(5,4)   =   1.0_dp
             el(5,5)   =   0.0_dp

             ! 2 : SECOND
             wrk = 0.5_dp * b2
             do l = niv+1 , nv-npv-nvv
                el ( 1 , l ) = wrk * psi_c ( l-niv )
             end do
             el ( 2   , niv+1:nv-npv-nvv ) = - el ( 1 , niv+1:nv-npv-nvv ) - el ( 1 , niv+1:nv-npv-nvv )
             el ( 3   , niv+1:nv-npv-nvv ) =   el ( 1 , niv+1:nv-npv-nvv )
             el ( 4:5 , niv+1:nv-npv-nvv ) =   0.0_dp
             if ( npv > 0 ) then
                el ( 1:5 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             ! 3 : THIRD
             el(6,1)   =   ( 1 + xi ) * q
             el(6,2)   = - ( 1 + xi ) * ux_c
             el(6,3)   = - ( 1 + xi ) * vy_c
             el(6,4)   = - ( 1 + xi ) * wz_c
             el(6,5)   =   ( 1 + xi )
             do l = niv+1 , nv-npv-nvv
                el ( niv+1 , l ) = xi * psi_c ( l-niv )
             end do
             if ( npv > 0 ) then
                el ( niv+1 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             if ( nv-npv-nvv >= 7 ) then
                ! 4 : FOURTH
                do l = 7 , nv-npv-nvv
                   el(l,1) = - b1        * Y_c (l-6)
                   el(l,2) =   b2 * ux_c * Y_c (l-6)
                   el(l,3) =   b2 * vy_c * Y_c (l-6)
                   el(l,4) =   b2 * wz_c * Y_c (l-6)
                   el(l,5) = - b2        * Y_c (l-6)
                end do
                if ( npv > 0 ) then
                   do l = nv-npv-nvv+1 , nv
                      el(l,1) = - b1        * Y_c (l-5)
                      el(l,2) =   b2 * ux_c * Y_c (l-5)
                      el(l,3) =   b2 * vy_c * Y_c (l-5)
                      el(l,4) =   b2 * wz_c * Y_c (l-5)
                      el(l,5) = - b2        * Y_c (l-5)
                   end do
                end if

                ! 5 : FIFTH
                do s1 = 1 , nv-npv-nvv-niv-1
                   do s2 = 1 , nv-npv-nvv-niv
                      el ( s1+6 , s2+5 ) = - b2 * Y_c (s1) * psi_c (s2)
                   end do
                end do
                do l = 1 , nv-npv-nvv-niv-1
                   el ( l+6 , l+5 ) = el ( l+6 , l+5 ) + 1.0_dp
                end do
                if ( npv > 0 ) then
                   el (7:nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , nv-npv-nvv-niv
                         el ( s1+nv-npv-nvv , s2+5 ) = - b2 * Y_c (s1+nv-npv-nvv-niv) * psi_c (s2)
                      end do
                   end do
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , npv+nvv
                         el ( s1+nv-npv-nvv , s2+nv-npv-nvv ) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      el ( l+nv-npv-nvv , l+nv-npv-nvv ) = 1.0_dp
                   end do
                end if
             end if


             ! right eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             er(1,1)   =  1.0_dp
             er(1,2)   =  1.0_dp
             er(1,3)   =  1.0_dp
             er(1,4)   =  0.0_dp
             er(1,5)   =  0.0_dp

             er(2,1)   =  ux_c
             er(2,2)   =  ux_c
             er(2,3)   =  ux_c
             er(2,4)   =  1.0_dp
             er(2,5)   =  0.0_dp

             er(3,1)   =  vy_c - cs_c
             er(3,2)   =  vy_c
             er(3,3)   =  vy_c + cs_c
             er(3,4)   =  0.0_dp
             er(3,5)   =  0.0_dp

             er(4,1)   =  wz_c
             er(4,2)   =  wz_c
             er(4,3)   =  wz_c
             er(4,4)   =  0.0_dp
             er(4,5)   =  1.0_dp

             er(5,1)   =  ht_c - vy_c * cs_c
             er(5,2)   =  q
             er(5,3)   =  ht_c + vy_c * cs_c
             er(5,4)   =  ux_c
             er(5,5)   =  wz_c

             ! 2 : SECOND
             er(1:4,niv+1:nv) = 0.0_dp

             ! 3 : THIRD
             do l = niv+1 , nv
                er(l,1) =  Y_c (l-5)
                er(l,2) =  0.0_dp
                er(l,3) =  Y_c (l-5)
                er(l,4) =  0.0_dp
                er(l,5) =  0.0_dp
             end do

             ! 4 : FOURTH
             do s1 = niv , nv-npv-nvv
                do s2 = niv+1 , nv
                   er(s1,s2) = 0.0_dp
                end do
             end do
             do l = niv+1 , nv-npv-nvv
                er(l-1,l) = 1.0_dp
             end do

             ! 5 : FIFTH
             er ( nv-npv-nvv , 6 ) = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             if ( nv-npv-nvv >= 7 ) then
                wrk = - 1.0_dp / psi_c (nv-npv-nvv-niv)
                do l = 1 , nv-npv-nvv-niv-1
                   er (nv-npv-nvv,l+niv+1) = wrk * psi_c (l)
                end do
                if ( npv > 0 ) then
                   er (nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = nv-npv-nvv+1 , nv
                      do s2 = niv+1 , nv
                         er(s1,s2) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      er(l+nv-npv-nvv,l+nv-npv-nvv) = 1.0_dp
                   end do
                end if
             end if





             do m = 1 , nv ! loop on the five char. fields

                eigenmax = -1.0_dp
                do st = 1 , stencil_m1 ! LLF
                   eigenmax = max ( L_s (st,m) , eigenmax )
                end do

                do st = 1 , stencil_m1 ! loop over the stencil centered at face i
                   wc(st) = 0.0_dp
                   gc(st) = 0.0_dp
                   j_s = j_l + st - ng
                   do mm = 1 , nv
                      wc(st) = wc(st) + el(m,mm) * v (i,j_s,k,mm)
                      gc(st) = gc(st) + el(m,mm) * f_s (st,mm)
                   end do
                   gplus (st) = 0.5_dp * ( gc(st) + eigenmax * wc(st) )
                   gminus(st) = gc(st) - gplus (st)
                end do

                ! Reconstruction of the '+' and '-' fluxes (page 32)
                if (any(shock_sensor(i,j_l-ng+1:j_l+ng,k))  .and. weno_avg ) then
                   call wenorec    ( gplus , gminus , gl , gr )
                else
                   call wenorec_nw ( gplus , gminus , gl , gr )
                end if

                ghat(m) = gl + gr ! char. flux

             end do


             ! Evaluation of fhat: the aim of this loop
             do m = 1 , nv
                fhat (j_l,m) = 0.0_dp
                do mm = 1 , nv
                   fhat (j_l,m) = fhat(j_l,m) + er(m,mm) * ghat(mm)
                end do
             end do


          end do ! end of loop on the cell faces


          ! evaluation of the flux


          do m = 1 , nv
             do j_c = sy , sy
                j_l = j_c - 1
                df = ( fhat(j_c,m) - fhat(j_l,m) ) * dy_i (j_c)
                fl (i,j_c,k,m) = fl (i,j_c,k,m) + df * fl_S
             end do
          end do

          do m = 1 , nv
             do j_c = sy+1 , ey-1 ! loop on the inner nodes
                j_l = j_c - 1
                df = ( fhat(j_c,m) - fhat(j_l,m) ) * dy_i (j_c)
                fl (i,j_c,k,m) = fl (i,j_c,k,m) + df
             end do
          end do

          do m = 1 , nv
             do j_c = ey , ey
                j_l = j_c - 1
                df = ( fhat(j_c,m) - fhat(j_l,m) ) * dy_i (j_c)
                fl (i,j_c,k,m) = fl (i,j_c,k,m) + df * fl_N
             end do
          end do


       end do ! end of j-loop
    end do ! end of k-loop


    if (neg_pres) then
       write (*,'(1X,A,I9,5(1X,I10))') 'WARNING: negative pressure at X-dir (rank,i,j,k)   =' , &
                    rank, neg_pres_coords(1) , neg_pres_coords(2) , neg_pres_coords(3)
       write (*,'(42X,A,9X,5(1X,1PE10.3))') '(x,y,z)   =' , &
                                      x (neg_pres_coords(1)) * adi % L_ref , &
                                      y (neg_pres_coords(2)) * adi % L_ref , &
                                      z (neg_pres_coords(3)) * adi % L_ref

       call abort_mpi (' ')
    end if


    deallocate ( fhat , vy , wz , ht , Ya , et , &
                 ux , cs , P , gam , rho_i )


  end subroutine euler_LLF_y


!> \brief WENO local Lax-Friedrich inviscid z-component flux.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine euler_LLF_z ( thd , adi , dz_i , x , y , z , T , W_i , cp , ha , shock_sensor , v , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i !< inverted dz array
    real (dp) , allocatable , dimension (:) , intent (in)          :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    logical , allocatable , dimension (:,:,:) , intent (inout)     :: shock_sensor !< shock criterion
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< z-component flux


    integer (ip) , parameter                  :: stencil_m1 = ng+ng

    integer (ip)                              :: ok

    integer (ip)                              :: st , k_s , k_c , k_l , k_r , i , j , l , m , mm , s1 , s2

    real (dp)                                 :: ux_c , vy_c , wz_c , T_c , ht_c , et_c , W_i_c , gam_c , cs_c , cp_c

    real (dp) , dimension (nrv)               :: psi_c , ha_c

    real (dp) , dimension (nrv+npv+nvv)       :: Y_c

    real (dp) , dimension (nv)                :: ghat

    real (dp)                                 :: drho , dpres , eigenmax , df , gl , gr

    real (dp)                                 :: r , r1 , cs_c_i , cs_c_2 , b1 , b2 , q , xi , fl_B , fl_F

    real (dp)                                 :: wrk

    real (dp) , dimension (nv,nv)             :: er , el

    logical                                   :: shk , neg_pres

    integer (ip) , dimension (ndimmax)        :: neg_pres_coords

    real (dp) , dimension (stencil_m1,nv)     :: f_s , L_s

    real (dp) , dimension (stencil_m1)        :: gplus , gminus , wc , gc

    real (dp) , dimension (:,:) , allocatable :: fhat , Ya

    real (dp) , dimension (:) , allocatable   :: rho_i , gam , ux , vy , wz , cs , P , ht , et


    ! real (dp) , dimension (nv,nv) :: matrix         , jacob     , &
    !                                  eigenvalmatrix , jacobtest


    allocate ( fhat  (sz-1:ez,nv)            , &
               ux    (sz-1:ez+1)             , &
               vy    (sz-1:ez+1)             , &
               et    (sz-1:ez+1)             , &
               ht    (sz-1:ez+1)             , &
               Ya    (sz-1:ez+1,nrv+npv+nvv) , &
               wz    (sz-ng:ez+ng)           , &
               cs    (sz-ng:ez+ng)           , &
               P     (sz-ng:ez+ng)           , &
               gam   (sz-ng:ez+ng)           , &
               rho_i (sz-ng:ez+ng)           , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate euler_LLF_z'


    fl_B = 1.0_dp
    if ( bc (B) == noreflection ) fl_B = 0.0_dp
    fl_F = 1.0_dp
    if ( bc (F) == noreflection ) fl_F = 0.0_dp


    neg_pres = .false.


    do j = sy , ey ! loop in the y-direction
       do i = sx , ex ! loop in the x-direction

          ! presolving some variables (1/2)
          do k_c = sz-ng , ez+ng

             rho_i (k_c) = 1.0_dp / v (i,j,k_c,1)
             wz (k_c)    = v (i,j,k_c,4) * rho_i (k_c)
             P (k_c)     = v (i,j,k_c,1) * T (i,j,k_c) * W_i (i,j,k_c)
             if ( P (k_c) <= 0.0_dp .and. .not. neg_pres ) then
                neg_pres = .true.
                neg_pres_coords (1) = i
                neg_pres_coords (2) = j
                neg_pres_coords (3) = k_c
             end if
             gam (k_c) = thd % gam2 * W_i (i,j,k_c) / cp (i,j,k_c)
             gam (k_c) = 1.0_dp / ( 1.0_dp - gam (k_c) )
             cs (k_c)  = sqrt ( gam (k_c) * P (k_c) * rho_i (k_c) )

          end do


          ! presolving some variables (2/2)
          do k_c = sz-1 , ez+1

             ux (k_c) = v (i,j,k_c,2) * rho_i (k_c)
             vy (k_c) = v (i,j,k_c,3) * rho_i (k_c)
             do l = 1 , nrv+npv+nvv
                Ya (k_c,l) = v (i,j,k_c,niv+l) * rho_i (k_c)
             end do
             ht (k_c) = 0.0_dp
             do l = 1 , nrv
                ht (k_c) = ht (k_c) + ha (i,j,k_c,l) * Ya (k_c,l)
             end do
             ht (k_c) = ht (k_c) + 0.5_dp * ( ux(k_c)*ux(k_c) + vy(k_c)*vy(k_c) + wz(k_c)*wz(k_c) )
             et (k_c) = v (i,j,k_c,5) * rho_i (k_c)

          end do


          do k_l = sz-1 , ez ! loop on the cell faces


             shk = .false.
             k_r = k_l + 1


             do st = 1 , stencil_m1

                k_s = k_l + st - ng

                f_s (st,1) = v (i,j,k_s,4)
                f_s (st,2) = v (i,j,k_s,2) * wz (k_s)
                f_s (st,3) = v (i,j,k_s,3) * wz (k_s)
                f_s (st,4) = v (i,j,k_s,4) * wz (k_s) + P (k_s)
                f_s (st,5) = wz (k_s) * ( v (i,j,k_s,5) + P (k_s) )
                do l = niv+1 , nv
                   f_s (st,l) = v (i,j,k_s,l) * wz (k_s)
                end do

                L_s (st,1) = abs ( wz (k_s) - cs (k_s) )
                L_s (st,2) = abs ( wz (k_s) )
                L_s (st,3) = abs ( wz (k_s) + cs (k_s) )
                L_s (st,4) = L_s (st,2)
                L_s (st,5) = L_s (st,2)
                do l = niv+1 , nv
                   L_s (st,l) = L_s (st,2)
                end do

                ! ! density criteria
                ! drho = abs ( v (i,j,min(k_s+1,ez),1) - v (i,j,k_s,1) ) * &
                !        min ( rho_i (k_s) , rho_i (min(k_s+1,ez)) )

                ! ! pressure criteria
                ! dpres = abs ( P (min(k_s+1,ez)) - P (k_s) ) / &
                !         min ( P (k_s) , P (min(k_s+1,ez)) )

                ! if ( drho > max_rel_weight .and. dpres > max_rel_weight ) shk = .true.

             end do


!             ! non-reflecting conditions
!              if ( k_l < sz+ng .and. bc (B) == noreflection ) shk = .true.
!              if ( k_l > ez-ng .and. bc (F) == noreflection ) shk = .true.


             ! ! activate 2D WENO at the boundaries (not 3D because of periodicity)
             ! if ( i < 1+ng )   shk = .true.
             ! if ( i > ntx-ng ) shk = .true.
             ! if ( j < 1+ng )   shk = .true.
             ! if ( j > nty-ng ) shk = .true.


             ! Roe's average state (by definition: velocities, mass fractios and total enthalpy)
             r     = sqrt ( v (i,j,k_r,1) * rho_i (k_l) )
             r1    = 1.0_dp / ( r + 1.0_dp )
             ux_c  = ( r * ux (k_r)  + ux (k_l)  ) * r1
             vy_c  = ( r * vy (k_r)  + vy (k_l)  ) * r1
             wz_c  = ( r * wz (k_r)  + wz (k_l)  ) * r1
             do l = 1 , nrv+npv+nvv
                Y_c (l) = ( r * ( Ya (k_r,l) ) + Ya (k_l,l) ) * r1
             end do
             ht_c  = ( r * ht (k_r)  + ht (k_l)  ) * r1
             et_c  = ( r * et (k_r)  + et (k_l)  ) * r1 ! but also total energy (Dieterding)


             ! bad average
             ! T_c   = ( r * T (i,j,k_r) + T (i,j,k_l) ) * r1
             ! gam_c = ( r * gam (k_r) + gam (k_l) ) * r1
             ! cs_c  = ( r * cs (k_r)  + cs (k_l)  ) * r1
             ! wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T_c
             ! do l = 1 , nrv
             !    psi_c (l) = - ( r * ha (i,j,k_r,l) + ha (i,j,k_l,l) ) * r1
             !    psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             ! end do


             ! average state (Roe approx. compatible with EOS)
             call Wmix_i_scalar ( thd , Y_c , W_i_c )
             T_c = ( ht_c - et_c ) / W_i_c ! averaged temperature compatible with EOS and previous Roe's averages
             call cp_scalar ( thd , T_c , Y_c , cp_c ) ! cp compatible with EOS and previous Roe's averages
             gam_c = thd % gam2 * W_i_c / cp_c
             gam_c = 1.0_dp / ( 1.0_dp - gam_c ) ! gamma compatible with EOS
             cs_c  = sqrt ( gam_c * T_c * W_i_c ) ! sound speed compatible with EOS
             call ha_scalar ( thd , T_c , ha_c )  ! species enthalpies compatible with EOS
             wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T_c
             do l = 1 , nrv
                psi_c (l) = - ( r * ha_c (l) + ha_c (l) ) * r1
                psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             end do


             cs_c_i = 1.0_dp / cs_c
             cs_c_2 = cs_c * cs_c
             q      = 0.5_dp * ( ux_c*ux_c + vy_c*vy_c + wz_c*wz_c )
             b2     = ( gam_c - 1.0_dp ) / cs_c_2
             b1     = b2 * q
             xi     = b1 - b2 * ht_c


             ! write (*,*) ux_c , vy_c , wz_c
             ! write (*,*)
             ! write (*,*) T_c , gam_c , cs_c  , ht_c , cs_c_i , cs_c_2 , q , b1 , b2 , xi
             ! write (*,*)
             ! write (*,*) Y_c (:)
             ! write (*,*)
             ! write (*,*) ha (i,j,k_l,:)
             ! write (*,*)
             ! write (*,*) psi_c
             ! write (*,*)
             ! write (*,*) f_s (1,:)
             ! write (*,*)
             ! read(*,*)


             ! left eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             el(1,1)   =   0.5_dp * ( b1        + wz_c * cs_c_i )
             el(1,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(1,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(1,4)   = - 0.5_dp * ( b2 * wz_c +        cs_c_i )
             el(1,5)   =   0.5_dp * b2

             el(2,1)   =   1.0_dp - b1
             el(2,2)   =   b2 * ux_c
             el(2,3)   =   b2 * vy_c
             el(2,4)   =   b2 * wz_c
             el(2,5)   = - b2

             el(3,1)   =   0.5_dp * ( b1        - wz_c * cs_c_i )
             el(3,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(3,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(3,4)   = - 0.5_dp * ( b2 * wz_c -        cs_c_i )
             el(3,5)   =   0.5_dp * b2

             el(4,1)   = - ux_c
             el(4,2)   =   1.0_dp
             el(4,3)   =   0.0_dp
             el(4,4)   =   0.0_dp
             el(4,5)   =   0.0_dp

             el(5,1)   = - vy_c
             el(5,2)   =   0.0_dp
             el(5,3)   =   1.0_dp
             el(5,4)   =   0.0_dp
             el(5,5)   =   0.0_dp

             ! 2 : SECOND
             wrk = 0.5_dp * b2
             do l = niv+1 , nv-npv-nvv
                el ( 1 , l ) = wrk * psi_c ( l-niv )
             end do
             el ( 2   , niv+1:nv-npv-nvv ) = - el ( 1 , niv+1:nv-npv-nvv ) - el ( 1 , niv+1:nv-npv-nvv )
             el ( 3   , niv+1:nv-npv-nvv ) =   el ( 1 , niv+1:nv-npv-nvv )
             el ( 4:5 , niv+1:nv-npv-nvv ) =   0.0_dp
             if ( npv > 0 ) then
                el ( 1:5 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             ! 3 : THIRD
             el(6,1)   =   ( 1 + xi ) * q
             el(6,2)   = - ( 1 + xi ) * ux_c
             el(6,3)   = - ( 1 + xi ) * vy_c
             el(6,4)   = - ( 1 + xi ) * wz_c
             el(6,5)   =   ( 1 + xi )
             do l = niv+1 , nv-npv-nvv
                el ( niv+1 , l ) = xi * psi_c ( l-niv )
             end do
             if ( npv > 0 ) then
                el ( niv+1 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             if ( nv-npv-nvv >= 7 ) then
                ! 4 : FOURTH
                do l = 7 , nv-npv-nvv
                   el(l,1) = - b1        * Y_c (l-6)
                   el(l,2) =   b2 * ux_c * Y_c (l-6)
                   el(l,3) =   b2 * vy_c * Y_c (l-6)
                   el(l,4) =   b2 * wz_c * Y_c (l-6)
                   el(l,5) = - b2        * Y_c (l-6)
                end do
                if ( npv > 0 ) then
                   do l = nv-npv-nvv+1 , nv
                      el(l,1) = - b1        * Y_c (l-5)
                      el(l,2) =   b2 * ux_c * Y_c (l-5)
                      el(l,3) =   b2 * vy_c * Y_c (l-5)
                      el(l,4) =   b2 * wz_c * Y_c (l-5)
                      el(l,5) = - b2        * Y_c (l-5)
                   end do
                end if

                ! 5 : FIFTH
                do s1 = 1 , nv-npv-nvv-niv-1
                   do s2 = 1 , nv-npv-nvv-niv
                      el ( s1+6 , s2+5 ) = - b2 * Y_c (s1) * psi_c (s2)
                   end do
                end do
                do l = 1 , nv-npv-nvv-niv-1
                   el ( l+6 , l+5 ) = el ( l+6 , l+5 ) + 1.0_dp
                end do
                if ( npv > 0 ) then
                   el (7:nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , nv-npv-nvv-niv
                         el ( s1+nv-npv-nvv , s2+5 ) = - b2 * Y_c (s1+nv-npv-nvv-niv) * psi_c (s2)
                      end do
                   end do
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , npv+nvv
                         el ( s1+nv-npv-nvv , s2+nv-npv-nvv ) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      el ( l+nv-npv-nvv , l+nv-npv-nvv ) = 1.0_dp
                   end do
                end if
             end if


             ! right eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             er(1,1)   =  1.0_dp
             er(1,2)   =  1.0_dp
             er(1,3)   =  1.0_dp
             er(1,4)   =  0.0_dp
             er(1,5)   =  0.0_dp

             er(2,1)   =  ux_c
             er(2,2)   =  ux_c
             er(2,3)   =  ux_c
             er(2,4)   =  1.0_dp
             er(2,5)   =  0.0_dp

             er(3,1)   =  vy_c
             er(3,2)   =  vy_c
             er(3,3)   =  vy_c
             er(3,4)   =  0.0_dp
             er(3,5)   =  1.0_dp

             er(4,1)   =  wz_c - cs_c
             er(4,2)   =  wz_c
             er(4,3)   =  wz_c + cs_c
             er(4,4)   =  0.0_dp
             er(4,5)   =  0.0_dp

             er(5,1)   =  ht_c - wz_c * cs_c
             er(5,2)   =  q
             er(5,3)   =  ht_c + wz_c * cs_c
             er(5,4)   =  ux_c
             er(5,5)   =  vy_c

             ! 2 : SECOND
             er(1:4,niv+1:nv) = 0.0_dp

             ! 3 : THIRD
             do l = niv+1 , nv
                er(l,1) =  Y_c (l-5)
                er(l,2) =  0.0_dp
                er(l,3) =  Y_c (l-5)
                er(l,4) =  0.0_dp
                er(l,5) =  0.0_dp
             end do

             ! 4 : FOURTH
             do s1 = niv , nv-npv-nvv
                do s2 = niv+1 , nv
                   er(s1,s2) = 0.0_dp
                end do
             end do
             do l = 6 , nv-npv-nvv
                er(l-1,l) = 1.0_dp
             end do

             ! 5 : FIFTH
             er ( nv-npv-nvv , 6 ) = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             if ( nv-npv-nvv >= 7 ) then
                wrk = - 1.0_dp / psi_c (nv-npv-nvv-niv)
                do l = 1 , nv-npv-nvv-niv-1
                   er (nv-npv-nvv,l+niv+1) = wrk * psi_c (l)
                end do
                if ( npv > 0 ) then
                   er (nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = nv-npv-nvv+1 , nv
                      do s2 = niv+1 , nv
                         er(s1,s2) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      er(l+nv-npv-nvv,l+nv-npv-nvv) = 1.0_dp
                   end do
                end if
             end if

             do m = 1 , nv ! loop on the five char. fields

                eigenmax = -1.0_dp
                do st = 1 , stencil_m1 ! LLF
                   eigenmax = max ( L_s (st,m) , eigenmax )
                end do

                do st = 1 , stencil_m1 ! loop over the stencil centered at face i
                   wc(st) = 0.0_dp
                   gc(st) = 0.0_dp
                   k_s = k_l + st - ng
                   do mm = 1 , nv
                      wc(st) = wc(st) + el(m,mm) * v (i,j,k_s,mm)
                      gc(st) = gc(st) + el(m,mm) * f_s (st,mm)
                   end do
                   gplus (st) = 0.5_dp * ( gc(st) + eigenmax * wc(st) )
                   gminus(st) = gc(st) - gplus (st)
                end do

                ! Reconstruction of the '+' and '-' fluxes (page 32)
                if (any(shock_sensor(i,j,k_l-ng+1:k_l+ng))  .and. weno_avg ) then
                   call wenorec    ( gplus , gminus , gl , gr )
                else
                   call wenorec_nw ( gplus , gminus , gl , gr )
                end if

                ghat(m) = gl + gr ! char. flux

             end do


             ! Evaluation of fhat: the aim of this loop
             do m = 1 , nv
                fhat (k_l,m) = 0.0_dp
                do mm = 1 , nv
                   fhat (k_l,m) = fhat(k_l,m) + er(m,mm) * ghat(mm)
                end do
             end do


          end do ! end of loop on the cell faces


          ! evaluation of the flux


          do m = 1 , nv
             do k_c = sz , sz
                k_l = k_c - 1
                df = ( fhat(k_c,m) - fhat(k_l,m) ) * dz_i (k_c)
                fl (i,j,k_c,m) = fl (i,j,k_c,m) + df * fl_B
             end do
          end do

          do m = 1 , nv
             do k_c = sz+1 , ez-1 ! loop on the inner nodes
                k_l = k_c - 1
                df = ( fhat(k_c,m) - fhat(k_l,m) ) * dz_i (k_c)
                fl (i,j,k_c,m) = fl (i,j,k_c,m) + df ! instead of fl = fl + df
             end do
          end do

          do m = 1 , nv
             do k_c = ez , ez
                k_l = k_c - 1
                df = ( fhat(k_c,m) - fhat(k_l,m) ) * dz_i (k_c)
                fl (i,j,k_c,m) = fl (i,j,k_c,m) + df * fl_F
             end do
          end do


       end do ! end of j-loop
    end do ! end of k-loop


    if (neg_pres) then
       write (*,'(1X,A,I9,5(1X,I10))') 'WARNING: negative pressure at X-dir (rank,i,j,k)   =' , &
                    rank, neg_pres_coords(1) , neg_pres_coords(2) , neg_pres_coords(3)
       write (*,'(42X,A,9X,5(1X,1PE10.3))') '(x,y,z)   =' , &
                                      x (neg_pres_coords(1)) * adi % L_ref , &
                                      y (neg_pres_coords(2)) * adi % L_ref , &
                                      z (neg_pres_coords(3)) * adi % L_ref

       call abort_mpi (' ')
    end if


    deallocate ( fhat , vy , wz , ht , Ya , et , &
                 ux , cs , P , gam , rho_i )


  end subroutine euler_LLF_z


!> \brief Modify WENO inviscid fluxes in a non-reflecting BC.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine bc_noreflectionweno ( face , thd , dx_i , dy_i , dz_i , T , W_i , cp , ha , v , fl )


    integer (ip) , intent (in)                                     :: face !< face of the BC (W, E...)
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< modified flux


    integer (ip) , parameter                  :: stencil_m1 = ng+ng , stencil_m2 = ng+ng-1

    integer (ip)                              :: st , i_s , j_s , k_s , i_c , j_c , k_c , i , j , k , l , m , mm , s1 , s2

    real (dp)                                 :: ux_c , vy_c , wz_c , ht_c , gam_c , cs_c , cs_c_i , cs_c_2

    real (dp) , dimension (nrv)               :: psi_c

    real (dp) , dimension (nrv+npv+nvv)       :: Y_c

    real (dp) , dimension (nv)                :: L_s , dw , dwc

    real (dp)                                 :: P , rho_i , b1 , b2 , q , xi , df

    real (dp)                                 :: wrk

    real (dp) , dimension (nv,nv)             :: er , el


    ! real (dp) , dimension (:,:) , allocatable :: fhat , Ya

    ! real (dp) , dimension (:) , allocatable   :: rho_i , gam , ux , vy , wz , cs , P , ht

    ! real (dp) , dimension (nv,nv) :: matrix         , jacob     , &
    !                                  eigenvalmatrix , jacobtest



    if ( face == W ) then


       i_c = sx


       do k = sz , ez ! loop in the z-direction
          do j = sy , ey ! loop in the y-direction


             ! presolving some variables (1/2)
             rho_i = 1.0_dp / v (i_c,j,k,1)
             ux_c  = v (i_c,j,k,2) * rho_i
             P     = v (i_c,j,k,1) * T (i_c,j,k) * W_i (i_c,j,k)
             gam_c = thd % gam2 * W_i (i_c,j,k) / cp (i_c,j,k)
             gam_c = 1.0_dp / ( 1.0_dp - gam_c )
             cs_c  = sqrt ( gam_c * P * rho_i )

             ! presolving some variables (2/2)
             vy_c = v (i_c,j,k,3) * rho_i
             wz_c = v (i_c,j,k,4) * rho_i
             do l = 1 , nrv+npv+nvv
                Y_c (l) = v (i_c,j,k,niv+l) * rho_i
             end do
             ht_c = 0.0_dp
             wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T (i_c,j,k)
             do l = 1 , nrv
                ht_c      = ht_c + ha (i_c,j,k,l) * Y_c (l)
                psi_c (l) = - ha (i_c,j,k,l)
                psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             end do


             cs_c_i = 1.0_dp / cs_c
             cs_c_2 = cs_c * cs_c
             q      = 0.5_dp * ( ux_c*ux_c + vy_c*vy_c + wz_c*wz_c )
             ht_c   = ht_c + q
             b2     = ( gam_c - 1.0_dp ) / cs_c_2
             b1     = b2 * q
             xi     = b1 - b2 * ht_c


             ! write (*,*) ux_c , vy_c , wz_c
             ! write (*,*)
             ! write (*,*) T_c , gam_c , cs_c  , ht_c , cs_c_i , cs_c_2 , q , b1 , b2 , xi
             ! write (*,*)
             ! write (*,*) Y_c (:)
             ! write (*,*)
             ! write (*,*) ha (i,j,k_l,:)
             ! write (*,*)
             ! write (*,*) psi_c
             ! write (*,*)
             ! write (*,*) f_s (1,:)
             ! write (*,*)
             ! read(*,*)


             ! left eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             el(1,1)   =   0.5_dp * ( b1        + ux_c * cs_c_i )
             el(1,2)   = - 0.5_dp * ( b2 * ux_c +        cs_c_i )
             el(1,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(1,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(1,5)   =   0.5_dp * b2

             el(2,1)   =   1.0_dp - b1
             el(2,2)   =   b2 * ux_c
             el(2,3)   =   b2 * vy_c
             el(2,4)   =   b2 * wz_c
             el(2,5)   = - b2

             el(3,1)   =   0.5_dp * ( b1        - ux_c * cs_c_i )
             el(3,2)   = - 0.5_dp * ( b2 * ux_c -        cs_c_i )
             el(3,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(3,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(3,5)   =   0.5_dp * b2

             el(4,1)   = - vy_c
             el(4,2)   =   0.0_dp
             el(4,3)   =   1.0_dp
             el(4,4)   =   0.0_dp
             el(4,5)   =   0.0_dp

             el(5,1)   = - wz_c
             el(5,2)   =   0.0_dp
             el(5,3)   =   0.0_dp
             el(5,4)   =   1.0_dp
             el(5,5)   =   0.0_dp

             ! 2 : SECOND
             wrk = 0.5_dp * b2
             do l = niv+1 , nv-npv-nvv
                el ( 1 , l ) = wrk * psi_c ( l-niv )
             end do
             el ( 2   , niv+1:nv-npv-nvv ) = - el ( 1 , niv+1:nv-npv-nvv ) - el ( 1 , niv+1:nv-npv-nvv )
             el ( 3   , niv+1:nv-npv-nvv ) =   el ( 1 , niv+1:nv-npv-nvv )
             el ( 4:5 , niv+1:nv-npv-nvv ) =   0.0_dp
             if ( npv > 0 ) then
                el ( 1:5 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             ! 3 : THIRD
             el(6,1)   =   ( 1 + xi ) * q
             el(6,2)   = - ( 1 + xi ) * ux_c
             el(6,3)   = - ( 1 + xi ) * vy_c
             el(6,4)   = - ( 1 + xi ) * wz_c
             el(6,5)   =   ( 1 + xi )
             do l = niv+1 , nv-npv-nvv
                el ( niv+1 , l ) = xi * psi_c ( l-niv )
             end do
             if ( npv > 0 ) then
                el ( niv+1 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             if ( nv-npv-nvv >= 7 ) then
                ! 4 : FOURTH
                do l = 7 , nv-npv-nvv
                   el(l,1) = - b1        * Y_c (l-6)
                   el(l,2) =   b2 * ux_c * Y_c (l-6)
                   el(l,3) =   b2 * vy_c * Y_c (l-6)
                   el(l,4) =   b2 * wz_c * Y_c (l-6)
                   el(l,5) = - b2        * Y_c (l-6)
                end do
                if ( npv > 0 ) then
                   do l = nv-npv-nvv+1 , nv
                      el(l,1) = - b1        * Y_c (l-5)
                      el(l,2) =   b2 * ux_c * Y_c (l-5)
                      el(l,3) =   b2 * vy_c * Y_c (l-5)
                      el(l,4) =   b2 * wz_c * Y_c (l-5)
                      el(l,5) = - b2        * Y_c (l-5)
                   end do
                end if

                ! 5 : FIFTH
                do s1 = 1 , nv-npv-nvv-niv-1
                   do s2 = 1 , nv-npv-nvv-niv
                      el ( s1+6 , s2+5 ) = - b2 * Y_c (s1) * psi_c (s2)
                   end do
                end do
                do l = 1 , nv-npv-nvv-niv-1
                   el ( l+6 , l+5 ) = el ( l+6 , l+5 ) + 1.0_dp
                end do
                if ( npv > 0 ) then
                   el (7:nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , nv-npv-nvv-niv
                         el ( s1+nv-npv-nvv , s2+5 ) = - b2 * Y_c (s1+nv-npv-nvv-niv) * psi_c (s2)
                      end do
                   end do
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , npv+nvv
                         el ( s1+nv-npv-nvv , s2+nv-npv-nvv ) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      el ( l+nv-npv-nvv , l+nv-npv-nvv ) = 1.0_dp
                   end do
                end if
             end if


             ! right eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             er(1,1)   =  1.0_dp
             er(1,2)   =  1.0_dp
             er(1,3)   =  1.0_dp
             er(1,4)   =  0.0_dp
             er(1,5)   =  0.0_dp

             er(2,1)   =  ux_c - cs_c
             er(2,2)   =  ux_c
             er(2,3)   =  ux_c + cs_c
             er(2,4)   =  0.0_dp
             er(2,5)   =  0.0_dp

             er(3,1)   =  vy_c
             er(3,2)   =  vy_c
             er(3,3)   =  vy_c
             er(3,4)   =  1.0_dp
             er(3,5)   =  0.0_dp

             er(4,1)   =  wz_c
             er(4,2)   =  wz_c
             er(4,3)   =  wz_c
             er(4,4)   =  0.0_dp
             er(4,5)   =  1.0_dp

             er(5,1)   =  ht_c - ux_c * cs_c
             er(5,2)   =  q
             er(5,3)   =  ht_c + ux_c * cs_c
             er(5,4)   =  vy_c
             er(5,5)   =  wz_c

             ! 2 : SECOND
             er(1:4,niv+1:nv) = 0.0_dp

             ! 3 : THIRD
             do l = niv+1 , nv
                er(l,1) =  Y_c (l-5)
                er(l,2) =  0.0_dp
                er(l,3) =  Y_c (l-5)
                er(l,4) =  0.0_dp
                er(l,5) =  0.0_dp
             end do

             ! 4 : FOURTH
             do s1 = niv , nv-npv-nvv
                do s2 = niv+1 , nv
                   er(s1,s2) = 0.0_dp
                end do
             end do
             do l = niv+1 , nv-npv-nvv
                er(l-1,l) = 1.0_dp
             end do

             ! 5 : FIFTH
             er ( nv-npv-nvv , 6 ) = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             wrk = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             if ( nv-npv-nvv >= 7 ) then
                do l = 1 , nv-npv-nvv-niv-1
                   er (nv-npv-nvv,l+niv+1) = wrk * psi_c (l)
                end do
                if ( npv > 0 ) then
                   er (nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = nv-npv-nvv+1 , nv
                      do s2 = niv+1 , nv
                         er(s1,s2) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      er(l+nv-npv-nvv,l+nv-npv-nvv) = 1.0_dp
                   end do
                end if
             end if



!          ! TESTING MATRICES
           !  if (rank==0) then
           !      write (*,*) 'euler noreflecting'
           !      matrix = matmul ( er , el )
           !      write (*,*) 'er*el'
           !      do l = 1 , nv
           !         write (*,'(50(1X,1PE15.6))') matrix ( l , : )
           !      end do
           ! end if
         ! do j1 = 1 , nvar
         !    do i1 = 1 , nvar
         !       if ( i1==j1 ) then
         !          if ( abs (matrix ( i1 , j1 ) - 1.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       else
         !          if ( abs (matrix ( i1 , j1 ) - 0.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       end if
         !    end do
         ! end do
!          eigenvalmatrix = 0.
!          do s1 = 1 , nvar
!             eigenvalmatrix ( s1 , s1 ) = uu
!          end do
!          eigenvalmatrix ( 1 , 1 ) = eigenvalmatrix ( 1 , 1 ) - c
!          eigenvalmatrix ( 3 , 3 ) = eigenvalmatrix ( 3 , 3 ) + c
! !
!          write (*,*) 'eigenvalues'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') eigenvalmatrix ( i1 , : )
!          end do
! !
! !
!          jacobtest = 0.
!          jacobtest = matmul ( eigenvalmatrix , el )
!          el        = matmul ( er , jacobtest )
!          jacobtest = el
!          write (*,*) 'calculated jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacobtest ( i1 , : )
!          end do
! !
! !
!          jacob = 0.
!          jacob ( 1,1 ) = 0.
!          jacob ( 1,2 ) = 1.
!          jacob ( 1 , 3:nvar ) = 0.
! !
!          jacob ( 2 , 1 ) = - ( ( 1. - gmm ) * qq + uu*uu)
!          jacob ( 2 , 2 ) =   ( 3. - gmm ) * uu
!          jacob ( 2 , 3 ) =   ( 1. - gmm ) * vv
!          jacob ( 2 , 4 ) =   ( 1. - gmm ) * ww
!          jacob ( 2 , 5 ) = - ( 1. - gmm )
!          do s1 = 1 , nvar-5
!             jacob ( 2 , s1+5 ) = - ( 1. - gmm ) * psim ( s1 )
!          end do
! !
!          jacob ( 3 , 1 ) = - uu * vv
!          jacob ( 3 , 2 ) = vv
!          jacob ( 3 , 3 ) = uu
!          jacob ( 3 , 4:nvar ) = 0.
! !
!          jacob ( 4 , 1 ) = - uu * ww
!          jacob ( 4 , 2 ) = ww
!          jacob ( 4 , 3 ) = 0.
!          jacob ( 4 , 4 ) = uu
!          jacob ( 4 , 5:nvar ) = 0.
! !
!          jacob ( 5 , 1 ) = - ( ( 1. - gmm ) * qq + h ) * uu
!          jacob ( 5 , 2 ) =   ( 1. - gmm ) * uu * uu + h
!          jacob ( 5 , 3 ) =   ( 1. - gmm ) * uu * vv
!          jacob ( 5 , 4 ) =   ( 1. - gmm ) * uu * ww
!          jacob ( 5 , 5 ) = gmm * uu
!          do s1 = 1 , nvar-5
!             jacob ( 5 , s1+5 ) = - ( 1. - gmm ) * uu * psim ( s1 )
!          end do
! !
!          do s1 = 1 , nvar-5
!             jacob ( 5 + s1 , 1 ) = - uu * chim ( s1 )
!             jacob ( 5 + s1 , 2 ) =        chim ( s1 )
!             jacob ( 5 + s1 , 3:5 ) = 0.
!          end do
! !
!          do s1 = 6 , nvar
!             jacob ( s1 , s1 ) = uu
!             ! the rest of the elements are set to zero in the
!             ! initialisation of the jacobian
!          end do
! !
!          write (*,*) 'analytical jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacob ( i1 , : )
!          end do
!          write (*,*) 'End of matrices'
!          read (*,*)


             ! Eigenvalues (killing the positive ones)
             L_s (1) =    min ( ux_c - cs_c , 0.0_dp )
             L_s (2) =    min ( ux_c        , 0.0_dp )
             L_s (3) =    min ( ux_c + cs_c , 0.0_dp )
             L_s (4) =    min ( ux_c        , 0.0_dp )
             L_s (5) =    min ( ux_c        , 0.0_dp )
             do m = niv+1 , nv
                L_s (m) = min ( ux_c        , 0.0_dp )
             end do


             ! Derivatives of conservative variables
             do m = 1 , nv
                dw (m) = 0.0_dp
                do st = 1 , stencil_m2
                   i_s = i_c + st - 1 ! be careful
                   dw (m) = dw (m) + cnr (st) * v (i_s,j,k,m)
                end do
             end do


             ! Projection on the eigendirections
             do m = 1 , nv
                dwc (m) = 0.0_dp
                do mm = 1 , nv
                   dwc (m) = dwc (m) + el (m,mm) * dw (mm)
                end do
             end do


             ! Returning to conservative variables
             do m = 1 , nv
                df = 0.0_dp
                do mm = 1 , nv
                   df = df + er (m,mm) * L_s (mm) * dwc (mm)
                end do
                fl (i_c,j,k,m) = fl (i_c,j,k,m) + df * dx_i (i_c)
             end do


          end do ! end of j-loop
       end do ! end of k-loop


    else if ( face == E ) then


       i_c = ex


       do k = sz , ez ! loop in the z-direction
          do j = sy , ey ! loop in the y-direction


             ! presolving some variables (1/2)
             rho_i = 1.0_dp / v (i_c,j,k,1)
             ux_c  = v (i_c,j,k,2) * rho_i
             P     = v (i_c,j,k,1) * T (i_c,j,k) * W_i (i_c,j,k)
             gam_c = thd % gam2 * W_i (i_c,j,k) / cp (i_c,j,k)
             gam_c = 1.0_dp / ( 1.0_dp - gam_c )
             cs_c  = sqrt ( gam_c * P * rho_i )

             ! presolving some variables (2/2)
             vy_c = v (i_c,j,k,3) * rho_i
             wz_c = v (i_c,j,k,4) * rho_i
             do l = 1 , nrv+npv+nvv
                Y_c (l) = v (i_c,j,k,niv+l) * rho_i
             end do
             ht_c = 0.0_dp
             wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T (i_c,j,k)
             do l = 1 , nrv
                ht_c = ht_c + ha (i_c,j,k,l) * Y_c (l)
                psi_c (l) = - ha (i_c,j,k,l)
                psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             end do


             cs_c_i = 1.0_dp / cs_c
             cs_c_2 = cs_c * cs_c
             q      = 0.5_dp * ( ux_c*ux_c + vy_c*vy_c + wz_c*wz_c )
             ht_c   = ht_c + q
             b2     = ( gam_c - 1.0_dp ) / cs_c_2
             b1     = b2 * q
             xi     = b1 - b2 * ht_c

             ! write (*,*)
             ! write (*,*) '***********************************'
             ! write (*,*)
             ! write (*,*) ux_c , vy_c , wz_c
             ! write (*,*)
!             write (*,*) T(i_c,j,k) , gam_c , cs_c  , ht_c , cs_c_i , cs_c_2 , q , b1 , b2 , xi
             ! write (*,*)
!              write (*,*) Y_c (:)
             ! write (*,*)
!             write (*,*) ha (i_c,j,k,:)
             ! write (*,*)
!             write (*,*) psi_c
             ! write (*,*)
             ! write (*,*) f_s (1,:)
             ! write (*,*)
             ! read(*,*)


             ! left eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             el(1,1)   =   0.5_dp * ( b1        + ux_c * cs_c_i )
             el(1,2)   = - 0.5_dp * ( b2 * ux_c +        cs_c_i )
             el(1,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(1,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(1,5)   =   0.5_dp * b2

             el(2,1)   =   1.0_dp - b1
             el(2,2)   =   b2 * ux_c
             el(2,3)   =   b2 * vy_c
             el(2,4)   =   b2 * wz_c
             el(2,5)   = - b2

             el(3,1)   =   0.5_dp * ( b1        - ux_c * cs_c_i )
             el(3,2)   = - 0.5_dp * ( b2 * ux_c -        cs_c_i )
             el(3,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(3,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(3,5)   =   0.5_dp * b2

             el(4,1)   = - vy_c
             el(4,2)   =   0.0_dp
             el(4,3)   =   1.0_dp
             el(4,4)   =   0.0_dp
             el(4,5)   =   0.0_dp

             el(5,1)   = - wz_c
             el(5,2)   =   0.0_dp
             el(5,3)   =   0.0_dp
             el(5,4)   =   1.0_dp
             el(5,5)   =   0.0_dp

             ! 2 : SECOND
             wrk = 0.5_dp * b2
             do l = niv+1 , nv-npv-nvv
                el ( 1 , l ) = wrk * psi_c ( l-niv )
             end do
             el ( 2   , niv+1:nv-npv-nvv ) = - el ( 1 , niv+1:nv-npv-nvv ) - el ( 1 , niv+1:nv-npv-nvv )
             el ( 3   , niv+1:nv-npv-nvv ) =   el ( 1 , niv+1:nv-npv-nvv )
             el ( 4:5 , niv+1:nv-npv-nvv ) =   0.0_dp
             if ( npv > 0 ) then
                el ( 1:5 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             ! 3 : THIRD
             el(6,1)   =   ( 1 + xi ) * q
             el(6,2)   = - ( 1 + xi ) * ux_c
             el(6,3)   = - ( 1 + xi ) * vy_c
             el(6,4)   = - ( 1 + xi ) * wz_c
             el(6,5)   =   ( 1 + xi )
             do l = niv+1 , nv-npv-nvv
                el ( niv+1 , l ) = xi * psi_c ( l-niv )
             end do
             if ( npv > 0 ) then
                el ( niv+1 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             if ( nv-npv-nvv >= 7 ) then
                ! 4 : FOURTH
                do l = 7 , nv-npv-nvv
                   el(l,1) = - b1        * Y_c (l-6)
                   el(l,2) =   b2 * ux_c * Y_c (l-6)
                   el(l,3) =   b2 * vy_c * Y_c (l-6)
                   el(l,4) =   b2 * wz_c * Y_c (l-6)
                   el(l,5) = - b2        * Y_c (l-6)
                end do
                if ( npv > 0 ) then
                   do l = nv-npv-nvv+1 , nv
                      el(l,1) = - b1        * Y_c (l-5)
                      el(l,2) =   b2 * ux_c * Y_c (l-5)
                      el(l,3) =   b2 * vy_c * Y_c (l-5)
                      el(l,4) =   b2 * wz_c * Y_c (l-5)
                      el(l,5) = - b2        * Y_c (l-5)
                   end do
                end if

                ! 5 : FIFTH
                do s1 = 1 , nv-npv-nvv-niv-1
                   do s2 = 1 , nv-npv-nvv-niv
                      el ( s1+6 , s2+5 ) = - b2 * Y_c (s1) * psi_c (s2)
                   end do
                end do
                do l = 1 , nv-npv-nvv-niv-1
                   el ( l+6 , l+5 ) = el ( l+6 , l+5 ) + 1.0_dp
                end do
                if ( npv > 0 ) then
                   el (7:nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , nv-npv-nvv-niv
                         el ( s1+nv-npv-nvv , s2+5 ) = - b2 * Y_c (s1+nv-npv-nvv-niv) * psi_c (s2)
                      end do
                   end do
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , npv+nvv
                         el ( s1+nv-npv-nvv , s2+nv-npv-nvv ) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      el ( l+nv-npv-nvv , l+nv-npv-nvv ) = 1.0_dp
                   end do
                end if
             end if


             ! right eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             er(1,1)   =  1.0_dp
             er(1,2)   =  1.0_dp
             er(1,3)   =  1.0_dp
             er(1,4)   =  0.0_dp
             er(1,5)   =  0.0_dp

             er(2,1)   =  ux_c - cs_c
             er(2,2)   =  ux_c
             er(2,3)   =  ux_c + cs_c
             er(2,4)   =  0.0_dp
             er(2,5)   =  0.0_dp

             er(3,1)   =  vy_c
             er(3,2)   =  vy_c
             er(3,3)   =  vy_c
             er(3,4)   =  1.0_dp
             er(3,5)   =  0.0_dp

             er(4,1)   =  wz_c
             er(4,2)   =  wz_c
             er(4,3)   =  wz_c
             er(4,4)   =  0.0_dp
             er(4,5)   =  1.0_dp

             er(5,1)   =  ht_c - ux_c * cs_c
             er(5,2)   =  q
             er(5,3)   =  ht_c + ux_c * cs_c
             er(5,4)   =  vy_c
             er(5,5)   =  wz_c

             ! 2 : SECOND
             er(1:4,niv+1:nv) = 0.0_dp

             ! 3 : THIRD
             do l = niv+1 , nv
                er(l,1) =  Y_c (l-5)
                er(l,2) =  0.0_dp
                er(l,3) =  Y_c (l-5)
                er(l,4) =  0.0_dp
                er(l,5) =  0.0_dp
             end do

             ! 4 : FOURTH
             do s1 = niv , nv-npv-nvv
                do s2 = niv+1 , nv
                   er(s1,s2) = 0.0_dp
                end do
             end do
             do l = niv+1 , nv-npv-nvv
                er(l-1,l) = 1.0_dp
             end do

             ! 5 : FIFTH
             er ( nv-npv-nvv , 6 ) = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             if ( nv-npv-nvv >= 7 ) then
                wrk = - 1.0_dp / psi_c (nv-npv-nvv-niv)
                do l = 1 , nv-npv-nvv-niv-1
                   er (nv-npv-nvv,l+niv+1) = wrk * psi_c (l)
                end do
                if ( npv > 0 ) then
                   er (nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = nv-npv-nvv+1 , nv
                      do s2 = niv+1 , nv
                         er(s1,s2) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      er(l+nv-npv-nvv,l+nv-npv-nvv) = 1.0_dp
                   end do
                end if
             end if



!          ! TESTING MATRICES
           !  if (rank==0) then
           !     write (*,*) 'euler noreflecting '
           !      matrix = matmul ( er , el )
           !      write (*,*) 'er*el'
           !      do l = 1 , nv
           !         write (*,'(50(1X,1PE15.6))') matrix ( l , : )
           !      end do
           ! end if
           ! stop
         ! do j1 = 1 , nvar
         !    do i1 = 1 , nvar
         !       if ( i1==j1 ) then
         !          if ( abs (matrix ( i1 , j1 ) - 1.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       else
         !          if ( abs (matrix ( i1 , j1 ) - 0.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       end if
         !    end do
         ! end do
!          eigenvalmatrix = 0.
!          do s1 = 1 , nvar
!             eigenvalmatrix ( s1 , s1 ) = uu
!          end do
!          eigenvalmatrix ( 1 , 1 ) = eigenvalmatrix ( 1 , 1 ) - c
!          eigenvalmatrix ( 3 , 3 ) = eigenvalmatrix ( 3 , 3 ) + c
! !
!          write (*,*) 'eigenvalues'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') eigenvalmatrix ( i1 , : )
!          end do
! !
! !
!          jacobtest = 0.
!          jacobtest = matmul ( eigenvalmatrix , el )
!          el        = matmul ( er , jacobtest )
!          jacobtest = el
!          write (*,*) 'calculated jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacobtest ( i1 , : )
!          end do
! !
! !
!          jacob = 0.
!          jacob ( 1,1 ) = 0.
!          jacob ( 1,2 ) = 1.
!          jacob ( 1 , 3:nvar ) = 0.
! !
!          jacob ( 2 , 1 ) = - ( ( 1. - gmm ) * qq + uu*uu)
!          jacob ( 2 , 2 ) =   ( 3. - gmm ) * uu
!          jacob ( 2 , 3 ) =   ( 1. - gmm ) * vv
!          jacob ( 2 , 4 ) =   ( 1. - gmm ) * ww
!          jacob ( 2 , 5 ) = - ( 1. - gmm )
!          do s1 = 1 , nvar-5
!             jacob ( 2 , s1+5 ) = - ( 1. - gmm ) * psim ( s1 )
!          end do
! !
!          jacob ( 3 , 1 ) = - uu * vv
!          jacob ( 3 , 2 ) = vv
!          jacob ( 3 , 3 ) = uu
!          jacob ( 3 , 4:nvar ) = 0.
! !
!          jacob ( 4 , 1 ) = - uu * ww
!          jacob ( 4 , 2 ) = ww
!          jacob ( 4 , 3 ) = 0.
!          jacob ( 4 , 4 ) = uu
!          jacob ( 4 , 5:nvar ) = 0.
! !
!          jacob ( 5 , 1 ) = - ( ( 1. - gmm ) * qq + h ) * uu
!          jacob ( 5 , 2 ) =   ( 1. - gmm ) * uu * uu + h
!          jacob ( 5 , 3 ) =   ( 1. - gmm ) * uu * vv
!          jacob ( 5 , 4 ) =   ( 1. - gmm ) * uu * ww
!          jacob ( 5 , 5 ) = gmm * uu
!          do s1 = 1 , nvar-5
!             jacob ( 5 , s1+5 ) = - ( 1. - gmm ) * uu * psim ( s1 )
!          end do
! !
!          do s1 = 1 , nvar-5
!             jacob ( 5 + s1 , 1 ) = - uu * chim ( s1 )
!             jacob ( 5 + s1 , 2 ) =        chim ( s1 )
!             jacob ( 5 + s1 , 3:5 ) = 0.
!          end do
! !
!          do s1 = 6 , nvar
!             jacob ( s1 , s1 ) = uu
!             ! the rest of the elements are set to zero in the
!             ! initialisation of the jacobian
!          end do
! !
!          write (*,*) 'analytical jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacob ( i1 , : )
!          end do
!          write (*,*) 'End of matrices'
!          read (*,*)


             ! Eigenvalues (killing the positive ones)
             L_s (1) =    max ( ux_c - cs_c , 0.0_dp )
             L_s (2) =    max ( ux_c        , 0.0_dp )
             L_s (3) =    max ( ux_c + cs_c , 0.0_dp )
             L_s (4) =    max ( ux_c        , 0.0_dp )
             L_s (5) =    max ( ux_c        , 0.0_dp )
             do m = niv+1 , nv
                L_s (m) = max ( ux_c        , 0.0_dp )
             end do


             ! Derivatives of conservative variables
             do m = 1 , nv
                dw (m) = 0.0_dp
                do st = 1 , stencil_m2
                   i_s = i_c - st + 1 ! be careful
                   dw (m) = dw (m) - cnr (st) * v (i_s,j,k,m) ! be careful
                end do
             end do

             ! Projection on the eigendirections
             do m = 1 , nv
                dwc (m) = 0.0_dp
                do mm = 1 , nv
                   dwc (m) = dwc (m) + el (m,mm) * dw (mm)
                end do
             end do


             ! Returning to conservative variables
             do m = 1 , nv
                df = 0.0_dp
                do mm = 1 , nv
                   df = df + er (m,mm) * L_s (mm) * dwc (mm)
                end do
                fl (i_c,j,k,m) = fl (i_c,j,k,m) + df * dx_i (i_c)
             end do


          end do ! end of j-loop
       end do ! end of k-loop


    else if ( face == S ) then


       j_c = sy


       do k = sz , ez ! loop in the z-direction
          do i = sx , ex ! loop in the x-direction


             ! presolving some variables (1/2)
             rho_i = 1.0_dp / v (i,j_c,k,1)
             ux_c  = v (i,j_c,k,2) * rho_i
             P     = v (i,j_c,k,1) * T (i,j_c,k) * W_i (i,j_c,k)
             gam_c = thd % gam2 * W_i (i,j_c,k) / cp (i,j_c,k)
             gam_c = 1.0_dp / ( 1.0_dp - gam_c )
             cs_c  = sqrt ( gam_c * P * rho_i )

             ! presolving some variables (2/2)
             vy_c = v (i,j_c,k,3) * rho_i
             wz_c = v (i,j_c,k,4) * rho_i
             do l = 1 , nrv+npv+nvv
                Y_c (l) = v (i,j_c,k,niv+l) * rho_i
             end do
             ht_c = 0.0_dp
             wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T (i,j_c,k)
             do l = 1 , nrv
                ht_c = ht_c + ha (i,j_c,k,l) * Y_c (l)
                psi_c (l) = - ha (i,j_c,k,l)
                psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             end do


             cs_c_i = 1.0_dp / cs_c
             cs_c_2 = cs_c * cs_c
             q      = 0.5_dp * ( ux_c*ux_c + vy_c*vy_c + wz_c*wz_c )
             ht_c   = ht_c + q
             b2     = ( gam_c - 1.0_dp ) / cs_c_2
             b1     = b2 * q
             xi     = b1 - b2 * ht_c


             ! write (*,*) ux_c , vy_c , wz_c
             ! write (*,*)
             ! write (*,*) T_c , gam_c , cs_c  , ht_c , cs_c_i , cs_c_2 , q , b1 , b2 , xi
             ! write (*,*)
             ! write (*,*) Y_c (:)
             ! write (*,*)
             ! write (*,*) ha (j_l,j,k,:)
             ! write (*,*)
             ! write (*,*) psi_c
             ! write (*,*)
             ! write (*,*) f_s (1,:)
             ! write (*,*)
             ! read(*,*)


             ! left eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             el(1,1)   =   0.5_dp * ( b1        + vy_c * cs_c_i )
             el(1,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(1,3)   = - 0.5_dp * ( b2 * vy_c +        cs_c_i )
             el(1,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(1,5)   =   0.5_dp * b2

             el(2,1)   =   1.0_dp - b1
             el(2,2)   =   b2 * ux_c
             el(2,3)   =   b2 * vy_c
             el(2,4)   =   b2 * wz_c
             el(2,5)   = - b2

             el(3,1)   =   0.5_dp * ( b1        - vy_c * cs_c_i )
             el(3,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(3,3)   = - 0.5_dp * ( b2 * vy_c -        cs_c_i )
             el(3,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(3,5)   =   0.5_dp * b2

             el(4,1)   = - ux_c
             el(4,2)   =   1.0_dp
             el(4,3)   =   0.0_dp
             el(4,4)   =   0.0_dp
             el(4,5)   =   0.0_dp

             el(5,1)   = - wz_c
             el(5,2)   =   0.0_dp
             el(5,3)   =   0.0_dp
             el(5,4)   =   1.0_dp
             el(5,5)   =   0.0_dp

             ! 2 : SECOND
             wrk = 0.5_dp * b2
             do l = niv+1 , nv-npv-nvv
                el ( 1 , l ) = wrk * psi_c ( l-niv )
             end do
             el ( 2   , niv+1:nv-npv-nvv ) = - el ( 1 , niv+1:nv-npv-nvv ) - el ( 1 , niv+1:nv-npv-nvv )
             el ( 3   , niv+1:nv-npv-nvv ) =   el ( 1 , niv+1:nv-npv-nvv )
             el ( 4:5 , niv+1:nv-npv-nvv ) =   0.0_dp
             if ( npv > 0 ) then
                el ( 1:5 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             ! 3 : THIRD
             el(6,1)   =   ( 1 + xi ) * q
             el(6,2)   = - ( 1 + xi ) * ux_c
             el(6,3)   = - ( 1 + xi ) * vy_c
             el(6,4)   = - ( 1 + xi ) * wz_c
             el(6,5)   =   ( 1 + xi )
             do l = niv+1 , nv-npv-nvv
                el ( niv+1 , l ) = xi * psi_c ( l-niv )
             end do
             if ( npv > 0 ) then
                el ( niv+1 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             if ( nv-npv-nvv >= 7 ) then
                ! 4 : FOURTH
                do l = 7 , nv-npv-nvv
                   el(l,1) = - b1        * Y_c (l-6)
                   el(l,2) =   b2 * ux_c * Y_c (l-6)
                   el(l,3) =   b2 * vy_c * Y_c (l-6)
                   el(l,4) =   b2 * wz_c * Y_c (l-6)
                   el(l,5) = - b2        * Y_c (l-6)
                end do
                if ( npv > 0 ) then
                   do l = nv-npv-nvv+1 , nv
                      el(l,1) = - b1        * Y_c (l-5)
                      el(l,2) =   b2 * ux_c * Y_c (l-5)
                      el(l,3) =   b2 * vy_c * Y_c (l-5)
                      el(l,4) =   b2 * wz_c * Y_c (l-5)
                      el(l,5) = - b2        * Y_c (l-5)
                   end do
                end if

                ! 5 : FIFTH
                do s1 = 1 , nv-npv-nvv-niv-1
                   do s2 = 1 , nv-npv-nvv-niv
                      el ( s1+6 , s2+5 ) = - b2 * Y_c (s1) * psi_c (s2)
                   end do
                end do
                do l = 1 , nv-npv-nvv-niv-1
                   el ( l+6 , l+5 ) = el ( l+6 , l+5 ) + 1.0_dp
                end do
                if ( npv > 0 ) then
                   el (7:nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , nv-npv-nvv-niv
                         el ( s1+nv-npv-nvv , s2+5 ) = - b2 * Y_c (s1+nv-npv-nvv-niv) * psi_c (s2)
                      end do
                   end do
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , npv+nvv
                         el ( s1+nv-npv-nvv , s2+nv-npv-nvv ) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      el ( l+nv-npv-nvv , l+nv-npv-nvv ) = 1.0_dp
                   end do
                end if
             end if


             ! right eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             er(1,1)   =  1.0_dp
             er(1,2)   =  1.0_dp
             er(1,3)   =  1.0_dp
             er(1,4)   =  0.0_dp
             er(1,5)   =  0.0_dp

             er(2,1)   =  ux_c
             er(2,2)   =  ux_c
             er(2,3)   =  ux_c
             er(2,4)   =  1.0_dp
             er(2,5)   =  0.0_dp

             er(3,1)   =  vy_c - cs_c
             er(3,2)   =  vy_c
             er(3,3)   =  vy_c + cs_c
             er(3,4)   =  0.0_dp
             er(3,5)   =  0.0_dp

             er(4,1)   =  wz_c
             er(4,2)   =  wz_c
             er(4,3)   =  wz_c
             er(4,4)   =  0.0_dp
             er(4,5)   =  1.0_dp

             er(5,1)   =  ht_c - vy_c * cs_c
             er(5,2)   =  q
             er(5,3)   =  ht_c + vy_c * cs_c
             er(5,4)   =  ux_c
             er(5,5)   =  wz_c

             ! 2 : SECOND
             er(1:4,niv+1:nv) = 0.0_dp

             ! 3 : THIRD
             do l = niv+1 , nv
                er(l,1) =  Y_c (l-5)
                er(l,2) =  0.0_dp
                er(l,3) =  Y_c (l-5)
                er(l,4) =  0.0_dp
                er(l,5) =  0.0_dp
             end do

             ! 4 : FOURTH
             do s1 = niv , nv-npv-nvv
                do s2 = niv+1 , nv
                   er(s1,s2) = 0.0_dp
                end do
             end do
             do l = niv+1 , nv-npv-nvv
                er(l-1,l) = 1.0_dp
             end do

             ! 5 : FIFTH
             er ( nv-npv-nvv , 6 ) = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             if ( nv-npv-nvv >= 7 ) then
                wrk = - 1.0_dp / psi_c (nv-npv-nvv-niv)
                do l = 1 , nv-npv-nvv-niv-1
                   er (nv-npv-nvv,l+niv+1) = wrk * psi_c (l)
                end do
                if ( npv > 0 ) then
                   er (nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = nv-npv-nvv+1 , nv
                      do s2 = niv+1 , nv
                         er(s1,s2) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      er(l+nv-npv-nvv,l+nv-npv-nvv) = 1.0_dp
                   end do
                end if
             end if



!          ! TESTING MATRICES
!             if (rank==0) then
                ! matrix = matmul ( er , el )
                ! write (*,*) 'er*el'
                ! do l = 1 , nv
                !    write (*,'(50(1X,1PE15.6))') matrix ( l , : )
                ! end do
                ! stop
!            end if
         ! do j1 = 1 , nvar
         !    do i1 = 1 , nvar
         !       if ( i1==j1 ) then
         !          if ( abs (matrix ( i1 , j1 ) - 1.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       else
         !          if ( abs (matrix ( i1 , j1 ) - 0.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       end if
         !    end do
         ! end do
!          eigenvalmatrix = 0.
!          do s1 = 1 , nvar
!             eigenvalmatrix ( s1 , s1 ) = uu
!          end do
!          eigenvalmatrix ( 1 , 1 ) = eigenvalmatrix ( 1 , 1 ) - c
!          eigenvalmatrix ( 3 , 3 ) = eigenvalmatrix ( 3 , 3 ) + c
! !
!          write (*,*) 'eigenvalues'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') eigenvalmatrix ( i1 , : )
!          end do
! !
! !
!          jacobtest = 0.
!          jacobtest = matmul ( eigenvalmatrix , el )
!          el        = matmul ( er , jacobtest )
!          jacobtest = el
!          write (*,*) 'calculated jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacobtest ( i1 , : )
!          end do
! !
! !
!          jacob = 0.
!          jacob ( 1,1 ) = 0.
!          jacob ( 1,2 ) = 1.
!          jacob ( 1 , 3:nvar ) = 0.
! !
!          jacob ( 2 , 1 ) = - ( ( 1. - gmm ) * qq + uu*uu)
!          jacob ( 2 , 2 ) =   ( 3. - gmm ) * uu
!          jacob ( 2 , 3 ) =   ( 1. - gmm ) * vv
!          jacob ( 2 , 4 ) =   ( 1. - gmm ) * ww
!          jacob ( 2 , 5 ) = - ( 1. - gmm )
!          do s1 = 1 , nvar-5
!             jacob ( 2 , s1+5 ) = - ( 1. - gmm ) * psim ( s1 )
!          end do
! !
!          jacob ( 3 , 1 ) = - uu * vv
!          jacob ( 3 , 2 ) = vv
!          jacob ( 3 , 3 ) = uu
!          jacob ( 3 , 4:nvar ) = 0.
! !
!          jacob ( 4 , 1 ) = - uu * ww
!          jacob ( 4 , 2 ) = ww
!          jacob ( 4 , 3 ) = 0.
!          jacob ( 4 , 4 ) = uu
!          jacob ( 4 , 5:nvar ) = 0.
! !
!          jacob ( 5 , 1 ) = - ( ( 1. - gmm ) * qq + h ) * uu
!          jacob ( 5 , 2 ) =   ( 1. - gmm ) * uu * uu + h
!          jacob ( 5 , 3 ) =   ( 1. - gmm ) * uu * vv
!          jacob ( 5 , 4 ) =   ( 1. - gmm ) * uu * ww
!          jacob ( 5 , 5 ) = gmm * uu
!          do s1 = 1 , nvar-5
!             jacob ( 5 , s1+5 ) = - ( 1. - gmm ) * uu * psim ( s1 )
!          end do
! !
!          do s1 = 1 , nvar-5
!             jacob ( 5 + s1 , 1 ) = - uu * chim ( s1 )
!             jacob ( 5 + s1 , 2 ) =        chim ( s1 )
!             jacob ( 5 + s1 , 3:5 ) = 0.
!          end do
! !
!          do s1 = 6 , nvar
!             jacob ( s1 , s1 ) = uu
!             ! the rest of the elements are set to zero in the
!             ! initialisation of the jacobian
!          end do
! !
!          write (*,*) 'analytical jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacob ( i1 , : )
!          end do
!          write (*,*) 'End of matrices'
!          read (*,*)


             ! Eigenvalues (killing the positive ones)
             L_s (1) =    min ( vy_c - cs_c , 0.0_dp )
             L_s (2) =    min ( vy_c        , 0.0_dp )
             L_s (3) =    min ( vy_c + cs_c , 0.0_dp )
             L_s (4) =    min ( vy_c        , 0.0_dp )
             L_s (5) =    min ( vy_c        , 0.0_dp )
             do m = niv+1 , nv
                L_s (m) = min ( vy_c        , 0.0_dp )
             end do


             ! Derivatives of conservative variables
             do m = 1 , nv
                dw (m) = 0.0_dp
                do st = 1 , stencil_m2
                   j_s = j_c + st - 1 ! be careful
                   dw (m) = dw (m) + cnr (st) * v (i,j_s,k,m) ! be careful
                end do
             end do


             ! Projection on the eigendirections
             do m = 1 , nv
                dwc (m) = 0.0_dp
                do mm = 1 , nv
                   dwc (m) = dwc (m) + el (m,mm) * dw (mm)
                end do
             end do


             ! Returning to conservative variables
             do m = 1 , nv
                df = 0.0_dp
                do mm = 1 , nv
                   df = df + er (m,mm) * L_s (mm) * dwc (mm)
                end do
                fl (i,j_c,k,m) = fl (i,j_c,k,m) + df * dy_i (j_c)
             end do


          end do ! end of j-loop
       end do ! end of k-loop


    else if ( face == N ) then


       j_c = ey


       do k = sz , ez ! loop in the z-direction
          do i = sx , ex ! loop in the x-direction


             ! presolving some variables (1/2)
             rho_i = 1.0_dp / v (i,j_c,k,1)
             ux_c  = v (i,j_c,k,2) * rho_i
             P     = v (i,j_c,k,1) * T (i,j_c,k) * W_i (i,j_c,k)
             gam_c = thd % gam2 * W_i (i,j_c,k) / cp (i,j_c,k)
             gam_c = 1.0_dp / ( 1.0_dp - gam_c )
             cs_c  = sqrt ( gam_c * P * rho_i )

             ! presolving some variables (2/2)
             vy_c = v (i,j_c,k,3) * rho_i
             wz_c = v (i,j_c,k,4) * rho_i
             do l = 1 , nrv+npv+nvv
                Y_c (l) = v (i,j_c,k,niv+l) * rho_i
             end do
             ht_c = 0.0_dp
             wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T (i,j_c,k)
             do l = 1 , nrv
                ht_c = ht_c + ha (i,j_c,k,l) * Y_c (l)
                psi_c (l) = - ha (i,j_c,k,l)
                psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             end do


             cs_c_i = 1.0_dp / cs_c
             cs_c_2 = cs_c * cs_c
             q      = 0.5_dp * ( ux_c*ux_c + vy_c*vy_c + wz_c*wz_c )
             ht_c   = ht_c + q
             b2     = ( gam_c - 1.0_dp ) / cs_c_2
             b1     = b2 * q
             xi     = b1 - b2 * ht_c


             ! write (*,*) ux_c , vy_c , wz_c
             ! write (*,*)
             ! write (*,*) T_c , gam_c , cs_c  , ht_c , cs_c_i , cs_c_2 , q , b1 , b2 , xi
             ! write (*,*)
             ! write (*,*) Y_c (:)
             ! write (*,*)
             ! write (*,*) ha (j_l,j,k,:)
             ! write (*,*)
             ! write (*,*) psi_c
             ! write (*,*)
             ! write (*,*) f_s (1,:)
             ! write (*,*)
             ! read(*,*)


             ! left eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             el(1,1)   =   0.5_dp * ( b1        + vy_c * cs_c_i )
             el(1,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(1,3)   = - 0.5_dp * ( b2 * vy_c +        cs_c_i )
             el(1,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(1,5)   =   0.5_dp * b2

             el(2,1)   =   1.0_dp - b1
             el(2,2)   =   b2 * ux_c
             el(2,3)   =   b2 * vy_c
             el(2,4)   =   b2 * wz_c
             el(2,5)   = - b2

             el(3,1)   =   0.5_dp * ( b1        - vy_c * cs_c_i )
             el(3,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(3,3)   = - 0.5_dp * ( b2 * vy_c -        cs_c_i )
             el(3,4)   = - 0.5_dp * ( b2 * wz_c                 )
             el(3,5)   =   0.5_dp * b2

             el(4,1)   = - ux_c
             el(4,2)   =   1.0_dp
             el(4,3)   =   0.0_dp
             el(4,4)   =   0.0_dp
             el(4,5)   =   0.0_dp

             el(5,1)   = - wz_c
             el(5,2)   =   0.0_dp
             el(5,3)   =   0.0_dp
             el(5,4)   =   1.0_dp
             el(5,5)   =   0.0_dp

             ! 2 : SECOND
             wrk = 0.5_dp * b2
             do l = niv+1 , nv-npv-nvv
                el ( 1 , l ) = wrk * psi_c ( l-niv )
             end do
             el ( 2   , niv+1:nv-npv-nvv ) = - el ( 1 , niv+1:nv-npv-nvv ) - el ( 1 , niv+1:nv-npv-nvv )
             el ( 3   , niv+1:nv-npv-nvv ) =   el ( 1 , niv+1:nv-npv-nvv )
             el ( 4:5 , niv+1:nv-npv-nvv ) =   0.0_dp
             if ( npv > 0 ) then
                el ( 1:5 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             ! 3 : THIRD
             el(6,1)   =   ( 1 + xi ) * q
             el(6,2)   = - ( 1 + xi ) * ux_c
             el(6,3)   = - ( 1 + xi ) * vy_c
             el(6,4)   = - ( 1 + xi ) * wz_c
             el(6,5)   =   ( 1 + xi )
             do l = niv+1 , nv-npv-nvv
                el ( niv+1 , l ) = xi * psi_c ( l-niv )
             end do
             if ( npv > 0 ) then
                el ( niv+1 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             if ( nv-npv-nvv >= 7 ) then
                ! 4 : FOURTH
                do l = 7 , nv-npv-nvv
                   el(l,1) = - b1        * Y_c (l-6)
                   el(l,2) =   b2 * ux_c * Y_c (l-6)
                   el(l,3) =   b2 * vy_c * Y_c (l-6)
                   el(l,4) =   b2 * wz_c * Y_c (l-6)
                   el(l,5) = - b2        * Y_c (l-6)
                end do
                if ( npv > 0 ) then
                   do l = nv-npv-nvv+1 , nv
                      el(l,1) = - b1        * Y_c (l-5)
                      el(l,2) =   b2 * ux_c * Y_c (l-5)
                      el(l,3) =   b2 * vy_c * Y_c (l-5)
                      el(l,4) =   b2 * wz_c * Y_c (l-5)
                      el(l,5) = - b2        * Y_c (l-5)
                   end do
                end if

                ! 5 : FIFTH
                do s1 = 1 , nv-npv-nvv-niv-1
                   do s2 = 1 , nv-npv-nvv-niv
                      el ( s1+6 , s2+5 ) = - b2 * Y_c (s1) * psi_c (s2)
                   end do
                end do
                do l = 1 , nv-npv-nvv-niv-1
                   el ( l+6 , l+5 ) = el ( l+6 , l+5 ) + 1.0_dp
                end do
                if ( npv > 0 ) then
                   el (7:nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , nv-npv-nvv-niv
                         el ( s1+nv-npv-nvv , s2+5 ) = - b2 * Y_c (s1+nv-npv-nvv-niv) * psi_c (s2)
                      end do
                   end do
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , npv+nvv
                         el ( s1+nv-npv-nvv , s2+nv-npv-nvv ) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      el ( l+nv-npv-nvv , l+nv-npv-nvv ) = 1.0_dp
                   end do
                end if
             end if


             ! right eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             er(1,1)   =  1.0_dp
             er(1,2)   =  1.0_dp
             er(1,3)   =  1.0_dp
             er(1,4)   =  0.0_dp
             er(1,5)   =  0.0_dp

             er(2,1)   =  ux_c
             er(2,2)   =  ux_c
             er(2,3)   =  ux_c
             er(2,4)   =  1.0_dp
             er(2,5)   =  0.0_dp

             er(3,1)   =  vy_c - cs_c
             er(3,2)   =  vy_c
             er(3,3)   =  vy_c + cs_c
             er(3,4)   =  0.0_dp
             er(3,5)   =  0.0_dp

             er(4,1)   =  wz_c
             er(4,2)   =  wz_c
             er(4,3)   =  wz_c
             er(4,4)   =  0.0_dp
             er(4,5)   =  1.0_dp

             er(5,1)   =  ht_c - vy_c * cs_c
             er(5,2)   =  q
             er(5,3)   =  ht_c + vy_c * cs_c
             er(5,4)   =  ux_c
             er(5,5)   =  wz_c

             ! 2 : SECOND
             er(1:4,niv+1:nv) = 0.0_dp

             ! 3 : THIRD
             do l = niv+1 , nv
                er(l,1) =  Y_c (l-5)
                er(l,2) =  0.0_dp
                er(l,3) =  Y_c (l-5)
                er(l,4) =  0.0_dp
                er(l,5) =  0.0_dp
             end do

             ! 4 : FOURTH
             do s1 = niv , nv-npv-nvv
                do s2 = niv+1 , nv
                   er(s1,s2) = 0.0_dp
                end do
             end do
             do l = niv+1 , nv-npv-nvv
                er(l-1,l) = 1.0_dp
             end do

             ! 5 : FIFTH
             er ( nv-npv-nvv , 6 ) = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             if ( nv-npv-nvv >= 7 ) then
                wrk = - 1.0_dp / psi_c (nv-npv-nvv-niv)
                do l = 1 , nv-npv-nvv-niv-1
                   er (nv-npv-nvv,l+niv+1) = wrk * psi_c (l)
                end do
                if ( npv > 0 ) then
                   er (nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = nv-npv-nvv+1 , nv
                      do s2 = niv+1 , nv
                         er(s1,s2) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      er(l+nv-npv-nvv,l+nv-npv-nvv) = 1.0_dp
                   end do
                end if
             end if



!          ! TESTING MATRICES
!             if (rank==0) then
                ! matrix = matmul ( er , el )
                ! write (*,*) 'er*el'
                ! do l = 1 , nv
                !    write (*,'(50(1X,1PE15.6))') matrix ( l , : )
                ! end do
                ! stop
!            end if
         ! do j1 = 1 , nvar
         !    do i1 = 1 , nvar
         !       if ( i1==j1 ) then
         !          if ( abs (matrix ( i1 , j1 ) - 1.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       else
         !          if ( abs (matrix ( i1 , j1 ) - 0.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       end if
         !    end do
         ! end do
!          eigenvalmatrix = 0.
!          do s1 = 1 , nvar
!             eigenvalmatrix ( s1 , s1 ) = uu
!          end do
!          eigenvalmatrix ( 1 , 1 ) = eigenvalmatrix ( 1 , 1 ) - c
!          eigenvalmatrix ( 3 , 3 ) = eigenvalmatrix ( 3 , 3 ) + c
! !
!          write (*,*) 'eigenvalues'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') eigenvalmatrix ( i1 , : )
!          end do
! !
! !
!          jacobtest = 0.
!          jacobtest = matmul ( eigenvalmatrix , el )
!          el        = matmul ( er , jacobtest )
!          jacobtest = el
!          write (*,*) 'calculated jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacobtest ( i1 , : )
!          end do
! !
! !
!          jacob = 0.
!          jacob ( 1,1 ) = 0.
!          jacob ( 1,2 ) = 1.
!          jacob ( 1 , 3:nvar ) = 0.
! !
!          jacob ( 2 , 1 ) = - ( ( 1. - gmm ) * qq + uu*uu)
!          jacob ( 2 , 2 ) =   ( 3. - gmm ) * uu
!          jacob ( 2 , 3 ) =   ( 1. - gmm ) * vv
!          jacob ( 2 , 4 ) =   ( 1. - gmm ) * ww
!          jacob ( 2 , 5 ) = - ( 1. - gmm )
!          do s1 = 1 , nvar-5
!             jacob ( 2 , s1+5 ) = - ( 1. - gmm ) * psim ( s1 )
!          end do
! !
!          jacob ( 3 , 1 ) = - uu * vv
!          jacob ( 3 , 2 ) = vv
!          jacob ( 3 , 3 ) = uu
!          jacob ( 3 , 4:nvar ) = 0.
! !
!          jacob ( 4 , 1 ) = - uu * ww
!          jacob ( 4 , 2 ) = ww
!          jacob ( 4 , 3 ) = 0.
!          jacob ( 4 , 4 ) = uu
!          jacob ( 4 , 5:nvar ) = 0.
! !
!          jacob ( 5 , 1 ) = - ( ( 1. - gmm ) * qq + h ) * uu
!          jacob ( 5 , 2 ) =   ( 1. - gmm ) * uu * uu + h
!          jacob ( 5 , 3 ) =   ( 1. - gmm ) * uu * vv
!          jacob ( 5 , 4 ) =   ( 1. - gmm ) * uu * ww
!          jacob ( 5 , 5 ) = gmm * uu
!          do s1 = 1 , nvar-5
!             jacob ( 5 , s1+5 ) = - ( 1. - gmm ) * uu * psim ( s1 )
!          end do
! !
!          do s1 = 1 , nvar-5
!             jacob ( 5 + s1 , 1 ) = - uu * chim ( s1 )
!             jacob ( 5 + s1 , 2 ) =        chim ( s1 )
!             jacob ( 5 + s1 , 3:5 ) = 0.
!          end do
! !
!          do s1 = 6 , nvar
!             jacob ( s1 , s1 ) = uu
!             ! the rest of the elements are set to zero in the
!             ! initialisation of the jacobian
!          end do
! !
!          write (*,*) 'analytical jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacob ( i1 , : )
!          end do
!          write (*,*) 'End of matrices'
!          read (*,*)


             ! Eigenvalues (killing the positive ones)
             L_s (1) =    max ( vy_c - cs_c , 0.0_dp )
             L_s (2) =    max ( vy_c        , 0.0_dp )
             L_s (3) =    max ( vy_c + cs_c , 0.0_dp )
             L_s (4) =    max ( vy_c        , 0.0_dp )
             L_s (5) =    max ( vy_c        , 0.0_dp )
             do m = niv+1 , nv
                L_s (m) = max ( vy_c        , 0.0_dp )
             end do


             ! Derivatives of conservative variables
             do m = 1 , nv
                dw (m) = 0.0_dp
                do st = 1 , stencil_m2
                   j_s = j_c - st + 1 ! be careful
                   dw (m) = dw (m) - cnr (st) * v (i,j_s,k,m) ! be careful
                end do
             end do


             ! Projection on the eigendirections
             do m = 1 , nv
                dwc (m) = 0.0_dp
                do mm = 1 , nv
                   dwc (m) = dwc (m) + el (m,mm) * dw (mm)
                end do
             end do


             ! Returning to conservative variables
             do m = 1 , nv
                df = 0.0_dp
                do mm = 1 , nv
                   df = df + er (m,mm) * L_s (mm) * dwc (mm)
                end do
                fl (i,j_c,k,m) = fl (i,j_c,k,m) + df * dy_i (j_c)
             end do


          end do ! end of j-loop
       end do ! end of k-loop


    else if ( face == B ) then


       k_c = sz


       do j = sy , ey ! loop in the y-direction
          do i = sx , ex ! loop in the x-direction


             ! presolving some variables (1/2)
             rho_i = 1.0_dp / v (i,j,k_c,1)
             ux_c  = v (i,j,k_c,2) * rho_i
             P     = v (i,j,k_c,1) * T (i,j,k_c) * W_i (i,j,k_c)
             gam_c = thd % gam2 * W_i (i,j,k_c) / cp (i,j,k_c)
             gam_c = 1.0_dp / ( 1.0_dp - gam_c )
             cs_c  = sqrt ( gam_c * P * rho_i )

             ! presolving some variables (2/2)
             vy_c = v (i,j,k_c,3) * rho_i
             wz_c = v (i,j,k_c,4) * rho_i
             do l = 1 , nrv+npv+nvv
                Y_c (l) = v (i,j,k_c,niv+l) * rho_i
             end do
             ht_c = 0.0_dp
             wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T (i,j,k_c)
             do l = 1 , nrv
                ht_c = ht_c + ha (i,j,k_c,l) * Y_c (l)
                psi_c (l) = - ha (i,j,k_c,l)
                psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             end do


             cs_c_i = 1.0_dp / cs_c
             cs_c_2 = cs_c * cs_c
             q      = 0.5_dp * ( ux_c*ux_c + vy_c*vy_c + wz_c*wz_c )
             ht_c   = ht_c + q
             b2     = ( gam_c - 1.0_dp ) / cs_c_2
             b1     = b2 * q
             xi     = b1 - b2 * ht_c


             ! write (*,*) ux_c , vy_c , wz_c
             ! write (*,*)
             ! write (*,*) T_c , gam_c , cs_c  , ht_c , cs_c_i , cs_c_2 , q , b1 , b2 , xi
             ! write (*,*)
             ! write (*,*) Y_c (:)
             ! write (*,*)
             ! write (*,*) ha (j_l,j,k,:)
             ! write (*,*)
             ! write (*,*) psi_c
             ! write (*,*)
             ! write (*,*) f_s (1,:)
             ! write (*,*)
             ! read(*,*)


             ! left eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             el(1,1)   =   0.5_dp * ( b1        + wz_c * cs_c_i )
             el(1,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(1,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(1,4)   = - 0.5_dp * ( b2 * wz_c +        cs_c_i )
             el(1,5)   =   0.5_dp * b2

             el(2,1)   =   1.0_dp - b1
             el(2,2)   =   b2 * ux_c
             el(2,3)   =   b2 * vy_c
             el(2,4)   =   b2 * wz_c
             el(2,5)   = - b2

             el(3,1)   =   0.5_dp * ( b1        - wz_c * cs_c_i )
             el(3,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(3,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(3,4)   = - 0.5_dp * ( b2 * wz_c -        cs_c_i )
             el(3,5)   =   0.5_dp * b2

             el(4,1)   = - ux_c
             el(4,2)   =   1.0_dp
             el(4,3)   =   0.0_dp
             el(4,4)   =   0.0_dp
             el(4,5)   =   0.0_dp

             el(5,1)   = - vy_c
             el(5,2)   =   0.0_dp
             el(5,3)   =   1.0_dp
             el(5,4)   =   0.0_dp
             el(5,5)   =   0.0_dp

             ! 2 : SECOND
             wrk = 0.5_dp * b2
             do l = niv+1 , nv-npv-nvv
                el ( 1 , l ) = wrk * psi_c ( l-niv )
             end do
             el ( 2   , niv+1:nv-npv-nvv ) = - el ( 1 , niv+1:nv-npv-nvv ) - el ( 1 , niv+1:nv-npv-nvv )
             el ( 3   , niv+1:nv-npv-nvv ) =   el ( 1 , niv+1:nv-npv-nvv )
             el ( 4:5 , niv+1:nv-npv-nvv ) =   0.0_dp
             if ( npv > 0 ) then
                el ( 1:5 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             ! 3 : THIRD
             el(6,1)   =   ( 1 + xi ) * q
             el(6,2)   = - ( 1 + xi ) * ux_c
             el(6,3)   = - ( 1 + xi ) * vy_c
             el(6,4)   = - ( 1 + xi ) * wz_c
             el(6,5)   =   ( 1 + xi )
             do l = niv+1 , nv-npv-nvv
                el ( niv+1 , l ) = xi * psi_c ( l-niv )
             end do
             if ( npv > 0 ) then
                el ( niv+1 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             if ( nv-npv-nvv >= 7 ) then
                ! 4 : FOURTH
                do l = 7 , nv-npv-nvv
                   el(l,1) = - b1        * Y_c (l-6)
                   el(l,2) =   b2 * ux_c * Y_c (l-6)
                   el(l,3) =   b2 * vy_c * Y_c (l-6)
                   el(l,4) =   b2 * wz_c * Y_c (l-6)
                   el(l,5) = - b2        * Y_c (l-6)
                end do
                if ( npv > 0 ) then
                   do l = nv-npv-nvv+1 , nv
                      el(l,1) = - b1        * Y_c (l-5)
                      el(l,2) =   b2 * ux_c * Y_c (l-5)
                      el(l,3) =   b2 * vy_c * Y_c (l-5)
                      el(l,4) =   b2 * wz_c * Y_c (l-5)
                      el(l,5) = - b2        * Y_c (l-5)
                   end do
                end if

                ! 5 : FIFTH
                do s1 = 1 , nv-npv-nvv-niv-1
                   do s2 = 1 , nv-npv-nvv-niv
                      el ( s1+6 , s2+5 ) = - b2 * Y_c (s1) * psi_c (s2)
                   end do
                end do
                do l = 1 , nv-npv-nvv-niv-1
                   el ( l+6 , l+5 ) = el ( l+6 , l+5 ) + 1.0_dp
                end do
                if ( npv > 0 ) then
                   el (7:nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , nv-npv-nvv-niv
                         el ( s1+nv-npv-nvv , s2+5 ) = - b2 * Y_c (s1+nv-npv-nvv-niv) * psi_c (s2)
                      end do
                   end do
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , npv+nvv
                         el ( s1+nv-npv-nvv , s2+nv-npv-nvv ) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      el ( l+nv-npv-nvv , l+nv-npv-nvv ) = 1.0_dp
                   end do
                end if
             end if


             ! right eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             er(1,1)   =  1.0_dp
             er(1,2)   =  1.0_dp
             er(1,3)   =  1.0_dp
             er(1,4)   =  0.0_dp
             er(1,5)   =  0.0_dp

             er(2,1)   =  ux_c
             er(2,2)   =  ux_c
             er(2,3)   =  ux_c
             er(2,4)   =  1.0_dp
             er(2,5)   =  0.0_dp

             er(3,1)   =  vy_c
             er(3,2)   =  vy_c
             er(3,3)   =  vy_c
             er(3,4)   =  0.0_dp
             er(3,5)   =  1.0_dp

             er(4,1)   =  wz_c - cs_c
             er(4,2)   =  wz_c
             er(4,3)   =  wz_c + cs_c
             er(4,4)   =  0.0_dp
             er(4,5)   =  0.0_dp

             er(5,1)   =  ht_c - wz_c * cs_c
             er(5,2)   =  q
             er(5,3)   =  ht_c + wz_c * cs_c
             er(5,4)   =  ux_c
             er(5,5)   =  vy_c

             ! 2 : SECOND
             er(1:4,niv+1:nv) = 0.0_dp

             ! 3 : THIRD
             do l = niv+1 , nv
                er(l,1) =  Y_c (l-5)
                er(l,2) =  0.0_dp
                er(l,3) =  Y_c (l-5)
                er(l,4) =  0.0_dp
                er(l,5) =  0.0_dp
             end do

             ! 4 : FOURTH
             do s1 = niv , nv-npv-nvv
                do s2 = niv+1 , nv
                   er(s1,s2) = 0.0_dp
                end do
             end do
             do l = niv+1 , nv-npv-nvv
                er(l-1,l) = 1.0_dp
             end do

             ! 5 : FIFTH
             er ( nv-npv-nvv , 6 ) = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             if ( nv-npv-nvv >= 7 ) then
                wrk = - 1.0_dp / psi_c (nv-npv-nvv-niv)
                do l = 1 , nv-npv-nvv-niv-1
                   er (nv-npv-nvv,l+niv+1) = wrk * psi_c (l)
                end do
                if ( npv > 0 ) then
                   er (nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = nv-npv-nvv+1 , nv
                      do s2 = niv+1 , nv
                         er(s1,s2) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      er(l+nv-npv-nvv,l+nv-npv-nvv) = 1.0_dp
                   end do
                end if
             end if



!          ! TESTING MATRICES
!             if (rank==0) then
                ! matrix = matmul ( er , el )
                ! write (*,*) 'er*el'
                ! do l = 1 , nv
                !    write (*,'(50(1X,1PE15.6))') matrix ( l , : )
                ! end do
                ! stop
!            end if
         ! do j1 = 1 , nvar
         !    do i1 = 1 , nvar
         !       if ( i1==j1 ) then
         !          if ( abs (matrix ( i1 , j1 ) - 1.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       else
         !          if ( abs (matrix ( i1 , j1 ) - 0.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       end if
         !    end do
         ! end do
!          eigenvalmatrix = 0.
!          do s1 = 1 , nvar
!             eigenvalmatrix ( s1 , s1 ) = uu
!          end do
!          eigenvalmatrix ( 1 , 1 ) = eigenvalmatrix ( 1 , 1 ) - c
!          eigenvalmatrix ( 3 , 3 ) = eigenvalmatrix ( 3 , 3 ) + c
! !
!          write (*,*) 'eigenvalues'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') eigenvalmatrix ( i1 , : )
!          end do
! !
! !
!          jacobtest = 0.
!          jacobtest = matmul ( eigenvalmatrix , el )
!          el        = matmul ( er , jacobtest )
!          jacobtest = el
!          write (*,*) 'calculated jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacobtest ( i1 , : )
!          end do
! !
! !
!          jacob = 0.
!          jacob ( 1,1 ) = 0.
!          jacob ( 1,2 ) = 1.
!          jacob ( 1 , 3:nvar ) = 0.
! !
!          jacob ( 2 , 1 ) = - ( ( 1. - gmm ) * qq + uu*uu)
!          jacob ( 2 , 2 ) =   ( 3. - gmm ) * uu
!          jacob ( 2 , 3 ) =   ( 1. - gmm ) * vv
!          jacob ( 2 , 4 ) =   ( 1. - gmm ) * ww
!          jacob ( 2 , 5 ) = - ( 1. - gmm )
!          do s1 = 1 , nvar-5
!             jacob ( 2 , s1+5 ) = - ( 1. - gmm ) * psim ( s1 )
!          end do
! !
!          jacob ( 3 , 1 ) = - uu * vv
!          jacob ( 3 , 2 ) = vv
!          jacob ( 3 , 3 ) = uu
!          jacob ( 3 , 4:nvar ) = 0.
! !
!          jacob ( 4 , 1 ) = - uu * ww
!          jacob ( 4 , 2 ) = ww
!          jacob ( 4 , 3 ) = 0.
!          jacob ( 4 , 4 ) = uu
!          jacob ( 4 , 5:nvar ) = 0.
! !
!          jacob ( 5 , 1 ) = - ( ( 1. - gmm ) * qq + h ) * uu
!          jacob ( 5 , 2 ) =   ( 1. - gmm ) * uu * uu + h
!          jacob ( 5 , 3 ) =   ( 1. - gmm ) * uu * vv
!          jacob ( 5 , 4 ) =   ( 1. - gmm ) * uu * ww
!          jacob ( 5 , 5 ) = gmm * uu
!          do s1 = 1 , nvar-5
!             jacob ( 5 , s1+5 ) = - ( 1. - gmm ) * uu * psim ( s1 )
!          end do
! !
!          do s1 = 1 , nvar-5
!             jacob ( 5 + s1 , 1 ) = - uu * chim ( s1 )
!             jacob ( 5 + s1 , 2 ) =        chim ( s1 )
!             jacob ( 5 + s1 , 3:5 ) = 0.
!          end do
! !
!          do s1 = 6 , nvar
!             jacob ( s1 , s1 ) = uu
!             ! the rest of the elements are set to zero in the
!             ! initialisation of the jacobian
!          end do
! !
!          write (*,*) 'analytical jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacob ( i1 , : )
!          end do
!          write (*,*) 'End of matrices'
!          read (*,*)


             ! Eigenvalues (killing the positive ones)
             L_s (1) =    min ( wz_c - cs_c , 0.0_dp )
             L_s (2) =    min ( wz_c        , 0.0_dp )
             L_s (3) =    min ( wz_c + cs_c , 0.0_dp )
             L_s (4) =    min ( wz_c        , 0.0_dp )
             L_s (5) =    min ( wz_c        , 0.0_dp )
             do m = niv+1 , nv
                L_s (m) = min ( wz_c        , 0.0_dp )
             end do


             ! Derivatives of conservative variables
             do m = 1 , nv
                dw (m) = 0.0_dp
                do st = 1 , stencil_m2
                   k_s = k_c + st - 1 ! be careful
                   dw (m) = dw (m) + cnr (st) * v (i,j,k_s,m) ! be careful
                end do
             end do


             ! Projection on the eigendirections
             do m = 1 , nv
                dwc (m) = 0.0_dp
                do mm = 1 , nv
                   dwc (m) = dwc (m) + el (m,mm) * dw (mm)
                end do
             end do


             ! Returning to conservative variables
             do m = 1 , nv
                df = 0.0_dp
                do mm = 1 , nv
                   df = df + er (m,mm) * L_s (mm) * dwc (mm)
                end do
                fl (i,j,k_c,m) = fl (i,j,k_c,m) + df * dz_i (k_c)
             end do


          end do ! end of j-loop
       end do ! end of k-loop


    else if ( face == F ) then


       k_c = ez


       do j = sy , ey ! loop in the y-direction
          do i = sx , ex ! loop in the x-direction


             ! presolving some variables (1/2)
             rho_i = 1.0_dp / v (i,j,k_c,1)
             ux_c  = v (i,j,k_c,2) * rho_i
             P     = v (i,j,k_c,1) * T (i,j,k_c) * W_i (i,j,k_c)
             gam_c = thd % gam2 * W_i (i,j,k_c) / cp (i,j,k_c)
             gam_c = 1.0_dp / ( 1.0_dp - gam_c )
             cs_c  = sqrt ( gam_c * P * rho_i )

             ! presolving some variables (2/2)
             vy_c = v (i,j,k_c,3) * rho_i
             wz_c = v (i,j,k_c,4) * rho_i
             do l = 1 , nrv+npv+nvv
                Y_c (l) = v (i,j,k_c,niv+l) * rho_i
             end do
             ht_c = 0.0_dp
             wrk = ( gam_c / ( gam_c - 1.0_dp ) ) * T (i,j,k_c)
             do l = 1 , nrv
                ht_c = ht_c + ha (i,j,k_c,l) * Y_c (l)
                psi_c (l) = - ha (i,j,k_c,l)
                psi_c (l) = psi_c (l) + wrk * thd % Wc_i (l)
             end do


             cs_c_i = 1.0_dp / cs_c
             cs_c_2 = cs_c * cs_c
             q      = 0.5_dp * ( ux_c*ux_c + vy_c*vy_c + wz_c*wz_c )
             ht_c   = ht_c + q
             b2     = ( gam_c - 1.0_dp ) / cs_c_2
             b1     = b2 * q
             xi     = b1 - b2 * ht_c


             ! write (*,*) ux_c , vy_c , wz_c
             ! write (*,*)
             ! write (*,*) T_c , gam_c , cs_c  , ht_c , cs_c_i , cs_c_2 , q , b1 , b2 , xi
             ! write (*,*)
             ! write (*,*) Y_c (:)
             ! write (*,*)
             ! write (*,*) ha (j_l,j,k,:)
             ! write (*,*)
             ! write (*,*) psi_c
             ! write (*,*)
             ! write (*,*) f_s (1,:)
             ! write (*,*)
             ! read(*,*)


             ! left eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             el(1,1)   =   0.5_dp * ( b1        + wz_c * cs_c_i )
             el(1,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(1,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(1,4)   = - 0.5_dp * ( b2 * wz_c +        cs_c_i )
             el(1,5)   =   0.5_dp * b2

             el(2,1)   =   1.0_dp - b1
             el(2,2)   =   b2 * ux_c
             el(2,3)   =   b2 * vy_c
             el(2,4)   =   b2 * wz_c
             el(2,5)   = - b2

             el(3,1)   =   0.5_dp * ( b1        - wz_c * cs_c_i )
             el(3,2)   = - 0.5_dp * ( b2 * ux_c                 )
             el(3,3)   = - 0.5_dp * ( b2 * vy_c                 )
             el(3,4)   = - 0.5_dp * ( b2 * wz_c -        cs_c_i )
             el(3,5)   =   0.5_dp * b2

             el(4,1)   = - ux_c
             el(4,2)   =   1.0_dp
             el(4,3)   =   0.0_dp
             el(4,4)   =   0.0_dp
             el(4,5)   =   0.0_dp

             el(5,1)   = - vy_c
             el(5,2)   =   0.0_dp
             el(5,3)   =   1.0_dp
             el(5,4)   =   0.0_dp
             el(5,5)   =   0.0_dp

             ! 2 : SECOND
             wrk = 0.5_dp * b2
             do l = niv+1 , nv-npv-nvv
                el ( 1 , l ) = wrk * psi_c ( l-niv )
             end do
             el ( 2   , niv+1:nv-npv-nvv ) = - el ( 1 , niv+1:nv-npv-nvv ) - el ( 1 , niv+1:nv-npv-nvv )
             el ( 3   , niv+1:nv-npv-nvv ) =   el ( 1 , niv+1:nv-npv-nvv )
             el ( 4:5 , niv+1:nv-npv-nvv ) =   0.0_dp
             if ( npv > 0 ) then
                el ( 1:5 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             ! 3 : THIRD
             el(6,1)   =   ( 1 + xi ) * q
             el(6,2)   = - ( 1 + xi ) * ux_c
             el(6,3)   = - ( 1 + xi ) * vy_c
             el(6,4)   = - ( 1 + xi ) * wz_c
             el(6,5)   =   ( 1 + xi )
             do l = niv+1 , nv-npv-nvv
                el ( niv+1 , l ) = xi * psi_c ( l-niv )
             end do
             if ( npv > 0 ) then
                el ( niv+1 , nv-npv-nvv+1:nv ) = 0.0_dp
             end if

             if ( nv-npv-nvv >= 7 ) then
                ! 4 : FOURTH
                do l = 7 , nv-npv-nvv
                   el(l,1) = - b1        * Y_c (l-6)
                   el(l,2) =   b2 * ux_c * Y_c (l-6)
                   el(l,3) =   b2 * vy_c * Y_c (l-6)
                   el(l,4) =   b2 * wz_c * Y_c (l-6)
                   el(l,5) = - b2        * Y_c (l-6)
                end do
                if ( npv > 0 ) then
                   do l = nv-npv-nvv+1 , nv
                      el(l,1) = - b1        * Y_c (l-5)
                      el(l,2) =   b2 * ux_c * Y_c (l-5)
                      el(l,3) =   b2 * vy_c * Y_c (l-5)
                      el(l,4) =   b2 * wz_c * Y_c (l-5)
                      el(l,5) = - b2        * Y_c (l-5)
                   end do
                end if

                ! 5 : FIFTH
                do s1 = 1 , nv-npv-nvv-niv-1
                   do s2 = 1 , nv-npv-nvv-niv
                      el ( s1+6 , s2+5 ) = - b2 * Y_c (s1) * psi_c (s2)
                   end do
                end do
                do l = 1 , nv-npv-nvv-niv-1
                   el ( l+6 , l+5 ) = el ( l+6 , l+5 ) + 1.0_dp
                end do
                if ( npv > 0 ) then
                   el (7:nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , nv-npv-nvv-niv
                         el ( s1+nv-npv-nvv , s2+5 ) = - b2 * Y_c (s1+nv-npv-nvv-niv) * psi_c (s2)
                      end do
                   end do
                   do s1 = 1 , npv+nvv
                      do s2 = 1 , npv+nvv
                         el ( s1+nv-npv-nvv , s2+nv-npv-nvv ) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      el ( l+nv-npv-nvv , l+nv-npv-nvv ) = 1.0_dp
                   end do
                end if
             end if


             ! right eigenvectors matrix (at Roe's state)


             ! 1 : FIRST
             er(1,1)   =  1.0_dp
             er(1,2)   =  1.0_dp
             er(1,3)   =  1.0_dp
             er(1,4)   =  0.0_dp
             er(1,5)   =  0.0_dp

             er(2,1)   =  ux_c
             er(2,2)   =  ux_c
             er(2,3)   =  ux_c
             er(2,4)   =  1.0_dp
             er(2,5)   =  0.0_dp

             er(3,1)   =  vy_c
             er(3,2)   =  vy_c
             er(3,3)   =  vy_c
             er(3,4)   =  0.0_dp
             er(3,5)   =  1.0_dp

             er(4,1)   =  wz_c - cs_c
             er(4,2)   =  wz_c
             er(4,3)   =  wz_c + cs_c
             er(4,4)   =  0.0_dp
             er(4,5)   =  0.0_dp

             er(5,1)   =  ht_c - wz_c * cs_c
             er(5,2)   =  q
             er(5,3)   =  ht_c + wz_c * cs_c
             er(5,4)   =  ux_c
             er(5,5)   =  vy_c

             ! 2 : SECOND
             er(1:4,niv+1:nv) = 0.0_dp

             ! 3 : THIRD
             do l = niv+1 , nv
                er(l,1) =  Y_c (l-5)
                er(l,2) =  0.0_dp
                er(l,3) =  Y_c (l-5)
                er(l,4) =  0.0_dp
                er(l,5) =  0.0_dp
             end do

             ! 4 : FOURTH
             do s1 = niv , nv-npv-nvv
                do s2 = niv+1 , nv
                   er(s1,s2) = 0.0_dp
                end do
             end do
             do l = niv+1 , nv-npv-nvv
                er(l-1,l) = 1.0_dp
             end do

             ! 5 : FIFTH
             er ( nv-npv-nvv , 6 ) = - 1.0_dp / psi_c (nv-npv-nvv-niv)
             if ( nv-npv-nvv >= 7 ) then
                wrk = - 1.0_dp / psi_c (nv-npv-nvv-niv)
                do l = 1 , nv-npv-nvv-niv-1
                   er (nv-npv-nvv,l+niv+1) = wrk * psi_c (l)
                end do
                if ( npv > 0 ) then
                   er (nv-npv-nvv,nv-npv-nvv+1:nv) = 0.0_dp
                   do s1 = nv-npv-nvv+1 , nv
                      do s2 = niv+1 , nv
                         er(s1,s2) = 0.0_dp
                      end do
                   end do
                   do l = 1 , npv+nvv
                      er(l+nv-npv-nvv,l+nv-npv-nvv) = 1.0_dp
                   end do
                end if
             end if



!          ! TESTING MATRICES
!             if (rank==0) then
                ! matrix = matmul ( er , el )
                ! write (*,*) 'er*el'
                ! do l = 1 , nv
                !    write (*,'(50(1X,1PE15.6))') matrix ( l , : )
                ! end do
                ! stop
!            end if
         ! do j1 = 1 , nvar
         !    do i1 = 1 , nvar
         !       if ( i1==j1 ) then
         !          if ( abs (matrix ( i1 , j1 ) - 1.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       else
         !          if ( abs (matrix ( i1 , j1 ) - 0.0e0 ) > 1.0e-8 ) then
         !             write (*,*) 'matmul error ' , i , i1 , j1
         !             stop
         !          end if
         !       end if
         !    end do
         ! end do
!          eigenvalmatrix = 0.
!          do s1 = 1 , nvar
!             eigenvalmatrix ( s1 , s1 ) = uu
!          end do
!          eigenvalmatrix ( 1 , 1 ) = eigenvalmatrix ( 1 , 1 ) - c
!          eigenvalmatrix ( 3 , 3 ) = eigenvalmatrix ( 3 , 3 ) + c
! !
!          write (*,*) 'eigenvalues'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') eigenvalmatrix ( i1 , : )
!          end do
! !
! !
!          jacobtest = 0.
!          jacobtest = matmul ( eigenvalmatrix , el )
!          el        = matmul ( er , jacobtest )
!          jacobtest = el
!          write (*,*) 'calculated jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacobtest ( i1 , : )
!          end do
! !
! !
!          jacob = 0.
!          jacob ( 1,1 ) = 0.
!          jacob ( 1,2 ) = 1.
!          jacob ( 1 , 3:nvar ) = 0.
! !
!          jacob ( 2 , 1 ) = - ( ( 1. - gmm ) * qq + uu*uu)
!          jacob ( 2 , 2 ) =   ( 3. - gmm ) * uu
!          jacob ( 2 , 3 ) =   ( 1. - gmm ) * vv
!          jacob ( 2 , 4 ) =   ( 1. - gmm ) * ww
!          jacob ( 2 , 5 ) = - ( 1. - gmm )
!          do s1 = 1 , nvar-5
!             jacob ( 2 , s1+5 ) = - ( 1. - gmm ) * psim ( s1 )
!          end do
! !
!          jacob ( 3 , 1 ) = - uu * vv
!          jacob ( 3 , 2 ) = vv
!          jacob ( 3 , 3 ) = uu
!          jacob ( 3 , 4:nvar ) = 0.
! !
!          jacob ( 4 , 1 ) = - uu * ww
!          jacob ( 4 , 2 ) = ww
!          jacob ( 4 , 3 ) = 0.
!          jacob ( 4 , 4 ) = uu
!          jacob ( 4 , 5:nvar ) = 0.
! !
!          jacob ( 5 , 1 ) = - ( ( 1. - gmm ) * qq + h ) * uu
!          jacob ( 5 , 2 ) =   ( 1. - gmm ) * uu * uu + h
!          jacob ( 5 , 3 ) =   ( 1. - gmm ) * uu * vv
!          jacob ( 5 , 4 ) =   ( 1. - gmm ) * uu * ww
!          jacob ( 5 , 5 ) = gmm * uu
!          do s1 = 1 , nvar-5
!             jacob ( 5 , s1+5 ) = - ( 1. - gmm ) * uu * psim ( s1 )
!          end do
! !
!          do s1 = 1 , nvar-5
!             jacob ( 5 + s1 , 1 ) = - uu * chim ( s1 )
!             jacob ( 5 + s1 , 2 ) =        chim ( s1 )
!             jacob ( 5 + s1 , 3:5 ) = 0.
!          end do
! !
!          do s1 = 6 , nvar
!             jacob ( s1 , s1 ) = uu
!             ! the rest of the elements are set to zero in the
!             ! initialisation of the jacobian
!          end do
! !
!          write (*,*) 'analytical jacobian'
!          do i1 = 1 , nvar
!             write (*,'(50(1X,1PE20.12))') jacob ( i1 , : )
!          end do
!          write (*,*) 'End of matrices'
!          read (*,*)


             ! Eigenvalues (killing the positive ones)
             L_s (1) =    max ( wz_c - cs_c , 0.0_dp )
             L_s (2) =    max ( wz_c        , 0.0_dp )
             L_s (3) =    max ( wz_c + cs_c , 0.0_dp )
             L_s (4) =    max ( wz_c        , 0.0_dp )
             L_s (5) =    max ( wz_c        , 0.0_dp )
             do m = niv+1 , nv
                L_s (m) = max ( wz_c        , 0.0_dp )
             end do


             ! Derivatives of conservative variables
             do m = 1 , nv
                dw (m) = 0.0_dp
                do st = 1 , stencil_m2
                   k_s = k_c - st + 1 ! be careful
                   dw (m) = dw (m) - cnr (st) * v (i,j,k_s,m) ! be careful
                end do
             end do


             ! Projection on the eigendirections
             do m = 1 , nv
                dwc (m) = 0.0_dp
                do mm = 1 , nv
                   dwc (m) = dwc (m) + el (m,mm) * dw (mm)
                end do
             end do


             ! Returning to conservative variables
             do m = 1 , nv
                df = 0.0_dp
                do mm = 1 , nv
                   df = df + er (m,mm) * L_s (mm) * dwc (mm)
                end do
                fl (i,j,k_c,m) = fl (i,j,k_c,m) + df * dz_i (k_c)
             end do


          end do ! end of j-loop
       end do ! end of k-loop


    end if


  end subroutine bc_noreflectionweno


end module weno
