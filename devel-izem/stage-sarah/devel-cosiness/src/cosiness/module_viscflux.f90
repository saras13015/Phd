!------------------------------------------------------------------------------
! MODULE: viscflux
!------------------------------------------------------------------------------
!> \brief Viscous fluxes.
!!
!! This module allows to calculate the viscous fluxes.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module viscflux


  use parameters
  use type_thd
  use adim
  use thermodynamics
  use deriv
  use SGS_models
  use tools , only : shock_det_slv , shock_det_ducros_slv

  implicit none

  real (dp) , parameter , private :: onethirds = 1.0_dp / 3.0_dp
  real (dp) , parameter , private :: twothirds = 2.0_dp * onethirds


contains


!> \brief 1D viscous flux.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_x ( thd , adi , grid , T , W_i , cp , ha , v , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< viscous flux


    integer (ip)                                  :: nfluid
    integer (ip)                                  :: ok , domain_id , i , j , k , l
    integer (ip)                                  :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                     :: rmu , rvd , rvd0 , rlam , rlam0 , rlamcp0 , Ya
    real (dp)                                     :: rho_i , uxs
    real (dp)                                     :: vdcx , vdtx , vdhx
    ! viscous variables
    real (dp) , dimension (:,:,:) , allocatable   :: mu , ct
    real (dp) , dimension (:,:,:,:) , allocatable :: dm
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable   :: ux
    real (dp) , dimension (:,:,:,:) , allocatable :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: q_x  , &
                                                     sx_x
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: wrk_x
    real (dp) , dimension (:,:,:,:) , allocatable :: wrka_x
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable :: fd


    nfluid  = ndim * ( ndim + 1 )
    rvd0    = adi % sqgmr / adi % Sc
    rlam0   = adi % sqgmr * adi % ggmopr
    rlamcp0 = adi % sqgmr / adi % Pr


    ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_x')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
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
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) ) ! nfluid = 2 in 1D

    do l = 1 , nrv+npv
       call dx_fixed1 ( grid % dx_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l) )
    end do


    deallocate (ux)


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species, you
    ! must specity the temperature and/or the velocity and/or the
    ! species normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       if ( bc (W) == symmetryplane ) then
          fd (sx,sy:ey,sz:ez,1) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx,sy:ey,sz:ez,2) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       if ( bc (E) == symmetryplane ) then
          fd (ex,sy:ey,sz:ez,1) = 0.0_dp ! normal velocity gradient
       end if
       fd (ex,sy:ey,sz:ez,2) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l) = 0.0_dp ! normal species gradient
       end do
    end if


    ! communicate the first derivative
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( dm   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv) , &
               mu   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               ct   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               q_x  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               sx_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_x 2')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call prim_vis_vars ( domain_id , thd , v , W_i , T , Xa , dm , mu , ct )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                uxs   = v (i,j,k,2) * rho_i

                rmu = - adi % sqgmr * mu (i,j,k)
                rvd = rvd0 * v (i,j,k,1)

                sx_x (i,j,k) = twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) )


                ! Correction velocity
                vdcx = 0.0_dp
                do l = 1 , nrv
                   vdcx = vdcx + dm (i,j,k, l )      *                         &
                          thd % Wc (l) * W_i (i,j,k) * fd (i,j,k,nfluid+ndim*l)
                end do


                vdhx = 0.0_dp
                do l = 1 , nrv

                   Ya = v (i,j,k, niv+l ) * rho_i
                   vdtx = vdcx * Ya

                   ! plus diffussion velocity
                   vdtx = vdtx - dm ( i,j,k, l )     *                         &
                          thd % Wc (l) * W_i (i,j,k) * fd (i,j,k,nfluid+ndim*l)

                   ! finally ... chi j = -cste * rho * Yalpha * V alpha j
                   fd (i,j,k,nfluid+ndim*l) = rvd * vdtx ! Xa_x

                   ! Enthalpy velocity diffusion term in energy equation
                   vdhx = vdhx + ha (i,j,k,l) * fd (i,j,k,nfluid+ndim*l)

                end do


                rlam = - rlam0 * ct (i,j,k)

                q_x (i,j,k) = uxs  * sx_x (i,j,k)   + &
                              rlam * fd   (i,j,k,2) + &
                              vdhx

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
                   mu (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) ! mu = - rlamcp
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l)
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( Xa , dm , mu , ct )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_x 3')


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
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_x 4')

    do l = 1 , nrv+npv
       call dx_fixed2 ( grid % dx_i                                                       , &
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


  end subroutine viscflux_x


!> \brief 1D viscous flux with LES models.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_LES_x ( inp , thd , adi , grid , T , W_i , cp , ha , v , mu_SGS , fl )


    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid   !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T      !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp     !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha     !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v      !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mu_SGS !< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl     !< viscous flux


    integer (ip)                                  :: nfluid , nSGS
    integer (ip)                                  :: ok , domain_id , i , j , k , l
    integer (ip)                                  :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                     :: wrk , rmu , rmu_SGS , rvd , rvd0 , rlam , rlam0 , rlamcp0 , &
                                                     Sc_SGS_i , Pr_SGS_i , gm1_g , ct_SGS , dm_SGS , dmrw
    real (dp)                                     :: vdcx , vdtx , vdhx , dxadx
    real (dp)                                     :: rho_i , uxs , Ya
    ! viscous variables
    real (dp) , dimension (:,:,:) , allocatable   :: mu , ct
    real (dp) , dimension (:,:,:,:) , allocatable :: dm
    ! SGS additional variables
    real (dp) , dimension (:,:,:) , allocatable   :: tau_iso_SGS
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable   :: ux
    real (dp) , dimension (:,:,:,:) , allocatable :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: q_x  , &
                                                     sx_x
    real (dp) , dimension (:,:,:) , allocatable   :: sx_x_SGS
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: wrk_x
    real (dp) , dimension (:,:,:,:) , allocatable :: wrka_x
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable :: fd


    nfluid  = ndim * ( ndim + 1 )
    nSGS    = nfluid + ndim * ( nrv + npv + nvv )

    rvd0    = adi % sqgmr / adi % Sc     ! M * sqrt(gamma) / ( Re * Sc )
    rlam0   = adi % sqgmr * adi % ggmopr ! M * sqrt(gamma) * gamma / ( Re * (gamma-1) * Pr )
    rlamcp0 = adi % sqgmr / adi % Pr     ! M * sqrt(gamma) / ( Re * Pr )

    Sc_SGS_i = 1.0_dp / inp % Sc_SGS
    Pr_SGS_i = 1.0_dp / inp % Pr_SGS
    gm1_g    = adi % gm1 / adi % gamma

    ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)             , &
               Xa          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv+nvv) , &
               fd          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)      , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_x')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
             end do
          end do
       end do

    end do


    ! passive scalar and variance equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do l = nrv+1,nrv+npv+nvv
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
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) ) ! dT/dx = fd (:,:,:,2) nfluid = 2 in 1D


    ! 1st derivative for molar fraction, passive scalar and variance
    do l = 1 , nrv+npv+nvv
       call dx_fixed1 ( grid % dx_i                                                   , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)             , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l) ) ! dXa/dx
    end do


    ! SGS: W_i derivatives need to be communicated
    call dx_fixed1 ( grid % dx_i , W_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nSGS+1) )    ! dW_i/dx

    ! SGS: shock detector criteria needs to be communicated
    fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) = 1.0_dp
    call shock_det_slv ( v , T , W_i , fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) )


    deallocate (ux)


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species, you
    ! must specity the temperature and/or the velocity and/or the
    ! species normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       if ( bc (W) == symmetryplane ) then
          fd (sx,sy:ey,sz:ez,1) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx,sy:ey,sz:ez,2) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv+nvv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l) = 0.0_dp ! normal species gradient
       end do
       fd (sx,sy:ey,sz:ez,nSGS+1)      = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       if ( bc (E) == symmetryplane ) then
          fd (ex,sy:ey,sz:ez,1) = 0.0_dp ! normal velocity gradient
       end if
       fd (ex,sy:ey,sz:ez,2) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv+nvv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l) = 0.0_dp ! normal species gradient
       end do
       fd (ex,sy:ey,sz:ez,nSGS+1)      = 0.0_dp ! normal mixture molecular weight gradient
    end if


    ! communicate the first derivatives
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( dm           (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)  , &
               mu           (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               ct           (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               q_x          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_x         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_x_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               tau_iso_SGS  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_x 2')


    ! SGS: viscosity does NOT need to be communicated
    do domain_id = -ndim , ndim
       call mu_SGS_selector ( inp , adi , domain_id , grid , v , fd , mu_SGS )
    end do

    ! SGS: tau_iso_SGS derivatives need to be COMMUNICATED
    if ( inp % tau_iso_SGS_switch ) then ! consider the isotropic part of SGS viscous stress tensor

       do domain_id = -ndim , ndim
          call tau_iso_SGS_selector ( inp , domain_id , grid , v , fd , tau_iso_SGS )
       end do

    else ! neglects the isotropic part of SGS viscous stress tensor

       tau_iso_SGS (:,:,:) = 0.0_dp

    end if


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call prim_vis_vars ( domain_id , thd , v , W_i , T , Xa , dm , mu , ct )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                uxs   = v (i,j,k,2) * rho_i

                ct_SGS = adi % Pr * mu_SGS (i,j,k) * cp (i,j,k) * Pr_SGS_i
                dm_SGS = adi % Sc * mu_SGS (i,j,k) * rho_i * Sc_SGS_i

                rmu          = - adi % sqgmr * mu (i,j,k)
                rmu_SGS      = - adi % sqgmr * mu_SGS (i,j,k) / rmu

                rvd = rvd0 * v (i,j,k,1)

                ! Viscous stress tensor: -M*\sqrt(gamma_inf)/Re * tau_ij
                sx_x (i,j,k) = twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) )
                sx_x_SGS (i,j,k) = sx_x (i,j,k) * ( 1.0_dp + rmu_SGS ) + onethirds * tau_iso_SGS (i,j,k)

                ! Correction velocity: sum_{beta=1,N} D_{beta m} W_beta/W \pd{Xbeta}{x_j}
                vdcx = 0.0_dp
                do l = 1 , nrv
                   dmrw = dm (i,j,k, l ) * thd % Wc (l) * W_i (i,j,k) ! No SGS term
                   vdcx = vdcx + dmrw * fd (i,j,k,nfluid+ndim*l)
                end do

                vdhx = 0.0_dp
                do l = 1 , nrv

                   Ya = v (i,j,k, niv+l ) * rho_i
                   vdtx = vdcx * Ya
                   dxadx = fd (i,j,k,nfluid+ndim*l)

                   dmrw = - dm ( i,j,k, l ) * thd % Wc (l) * W_i (i,j,k)

                   ! Enthalpy velocity diffusion term in energy equation (No SGS term)
                   ! - M*\sqrt(gamma_inf)/(Re*Le*Pr) * sum_{alpha=1,N} rho * Yalpha * V alpha j * halpha
                   vdhx = vdhx + ha (i,j,k,l) * rvd * ( vdtx + dmrw * dxadx )

                   wrk = dmrw - dm_SGS * thd % Wc (l) * W_i (i,j,k)

                   ! plus diffussion velocity
                   vdtx = vdtx + wrk * dxadx

                   wrk = - dm_SGS * thd % Wc (l) * Xa (i,j,k,l)

                   ! plus an extra SGS term
                   vdtx = vdtx + wrk * fd (i,j,k,nSGS+1)

                   ! finally ... chi j = -cste * rho * Yalpha * V alpha j
                   fd (i,j,k,nfluid+ndim*l) = rvd * vdtx ! Xa_x

                end do

                rlam = - rlam0 * ( ct (i,j,k) + ct_SGS )

                q_x (i,j,k) = uxs  * sx_x_SGS (i,j,k)       + &
                              rlam * fd   (i,j,k,2)         + &
                              vdhx

             end do
          end do
       end do

    end do

    deallocate ( tau_iso_SGS )
    deallocate ( sx_x )

    ! passive scalar and variance equations: viscous term
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   dm_SGS = adi % Sc * mu_SGS (i,j,k) / v (i,j,k,1) * Sc_SGS_i
                   mu (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) & ! > laminar/molecular term
                                - rvd0 * v (i,j,k,1) * dm_sgs         ! > SGS term for LES
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv+nvv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l)
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( Xa , dm , mu , ct )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_x 3')


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


    allocate ( wrka_x (sx:ex,sy:ey,sz:ez,nrv+npv+nvv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_x 4')


    ! viscous flux in rho*Y_alpha equation -> fl (:,:,:,niv+1:niv+npv+nvv)
    do l = 1 , nrv+npv+nvv
       call dx_fixed2 ( grid % dx_i                                                       , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                               )
    end do

    do l = 1 , nrv+npv+nvv
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


  end subroutine viscflux_LES_x


!> \brief 2D viscous flux.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_xy ( thd , adi , grid , T , W_i , cp , ha , v , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< viscous flux


    integer (ip)                                  :: nfluid
    integer (ip)                                  :: ok , domain_id , i , j , k , l
    integer (ip)                                  :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                     :: wrk , rmu , rvd , rvd0 , rlam , rlam0 , rlamcp0 , Ya
    real (dp)                                     :: rho_i , uxs , vys
    real (dp)                                     :: vdcx , vdtx , vdhx , &
                                                     vdcy , vdty , vdhy
    ! viscous variables
    real (dp) , dimension (:,:,:) , allocatable   :: mu , ct
    real (dp) , dimension (:,:,:,:) , allocatable :: dm
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable   :: ux , vy
    real (dp) , dimension (:,:,:,:) , allocatable :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: q_x  , q_y  , &
                                                     sx_x , sx_y , &
                                                     sy_y
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: wrk_x , wrk_y
    real (dp) , dimension (:,:,:,:) , allocatable :: wrka_x , wrka_y
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable :: fd

    ! wallinjection variables
    real (dp)                       :: rinj , x0 , rrr

    nfluid  = ndim * ( ndim + 1 )
    rvd0    = adi % sqgmr / adi % Sc
    rlam0   = adi % sqgmr * adi % ggmopr
    rlamcp0 = adi % sqgmr / adi % Pr


    ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_xy')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
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

    call dx_fixed1 ( grid % dx_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,3) ) ! dv/dx = fd (:,:,:,3)
    call dy_fixed1 ( grid % dy_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,4) ) ! dv/dy = fd (:,:,:,4)

    call dx_fixed1 ( grid % dx_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,5) ) ! dT/dx = fd (:,:,:,5)
    call dy_fixed1 ( grid % dy_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,6) ) ! dT/dy = fd (:,:,:,6) nfluid = 6 in 2D


    do l = 1 , nrv+npv
       call dx_fixed1 ( grid % dx_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) ) ! dXa/dx
       call dy_fixed1 ( grid % dy_i                                       , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l) , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) ) ! dXa/dy
    end do


    deallocate ( ux , vy )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species, you
    ! must specity the temperature and/or the velocity and/or the
    ! species normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       if ( bc (W) == symmetryplane ) then
          fd (sx,sy:ey,sz:ez,1) = 0.0_dp
          fd (sx,sy:ey,sz:ez,3) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx,sy:ey,sz:ez,5) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       if ( bc (E) == symmetryplane ) then
          fd (ex,sy:ey,sz:ez,1) = 0.0_dp
          fd (ex,sy:ey,sz:ez,3) = 0.0_dp ! normal velocity gradient
       end if
       fd (ex,sy:ey,sz:ez,5) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
    end if

    if ( bc (S) == symmetryplane .or. bc (S) == adiabaticwall ) then
       if ( bc (S) == symmetryplane ) then
          fd (sx:ex,sy,sz:ez,2) = 0.0_dp
          fd (sx:ex,sy,sz:ez,4) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,sy,sz:ez,6) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy,sz:ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (N) == symmetryplane .or. bc (N) == adiabaticwall ) then
       if ( bc (N) == symmetryplane ) then
          fd (sx:ex,ey,sz:ez,2) = 0.0_dp
          fd (sx:ex,ey,sz:ez,4) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,ey,sz:ez,6) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (sx:ex,ey,sz:ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
    end if

    if ( bc (N) == wallinjection ) then
       rinj = 0.0010_dp / adi % L_ref ! if it modified, changed also in module_BCs.f90 !!!
       x0   = 0.2000_dp / adi % L_ref
       do i = sx , ex
          rrr = sqrt ( (grid % x(i)-x0)*(grid % x(i)-x0) )
          if ( rrr > rinj ) then ! wall/slipwall
!             fd (i,ey,sz:ez,2) = 0.0_dp
!             fd (i,ey,sz:ez,4) = 0.0_dp ! normal velocity gradient
             fd (i,ey,sz:ez,6) = 0.0_dp ! normal temperature gradient
             do l = 1 , nrv+npv
                fd (i,ey,sz:ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
             end do
          end if
       end do
    end if

    if ( bc (S) == wallinjection ) then
       rinj = 0.0010_dp / adi % L_ref ! if it modified, changed also in module_BCs.f90 !!!
       x0   = 0.2000_dp / adi % L_ref
       do i = sx , ex
          rrr = sqrt ( (grid % x(i)-x0)*(grid % x(i)-x0) )
          if ( rrr > rinj ) then ! wall/slipwall
!             fd (i,sy,sz:ez,2) = 0.0_dp
!             fd (i,sy,sz:ez,4) = 0.0_dp ! normal velocity gradient
             fd (i,sy,sz:ez,6) = 0.0_dp ! normal temperature gradient
             do l = 1 , nrv+npv
                fd (i,sy,sz:ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
             end do
          end if
       end do
    end if

    ! communicate the first derivative
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( dm   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)     , &
               mu   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               ct   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               q_x  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               q_y  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               sx_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               sx_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               sy_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_xy 2')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call prim_vis_vars ( domain_id , thd , v , W_i , T , Xa , dm , mu , ct )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                uxs   = v (i,j,k,2) * rho_i
                vys   = v (i,j,k,3) * rho_i

                rmu = - adi % sqgmr * mu (i,j,k)
                rvd = rvd0 * v (i,j,k,1)

                sx_x (i,j,k) = twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) - &
                                                   fd (i,j,k,4) )

                sx_y (i,j,k) = rmu * ( fd (i,j,k,3) + fd (i,j,k,2) )

                sy_y (i,j,k) = twothirds * rmu * ( fd (i,j,k,4) + fd (i,j,k,4) - &
                                                   fd (i,j,k,1) )


                ! Correction velocity
                vdcx = 0.0_dp ; vdcy = 0.0_dp
                do l = 1 , nrv
                   wrk  = dm (i,j,k, l ) * thd % Wc (l) * W_i (i,j,k)
                   vdcx = vdcx + wrk * fd (i,j,k,nfluid+ndim*l-1)
                   vdcy = vdcy + wrk * fd (i,j,k,nfluid+ndim*l  )
                end do


                vdhx = 0.0_dp ; vdhy = 0.0_dp
                do l = 1 , nrv

                   Ya = v (i,j,k, niv+l ) * rho_i
                   vdtx = vdcx * Ya ; vdty = vdcy * Ya

                   wrk = - dm ( i,j,k, l ) * thd % Wc (l) * W_i (i,j,k)

                   ! plus diffussion velocity
                   vdtx = vdtx + wrk * fd (i,j,k,nfluid+ndim*l-1)
                   vdty = vdty + wrk * fd (i,j,k,nfluid+ndim*l  )

                   ! finally ... chi j = -cste * rho * Yalpha * V alpha j
                   fd (i,j,k,nfluid+ndim*l-1) = rvd * vdtx ! Xa_x
                   fd (i,j,k,nfluid+ndim*l  ) = rvd * vdty ! Xa_y

                   ! Enthalpy velocity diffusion term in energy equation
                   vdhx = vdhx + ha (i,j,k,l) * fd (i,j,k,nfluid+ndim*l-1)
                   vdhy = vdhy + ha (i,j,k,l) * fd (i,j,k,nfluid+ndim*l  )

                end do


                rlam = - rlam0 * ct (i,j,k)

                q_x (i,j,k) = uxs  * sx_x (i,j,k)   + &
                              vys  * sx_y (i,j,k)   + &
                              rlam * fd   (i,j,k,5) + &
                              vdhx

                q_y (i,j,k) = uxs  * sx_y (i,j,k)   + &
                              vys  * sy_y (i,j,k)   + &
                              rlam * fd   (i,j,k,6) + &
                              vdhy

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
                   mu (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) ! mu = -rlamcp
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l-1) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l-1)
                      fd (i,j,k,nfluid+ndim*l  ) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l  )
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( Xa , dm , mu , ct )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               wrk_y (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_xy 3')


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
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_xy 4')

    do l = 1 , nrv+npv
       call dx_fixed2 ( grid % dx_i                                                         , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                                 )
       call dy_fixed2 ( grid % dy_i                                                         , &
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


  end subroutine viscflux_xy


!> \brief 2D viscous flux with LES models.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_LES_xy ( inp , thd , adi , grid , T , W_i , cp , ha , v , mu_SGS , fl )


    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid   !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T      !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp     !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha     !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v      !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mu_SGS !< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl     !< viscous flux


    integer (ip)                                  :: nfluid , nSGS
    integer (ip)                                  :: ok , domain_id , i , j , k , l
    integer (ip)                                  :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                     :: wrk , rmu , rmu_SGS , rvd , rvd0 , rlam , rlam0 , rlamcp0 , &
                                                     Sc_SGS_i , Pr_SGS_i , gm1_g , ct_SGS , dm_SGS , dmrw
    real (dp)                                     :: vdcx , vdtx , vdhx , dxadx , &
                                                     vdcy , vdty , vdhy , dxady
    real (dp)                                     :: rho_i , uxs , vys , Ya
    ! viscous variables
    real (dp) , dimension (:,:,:) , allocatable   :: mu , ct
    real (dp) , dimension (:,:,:,:) , allocatable :: dm
    ! SGS additional variables
    real (dp) , dimension (:,:,:) , allocatable   :: tau_iso_SGS
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable   :: ux , vy
    real (dp) , dimension (:,:,:,:) , allocatable :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: q_x  , q_y  , &
                                                     sx_x , sx_y , &
                                                     sy_y
    real (dp) , dimension (:,:,:) , allocatable   :: sx_x_SGS , sx_y_SGS , &
                                                     sy_y_SGS
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: wrk_x , wrk_y
    real (dp) , dimension (:,:,:,:) , allocatable :: wrka_x , wrka_y
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable :: fd

    ! BC_wallinjection variables
    real (dp)                       :: rinj , x0 , rrr

    ! Constants
    real (dp) , parameter         :: C_xi      = 1.0_dp          , &
                                     eps       = 1.0e-15_dp


    nfluid  = ndim * ( ndim + 1 )
    nSGS    = nfluid + ndim * ( nrv + npv + nvv )

    rvd0    = adi % sqgmr / adi % Sc     ! M * sqrt(gamma) / ( Re * Sc )
    rlam0   = adi % sqgmr * adi % ggmopr ! M * sqrt(gamma) * gamma / ( Re * (gamma-1) * Pr )
    rlamcp0 = adi % sqgmr / adi % Pr     ! M * sqrt(gamma) / ( Re * Pr )

    Sc_SGS_i = 1.0_dp / inp % Sc_SGS
    Pr_SGS_i = 1.0_dp / inp % Pr_SGS
    gm1_g    = adi % gm1 / adi % gamma

    ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
    allocate ( ux          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)             , &
               vy          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)             , &
               Xa          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv+nvv) , &
               fd          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)      , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_xy')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call molefrac ( domain_id , thd , W_i , v , Xa )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
             end do
          end do
       end do

    end do


    ! passive scalar and variance equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do l = nrv+1,nrv+npv+nvv
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
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,6) ) ! dT/dy = fd (:,:,:,6) nfluid = 6 in 2D


    ! 1st derivative for molar fraction, passive scalar and variance
    do l = 1 , nrv+npv+nvv
       call dx_fixed1 ( grid % dx_i                                                     , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)               , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) ) ! dXa/dx
       call dy_fixed1 ( grid % dy_i                                                     , &
                        Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)               , &
                        fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) ) ! dXa/dy
    end do


    ! SGS: W_i derivatives need to be communicated
    call dx_fixed1 ( grid % dx_i , W_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nSGS+1) ) ! dW_i/dx = fd (:,:,:,nSGS+1)
    call dy_fixed1 ( grid % dy_i , W_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nSGS+2) ) ! dW_i/dy = fd (:,:,:,nSGS+2)

    ! SGS: shock detector criteria needs to be communicated
    fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) = 1.0_dp
!    call shock_det_slv ( v , T , W_i , fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) )
    call shock_det_ducros_slv ( 0 , fd , fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) )


    deallocate ( ux , vy )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species, you
    ! must specity the temperature and/or the velocity and/or the
    ! species normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       if ( bc (W) == symmetryplane ) then
          fd (sx,sy:ey,sz:ez,1) = 0.0_dp
          fd (sx,sy:ey,sz:ez,3) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx,sy:ey,sz:ez,5) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv+nvv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
       fd (sx,sy:ey,sz:ez,nSGS+1)      = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       if ( bc (E) == symmetryplane ) then
          fd (ex,sy:ey,sz:ez,1) = 0.0_dp
          fd (ex,sy:ey,sz:ez,3) = 0.0_dp ! normal velocity gradient
       end if
       fd (ex,sy:ey,sz:ez,5) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv+nvv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
       fd (ex,sy:ey,sz:ez,nSGS+1)      = 0.0_dp ! normal mixture molecular weight gradient
    end if

    if ( bc (S) == symmetryplane .or. bc (S) == adiabaticwall ) then
       if ( bc (S) == symmetryplane ) then
          fd (sx:ex,sy,sz:ez,2) = 0.0_dp
          fd (sx:ex,sy,sz:ez,4) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,sy,sz:ez,6) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv+nvv
          fd (sx:ex,sy,sz:ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
       fd (sx:ex,sy,sz:ez,nSGS+2)      = 0.0_dp ! normal mixture molecular weight gradient
    end if
    if ( bc (N) == symmetryplane .or. bc (N) == adiabaticwall ) then
       if ( bc (N) == symmetryplane ) then
          fd (sx:ex,ey,sz:ez,2) = 0.0_dp
          fd (sx:ex,ey,sz:ez,4) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,ey,sz:ez,6) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv+nvv
          fd (sx:ex,ey,sz:ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
       fd (sx:ex,ey,sz:ez,nSGS+2)      = 0.0_dp ! normal mixture molecular weight gradient
    end if

    if ( bc (N) == wallinjection ) then
       rinj = 0.0010_dp / adi % L_ref ! if it modified, changed also in module_BCs.f90 !!!
       x0   = 0.2000_dp / adi % L_ref
       do i = sx , ex
          rrr = sqrt ( (grid % x(i)-x0)*(grid % x(i)-x0) )
          if ( rrr > rinj ) then ! wall/slipwall
!             fd (i,ey,sz:ez,2) = 0.0_dp
!             fd (i,ey,sz:ez,4) = 0.0_dp ! normal velocity gradient
             fd (i,ey,sz:ez,6) = 0.0_dp ! normal temperature gradient
             do l = 1 , nrv+npv+nvv
                fd (i,ey,sz:ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
             end do
             fd (i,ey,sz:ez,nSGS+2)      = 0.0_dp ! normal mixture molecular weight gradient
          end if
       end do
    end if

    if ( bc (S) == wallinjection ) then
       rinj = 0.0010_dp / adi % L_ref ! if it modified, changed also in module_BCs.f90 !!!
       x0   = 0.2000_dp / adi % L_ref
       do i = sx , ex
          rrr = sqrt ( (grid % x(i)-x0)*(grid % x(i)-x0) )
          if ( rrr > rinj ) then ! wall/slipwall
!             fd (i,sy,sz:ez,2) = 0.0_dp
!             fd (i,sy,sz:ez,4) = 0.0_dp ! normal velocity gradient
             fd (i,sy,sz:ez,6) = 0.0_dp ! normal temperature gradient
             do l = 1 , nrv+npv+nvv
                fd (i,sy,sz:ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
             end do
             fd (i,sy,sz:ez,nSGS+2)      = 0.0_dp ! normal mixture molecular weight gradient
          end if
       end do
    end if


    ! communicate the first derivatives
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( dm           (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)  , &
               mu           (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               ct           (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               q_x          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               q_y          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_x         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_y         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sy_y         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_x_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_y_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sy_y_SGS     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               tau_iso_SGS  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_xy 2')


    ! SGS: viscosity does NOT need to be communicated
    do domain_id = -ndim , ndim
       call mu_SGS_selector ( inp , adi , domain_id , grid , v , fd , mu_SGS )
    end do

    ! SGS: tau_iso_SGS derivatives need to be COMMUNICATED
    if ( inp % tau_iso_SGS_switch ) then ! consider the isotropic part of SGS viscous stress tensor

       do domain_id = -ndim , ndim
          call tau_iso_SGS_selector ( inp , domain_id , grid , v , fd , tau_iso_SGS )
       end do

    else ! neglects the isotropic part of SGS viscous stress tensor

       tau_iso_SGS (:,:,:) = 0.0_dp

    end if


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call prim_vis_vars ( domain_id , thd , v , W_i , T , Xa , dm , mu , ct )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                uxs   = v (i,j,k,2) * rho_i
                vys   = v (i,j,k,3) * rho_i

                ct_SGS = adi % Pr * mu_SGS (i,j,k) * cp (i,j,k) * Pr_SGS_i
                dm_SGS = adi % Sc * mu_SGS (i,j,k) * rho_i * Sc_SGS_i

                rmu          = - adi % sqgmr * mu (i,j,k)
                rmu_SGS      = - adi % sqgmr * mu_SGS (i,j,k) / rmu

                rvd = rvd0 * v (i,j,k,1)

                ! Viscous stress tensor: -M*\sqrt(gamma_inf)/Re * tau_ij
                sx_x (i,j,k) = twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) - &
                                                   fd (i,j,k,4) )
                sx_x_SGS (i,j,k) = sx_x (i,j,k) * ( 1.0_dp + rmu_SGS ) + onethirds * tau_iso_SGS (i,j,k)

                sx_y (i,j,k) = rmu * ( fd (i,j,k,3) + fd (i,j,k,2) )
                sx_y_SGS (i,j,k) = sx_y (i,j,k) * ( 1.0_dp + rmu_SGS )

                sy_y (i,j,k) = twothirds * rmu * ( fd (i,j,k,4) + fd (i,j,k,4) - &
                                                   fd (i,j,k,1) )
                sy_y_SGS (i,j,k) = sy_y (i,j,k) * ( 1.0_dp + rmu_SGS ) + onethirds * tau_iso_SGS (i,j,k)

                ! Correction velocity: sum_{beta=1,N} D_{beta m} W_beta/W \pd{Xbeta}{x_j}
                vdcx = 0.0_dp ; vdcy = 0.0_dp
                do l = 1 , nrv
                   dmrw = dm (i,j,k, l ) * thd % Wc (l) * W_i (i,j,k) ! No SGS term
                   vdcx = vdcx + dmrw * fd (i,j,k,nfluid+ndim*l-1)
                   vdcy = vdcy + dmrw * fd (i,j,k,nfluid+ndim*l  )
                end do

                vdhx = 0.0_dp ; vdhy = 0.0_dp
                do l = 1 , nrv

                   Ya = v (i,j,k, niv+l ) * rho_i
                   vdtx = vdcx * Ya ; vdty = vdcy * Ya
                   dxadx = fd (i,j,k,nfluid+ndim*l-1)
                   dxady = fd (i,j,k,nfluid+ndim*l  )

                   dmrw = - dm ( i,j,k, l ) * thd % Wc (l) * W_i (i,j,k)

                   ! Enthalpy velocity diffusion term in energy equation (No SGS term)
                   ! - M*\sqrt(gamma_inf)/(Re*Le*Pr) * sum_{alpha=1,N} rho * Yalpha * V alpha j * halpha
                   vdhx = vdhx + ha (i,j,k,l) * rvd * ( vdtx + dmrw * dxadx )
                   vdhy = vdhy + ha (i,j,k,l) * rvd * ( vdty + dmrw * dxady )

                   wrk = dmrw - dm_SGS * thd % Wc (l) * W_i (i,j,k)

                   ! plus diffussion velocity
                   vdtx = vdtx + wrk * dxadx
                   vdty = vdty + wrk * dxady

                   wrk = - dm_SGS * thd % Wc (l) * Xa (i,j,k,l)

                   ! plus an extra SGS term
                   vdtx = vdtx + wrk * fd (i,j,k,nSGS+1)
                   vdty = vdty + wrk * fd (i,j,k,nSGS+2)

                   ! finally ... chi j = -cste * rho * Yalpha * V alpha j
                   fd (i,j,k,nfluid+ndim*l-1) = rvd * vdtx ! Xa_x
                   fd (i,j,k,nfluid+ndim*l  ) = rvd * vdty ! Xa_y

                end do

                rlam = - rlam0 * ( ct (i,j,k) + ct_SGS )

                q_x (i,j,k) = uxs  * sx_x_SGS (i,j,k)       + &
                              vys  * sx_y_SGS (i,j,k)       + &
                              rlam * fd (i,j,k,5)           + &
                              vdhx

                q_y (i,j,k) = uxs  * sx_y_SGS (i,j,k)       + &
                              vys  * sy_y_SGS (i,j,k)       + &
                              rlam * fd (i,j,k,6)           + &
                              vdhy

             end do
          end do
       end do

    end do

    deallocate ( tau_iso_SGS )
    deallocate ( sx_x , sx_y , &
                        sy_y   )

    ! passive scalar and variance equations: viscous term
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   dm_SGS = adi % Sc * mu_SGS (i,j,k) / v (i,j,k,1) * Sc_SGS_i
                   mu (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) & ! > laminar/molecular term
                                - rvd0 * v (i,j,k,1) * dm_sgs         ! > SGS term for LES
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv+nvv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l-1) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l-1)
                      fd (i,j,k,nfluid+ndim*l  ) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l  )
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( Xa , dm , mu , ct )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               wrk_y (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_xy 3')


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


    ! viscous flux in rho*et equation -> fl (:,:,:,5)
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


    allocate ( wrka_x (sx:ex,sy:ey,sz:ez,nrv+npv+nvv) , &
               wrka_y (sx:ex,sy:ey,sz:ez,nrv+npv+nvv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_xy 4')


    ! viscous flux in rho*Y_alpha equation -> fl (:,:,:,niv+1:niv+npv+nvv)
    do l = 1 , nrv+npv+nvv
       call dx_fixed2 ( grid % dx_i                                                         , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                                 )
       call dy_fixed2 ( grid % dy_i                                                         , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) , &
                        wrka_y (sx:ex,sy:ey,sz:ez,l)                                 )
    end do

    do l = 1 , nrv+npv+nvv
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


  end subroutine viscflux_LES_xy


!> \brief 3D viscous flux.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_xyz ( thd , adi , grid , T , W_i , cp , ha , v , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< viscous flux


    integer (ip)                                  :: nfluid
    integer (ip)                                  :: ok , domain_id , i , j , k , l
    integer (ip)                                  :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                     :: wrk , rmu , rvd , rvd0 , rlam , rlam0 , rlamcp0 , Ya
    real (dp)                                     :: vdcx , vdtx , vdhx , &
                                                     vdcy , vdty , vdhy , &
                                                     vdcz , vdtz , vdhz
    real (dp)                                     :: rho_i , uxs , vys , wzs
    ! viscous variables
    real (dp) , dimension (:,:,:) , allocatable   :: mu , ct
    real (dp) , dimension (:,:,:,:) , allocatable :: dm
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable   :: ux , vy , wz
    real (dp) , dimension (:,:,:,:) , allocatable :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: q_x  , q_y  , q_z  , &
                                                     sx_x , sx_y , sx_z , &
                                                     sy_y , sy_z , sz_z
    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: wrk_x , wrk_y , wrk_z
    real (dp) , dimension (:,:,:,:) , allocatable :: wrka_x , wrka_y , wrka_z
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable :: fd

    ! BC_wallinjection variables
    real (dp)                       :: rinj , x0 , z0 , rrr


    nfluid  = ndim * ( ndim + 1 )
    rvd0    = adi % sqgmr / adi % Sc     ! = Ma sqrt(gamma) / ( Re Sc )
    rlam0   = adi % sqgmr * adi % ggmopr ! = Ma sqrt(gamma) gamma / ( Re (gamma-1) Pr )
    rlamcp0 = adi % sqgmr / adi % Pr     ! = Ma sqrt(gamma) / ( Re Pr )


    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
               Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_xyz')


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
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,12) ) ! dT/dz = fd (:,:,:,12) nfluid = 12 in 3D


    ! 1st derivative for molar fraction
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


    deallocate ( ux , vy , wz )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species, you
    ! must specity the temperature and/or the velocity and/or the
    ! species normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       if ( bc (W) == symmetryplane ) then
          fd (sx,sy:ey,sz:ez,1) = 0.0_dp
          fd (sx,sy:ey,sz:ez,4) = 0.0_dp
          fd (sx,sy:ey,sz:ez,7) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx,sy:ey,sz:ez,10) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (sx,sy:ey,sz:ez,nfluid+ndim*l-2) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (E) == symmetryplane .or. bc (E) == adiabaticwall ) then
       if ( bc (E) == symmetryplane ) then
          fd (ex,sy:ey,sz:ez,1) = 0.0_dp
          fd (ex,sy:ey,sz:ez,4) = 0.0_dp
          fd (ex,sy:ey,sz:ez,7) = 0.0_dp ! normal velocity gradient
       end if
       fd (ex,sy:ey,sz:ez,10) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (ex,sy:ey,sz:ez,nfluid+ndim*l-2) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (S) == symmetryplane .or. bc (S) == adiabaticwall ) then
       if ( bc (S) == symmetryplane ) then
          fd (sx:ex,sy,sz:ez,2) = 0.0_dp
          fd (sx:ex,sy,sz:ez,5) = 0.0_dp
          fd (sx:ex,sy,sz:ez,8) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,sy,sz:ez,11) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (N) == symmetryplane .or. bc (N) == adiabaticwall ) then
       if ( bc (N) == symmetryplane ) then
          fd (sx:ex,ey,sz:ez,2) = 0.0_dp
          fd (sx:ex,ey,sz:ez,5) = 0.0_dp
          fd (sx:ex,ey,sz:ez,8) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,ey,sz:ez,11) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (sx:ex,ey,sz:ez,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (B) == symmetryplane .or. bc (B) == adiabaticwall ) then
       if ( bc (B) == symmetryplane ) then
          fd (sx:ex,sy:ey,sz,3) = 0.0_dp
          fd (sx:ex,sy:ey,sz,6) = 0.0_dp
          fd (sx:ex,sy:ey,sz,9) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,sy:ey,sz,12) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy:ey,sz,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
    end if
    if ( bc (F) == symmetryplane .or. bc (F) == adiabaticwall ) then
       if ( bc (F) == symmetryplane ) then
          fd (sx:ex,sy:ey,ez,3) = 0.0_dp
          fd (sx:ex,sy:ey,ez,6) = 0.0_dp
          fd (sx:ex,sy:ey,ez,9) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx:ex,sy:ey,ez,12) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv
          fd (sx:ex,sy:ey,ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
    end if

    if ( bc (N) == wallinjection ) then
       rinj = 0.0010_dp / adi % L_ref ! if it modified, changed also in module_BCs.f90 !!!
       x0   = 0.2000_dp / adi % L_ref
       z0   = 0.0000_dp / adi % L_ref
       do k = sz , ez
          do i = sx , ex
             rrr = sqrt ( (grid % x(i)-x0)*(grid % x(i)-x0) + (grid % z(k)-z0)*(grid % z(k)-z0) )
             if ( rrr > rinj ) then ! wall/slipwall
!                fd (i,ey,k,2) = 0.0_dp
!                fd (i,ey,k,5) = 0.0_dp
!                fd (i,ey,k,8) = 0.0_dp ! normal velocity gradient
                fd (i,ey,k,11) = 0.0_dp ! normal temperature gradient
                do l = 1 , nrv+npv
                   fd (i,ey,k,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
                end do
             end if
          end do
       end do
    end if

    if ( bc (S) == wallinjection ) then
       rinj = 0.0010_dp / adi % L_ref ! if it modified, changed also in module_BCs.f90 !!!
       x0   = 0.2000_dp / adi % L_ref
       z0   = 0.0000_dp / adi % L_ref
       do k = sz , ez
          do i = sx , ex
             rrr = sqrt ( (grid % x(i)-x0)*(grid % x(i)-x0) + (grid % z(k)-z0)*(grid % z(k)-z0) )
             if ( rrr > rinj ) then ! wall/slipwall
!                fd (i,sy,k,2) = 0.0_dp
!                fd (i,sy,k,5) = 0.0_dp
!                fd (i,sy,k,8) = 0.0_dp ! normal velocity gradient
                fd (i,sy,k,11) = 0.0_dp ! normal temperature gradient
                do l = 1 , nrv+npv
                   fd (i,sy,k,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
                end do
             end if
          end do
       end do
    end if


    ! communicate the first derivative
    call comm_deriv (fd)


    ! components of the viscous stresses
    allocate ( dm   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv) , &
               mu   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               ct   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               q_x  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               q_y  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               q_z  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               sx_x (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               sx_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               sx_z (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               sy_y (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               sy_z (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               sz_z (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_xyz 2')


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call prim_vis_vars ( domain_id , thd , v , W_i , T , Xa , dm , mu , ct )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                uxs   = v (i,j,k,2) * rho_i
                vys   = v (i,j,k,3) * rho_i
                wzs   = v (i,j,k,4) * rho_i

                rmu = - adi % sqgmr * mu (i,j,k)
                rvd = rvd0 * v (i,j,k,1)


                ! Viscous stress tensor: -M*\sqrt(gamma_inf)/Re * tau_ij
                sx_x (i,j,k) = twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1) - &
                                                 ( fd (i,j,k,5) + fd (i,j,k,9) ) )

                sx_y (i,j,k) = rmu * ( fd (i,j,k,4) + fd (i,j,k,2) )

                sx_z (i,j,k) = rmu * ( fd (i,j,k,7) + fd (i,j,k,3) )

                sy_y (i,j,k) = twothirds * rmu * ( fd (i,j,k,5) + fd (i,j,k,5) - &
                                                 ( fd (i,j,k,1) + fd (i,j,k,9) ) )

                sy_z (i,j,k) = rmu * ( fd (i,j,k,8) + fd (i,j,k,6) )

                sz_z (i,j,k) = twothirds * rmu * ( fd (i,j,k,9) + fd (i,j,k,9) - &
                                                 ( fd (i,j,k,1) + fd (i,j,k,5) ) )


                ! Correction velocity: sum_{beta=1,N} D_{beta m} W_beta/W \pd{Ybeta}{x_j}
                vdcx = 0.0_dp ; vdcy = 0.0_dp ; vdcz = 0.0_dp
                do l = 1 , nrv
                   wrk  = dm (i,j,k, l ) * thd % Wc (l) * W_i (i,j,k)
                   vdcx = vdcx + wrk * fd (i,j,k,nfluid+ndim*l-2)
                   vdcy = vdcy + wrk * fd (i,j,k,nfluid+ndim*l-1)
                   vdcz = vdcz + wrk * fd (i,j,k,nfluid+ndim*l  )
                end do


                vdhx = 0.0_dp ; vdhy = 0.0_dp ; vdhz = 0.0_dp
                do l = 1 , nrv

                   Ya = v (i,j,k, niv+l ) * rho_i
                   vdtx = vdcx * Ya ; vdty = vdcy * Ya ; vdtz = vdcz * Ya

                   wrk  = - dm ( i,j,k, l ) * thd % Wc (l) * W_i (i,j,k)

                   ! plus diffussion velocity: Yalpha * V alpha j
                   vdtx = vdtx + wrk * fd (i,j,k,nfluid+ndim*l-2)
                   vdty = vdty + wrk * fd (i,j,k,nfluid+ndim*l-1)
                   vdtz = vdtz + wrk * fd (i,j,k,nfluid+ndim*l  )

                   ! Visous term in rho*Yalpha equation: 
                   ! - M*\sqrt(gamma_inf)/(Re*Le*Pr) * rho * Yalpha * V alpha j
                   fd (i,j,k,nfluid+ndim*l-2) = rvd * vdtx ! Xa_x
                   fd (i,j,k,nfluid+ndim*l-1) = rvd * vdty ! Xa_y
                   fd (i,j,k,nfluid+ndim*l  ) = rvd * vdtz ! Xa_z

                   ! Enthalpy velocity diffusion term in energy equation
                   ! - M*\sqrt(gamma_inf)/(Re*Le*Pr) * rho * Yalpha * V alpha j * halpha
                   vdhx = vdhx + ha (i,j,k,l) * fd (i,j,k,nfluid+ndim*l-2)
                   vdhy = vdhy + ha (i,j,k,l) * fd (i,j,k,nfluid+ndim*l-1)
                   vdhz = vdhz + ha (i,j,k,l) * fd (i,j,k,nfluid+ndim*l  )

                end do


                ! All visous terms in rho*e_tot equation: - u_i * tau_ij + q_j
                rlam = - rlam0 * ct (i,j,k)

                q_x (i,j,k) = uxs  * sx_x (i,j,k)    + &
                              vys  * sx_y (i,j,k)    + &
                              wzs  * sx_z (i,j,k)    + &
                              rlam * fd   (i,j,k,10) + &
                              vdhx

                q_y (i,j,k) = uxs  * sx_y (i,j,k)    + &
                              vys  * sy_y (i,j,k)    + &
                              wzs  * sy_z (i,j,k)    + &
                              rlam * fd   (i,j,k,11) + &
                              vdhy

                q_z (i,j,k) = uxs  * sx_z (i,j,k)    + &
                              vys  * sy_z (i,j,k)    + &
                              wzs  * sz_z (i,j,k)    + &
                              rlam * fd   (i,j,k,12) + &
                              vdhz

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
                   mu (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) ! mu = -rlamcp
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l-2) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l-2)
                      fd (i,j,k,nfluid+ndim*l-1) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l-1)
                      fd (i,j,k,nfluid+ndim*l  ) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l  )
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( Xa , dm , mu , ct )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               wrk_y (sx:ex,sy:ey,sz:ez) , &
               wrk_z (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_xyz 3')


    ! viscous flux in rho*u equation -> fl (:,:,:,2)
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


    ! viscous flux in rho*v equation -> fl (:,:,:,3)
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


    ! viscous flux in rho*w equation -> fl (:,:,:,4)
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
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_xyz 4')

    ! viscous flux in rho*Y_alpha equation -> fl (:,:,:,niv+1:niv+npv)
    do l = 1 , nrv+npv
       call dx_fixed2 ( grid % dx_i                                                         , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-2) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                                 )
       call dy_fixed2 ( grid % dy_i                                                         , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) , &
                        wrka_y (sx:ex,sy:ey,sz:ez,l)                                 )
       call dz_fixed2 ( grid % dz_i                                                         , &
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


  end subroutine viscflux_xyz


!> \brief 3D viscous flux with LES models.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine viscflux_LES_xyz ( inp , thd , adi , grid , T , W_i , cp , ha , v , mu_SGS , fl )


    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid   !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T      !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp     !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha     !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v      !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mu_SGS !< turbulent viscosity (for LES)
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl     !< viscous flux


    integer (ip)                                  :: nfluid , nSGS
    integer (ip)                                  :: ok , domain_id , i , j , k , l
    integer (ip)                                  :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                     :: wrk , rmu , rmu_SGS , rvd , rvd0 , rlam , rlam0 , rlamcp0 , &
                                                     Sc_SGS_i , Pr_SGS_i , gm1_g , ct_SGS , dm_SGS , dmrw
    real (dp)                                     :: vdcx , vdtx , vdhx , dxadx , &
                                                     vdcy , vdty , vdhy , dxady , &
                                                     vdcz , vdtz , vdhz , dxadz
    real (dp)                                     :: rho_i , uxs , vys , wzs , Ya
    ! viscous variables
    real (dp) , dimension (:,:,:) , allocatable   :: mu , ct
    real (dp) , dimension (:,:,:,:) , allocatable :: dm
    ! SGS additional variables
    real (dp) , dimension (:,:,:) , allocatable   :: tau_iso_SGS
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable   :: ux , vy , wz
    real (dp) , dimension (:,:,:,:) , allocatable :: Xa
    ! additional first derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: q_x  , q_y  , q_z  , &
                                                     sx_x , sx_y , sx_z , &
                                                     sy_y , sy_z , sz_z
    real (dp) , dimension (:,:,:) , allocatable   :: sx_x_SGS , sx_y_SGS , sx_z_SGS , &
                                                     sy_y_SGS , sy_z_SGS , sz_z_SGS

    ! additional second derivatives
    real (dp) , dimension (:,:,:) , allocatable   :: wrk_x , wrk_y , wrk_z
    real (dp) , dimension (:,:,:,:) , allocatable :: wrka_x , wrka_y , wrka_z
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable :: fd

    ! BC_wallinjection variables
    real (dp)                       :: rinj , x0 , z0 , rrr , xpa , xpb

    ! Integration of the source_sink terms in SGS variance equations
    real (dp) , dimension (:,:,:,:) , allocatable :: grad2
    real (dp)                                     :: grad2max , tau , wrktau , prod , dissip
    ! Constants
    real (dp) , parameter         :: C_xi      = 1.0_dp          , &
                                     eps       = 1.0e-15_dp


    nfluid  = ndim * ( ndim + 1 )
    nSGS    = nfluid + ndim * ( nrv + npv + nvv )

    rvd0    = adi % sqgmr / adi % Sc     ! M * sqrt(gamma) / ( Re * Sc )
    rlam0   = adi % sqgmr * adi % ggmopr ! M * sqrt(gamma) * gamma / ( Re * (gamma-1) * Pr )
    rlamcp0 = adi % sqgmr / adi % Pr     ! M * sqrt(gamma) / ( Re * Pr )

    Sc_SGS_i = 1.0_dp / inp % Sc_SGS
    Pr_SGS_i = 1.0_dp / inp % Pr_SGS
    gm1_g    = adi % gm1 / adi % gamma

    allocate ( ux          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)              , &
               vy          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)              , &
               wz          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)              , &
               Xa          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv+nvv)  , &
               fd          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)       , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_xyz')


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
             end do
          end do
       end do

    end do


    ! passive scalar and variance equations
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do l = nrv+1,nrv+npv+nvv
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
                            fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,12) ) ! dT/dz = fd (:,:,:,12) nfluid = 12 in 3D


    ! 1st derivative for molar fraction, passive scalar and variance
    do l = 1 , nrv+npv+nvv
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
    fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) = 1.0_dp
!    call shock_det_slv ( v , T , W_i , fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) )
    call shock_det_ducros_slv ( 0 , fd , fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv) )


    deallocate ( ux , vy , wz )


    ! symmetry plane and adiabatic wall boundary conditions. For
    ! non-adiabatic and/or non-zero-shear and/or non-zero-species, you
    ! must specity the temperature and/or the velocity and/or the
    ! species normal gradients
    if ( bc (W) == symmetryplane .or. bc (W) == adiabaticwall ) then
       if ( bc (W) == symmetryplane ) then
          fd (sx,sy:ey,sz:ez,1) = 0.0_dp
          fd (sx,sy:ey,sz:ez,4) = 0.0_dp
          fd (sx,sy:ey,sz:ez,7) = 0.0_dp ! normal velocity gradient
       end if
       fd (sx,sy:ey,sz:ez,10) = 0.0_dp ! normal temperature gradient
       do l = 1 , nrv+npv+nvv
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
       do l = 1 , nrv+npv+nvv
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
       do l = 1 , nrv+npv+nvv
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
       do l = 1 , nrv+npv+nvv
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
       do l = 1 , nrv+npv+nvv
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
       do l = 1 , nrv+npv+nvv
          fd (sx:ex,sy:ey,ez,nfluid+ndim*l  ) = 0.0_dp ! normal species gradient
       end do
       fd (sx:ex,sy:ey,ez,nSGS+3) = 0.0_dp ! normal mixture molecular weight gradient
    end if

    if ( bc (N) == wallinjection ) then
       rinj = 0.0010_dp / adi % L_ref ! if it modified, changed also in module_BCs.f90 !!!
       x0   = 0.2000_dp / adi % L_ref
       z0   = 0.0000_dp / adi % L_ref
       do k = sz , ez
          do i = sx , ex
             rrr = sqrt ( (grid % x(i)-x0)*(grid % x(i)-x0) + (grid % z(k)-z0)*(grid % z(k)-z0) )
             if ( rrr > rinj ) then ! wall/slipwall
!                fd (i,ey,k,2) = 0.0_dp
!                fd (i,ey,k,5) = 0.0_dp
!                fd (i,ey,k,8) = 0.0_dp ! normal velocity gradient
                fd (i,ey,k,11) = 0.0_dp ! normal temperature gradient
                do l = 1 , nrv+npv+nvv
                   fd (i,ey,k,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
                end do
                fd (i,ey,k,nSGS+2)      = 0.0_dp ! normal mixture molecular weight gradient
             end if
          end do
       end do
    end if

    if ( bc (S) == wallinjection ) then
       rinj = 0.0010_dp / adi % L_ref ! if it modified, changed also in module_BCs.f90 !!!
       x0   = 0.2000_dp / adi % L_ref
       z0   = 0.0000_dp / adi % L_ref
       xpa  = 0.122_dp  / adi % L_ref
       xpb  = 0.122_dp  / adi % L_ref + rinj + rinj
       do k = sz , ez
          do i = sx , ex
             rrr = sqrt ( (grid % x(i)-x0)*(grid % x(i)-x0) + (grid % z(k)-z0)*(grid % z(k)-z0) )
             if ( rrr > rinj ) then !.and. ( x(i) < xpa .or. xpb < x(i) ) ) then ! wall/slipwall
!                fd (i,sy,k,2) = 0.0_dp
!                fd (i,sy,k,5) = 0.0_dp
!                fd (i,sy,k,8) = 0.0_dp ! normal velocity gradient
                fd (i,sy,k,11) = 0.0_dp ! normal temperature gradient
                do l = 1 , nrv+npv+nvv
                   fd (i,sy,k,nfluid+ndim*l-1) = 0.0_dp ! normal species gradient
                end do
                fd (i,sy,k,nSGS+2)      = 0.0_dp ! normal mixture molecular weight gradient
             end if
          end do
       end do
    end if


    ! communicate the first derivatives
    call comm_deriv (fd)


    if ( nvv > 0 ) then

       ! Allocate variable for SGS variances equations. (A.Techer)
       allocate ( grad2       (sx:ex, sy:ey, sz:ez, 1:npv) , & ! not need in ghost points +/-ng
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_xyz SGS variance variables')

       ! Squared scalar gradient
       do l = nrv+1 , nrv+npv
          do k = sz , ez
             do j = sy , ey
                do i = sx , ex
                   grad2 (i,j,k,l-nrv) = fd (i,j,k,nfluid+ndim*l-2) * fd (i,j,k,nfluid+ndim*l-2) &
                                       + fd (i,j,k,nfluid+ndim*l-1) * fd (i,j,k,nfluid+ndim*l-1) &
                                       + fd (i,j,k,nfluid+ndim*l  ) * fd (i,j,k,nfluid+ndim*l  )

                   grad2max = grid % dx_i (i) * grid % dx_i (i) &
                            + grid % dy_i (j) * grid % dy_i (j) &
                            + grid % dz_i (k) + grid % dz_i (k)

                   grad2 (i,j,k,l-nrv) = min ( grad2 (i,j,k,l-nrv) , grad2max )
                end do
             end do
          end do
       end do

    end if


    ! components of the viscous stresses
    allocate ( dm           (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv)  , &
               mu           (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               ct           (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               q_x          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               q_y          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               q_z          (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
               sx_x         (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
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
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_xyz 2')


    ! SGS: viscosity
    do domain_id = -ndim , ndim
       call mu_SGS_selector ( inp , adi , domain_id , grid , v , fd , mu_SGS )
    end do

    ! SGS: tau_iso_SGS derivatives need to be COMMUNICATED
    if ( inp % tau_iso_SGS_switch ) then ! consider the isotropic part of SGS viscous stress tensor

       do domain_id = -ndim , ndim
          call tau_iso_SGS_selector ( inp , domain_id , grid , v , fd , tau_iso_SGS )
       end do

    else ! neglects the isotropic part of SGS viscous stress tensor

       tau_iso_SGS (:,:,:) = 0.0_dp

    end if


    do domain_id = -ndim , ndim

       call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
       call prim_vis_vars ( domain_id , thd , v , W_i , T , Xa , dm , mu , ct )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)
                uxs   = v (i,j,k,2) * rho_i
                vys   = v (i,j,k,3) * rho_i
                wzs   = v (i,j,k,4) * rho_i

                ct_SGS = adi % Pr * mu_SGS (i,j,k) * cp (i,j,k) * Pr_SGS_i
                dm_SGS = adi % Sc * mu_SGS (i,j,k) * rho_i * Sc_SGS_i

                rmu          = - adi % sqgmr * mu (i,j,k)
                rmu_SGS      = - adi % sqgmr * mu_SGS (i,j,k) / rmu

                rvd = rvd0 * v (i,j,k,1)

                ! Viscous stress tensor: -M*\sqrt(gamma_inf)/Re * tau_ij
                sx_x (i,j,k) = twothirds * rmu * ( fd (i,j,k,1) + fd (i,j,k,1)   - &
                                                 ( fd (i,j,k,5) + fd (i,j,k,9) ) )
                sx_x_SGS (i,j,k) = sx_x (i,j,k) * ( 1.0_dp + rmu_SGS ) + onethirds * tau_iso_SGS (i,j,k)

                sx_y (i,j,k) = rmu * ( fd (i,j,k,4) + fd (i,j,k,2) )
                sx_y_SGS (i,j,k) = sx_y (i,j,k) * ( 1.0_dp + rmu_SGS )

                sx_z (i,j,k) = rmu * ( fd (i,j,k,7) + fd (i,j,k,3) )
                sx_z_SGS (i,j,k) = sx_z (i,j,k) * ( 1.0_dp + rmu_SGS )

                sy_y (i,j,k) = twothirds * rmu * ( fd (i,j,k,5) + fd (i,j,k,5)   - &
                                                 ( fd (i,j,k,1) + fd (i,j,k,9) ) )
                sy_y_SGS (i,j,k) = sy_y (i,j,k) * ( 1.0_dp + rmu_SGS ) + onethirds * tau_iso_SGS (i,j,k)

                sy_z (i,j,k) = rmu * ( fd (i,j,k,8) + fd (i,j,k,6) )
                sy_z_SGS (i,j,k) = sy_z (i,j,k) * ( 1.0_dp + rmu_SGS )

                sz_z (i,j,k) = twothirds * rmu * ( fd (i,j,k,9) + fd (i,j,k,9)   - &
                                                 ( fd (i,j,k,1) + fd (i,j,k,5) ) )
                sz_z_SGS (i,j,k) = sz_z (i,j,k) * ( 1.0_dp + rmu_SGS ) + onethirds * tau_iso_SGS (i,j,k)

                ! Correction velocity: sum_{beta=1,N} D_{beta m} W_beta/W \pd{Xbeta}{x_j}
                vdcx = 0.0_dp ; vdcy = 0.0_dp ; vdcz = 0.0_dp
                do l = 1 , nrv
                   dmrw = dm (i,j,k, l ) * thd % Wc (l) * W_i (i,j,k) ! No SGS term
                   vdcx = vdcx + dmrw * fd (i,j,k,nfluid+ndim*l-2)
                   vdcy = vdcy + dmrw * fd (i,j,k,nfluid+ndim*l-1)
                   vdcz = vdcz + dmrw * fd (i,j,k,nfluid+ndim*l  )
                end do

                vdhx = 0.0_dp ; vdhy = 0.0_dp ; vdhz = 0.0_dp
                do l = 1 , nrv

                   Ya = v (i,j,k, niv+l ) * rho_i
                   vdtx = vdcx * Ya ; vdty = vdcy * Ya ; vdtz = vdcz * Ya
                   dxadx = fd (i,j,k,nfluid+ndim*l-2)
                   dxady = fd (i,j,k,nfluid+ndim*l-1)
                   dxadz = fd (i,j,k,nfluid+ndim*l  )

                   dmrw = - dm ( i,j,k, l ) * thd % Wc (l) * W_i (i,j,k)

                   ! Enthalpy velocity diffusion term in energy equation (No SGS term)
                   ! - M*\sqrt(gamma_inf)/(Re*Le*Pr) * sum_{alpha=1,N} rho * Yalpha * V alpha j * halpha
                   vdhx = vdhx + ha (i,j,k,l) * rvd * ( vdtx + dmrw * dxadx )
                   vdhy = vdhy + ha (i,j,k,l) * rvd * ( vdty + dmrw * dxady )
                   vdhz = vdhz + ha (i,j,k,l) * rvd * ( vdtz + dmrw * dxadz )

                   wrk = dmrw - dm_SGS * thd % Wc (l) * W_i (i,j,k)

                   ! plus diffussion velocity
                   vdtx = vdtx + wrk * dxadx
                   vdty = vdty + wrk * dxady
                   vdtz = vdtz + wrk * dxadz

                   wrk = - dm_SGS * thd % Wc (l) * Xa (i,j,k,l)

                   ! plus an extra SGS term
                   vdtx = vdtx + wrk * fd (i,j,k,nSGS+1)
                   vdty = vdty + wrk * fd (i,j,k,nSGS+2)
                   vdtz = vdtz + wrk * fd (i,j,k,nSGS+3)

                   ! finally ... chi j = -cste * rho * Yalpha * V alpha j
                   fd (i,j,k,nfluid+ndim*l-2) = rvd * vdtx ! Xa_x
                   fd (i,j,k,nfluid+ndim*l-1) = rvd * vdty ! Xa_y
                   fd (i,j,k,nfluid+ndim*l  ) = rvd * vdtz ! Xa_z

                end do

                rlam = - rlam0 * ( ct (i,j,k) + ct_SGS )

                q_x (i,j,k) = uxs  * sx_x_SGS (i,j,k)       + &
                              vys  * sx_y_SGS (i,j,k)       + &
                              wzs  * sx_z_SGS (i,j,k)       + &
                              rlam * fd (i,j,k,10)          + &
                              vdhx

                q_y (i,j,k) = uxs  * sx_y_SGS (i,j,k)       + &
                              vys  * sy_y_SGS (i,j,k)       + &
                              wzs  * sy_z_SGS (i,j,k)       + &
                              rlam * fd (i,j,k,11)          + &
                              vdhy

                q_z (i,j,k) = uxs  * sx_z_SGS (i,j,k)       + &
                              vys  * sy_z_SGS (i,j,k)       + &
                              wzs  * sz_z_SGS (i,j,k)       + &
                              rlam * fd (i,j,k,12)          + &
                              vdhz

             end do
          end do
       end do

    end do

    deallocate ( tau_iso_SGS )
    deallocate ( sx_x , sx_y , sx_z , &
                        sy_y , sy_z , &
                               sz_z   )


    ! passive scalar and variance equations: viscous term
    if ( npv > 0 ) then

       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   dm_SGS = adi % Sc * mu_SGS (i,j,k) / v (i,j,k,1) * Sc_SGS_i
                   mu (i,j,k) = - rlamcp0 * ct (i,j,k) / cp (i,j,k) & ! > laminar/molecular term
                                - rvd0 * v (i,j,k,1) * dm_sgs         ! > SGS term for LES
                end do
             end do
          end do

          do l = nrv+1 , nrv+npv+nvv
             do k = k0 , k1
                do j = j0 , j1
                   do i = i0 , i1
                      fd (i,j,k,nfluid+ndim*l-2) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l-2)
                      fd (i,j,k,nfluid+ndim*l-1) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l-1)
                      fd (i,j,k,nfluid+ndim*l  ) = mu (i,j,k) * fd (i,j,k,nfluid+ndim*l  )
                   end do
                end do
             end do
          end do

       end do

    end if


    deallocate ( Xa , dm , mu )
! A.Techer : for integration of source_sink terms in SGS variance equations
!    deallocate ( ct )


    allocate ( wrk_x (sx:ex,sy:ey,sz:ez) , &
               wrk_y (sx:ex,sy:ey,sz:ez) , &
               wrk_z (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_xyz 3')


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


    allocate ( wrka_x (sx:ex,sy:ey,sz:ez,nrv+npv+nvv) , &
               wrka_y (sx:ex,sy:ey,sz:ez,nrv+npv+nvv) , &
               wrka_z (sx:ex,sy:ey,sz:ez,nrv+npv+nvv) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate viscflux_LES_xyz 4')


    ! viscous flux in rho*Y_alpha equation -> fl (:,:,:,niv+1:niv+npv+nvv)
    do l = 1 , nrv+npv+nvv
       call dx_fixed2 ( grid % dx_i                                                         , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-2) , &
                        wrka_x (sx:ex,sy:ey,sz:ez,l)                                 )
       call dy_fixed2 ( grid % dy_i                                                         , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) , &
                        wrka_y (sx:ex,sy:ey,sz:ez,l)                                 )
       call dz_fixed2 ( grid % dz_i                                                         , &
                        fd     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) , &
                        wrka_z (sx:ex,sy:ey,sz:ez,l)                                 )
    end do

    do l = 1 , nrv+npv+nvv
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


    ! Integration the source_sink terms in SGS variance equations (A.Techer)
    if ( nvv > 0 ) then

    do l = 1 , nvv
       do k = sz , ez
          do j = sy , ey
             do i = sx , ex

                rho_i = 1.0_dp / v (i,j,k,1)

                dm_SGS = adi % Sc * mu_SGS (i,j,k) * rho_i * Sc_SGS_i ! (-)

                ! Relaxation time ( + eps to avoid divisions by zero )
                tau = grid % delta (i,j,k) * grid % delta (i,j,k) / ( rvd0 * dm_SGS + eps ) ! (-)

                wrktau = - 2.0_dp * C_xi / tau ! independent to SGS variance => necessarily < 0

               !=======================================================
               ! SGS Variance equation
               !=======================================================

!                l = 6

                prod   = 0.0_dp
                dissip = 0.0_dp

                ! Production term
                prod = rvd0 * 2.0_dp * adi % Sc * Sc_sgs_i * mu_SGS (i,j,k) * grad2 (i,j,k,l)

                ! Dissipation term
                dissip = wrktau * v (i,j,k,niv+nrv+npv+l)

                fl (i,j,k,niv+nrv+npv+l) = fl (i,j,k,niv+nrv+npv+l) - ( prod + dissip )

!                !=======================================================
!                ! Squared filtered passive scalar equation \tilde{\xi}^2
!                !=======================================================

!                l = 8

!                dissip = 0.0_dp

!                ! Dissipation term
!                dissip = - 2.0_dp * ( rlamcp0 * ct (i,j,k) / cp (i,j,k) &
!                                    + rvd0 * adi % Sc * Sc_sgs_i * mu_SGS (i,j,k) ) * grad2 (i,j,k,l)

!                fl (i,j,k,niv+l) = fl (i,j,k,niv+l) - dissip

!                !=======================================================
!                ! Filtered squared passive scalar equation \tilde{\xi^2}
!                !=======================================================

!                l = 7

!                dissip = 0.0_dp

!                ! Dissipation term
!                dissip = - rlamcp0 * 2.0_dp * ct (i,j,k) / cp (i,j,k) * grad2 (i,j,k,l)

!                wrk = wrktau * ( v (i,j,k,niv+7) - ( v (i,j,k,niv+5) * v (i,j,k,niv+5) ) * rho_i )
!                dissip = dissip + wrk

!                fl (i,j,k,niv+l) = fl (i,j,k,niv+l) - dissip

!                !=======================================================
!                ! Departure from maximal variance equation \Delta_\xi = \xi*(1-\xi) - variance
!                !=======================================================

!                l = 9

!                prod = 0.0_dp

!                ! Production term
!                prod = rlamcp0 * 2.0_dp * ct (i,j,k) / cp (i,j,k) * grad2 (i,j,k,l)

!                wrk = - wrktau * ( v (i,j,k,niv+5) * ( 1.0_dp - v (i,j,k,niv+5) * rho_i ) - v (i,j,k,niv+9) )
!                prod = prod + wrk

!                fl (i,j,k,niv+l) = fl (i,j,k,niv+l) - prod

             end do
          end do
       end do
    end do

    ! Deallocate for SGS variance equations
    if ( allocated (grad2) ) deallocate (grad2)

    end if ! ( nvv > 0 )

    deallocate ( ct )


  end subroutine viscflux_LES_xyz


end module viscflux
