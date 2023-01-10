!------------------------------------------------------------------------------
! MODULE: sgs_models
!------------------------------------------------------------------------------
!> \brief Subgrid-scale (SGS) models.
!!
!! Calculates the dynamic _turbulent_ viscosity using:\n
!!   * the Smagorinsky model,\n
!!   * the structure function.\n
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module SGS_models


  use parameters
  use parallel
  use input
  use adim
  use deriv
  use tools , only : shock_det_slv , shock_det_ducros_slv , shock_det_ducros_post , filter_width

  implicit none

  real (dp) , parameter , private :: Cst_Kolmogorov = 1.4_dp   !< Kolmogorov constant
  real (dp) , parameter , private :: Cst_Yoshizawa = 6.6e-3_dp !< Yoshizawa constant


contains


!> \brief Selector for the viscosity SGS models
!!
!! Reads the keyword to select the corresponding viscosity SGS model.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine mu_SGS_selector ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< adi derived type
    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i      !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i      !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd        !< array of derivative variables
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< dynamic SGS viscosity


    if ( inp % mu_SGS_model == Smagorinsky ) then
       call mu_SGS_smagorinsky ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , mu_SGS )
    else if ( inp % mu_SGS_model == StructureFunction ) then
       call mu_SGS_structure_function ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , mu_SGS )
    else if ( inp % mu_SGS_model == WALE ) then
       call mu_SGS_WALE ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , mu_SGS )
    else
       call abort_mpi ('SGS model ' // trim (inp % mu_SGS_model) // ' not defined')
    end if


  end subroutine mu_SGS_selector


!> \brief Selector for the viscosity SGS models for post-processing with POSTizem
!!
!! Reads the keyword to select the corresponding viscosity SGS model.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine mu_SGS_selector_post ( inp , adi , dx_i , dy_i , dz_i , v , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< adi derived type
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i      !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i      !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< dynamic SGS viscosity

    real (dp) , dimension (:,:,:) , allocatable         :: ux , vy , wz
    real (dp)                                           :: rho_i
    integer (ip)                                        :: ok , i , j , k

    allocate ( ux       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate mu_SGS_selector_post'

    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             rho_i      = 1.0_dp / v (i,j,k,1)
             ux (i,j,k) = v (i,j,k,2) * rho_i
             vy (i,j,k) = v (i,j,k,3) * rho_i
             wz (i,j,k) = v (i,j,k,4) * rho_i
          end do
       end do
    end do


    if ( inp % mu_SGS_model == Smagorinsky ) then
       call mu_SGS_smagorinsky_post ( inp , adi , dx_i , dy_i , dz_i , v , ux , vy , wz , mu_SGS )
    else if ( inp % mu_SGS_model == StructureFunction ) then
       call mu_SGS_structure_function_post ( inp , adi , dx_i , dy_i , dz_i , v , ux , vy , wz , mu_SGS )
    else if ( inp % mu_SGS_model == WALE ) then
       call mu_SGS_WALE_post ( inp , adi , dx_i , dy_i , dz_i , v , ux , vy , wz , mu_SGS )
    else
       call abort_mpi ('SGS model_post ' // trim (inp % mu_SGS_model) // ' not defined')
    end if

    deallocate ( ux , vy , wz )


  end subroutine mu_SGS_selector_post


!! !> \brief Selector for the SGS tensor model
!!
!! Reads the keyword to select the corresponding SGS tensor model.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine tau_iso_SGS_selector ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , tau_iso_SGS )


    type (inp_type) , intent (in)                                  :: inp         !< input derived type
    type (adi_type) , intent (in)                                  :: adi         !< non-dimensional derived type
    integer (ip) , intent (in)                                     :: domain_id   !< subdomain selection
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i        !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i        !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i        !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v           !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd          !< array of derivative variables
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: tau_iso_SGS !< subgrid isotropic tensor


    if ( inp % tau_iso_SGS_model == Yoshizawa ) then
       call tau_iso_SGS_yoshizawa ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , tau_iso_SGS )
    else
       call abort_mpi ('SGS isotropic tensor model ' // trim (inp % tau_iso_SGS_model) // ' not defined')
    end if


  end subroutine tau_iso_SGS_selector


!> \brief Smagorinsky model
!!
!! Calculates the dynamic _turbulent_ viscosity using the Smagorinsky's model.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine mu_SGS_smagorinsky ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i      !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i      !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd        !< array of derivative variables
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< dynamic SGS viscosity


    integer (ip)                                      :: ok , i , j , k
    integer (ip)                                      :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                         :: Cs_2 ! Smagorinsky "square" constant
    real (dp)                                         :: Sbar
    real (dp) , allocatable , dimension (:,:,:)       :: delta_2
    real (dp)                                         :: s11 , s12 , s13 , &
                                                               s22 , s23 , &
                                                                     s33


    Cs_2 = inp % mu_SGS_factor
    Cs_2 = Cs_2 * Cs_2 / adi % sqgmr

    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    allocate  ( delta_2 ( i0:i1 , j0:j1 , k0:k1 )                   , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_smagorinsky')

    call filter_width ( domain_id , dx_i , dy_i , dz_i , delta_2 )
    delta_2 = delta_2 * delta_2


    if ( ndim == 1 ) then ! 1D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                s11 = fd (i,j,k,1)

                Sbar = s11*s11

                Sbar = sqrt ( Sbar + Sbar )

                mu_SGS (i,j,k) = v (i,j,k,1) * Cs_2 * delta_2 (i,j,k) * Sbar * fd (i,j,k,nderiv)

             end do
          end do
       end do


    else if ( ndim == 2 ) then ! 2D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                s11 = fd (i,j,k,1)
                s12 = 0.5_dp * ( fd (i,j,k,2) + fd (i,j,k,3) )

                s22 = fd (i,j,k,4)

                Sbar = s11*s11 + s12*s12 + &
                       s12*s12 + s22*s22

                Sbar = sqrt ( Sbar + Sbar )

                mu_SGS (i,j,k) = v (i,j,k,1) * Cs_2 * delta_2 (i,j,k) * Sbar * fd (i,j,k,nderiv)

             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                s11 = fd (i,j,k,1)
                s12 = 0.5_dp * ( fd (i,j,k,2) + fd (i,j,k,4) )
                s13 = 0.5_dp * ( fd (i,j,k,3) + fd (i,j,k,7) )

                s22 = fd (i,j,k,5)
                s23 = 0.5_dp * ( fd (i,j,k,6) + fd (i,j,k,8) )

                s33 = fd (i,j,k,9)

                Sbar = s11*s11 + s12*s12 + s13*s13 + &
                       s12*s12 + s22*s22 + s23*s23 + &
                       s13*s13 + s23*s23 + s33*s33

                Sbar = sqrt ( Sbar + Sbar )

                mu_SGS (i,j,k) = v (i,j,k,1) * Cs_2 * delta_2 (i,j,k) * Sbar * fd (i,j,k,nderiv)

             end do
          end do
       end do


    end if


    deallocate  ( delta_2 )


  end subroutine mu_SGS_smagorinsky


!> \brief Smagorinsky model (for post-processing)
!!
!! Calculates the dynamic _turbulent_ viscosity using the Smagorinsky's model.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine mu_SGS_smagorinsky_post ( inp , adi , dx_i , dy_i , dz_i , v , ux , vy , wz , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i      !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i      !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: ux        !< x-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: vy        !< y-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: wz        !< z-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< dynamic SGS viscosity


    integer (ip)                                      :: ok , i , j , k
    real (dp)                                         :: Cs_2 ! Smagorinsky "square" constant
    real (dp)                                         :: Sbar
    real (dp) , allocatable , dimension (:,:,:)       :: delta_2
    real (dp)                                         :: s11 , s12 , s13 , &
                                                               s22 , s23 , &
                                                                     s33
    real (dp) , dimension (:,:,:) , allocatable       :: shk

    real (dp) , dimension (:,:,:) , allocatable       :: dudx , dudy , dudz , &
                                                         dvdx , dvdy , dvdz , &
                                                         dwdx , dwdy , dwdz

    Cs_2 = inp % mu_SGS_factor
    Cs_2 = Cs_2 * Cs_2 / adi % sqgmr

    allocate  ( delta_2 ( sx:ex , sy:ey , sz:ez )                   , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_Smago_post 0')

    call filter_width ( 0 , dx_i , dy_i , dz_i , delta_2 )
    delta_2 = delta_2 * delta_2


    if ( ndim == 1 ) then ! 1D problem


       allocate ( dudx   (sx:ex,sy:ey,sz:ez) , &
                  shk    (sx:ex,sy:ey,sz:ez) , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_Smago_post 1')

       call shock_det_ducros_post ( dx_i , dx_i , dx_i , ux , ux , ux , shk )

       call comm_one (ux)

       call dx ( dx_i , ux , dudx )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                s11 = dudx (i,j,k)

                Sbar = s11*s11

                Sbar = sqrt ( Sbar + Sbar )

                mu_SGS (i,j,k) = v (i,j,k,1) * Cs_2 * delta_2 (i,j,k) * Sbar * shk (i,j,k)

             end do
          end do
       end do

       deallocate ( dudx )


    else if ( ndim == 2 ) then ! 2D problem


       allocate ( dudx   (sx:ex,sy:ey,sz:ez) , &
                  dudy   (sx:ex,sy:ey,sz:ez) , &
                  dvdx   (sx:ex,sy:ey,sz:ez) , &
                  dvdy   (sx:ex,sy:ey,sz:ez) , &
                  shk    (sx:ex,sy:ey,sz:ez) , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_Smago_post 2')

       call shock_det_ducros_post ( dx_i , dy_i , dy_i , ux , vy , vy , shk )

       call comm_one (ux) ; call comm_one (vy)

       call dx ( dx_i , ux , dudx )
       call dx ( dx_i , vy , dvdx )

       call dy ( dy_i , ux , dudy )
       call dy ( dy_i , vy , dvdy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                s11 = dudx (i,j,k)
                s12 = 0.5_dp * ( dudy (i,j,k) + dvdx (i,j,k) )

                s22 = dvdy (i,j,k)

                Sbar = s11*s11 + s12*s12 + &
                       s12*s12 + s22*s22

                Sbar = sqrt ( Sbar + Sbar )

                mu_SGS (i,j,k) = v (i,j,k,1) * Cs_2 * delta_2 (i,j,k) * Sbar * shk (i,j,k)

             end do
          end do
       end do

       deallocate ( dudx , dudy , &
                    dvdx , dvdy )


    else if ( ndim == 3 ) then ! 3D problem


       allocate ( dudx   (sx:ex,sy:ey,sz:ez) , &
                  dudy   (sx:ex,sy:ey,sz:ez) , &
                  dudz   (sx:ex,sy:ey,sz:ez) , &
                  dvdx   (sx:ex,sy:ey,sz:ez) , &
                  dvdy   (sx:ex,sy:ey,sz:ez) , &
                  dvdz   (sx:ex,sy:ey,sz:ez) , &
                  dwdx   (sx:ex,sy:ey,sz:ez) , &
                  dwdy   (sx:ex,sy:ey,sz:ez) , &
                  dwdz   (sx:ex,sy:ey,sz:ez) , &
                  shk    (sx:ex,sy:ey,sz:ez) , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_Smago_post 3')

       call shock_det_ducros_post ( dx_i , dy_i , dz_i , ux , vy , wz , shk )

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( dx_i , ux , dudx )
       call dx ( dx_i , vy , dvdx )
       call dx ( dx_i , wz , dwdx )

       call dy ( dy_i , ux , dudy )
       call dy ( dy_i , vy , dvdy )
       call dy ( dy_i , wz , dwdy )

       call dz ( dz_i , ux , dudz )
       call dz ( dz_i , vy , dvdz )
       call dz ( dz_i , wz , dwdz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                s11 = dudx (i,j,k)
                s12 = 0.5_dp * ( dudy (i,j,k) + dvdx (i,j,k) )
                s13 = 0.5_dp * ( dwdz (i,j,k) + dudz (i,j,k) )

                s22 = dvdy (i,j,k)
                s23 = 0.5_dp * ( dvdz (i,j,k) + dwdz (i,j,k) )

                s33 = dwdz (i,j,k)

                Sbar = s11*s11 + s12*s12 + s13*s13 + &
                       s12*s12 + s22*s22 + s23*s23 + &
                       s13*s13 + s23*s23 + s33*s33

                Sbar = sqrt ( Sbar + Sbar )

                mu_SGS (i,j,k) = v (i,j,k,1) * Cs_2 * delta_2 (i,j,k) * Sbar * shk (i,j,k)

             end do
          end do
       end do

       deallocate ( dudx , dudy , dudz , &
                    dvdx , dvdy , dvdz , &
                    dwdx , dwdy , dwdz )


    end if


    deallocate ( shk , delta_2 )


  end subroutine mu_SGS_smagorinsky_post


!> \brief Structure function model
!!
!! Calculates the dynamic _turbulent_ viscosity using the structure function.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine mu_SGS_structure_function ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    integer (ip)                                  , intent (in)    :: domain_id !< subdomain selection
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i      !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i      !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd        !< array of derivative variables
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< dynamic SGS viscosity


    integer (ip)                                                   :: ok , i , j , k 
    integer (ip)                                                   :: i0 , i1 , j0 , j1 , k0 , k1
    real    (dp) , parameter              :: onethird = 1.0_dp / 3.0_dp , power = -3.0_dp / 2.0_dp
    real    (dp)                                                   :: A_W , b_E , c_S , d_N , e_B , f_F
    real    (dp)                                                   :: cst , rho_i , flux
    real    (dp) , allocatable , dimension (:,:,:)                 :: delta



    cst = 0.105_dp * ( inp % mu_SGS_factor ) ** power
    cst = cst / adi % sqgmr


    ! allocate  ( delta ( i0:i1 , j0:j1 , k0:k1 )  ,  stat = ok )
    ! if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_structure_function')

    ! call filter_width ( domain_id , dx_i , dy_i , dz_i , delta )
    ! delta = delta * delta


    ! call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    ! if ( ndim == 1 ) then ! 1D problem


    !    do k = k0 , k1
    !       do j = j0 , j1
    !          do i = i0 , i1

    !             rho_i = 1.0_dp / v (i,j,k,1)

    !             a_W = 1.0_dp ; if ( i==1   ) a_W = 0.0_dp
    !             b_E = 1.0_dp ; if ( i==ntx ) b_E = 0.0_dp

    !             flux = 0.0_dp

    !             flux = flux +                                                                &
    !                  ( a_W * (v (i,j,k,2) - v (i-1,j,k,2)) * (v (i,j,k,2) - v (i-1,j,k,2)) + &
    !                    b_E * (v (i,j,k,2) - v (i+1,j,k,2)) * (v (i,j,k,2) - v (i+1,j,k,2))   )

    !             flux = flux * rho_i * rho_i * 0.5_dp

    !             mu_SGS (i,j,k) = v (i,j,k,1) * cst * delta (i,j,k) * sqrt (flux) * fd (i,j,k,nderiv)

    !          end do
    !       end do
    !    end do


    ! else if ( ndim == 2 ) then ! 2D problem


    !    do k = k0 , k1
    !       do j = j0 , j1
    !          do i = i0 , i1

    !             rho_i = 1.0_dp / v (i,j,k,1)

    !             a_W = 1.0_dp ; if ( i==1   ) a_W = 0.0_dp
    !             b_E = 1.0_dp ; if ( i==ntx ) b_E = 0.0_dp
    !             c_S = 1.0_dp ; if ( j==1   ) c_S = 0.0_dp
    !             d_N = 1.0_dp ; if ( j==nty ) d_N = 0.0_dp

    !             flux = 0.0_dp

    !             flux = flux +                                                                &
    !                  ( a_W * (v (i,j,k,2) - v (i-1,j,k,2)) * (v (i,j,k,2) - v (i-1,j,k,2)) + &
    !                    b_E * (v (i,j,k,2) - v (i+1,j,k,2)) * (v (i,j,k,2) - v (i+1,j,k,2)) + &
    !                    c_S * (v (i,j,k,2) - v (i,j-1,k,2)) * (v (i,j,k,2) - v (i,j-1,k,2)) + &
    !                    d_N * (v (i,j,k,2) - v (i,j+1,k,2)) * (v (i,j,k,2) - v (i,j+1,k,2))   )

    !             flux = flux +                                                                &
    !                  ( a_W * (v (i,j,k,3) - v (i-1,j,k,3)) * (v (i,j,k,3) - v (i-1,j,k,3)) + &
    !                    b_E * (v (i,j,k,3) - v (i+1,j,k,3)) * (v (i,j,k,3) - v (i+1,j,k,3)) + &
    !                    c_S * (v (i,j,k,3) - v (i,j-1,k,3)) * (v (i,j,k,3) - v (i,j-1,k,3)) + &
    !                    d_N * (v (i,j,k,3) - v (i,j+1,k,3)) * (v (i,j,k,3) - v (i,j+1,k,3))   )

    !             flux = flux * rho_i * rho_i * 0.25_dp

    !             mu_SGS (i,j,k) = v (i,j,k,1) * cst * delta (i,j,k) * sqrt (flux) * fd (i,j,k,nderiv)

    !          end do
    !       end do
    !    end do


    ! else if ( ndim == 3 ) then ! 3D problem


    !    do k = k0 , k1
    !       do j = j0 , j1
    !          do i = i0 , i1

    !             rho_i = 1.0_dp / v (i,j,k,1)

    !             a_W = 1.0_dp ; if ( i==1   ) a_W = 0.0_dp
    !             b_E = 1.0_dp ; if ( i==ntx ) b_E = 0.0_dp
    !             c_S = 1.0_dp ; if ( j==1   ) c_S = 0.0_dp
    !             d_N = 1.0_dp ; if ( j==nty ) d_N = 0.0_dp
    !             e_B = 1.0_dp ; if ( k==1   ) e_B = 0.0_dp
    !             f_F = 1.0_dp ; if ( k==ntz ) f_F = 0.0_dp

    !             flux = 0.0_dp

    !             flux = flux +                                                                &
    !                  ( a_W * (v (i,j,k,2) - v (i-1,j,k,2)) * (v (i,j,k,2) - v (i-1,j,k,2)) + &
    !                    b_E * (v (i,j,k,2) - v (i+1,j,k,2)) * (v (i,j,k,2) - v (i+1,j,k,2)) + &
    !                    c_S * (v (i,j,k,2) - v (i,j-1,k,2)) * (v (i,j,k,2) - v (i,j-1,k,2)) + &
    !                    d_N * (v (i,j,k,2) - v (i,j+1,k,2)) * (v (i,j,k,2) - v (i,j+1,k,2)) + &
    !                    e_B * (v (i,j,k,2) - v (i,j,k-1,2)) * (v (i,j,k,2) - v (i,j,k-1,2)) + &
    !                    f_F * (v (i,j,k,2) - v (i,j,k+1,2)) * (v (i,j,k,2) - v (i,j,k+1,2))   )

    !             flux = flux +                                                                &
    !                  ( a_W * (v (i,j,k,3) - v (i-1,j,k,3)) * (v (i,j,k,3) - v (i-1,j,k,3)) + &
    !                    b_E * (v (i,j,k,3) - v (i+1,j,k,3)) * (v (i,j,k,3) - v (i+1,j,k,3)) + &
    !                    c_S * (v (i,j,k,3) - v (i,j-1,k,3)) * (v (i,j,k,3) - v (i,j-1,k,3)) + &
    !                    d_N * (v (i,j,k,3) - v (i,j+1,k,3)) * (v (i,j,k,3) - v (i,j+1,k,3)) + &
    !                    e_B * (v (i,j,k,3) - v (i,j,k-1,3)) * (v (i,j,k,3) - v (i,j,k-1,3)) + &
    !                    f_F * (v (i,j,k,3) - v (i,j,k+1,3)) * (v (i,j,k,3) - v (i,j,k+1,3))   )

    !             flux = flux +                                                                &
    !                  ( a_W * (v (i,j,k,4) - v (i-1,j,k,4)) * (v (i,j,k,4) - v (i-1,j,k,4)) + &
    !                    b_E * (v (i,j,k,4) - v (i+1,j,k,4)) * (v (i,j,k,4) - v (i+1,j,k,4)) + &
    !                    c_S * (v (i,j,k,4) - v (i,j-1,k,4)) * (v (i,j,k,4) - v (i,j-1,k,4)) + &
    !                    d_N * (v (i,j,k,4) - v (i,j+1,k,4)) * (v (i,j,k,4) - v (i,j+1,k,4)) + &
    !                    e_B * (v (i,j,k,4) - v (i,j,k-1,4)) * (v (i,j,k,4) - v (i,j,k-1,4)) + &
    !                    f_F * (v (i,j,k,4) - v (i,j,k+1,4)) * (v (i,j,k,4) - v (i,j,k+1,4))   )

    !             flux = flux * rho_i * rho_i / 6.0_dp

    !             mu_SGS (i,j,k) = v (i,j,k,1) * cst * delta (i,j,k) * sqrt (flux) * fd (i,j,k,nderiv)

    !          end do
    !       end do
    !    end do


    ! end if


    deallocate ( delta )


  end subroutine mu_SGS_structure_function


!> \brief Structure function model (for post-processing)
!!
!! Calculates the dynamic _turbulent_ viscosity using the structure function.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine mu_SGS_structure_function_post ( inp , adi , dx_i , dy_i , dz_i , v , ux , vy , wz , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i   !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i   !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i   !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v      !< conserved variables array
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: ux     !< x-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: vy     !< y-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: wz     !< z-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS !< dynamic SGS viscosity


    integer (ip)                                                   :: ok , i , j , k
    integer (ip)                                                   :: i0 , i1 , j0 , j1 , k0 , k1
    real    (dp) , parameter                                       :: onethird = 1.0_dp / 3.0_dp , power = -3.0_dp / 2.0_dp
    real    (dp)                                                   :: A_W , b_E , c_S , d_N , e_B , f_F
    real    (dp)                                                   :: cst , rho_i , flux
    real    (dp) , dimension (:,:,:) , allocatable                 :: shk , delta


    allocate ( shk   (sx:ex,sy:ey,sz:ez) , &
               delta (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate mu_SGS_structure_function_post'

    call shock_det_ducros_post ( dx_i , dy_i , dz_i , ux , vy , wz , shk )

    call filter_width ( 0 , dx_i , dy_i , dz_i , delta )
    delta = delta * delta


    cst = 0.105_dp * ( inp % mu_SGS_factor ) ** power
    cst = cst / adi % sqgmr


    call domain_select ( 0 , i0 , i1 , j0 , j1 , k0 , k1 )


    if ( ndim == 1 ) then ! 1D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)

                a_W = 1.0_dp ; if ( i==1   ) a_W = 0.0_dp
                b_E = 1.0_dp ; if ( i==ntx ) b_E = 0.0_dp

                flux = 0.0_dp

                flux = flux +                                                                &
                     ( a_W * (v (i,j,k,2) - v (i-1,j,k,2)) * (v (i,j,k,2) - v (i-1,j,k,2)) + &
                       b_E * (v (i,j,k,2) - v (i+1,j,k,2)) * (v (i,j,k,2) - v (i+1,j,k,2))   )

                flux = flux * rho_i * rho_i * 0.5_dp

                mu_SGS (i,j,k) = v (i,j,k,1) * cst * delta (i,j,k) * sqrt (flux) * shk (i,j,k)

             end do
          end do
       end do


    else if ( ndim == 2 ) then ! 2D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)

                a_W = 1.0_dp ; if ( i==1   ) a_W = 0.0_dp
                b_E = 1.0_dp ; if ( i==ntx ) b_E = 0.0_dp
                c_S = 1.0_dp ; if ( j==1   ) c_S = 0.0_dp
                d_N = 1.0_dp ; if ( j==nty ) d_N = 0.0_dp

                flux = 0.0_dp

                flux = flux +                                                                &
                     ( a_W * (v (i,j,k,2) - v (i-1,j,k,2)) * (v (i,j,k,2) - v (i-1,j,k,2)) + &
                       b_E * (v (i,j,k,2) - v (i+1,j,k,2)) * (v (i,j,k,2) - v (i+1,j,k,2)) + &
                       c_S * (v (i,j,k,2) - v (i,j-1,k,2)) * (v (i,j,k,2) - v (i,j-1,k,2)) + &
                       d_N * (v (i,j,k,2) - v (i,j+1,k,2)) * (v (i,j,k,2) - v (i,j+1,k,2))   )

                flux = flux +                                                                &
                     ( a_W * (v (i,j,k,3) - v (i-1,j,k,3)) * (v (i,j,k,3) - v (i-1,j,k,3)) + &
                       b_E * (v (i,j,k,3) - v (i+1,j,k,3)) * (v (i,j,k,3) - v (i+1,j,k,3)) + &
                       c_S * (v (i,j,k,3) - v (i,j-1,k,3)) * (v (i,j,k,3) - v (i,j-1,k,3)) + &
                       d_N * (v (i,j,k,3) - v (i,j+1,k,3)) * (v (i,j,k,3) - v (i,j+1,k,3))   )

                flux = flux * rho_i * rho_i * 0.25_dp

                mu_SGS (i,j,k) = v (i,j,k,1) * cst * delta (i,j,k) * sqrt (flux) * shk (i,j,k)

             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho_i = 1.0_dp / v (i,j,k,1)

                a_W = 1.0_dp ; if ( i==1   ) a_W = 0.0_dp
                b_E = 1.0_dp ; if ( i==ntx ) b_E = 0.0_dp
                c_S = 1.0_dp ; if ( j==1   ) c_S = 0.0_dp
                d_N = 1.0_dp ; if ( j==nty ) d_N = 0.0_dp
                e_B = 1.0_dp ; if ( k==1   ) e_B = 0.0_dp
                f_F = 1.0_dp ; if ( k==ntz ) f_F = 0.0_dp

                flux = 0.0_dp

                flux = flux +                                                                &
                     ( a_W * (v (i,j,k,2) - v (i-1,j,k,2)) * (v (i,j,k,2) - v (i-1,j,k,2)) + &
                       b_E * (v (i,j,k,2) - v (i+1,j,k,2)) * (v (i,j,k,2) - v (i+1,j,k,2)) + &
                       c_S * (v (i,j,k,2) - v (i,j-1,k,2)) * (v (i,j,k,2) - v (i,j-1,k,2)) + &
                       d_N * (v (i,j,k,2) - v (i,j+1,k,2)) * (v (i,j,k,2) - v (i,j+1,k,2)) + &
                       e_B * (v (i,j,k,2) - v (i,j,k-1,2)) * (v (i,j,k,2) - v (i,j,k-1,2)) + &
                       f_F * (v (i,j,k,2) - v (i,j,k+1,2)) * (v (i,j,k,2) - v (i,j,k+1,2))   )

                flux = flux +                                                                &
                     ( a_W * (v (i,j,k,3) - v (i-1,j,k,3)) * (v (i,j,k,3) - v (i-1,j,k,3)) + &
                       b_E * (v (i,j,k,3) - v (i+1,j,k,3)) * (v (i,j,k,3) - v (i+1,j,k,3)) + &
                       c_S * (v (i,j,k,3) - v (i,j-1,k,3)) * (v (i,j,k,3) - v (i,j-1,k,3)) + &
                       d_N * (v (i,j,k,3) - v (i,j+1,k,3)) * (v (i,j,k,3) - v (i,j+1,k,3)) + &
                       e_B * (v (i,j,k,3) - v (i,j,k-1,3)) * (v (i,j,k,3) - v (i,j,k-1,3)) + &
                       f_F * (v (i,j,k,3) - v (i,j,k+1,3)) * (v (i,j,k,3) - v (i,j,k+1,3))   )

                flux = flux +                                                                &
                     ( a_W * (v (i,j,k,4) - v (i-1,j,k,4)) * (v (i,j,k,4) - v (i-1,j,k,4)) + &
                       b_E * (v (i,j,k,4) - v (i+1,j,k,4)) * (v (i,j,k,4) - v (i+1,j,k,4)) + &
                       c_S * (v (i,j,k,4) - v (i,j-1,k,4)) * (v (i,j,k,4) - v (i,j-1,k,4)) + &
                       d_N * (v (i,j,k,4) - v (i,j+1,k,4)) * (v (i,j,k,4) - v (i,j+1,k,4)) + &
                       e_B * (v (i,j,k,4) - v (i,j,k-1,4)) * (v (i,j,k,4) - v (i,j,k-1,4)) + &
                       f_F * (v (i,j,k,4) - v (i,j,k+1,4)) * (v (i,j,k,4) - v (i,j,k+1,4))   )

                flux = flux * rho_i * rho_i / 6.0_dp

                mu_SGS (i,j,k) = v (i,j,k,1) * cst * delta (i,j,k) * sqrt (flux) * shk (i,j,k)

             end do
          end do
       end do


    end if


    deallocate (shk)


  end subroutine mu_SGS_structure_function_post


!> \brief WALE model
!!
!! Calculates the dynamic _turbulent_ viscosity using the WALE (Wall Adaptative Local Eddy) model.
!!
!! References:
!! -# F. Nicoud, F. Ducros. _Subgrid-scale stress modelling 
!!    based on the square of the velocity gradient tensor_
!!    Flow, Turbulence and Combustion, (1999)
!!
!! -# E. Garnier et al., Large Eddy Simulation for Compressible Flows,
!!    Scientific Computation, p. 88, (2009).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco


  subroutine mu_SGS_WALE ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i      !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i      !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd        !< array of derivative variables
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< turbulent dynamic viscosity


    integer (ip)                                      :: ok , i , j , k
    integer (ip)                                      :: i0 , i1 , j0 , j1 , k0 , k1
    integer (ip)                                      :: is , js , ls , ms
    real (dp)                                         :: Cw_2 ! WALE "square" constant
    real (dp) , parameter                             :: onethird = 1.0_dp / 3.0_dp
    real (dp)                                         :: Sbar , Sd , Sra , sij , sijd
    real (dp) , dimension (:,:)   , allocatable       :: dudx , kdelta
    real (dp) , dimension (:,:,:) , allocatable       :: delta_2


    Cw_2 = inp % mu_SGS_factor * sqrt ( 10.6_dp )
    Cw_2 = Cw_2 * Cw_2 / adi % sqgmr


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    allocate ( dudx    ( ndim , ndim )             , &
               kdelta  ( ndim , ndim )             , &
               delta_2 ( i0:i1 , j0:j1 , k0:k1 )   , &
               stat = ok )

    if ( ok > 0 ) stop 'error allocate mu_SGS_WALE'


    ! Kronecker delta Dij
    kdelta = 0.0_dp
    do i = 1 , ndim
       kdelta (i,i) = 1.0_dp
    end do

    ! Square filter-width
    call filter_width ( domain_id , dx_i , dy_i , dz_i , delta_2 )
    delta_2 = delta_2 * delta_2


    do k = k0 , k1
      do j = j0 , j1
         do i = i0 , i1

            ! Definition of dudx
            ls = 1

            do is = 1 , ndim
               do js = 1 , ndim
                  dudx (is,js) = fd (i,j,k,ls)
                  ls = ls + 1
               end do
            end do

            ! Initialization of Sbar , Sd , Sra
            Sbar   = 0.0_dp               ! Sij
            Sd     = 0.0_dp               ! Sijd
            Sra    = 0.0_dp               ! Sratio

            do is = 1 , ndim
               do js = 1 , ndim

                  ! Sij = 1/2 * ( dudx(i,j) + dudx(j,i) )
                  sij  = 0.5_dp * ( dudx (is,js) + dudx (js,is) )
                  Sbar = Sbar + sij * sij

                  ! Traceless symmetric part of the square of the velocity gradient tensor
                  !   Sij^d = 1/2 * ( dudx(i,l) * dudx(l,j) + dudx(j,l) * dudx(l,i) )
                  !         - 1/3 * ( dudx(m,l) * dudx(l,m) * kdelta(i,j) )
                  do ls = 1 , ndim
                     do ms = 1 , ndim

                        sijd = 0.5_dp * ( dudx (is,ls) * dudx (ls,js) + dudx (js,ls) * dudx (ls,is) ) - &
                               onethird * dudx (ms,ls) * dudx (ls,ms) * kdelta(is,js)
                        Sd   = Sd + sijd * sijd

                     end do
                  end do

               end do
            end do

            ! Turbulent inverse time scale
            !   Sra = (Sijd * Sijd)^3/2 / [ (Sij * Sij)^5/2 + (Sijd * Sijd)^5/4 ]
            Sra = ( Sbar**2.5_dp + Sd**1.25_dp )

            if ( Sra > 0.0_dp ) then 
               Sra = Sd**1.5_dp / Sra
            else
               Sra = 0.0_dp
            end if

            mu_SGS (i,j,k) = v (i,j,k,1) * Cw_2 * delta_2 (i,j,k) * Sra * fd (i,j,k,nderiv)

         end do
      end do
    end do


    deallocate ( dudx , kdelta , delta_2 )


  end subroutine mu_SGS_WALE


!> \brief WALE model (for post-processing)
!!
!! Calculates the dynamic _turbulent_ viscosity using the WALE (Wall Adaptative Local Eddy) model.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco


  subroutine mu_SGS_WALE_post ( inp , adi , dx_i , dy_i , dz_i , v , ux , vy , wz , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i      !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i      !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: ux        !< x-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: vy        !< y-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: wz        !< z-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< turbulent dynamic viscosity


    integer (ip)                                      :: ok , i , j , k
    integer (ip)                                      :: is , js , ls , ms
    real (dp)                                         :: Cw_2 ! WALE "square" constant
    real (dp) , parameter                             :: onethird  = 1.0_dp / 3.0_dp
    real (dp)                                         :: Sbar , Sd , Sra , sij , sijd
    real (dp) , dimension (:,:)   , allocatable       :: duidxj , kdelta
    real (dp) , dimension (:,:,:) , allocatable       :: delta_2

    real (dp) , dimension (:,:,:) , allocatable       :: dudx , dudy , dudz , &
                                                         dvdx , dvdy , dvdz , &
                                                         dwdx , dwdy , dwdz

    real (dp) , dimension (:,:,:) , allocatable       :: shk


    Cw_2 = inp % mu_SGS_factor * sqrt ( 10.6_dp )
    Cw_2 = Cw_2 * Cw_2 / adi % sqgmr

    allocate ( duidxj  (ndim,ndim)         , &
               kdelta  (ndim,ndim)         , &
               delta_2 (sx:ex,sy:ey,sz:ez) , &
               shk     (sx:ex,sy:ey,sz:ez) , &
               dudx    (sx:ex,sy:ey,sz:ez) , &
               dudy    (sx:ex,sy:ey,sz:ez) , &
               dudz    (sx:ex,sy:ey,sz:ez) , &
               dvdx    (sx:ex,sy:ey,sz:ez) , &
               dvdy    (sx:ex,sy:ey,sz:ez) , &
               dvdz    (sx:ex,sy:ey,sz:ez) , &
               dwdx    (sx:ex,sy:ey,sz:ez) , &
               dwdy    (sx:ex,sy:ey,sz:ez) , &
               dwdz    (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) then
       write (*,*) 'error allocate mu_SGS_WALE_post'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if


    call shock_det_ducros_post ( dx_i , dy_i , dz_i , ux , vy , wz , shk )

    ! Kronecker delta Dij
    kdelta = 0.0_dp
    do i = 1 , ndim
       kdelta (i,i) = 1.0_dp
    end do

    ! Square filter-width & velocity components 1st derivative
    call filter_width ( 0 , dx_i , dy_i , dz_i , delta_2 )
    delta_2 = delta_2 * delta_2

    if      ( ndim == 1 ) then ! 1D problem

       call comm_one (ux)

       call dx ( dx_i , ux , dudx )

    else if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux) ; call comm_one (vy)

       call dx ( dx_i , ux , dudx )
       call dx ( dx_i , vy , dvdx )

       call dy ( dy_i , ux , dudy )
       call dy ( dy_i , vy , dvdy )

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( dx_i , ux , dudx )
       call dx ( dx_i , vy , dvdx )
       call dx ( dx_i , wz , dwdx )

       call dy ( dy_i , ux , dudy )
       call dy ( dy_i , vy , dvdy )
       call dy ( dy_i , wz , dwdy )

       call dz ( dz_i , ux , dudz )
       call dz ( dz_i , vy , dvdz )
       call dz ( dz_i , wz , dwdz )

    end if


    do k = sz,ez
       do j = sy,ey
          do i = sx,ex

            ! Definition of duidxj
            if      ( ndim == 1 ) then

               duidxj (1,1) = dudx (i,j,k)

            else if ( ndim == 2 ) then

               duidxj (1,1) = dudx (i,j,k)
               duidxj (2,1) = dvdx (i,j,k)

               duidxj (1,2) = dudy (i,j,k)
               duidxj (2,2) = dvdy (i,j,k)

            else if ( ndim == 3 ) then

               duidxj (1,1) = dudx (i,j,k)
               duidxj (2,1) = dvdx (i,j,k)
               duidxj (3,1) = dwdx (i,j,k)

               duidxj (1,2) = dudy (i,j,k)
               duidxj (2,2) = dvdy (i,j,k)
               duidxj (3,2) = dwdy (i,j,k)

               duidxj (1,3) = dudz (i,j,k)
               duidxj (2,3) = dvdz (i,j,k)
               duidxj (3,3) = dwdz (i,j,k)

            end if

            ! Initialization of Sbar , Sd , Sra
            Sbar   = 0.0_dp               ! Sij
            Sd     = 0.0_dp               ! Sijd
            Sra    = 0.0_dp               ! Sratio

            do is = 1 , ndim
               do js = 1 , ndim

                  ! Sij = 1/2 * ( dudx(i,j) + dudx(j,i) )
                  sij  = 0.5_dp * ( duidxj (is,js) + duidxj (js,is) )
                  Sbar = Sbar + sij * sij

                  ! Traceless symmetric part of the square of the velocity gradient tensor
                  !   Sij^d = 1/2 * ( dudx(i,l) * dudx(l,j) + dudx(j,l) * dudx(l,i) )
                  !         - 1/3 * ( dudx(m,l) * dudx(l,m) * kdelta(i,j) )
                  do ls = 1 , ndim
                     do ms = 1 , ndim

                        sijd = 0.5_dp * ( duidxj (is,ls) * duidxj (ls,js) + duidxj (js,ls) * duidxj (ls,is) ) - &
                               onethird * duidxj (ms,ls) * duidxj (ls,ms) * kdelta(is,js)
                        Sd   = Sd + sijd * sijd

                     end do
                  end do

               end do
            end do

            ! Turbulent inverse time scale
            !   Sra = (Sijd * Sijd)^3/2 / [ (Sij * Sij)^5/2 + (Sijd * Sijd)^5/4 ]
            Sra = ( Sbar**2.5_dp + Sd**1.25_dp )

            if ( Sra > 0.0_dp ) then 
               Sra = Sd**1.5_dp / Sra
            else
               Sra = 0.0_dp
            end if

            mu_SGS (i,j,k) = v (i,j,k,1) * Cw_2 * delta_2 (i,j,k) * Sra * shk (i,j,k)

         end do
      end do
    end do


    deallocate ( duidxj , kdelta , delta_2 , shk )
    deallocate ( dudx , dudy , dudz , &
                 dvdx , dvdy , dvdz , &
                 dwdx , dwdy , dwdz )


  end subroutine mu_SGS_WALE_post


!> \brief Yoshizawa (1986) isotropic tensor model
!!
!! Calculates the isotropic tensor (turbulent pressure) using the Yoshizawa's model.
!!
!! References: -# A. Yoshizawa, ”Statistical theory for compressible
!! turbulent shear flows, with the application to subgrid modeling”,
!! Phys. Fluids A29, 2152 (1986).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine tau_iso_SGS_yoshizawa ( inp , adi , domain_id , dx_i , dy_i , dz_i , v , fd , tau_iso_SGS )


    type (inp_type) , intent (in)                                  :: inp         !< input derived type
    type (adi_type) , intent (in)                                  :: adi         !< non-dimensional derived type
    integer (ip)                                  , intent (in)    :: domain_id   !< subdomain selection
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i        !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i        !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i        !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v           !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd          !< strain rate array
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: tau_iso_SGS !< isotropic tensor array



    integer (ip)                                        :: ok , i , j , k
    integer (ip)                                        :: i0 , i1 , j0 , j1 , k0 , k1
    real    (dp)                                        :: cst , Sbar2
    real    (dp) , allocatable , dimension (:,:,:)      :: delta_2
    real    (dp)                                        :: s11 , s12 , s13 , &
                                                                 s22 , s23 , &
                                                                       s33


    cst = inp % tau_iso_SGS_factor


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    allocate  ( delta_2 ( i0:i1 , j0:j1 , k0:k1 )                   , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_smagorinsky')

    call filter_width ( domain_id , dx_i , dy_i , dz_i , delta_2 )
    delta_2 = delta_2 * delta_2


    if ( ndim == 1 ) then ! 1D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                s11 = fd (i,j,k,1)

                Sbar2 = s11*s11

                Sbar2 = Sbar2 + Sbar2

                tau_iso_SGS (i,j,k) = v (i,j,k,1) * cst * delta_2 (i,j,k) * Sbar2
                tau_iso_SGS (i,j,k) = ( tau_iso_SGS (i,j,k) + tau_iso_SGS (i,j,k) ) * fd (i,j,k,nderiv)

             end do
          end do
       end do


    else if ( ndim == 2 ) then ! 2D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                s11 = fd (i,j,k,1)
                s12 = 0.5_dp * ( fd (i,j,k,2) + fd (i,j,k,3) )

                s22 = fd (i,j,k,4)

                Sbar2 = s11*s11 + s12*s12 + &
                        s12*s12 + s22*s22

                Sbar2 = Sbar2 + Sbar2

                tau_iso_SGS (i,j,k) = v (i,j,k,1) * cst * delta_2 (i,j,k) * Sbar2
                tau_iso_SGS (i,j,k) = ( tau_iso_SGS (i,j,k) + tau_iso_SGS (i,j,k) ) * fd (i,j,k,nderiv)

             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                s11 = fd (i,j,k,1)
                s12 = 0.5_dp * ( fd (i,j,k,2) + fd (i,j,k,4) )
                s13 = 0.5_dp * ( fd (i,j,k,3) + fd (i,j,k,7) )

                s22 = fd (i,j,k,5)
                s23 = 0.5_dp * ( fd (i,j,k,6) + fd (i,j,k,8) )

                s33 = fd (i,j,k,9)

                Sbar2 = s11*s11 + s12*s12 + s13*s13 + &
                        s12*s12 + s22*s22 + s23*s23 + &
                        s13*s13 + s23*s23 + s33*s33

                Sbar2 = Sbar2 + Sbar2

                tau_iso_SGS (i,j,k) = v (i,j,k,1) * cst * delta_2 (i,j,k) * Sbar2
                tau_iso_SGS (i,j,k) = ( tau_iso_SGS (i,j,k) + tau_iso_SGS (i,j,k) ) * fd (i,j,k,nderiv)

             end do
          end do
       end do


    end if


    deallocate  ( delta_2 )


  end subroutine tau_iso_SGS_yoshizawa


end module SGS_models
