!------------------------------------------------------------------------------
! MODULE: sgs_models
!------------------------------------------------------------------------------
!> \brief Subgrid-scale (SGS) models.
!!
!! Calculates the dynamic _turbulent_ viscosity using:\n
!!   * the Smagorinsky model,\n
!!   * the structure function.\n
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module SGS_models


  use parameters
  use parallel
  use input
  use adim
  use deriv
  use tools , only : shock_det_slv , shock_det_ducros_slv , shock_det_ducros_post

  implicit none

  real (dp) , parameter , private :: Cst_Kolmogorov = 1.4_dp   !< Kolmogorov constant
  real (dp) , parameter , private :: Cst_Yoshizawa = 6.6e-3_dp !< Yoshizawa constant for SGS mu
  real (dp) , parameter , private :: Cm_Yoshizawa = 0.069_dp   !< Yoshizawa constant for SGS tke
  real (dp) , parameter , private :: epsi = 1.0e-5_dp          !< to avoid divisions by zero



contains


!> \brief Selector for the viscosity SGS models
!!
!! Reads the keyword to select the corresponding viscosity SGS model.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine mu_SGS_selector ( inp , adi , domain_id , grid , v , fd , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< adi derived type
    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (inp_grid) , intent (in)                                  :: grid      !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd        !< array of derivative variables
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< dynamic SGS viscosity


    if ( inp % mu_SGS_model == Smagorinsky ) then
       call mu_SGS_smagorinsky ( inp , adi , domain_id , grid , v , fd , mu_SGS )
    else if ( inp % mu_SGS_model == StructureFunction ) then
       call mu_SGS_structure_function ( inp , adi , domain_id , grid , v , fd , mu_SGS )
    else if ( inp % mu_SGS_model == WALE ) then
       call mu_SGS_WALE ( inp , adi , domain_id , grid , v , fd , mu_SGS )
    else
       call abort_mpi ('SGS model ' // trim (inp % mu_SGS_model) // ' not defined')
    end if


  end subroutine mu_SGS_selector


!> \brief Selector for the viscosity SGS models for post-processing with POSTCOSINESS
!!
!! Reads the keyword to select the corresponding viscosity SGS model.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine mu_SGS_selector_post ( inp , adi , grid , v , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< adi derived type
    type (inp_grid) , intent (in)                                  :: grid      !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< dynamic SGS viscosity

    real (dp) , dimension (:,:,:) , allocatable         :: ux , vy , wz
    real (dp)                                           :: rho_i
    integer (ip)                                        :: ok , i , j , k

    allocate ( ux       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_selector_post')

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
       call mu_SGS_smagorinsky_post ( inp , adi , grid , v , ux , vy , wz , mu_SGS )
    else if ( inp % mu_SGS_model == StructureFunction ) then
       call mu_SGS_structure_function_post ( inp , adi , grid , v , ux , vy , wz , mu_SGS )
    else if ( inp % mu_SGS_model == WALE ) then
       call mu_SGS_WALE_post ( inp , adi , grid , v , ux , vy , wz , mu_SGS )
    else
       call abort_mpi ('SGS model_post ' // trim (inp % mu_SGS_model) // ' not defined')
    end if

    deallocate ( ux , vy , wz )


  end subroutine mu_SGS_selector_post


!! !> \brief Selector for the SGS tensor model
!!
!! Reads the keyword to select the corresponding SGS tensor model.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tau_iso_SGS_selector ( inp , domain_id , grid , v , fd , tau_iso_SGS )


    type (inp_type) , intent (in)                                  :: inp         !< input derived type
    integer (ip) , intent (in)                                     :: domain_id   !< subdomain selection
    type (inp_grid) , intent (in)                                  :: grid        !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v           !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd          !< array of derivative variables
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: tau_iso_SGS !< subgrid isotropic tensor


    if ( inp % tau_iso_SGS_model == Yoshizawa ) then
       call tau_iso_SGS_yoshizawa ( inp , domain_id , grid , v , fd , tau_iso_SGS )
    else
       call abort_mpi ('SGS isotropic tensor model ' // trim (inp % tau_iso_SGS_model) // ' not defined')
    end if


  end subroutine tau_iso_SGS_selector


!! !> \brief Selector for the SGS tensor model for post-processing with POSTCOSINESS
!!
!! Reads the keyword to select the corresponding SGS tensor model.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tau_iso_SGS_selector_post ( inp , grid , v , tau_iso_SGS )


    type (inp_type) , intent (in)                                  :: inp         !< input derived type
    type (inp_grid) , intent (in)                                  :: grid        !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v           !< conserved variables array
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: tau_iso_SGS !< subgrid isotropic tensor

    real (dp) , dimension (:,:,:) , allocatable         :: ux , vy , wz
    real (dp)                                           :: rho_i
    integer (ip)                                        :: ok , i , j , k

    allocate ( ux       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate tau_iso_SGS_selector_post')

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


    if ( inp % tau_iso_SGS_model == Yoshizawa ) then
       call tau_iso_SGS_yoshizawa_post ( inp , grid , v , ux , vy , wz , tau_iso_SGS )
    else
       call abort_mpi ('SGS isotropic tensor model ' // trim (inp % tau_iso_SGS_model) // ' not defined')
    end if

    deallocate ( ux , vy , wz )


  end subroutine tau_iso_SGS_selector_post


!> \brief Smagorinsky model
!!
!! Calculates the dynamic _turbulent_ viscosity using the Smagorinsky's model.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine mu_SGS_smagorinsky ( inp , adi , domain_id , grid , v , fd , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (inp_grid) , intent (in)                                  :: grid      !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd        !< array of derivative variables and shock sensor
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

    delta_2 (i0:i1,j0:j1,k0:k1) = grid % delta (i0:i1,j0:j1,k0:k1)
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
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine mu_SGS_smagorinsky_post ( inp , adi , grid , v , ux , vy , wz , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid      !< grid derived type
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

    delta_2 (sx:ex,sy:ey,sz:ez) = grid % delta (sx:ex,sy:ey,sz:ez)
    delta_2 = delta_2 * delta_2


    if ( ndim == 1 ) then ! 1D problem


       allocate ( dudx   (sx:ex,sy:ey,sz:ez) , &
                  shk    (sx:ex,sy:ey,sz:ez) , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_Smago_post 1')

       call shock_det_ducros_post ( grid % dx_i , grid % dy_i , grid % dz_i , ux , vy , wz , shk )

       call comm_one (ux)

       call dx ( grid % dx_i , ux , dudx )

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

       call shock_det_ducros_post ( grid % dx_i , grid % dy_i , grid % dz_i , ux , vy , wz , shk )

       call comm_one (ux) ; call comm_one (vy)

       call dx ( grid % dx_i , ux , dudx )
       call dx ( grid % dx_i , vy , dvdx )

       call dy ( grid % dy_i , ux , dudy )
       call dy ( grid % dy_i , vy , dvdy )

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

       call shock_det_ducros_post ( grid % dx_i , grid % dy_i , grid % dz_i , ux , vy , wz , shk )

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( grid % dx_i , ux , dudx )
       call dx ( grid % dx_i , vy , dvdx )
       call dx ( grid % dx_i , wz , dwdx )

       call dy ( grid % dy_i , ux , dudy )
       call dy ( grid % dy_i , vy , dvdy )
       call dy ( grid % dy_i , wz , dwdy )

       call dz ( grid % dz_i , ux , dudz )
       call dz ( grid % dz_i , vy , dvdz )
       call dz ( grid % dz_i , wz , dwdz )

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
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine mu_SGS_structure_function ( inp , adi , domain_id , grid , v , fd , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    integer (ip)                                  , intent (in)    :: domain_id !< subdomain selection
    type (inp_grid) , intent (in)                                  :: grid      !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd        !< array of derivative variables and shock sensor
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< dynamic SGS viscosity


    integer (ip)                                                   :: i , j , k
    integer (ip)                                                   :: i0 , i1 , j0 , j1 , k0 , k1
    real    (dp) , parameter                                       :: onethird = 1.0_dp / 3.0_dp , power = -3.0_dp / 2.0_dp
    real    (dp)                                                   :: A_W , b_E , c_S , d_N , e_B , f_F
    real    (dp)                                                   :: cst , rho_i , flux



    cst = 0.105_dp * ( inp % mu_SGS_factor ) ** power
    cst = cst / adi % sqgmr


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


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

                mu_SGS (i,j,k) = v (i,j,k,1) * cst * grid % delta (i,j,k) * sqrt (flux) * fd (i,j,k,nderiv)

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

                mu_SGS (i,j,k) = v (i,j,k,1) * cst * grid % delta (i,j,k) * sqrt (flux) * fd (i,j,k,nderiv)

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

                mu_SGS (i,j,k) = v (i,j,k,1) * cst * grid % delta (i,j,k) * sqrt (flux) * fd (i,j,k,nderiv)

             end do
          end do
       end do


    end if


  end subroutine mu_SGS_structure_function


!> \brief Structure function model (for post-processing)
!!
!! Calculates the dynamic _turbulent_ viscosity using the structure function.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine mu_SGS_structure_function_post ( inp , adi , grid , v , ux , vy , wz , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid   !< grid derived type
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
    real    (dp) , dimension (:,:,:) , allocatable                 :: shk


    allocate ( shk   (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_structure_function_post')

    call shock_det_ducros_post ( grid % dx_i , grid % dy_i , grid % dz_i , ux , vy , wz , shk )

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

                mu_SGS (i,j,k) = v (i,j,k,1) * cst * grid % delta (i,j,k) * sqrt (flux) * shk (i,j,k)

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

                mu_SGS (i,j,k) = v (i,j,k,1) * cst * grid % delta (i,j,k) * sqrt (flux) * shk (i,j,k)

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

                mu_SGS (i,j,k) = v (i,j,k,1) * cst * grid % delta (i,j,k) * sqrt (flux) * shk (i,j,k)

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
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine mu_SGS_WALE ( inp , adi , domain_id , grid , v , fd , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (inp_grid) , intent (in)                                  :: grid      !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd        !< array of derivative variables and shock sensor
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< turbulent dynamic viscosity


    integer (ip)                                      :: ok , i , j , k
    integer (ip)                                      :: i0 , i1 , j0 , j1 , k0 , k1
    integer (ip)                                      :: is , js , ls , ms
    real (dp)                                         :: Cw_2 ! WALE "square" constant
    real (dp) , parameter                             :: onethird = 1.0_dp / 3.0_dp
    real (dp)                                         :: Sbar , Sd , Sra , sij , sijd , sum_gmlglm
    real (dp) , dimension (:,:,:) , allocatable       :: delta_2


    Cw_2 = inp % mu_SGS_factor * sqrt ( 10.6_dp )
    Cw_2 = Cw_2 * Cw_2 / adi % sqgmr


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    allocate ( delta_2 ( i0:i1 , j0:j1 , k0:k1 )   , &
               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_WALE')


    ! Kronecker delta Dij
    kdelta = 0.0_dp
    do i = 1 , ndim
       kdelta (i,i) = 1.0_dp
    end do

    ! Square filter-width
    delta_2 (i0:i1,j0:j1,k0:k1) = grid % delta (i0:i1,j0:j1,k0:k1)
    delta_2 = delta_2 * delta_2


    do k = k0 , k1
      do j = j0 , j1
         do i = i0 , i1

            ! Definition of duidxj
            ls = 1

            do is = 1 , ndim
               do js = 1 , ndim
                  duidxj (is,js) = fd (i,j,k,ls)
                  ls = ls + 1
               end do
            end do

            ! Initialization of Sbar , Sd , Sra
            Sbar   = 0.0_dp               ! Sij
            Sd     = 0.0_dp               ! Sijd
            Sra    = 0.0_dp               ! Sratio

            ! sum_gmlglm = duidxj(m,l) * duidxj(l,m)
            sum_gmlglm = 0
            do ls = 1 , ndim
               do ms = 1 , ndim
                  sum_gmlglm = sum_gmlglm + duidxj (ls,ms) * duidxj (ms,ls)
               end do
            end do

            do is = 1 , ndim
               do js = 1 , ndim

                  ! Sij = 1/2 * ( duidxj(i,j) + duidxj(j,i) )
                  sij  = 0.5_dp * ( duidxj (is,js) + duidxj (js,is) )
                  Sbar = Sbar + sij * sij

                  ! Traceless symmetric part of the square of the velocity gradient tensor
                  !   Sij^d = 1/2 * ( duidxj(i,l) * duidxj(l,j) + duidxj(j,l) * duidxj(l,i) )
                  !         - 1/3 * sum_gmlglm * kdelta(i,j)
                  sijd = 0.0_dp ! Initialization of sijd
                  do ls = 1 , ndim
                     sijd = sijd + 0.5_dp * ( duidxj (is,ls) * duidxj (ls,js) + duidxj (js,ls) * duidxj (ls,is) )
                  end do
                  sijd = sijd - onethird * sum_gmlglm * kdelta(is,js)
                  Sd   = Sd + sijd * sijd

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


    deallocate ( delta_2 )


  end subroutine mu_SGS_WALE


!> \brief WALE model (for post-processing)
!!
!! Calculates the dynamic _turbulent_ viscosity using the WALE (Wall Adaptative Local Eddy) model.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine mu_SGS_WALE_post ( inp , adi , grid , v , ux , vy , wz , mu_SGS )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (adi_type) , intent (in)                                  :: adi       !< non-dimensional derived type
    type (inp_grid) , intent (in)                                  :: grid      !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v         !< conserved variables array
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: ux        !< x-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: vy        !< y-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: wz        !< z-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: mu_SGS    !< turbulent dynamic viscosity


    integer (ip)                                      :: ok , i , j , k
    integer (ip)                                      :: is , js , ls , ms
    real (dp)                                         :: Cw_2 ! WALE "square" constant
    real (dp) , parameter                             :: onethird = 1.0_dp / 3.0_dp
    real (dp)                                         :: Sbar , Sd , Sra , sij , sijd , sum_gmlglm
    real (dp) , dimension (:,:,:) , allocatable       :: delta_2

    real (dp) , dimension (:,:,:) , allocatable       :: dudx , dudy , dudz , &
                                                         dvdx , dvdy , dvdz , &
                                                         dwdx , dwdy , dwdz

    real (dp) , dimension (:,:,:) , allocatable       :: shk


    Cw_2 = inp % mu_SGS_factor * sqrt ( 10.6_dp )
    Cw_2 = Cw_2 * Cw_2 / adi % sqgmr

    allocate ( delta_2 (sx:ex,sy:ey,sz:ez) , &
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

    if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS_WALE_post')

    call shock_det_ducros_post ( grid % dx_i , grid % dy_i , grid % dz_i , ux , vy , wz , shk )

    ! Square filter-width & velocity components 1st derivative
    delta_2 (sx:ex,sy:ey,sz:ez) = grid % delta (sx:ex,sy:ey,sz:ez)
    delta_2 = delta_2 * delta_2

    if      ( ndim == 1 ) then ! 1D problem

       call comm_one (ux)

       call dx ( grid % dx_i , ux , dudx )

    else if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux) ; call comm_one (vy)

       call dx ( grid % dx_i , ux , dudx )
       call dx ( grid % dx_i , vy , dvdx )

       call dy ( grid % dy_i , ux , dudy )
       call dy ( grid % dy_i , vy , dvdy )

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( grid % dx_i , ux , dudx )
       call dx ( grid % dx_i , vy , dvdx )
       call dx ( grid % dx_i , wz , dwdx )

       call dy ( grid % dy_i , ux , dudy )
       call dy ( grid % dy_i , vy , dvdy )
       call dy ( grid % dy_i , wz , dwdy )

       call dz ( grid % dz_i , ux , dudz )
       call dz ( grid % dz_i , vy , dvdz )
       call dz ( grid % dz_i , wz , dwdz )

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

            ! sum_gmlglm = duidxj(m,l) * dudx(l,m)
            sum_gmlglm = 0
            do ls = 1 , ndim
               do ms = 1 , ndim
                  sum_gmlglm = sum_gmlglm + duidxj (ls,ms) * duidxj (ms,ls)
               end do
            end do

            do is = 1 , ndim
               do js = 1 , ndim

                  ! Sij = 1/2 * ( dudx(i,j) + dudx(j,i) )
                  sij  = 0.5_dp * ( duidxj (is,js) + duidxj (js,is) )
                  Sbar = Sbar + sij * sij

                  ! Traceless symmetric part of the square of the velocity gradient tensor
                  !   Sij^d = 1/2 * ( dudx(i,l) * dudx(l,j) + dudx(j,l) * dudx(l,i) )
                  !         - 1/3 * sum_gmlglm * kdelta(i,j)
                  sijd = 0.0_dp ! Initialization of sijd
                  do ls = 1 , ndim
                     sijd = sijd + 0.5_dp * ( duidxj (is,ls) * duidxj (ls,js) + duidxj (js,ls) * duidxj (ls,is) )
                  end do
                  sijd = sijd - onethird * sum_gmlglm * kdelta(is,js)
                  Sd   = Sd + sijd * sijd

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


    deallocate ( delta_2 , shk )
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
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tau_iso_SGS_yoshizawa ( inp , domain_id , grid , v , fd , tau_iso_SGS )


    type (inp_type) , intent (in)                                  :: inp         !< input derived type
    integer (ip)                                  , intent (in)    :: domain_id   !< subdomain selection
    type (inp_grid) , intent (in)                                  :: grid        !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v           !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: fd          !< strain rate array and shock sensor
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: tau_iso_SGS !< isotropic tensor array



    integer (ip)                                        :: ok , i , j , k
    integer (ip)                                        :: i0 , i1 , j0 , j1 , k0 , k1
    real    (dp)                                        :: cst , Sbar2
    real    (dp) , allocatable , dimension (:,:,:)      :: delta_2
    real    (dp)                                        :: s11 , s12 , s13 , &
                                                                 s22 , s23 , &
                                                                       s33


    cst = inp % tau_iso_SGS_factor
    cst = cst + cst

    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    allocate  ( delta_2 ( i0:i1 , j0:j1 , k0:k1 )                   , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate tau_iso_SGS_yoshizawa')

    delta_2 (i0:i1,j0:j1,k0:k1) = grid % delta (i0:i1,j0:j1,k0:k1)
    delta_2 = delta_2 * delta_2


    if ( ndim == 1 ) then ! 1D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                s11 = fd (i,j,k,1)

                Sbar2 = s11*s11

                Sbar2 = Sbar2 + Sbar2

                tau_iso_SGS (i,j,k) = v (i,j,k,1) * cst * delta_2 (i,j,k) * Sbar2 * fd (i,j,k,nderiv)

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

                tau_iso_SGS (i,j,k) = v (i,j,k,1) * cst * delta_2 (i,j,k) * Sbar2 * fd (i,j,k,nderiv)

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

                tau_iso_SGS (i,j,k) = v (i,j,k,1) * cst * delta_2 (i,j,k) * Sbar2 * fd (i,j,k,nderiv)

             end do
          end do
       end do


    end if


    deallocate  ( delta_2 )


  end subroutine tau_iso_SGS_yoshizawa


!> \brief Yoshizawa (1986) isotropic tensor model (for post-processing)
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tau_iso_SGS_yoshizawa_post ( inp , grid , v , ux , vy , wz , tau_iso_SGS )


    type (inp_type) , intent (in)                                  :: inp         !< input derived type
    type (inp_grid) , intent (in)                                  :: grid        !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)    :: v           !< conserved variables array
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: ux          !< x-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: vy          !< y-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: wz          !< z-component of the velocity
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: tau_iso_SGS !< isotropic tensor array



    integer (ip)                                        :: ok , i , j , k
    real    (dp)                                        :: cst , Sbar2

    integer (ip)                                      :: is , js
    real (dp)                                         :: Sbar , sij
    real (dp) , dimension (:,:,:) , allocatable       :: delta_2

    real (dp) , dimension (:,:,:) , allocatable       :: dudx , dudy , dudz , &
                                                         dvdx , dvdy , dvdz , &
                                                         dwdx , dwdy , dwdz

    real (dp) , dimension (:,:,:) , allocatable       :: shk


    cst = inp % tau_iso_SGS_factor
    cst = cst + cst

    allocate ( delta_2 (sx:ex,sy:ey,sz:ez) , &
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

    if ( ok > 0 ) call abort_mpi ('error allocate tau_iso_SGS_yoshizawa_post')

    call shock_det_ducros_post ( grid % dx_i , grid % dy_i , grid % dz_i , ux , vy , wz , shk )

    ! Square filter-width & velocity components 1st derivative
    delta_2 (sx:ex,sy:ey,sz:ez) = grid % delta (sx:ex,sy:ey,sz:ez)
    delta_2 = delta_2 * delta_2


    if      ( ndim == 1 ) then ! 1D problem

       call comm_one (ux)

       call dx ( grid % dx_i , ux , dudx )

    else if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux) ; call comm_one (vy)

       call dx ( grid % dx_i , ux , dudx )
       call dx ( grid % dx_i , vy , dvdx )

       call dy ( grid % dy_i , ux , dudy )
       call dy ( grid % dy_i , vy , dvdy )

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( grid % dx_i , ux , dudx )
       call dx ( grid % dx_i , vy , dvdx )
       call dx ( grid % dx_i , wz , dwdx )

       call dy ( grid % dy_i , ux , dudy )
       call dy ( grid % dy_i , vy , dvdy )
       call dy ( grid % dy_i , wz , dwdy )

       call dz ( grid % dz_i , ux , dudz )
       call dz ( grid % dz_i , vy , dvdz )
       call dz ( grid % dz_i , wz , dwdz )

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

            ! Initialization of Sbar
            Sbar   = 0.0_dp               ! Sij
            do is = 1 , ndim
               do js = 1 , ndim

                  ! Sij = 1/2 * ( dudx(i,j) + dudx(j,i) )
                  sij  = 0.5_dp * ( duidxj (is,js) + duidxj (js,is) )
                  Sbar = Sbar + sij * sij

               end do
            end do

            Sbar2 = Sbar + Sbar

            tau_iso_SGS (i,j,k) = v (i,j,k,1) * cst * delta_2 (i,j,k) * Sbar2 * shk (i,j,k)

         end do
      end do
    end do


    deallocate ( delta_2 , shk )
    deallocate ( dudx , dudy , dudz , &
                 dvdx , dvdy , dvdz , &
                 dwdx , dwdy , dwdz )


  end subroutine tau_iso_SGS_yoshizawa_post


!> \brief SGS turbulent kinetic energy model
!!
!! Calculates the subgrid-scale turbulent kinetic energy using the 
!! Yoshizawa model.
!!
!! References: 
!! -# A. Yoshizawa, ”Statistical theory for compressible
!!    turbulent shear flows, with the application to subgrid modeling”,
!!    Phys. Fluids A29, 2152 (1986).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tke_SGS_yoshizawa ( adi , mu_SGS , rho , delta , tke_SGS )


    type (adi_type) , intent (in)                                  :: adi     !< non-dimensional derived type
    real (dp) , dimension (:,:,:)   , allocatable , intent (in)    :: mu_SGS  !< dynamic SGS viscosity
    real (dp) , dimension (:,:,:)   , allocatable , intent (in)    :: rho     !< density
    real (dp) , dimension (:,:,:)   , allocatable , intent (in)    :: delta   !< filter width = mesh size
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: tke_SGS !< SGS turbulent kinetic energy

    integer (ip)                                        :: i , j , k

    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             tke_SGS (i,j,k) = adi % sqgmr * mu_SGS (i,j,k) &
                           / ( Cm_Yoshizawa * rho (i,j,k) * delta (i,j,k) )
             tke_SGS (i,j,k) = tke_SGS (i,j,k) * tke_SGS (i,j,k)
          end do
       end do
    end do


  end subroutine tke_SGS_yoshizawa


!> \brief Equilibrium SGS variance
!!
!! Calculates the subgrid-scale variance using the equilibrium assumption
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine equilibrium_SGS_variance ( inp , grid , Zm , eqSGS_var )


    type (inp_type) , intent (in)                                  :: inp       !< input derived type
    type (inp_grid) , intent (in)                                  :: grid      !< grid derived type
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: Zm        !< numerical mixture fraction
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: eqSGS_var !< equilibrium SGS variance

    integer (ip)                                        :: ok , i , j , k

    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3
    real (dp)                                           :: Zvmin , Zvmax
    real (dp) , parameter                               :: C_xi_i = 1.0_dp / 1.0_dp


    call comm_one (Zm)

    allocate ( wrk1 (sx:ex,sy:ey,sz:ez) , &
               wrk2 (sx:ex,sy:ey,sz:ez) , &
               wrk3 (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate equilibrium_SGS_variance')

    ! Zm gradient
    wrk1 = 0.0_dp ; wrk2 = 0.0_dp ; wrk3 = 0.0_dp
    if ( ndim >= 1 ) call dx ( grid % dx_i , Zm , wrk1 )
    if ( ndim >= 2 ) call dy ( grid % dy_i , Zm , wrk2 )
    if ( ndim == 3 ) call dz ( grid % dz_i , Zm , wrk3 )

    ! Equilibrium SGS variance
    do k = sz , ez
       do j = sy , ey
          do i = sx , ex

             eqSGS_var (i,j,k) = wrk1 (i,j,k) * wrk1 (i,j,k) &
                               + wrk2 (i,j,k) * wrk2 (i,j,k) &
                               + wrk3 (i,j,k) * wrk3 (i,j,k)
             eqSGS_var (i,j,k) = grid % delta (i,j,k) * grid % delta (i,j,k) &
                               * C_xi_i * eqSGS_var (i,j,k)

             if ( inp % forcesum ) then
                Zvmin             = 0.0_dp
                Zvmax             = Zm (i,j,k) * ( 1.0_dp - Zm (i,j,k) )
                eqSGS_var (i,j,k) = min ( max ( eqSGS_var (i,j,k) , Zvmin ) , Zvmax )
             end if

          end do
       end do
    end do


    deallocate ( wrk1 , wrk2 , wrk3 )


  end subroutine equilibrium_SGS_variance


!> \brief Segregation rate
!!
!! Calculates the segregation rate = ratio between the subgrid-scale variance 
!! and the subgrid-scale variance maxima.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine segregation_rate ( Zm , SGS_var , seg )


    real (dp) , dimension (:,:,:)   , allocatable , intent (in)    :: Zm        !< numerical mixture fraction
    real (dp) , dimension (:,:,:)   , allocatable , intent (in)    :: SGS_var   !< SGS variance
    real (dp) , dimension (:,:,:)   , allocatable , intent (inout) :: seg       !< segregation rate

    integer (ip)                                        :: i , j , k

    real (dp) , parameter                               :: Smin = 0.0_dp , Smax = 1.0_dp


    do k = sz , ez
       do j = sy , ey
          do i = sx , ex

             seg (i,j,k) = ( Zm (i,j,k) + epsi ) * ( 1.0_dp - Zm (i,j,k) - epsi )
             seg (i,j,k) = SGS_var (i,j,k) / seg (i,j,k)

             seg (i,j,k) = max ( min ( seg (i,j,k) , Smax ) , Smin )

          end do
       end do
    end do


  end subroutine segregation_rate


!> \brief Calculate SGS viscous stress tensor.
!!
!! \tau_{ij}^\sgs
!!  = \bar{\rho u_i u_j} - \bar{\rho} \tilde{u}_i \tilde{u}_j
!!  = - 2 \mu_\sgs \left( \tilde{S}_{ij} - \frac13 \tilde{S}_{kk}\delta_{ij} \right)
!!    + \frac13 \tau_{kk}^\sgs \delta_{ij}
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tau_ij_SGS ( inp , adi , grid , v , ux , vy , wz , mu_SGS , tau_sgs )


    type (inp_type) , intent (in)                                        :: inp     !< input derived type
    type (adi_type) , intent (in)                                        :: adi     !< non-dimensional derived type
    type (inp_grid) , intent (in)                                        :: grid    !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v       !< conserved variables array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux      !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy      !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz      !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: mu_SGS  !< SGS viscosity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (inout)     :: tau_SGS !< Reynolds stress tensor


    integer (ip)                                        :: ok , i , j , k
    real (dp) , parameter                               :: onethird = 1.0_dp / 3.0_dp
    real (dp)                                           :: divu , wrk , cst

    integer (ip)                                      :: is , js
    real (dp)                                         :: sij

    real (dp) , dimension (:,:,:) , allocatable       :: dudx , dudy , dudz , &
                                                         dvdx , dvdy , dvdz , &
                                                         dwdx , dwdy , dwdz , &
                                                         tau_iso_SGS

    allocate ( dudx    (sx:ex,sy:ey,sz:ez) , &
               dudy    (sx:ex,sy:ey,sz:ez) , &
               dudz    (sx:ex,sy:ey,sz:ez) , &
               dvdx    (sx:ex,sy:ey,sz:ez) , &
               dvdy    (sx:ex,sy:ey,sz:ez) , &
               dvdz    (sx:ex,sy:ey,sz:ez) , &
               dwdx    (sx:ex,sy:ey,sz:ez) , &
               dwdy    (sx:ex,sy:ey,sz:ez) , &
               dwdz    (sx:ex,sy:ey,sz:ez) , &
               tau_iso_SGS (sx:ex,sy:ey,sz:ez) , &
               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate tau_ij_SGS')


    if      ( ndim == 1 ) then ! 1D problem

       call comm_one (ux)

       call dx ( grid % dx_i , ux , dudx )

    else if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux) ; call comm_one (vy)

       call dx ( grid % dx_i , ux , dudx )
       call dx ( grid % dx_i , vy , dvdx )

       call dy ( grid % dy_i , ux , dudy )
       call dy ( grid % dy_i , vy , dvdy )

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( grid % dx_i , ux , dudx )
       call dx ( grid % dx_i , vy , dvdx )
       call dx ( grid % dx_i , wz , dwdx )

       call dy ( grid % dy_i , ux , dudy )
       call dy ( grid % dy_i , vy , dvdy )
       call dy ( grid % dy_i , wz , dwdy )

       call dz ( grid % dz_i , ux , dudz )
       call dz ( grid % dz_i , vy , dvdz )
       call dz ( grid % dz_i , wz , dwdz )

    end if


    ! Isotropic part of the SGS viscous stress tensor
    if ( inp % tau_iso_SGS_switch ) then ! consider it
       call tau_iso_SGS_selector_post ( inp , grid , v , tau_iso_SGS )
    else ! neglect it
       tau_iso_SGS (:,:,:) = 0.0_dp
    end if

    cst = - adi % sqgmr * 2.0_dp

!    do k = sz,ez
!       do j = sy,ey
!          do i = sx,ex

!            ! Definition of duidxj
!            if      ( ndim == 1 ) then

!               duidxj (1,1) = dudx (i,j,k)

!            else if ( ndim == 2 ) then

!               duidxj (1,1) = dudx (i,j,k)
!               duidxj (2,1) = dvdx (i,j,k)

!               duidxj (1,2) = dudy (i,j,k)
!               duidxj (2,2) = dvdy (i,j,k)

!            else if ( ndim == 3 ) then

!               duidxj (1,1) = dudx (i,j,k)
!               duidxj (2,1) = dvdx (i,j,k)
!               duidxj (3,1) = dwdx (i,j,k)

!               duidxj (1,2) = dudy (i,j,k)
!               duidxj (2,2) = dvdy (i,j,k)
!               duidxj (3,2) = dwdy (i,j,k)

!               duidxj (1,3) = dudz (i,j,k)
!               duidxj (2,3) = dvdz (i,j,k)
!               duidxj (3,3) = dwdz (i,j,k)

!            end if

!            divu = 0
!            do is = 1 , ndim
!               divu = divu + duidxj (is,is)
!            end do
!            wrk = onethird * ( - cst * mu_SGS (i,j,k) * divu + tau_iso_SGS (i,j,k) )

!            do is = 1 , ndim
!               do js = 1 , ndim

!                  sij  = 0.5_dp * ( duidxj (is,js) + duidxj (js,is) )
!                  tau (i,j,k,is,js) = cst * mu_SGS (i,j,k) * sij + kdelta(is,js) * wrk

!               end do
!            end do

!         end do
!      end do
!    end do


    if      ( ndim == 1 ) then

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

               duidxj (1,1) = dudx (i,j,k)

               divu = 0
               do is = 1 , ndim
                  divu = divu + duidxj (is,is)
               end do
               wrk = onethird * ( - cst * mu_SGS (i,j,k) * divu + tau_iso_SGS (i,j,k) )

               do is = 1 , ndim
                  do js = 1 , ndim

                     sij  = 0.5_dp * ( duidxj (is,js) + duidxj (js,is) )
                     tau_SGS (i,j,k,is,js) = cst * mu_SGS (i,j,k) * sij + kdelta(is,js) * wrk

                  end do
               end do

            end do
         end do
       end do

    else if ( ndim == 2 ) then

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

               duidxj (1,1) = dudx (i,j,k)
               duidxj (2,1) = dvdx (i,j,k)

               duidxj (1,2) = dudy (i,j,k)
               duidxj (2,2) = dvdy (i,j,k)

               divu = 0
               do is = 1 , ndim
                  divu = divu + duidxj (is,is)
               end do
               wrk = onethird * ( - cst * mu_SGS (i,j,k) * divu + tau_iso_SGS (i,j,k) )

               do is = 1 , ndim
                  do js = 1 , ndim

                     sij  = 0.5_dp * ( duidxj (is,js) + duidxj (js,is) )
                     tau_SGS (i,j,k,is,js) = cst * mu_SGS (i,j,k) * sij + kdelta(is,js) * wrk

                  end do
               end do

            end do
         end do
       end do

    else if ( ndim == 3 ) then

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

               duidxj (1,1) = dudx (i,j,k)
               duidxj (2,1) = dvdx (i,j,k)
               duidxj (3,1) = dwdx (i,j,k)

               duidxj (1,2) = dudy (i,j,k)
               duidxj (2,2) = dvdy (i,j,k)
               duidxj (3,2) = dwdy (i,j,k)

               duidxj (1,3) = dudz (i,j,k)
               duidxj (2,3) = dvdz (i,j,k)
               duidxj (3,3) = dwdz (i,j,k)

               divu = 0
               do is = 1 , ndim
                  divu = divu + duidxj (is,is)
               end do
               wrk = onethird * ( - cst * mu_SGS (i,j,k) * divu + tau_iso_SGS (i,j,k) )

               do is = 1 , ndim
                  do js = 1 , ndim

                     sij  = 0.5_dp * ( duidxj (is,js) + duidxj (js,is) )
                     tau_SGS (i,j,k,is,js) = cst * mu_SGS (i,j,k) * sij + kdelta(is,js) * wrk

                  end do
               end do

            end do
         end do
       end do

    end if


    deallocate ( dudx , dudy , dudz , &
                 dvdx , dvdy , dvdz , &
                 dwdx , dwdy , dwdz )
    deallocate ( tau_iso_SGS )


  end subroutine tau_ij_SGS


end module SGS_models
