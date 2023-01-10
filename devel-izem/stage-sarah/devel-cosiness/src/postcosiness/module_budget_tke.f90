!------------------------------------------------------------------------------
! MODULE: budget_tke
!------------------------------------------------------------------------------
!> \brief Budget equations for the turbulent kinetic energy (TKE) transport equation.
!!
!! ### Each term has its own sign!
!!
!! In module_stat_plot.f90 specific signs and other customizations are
!! applied to each term to make the _addition_ of all the terms in the
!! equation equal to zero.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module budget_tke


  use parameters

  use parallel

  use adim

  use deriv


  implicit none


contains


!> \brief Calculate the dissipation term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tke_dissipation ( dx_i , dy_i , dz_i , ux , vy , wz , tau , dissip )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i   !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i   !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i   !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux     !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy     !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz     !< z-component of the velocity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau    !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: dissip !< dissipation term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_ux , dz_ux , &
                                                           dx_vy , dy_vy , dz_vy , &
                                                           dx_wz , dy_wz , dz_wz


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               dx_wz (sx:ex,sy:ey,sz:ez) , &
               dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate tke_dissipation'


    if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,1,1) * dx_ux (i,j,k) + tau (i,j,k,1,2) * dy_ux (i,j,k) + &
                                 tau (i,j,k,2,1) * dx_vy (i,j,k) + tau (i,j,k,2,2) * dy_vy (i,j,k)
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,1,1) * dx_ux (i,j,k) + &
                                 tau (i,j,k,1,2) * dy_ux (i,j,k) + &
                                 tau (i,j,k,1,3) * dz_ux (i,j,k) + &
                                 tau (i,j,k,2,1) * dx_vy (i,j,k) + &
                                 tau (i,j,k,2,2) * dy_vy (i,j,k) + &
                                 tau (i,j,k,2,3) * dz_vy (i,j,k) + &
                                 tau (i,j,k,3,1) * dx_wz (i,j,k) + &
                                 tau (i,j,k,3,2) * dy_wz (i,j,k) + &
                                 tau (i,j,k,3,3) * dz_wz (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_ux , dy_ux , dz_ux , &
                 dx_vy , dy_vy , dz_vy , &
                 dx_wz , dy_wz , dz_wz )


  end subroutine tke_dissipation


!> \brief Calculate the pressure-strain term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tke_pressure_strain ( dx_i , dy_i , dz_i , ux , vy , wz , pres , ps )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ps   !< pressure-strain term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_vy , dz_wz


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate tke_pressure_strain'


    if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux )
       call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                ps (i,j,k) = pres (i,j,k) * ( dx_ux (i,j,k) + dy_vy (i,j,k) )
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)
       call dx ( dx_i , ux , dx_ux )
       call dy ( dy_i , vy , dy_vy )
       call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                ps (i,j,k) = pres (i,j,k) * ( dx_ux (i,j,k) + dy_vy (i,j,k) + dz_wz (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( dx_ux , dy_vy , dz_wz )


  end subroutine tke_pressure_strain


!> \brief Calculate the transport term (first call).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tke_transport1 ( adi , rux , rvy , rwz , fux , fvy , fwz , rho , pres , tau , tr )


    type (adi_type) , intent (in)                                        :: adi  !< non-dimensional derived type
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rux  !< x-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rvy  !< y-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rwz  !< z-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fux  !< x-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fvy  !< y-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fwz  !< z-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau  !< Reynolds stress tensor
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: tr   !< transport term


    integer (ip) :: i , j , k
    real (dp)    :: sqgmr


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fux (i,j,k) * fux (i,j,k) * fux (i,j,k)        + &
                               pres (i,j,k) * ( rux (i,j,k) + rux (i,j,k) )                 - &
                               sqgmr * tau (i,j,k,1,1) * ( fux (i,j,k) + fux (i,j,k) )      + &

                               rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fux (i,j,k)        - &
                               sqgmr * tau (i,j,k,2,1) * ( fvy (i,j,k) + fvy (i,j,k) )


                tr (i,j,k,2) = rho (i,j,k) * fux (i,j,k) * fux (i,j,k) * fvy (i,j,k)        - &
                               sqgmr * tau (i,j,k,1,2) * ( fux (i,j,k) + fux (i,j,k) )      + &

                               rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fvy (i,j,k)        + &
                               pres (i,j,k) * ( rvy (i,j,k) + rvy (i,j,k) )                 - &
                               sqgmr * tau (i,j,k,2,2) * ( fvy (i,j,k) + fvy (i,j,k) )


                tr (i,j,k,3) = 0.0_dp


             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

      do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fux (i,j,k) * fux (i,j,k) * fux (i,j,k)        + &
                               pres (i,j,k) * ( rux (i,j,k) + rux (i,j,k) )                 - &
                               sqgmr * tau (i,j,k,1,1) * ( fux (i,j,k) + fux (i,j,k) )      + &

                               rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fux (i,j,k)        - &
                               sqgmr * tau (i,j,k,2,1) * ( fvy (i,j,k) + fvy (i,j,k) )      + &

                               rho (i,j,k) * fwz (i,j,k) * fwz (i,j,k) * fux (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,1) * ( fwz (i,j,k) + fwz (i,j,k) )


                tr (i,j,k,2) = rho (i,j,k) * fux (i,j,k) * fux (i,j,k) * fvy (i,j,k)        - &
                               sqgmr * tau (i,j,k,1,2) * ( fux (i,j,k) + fux (i,j,k) )      + &

                               rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fvy (i,j,k)        + &
                               pres (i,j,k) * ( rvy (i,j,k) + rvy (i,j,k) )                 - &
                               sqgmr * tau (i,j,k,2,2) * ( fvy (i,j,k) + fvy (i,j,k) )      + &

                               rho (i,j,k) * fwz (i,j,k) * fwz (i,j,k) * fvy (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,2) * ( fwz (i,j,k) + fwz (i,j,k) )


                tr (i,j,k,3) = rho (i,j,k) * fux (i,j,k) * fux (i,j,k) * fwz (i,j,k)        - &
                               sqgmr * tau (i,j,k,1,3) * ( fux (i,j,k) + fux (i,j,k) )      + &

                               rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fwz (i,j,k)        - &
                               sqgmr * tau (i,j,k,2,3) * ( fvy (i,j,k) + fvy (i,j,k) )      + &

                               rho (i,j,k) * fwz (i,j,k) * fwz (i,j,k) * fwz (i,j,k)        + &
                               pres (i,j,k) * ( rwz (i,j,k) + rwz (i,j,k) )                 - &
                               sqgmr * tau (i,j,k,3,3) * ( fwz (i,j,k) + fwz (i,j,k) )


             end do
          end do
       end do

    end if


  end subroutine tke_transport1


!> \brief Calculate the transport term (second call).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tke_transport2 ( dx_i , dy_i , dz_i , tr_k , tr )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: tr_k !< previous transport term
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: tr   !< final transport term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: tr_1  , tr_2  , tr_3 , &
                                                           dx_tr , dy_tr , dz_tr


    allocate ( tr_1 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tr_2 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tr_3 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_tr (sx:ex,sy:ey,sz:ez)                  , &
               dy_tr (sx:ex,sy:ey,sz:ez)                  , &
               dz_tr (sx:ex,sy:ey,sz:ez)                  , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate tke_transport2'


    if ( ndim == 2 ) then ! 2D problem

       tr_1 (sx:ex,sy:ey,sz:ez) = tr_k (sx:ex,sy:ey,sz:ez,1)
       tr_2 (sx:ex,sy:ey,sz:ez) = tr_k (sx:ex,sy:ey,sz:ez,2)

       call comm_one (tr_1) ; call comm_one (tr_2)
       call dx ( dx_i , tr_1 , dx_tr )
       call dy ( dy_i , tr_2 , dy_tr )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                tr (i,j,k) = 0.5_dp * ( dx_tr (i,j,k) + dy_tr (i,j,k) )
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       tr_1 (sx:ex,sy:ey,sz:ez) = tr_k (sx:ex,sy:ey,sz:ez,1)
       tr_2 (sx:ex,sy:ey,sz:ez) = tr_k (sx:ex,sy:ey,sz:ez,2)
       tr_3 (sx:ex,sy:ey,sz:ez) = tr_k (sx:ex,sy:ey,sz:ez,3)

       call comm_one (tr_1) ; call comm_one (tr_2) ; call comm_one (tr_3)
       call dx ( dx_i , tr_1 , dx_tr )
       call dy ( dy_i , tr_2 , dy_tr )
       call dz ( dz_i , tr_3 , dz_tr )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                tr (i,j,k) = 0.5_dp * ( dx_tr (i,j,k) + dy_tr (i,j,k) + dz_tr (i,j,k) )
             end do
          end do
       end do

    end if


  end subroutine tke_transport2


!> \brief Calculate the production term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tke_production ( dx_i , dy_i , dz_i , ux , vy , wz , var_ux , var_vy , var_wz , &
                              var_uxvy , var_uxwz , var_vywz , prod )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i     !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i     !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i     !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux       !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy       !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz       !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_ux   !< ux variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vy   !< vy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_wz   !< wz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxvy !< uxvy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxwz !< uxwz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vywz !< vywz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: prod     !< production term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_ux , dz_ux , &
                                                           dx_vy , dy_vy , dz_vy , &
                                                           dx_wz , dy_wz , dz_wz


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               dx_wz (sx:ex,sy:ey,sz:ez) , &
               dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate tke_production'


    if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_ux (i,j,k) * dx_ux (i,j,k)                       + &
                               var_vy (i,j,k) * dy_vy (i,j,k)                       + &
                               var_uxvy (i,j,k) * ( dy_ux (i,j,k) + dx_vy (i,j,k) )
                prod (i,j,k) = - prod (i,j,k)
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_ux (i,j,k) * dx_ux (i,j,k)                       + &
                               var_vy (i,j,k) * dy_vy (i,j,k)                       + &
                               var_wz (i,j,k) * dz_wz (i,j,k)                       + &
                               var_uxvy (i,j,k) * ( dy_ux (i,j,k) + dx_vy (i,j,k) ) + &
                               var_uxwz (i,j,k) * ( dz_ux (i,j,k) + dx_wz (i,j,k) ) + &
                               var_vywz (i,j,k) * ( dz_vy (i,j,k) + dy_wz (i,j,k) )
                prod (i,j,k) = - prod (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_ux , dy_ux , dz_ux , &
                 dx_vy , dy_vy , dz_vy , &
                 dx_wz , dy_wz , dz_wz )


  end subroutine tke_production


!> \brief Calculate the convection term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tke_convection ( dx_i , dy_i , dz_i , rho , ux , vy , wz , K_ , cv )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: K_   !< square of the Favre variance velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: cv   !< convection term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: cvx   , cvy   , cvz   , &
                                                           dx_cv , dy_cv , dz_cv


    allocate ( cvx   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               cvy   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               cvz   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_cv (sx:ex,sy:ey,sz:ez)                   , &
               dy_cv (sx:ex,sy:ey,sz:ez)                   , &
               dz_cv (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate tke_convection'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                cvx (i,j,k) = rho (i,j,k) * ux (i,j,k) * K_ (i,j,k)
                cvy (i,j,k) = rho (i,j,k) * vy (i,j,k) * K_ (i,j,k)
             end do
          end do
       end do

       call comm_one (cvx) ; call comm_one (cvy)
       call dx ( dx_i , cvx , dx_cv )
       call dy ( dy_i , cvy , dy_cv )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                cv (i,j,k) = dx_cv (i,j,k) + dy_cv (i,j,k)
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                cvx (i,j,k) = rho (i,j,k) * ux (i,j,k) * K_ (i,j,k)
                cvy (i,j,k) = rho (i,j,k) * vy (i,j,k) * K_ (i,j,k)
                cvz (i,j,k) = rho (i,j,k) * wz (i,j,k) * K_ (i,j,k)
             end do
          end do
       end do

       call comm_one (cvx) ; call comm_one (cvy) ; call comm_one (cvz)
       call dx ( dx_i , cvx , dx_cv )
       call dy ( dy_i , cvy , dy_cv )
       call dz ( dz_i , cvz , dz_cv )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                cv (i,j,k) = dx_cv (i,j,k) + dy_cv (i,j,k) + dz_cv (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( cvx   , cvy   , cvz   , &
                 dx_cv , dy_cv , dz_cv )


  end subroutine tke_convection



!> \brief Calculate the mass-flux coupling term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine tke_mass_flux_coupl ( adi , dx_i , dy_i , dz_i , ux , vy , wz , pres , tau , sigma )


    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i  !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux    !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vy    !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: wz    !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: pres  !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau   !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: sigma !< mass-flux coupling term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: sqgmr
    real (dp) , dimension (:,:,:) , allocatable         :: tau_11    , tau_12    , tau_13    , &
                                                           tau_21    , tau_22    , tau_23    , &
                                                           tau_31    , tau_32    , tau_33    , &
                                                           dx_tau_11 , dx_tau_21 , dx_tau_31 , &
                                                           dy_tau_12 , dy_tau_22 , dy_tau_32 , &
                                                           dz_tau_13 , dz_tau_23 , dz_tau_33 , &
                                                           dx_pres   , dy_pres   , dz_pres


    allocate ( tau_11 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_12 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_13 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_21 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_22 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_23 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_31 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_32 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_33 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_tau_11 (sx:ex,sy:ey,sz:ez)                , &
               dx_tau_21 (sx:ex,sy:ey,sz:ez)                , &
               dx_tau_31 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_12 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_22 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_32 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_13 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_23 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_33 (sx:ex,sy:ey,sz:ez)                , &
               dx_pres   (sx:ex,sy:ey,sz:ez)                , &
               dy_pres   (sx:ex,sy:ey,sz:ez)                , &
               dz_pres   (sx:ex,sy:ey,sz:ez)                , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate tke_mass_flux_coupl'


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       tau_11 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,1)
       tau_12 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,2)
       tau_21 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,1)
       tau_22 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,2)

       call comm_one (pres)
       call comm_one (tau_11) ; call comm_one (tau_12)
       call comm_one (tau_21) ; call comm_one (tau_22)
       call dx ( dx_i , pres , dx_pres )
       call dx ( dx_i , tau_11 , dx_tau_11 ) ; call dx ( dx_i , tau_21 , dx_tau_21 )
       call dy ( dy_i , pres , dy_pres )
       call dy ( dy_i , tau_12 , dy_tau_12 ) ; call dy ( dy_i , tau_22 , dy_tau_22 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ux (i,j,k) * ( dx_tau_11 (i,j,k) + dy_tau_12 (i,j,k) - dx_pres (i,j,k) ) + &
                                vy (i,j,k) * ( dx_tau_21 (i,j,k) + dy_tau_22 (i,j,k) - dy_pres (i,j,k) )
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       tau_11 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,1)
       tau_12 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,2)
       tau_13 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,3)
       tau_21 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,1)
       tau_22 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,2)
       tau_23 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,3)
       tau_31 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,1)
       tau_32 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,2)
       tau_33 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,3)

       call comm_one (pres)
       call comm_one (tau_11) ; call comm_one (tau_12) ; call comm_one (tau_13)
       call comm_one (tau_21) ; call comm_one (tau_22) ; call comm_one (tau_23)
       call comm_one (tau_31) ; call comm_one (tau_32) ; call comm_one (tau_33)
       call dx ( dx_i , pres , dx_pres )
       call dx ( dx_i , tau_11 , dx_tau_11 ) ; call dx ( dx_i , tau_21 , dx_tau_21 ) ; call dx ( dx_i , tau_31 , dx_tau_31 )
       call dy ( dy_i , pres , dy_pres )
       call dy ( dy_i , tau_12 , dy_tau_12 ) ; call dy ( dy_i , tau_22 , dy_tau_22 ) ; call dy ( dy_i , tau_32 , dy_tau_32 )
       call dz ( dz_i , pres , dz_pres )
       call dz ( dz_i , tau_13 , dz_tau_13 ) ; call dz ( dz_i , tau_23 , dz_tau_23 ) ; call dz ( dz_i , tau_33 , dz_tau_33 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ux (i,j,k) * ( dx_tau_11 (i,j,k) + dy_tau_12 (i,j,k) + dz_tau_13 (i,j,k) - dx_pres (i,j,k) ) + &
                                vy (i,j,k) * ( dx_tau_21 (i,j,k) + dy_tau_22 (i,j,k) + dz_tau_23 (i,j,k) - dy_pres (i,j,k) ) + &
                                wz (i,j,k) * ( dx_tau_31 (i,j,k) + dy_tau_32 (i,j,k) + dz_tau_33 (i,j,k) - dz_pres (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( tau_11    , tau_12    , tau_13    , &
                 tau_21    , tau_22    , tau_23    , &
                 tau_31    , tau_32    , tau_33    , &
                 dx_tau_11 , dx_tau_21 , dx_tau_31 , &
                 dy_tau_12 , dy_tau_22 , dy_tau_32 , &
                 dz_tau_13 , dz_tau_23 , dz_tau_33 , &
                 dx_pres   , dy_pres   , dz_pres )


  end subroutine tke_mass_flux_coupl


end module budget_tke
