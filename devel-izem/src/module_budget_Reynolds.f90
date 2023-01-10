!------------------------------------------------------------------------------
! MODULE: budget_Reynolds
!------------------------------------------------------------------------------
!> \brief Budget equations for the Reynolds stress transport equations.
!!
!! ### Each term has its own sign!
!!
!! In module_stat_plot.f90 specific signs and other customizations are
!! applied to each term to make the _addition_ of all the terms in the
!! equation equal to zero.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module budget_Reynolds


  use parameters

  use parallel

  use adim

  use deriv


  implicit none


contains


!> \brief Calculate the 11-component of the dissipation term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_dissipation_11 ( dx_i , dy_i , dz_i , ux , tau , dissip )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i   !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i   !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i   !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux     !< x-component of the velocity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau    !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: dissip !< component 11 of the dissipation term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_ux , dz_ux


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_dissipation_11'


    if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,1,1) * dx_ux (i,j,k) + tau (i,j,k,1,2) * dy_ux (i,j,k)
                dissip (i,j,k) = dissip (i,j,k) + dissip (i,j,k)
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,1,1) * dx_ux (i,j,k) + &
                                 tau (i,j,k,1,2) * dy_ux (i,j,k) + &
                                 tau (i,j,k,1,3) * dz_ux (i,j,k)
                dissip (i,j,k) = dissip (i,j,k) + dissip (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_ux , dy_ux , dz_ux )


  end subroutine rey_dissipation_11


!> \brief Calculate the 22-component of the dissipation term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_dissipation_22 ( dx_i , dy_i , dz_i , vy , tau , dissip )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i   !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i   !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i   !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy     !< y-component of the velocity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau    !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: dissip !< component 22 of the dissipation term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_vy , dy_vy , dz_vy


    allocate ( dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_dissipation_22'


    if ( ndim == 2 ) then ! 2D problem

       call comm_one (vy)
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,2,1) * dx_vy (i,j,k) + tau (i,j,k,2,2) * dy_vy (i,j,k)
                dissip (i,j,k) = dissip (i,j,k) + dissip (i,j,k)
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (vy)
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,2,1) * dx_vy (i,j,k) + &
                                 tau (i,j,k,2,2) * dy_vy (i,j,k) + &
                                 tau (i,j,k,2,3) * dz_vy (i,j,k)
                dissip (i,j,k) = dissip (i,j,k) + dissip (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_vy , dy_vy , dz_vy )


  end subroutine rey_dissipation_22


!> \brief Calculate the 33-component of the dissipation term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_dissipation_33 ( dx_i , dy_i , dz_i , wz , tau , dissip )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i    !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i    !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i    !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz      !< z-component of the velocity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau     !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: dissip  !< component 33 of the dissipation term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_wz , dy_wz , dz_wz


    allocate ( dx_wz (sx:ex,sy:ey,sz:ez) , &
               dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_dissipation_33'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (wz)
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,3,1) * dx_wz (i,j,k) + &
                                 tau (i,j,k,3,2) * dy_wz (i,j,k) + &
                                 tau (i,j,k,3,3) * dz_wz (i,j,k)
                dissip (i,j,k) = dissip (i,j,k) + dissip (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_wz , dy_wz , dz_wz )


  end subroutine rey_dissipation_33


!> \brief Calculate the 12-component of the dissipation term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_dissipation_12 ( dx_i , dy_i , dz_i , ux , vy , tau , dissip )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i   !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i   !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i   !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux     !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy     !< y-component of the velocity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau    !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: dissip !< component 12 of the dissipation term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_ux , dz_ux , &
                                                           dx_vy , dy_vy , dz_vy


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_dissipation_12'


    if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,2,1) * dx_ux (i,j,k) + tau (i,j,k,2,2) * dy_ux (i,j,k) + &
                                 tau (i,j,k,1,1) * dx_vy (i,j,k) + tau (i,j,k,1,2) * dy_vy (i,j,k)
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,2,1) * dx_ux (i,j,k) + &
                                 tau (i,j,k,2,2) * dy_ux (i,j,k) + &
                                 tau (i,j,k,2,3) * dz_ux (i,j,k) + &
                                 tau (i,j,k,1,1) * dx_vy (i,j,k) + &
                                 tau (i,j,k,1,2) * dy_vy (i,j,k) + &
                                 tau (i,j,k,1,3) * dz_vy (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_ux , dy_ux , dz_ux , &
                 dx_vy , dy_vy , dz_vy )


  end subroutine rey_dissipation_12


!> \brief Calculate the 13-component of the dissipation term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_dissipation_13 ( dx_i , dy_i , dz_i , ux , wz , tau , dissip )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i   !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i   !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i   !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux     !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz     !< z-component of the velocity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau    !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: dissip !< component 13 of the dissipation term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_ux , dz_ux , &
                                                           dx_wz , dy_wz , dz_wz


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               dx_wz (sx:ex,sy:ey,sz:ez) , &
               dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_dissipation_13'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (wz)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,3,1) * dx_ux (i,j,k) + &
                                 tau (i,j,k,3,2) * dy_ux (i,j,k) + &
                                 tau (i,j,k,3,3) * dz_ux (i,j,k) + &
                                 tau (i,j,k,1,1) * dx_wz (i,j,k) + &
                                 tau (i,j,k,1,2) * dy_wz (i,j,k) + &
                                 tau (i,j,k,1,3) * dz_wz (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_ux , dy_ux , dz_ux , &
                 dx_wz , dy_wz , dz_wz )


  end subroutine rey_dissipation_13


!> \brief Calculate the 23-component of the dissipation term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_dissipation_23 ( dx_i , dy_i , dz_i , vy , wz , tau , dissip )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i   !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i   !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i   !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy     !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz     !< z-component of the velocity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau    !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: dissip !< component 23 of the dissipation term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_vy , dy_vy , dz_vy , &
                                                           dx_wz , dy_wz , dz_wz


    allocate ( dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               dx_wz (sx:ex,sy:ey,sz:ez) , &
               dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_dissipation_13'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (vy) ; call comm_one (wz)
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dissip (i,j,k) = tau (i,j,k,3,1) * dx_vy (i,j,k) + &
                                 tau (i,j,k,3,2) * dy_vy (i,j,k) + &
                                 tau (i,j,k,3,3) * dz_vy (i,j,k) + &
                                 tau (i,j,k,1,1) * dx_wz (i,j,k) + &
                                 tau (i,j,k,1,2) * dy_wz (i,j,k) + &
                                 tau (i,j,k,1,3) * dz_wz (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_vy , dy_vy , dz_vy , &
                 dx_wz , dy_wz , dz_wz )


  end subroutine rey_dissipation_23


!> \brief Calculate the 11-component of the pressure-strain term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_pressure_strain_11 ( dx_i , ux , pres , ps )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ps   !< component 11 of the pressure-strain


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_pressure_strain_11'


    call comm_one (ux)
    call dx ( dx_i , ux , dx_ux )

    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             ps (i,j,k) = pres (i,j,k) * dx_ux (i,j,k)
             ps (i,j,k) = ps (i,j,k) + ps (i,j,k)
          end do
       end do
    end do

    deallocate (dx_ux)


  end subroutine rey_pressure_strain_11


!> \brief Calculate the 22-component of the pressure-strain term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_pressure_strain_22 ( dy_i , vy , pres , ps )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ps   !< component 22 of the pressure-strain


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dy_vy


    allocate ( dy_vy (sx:ex,sy:ey,sz:ez) , stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_pressure_strain_22'


    call comm_one (vy)
    call dy ( dy_i , vy , dy_vy )

    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             ps (i,j,k) = pres (i,j,k) * dy_vy (i,j,k)
             ps (i,j,k) = ps (i,j,k) + ps (i,j,k)
          end do
       end do
    end do

    deallocate (dy_vy)


  end subroutine rey_pressure_strain_22


!> \brief Calculate the 33-component of the pressure-strain term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_pressure_strain_33 ( dz_i , wz , pres , ps )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ps   !< component 33 of the pressure-strain


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dz_wz


    allocate ( dz_wz (sx:ex,sy:ey,sz:ez) , stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_pressure_strain_33'


    call comm_one (wz)
    call dz ( dz_i , wz , dz_wz )

    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             ps (i,j,k) = pres (i,j,k) * dz_wz (i,j,k)
             ps (i,j,k) = ps (i,j,k) + ps (i,j,k)
          end do
       end do
    end do

    deallocate (dz_wz)


  end subroutine rey_pressure_strain_33


!> \brief Calculate the 12-component of the pressure-strain term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_pressure_strain_12 ( dx_i , dy_i , ux , vy , pres , ps )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ps   !< component 12 of the pressure-strain


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_vy , dy_ux


    allocate ( dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_pressure_strain_12'


    call comm_one (ux)            ; call comm_one (vy)
    call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , ux , dy_ux )

    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             ps (i,j,k) = pres (i,j,k) * ( dy_ux (i,j,k) + dx_vy (i,j,k) )
          end do
       end do
    end do


    deallocate ( dx_vy , dy_ux )


  end subroutine rey_pressure_strain_12


!> \brief Calculate the 13-component of the pressure-strain term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_pressure_strain_13 ( dx_i , dz_i , ux , wz , pres , ps )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ps   !< component 13 of the pressure-strain


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_wz , dz_ux


    allocate ( dx_wz (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_pressure_strain_13'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                ps (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux)            ; call comm_one (wz)
       call dx ( dx_i , wz , dx_wz ) ; call dz ( dz_i , ux , dz_ux )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                ps (i,j,k) = pres (i,j,k) * ( dz_ux (i,j,k) + dx_wz (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( dz_ux , dx_wz )


  end subroutine rey_pressure_strain_13


!> \brief Calculate the 23-component of the pressure-strain term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_pressure_strain_23 ( dy_i , dz_i , vy , wz , pres , ps )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ps   !< component 23 of the pressure-strain


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dy_wz , dz_vy


    allocate ( dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_pressure_strain_23'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                ps (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (vy)            ; call comm_one (wz)
       call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , vy , dz_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                ps (i,j,k) = pres (i,j,k) * ( dz_vy (i,j,k) + dy_wz (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( dz_vy , dy_wz )


  end subroutine rey_pressure_strain_23


!> \brief Calculate the 11-component of the transport term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_transport1_11 ( adi , rux , fux , fvy , fwz , rho , pres , tau , tr )


    type (adi_type) , intent (in)                                        :: adi  !< non-dimensional derived type
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rux  !< x-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fux  !< x-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fvy  !< y-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fwz  !< z-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau  !< Reynolds stress tensor
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: tr   !< component 11 of the transport term


    integer (ip) :: i , j , k
    real (dp)    :: sqgmr


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fux (i,j,k) * fux (i,j,k) * fux (i,j,k)        + &
                               pres (i,j,k) * ( rux (i,j,k) + rux (i,j,k) )                 - &
                               sqgmr * tau (i,j,k,1,1) * ( fux (i,j,k) + fux (i,j,k) )

                tr (i,j,k,2) = rho (i,j,k) * fux (i,j,k) * fux (i,j,k) * fvy (i,j,k)        - &
                               sqgmr * tau (i,j,k,1,2) * ( fux (i,j,k) + fux (i,j,k) )

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
                               sqgmr * tau (i,j,k,1,1) * ( fux (i,j,k) + fux (i,j,k) )

                tr (i,j,k,2) = rho (i,j,k) * fux (i,j,k) * fux (i,j,k) * fvy (i,j,k)        - &
                               sqgmr * tau (i,j,k,1,2) * ( fux (i,j,k) + fux (i,j,k) )

                tr (i,j,k,3) = rho (i,j,k) * fux (i,j,k) * fux (i,j,k) * fwz (i,j,k)        - &
                               sqgmr * tau (i,j,k,1,3) * ( fux (i,j,k) + fux (i,j,k) )


             end do
          end do
       end do

    end if


  end subroutine rey_transport1_11


!> \brief Calculate the 22-component of the transport term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_transport1_22 ( adi , rvy , fux , fvy , fwz , rho , pres , tau , tr )


    type (adi_type) , intent (in)                                        :: adi  !< non-dimensional derived type
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rvy  !< y-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fux  !< x-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fvy  !< y-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fwz  !< z-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau  !< Reynolds stress tensor
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: tr   !< component 22 of the transport term


    integer (ip) :: i , j , k
    real (dp)    :: sqgmr


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fux (i,j,k)        - &
                               sqgmr * tau (i,j,k,2,1) * ( fvy (i,j,k) + fvy (i,j,k) )

                tr (i,j,k,2) = rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fvy (i,j,k)        + &
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


                tr (i,j,k,1) = rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fux (i,j,k)        - &
                               sqgmr * tau (i,j,k,2,1) * ( fvy (i,j,k) + fvy (i,j,k) )

                tr (i,j,k,2) = rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fvy (i,j,k)        + &
                               pres (i,j,k) * ( rvy (i,j,k) + rvy (i,j,k) )                 - &
                               sqgmr * tau (i,j,k,2,2) * ( fvy (i,j,k) + fvy (i,j,k) )

                tr (i,j,k,3) = rho (i,j,k) * fvy (i,j,k) * fvy (i,j,k) * fwz (i,j,k)        - &
                               sqgmr * tau (i,j,k,2,3) * ( fvy (i,j,k) + fvy (i,j,k) )


             end do
          end do
       end do

    end if


  end subroutine rey_transport1_22


!> \brief Calculate the 33-component of the transport term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_transport1_33 ( adi , rwz , fux , fvy , fwz , rho , pres , tau , tr )


    type (adi_type) , intent (in)                                        :: adi  !< non-dimensional derived type
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rwz  !< x-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fux  !< x-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fvy  !< y-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fwz  !< z-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau  !< Reynolds stress tensor
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: tr   !< component 11 of the transport term


    integer (ip) :: i , j , k
    real (dp)    :: sqgmr


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fwz (i,j,k) * fwz (i,j,k) * fux (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,1) * ( fwz (i,j,k) + fwz (i,j,k) )

                tr (i,j,k,2) = rho (i,j,k) * fwz (i,j,k) * fwz (i,j,k) * fvy (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,2) * ( fwz (i,j,k) + fwz (i,j,k) )

                tr (i,j,k,3) = 0.0_dp


             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

      do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fwz (i,j,k) * fwz (i,j,k) * fux (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,1) * ( fwz (i,j,k) + fwz (i,j,k) )

                tr (i,j,k,2) = rho (i,j,k) * fwz (i,j,k) * fwz (i,j,k) * fvy (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,2) * ( fwz (i,j,k) + fwz (i,j,k) )

                tr (i,j,k,3) = rho (i,j,k) * fwz (i,j,k) * fwz (i,j,k) * fwz (i,j,k)        + &
                               pres (i,j,k) * ( rwz (i,j,k) + rwz (i,j,k) )                 - &
                               sqgmr * tau (i,j,k,3,3) * ( fwz (i,j,k) + fwz (i,j,k) )


             end do
          end do
       end do

    end if


  end subroutine rey_transport1_33


!> \brief Calculate the 12-component of the transport term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_transport1_12 ( adi , rux , rvy , fux , fvy , fwz , rho , pres , tau , tr )


    type (adi_type) , intent (in)                                        :: adi  !< non-dimensional derived type
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rux  !< x-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rvy  !< y-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fux  !< x-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fvy  !< y-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fwz  !< z-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau  !< Reynolds stress tensor
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: tr   !< component 12 of the transport term


    integer (ip) :: i , j , k
    real (dp)    :: sqgmr


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fux (i,j,k) * fvy (i,j,k) * fux (i,j,k)        + &
                               pres (i,j,k) * rvy (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,2,1) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,1) * fvy (i,j,k)

                tr (i,j,k,2) = rho (i,j,k) * fux (i,j,k) * fvy (i,j,k) * fvy (i,j,k)        + &
                               pres (i,j,k) * rux (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,2,2) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,2) * fvy (i,j,k)

                tr (i,j,k,3) = 0.0_dp


             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

      do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fux (i,j,k) * fvy (i,j,k) * fux (i,j,k)        + &
                               pres (i,j,k) * rvy (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,2,1) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,1) * fvy (i,j,k)

                tr (i,j,k,2) = rho (i,j,k) * fux (i,j,k) * fvy (i,j,k) * fvy (i,j,k)        + &
                               pres (i,j,k) * rux (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,2,2) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,2) * fvy (i,j,k)

                tr (i,j,k,3) = rho (i,j,k) * fux (i,j,k) * fvy (i,j,k) * fwz (i,j,k)        - &
                               sqgmr * tau (i,j,k,2,3) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,3) * fvy (i,j,k)


             end do
          end do
       end do

    end if


  end subroutine rey_transport1_12


!> \brief Calculate the 13-component of the transport term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_transport1_13 ( adi , rux , rwz , fux , fvy , fwz , rho , pres , tau , tr )


    type (adi_type) , intent (in)                                        :: adi  !< non-dimensional derived type
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rux  !< x-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rwz  !< z-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fux  !< x-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fvy  !< y-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fwz  !< z-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau  !< Reynolds stress tensor
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: tr   !< component 12 of the transport term


    integer (ip) :: i , j , k
    real (dp)    :: sqgmr


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fux (i,j,k) * fwz (i,j,k) * fux (i,j,k)        + &
                               pres (i,j,k) * rwz (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,3,1) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,1) * fwz (i,j,k)

                tr (i,j,k,2) = rho (i,j,k) * fux (i,j,k) * fwz (i,j,k) * fvy (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,2) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,2) * fvy (i,j,k)

                tr (i,j,k,3) = 0.0_dp


             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

      do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fux (i,j,k) * fwz (i,j,k) * fux (i,j,k)        + &
                               pres (i,j,k) * rwz (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,3,1) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,1) * fwz (i,j,k)

                tr (i,j,k,2) = rho (i,j,k) * fux (i,j,k) * fwz (i,j,k) * fvy (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,2) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,2) * fvy (i,j,k)

                tr (i,j,k,3) = rho (i,j,k) * fux (i,j,k) * fwz (i,j,k) * fwz (i,j,k)        + &
                               pres (i,j,k) * rux (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,3,3) * fux (i,j,k)                        - &
                               sqgmr * tau (i,j,k,1,3) * fwz (i,j,k)


             end do
          end do
       end do

    end if


  end subroutine rey_transport1_13


!> \brief Calculate the 23-component of the transport term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_transport1_23 ( adi , rvy , rwz , fux , fvy , fwz , rho , pres , tau , tr )


    type (adi_type) , intent (in)                                        :: adi  !< non-dimensional derived type
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rvy  !< y-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rwz  !< z-component of the Reynolds fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fux  !< x-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fvy  !< y-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fwz  !< z-component of the Favre fluctuating velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: pres !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau  !< Reynolds stress tensor
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: tr   !< component 23 of the transport term


    integer (ip) :: i , j , k
    real (dp)    :: sqgmr


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fvy (i,j,k) * fwz (i,j,k) * fux (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,1) * fvy (i,j,k)                        - &
                               sqgmr * tau (i,j,k,2,1) * fwz (i,j,k)

                tr (i,j,k,2) = rho (i,j,k) * fvy (i,j,k) * fwz (i,j,k) * fvy (i,j,k)        + &
                               pres (i,j,k) * rwz (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,3,2) * fvy (i,j,k)                        - &
                               sqgmr * tau (i,j,k,2,2) * fwz (i,j,k)

                tr (i,j,k,3) = 0.0_dp


             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

      do k = sz,ez
          do j = sy,ey
             do i = sx,ex


                tr (i,j,k,1) = rho (i,j,k) * fvy (i,j,k) * fwz (i,j,k) * fux (i,j,k)        - &
                               sqgmr * tau (i,j,k,3,1) * fvy (i,j,k)                        - &
                               sqgmr * tau (i,j,k,2,1) * fwz (i,j,k)

                tr (i,j,k,2) = rho (i,j,k) * fvy (i,j,k) * fwz (i,j,k) * fvy (i,j,k)        + &
                               pres (i,j,k) * rwz (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,3,2) * fvy (i,j,k)                        - &
                               sqgmr * tau (i,j,k,2,2) * fwz (i,j,k)

                tr (i,j,k,3) = rho (i,j,k) * fvy (i,j,k) * fwz (i,j,k) * fwz (i,j,k)        + &
                               pres (i,j,k) * rvy (i,j,k)                                   - &
                               sqgmr * tau (i,j,k,3,3) * fvy (i,j,k)                        - &
                               sqgmr * tau (i,j,k,2,3) * fwz (i,j,k)


             end do
          end do
       end do

    end if


  end subroutine rey_transport1_23


!> \brief Calculate the transport term (second call).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_transport2 ( dx_i , dy_i , dz_i , tr_k , tr )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: tr_k !< TKE transport term
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: tr   !< transport term


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
                tr (i,j,k) = dx_tr (i,j,k) + dy_tr (i,j,k)
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
                tr (i,j,k) = dx_tr (i,j,k) + dy_tr (i,j,k) + dz_tr (i,j,k)
             end do
          end do
       end do

    end if


  end subroutine rey_transport2


!> \brief Calculate the 11-component of the production term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_production_11 ( dx_i , dy_i , dz_i , ux , var_uxux , &
                                 var_uxvy , var_uxwz , prod )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i     !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i     !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i     !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux       !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxux !< uxux variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxvy !< uxvy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxwz !< uxwz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: prod     !< component 11 of the production term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_ux , dz_ux


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_production_11'


    if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_uxux (i,j,k) * dx_ux (i,j,k) + &
                               var_uxvy (i,j,k) * dy_ux (i,j,k)
                prod (i,j,k) = - ( prod (i,j,k) + prod (i,j,k) )
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_uxux (i,j,k) * dx_ux (i,j,k) + &
                               var_uxvy (i,j,k) * dy_ux (i,j,k) + &
                               var_uxwz (i,j,k) * dz_ux (i,j,k)
                prod (i,j,k) = - ( prod (i,j,k) + prod (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( dx_ux , dy_ux , dz_ux )


  end subroutine rey_production_11


!> \brief Calculate the 11-component of the production term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_production_22 ( dx_i , dy_i , dz_i , vy , var_vyvy , &
                                 var_uxvy , var_vywz , prod )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i     !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i     !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i     !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy       !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vyvy !< vyvy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxvy !< uxvy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vywz !< vywz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: prod     !< component 22 of the production term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_vy , dy_vy , dz_vy


    allocate ( dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_production_22'


    if ( ndim == 2 ) then ! 2D problem

       call comm_one (vy)
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_uxvy (i,j,k) * dx_vy (i,j,k) + &
                               var_vyvy (i,j,k) * dy_vy (i,j,k)
                prod (i,j,k) = - ( prod (i,j,k) + prod (i,j,k) )
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (vy)
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_uxvy (i,j,k) * dx_vy (i,j,k) + &
                               var_vyvy (i,j,k) * dy_vy (i,j,k) + &
                               var_vywz (i,j,k) * dz_vy (i,j,k)
                prod (i,j,k) = - ( prod (i,j,k) + prod (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( dx_vy , dy_vy , dz_vy )


  end subroutine rey_production_22


!> \brief Calculate the 33-component of the production term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_production_33 ( dx_i , dy_i , dz_i , wz , var_wzwz , &
                                 var_uxwz , var_vywz , prod )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i     !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i     !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i     !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz       !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_wzwz !< wzwz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxwz !< uxwz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vywz !< vywz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: prod     !< component 33 of the production term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_wz , dy_wz , dz_wz


    allocate ( dx_wz (sx:ex,sy:ey,sz:ez) , &
               dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_production_33'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (wz)
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_uxwz (i,j,k) * dx_wz (i,j,k) + &
                               var_vywz (i,j,k) * dy_wz (i,j,k) + &
                               var_wzwz (i,j,k) * dz_wz (i,j,k)
                prod (i,j,k) = - ( prod (i,j,k) + prod (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( dx_wz , dy_wz , dz_wz )


  end subroutine rey_production_33


!> \brief Calculate the 12-component of the production term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_production_12 ( dx_i , dy_i , dz_i , ux , vy , var_uxux , var_vyvy , &
                                 var_uxvy , var_uxwz , var_vywz , prod )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i     !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i     !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i     !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux       !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy       !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxux !< uxux variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vyvy !< vyvy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxvy !< uxvy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxwz !< uxwz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vywz !< vywz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: prod     !< component 12 of the production term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_ux , dz_ux , &
                                                           dx_vy , dy_vy , dz_vy


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_production_12'


    if ( ndim == 2 ) then ! 2D problem

       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_uxux (i,j,k) * dx_vy (i,j,k) + var_uxvy (i,j,k) * dx_ux (i,j,k) + &
                               var_uxvy (i,j,k) * dy_vy (i,j,k) + var_vyvy (i,j,k) * dy_ux (i,j,k)
                prod (i,j,k) = - prod (i,j,k)
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_uxux (i,j,k) * dx_vy (i,j,k) + var_uxvy (i,j,k) * dx_ux (i,j,k) + &
                               var_uxvy (i,j,k) * dy_vy (i,j,k) + var_vyvy (i,j,k) * dy_ux (i,j,k) + &
                               var_uxwz (i,j,k) * dz_vy (i,j,k) + var_vywz (i,j,k) * dz_ux (i,j,k)
                prod (i,j,k) = - prod (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_ux , dy_ux , dz_ux , &
                 dx_vy , dy_vy , dz_vy )


  end subroutine rey_production_12


!> \brief Calculate the 13-component of the production term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_production_13 ( dx_i , dy_i , dz_i , ux , wz , var_uxux , var_wzwz , &
                                 var_uxvy , var_uxwz , var_vywz , prod )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i     !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i     !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i     !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux       !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz       !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxux !< uxux variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_wzwz !< wzwz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxvy !< uxvy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxwz !< uxwz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vywz !< vywz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: prod     !< component 13 of the production term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_ux , dz_ux , &
                                                           dx_wz , dy_wz , dz_wz


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               dx_wz (sx:ex,sy:ey,sz:ez) , &
               dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_production_13'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (ux) ; call comm_one (wz)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_uxux (i,j,k) * dx_wz (i,j,k) + var_uxwz (i,j,k) * dx_ux (i,j,k) + &
                               var_uxvy (i,j,k) * dy_wz (i,j,k) + var_vywz (i,j,k) * dy_ux (i,j,k) + &
                               var_uxwz (i,j,k) * dz_wz (i,j,k) + var_wzwz (i,j,k) * dz_ux (i,j,k)
                prod (i,j,k) = - prod (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_ux , dy_ux , dz_ux , &
                 dx_wz , dy_wz , dz_wz )


  end subroutine rey_production_13


!> \brief Calculate the 23-component of the production term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_production_23 ( dx_i , dy_i , dz_i , vy , wz , var_vyvy , var_wzwz , &
                                 var_uxvy , var_uxwz , var_vywz , prod )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i     !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i     !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i     !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy       !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz       !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vyvy !< vyvy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_wzwz !< wzwz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxvy !< uxvy variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_uxwz !< uxwz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var_vywz !< vywz variance
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: prod     !< component 23 of the production term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dx_vy , dy_vy , dz_vy , &
                                                           dx_wz , dy_wz , dz_wz


    allocate ( dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               dx_wz (sx:ex,sy:ey,sz:ez) , &
               dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_production_23'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       call comm_one (vy) ; call comm_one (wz)
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                prod (i,j,k) = var_uxvy (i,j,k) * dx_wz (i,j,k) + var_uxwz (i,j,k) * dx_vy (i,j,k) + &
                               var_vyvy (i,j,k) * dy_wz (i,j,k) + var_vywz (i,j,k) * dy_vy (i,j,k) + &
                               var_vywz (i,j,k) * dz_wz (i,j,k) + var_wzwz (i,j,k) * dz_vy (i,j,k)
                prod (i,j,k) = - prod (i,j,k)
             end do
          end do
       end do

    end if


    deallocate ( dx_vy , dy_vy , dz_vy , &
                 dx_wz , dy_wz , dz_wz )


  end subroutine rey_production_23


!> \brief Calculate the convection term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_convection ( dx_i , dy_i , dz_i , rho , ux , vy , wz , var , cv )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dx array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var  !< velocity component variance
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
    if ( ok > 0 ) stop 'error allocate rey_convection'


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                cvx (i,j,k) = rho (i,j,k) * ux (i,j,k) * var (i,j,k)
                cvy (i,j,k) = rho (i,j,k) * vy (i,j,k) * var (i,j,k)
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
                cvx (i,j,k) = rho (i,j,k) * ux (i,j,k) * var (i,j,k)
                cvy (i,j,k) = rho (i,j,k) * vy (i,j,k) * var (i,j,k)
                cvz (i,j,k) = rho (i,j,k) * wz (i,j,k) * var (i,j,k)
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


  end subroutine rey_convection


!> \brief Calculate the 11-component of the mass-flux coupling term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_mass_flux_coupl_11 ( adi , dx_i , dy_i , dz_i , ux , pres , tau , sigma )


    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i  !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux    !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: pres  !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau   !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: sigma !< component 11 of the mass-flux coupling term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: sqgmr
    real (dp) , dimension (:,:,:) , allocatable         :: tau_11    , tau_12    , tau_13    , &
                                                           dx_tau_11 , dy_tau_12 , dz_tau_13 , &
                                                           dx_pres


    allocate ( tau_11 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_12 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_13 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_tau_11 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_12 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_13 (sx:ex,sy:ey,sz:ez)                , &
               dx_pres   (sx:ex,sy:ey,sz:ez)                , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_mass_flux_coupl_11'


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       tau_11 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,1)
       tau_12 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,2)

       call comm_one (pres)
       call comm_one (tau_11) ; call comm_one (tau_12)
       call dx ( dx_i , pres , dx_pres )
       call dx ( dx_i , tau_11 , dx_tau_11 ) ; call dy ( dy_i , tau_12 , dy_tau_12 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ( ux (i,j,k) + ux (i,j,k) ) *                               &
                                ( dx_tau_11 (i,j,k) + dy_tau_12 (i,j,k) - dx_pres (i,j,k) )
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       tau_11 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,1)
       tau_12 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,2)
       tau_13 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,3)

       call comm_one (pres)
       call comm_one (tau_11) ; call comm_one (tau_12) ; call comm_one (tau_13)
       call dx ( dx_i , pres , dx_pres )
       call dx ( dx_i , tau_11 , dx_tau_11 ) ; call dy ( dy_i , tau_12 , dy_tau_12 ) ; call dz ( dz_i , tau_13 , dz_tau_13 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ( ux (i,j,k) + ux (i,j,k) ) *                                                   &
                                ( dx_tau_11 (i,j,k) + dy_tau_12 (i,j,k) + dz_tau_13 (i,j,k) - dx_pres (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( tau_11    , tau_12    , tau_13    , &
                 dx_tau_11 , dy_tau_12 , dz_tau_13 , &
                 dx_pres )


  end subroutine rey_mass_flux_coupl_11


!> \brief Calculate the 22-component of the mass-flux coupling term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_mass_flux_coupl_22 ( adi , dx_i , dy_i , dz_i , vy , pres , tau , sigma )


    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i  !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vy    !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: pres  !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau   !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: sigma !< component 22 of the mass-flux coupling term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: sqgmr
    real (dp) , dimension (:,:,:) , allocatable         :: tau_21    , tau_22    , tau_23    , &
                                                           dx_tau_21 , dy_tau_22 , dz_tau_23 , &
                                                           dy_pres


    allocate ( tau_21 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_22 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_23 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_tau_21 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_22 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_23 (sx:ex,sy:ey,sz:ez)                , &
               dy_pres   (sx:ex,sy:ey,sz:ez)                , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_mass_flux_coupl_22'


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       tau_21 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,1)
       tau_22 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,2)

       call comm_one (pres)
       call comm_one (tau_21) ; call comm_one (tau_22)
       call dy ( dy_i , pres , dy_pres )
       call dx ( dx_i , tau_21 , dx_tau_21 ) ; call dy ( dy_i , tau_22 , dy_tau_22 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ( vy (i,j,k) + vy (i,j,k) ) *                               &
                                ( dx_tau_21 (i,j,k) + dy_tau_22 (i,j,k) - dy_pres (i,j,k) )
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       tau_21 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,1)
       tau_22 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,2)
       tau_23 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,3)

       call comm_one (pres)
       call comm_one (tau_21) ; call comm_one (tau_22) ; call comm_one (tau_23)
       call dy ( dy_i , pres , dy_pres )
       call dx ( dx_i , tau_21 , dx_tau_21 ) ; call dy ( dy_i , tau_22 , dy_tau_22 ) ; call dz ( dz_i , tau_23 , dz_tau_23 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ( vy (i,j,k) + vy (i,j,k) ) *                                                   &
                                ( dx_tau_21 (i,j,k) + dy_tau_22 (i,j,k) + dz_tau_23 (i,j,k) - dy_pres (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( tau_21    , tau_22    , tau_23    , &
                 dx_tau_21 , dy_tau_22 , dz_tau_23 , &
                 dy_pres )


  end subroutine rey_mass_flux_coupl_22


!> \brief Calculate the 33-component of the mass-flux coupling term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_mass_flux_coupl_33 ( adi , dx_i , dy_i , dz_i , wz , pres , tau , sigma )


    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i  !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: wz    !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: pres  !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau   !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: sigma !< component 33 of the mass-flux coupling term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: sqgmr
    real (dp) , dimension (:,:,:) , allocatable         :: tau_31    , tau_32    , tau_33    , &
                                                           dx_tau_31 , dy_tau_32 , dz_tau_33 , &
                                                           dz_pres


    allocate ( tau_31 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_32 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_33 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_tau_31 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_32 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_33 (sx:ex,sy:ey,sz:ez)                , &
               dz_pres   (sx:ex,sy:ey,sz:ez)                , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_mass_flux_coupl_33'


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       tau_31 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,1)
       tau_32 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,2)
       tau_33 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,3)

       call comm_one (pres)
       call comm_one (tau_31) ; call comm_one (tau_32) ; call comm_one (tau_33)
       call dz ( dz_i , pres , dz_pres )
       call dx ( dx_i , tau_31 , dx_tau_31 ) ; call dy ( dy_i , tau_32 , dy_tau_32 ) ; call dz ( dz_i , tau_33 , dz_tau_33 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ( wz (i,j,k) + wz (i,j,k) ) *                                                   &
                                ( dx_tau_31 (i,j,k) + dy_tau_32 (i,j,k) + dz_tau_33 (i,j,k) - dz_pres (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( tau_31    , tau_32    , tau_33    , &
                 dx_tau_31 , dy_tau_32 , dz_tau_33 , &
                 dz_pres )


  end subroutine rey_mass_flux_coupl_33


!> \brief Calculate the 12-component of the mass-flux coupling term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_mass_flux_coupl_12 ( adi , dx_i , dy_i , dz_i , ux , vy , pres , tau , sigma )


    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i  !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux    !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vy    !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: pres  !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau   !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: sigma !< component 12 of the mass-flux coupling term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: sqgmr
    real (dp) , dimension (:,:,:) , allocatable         :: tau_11    , tau_12    , tau_13    , &
                                                           tau_21    , tau_22    , tau_23    , &
                                                           dx_tau_11 , dy_tau_12 , dz_tau_13 , &
                                                           dx_tau_21 , dy_tau_22 , dz_tau_23 , &
                                                           dx_pres   , dy_pres


    allocate ( tau_11 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_12 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_13 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_21 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_22 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_23 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_tau_11 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_12 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_13 (sx:ex,sy:ey,sz:ez)                , &
               dx_tau_21 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_22 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_23 (sx:ex,sy:ey,sz:ez)                , &
               dx_pres   (sx:ex,sy:ey,sz:ez)                , &
               dy_pres   (sx:ex,sy:ey,sz:ez)                , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_mass_flux_coupl_12'


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       tau_11 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,1)
       tau_12 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,2)
       tau_21 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,1)
       tau_22 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,2)

       call comm_one (pres)
       call comm_one (tau_11) ; call comm_one (tau_12)
       call comm_one (tau_21) ; call comm_one (tau_22)
       call dx ( dx_i , pres , dx_pres )     ; call dy ( dy_i , pres , dy_pres )
       call dx ( dx_i , tau_11 , dx_tau_11 ) ; call dy ( dy_i , tau_12 , dy_tau_12 )
       call dx ( dx_i , tau_21 , dx_tau_21 ) ; call dy ( dy_i , tau_22 , dy_tau_22 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ux (i,j,k) * ( dx_tau_21 (i,j,k) + dy_tau_22 (i,j,k) - dy_pres (i,j,k) ) + &
                                vy (i,j,k) * ( dx_tau_11 (i,j,k) + dy_tau_12 (i,j,k) - dx_pres (i,j,k) )
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

       call comm_one (pres)
       call comm_one (tau_11) ; call comm_one (tau_12) ; call comm_one (tau_13)
       call comm_one (tau_21) ; call comm_one (tau_22) ; call comm_one (tau_23)
       call dx ( dx_i , pres , dx_pres )     ; call dy ( dy_i , pres , dy_pres )
       call dx ( dx_i , tau_11 , dx_tau_11 ) ; call dy ( dy_i , tau_12 , dy_tau_12 ) ; call dz ( dz_i , tau_13 , dz_tau_13 )
       call dx ( dx_i , tau_21 , dx_tau_21 ) ; call dy ( dy_i , tau_22 , dy_tau_22 ) ; call dz ( dz_i , tau_23 , dz_tau_23 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ux (i,j,k) * ( dx_tau_21 (i,j,k) + dy_tau_22 (i,j,k) + dz_tau_23 (i,j,k) - dy_pres (i,j,k) ) + &
                                vy (i,j,k) * ( dx_tau_11 (i,j,k) + dy_tau_12 (i,j,k) + dz_tau_13 (i,j,k) - dx_pres (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( tau_11    , tau_12    , tau_13    , &
                 tau_21    , tau_22    , tau_23    , &
                 dx_tau_11 , dy_tau_12 , dz_tau_13 , &
                 dx_tau_21 , dy_tau_22 , dz_tau_23 , &
                 dx_pres   , dy_pres )


  end subroutine rey_mass_flux_coupl_12


!> \brief Calculate the 13-component of the mass-flux coupling term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_mass_flux_coupl_13 ( adi , dx_i , dy_i , dz_i , ux , wz , pres , tau , sigma )


    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i  !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux    !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: wz    !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: pres  !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau   !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: sigma !< component 13 of the mass-flux coupling term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: sqgmr
    real (dp) , dimension (:,:,:) , allocatable         :: tau_11    , tau_12    , tau_13    , &
                                                           tau_31    , tau_32    , tau_33    , &
                                                           dx_tau_11 , dy_tau_12 , dz_tau_13 , &
                                                           dx_tau_31 , dy_tau_32 , dz_tau_33 , &
                                                           dx_pres   , dz_pres


    allocate ( tau_11 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_12 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_13 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_31 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_32 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_33 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_tau_11 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_12 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_13 (sx:ex,sy:ey,sz:ez)                , &
               dx_tau_31 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_32 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_33 (sx:ex,sy:ey,sz:ez)                , &
               dx_pres   (sx:ex,sy:ey,sz:ez)                , &
               dz_pres   (sx:ex,sy:ey,sz:ez)                , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_mass_flux_coupl_13'


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       tau_11 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,1)
       tau_12 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,2)
       tau_13 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,1,3)
       tau_31 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,1)
       tau_32 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,2)
       tau_33 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,3)

       call comm_one (pres)
       call comm_one (tau_11) ; call comm_one (tau_12) ; call comm_one (tau_13)
       call comm_one (tau_31) ; call comm_one (tau_32) ; call comm_one (tau_33)
       call dx ( dx_i , pres , dx_pres )     ; call dz ( dz_i , pres , dz_pres )
       call dx ( dx_i , tau_11 , dx_tau_11 ) ; call dy ( dy_i , tau_12 , dy_tau_12 ) ; call dz ( dz_i , tau_13 , dz_tau_13 )
       call dx ( dx_i , tau_31 , dx_tau_31 ) ; call dy ( dy_i , tau_32 , dy_tau_32 ) ; call dz ( dz_i , tau_33 , dz_tau_33 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = ux (i,j,k) * ( dx_tau_31 (i,j,k) + dy_tau_32 (i,j,k) + dz_tau_33 (i,j,k) - dz_pres (i,j,k) ) + &
                                wz (i,j,k) * ( dx_tau_11 (i,j,k) + dy_tau_12 (i,j,k) + dz_tau_13 (i,j,k) - dx_pres (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( tau_11    , tau_12    , tau_13    , &
                 tau_31    , tau_32    , tau_33    , &
                 dx_tau_11 , dy_tau_12 , dz_tau_13 , &
                 dx_tau_31 , dy_tau_32 , dz_tau_33 , &
                 dx_pres   , dz_pres )


  end subroutine rey_mass_flux_coupl_13


!> \brief Calculate the 23-component of the mass-flux coupling term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rey_mass_flux_coupl_23 ( adi , dx_i , dy_i , dz_i , vy , wz , pres , tau , sigma )


    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i  !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vy    !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: wz    !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: pres  !< pressure
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau   !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: sigma !< component 23 of the mass-flux coupling term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: sqgmr
    real (dp) , dimension (:,:,:) , allocatable         :: tau_21    , tau_22    , tau_23    , &
                                                           tau_31    , tau_32    , tau_33    , &
                                                           dx_tau_21 , dy_tau_22 , dz_tau_23 , &
                                                           dx_tau_31 , dy_tau_32 , dz_tau_33 , &
                                                           dy_pres   , dz_pres


    allocate ( tau_21 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_22 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_23 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_31 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_32 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               tau_33 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_tau_21 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_22 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_23 (sx:ex,sy:ey,sz:ez)                , &
               dx_tau_31 (sx:ex,sy:ey,sz:ez)                , &
               dy_tau_32 (sx:ex,sy:ey,sz:ez)                , &
               dz_tau_33 (sx:ex,sy:ey,sz:ez)                , &
               dy_pres   (sx:ex,sy:ey,sz:ez)                , &
               dz_pres   (sx:ex,sy:ey,sz:ez)                , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate rey_mass_flux_coupl_23'


    sqgmr = adi % sqgmr


    if ( ndim == 2 ) then ! 2D problem

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = 0.0_dp
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       tau_21 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,1)
       tau_22 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,2)
       tau_23 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,2,3)
       tau_31 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,1)
       tau_32 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,2)
       tau_33 (sx:ex,sy:ey,sz:ez) = sqgmr * tau (sx:ex,sy:ey,sz:ez,3,3)

       call comm_one (pres)
       call comm_one (tau_21) ; call comm_one (tau_22) ; call comm_one (tau_23)
       call comm_one (tau_31) ; call comm_one (tau_32) ; call comm_one (tau_33)
       call dy ( dy_i , pres , dy_pres )     ; call dz ( dz_i , pres , dz_pres )
       call dx ( dx_i , tau_21 , dx_tau_21 ) ; call dy ( dy_i , tau_22 , dy_tau_22 ) ; call dz ( dz_i , tau_23 , dz_tau_23 )
       call dx ( dx_i , tau_31 , dx_tau_31 ) ; call dy ( dy_i , tau_32 , dy_tau_32 ) ; call dz ( dz_i , tau_33 , dz_tau_33 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                sigma (i,j,k) = vy (i,j,k) * ( dx_tau_31 (i,j,k) + dy_tau_32 (i,j,k) + dz_tau_33 (i,j,k) - dz_pres (i,j,k) ) + &
                                wz (i,j,k) * ( dx_tau_21 (i,j,k) + dy_tau_22 (i,j,k) + dz_tau_23 (i,j,k) - dy_pres (i,j,k) )
             end do
          end do
       end do

    end if


    deallocate ( tau_21    , tau_22    , tau_23    , &
                 tau_31    , tau_32    , tau_33    , &
                 dx_tau_21 , dy_tau_22 , dz_tau_23 , &
                 dx_tau_31 , dy_tau_32 , dz_tau_33 , &
                 dy_pres   , dz_pres )


  end subroutine rey_mass_flux_coupl_23


end module budget_Reynolds
