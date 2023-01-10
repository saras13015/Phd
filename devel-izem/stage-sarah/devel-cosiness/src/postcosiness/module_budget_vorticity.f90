!------------------------------------------------------------------------------
! MODULE: budget_vorticity
!------------------------------------------------------------------------------
!> \brief Budget equations for the vorticity transport equation.
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
module budget_vorticity


  use parameters

  use parallel

  use deriv


  implicit none


contains


!> \brief Calculate the vorticity.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine vorticity_i ( dx_i , dy_i , dz_i , v , vort_i )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i   !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i   !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i   !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v      !< velocity array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: vort_i !< vorticity array


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: rho_i
    real (dp) , dimension (:,:,:) , allocatable         :: ux , vy , wz
    real (dp) , dimension (:,:,:) , allocatable         :: dxvary , dxvarz , &
                                                           dyvarx , dyvarz , &
                                                           dzvarx , dzvary


    allocate ( ux     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxvary (sx:ex,sy:ey,sz:ez)                   , &
               dxvarz (sx:ex,sy:ey,sz:ez)                   , &
               dyvarx (sx:ex,sy:ey,sz:ez)                   , &
               dyvarz (sx:ex,sy:ey,sz:ez)                   , &
               dzvarx (sx:ex,sy:ey,sz:ez)                   , &
               dzvary (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate vorticity_i'


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
             end do
          end do
       end do

       call comm_one (ux) ; call comm_one (vy)

       call dx ( dx_i , vy , dxvary )
       call dy ( dy_i , ux , dyvarx )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                vort_i (i,j,k,1) = 0.0_dp
                vort_i (i,j,k,2) = 0.0_dp
                vort_i (i,j,k,3) = dxvary (i,j,k) - dyvarx (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


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

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( dx_i , vy , dxvary )
       call dx ( dx_i , wz , dxvarz )
       call dy ( dy_i , ux , dyvarx )
       call dy ( dy_i , wz , dyvarz )
       call dz ( dz_i , ux , dzvarx )
       call dz ( dz_i , vy , dzvary )


       ! this is the vorticity vector
       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                vort_i (i,j,k,1) = dyvarz (i,j,k) - dzvary (i,j,k)
                vort_i (i,j,k,2) = dzvarx (i,j,k) - dxvarz (i,j,k)
                vort_i (i,j,k,3) = dxvary (i,j,k) - dyvarx (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( ux , vy , wz , dxvary , dxvarz , dyvarx , dyvarz , dzvarx , dzvary )


  end subroutine vorticity_i


!> \brief Calculate the vorticity convection term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine vorticity_convection ( dx_i , dy_i , dz_i , v , vort_i , vort_conv )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: vort_i    !< vorticity array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vort_conv !< vorticity convection term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: rho_i , ux , vy , wz
    real (dp) , dimension (:,:,:) , allocatable         :: w3 , dxvortz , dyvortz , dzvortz


    allocate ( w3      (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxvortz (sx:ex,sy:ey,sz:ez)                   , &
               dyvortz (sx:ex,sy:ey,sz:ez)                   , &
               dzvortz (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate vorticity_convection'


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                w3 (i,j,k) = vort_i (i,j,k,3)
             end do
          end do
       end do

       call comm_one (w3)

       call dx ( dx_i , w3 , dxvortz )
       call dy ( dy_i , w3 , dyvortz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i             = 1.0_dp / v (i,j,k,1)
                ux                = v (i,j,k,2) * rho_i
                vy                = v (i,j,k,3) * rho_i
                vort_conv (i,j,k) = ux * dxvortz (i,j,k) + vy * dyvortz (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                w3 (i,j,k) = vort_i (i,j,k,3)
             end do
          end do
       end do

       call comm_one (w3)

       call dx ( dx_i , w3 , dxvortz )
       call dy ( dy_i , w3 , dyvortz )
       call dz ( dz_i , w3 , dzvortz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i             = 1.0_dp / v (i,j,k,1)
                ux                = v (i,j,k,2) * rho_i
                vy                = v (i,j,k,3) * rho_i
                wz                = v (i,j,k,4) * rho_i
                vort_conv (i,j,k) = ux * dxvortz (i,j,k) + vy * dyvortz (i,j,k) + wz * dzvortz (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( w3 , dxvortz , dyvortz , dzvortz )


  end subroutine vorticity_convection


!> \brief Calculate the enstrophy convection term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine enstrophy_convection ( dx_i , dy_i , dz_i , v , vort_i , enst_conv )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: vort_i    !< vorticity array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: enst_conv !< enstrophy convection term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: rho_i , ux , vy , wz
    real (dp) , dimension (:,:,:) , allocatable         :: w1     , w2     , w3     , &
                                                           dxvarx , dyvarx , dzvarx , &
                                                           dxvary , dyvary , dzvary , &
                                                           dxvarz , dyvarz , dzvarz


    allocate ( w1      (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               w2      (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               w3      (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxvarx (sx:ex,sy:ey,sz:ez)                    , &
               dyvarx (sx:ex,sy:ey,sz:ez)                    , &
               dzvarx (sx:ex,sy:ey,sz:ez)                    , &
               dxvary (sx:ex,sy:ey,sz:ez)                    , &
               dyvary (sx:ex,sy:ey,sz:ez)                    , &
               dzvary (sx:ex,sy:ey,sz:ez)                    , &
               dxvarz (sx:ex,sy:ey,sz:ez)                    , &
               dyvarz (sx:ex,sy:ey,sz:ez)                    , &
               dzvarz (sx:ex,sy:ey,sz:ez)                    , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate enstrophy_convection'


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                w3 (i,j,k) = vort_i (i,j,k,3)
             end do
          end do
       end do

       call comm_one (w3)

       call dx ( dx_i , w3 , dxvarz )
       call dy ( dy_i , w3 , dyvarz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i             = 1.0_dp / v (i,j,k,1)
                ux                = v (i,j,k,2) * rho_i
                vy                = v (i,j,k,3) * rho_i
                enst_conv (i,j,k) = ux * dxvarz (i,j,k) + vy * dyvarz (i,j,k)
                enst_conv (i,j,k) = enst_conv (i,j,k) * ( w3 (i,j,k) + w3 (i,j,k) )
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                w1 (i,j,k) = vort_i (i,j,k,1)
                w2 (i,j,k) = vort_i (i,j,k,2)
                w3 (i,j,k) = vort_i (i,j,k,3)
             end do
          end do
       end do

       call comm_one (w1) ; call comm_one (w2) ; call comm_one (w3)

       call dx ( dx_i , w1 , dxvarx ) ; call dx ( dx_i , w2 , dxvary ) ; call dx ( dx_i , w3 , dxvarz )
       call dy ( dy_i , w1 , dyvarx ) ; call dy ( dy_i , w2 , dyvary ) ; call dy ( dy_i , w3 , dyvarz )
       call dz ( dz_i , w1 , dzvarx ) ; call dz ( dz_i , w2 , dzvary ) ; call dz ( dz_i , w3 , dzvarz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                rho_i = 1.0_dp / v (i,j,k,1)
                ux    = v (i,j,k,2) * rho_i
                vy    = v (i,j,k,3) * rho_i
                wz    = v (i,j,k,4) * rho_i

                enst_conv (i,j,k) = w1 (i,j,k) * ( dxvarx (i,j,k) * ux + dyvarx (i,j,k) * vy + dzvarx (i,j,k) * wz ) + &
                                    w2 (i,j,k) * ( dxvary (i,j,k) * ux + dyvary (i,j,k) * vy + dzvary (i,j,k) * wz ) + &
                                    w3 (i,j,k) * ( dxvarz (i,j,k) * ux + dyvarz (i,j,k) * vy + dzvarz (i,j,k) * wz )

                enst_conv (i,j,k) = enst_conv (i,j,k) + enst_conv (i,j,k)

             end do
          end do
       end do


    end if


    deallocate ( w1     , w2     , w3     , &
                 dxvarx , dyvarx , dzvarx , &
                 dxvary , dyvary , dzvary , &
                 dxvarz , dyvarz , dzvarz )


  end subroutine enstrophy_convection


!> \brief Calculate the vorticity stretching term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine vorticity_stretching ( dx_i , dy_i , dz_i , v , vort_i , vort_stret )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i       !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i       !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i       !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v          !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: vort_i     !< vorticity array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vort_stret !< vorticity stretching term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: wz
    real (dp) , dimension (:,:,:) , allocatable         :: dxvarz , dyvarz , dzvarz


    allocate ( wz     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxvarz (sx:ex,sy:ey,sz:ez)                   , &
               dyvarz (sx:ex,sy:ey,sz:ez)                   , &
               dzvarz (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate vorticity_stretching'


    if ( ndim == 2 ) then ! 2D problem


       vort_stret (sx:ex,sy:ey,sz:ez) = 0.0_dp


    else if ( ndim == 3 ) then ! 3D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wz (i,j,k) = v (i,j,k,4) / v (i,j,k,1)
             end do
          end do
       end do

       call comm_one (wz)

       call dx ( dx_i , wz , dxvarz )
       call dy ( dy_i , wz , dyvarz )
       call dz ( dz_i , wz , dzvarz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                vort_stret (i,j,k) = vort_i (i,j,k,1) * dxvarz (i,j,k) + &
                                     vort_i (i,j,k,2) * dyvarz (i,j,k) + &
                                     vort_i (i,j,k,3) * dzvarz (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( wz , dxvarz , dyvarz , dzvarz )


  end subroutine vorticity_stretching


!> \brief Calculate the enstrophy stretching term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine enstrophy_stretching ( dx_i , dy_i , dz_i , v , vort_i , enst_stret )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i       !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i       !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i       !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v          !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: vort_i     !< vorticity array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: enst_stret !< enstrophy stretching term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: rho_i , w1 , w2 , w3
    real (dp) , dimension (:,:,:) , allocatable         :: ux , vy , wz
    real (dp) , dimension (:,:,:) , allocatable         :: dxvarx , dyvarx , dzvarx , &
                                                           dxvary , dyvary , dzvary , &
                                                           dxvarz , dyvarz , dzvarz


    allocate ( ux     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxvarx (sx:ex,sy:ey,sz:ez)                   , &
               dyvarx (sx:ex,sy:ey,sz:ez)                   , &
               dzvarx (sx:ex,sy:ey,sz:ez)                   , &
               dxvary (sx:ex,sy:ey,sz:ez)                   , &
               dyvary (sx:ex,sy:ey,sz:ez)                   , &
               dzvary (sx:ex,sy:ey,sz:ez)                   , &
               dxvarz (sx:ex,sy:ey,sz:ez)                   , &
               dyvarz (sx:ex,sy:ey,sz:ez)                   , &
               dzvarz (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate enstrophy_stretching'


    if ( ndim == 2 ) then ! 2D problem


       enst_stret (sx:ex,sy:ey,sz:ez) = 0.0_dp


    else if ( ndim == 3 ) then ! 3D problem


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

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( dx_i , ux , dxvarx ) ; call dx ( dx_i , vy , dxvary ) ; call dx ( dx_i , wz , dxvarz )
       call dy ( dy_i , ux , dyvarx ) ; call dy ( dy_i , vy , dyvary ) ; call dy ( dy_i , wz , dyvarz )
       call dz ( dz_i , ux , dzvarx ) ; call dz ( dz_i , vy , dzvary ) ; call dz ( dz_i , wz , dzvarz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                w1 = vort_i (i,j,k,1)
                w2 = vort_i (i,j,k,2)
                w3 = vort_i (i,j,k,3)

                enst_stret (i,j,k) = w1 * ( dxvarx (i,j,k) * w1 + dyvarx (i,j,k) * w2 + dzvarx (i,j,k) * w3 ) + &
                                     w2 * ( dxvary (i,j,k) * w1 + dyvary (i,j,k) * w2 + dzvary (i,j,k) * w3 ) + &
                                     w3 * ( dxvarz (i,j,k) * w1 + dyvarz (i,j,k) * w2 + dzvarz (i,j,k) * w3 )

                enst_stret (i,j,k) = enst_stret (i,j,k) + enst_stret (i,j,k)

             end do
          end do
       end do


    end if


    deallocate ( ux     , vy     , wz     , &
                 dxvarx , dyvarx , dzvarx , &
                 dxvary , dyvary , dzvary , &
                 dxvarz , dyvarz , dzvarz )


  end subroutine enstrophy_stretching


!> \brief Calculate the vorticity dilatation term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine vorticity_dilatation ( dx_i , dy_i , dz_i , v , vort_i , vort_dila )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: vort_i    !< vorticity array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vort_dila !< vorticity dilatation term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: rho_i
    real (dp) , dimension (:,:,:) , allocatable         :: ux , vy , wz
    real (dp) , dimension (:,:,:) , allocatable         :: dxvarx , dyvary , dzvarz


    allocate ( ux     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxvarx (sx:ex,sy:ey,sz:ez)                   , &
               dyvary (sx:ex,sy:ey,sz:ez)                   , &
               dzvarz (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate vorticity_dilatation'


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
             end do
          end do
       end do

       call comm_one (ux) ; call comm_one (vy)

       call dx ( dx_i , ux , dxvarx )
       call dy ( dy_i , vy , dyvary )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                vort_dila (i,j,k) = vort_i (i,j,k,3) * ( dxvarx (i,j,k) + dyvary (i,j,k) )
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


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

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( dx_i , ux , dxvarx )
       call dy ( dy_i , vy , dyvary )
       call dz ( dz_i , wz , dzvarz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                vort_dila (i,j,k) = vort_i (i,j,k,3) * ( dxvarx (i,j,k) + dyvary (i,j,k) + dzvarz (i,j,k) )
             end do
          end do
       end do


    end if


    deallocate ( ux , vy , wz , dxvarx , dyvary , dzvarz )


  end subroutine vorticity_dilatation


!> \brief Calculate the enstrophy dilatation term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine enstrophy_dilatation ( dx_i , dy_i , dz_i , v , vort_i , enst_dila )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i      !< inverted dx array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: vort_i    !< vorticity array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: enst_dila !< enstrophy dilatation term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: rho_i , w2
    real (dp) , dimension (:,:,:) , allocatable         :: ux , vy , wz
    real (dp) , dimension (:,:,:) , allocatable         :: dxvarx , dyvary , dzvarz


    allocate ( ux     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxvarx (sx:ex,sy:ey,sz:ez)                   , &
               dyvary (sx:ex,sy:ey,sz:ez)                   , &
               dzvarz (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate enstrophy_dilatation'


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
             end do
          end do
       end do

       call comm_one (ux) ; call comm_one (vy)

       call dx ( dx_i , ux , dxvarx )
       call dy ( dy_i , vy , dyvary )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                w2 = vort_i (i,j,k,3) * vort_i (i,j,k,3)
                enst_dila (i,j,k) = w2 * ( dxvarx (i,j,k) + dyvary (i,j,k) )
                enst_dila (i,j,k) = enst_dila (i,j,k) + enst_dila (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


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

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( dx_i , ux , dxvarx )
       call dy ( dy_i , vy , dyvary )
       call dz ( dz_i , wz , dzvarz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                w2 = vort_i (i,j,k,1) * vort_i (i,j,k,1) + &
                     vort_i (i,j,k,2) * vort_i (i,j,k,2) + &
                     vort_i (i,j,k,3) * vort_i (i,j,k,3)
                enst_dila (i,j,k) = w2 * ( dxvarx (i,j,k) + dyvary (i,j,k) + dzvarz (i,j,k) )
                enst_dila (i,j,k) = enst_dila (i,j,k) + enst_dila (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( ux , vy , wz , dxvarx , dyvary , dzvarz )


  end subroutine enstrophy_dilatation


!> \brief Calculate the vorticity baroclininc term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine vorticity_baroclinic ( dx_i , dy_i , W_i , T , v , vort_baro )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: W_i       !< inverted molar mass
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: T         !< temperature
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v         !< conserved variables array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vort_baro !< vorticity baroclinic term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: rho , pres
    real (dp) , dimension (:,:,:) , allocatable         :: dxvarrho  , dyvarrho  , &
                                                           dxvarpres , dyvarpres


    allocate ( rho       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               pres      (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxvarrho  (sx:ex,sy:ey,sz:ez)                   , &
               dyvarrho  (sx:ex,sy:ey,sz:ez)                   , &
               dxvarpres (sx:ex,sy:ey,sz:ez)                   , &
               dyvarpres (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate vorticity_baroclinic'


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho (i,j,k)  = v (i,j,k,1)
                pres (i,j,k) = rho (i,j,k) * T (i,j,k) * W_i (i,j,k)
             end do
          end do
       end do

       call comm_one (rho) ; call comm_one (pres)

       call dx ( dx_i , rho , dxvarrho )   ; call dy ( dy_i , rho , dyvarrho )
       call dx ( dx_i , pres , dxvarpres ) ; call dy ( dy_i , pres , dyvarpres )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                vort_baro (i,j,k) = dxvarrho (i,j,k) * dyvarpres (i,j,k) - &
                                    dyvarrho (i,j,k) * dxvarpres (i,j,k)
                vort_baro (i,j,k) = vort_baro (i,j,k) / ( rho (i,j,k) * rho (i,j,k) )
             end do
          end do
       end do


    deallocate ( rho , pres , dxvarrho , dyvarrho , dxvarpres , dyvarpres )


  end subroutine vorticity_baroclinic


!> \brief Calculate the enstrophy baroclininc term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine enstrophy_baroclinic ( dx_i , dy_i , dz_i , W_i , T , v , vort_i , enst_baro )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: W_i       !< inverted molar mass
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: T         !< temperature
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: vort_i    !< vorticity array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: enst_baro !< enstrophy baroclinic term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: rho , pres
    real (dp) , dimension (:,:,:) , allocatable         :: dxvarrho  , dyvarrho  , dzvarrho  , &
                                                           dxvarpres , dyvarpres , dzvarpres


    allocate ( rho       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               pres      (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxvarrho  (sx:ex,sy:ey,sz:ez)                   , &
               dyvarrho  (sx:ex,sy:ey,sz:ez)                   , &
               dzvarrho  (sx:ex,sy:ey,sz:ez)                   , &
               dxvarpres (sx:ex,sy:ey,sz:ez)                   , &
               dyvarpres (sx:ex,sy:ey,sz:ez)                   , &
               dzvarpres (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate enstrophy_baroclinic'


    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             rho (i,j,k)  = v (i,j,k,1)
             pres (i,j,k) = rho (i,j,k) * T (i,j,k) * W_i (i,j,k)
          end do
       end do
    end do

    call comm_one (rho) ; call comm_one (pres)


    if ( ndim == 2 ) then ! 2D problem


       call dx ( dx_i , rho , dxvarrho )   ; call dy ( dy_i , rho , dyvarrho )
       call dx ( dx_i , pres , dxvarpres ) ; call dy ( dy_i , pres , dyvarpres )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                enst_baro (i,j,k) = ( dxvarrho (i,j,k) * dyvarpres (i,j,k) -   &
                                      dyvarrho (i,j,k) * dxvarpres (i,j,k) ) * &
                                      vort_i (i,j,k,3)
                enst_baro (i,j,k) = enst_baro (i,j,k) / ( rho (i,j,k) * rho (i,j,k) )
                enst_baro (i,j,k) = enst_baro (i,j,k) + enst_baro (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       call dx ( dx_i , rho , dxvarrho )   ; call dy ( dy_i , rho , dyvarrho )   ; call dz ( dz_i , rho , dzvarrho )
       call dx ( dx_i , pres , dxvarpres ) ; call dy ( dy_i , pres , dyvarpres ) ; call dz ( dz_i , pres , dzvarpres )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                enst_baro (i,j,k) = ( dyvarrho (i,j,k) * dzvarpres (i,j,k) - dzvarrho (i,j,k) * dyvarpres (i,j,k) ) * &
                                      vort_i (i,j,k,1)                                                              + &
                                    ( dzvarrho (i,j,k) * dxvarpres (i,j,k) - dxvarrho (i,j,k) * dzvarpres (i,j,k) ) * &
                                      vort_i (i,j,k,2)                                                              + &
                                    ( dxvarrho (i,j,k) * dyvarpres (i,j,k) - dyvarrho (i,j,k) * dxvarpres (i,j,k) ) * &
                                      vort_i (i,j,k,3)

                enst_baro (i,j,k) = enst_baro (i,j,k) / ( rho (i,j,k) * rho (i,j,k) )

                enst_baro (i,j,k) = enst_baro (i,j,k) + enst_baro (i,j,k)

             end do
          end do
       end do


    end if


    deallocate ( rho , pres , dxvarrho , dyvarrho , dzvarrho , dxvarpres , dyvarpres , dzvarpres )


  end subroutine enstrophy_baroclinic


!> \brief Calculate the vorticity viscous term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine vorticity_viscous ( dx_i , dy_i , dz_i , v , tau , vort_visc )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau       !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vort_visc !< vorticity viscous term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: rho_i
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3
    real (dp) , dimension (:,:,:) , allocatable         :: dxtau11 , dytau12 , dztau13 , &
                                                           dxtau21 , dytau22 , dztau23


    allocate ( wrk1    (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk2    (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk3    (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxtau11 (sx:ex,sy:ey,sz:ez)                   , &
               dytau12 (sx:ex,sy:ey,sz:ez)                   , &
               dztau13 (sx:ex,sy:ey,sz:ez)                   , &
               dxtau21 (sx:ex,sy:ey,sz:ez)                   , &
               dytau22 (sx:ex,sy:ey,sz:ez)                   , &
               dztau23 (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate vorticity_viscous'


    if ( ndim == 2 ) then ! 2D problem


       wrk1 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,2)

       call comm_one (wrk1) ; call comm_one (wrk2)

       call dx ( dx_i , wrk1 , dxtau11 )
       call dy ( dy_i , wrk2 , dytau12 )

       wrk1 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,2)

       call comm_one (wrk1) ; call comm_one (wrk2)

       call dx ( dx_i , wrk1 , dxtau21 )
       call dy ( dy_i , wrk2 , dytau22 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i        = 1.0_dp / v (i,j,k,1)
                wrk1 (i,j,k) = rho_i * ( dxtau11 (i,j,k) + dytau12 (i,j,k) )
                wrk2 (i,j,k) = rho_i * ( dxtau21 (i,j,k) + dytau22 (i,j,k) )
             end do
          end do
       end do

       call comm_one (wrk1) ; call comm_one (wrk2)

       call dx ( dx_i , wrk2 , dxtau11 )
       call dy ( dy_i , wrk1 , dytau12 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                vort_visc (i,j,k) = dxtau11 (i,j,k) - dytau12 (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       wrk1 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,2)
       wrk3 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,3)

       call comm_one (wrk1) ; call comm_one (wrk2) ; call comm_one (wrk3)

       call dx ( dx_i , wrk1 , dxtau11 )
       call dy ( dy_i , wrk2 , dytau12 )
       call dz ( dz_i , wrk3 , dztau13 )

       wrk1 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,2)
       wrk3 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,3)

       call comm_one (wrk1) ; call comm_one (wrk2) ; call comm_one (wrk3)

       call dx ( dx_i , wrk1 , dxtau21 )
       call dy ( dy_i , wrk2 , dytau22 )
       call dz ( dz_i , wrk3 , dztau23 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i        = 1.0_dp / v (i,j,k,1)
                wrk1 (i,j,k) = rho_i * ( dxtau11 (i,j,k) + dytau12 (i,j,k) + dztau13 (i,j,k) )
                wrk2 (i,j,k) = rho_i * ( dxtau21 (i,j,k) + dytau22 (i,j,k) + dztau23 (i,j,k) )
             end do
          end do
       end do

       call comm_one (wrk1) ; call comm_one (wrk2)

       call dx ( dx_i , wrk2 , dxtau11 )
       call dy ( dy_i , wrk1 , dytau12 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                vort_visc (i,j,k) = dxtau11 (i,j,k) - dytau12 (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( wrk1    , wrk2    , wrk3    , &
                 dxtau21 , dytau22 , dztau23 , &
                 dxtau11 , dytau12 , dztau13 )


  end subroutine vorticity_viscous


!> \brief Calculate the enstrophy viscous term.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine enstrophy_viscous ( dx_i , dy_i , dz_i , v , vort_i , tau , enst_visc )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v         !< conserved variables array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: vort_i    !< vorticity array
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: tau       !< Reynolds stress tensor
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: enst_visc !< enstrophy viscous term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: rho_i
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3
    real (dp) , dimension (:,:,:) , allocatable         :: dxtau11 , dytau12 , dztau13 , &
                                                           dxtau21 , dytau22 , dztau23 , &
                                                           dxtau31 , dytau32 , dztau33


    allocate ( wrk1    (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk2    (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk3    (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dxtau21 (sx:ex,sy:ey,sz:ez)                   , &
               dytau22 (sx:ex,sy:ey,sz:ez)                   , &
               dztau23 (sx:ex,sy:ey,sz:ez)                   , &
               dxtau11 (sx:ex,sy:ey,sz:ez)                   , &
               dytau12 (sx:ex,sy:ey,sz:ez)                   , &
               dztau13 (sx:ex,sy:ey,sz:ez)                   , &
               dxtau31 (sx:ex,sy:ey,sz:ez)                   , &
               dytau32 (sx:ex,sy:ey,sz:ez)                   , &
               dztau33 (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate vorticity_viscous'


    if ( ndim == 2 ) then ! 2D problem


       wrk1 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,2)

       call comm_one (wrk1) ; call comm_one (wrk2)

       call dx ( dx_i , wrk1 , dxtau11 )
       call dy ( dy_i , wrk2 , dytau12 )

       wrk1 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,2)

       call comm_one (wrk1) ; call comm_one (wrk2)

       call dx ( dx_i , wrk1 , dxtau21 )
       call dy ( dy_i , wrk2 , dytau22 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i        = 1.0_dp / v (i,j,k,1)
                wrk1 (i,j,k) = rho_i * ( dxtau11 (i,j,k) + dytau12 (i,j,k) )
                wrk2 (i,j,k) = rho_i * ( dxtau21 (i,j,k) + dytau22 (i,j,k) )
             end do
          end do
       end do

       call comm_one (wrk1) ; call comm_one (wrk2)

       call dx ( dx_i , wrk2 , dxtau11 )
       call dy ( dy_i , wrk1 , dytau12 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                enst_visc (i,j,k) = ( vort_i (i,j,k,3) + vort_i (i,j,k,3) ) * &
                                    ( dxtau11 (i,j,k) - dytau12 (i,j,k) )
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       wrk1 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,2)
       wrk3 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,1,3)

       call comm_one (wrk1) ; call comm_one (wrk2) ; call comm_one (wrk3)

       call dx ( dx_i , wrk1 , dxtau11 )
       call dy ( dy_i , wrk2 , dytau12 )
       call dz ( dz_i , wrk3 , dztau13 )

       wrk1 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,2)
       wrk3 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,2,3)

       call comm_one (wrk1) ; call comm_one (wrk2) ; call comm_one (wrk3)

       call dx ( dx_i , wrk1 , dxtau21 )
       call dy ( dy_i , wrk2 , dytau22 )
       call dz ( dz_i , wrk3 , dztau23 )

       wrk1 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,3,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,3,2)
       wrk3 (sx:ex,sy:ey,sz:ez) = tau (sx:ex,sy:ey,sz:ez,3,3)

       call comm_one (wrk1) ; call comm_one (wrk2) ; call comm_one (wrk3)

       call dx ( dx_i , wrk1 , dxtau31 )
       call dy ( dy_i , wrk2 , dytau32 )
       call dz ( dz_i , wrk3 , dztau33 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i        = 1.0_dp / v (i,j,k,1)
                wrk1 (i,j,k) = rho_i * ( dxtau11 (i,j,k) + dytau12 (i,j,k) + dztau13 (i,j,k) )
                wrk2 (i,j,k) = rho_i * ( dxtau21 (i,j,k) + dytau22 (i,j,k) + dztau23 (i,j,k) )
                wrk3 (i,j,k) = rho_i * ( dxtau31 (i,j,k) + dytau32 (i,j,k) + dztau33 (i,j,k) )
             end do
          end do
       end do

       call comm_one (wrk1) ; call comm_one (wrk2) ; call comm_one (wrk3)

       call dy ( dy_i , wrk3 , dxtau31 ) ! +
       call dz ( dz_i , wrk2 , dytau32 ) ! -

       call dx ( dx_i , wrk3 , dxtau21 ) ! -
       call dz ( dz_i , wrk1 , dytau22 ) ! +

       call dx ( dx_i , wrk2 , dxtau11 ) ! +
       call dy ( dy_i , wrk1 , dytau12 ) ! -

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                enst_visc (i,j,k) = vort_i (i,j,k,1)                    * &
                                  ( dxtau31 (i,j,k) - dytau32 (i,j,k) ) + &
                                    vort_i (i,j,k,2)                    * &
                                  ( dytau22 (i,j,k) - dxtau21 (i,j,k) ) + &
                                    vort_i (i,j,k,3)                    * &
                                  ( dxtau11 (i,j,k) - dytau12 (i,j,k) )
                enst_visc (i,j,k) = enst_visc (i,j,k) + enst_visc (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( wrk1    , wrk2    , wrk3    , &
                 dxtau11 , dytau12 , dztau13 , &
                 dxtau21 , dytau22 , dztau23 , &
                 dxtau31 , dytau32 , dztau33 )


  end subroutine enstrophy_viscous


end module budget_vorticity
