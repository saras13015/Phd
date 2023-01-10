!------------------------------------------------------------------------------
! MODULE: deriv
!------------------------------------------------------------------------------
!> \brief Derivatives definition.
!!
!! This module contains the definitions of derivatives, up to
!! 8th-order of accuracy, applied mainly for the viscous terms.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module deriv


  use parameters
  use parallel


  implicit none


  real (dp) , dimension (-4:4) , parameter , private :: &

       !> fourth-order one-side finite difference (LEFT)
       alph4ol = (/  0.0e0_dp               ,   0.0e0_dp              ,   0.0e0_dp              , &
                     0.0e0_dp               , -25.0e0_dp / 12.0e0_dp  ,  48.0e0_dp / 12.0e0_dp  , &
                   -36.0e0_dp / 12.0e0_dp   ,  16.0e0_dp / 12.0e0_dp  ,  -3.0e0_dp / 12.0e0_dp    /) , &

       !> fourth-order upwind finite difference (LEFT)
       alph4ul = (/   0.0e0_dp              ,   0.0e0_dp              ,   0.0e0_dp              , &
                     -3.0e0_dp / 12.0e0_dp  , -10.0e0_dp / 12.0e0_dp  ,  18.0e0_dp / 12.0e0_dp  , &
                     -6.0e0_dp / 12.0e0_dp  ,   1.0e0_dp / 12.0e0_dp  ,   0.0e0_dp                /) , &

       !> fourth-order central finite difference
       alph4   = (/   0.0e0_dp              ,   0.0e0_dp              ,   1.0e0_dp / 12.0e0_dp  , &
                     -8.0e0_dp / 12.0e0_dp  ,   0.0_dp                ,   8.0e0_dp / 12.0e0_dp  , &
                     -1.0e0_dp / 12.0e0_dp  ,   0.0e0_dp              ,   0.0e0_dp                /) , &

       !> sixth-order central finite difference
       alph6   = (/   0.0e0_dp              ,  -1.0e0_dp / 60.0e0_dp  ,   9.0e0_dp / 60.0e0_dp  , &
                    -45.0e0_dp / 60.0e0_dp  ,   0.0_dp                ,  45.0e0_dp / 60.0e0_dp  , &
                     -9.0e0_dp / 60.0e0_dp  ,   1.0e0_dp / 60.0e0_dp  ,   0.0e0_dp                /) , &

       !> eighth-order central finite difference
       alph8   = (/   3.0e0_dp / 840.0e0_dp , -32.0e0_dp / 840.0e0_dp , 168.0e0_dp / 840.0e0_dp , &
                   -672.0e0_dp / 840.0e0_dp ,   0.0_dp                , 672.0e0_dp / 840.0e0_dp , &
                   -168.0e0_dp / 840.0e0_dp ,  32.0e0_dp / 840.0e0_dp ,  -3.0e0_dp / 840.0e0_dp   /) , &

       !> fourth-order upwind finite difference (RIGHT)
       alph4ur = (/   0.0e0_dp              ,   1.0e0_dp / 12.0e0_dp  ,  -6.0e0_dp / 12.0e0_dp  , &
                     18.0e0_dp / 12.0e0_dp  , -10.0e0_dp / 12.0e0_dp  ,  -3.0e0_dp / 12.0e0_dp  , &
                      0.0e0_dp              ,   0.0e0_dp              ,   0.0e0_dp                /) , &

       !> fourth-order one-side finite difference (RIGHT)
       alph4or = (/ -3.0e0_dp / 12.0e0_dp   ,  16.0e0_dp / 12.0e0_dp  , -36.0e0_dp / 12.0e0_dp  , &
                    48.0e0_dp / 12.0e0_dp   , -25.0e0_dp / 12.0e0_dp  ,   0.0e0_dp              , &
                     0.0e0_dp               ,   0.0e0_dp              ,   0.0e0_dp              /)


contains


!> \brief Calculate the derivative of a generic array.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dscalar ( start , end , d_i , v , dv )


    integer (ip) , intent (in)                                     :: start , end !< start (sx, sy or sz) and end (ex, ey or ez) points of the subdomain
    real (dp) , allocatable , dimension (:) , intent (in)          :: d_i         !< 1/dx denominator
    real (dp) , allocatable , dimension (:) , intent (in)          :: v           !< generic array

    real (dp) , allocatable , dimension (:) , intent (inout)       :: dv          !< derivative of the generic array : v * (1/dx)


    integer (ip)  :: m , l

    real (dp)  :: dd


    ! eighth-order central finite difference
    do m = start+ng , end-ng
       dd = 0.0_dp
       do l = -4 , 4
          dd = dd + alph8(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
    end do

    ! fourth-order one-side finite difference (LEFT)
    do m = start , start
       dd = 0.0_dp
       do l = 0 , 4
          dd = dd + alph4ol(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
    end do

    ! fourth-order upwind finite difference (LEFT)
    do m = start+1 , start+1
       dd = 0.0_dp
       do l = -1 , 3
          dd = dd + alph4ul(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
    end do

    ! fourth-order central finite difference
    do m = start+2 , start+2
       dd = 0.0_dp
       do l = -2 , 2
          dd = dd + alph4(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
    end do

    ! sixth-order central finite difference
    do m = start+3 , start+3
       dd = 0.0_dp
       do l = -3 , 3
          dd = dd + alph6(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
    end do


    ! sixth-order central finite difference
    do m = end-3 , end-3
       dd = 0.0_dp
       do l = -3 , 3
          dd = dd + alph6(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
    end do

    ! fourth-order central finite difference
    do m = end-2 , end-2
       dd = 0.0_dp
       do l = -2 , 2
          dd = dd + alph4(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
    end do

    ! fourth-order upwind finite difference (RIGHT)
    do m = end-1 , end-1
       dd = 0.0_dp
       do l = -3 , 1
          dd = dd + alph4ur(l) * v ( m+l )
       end do
       dv (m) = - dd * d_i (m)
    end do

    ! fourth-order one-side finite difference (RIGHT)
    do m = end , end
       dd = 0.0_dp
       do l = -4 , 0
          dd = dd + alph4or(l) * v ( m+l )
       end do
       dv (m) = - dd * d_i (m)
    end do


  end subroutine dscalar


!> \brief Calculate a generic 3D array derivative in x-direction.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dx ( dx_i , v , dv )


    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i !< inverted dx array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: v    !< generic 3D array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: dv   !< derivative


    integer (ip)  :: i , j , k , l

    real (dp) :: dd


    ! common
    ! eighth-order central finite difference
    do k = sz , ez
       do j = sy , ey

          do i = sx+ng , ex-ng
             dd = 0.0_dp
             do l = -4 , 4
                dd = dd + alph8(l) * v ( i+l , j , k )
             end do
             dv (i,j,k) = dd * dx_i (i)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (W) == MPI_PROC_NULL ) then


       do k = sz , ez
          do j = sy , ey

             ! fourth-order one-side finite difference (LEFT)
             do i = sx , sx
                dd = 0.0_dp
                do l = 0 , 4
                   dd = dd + alph4ol(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order upwind finite difference (LEFT)
             do i = sx+1 , sx+1
                dd = 0.0_dp
                do l = -1 , 3
                   dd = dd + alph4ul(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order central finite difference
             do i = sx+2 , sx+2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! sixth-order central finite difference
             do i = sx+3 , sx+3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do j = sy , ey

             do i = sx , sx+ng-1
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    end if



    ! right boundary
    if ( neigh (E) == MPI_PROC_NULL ) then


       do k = sz , ez
          do j = sy , ey


             ! sixth-order central finite difference
             do i = ex-3 , ex-3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order central finite difference
             do i = ex-2 , ex-2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order upwind finite difference (RIGHT)
             do i = ex-1 , ex-1
                dd = 0.0_dp
                do l = -3 , 1
                   dd = dd + alph4ur(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = - dd * dx_i (i)
             end do
!
             ! fourth-order one-side finite difference (RIGHT)
             do i = ex , ex
                dd = 0.0_dp
                do l = -4 , 0
                   dd = dd + alph4or(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = - dd * dx_i (i)
             end do


          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do j = sy , ey

             do i = ex-ng+1 , ex
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    end if


  end subroutine dx


!> \brief Calculate a generic 3D array derivative in y-direction.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dy ( dy_i , v , dv )


    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i !< inverted dy array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: v    !< generic 3D array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: dv   !< derivative


    integer (ip)  :: i , j , k , l

    real (dp) :: dd


    ! common
    ! eighth-order central finite difference
    do k = sz , ez
       do i = sx , ex

          do j = sy+ng , ey-ng
             dd = 0.0_dp
             do l = -4 , 4
                dd = dd + alph8(l) * v ( i , j+l , k )
             end do
             dv (i,j,k) = dd * dy_i (j)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (S) == MPI_PROC_NULL ) then


       do k = sz , ez
          do i = sx , ex

             ! fourth-order one-side finite difference (LEFT)
             do j = sy , sy
                dd = 0.0_dp
                do l = 0 , 4
                   dd = dd + alph4ol(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order upwind finite difference (LEFT)
             do j = sy+1 , sy+1
                dd = 0.0_dp
                do l = -1 , 3
                   dd = dd + alph4ul(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order central finite difference
             do j = sy+2 , sy+2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! siyth-order central finite difference
             do j = sy+3 , sy+3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do i = sx , ex

             do j = sy , sy+ng-1
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    end if



    ! right boundary
    if ( neigh (N) == MPI_PROC_NULL ) then


       do k = sz , ez
          do i = sx , ex


             ! sixth-order central finite difference
             do j = ey-3 , ey-3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order central finite difference
             do j = ey-2 , ey-2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order upwind finite difference (RIGHT)
             do j = ey-1 , ey-1
                dd = 0.0_dp
                do l = -3 , 1
                   dd = dd + alph4ur(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = - dd * dy_i (j)
             end do
!
             ! fourth-order one-side finite difference (RIGHT)
             do j = ey , ey
                dd = 0.0_dp
                do l = -4 , 0
                   dd = dd + alph4or(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = - dd * dy_i (j)
             end do


          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do i = sx , ex

             do j = ey-ng+1 , ey
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    end if


  end subroutine dy


!> \brief Calculate a generic 3D array derivative in z-direction.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dz ( dz_i , v , dv )


    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: v    !< generic 3D array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: dv   !< derivative

    integer (ip)  :: i , j , k , l

    real (dp) :: dd


    ! common
    ! eighth-order central finite difference
    do j = sy , ey
       do i = sx , ex

          do k = sz+ng , ez-ng
             dd = 0.0_dp
             do l = -4 , 4
                dd = dd + alph8(l) * v ( i , j , k+l )
             end do
             dv (i,j,k) = dd * dz_i (k)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (B) == MPI_PROC_NULL ) then


       do j = sy , ey
          do i = sx , ex

             ! fourth-order one-side finite difference (LEFT)
             do k = sz , sz
                dd = 0.0_dp
                do l = 0 , 4
                   dd = dd + alph4ol(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order upwind finite difference (LEFT)
             do k = sz+1 , sz+1
                dd = 0.0_dp
                do l = -1 , 3
                   dd = dd + alph4ul(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order central finite difference
             do k = sz+2 , sz+2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! sixth-order central finite difference
             do k = sz+3 , sz+3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    else


       ! eighth-order central finite difference
       do j = sy , ey
          do i = sx , ex

             do k = sz , sz+ng-1
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    end if



    ! right boundary
    if ( neigh (F) == MPI_PROC_NULL ) then


       do j = sy , ey
          do i = sx , ex


             ! sixth-order central finite difference
             do k = ez-3 , ez-3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order central finite difference
             do k = ez-2 , ez-2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order upwind finite difference (RIGHT)
             do k = ez-1 , ez-1
                dd = 0.0_dp
                do l = -3 , 1
                   dd = dd + alph4ur(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = - dd * dz_i (k)
             end do
!
             ! fourth-order one-side finite difference (RIGHT)
             do k = ez , ez
                dd = 0.0_dp
                do l = -4 , 0
                   dd = dd + alph4or(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = - dd * dz_i (k)
             end do


          end do
       end do


    else


       ! eighth-order central finite difference
       do j = sy , ey
          do i = sx , ex

             do k = ez-ng+1 , ez
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    end if


  end subroutine dz


!> \brief Calculate a fixed (type 1) 3D array derivative in x-direction.
!!
!! This subroutine is exactly the same that the subroutine dx.
!! But it changes ONLY in dimension of the variables \f$ v \f$ and \f$ dv \f$ which have 
!! \f$ (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) \f$ dimension.
!! It must be followed by a com_deriv to calculate the second derivative by using dx_fixed2
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dx_fixed1 ( dx_i , v , dv )


    real (dp) , allocatable , dimension (:) , intent (in)                        :: dx_i !< inverted dx array
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (in)    :: v    !< fixed 3D array
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (inout) :: dv   !< fixed derivative


    integer (ip)  :: i , j , k , l

    real (dp) :: dd


    ! common
    ! eighth-order central finite difference
    do k = sz , ez
       do j = sy , ey

          do i = sx+ng , ex-ng
             dd = 0.0_dp
             do l = -4 , 4
                dd = dd + alph8(l) * v ( i+l , j , k )
             end do
             dv (i,j,k) = dd * dx_i (i)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (W) == MPI_PROC_NULL ) then


       do k = sz , ez
          do j = sy , ey

             ! fourth-order one-side finite difference (LEFT)
             do i = sx , sx
                dd = 0.0_dp
                do l = 0 , 4
                   dd = dd + alph4ol(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order upwind finite difference (LEFT)
             do i = sx+1 , sx+1
                dd = 0.0_dp
                do l = -1 , 3
                   dd = dd + alph4ul(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order central finite difference
             do i = sx+2 , sx+2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! sixth-order central finite difference
             do i = sx+3 , sx+3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do j = sy , ey

             do i = sx , sx+ng-1
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    end if



    ! right boundary
    if ( neigh (E) == MPI_PROC_NULL ) then


       do k = sz , ez
          do j = sy , ey


             ! sixth-order central finite difference
             do i = ex-3 , ex-3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order central finite difference
             do i = ex-2 , ex-2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order upwind finite difference (RIGHT)
             do i = ex-1 , ex-1
                dd = 0.0_dp
                do l = -3 , 1
                   dd = dd + alph4ur(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = - dd * dx_i (i)
             end do
!
             ! fourth-order one-side finite difference (RIGHT)
             do i = ex , ex
                dd = 0.0_dp
                do l = -4 , 0
                   dd = dd + alph4or(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = - dd * dx_i (i)
             end do


          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do j = sy , ey

             do i = ex-ng+1 , ex
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    end if


  end subroutine dx_fixed1


!> \brief Calculate a fixed (type 1) 3D array derivative in y-direction.
!!
!! This subroutine is exactly the same that the subroutine dx.
!! But it change ONLY in dimension of the variables \f$ v \f$ and \f$ dv \f$ which have 
!! \f$ (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) \f$ dimension.
!! It must be followed by a com_deriv to calculate the second derivative by using dy_fixed2
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dy_fixed1 ( dy_i , v , dv )


    real (dp) , allocatable , dimension (:) , intent (in)                        :: dy_i !< inverted dy array
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (in)    :: v    !< fixed 3D array
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (inout) :: dv   !< fixed derivative


    integer (ip)  :: i , j , k , l

    real (dp) :: dd


    ! common
    ! eighth-order central finite difference
    do k = sz , ez
       do i = sx , ex

          do j = sy+ng , ey-ng
             dd = 0.0_dp
             do l = -4 , 4
                dd = dd + alph8(l) * v ( i , j+l , k )
             end do
             dv (i,j,k) = dd * dy_i (j)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (S) == MPI_PROC_NULL ) then


       do k = sz , ez
          do i = sx , ex

             ! fourth-order one-side finite difference (LEFT)
             do j = sy , sy
                dd = 0.0_dp
                do l = 0 , 4
                   dd = dd + alph4ol(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order upwind finite difference (LEFT)
             do j = sy+1 , sy+1
                dd = 0.0_dp
                do l = -1 , 3
                   dd = dd + alph4ul(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order central finite difference
             do j = sy+2 , sy+2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! siyth-order central finite difference
             do j = sy+3 , sy+3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do i = sx , ex

             do j = sy , sy+ng-1
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    end if



    ! right boundary
    if ( neigh (N) == MPI_PROC_NULL ) then


       do k = sz , ez
          do i = sx , ex


             ! sixth-order central finite difference
             do j = ey-3 , ey-3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order central finite difference
             do j = ey-2 , ey-2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order upwind finite difference (RIGHT)
             do j = ey-1 , ey-1
                dd = 0.0_dp
                do l = -3 , 1
                   dd = dd + alph4ur(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = - dd * dy_i (j)
             end do
!
             ! fourth-order one-side finite difference (RIGHT)
             do j = ey , ey
                dd = 0.0_dp
                do l = -4 , 0
                   dd = dd + alph4or(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = - dd * dy_i (j)
             end do


          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do i = sx , ex

             do j = ey-ng+1 , ey
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    end if


  end subroutine dy_fixed1


!> \brief Calculate a fixed (type 1) 3D array derivative in z-direction.
!!
!! This subroutine is exactly the same that the subroutine dx.
!! But it change ONLY in dimension of the variables \f$ v \f$ and \f$ dv \f$ which have 
!! \f$ (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) \f$ dimension.
!! It must be followed by a com_deriv to calculate the second derivative by using dz_fixed2
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dz_fixed1 ( dz_i , v , dv )


    real (dp) , allocatable , dimension (:) , intent (in)                        :: dz_i !< inverted dz array
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (in)    :: v    !< fixed 3D array
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (inout) :: dv   !< fixed derivative


    integer (ip)  :: i , j , k , l

    real (dp) :: dd


    ! common
    ! eighth-order central finite difference
    do j = sy , ey
       do i = sx , ex

          do k = sz+ng , ez-ng
             dd = 0.0_dp
             do l = -4 , 4
                dd = dd + alph8(l) * v ( i , j , k+l )
             end do
             dv (i,j,k) = dd * dz_i (k)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (B) == MPI_PROC_NULL ) then


       do j = sy , ey
          do i = sx , ex

             ! fourth-order one-side finite difference (LEFT)
             do k = sz , sz
                dd = 0.0_dp
                do l = 0 , 4
                   dd = dd + alph4ol(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order upwind finite difference (LEFT)
             do k = sz+1 , sz+1
                dd = 0.0_dp
                do l = -1 , 3
                   dd = dd + alph4ul(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order central finite difference
             do k = sz+2 , sz+2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! sixth-order central finite difference
             do k = sz+3 , sz+3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    else


       ! eighth-order central finite difference
       do j = sy , ey
          do i = sx , ex

             do k = sz , sz+ng-1
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    end if



    ! right boundary
    if ( neigh (F) == MPI_PROC_NULL ) then


       do j = sy , ey
          do i = sx , ex


             ! sixth-order central finite difference
             do k = ez-3 , ez-3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order central finite difference
             do k = ez-2 , ez-2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order upwind finite difference (RIGHT)
             do k = ez-1 , ez-1
                dd = 0.0_dp
                do l = -3 , 1
                   dd = dd + alph4ur(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = - dd * dz_i (k)
             end do
!
             ! fourth-order one-side finite difference (RIGHT)
             do k = ez , ez
                dd = 0.0_dp
                do l = -4 , 0
                   dd = dd + alph4or(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = - dd * dz_i (k)
             end do


          end do
       end do


    else


       ! eighth-order central finite difference
       do j = sy , ey
          do i = sx , ex

             do k = ez-ng+1 , ez
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    end if


  end subroutine dz_fixed1


!> \brief Calculate a fixed (type 2) 3D array derivative in x-direction.
!!
!! This subroutine is exactly the same that the subroutine dx_fixed1.
!! But it changes ONLY in dimension of the variables \f$ v \f$ which has 
!! \f$ (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) \f$ dimension and \f$ dv \f$ which has 
!! \f$ (sx:ex,sy:ey,sz:ez) \f$  dimension.\n
!!
!! This subroutine is used to calculate the second derivative of viscous flux.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dx_fixed2 ( dx_i , v , dv )


    real (dp) , allocatable , dimension (:) , intent (in)                     :: dx_i !< inverted dx array
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (in) :: v    !< fixed 3D array
    real (dp) , dimension (sx:ex,sy:ey,sz:ez) , intent (inout)                :: dv   !< fixed derivative


    integer (ip)  :: i , j , k , l

    real (dp) :: dd


    ! common
    ! eighth-order central finite difference
    do k = sz , ez
       do j = sy , ey

          do i = sx+ng , ex-ng
             dd = 0.0_dp
             do l = -4 , 4
                dd = dd + alph8(l) * v ( i+l , j , k )
             end do
             dv (i,j,k) = dd * dx_i (i)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (W) == MPI_PROC_NULL ) then


       do k = sz , ez
          do j = sy , ey

             ! fourth-order one-side finite difference (LEFT)
             do i = sx , sx
                dd = 0.0_dp
                do l = 0 , 4
                   dd = dd + alph4ol(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order upwind finite difference (LEFT)
             do i = sx+1 , sx+1
                dd = 0.0_dp
                do l = -1 , 3
                   dd = dd + alph4ul(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order central finite difference
             do i = sx+2 , sx+2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! sixth-order central finite difference
             do i = sx+3 , sx+3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do j = sy , ey

             do i = sx , sx+ng-1
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    end if



    ! right boundary
    if ( neigh (E) == MPI_PROC_NULL ) then


       do k = sz , ez
          do j = sy , ey


             ! sixth-order central finite difference
             do i = ex-3 , ex-3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order central finite difference
             do i = ex-2 , ex-2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do
!
             ! fourth-order upwind finite difference (RIGHT)
             do i = ex-1 , ex-1
                dd = 0.0_dp
                do l = -3 , 1
                   dd = dd + alph4ur(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = - dd * dx_i (i)
             end do
!
             ! fourth-order one-side finite difference (RIGHT)
             do i = ex , ex
                dd = 0.0_dp
                do l = -4 , 0
                   dd = dd + alph4or(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = - dd * dx_i (i)
             end do


          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do j = sy , ey

             do i = ex-ng+1 , ex
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    end if


  end subroutine dx_fixed2


!> \brief Calculate a fixed (type 2) 3D array derivative in y-direction.
!!
!! This subroutine is exactly the same that the subroutine dx_fixed1.
!! But it changes ONLY in dimension of the variables \f$ v \f$ which has 
!! \f$ (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) \f$ dimension and \f$ dv \f$ which has 
!! \f$ (sx:ex,sy:ey,sz:ez) \f$  dimension.
!!
!! This subroutine is used to calculate the second derivative of viscous flux.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dy_fixed2 ( dy_i , v , dv )


    real (dp) , allocatable , dimension (:) , intent (in)                     :: dy_i !< inverted dy array
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (in) :: v    !< fixed 3D array
    real (dp) , dimension (sx:ex,sy:ey,sz:ez) , intent (inout)                :: dv   !< fixed derivative


    integer (ip)  :: i , j , k , l

    real (dp) :: dd


    ! common
    ! eighth-order central finite difference
    do k = sz , ez
       do i = sx , ex

          do j = sy+ng , ey-ng
             dd = 0.0_dp
             do l = -4 , 4
                dd = dd + alph8(l) * v ( i , j+l , k )
             end do
             dv (i,j,k) = dd * dy_i (j)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (S) == MPI_PROC_NULL ) then


       do k = sz , ez
          do i = sx , ex

             ! fourth-order one-side finite difference (LEFT)
             do j = sy , sy
                dd = 0.0_dp
                do l = 0 , 4
                   dd = dd + alph4ol(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order upwind finite difference (LEFT)
             do j = sy+1 , sy+1
                dd = 0.0_dp
                do l = -1 , 3
                   dd = dd + alph4ul(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order central finite difference
             do j = sy+2 , sy+2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! siyth-order central finite difference
             do j = sy+3 , sy+3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do i = sx , ex

             do j = sy , sy+ng-1
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    end if



    ! right boundary
    if ( neigh (N) == MPI_PROC_NULL ) then


       do k = sz , ez
          do i = sx , ex


             ! sixth-order central finite difference
             do j = ey-3 , ey-3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order central finite difference
             do j = ey-2 , ey-2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do
!
             ! fourth-order upwind finite difference (RIGHT)
             do j = ey-1 , ey-1
                dd = 0.0_dp
                do l = -3 , 1
                   dd = dd + alph4ur(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = - dd * dy_i (j)
             end do
!
             ! fourth-order one-side finite difference (RIGHT)
             do j = ey , ey
                dd = 0.0_dp
                do l = -4 , 0
                   dd = dd + alph4or(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = - dd * dy_i (j)
             end do


          end do
       end do


    else


       ! eighth-order central finite difference
       do k = sz , ez
          do i = sx , ex

             do j = ey-ng+1 , ey
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    end if


  end subroutine dy_fixed2


!> \brief Calculate a fixed (type 2) 3D array derivative in z-direction.
!!
!! This subroutine is exactly the same that the subroutine dx_fixed1.
!! But it changes ONLY in dimension of the variables \f$ v \f$ which has 
!! \f$ (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) \f$ dimension and \f$ dv \f$ which has 
!! \f$ (sx:ex,sy:ey,sz:ez) \f$  dimension.
!!
!! This subroutine is used to calculate the second derivative of viscous flux.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dz_fixed2 ( dz_i , v , dv )


    real (dp) , allocatable , dimension (:) , intent (in)                     :: dz_i !< inverted dz array
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (in) :: v    !< fixed 3D array
    real (dp) , dimension (sx:ex,sy:ey,sz:ez) , intent (inout)                :: dv   !< fixed derivative


    integer (ip)  :: i , j , k , l

    real (dp) :: dd


    ! common
    ! eighth-order central finite difference
    do j = sy , ey
       do i = sx , ex

          do k = sz+ng , ez-ng
             dd = 0.0_dp
             do l = -4 , 4
                dd = dd + alph8(l) * v ( i , j , k+l )
             end do
             dv (i,j,k) = dd * dz_i (k)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (B) == MPI_PROC_NULL ) then


       do j = sy , ey
          do i = sx , ex

             ! fourth-order one-side finite difference (LEFT)
             do k = sz , sz
                dd = 0.0_dp
                do l = 0 , 4
                   dd = dd + alph4ol(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order upwind finite difference (LEFT)
             do k = sz+1 , sz+1
                dd = 0.0_dp
                do l = -1 , 3
                   dd = dd + alph4ul(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order central finite difference
             do k = sz+2 , sz+2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! sixth-order central finite difference
             do k = sz+3 , sz+3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    else


       ! eighth-order central finite difference
       do j = sy , ey
          do i = sx , ex

             do k = sz , sz+ng-1
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    end if



    ! right boundary
    if ( neigh (F) == MPI_PROC_NULL ) then


       do j = sy , ey
          do i = sx , ex


             ! sixth-order central finite difference
             do k = ez-3 , ez-3
                dd = 0.0_dp
                do l = -3 , 3
                   dd = dd + alph6(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order central finite difference
             do k = ez-2 , ez-2
                dd = 0.0_dp
                do l = -2 , 2
                   dd = dd + alph4(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do
!
             ! fourth-order upwind finite difference (RIGHT)
             do k = ez-1 , ez-1
                dd = 0.0_dp
                do l = -3 , 1
                   dd = dd + alph4ur(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = - dd * dz_i (k)
             end do
!
             ! fourth-order one-side finite difference (RIGHT)
             do k = ez , ez
                dd = 0.0_dp
                do l = -4 , 0
                   dd = dd + alph4or(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = - dd * dz_i (k)
             end do


          end do
       end do


    else


       ! eighth-order central finite difference
       do j = sy , ey
          do i = sx , ex

             do k = ez-ng+1 , ez
                dd = 0.0_dp
                do l = -4 , 4
                   dd = dd + alph8(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    end if


  end subroutine dz_fixed2


end module deriv
