!------------------------------------------------------------------------------
! MODULE: deriv_post
!------------------------------------------------------------------------------
!> \brief Derivatives definition forced to low accuracy.
!!
!! This module is a modified version of the original module_deriv.f90.
!! The interest of such module, refered as _deriv_ and not
!! _deriv_post_, is to use the same functionality as the original
!! one. Nevertheless, here only 2nd-order derivatives are evaluated
!! for post-treatment calculations.
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


  real (dp) , dimension (-1:1) , parameter , private :: &

       !> second-order centered finite difference
       alph2c = (/ -0.5_dp , 0.0_dp , 0.5_dp /) , &

       !> first-order forward finite difference
       alph1f = (/ 0.0_dp , -1.0_dp , 1.0_dp /) , &

       !> first-order backward finite difference
       alph1b = (/ -1.0_dp , 1.0_dp , 0.0_dp /)


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

    real (dp) :: dd


    ! second-order centered finite difference
    do m = start+1 , end-1
       dd = 0.0_dp
       do l = -1 , 1
          dd = dd + alph2c(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
    end do

    ! first-order forward finite difference (LEFT)
    do m = start , start
       dd = 0.0_dp
       do l = 0 , 1
          dd = dd + alph1f(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
    end do

    ! first-order backward finite difference
    do m = end , end
       dd = 0.0_dp
       do l = -1 , 0
          dd = dd + alph1b(l) * v ( m+l )
       end do
       dv (m) = dd * d_i (m)
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
    ! second-order centered finite difference
    do k = sz , ez
       do j = sy , ey

          do i = sx+1 , ex-1
             dd = 0.0_dp
             do l = -1 , 1
                dd = dd + alph2c(l) * v ( i+l , j , k )
             end do
             dv (i,j,k) = dd * dx_i (i)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (W) == MPI_PROC_NULL ) then


       do k = sz , ez
          do j = sy , ey

             ! first-order forward finite difference (LEFT)
             do i = sx , sx
                dd = 0.0_dp
                do l = 0 , 1
                   dd = dd + alph1f(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do j = sy , ey

             do i = sx , sx
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i+l , j , k )
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


             ! first-order backward finite difference
             do i = ex , ex
                dd = 0.0_dp
                do l = -1 , 0
                   dd = dd + alph1b(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do


          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do j = sy , ey

             do i = ex , ex
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i+l , j , k )
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
    ! second-order centered finite difference
    do k = sz , ez
       do i = sx , ex

          do j = sy+1 , ey-1
             dd = 0.0_dp
             do l = -1 , 1
                dd = dd + alph2c(l) * v ( i , j+l , k )
             end do
             dv (i,j,k) = dd * dy_i (j)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (S) == MPI_PROC_NULL ) then


       do k = sz , ez
          do i = sx , ex

             ! first-order forward finite difference (LEFT)
             do j = sy , sy
                dd = 0.0_dp
                do l = 0 , 1
                   dd = dd + alph1f(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do i = sx , ex

             do j = sy , sy
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j+l , k )
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


             ! fisrt-order backward finite difference
             do j = ey , ey
                dd = 0.0_dp
                do l = -1 , 0
                   dd = dd + alph1b(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do


          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do i = sx , ex

             do j = ey , ey
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j+l , k )
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
    ! second-order centered finite difference
    do j = sy , ey
       do i = sx , ex

          do k = sz+1 , ez-1
             dd = 0.0_dp
             do l = -1 , 1
                dd = dd + alph2c(l) * v ( i , j , k+l )
             end do
             dv (i,j,k) = dd * dz_i (k)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (B) == MPI_PROC_NULL ) then


       do j = sy , ey
          do i = sx , ex

             ! first-order forward finite difference (LEFT)
             do k = sz , sz
                dd = 0.0_dp
                do l = 0 , 1
                   dd = dd + alph1f(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    else


       ! second-order centered finite difference
       do j = sy , ey
          do i = sx , ex

             do k = sz , sz
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j , k+l )
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


             ! first-order backward finite difference
             do k = ez , ez
                dd = 0.0_dp
                do l = -1 , 0
                   dd = dd + alph1b(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do


          end do
       end do


    else


       ! second-order centered finite difference
       do j = sy , ey
          do i = sx , ex

             do k = ez , ez
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    end if


  end subroutine dz


!> \brief Calculate a fixed (type 1) 3D array derivative in x-direction.
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
    ! second-order centered finite difference
    do k = sz , ez
       do j = sy , ey

          do i = sx+1 , ex-1
             dd = 0.0_dp
             do l = -1 , 1
                dd = dd + alph2c(l) * v ( i+l , j , k )
             end do
             dv (i,j,k) = dd * dx_i (i)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (W) == MPI_PROC_NULL ) then


       do k = sz , ez
          do j = sy , ey

             ! first-order forward finite difference (LEFT)
             do i = sx , sx
                dd = 0.0_dp
                do l = 0 , 1
                   dd = dd + alph1f(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do j = sy , ey

             do i = sx , sx
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i+l , j , k )
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


             ! first-order backward finite difference
             do i = ex , ex
                dd = 0.0_dp
                do l = -1 , 0
                   dd = dd + alph1b(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do


          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do j = sy , ey

             do i = ex , ex
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    end if


  end subroutine dx_fixed1


!> \brief Calculate a fixed (type 1) 3D array derivative in y-direction.
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
    ! second-order centered finite difference
    do k = sz , ez
       do i = sx , ex

          do j = sy+1 , ey-1
             dd = 0.0_dp
             do l = -1 , 1
                dd = dd + alph2c(l) * v ( i , j+l , k )
             end do
             dv (i,j,k) = dd * dy_i (j)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (S) == MPI_PROC_NULL ) then


       do k = sz , ez
          do i = sx , ex

             ! first-order forward finite difference (LEFT)
             do j = sy , sy
                dd = 0.0_dp
                do l = 0 , 1
                   dd = dd + alph1f(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do i = sx , ex

             do j = sy , sy
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j+l , k )
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


             ! fisrt-order backward finite difference
             do j = ey , ey
                dd = 0.0_dp
                do l = -1 , 0
                   dd = dd + alph1b(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do


          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do i = sx , ex

             do j = ey , ey
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    end if


  end subroutine dy_fixed1


!> \brief Calculate a fixed (type 1) 3D array derivative in z-direction.
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
    ! second-order centered finite difference
    do j = sy , ey
       do i = sx , ex

          do k = sz+1 , ez-1
             dd = 0.0_dp
             do l = -1 , 1
                dd = dd + alph2c(l) * v ( i , j , k+l )
             end do
             dv (i,j,k) = dd * dz_i (k)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (B) == MPI_PROC_NULL ) then


       do j = sy , ey
          do i = sx , ex

             ! first-order forward finite difference (LEFT)
             do k = sz , sz
                dd = 0.0_dp
                do l = 0 , 1
                   dd = dd + alph1f(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    else


       ! second-order centered finite difference
       do j = sy , ey
          do i = sx , ex

             do k = sz , sz
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j , k+l )
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


             ! first-order backward finite difference
             do k = ez , ez
                dd = 0.0_dp
                do l = -1 , 0
                   dd = dd + alph1b(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do


          end do
       end do


    else


       ! second-order centered finite difference
       do j = sy , ey
          do i = sx , ex

             do k = ez , ez
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    end if


  end subroutine dz_fixed1


!> \brief Calculate a fixed (type 2) 3D array derivative in x-direction.
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
    ! second-order centered finite difference
    do k = sz , ez
       do j = sy , ey

          do i = sx+1 , ex-1
             dd = 0.0_dp
             do l = -1 , 1
                dd = dd + alph2c(l) * v ( i+l , j , k )
             end do
             dv (i,j,k) = dd * dx_i (i)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (W) == MPI_PROC_NULL ) then


       do k = sz , ez
          do j = sy , ey

             ! first-order forward finite difference (LEFT)
             do i = sx , sx
                dd = 0.0_dp
                do l = 0 , 1
                   dd = dd + alph1f(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do j = sy , ey

             do i = sx , sx
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i+l , j , k )
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


             ! first-order backward finite difference
             do i = ex , ex
                dd = 0.0_dp
                do l = -1 , 0
                   dd = dd + alph1b(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do


          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do j = sy , ey

             do i = ex , ex
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i+l , j , k )
                end do
                dv (i,j,k) = dd * dx_i (i)
             end do

          end do
       end do


    end if


  end subroutine dx_fixed2


!> \brief Calculate a fixed (type 2) 3D array derivative in y-direction.
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
    ! second-order centered finite difference
    do k = sz , ez
       do i = sx , ex

          do j = sy+1 , ey-1
             dd = 0.0_dp
             do l = -1 , 1
                dd = dd + alph2c(l) * v ( i , j+l , k )
             end do
             dv (i,j,k) = dd * dy_i (j)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (S) == MPI_PROC_NULL ) then


       do k = sz , ez
          do i = sx , ex

             ! first-order forward finite difference (LEFT)
             do j = sy , sy
                dd = 0.0_dp
                do l = 0 , 1
                   dd = dd + alph1f(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do i = sx , ex

             do j = sy , sy
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j+l , k )
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


             ! fisrt-order backward finite difference
             do j = ey , ey
                dd = 0.0_dp
                do l = -1 , 0
                   dd = dd + alph1b(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do


          end do
       end do


    else


       ! second-order centered finite difference
       do k = sz , ez
          do i = sx , ex

             do j = ey , ey
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j+l , k )
                end do
                dv (i,j,k) = dd * dy_i (j)
             end do

          end do
       end do


    end if


  end subroutine dy_fixed2


!> \brief Calculate a fixed (type 2) 3D array derivative in z-direction.
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
    ! second-order centered finite difference
    do j = sy , ey
       do i = sx , ex

          do k = sz+1 , ez-1
             dd = 0.0_dp
             do l = -1 , 1
                dd = dd + alph2c(l) * v ( i , j , k+l )
             end do
             dv (i,j,k) = dd * dz_i (k)
          end do

       end do
    end do


    ! left boundary
    if ( neigh (B) == MPI_PROC_NULL ) then


       do j = sy , ey
          do i = sx , ex

             ! first-order forward finite difference (LEFT)
             do k = sz , sz
                dd = 0.0_dp
                do l = 0 , 1
                   dd = dd + alph1f(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    else


       ! second-order centered finite difference
       do j = sy , ey
          do i = sx , ex

             do k = sz , sz
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j , k+l )
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


             ! first-order backward finite difference
             do k = ez , ez
                dd = 0.0_dp
                do l = -1 , 0
                   dd = dd + alph1b(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do


          end do
       end do


    else


       ! second-order centered finite difference
       do j = sy , ey
          do i = sx , ex

             do k = ez , ez
                dd = 0.0_dp
                do l = -1 , 1
                   dd = dd + alph2c(l) * v ( i , j , k+l )
                end do
                dv (i,j,k) = dd * dz_i (k)
             end do

          end do
       end do


    end if


  end subroutine dz_fixed2


end module deriv
