!------------------------------------------------------------------------------
! MODULE: BCs
!------------------------------------------------------------------------------
!> \brief Boundary conditions.
!!
!! This module contains the definition of all the boundary conditions.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module BCs

  use parameters
  use parallel
  use input
  use adim
  use thermodynamics
  use weno
  use Rankine_Hugoniot

  implicit none

  real (dp) , parameter , private :: pi2 = pi + pi                      , &
                                     sqrt27_i = 1.0_dp / sqrt (27.0_dp)


contains


!> \brief Selector to update the boundary.
!!
!! This subroutine reads the keyword to select the corresponding boundary condition.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine upd_boundaries ( inp , adi , thd , time , dtime , grid , T , W_i , cp , ha , v )


    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type
    real (dp) , intent (in)                                        :: time  !< time
    real (dp) , intent (in)                                        :: dtime !< time step
    type (inp_grid) , intent (in)                                  :: grid  !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i   !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp    !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha    !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v     !< conserved variables array


    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -10 , 10 , -20 , 20 , -30 , 30 /)

    integer (ip) :: l

    logical      :: dirichlet


    do l = 1 , ndim+ndim

       dirichlet = .false.

       if ( neigh (l) == MPI_PROC_NULL ) then

          if ( bc (l) == extrapolation ) then
             call bc_extrapolation ( l , v )
          else if ( bc (l) == noreflection ) then
             call bc_noreflection ( l , v )
          else if ( bc (l) == symmetryplane ) then
             dirichlet = .true.
             call bc_symmetryplane ( l , v )
          else if ( bc (l) == adiabaticwall ) then
             dirichlet = .true.
             call bc_adiabaticwall ( l , v )
          else if ( bc (l) == syntheticturb ) then
             dirichlet = .true.
             call bc_syntheticturb ( l , inp , adi , thd , grid % y , grid % z , v )
          else if ( bc (l) == shock ) then
             dirichlet = .true.
             call bc_shock ( l , inp , adi , thd , v )
          else if ( bc (l) == postshock_slipwall ) then
             dirichlet = .true.
             call bc_postshock_slipwall ( l , adi , thd , grid % x , v )
          else if ( bc (l) == inflowcheng ) then
             dirichlet = .true.
             call bc_inflowcheng ( l , adi , thd , grid % y , v )
          else if ( bc (l) == inflowmiller ) then
             dirichlet = .true.
             call bc_inflowmiller ( l , adi , thd , grid % y , v )
          else if ( bc (l) == inflowmixchengmiller ) then
             dirichlet = .true.
             call bc_inflowmixchengmiller ( l , inp , adi , thd , grid % y , v )
          else if ( bc (l) == inflowbogey ) then
             dirichlet = .true.
             call bc_inflowbogey ( l , adi , thd , grid % y , v )
          else if ( bc (l) == inflowjet ) then
             dirichlet = .true.
             call bc_inflowjet ( l , adi , thd , grid % y , grid % z , v , time )
          else if ( bc (l) == wallinjection ) then
             dirichlet = .true.
             call bc_wallinjection ( l , inp , adi , thd , grid % x , grid % z , time , v )
          else if ( bc (l) == wallperturbation ) then
             dirichlet = .true.
             call bc_wallperturbation ( l , inp , adi , thd , grid % x , v , time )
          else if ( bc (l) == premixfix ) then
             dirichlet = .true.
             call bc_premix ( l , v )
          else if ( bc (l) == nscbc_out ) then
             dirichlet = .true.
             call bc_nscbc_out ( l , adi , thd , dtime , grid , T , W_i , cp , v )
          else if ( bc (l) == supersonicflow ) then
             dirichlet = .true.
             call bc_supersonicflow ( l , inp , adi , thd , v )

          else
             call abort_mpi ('boundary condition ' // trim (bc(l)) // ' not defined')
          end if

          if (dirichlet) call upd_prim_var_dirichlet ( thd , face_domain (l) , T , W_i , cp , ha , v )

       end if

    end do


  end subroutine upd_boundaries


!> \brief Extrapolation pseudo-BC (to initialize solutions).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine upd_extrapolation (v)


    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v !< conserved variables array


    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -10 , 10 , -20 , 20 , -30 , 30 /)

    integer (ip) :: l


    do l = 1 , ndim+ndim

       if ( neigh (l) == MPI_PROC_NULL .or. bc (l) == periodic ) then
          call bc_extrapolation ( l , v )
       end if

    end do


  end subroutine upd_extrapolation


!> \brief Non-reflecting BC (only for WENO).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine NR_boundaries ( thd , grid , T , W_i , cp , ha , v , fl )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    type (inp_grid) , intent (in)                                  :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: fl   !< inviscid flux array


    integer (ip) :: l


    do l = 1 , ndim + ndim

       if ( bc (l) == noreflection ) then
          call bc_noreflectionweno ( l , thd , grid , T , W_i , cp , ha , v , fl )
       end if

    end do


  end subroutine NR_boundaries


!> \brief Conserved variables pseudo-BC.
!!
!! Updates the array of conserved variables in the inner domain.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine upd_prim_var_domain ( thd , T , W_i , cp , ha , v )


    type (thd_type) , intent (in)                                  :: thd !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T   !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp  !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha  !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v   !< conserved variables array


    if (forcesum) call force_sumY ( 0 , v )

    call prim_inv_var ( 0 , thd , v , W_i , T , cp , ha )


  end subroutine upd_prim_var_domain


!> \brief Ghost BC.
!!
!! Updates the array of conserved variables in the ghost points.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine upd_prim_var_ghost ( thd , T , W_i , cp , ha , v )


    type (thd_type) , intent (in)                                  :: thd !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T   !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp  !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha  !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v   !< conserved variables array

    integer (ip) , dimension (ndimmax+ndimmax) , parameter :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip) :: l


    if (forcesum) then

       do l = 1 , ndim+ndim
          call force_sumY ( face_domain (l) , v )
          call prim_inv_var ( face_domain (l) , thd , v , W_i , T , cp , ha )
       end do

    else

       do l = 1 , ndim+ndim
          call prim_inv_var ( face_domain (l) , thd , v , W_i , T , cp , ha )
       end do

    end if


  end subroutine upd_prim_var_ghost


!> \brief Dirichlet BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine upd_prim_var_dirichlet ( thd , face , T , W_i , cp , ha , v )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    integer (ip) , intent (in)                                     :: face !< face domain
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    if (forcesum) call force_sumY ( face , v )

    call prim_inv_var ( face , thd , v , W_i , T , cp , ha )


  end subroutine upd_prim_var_dirichlet


!> \brief Extrapolation BC.
!!
!! Symmetry of all conservative variables (on \f$ ng \f$ ghost points) on both sides of the boundary
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_extrapolation ( face , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip) :: i , j , k , l


    if ( face == W ) then ! ghost point at x beginning

       do l = 1 , nv
          do k = sz , ez
             do j = sy , ey
                do i = 1 , ng
                   v (sx-i,j,k,l) = v (sx,j,k,l)
                end do
             end do
          end do
       end do

    else if ( face == E ) then ! ghost point at x ending

       do l = 1 , nv
          do k = sz , ez
             do j = sy , ey
                do i = 1 , ng
                   v (ex+i,j,k,l) = v (ex,j,k,l)
                end do
             end do
          end do
       end do

    else if ( face == S ) then ! ghost point at y beginning

       do l = 1 , nv
          do k = sz , ez
             do j = 1 , ng
                do i = sx , ex
                   v (i,sy-j,k,l) = v (i,sy,k,l)
                end do
             end do
          end do
       end do

    else if ( face == N ) then ! ghost point at y ending

       do l = 1 , nv
          do k = sz , ez
             do j = 1 , ng
                do i = sx , ex
                   v (i,ey+j,k,l) = v (i,ey,k,l)
                end do
             end do
          end do
       end do

    else if ( face == B ) then ! ghost point at z beginning

       do l = 1 , nv
          do k = 1 , ng
             do j = sy , ey
                do i = sx , ex
                   v (i,j,sz-k,l) = v (i,j,sz,l)
                end do
             end do
          end do
       end do

    else if ( face == F ) then ! ghost point at z ending

       do l = 1 , nv
          do k = 1 , ng
             do j = sy , ey
                do i = sx , ex
                   v (i,j,ez+k,l) = v (i,j,ez,l)
                end do
             end do
          end do
       end do

    end if


  end subroutine bc_extrapolation


!> \brief Non-reflecting BC.
!!
!! Same subroutine than "Extrapolation" but with one garbage coefficient to artificially 
!! hoover the flow.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_noreflection ( face , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip) :: i , j , k , l
    real (dp) , parameter :: garbage = 10.0_dp


    if ( face == W ) then

       do l = 1 , nv
          do k = sz , ez
             do j = sy , ey
                do i = 1 , ng
                   v (sx-i,j,k,l) = v (sx,j,k,l) * garbage
                end do
             end do
          end do
       end do

    else if ( face == E ) then

       do l = 1 , nv
          do k = sz , ez
             do j = sy , ey
                do i = 1 , ng
                   v (ex+i,j,k,l) = v (ex,j,k,l) * garbage
                end do
             end do
          end do
       end do

    else if ( face == S ) then

       do l = 1 , nv
          do k = sz , ez
             do j = 1 , ng
                do i = sx , ex
                   v (i,sy-j,k,l) = v (i,sy,k,l) * garbage
                end do
             end do
          end do
       end do

    else if ( face == N ) then

       do l = 1 , nv
          do k = sz , ez
             do j = 1 , ng
                do i = sx , ex
                   v (i,ey+j,k,l) = v (i,ey,k,l) * garbage
                end do
             end do
          end do
       end do

    else if ( face == B ) then

       do l = 1 , nv
          do k = 1 , ng
             do j = sy , ey
                do i = sx , ex
                   v (i,j,sz-k,l) = v (i,j,sz,l) * garbage
                end do
             end do
          end do
       end do

    else if ( face == F ) then

       do l = 1 , nv
          do k = 1 , ng
             do j = sy , ey
                do i = sx , ex
                   v (i,j,ez+k,l) = v (i,j,ez,l) * garbage
                end do
             end do
          end do
       end do

    end if


  end subroutine bc_noreflection


!> \brief Symmetry plane BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_symmetryplane ( face , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip) :: i , j , k , l


    if ( face == W ) then


       do l = 1 , nv
          do k = sz , ez
             do j = sy , ey
                do i = 1 , ng
                   v (sx-i,j,k,l) = v (sx+i,j,k,l)
                end do
             end do
          end do
       end do

       do k = sz , ez
          do j = sy , ey
             do i = sx , sx
                v (i,j,k,2) = 0.0_dp
             end do
          end do
       end do


    else if ( face == E ) then


       do l = 1 , nv
          do k = sz , ez
             do j = sy , ey
                do i = 1 , ng
                   v (ex+i,j,k,l) = v (ex-i,j,k,l)
                end do
             end do
          end do
       end do

       do k = sz , ez
          do j = sy , ey
             do i = ex , ex
                v (i,j,k,2) = 0.0_dp
             end do
          end do
       end do


    else if ( face == S ) then


       do l = 1 , nv
          do k = sz , ez
             do j = 1 , ng
                do i = sx , ex
                   v (i,sy-j,k,l) = v (i,sy+j,k,l)
                end do
             end do
          end do
       end do

       do k = sz , ez
          do j = sy , sy
             do i = sx , ex
                v (i,j,k,3) = 0.0_dp
             end do
          end do
       end do


    else if ( face == N ) then


       do l = 1 , nv
          do k = sz , ez
             do j = 1 , ng
                do i = sx , ex
                   v (i,ey+j,k,l) = v (i,ey-j,k,l)
                end do
             end do
          end do
       end do

       do k = sz , ez
          do j = ey , ey
             do i = sx , ex
                v (i,j,k,3) = 0.0_dp
             end do
          end do
       end do


    else if ( face == B ) then


       do l = 1 , nv
          do k = 1 , ng
             do j = sy , ey
                do i = sx , ex
                   v (i,j,sz-k,l) = v (i,j,sz+k,l)
                end do
             end do
          end do
       end do

       do k = sz , sz
          do j = sy , ey
             do i = sx , ex
                v (i,j,k,4) = 0.0_dp
             end do
          end do
       end do


    else if ( face == F ) then


       do l = 1 , nv
          do j = sy , ey
             do i = sx , ex
                do k = 1 , ng
                   v (i,j,ez+k,l) = v (i,j,ez-k,l)
                end do
             end do
          end do
       end do

       do k = ez , ez
          do j = sy , ey
             do i = sx , ex
                v (i,j,k,4) = 0.0_dp
             end do
          end do
       end do


    end if


  end subroutine bc_symmetryplane


!> \brief Adiabatic wall BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_adiabaticwall ( face , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip) :: i , j , k , l


    if ( face == W ) then


       do l = 1 , nv
          do k = sz , ez
             do j = sy , ey
                do i = 1 , ng
                   v (sx-i,j,k,l) = v (sx+i,j,k,l)
                end do
             end do
          end do
       end do

       do l = 2 , 4
          do k = sz , ez
             do j = sy , ey
                do i = sx , sx
                   v (i,j,k,l) = 0.0_dp
                end do
             end do
          end do
       end do


    else if ( face == E ) then


       do l = 1 , nv
          do k = sz , ez
             do j = sy , ey
                do i = 1 , ng
                   v (ex+i,j,k,l) = v (ex-i,j,k,l)
                end do
             end do
          end do
       end do

       do l = 2 , 4
          do k = sz , ez
             do j = sy , ey
                do i = ex , ex
                   v (i,j,k,l) = 0.0_dp
                end do
             end do
          end do
       end do


    else if ( face == S ) then


       do l = 1 , nv
          do k = sz , ez
             do j = 1 , ng
                do i = sx , ex
                   v (i,sy-j,k,l) = v (i,sy+j,k,l)
                end do
             end do
          end do
       end do

       do l = 2 , 4
          do k = sz , ez
             do j = sy , sy
                do i = sx , ex
                   v (i,j,k,l) = 0.0_dp
                end do
             end do
          end do
       end do


    else if ( face == N ) then


       do l = 1 , nv
          do k = sz , ez
             do j = 1 , ng
                do i = sx , ex
                   v (i,ey+j,k,l) = v (i,ey-j,k,l)
                end do
             end do
          end do
       end do

       do l = 2 , 4
          do k = sz , ez
             do j = ey , ey
                do i = sx , ex
                   v (i,j,k,l) = 0.0_dp
                end do
             end do
          end do
       end do


    else if ( face == B ) then


       do l = 1 , nv
          do k = 1 , ng
             do j = sy , ey
                do i = sx , ex
                   v (i,j,sz-k,l) = v (i,j,sz+k,l)
                end do
             end do
          end do
       end do

       do l = 2 , 4
          do k = sz , sz
             do j = sy , ey
                do i = sx , ex
                   v (i,j,k,l) = 0.0_dp
                end do
             end do
          end do
       end do


    else if ( face == F ) then


       do l = 1 , nv
          do k = 1 , ng
             do j = sy , ey
                do i = sx , ex
                   v (i,j,ez+k,l) = v (i,j,ez-k,l)
                end do
             end do
          end do
       end do

       do l = 2 , 4
          do k = ez , ez
             do j = sy , ey
                do i = sx , ex
                   v (i,j,k,l) = 0.0_dp
                end do
             end do
          end do
       end do


    end if


  end subroutine bc_adiabaticwall


!> \brief "Synthetic turbulence generators" BC.
!!
!!  Generate turbulence, synthetically, at the inlet of a computational domain.
!!  Available methods:
!!   - the white noise (isotropic, anispotropic)
!!   - the Smirnov method, based on Fourier techniques
!!   - the Klein method, based on Digital filters
!!
!! References:
!! -# J.M. Vedovoto, A.S. Neto, L.F.F. Silva, A. Mura, 
!!    Influence of synthetic inlet turbulence on the prediction of low Mach number flows,
!!    Computers and Fluids, 106, 135-153 (2015)
!! -# A. Smirnov, S. Shi, I. Celik, 
!!    Random flow generation technique for Large Eddy Simulation and particle-dynamics modeling,
!!    J Fluids Eng 2001, 123:359-71
!! -# M. Klein, A. Sadiki, and J. Janicka. 
!!    A Digital Filter Based Generation of Inflow Data for Spatially Developing Direct Numerical or Large Eddy Simulations,
!!    J of Computational Physics, 186, 2, 2003
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_syntheticturb ( face , inp , adi , thd , y , z , v )


    integer (ip)                                  , intent (in)    :: face !< face domain
    type (inp_type)                               , intent (in)    :: inp  !< external input data derived type
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    type (thd_type)                               , intent (in)    :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:)       , intent (in)    :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip)                               :: ok , i , j , k , l
    integer (ip)                               :: ipr , jpr , kpr
    real (dp)                                  :: P , T , rho , ux , vy , wz
    real (dp)                                  :: Ws_i , hm , Ma
    real (dp)                                  :: Ya (nrv+npv+nvv) , has (nrv)

    real (dp) , parameter                      :: twothirds = 2.0_dp / 3.0_dp
    real (dp)                                  :: fluct , fluct_max_i , & ! fluctuations
                                                  factor , vin
    real (dp)                                  :: epsv (ndim) , noise (ndim)

    real (dp) , allocatable , dimension (:,:,:) :: randu

    character (len_default)                    :: stg_type


    !======== Parameters of the problem ========
!       X  =  X/Xtot   *   Xtot   / adi % X_ref ! where X/Xtot is function of gamma and Mach number
        ! Values of P and T for a
        ! Mach number = 2 and T0 = 1695 K, P0 = 0.41e6 Pa and gamma = 1.272 (vitiated air)
        Ma = 2.0_dp
        P  = 0.056e6_dp / adi % P_ref
        T  = 1108.0_dp  / adi % T_ref
        vin = - 13.0_dp / adi % u_ref

        vy = 0.0_dp
        wz = 0.0_dp

        Ya(:) = 0.0_dp
        ! Species in chem.inp
        ! air only
!        Ya (inp % index_O2) = 0.233_dp ! O2 (0.233_dp for air)
!        Ya (inp % index_N2) = 0.767_dp ! N2 (0.767_dp for air)

        ! vitiated air for H2/air combustion
        Ya (inp % index_O2)  = 0.2527_dp ! O2
        Ya (inp % index_H2O) = 0.1631_dp ! H2O
        Ya (inp % index_N2)  = 0.5842_dp ! N2

        stg_type = 'White_Noise_isotropic'
!        stg_type = 'White_Noise_anisotropic'
!        stg_type = 'Klein_Digital_Filter'
        factor = 10.0_dp
        noise (1) = 20.0_dp ! Noise amplitude in percent
        noise (2) = 10.0_dp
        fluct_max_i = 80.0_dp / adi % u_ref
    !===========================================

    noise = noise * 1.0e-2_dp
    fluct_max_i = 1.0_dp / fluct_max_i

    call Wmix_i_scalar ( thd , Ya , Ws_i )

    if ( face == W ) then

       if ( stg_type == 'White_Noise_isotropic' ) then
       !================================================================

          do k = sz , ez
             do j = sy , ey

                fluct = sqrt ( twothirds * inp % ext % tke (j) )

                T = inp % ext % T (j)

                rho = P / ( T * Ws_i )

                call ha_scalar ( thd , T , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * Ya (l)
                end do

!                do l = 1 , ndim
!                   call random_number ( epsv (l) )
!                   epsv (l) = epsv (l) + epsv (l) - 1.0_dp
!                end do

                epsv (1) = uperturb ( inp , y(j) , z(k) )
                epsv (1) = noise (1) * ( epsv (1) + epsv (1) - 1.0_dp )
                call random_number ( epsv (2) )
                epsv (2) = noise (2) * ( epsv (2) + epsv (2) - 1.0_dp )

!                ux = inp % ext % uin (j) + fluct * epsv (1) * factor
!                vy = vin * ( 1.0_dp + epsv (2) ) !+ fluct * epsv (2) !* factor
!                wz = 0.0_dp              !+ fluct * epsv (3) !* factor

                ux = inp % ext % uin (j) * ( 1.0_dp + fluct * fluct_max_i * epsv (1) )
                vy = 0.0_dp + ux * epsv (2)
                wz = 0.0_dp


                v (sx-ng:sx,j,k,1) = rho
                v (sx-ng:sx,j,k,2) = rho * ux
                v (sx-ng:sx,j,k,3) = rho * vy
                v (sx-ng:sx,j,k,4) = rho * wz
                v (sx-ng:sx,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
                do l = 1 , nrv+npv+nvv
                   v (sx-ng:sx,j,k,niv+l) = rho * Ya (l)
                end do

             end do
          end do

       else if ( stg_type == 'White_Noise_anisotropic' ) then
       !================================================================

          do k = sz , ez
             do j = sy , ey

                T = inp % ext % T (j)

                rho = P / ( T * Ws_i )

                call ha_scalar ( thd , T , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * Ya (l)
                end do

                do i = sx-ng , sx

                   do l = 1 , ndim
                      call random_number ( epsv (l) )
                      epsv (l) = epsv (l) + epsv (l) - 1.0_dp
                   end do

                   ux = inp % ext % uin (j) + ( inp % ext % A11 (j) * epsv (1) ) * factor
                   vy = 0.0_dp              + ( inp % ext % A21 (j) * epsv (1) + &
                                                inp % ext % A22 (j) * epsv (2) ) * factor
                   wz = 0.0_dp              + ( inp % ext % A31 (j) * epsv (1) + &
                                                inp % ext % A32 (j) * epsv (2) + &
                                                inp % ext % A33 (j) * epsv (3) ) * factor

                   v (i,j,k,1) = rho
                   v (i,j,k,2) = rho * ux
                   v (i,j,k,3) = rho * vy
                   v (i,j,k,4) = rho * wz
                   v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
                   do l = 1 , nrv+npv+nvv
                      v (i,j,k,niv+l) = rho * Ya (l)
                   end do

                end do

             end do
          end do

       else if ( stg_type == 'Klein_Digital_Filter' ) then
       !================================================================

          allocate( randu ( sy:ey , sz:ez , ndim ) , &
                    stat=ok )
          if ( ok > 0 ) call abort_mpi ('error allocate bc_syntheticturb randu')

          randu = 0.0_dp ! Initialization

          do k = sz , ez
             do j = sy , ey

                ! Filtering with convolution (between i and (j,k) loop !)
                do kpr = -nzfw , nzfw
                   do jpr = -nyfw , nyfw
                      do ipr = -nxfw , nxfw
                         do l = 1 , ndim
                            randu (j,k,l) = randu (j,k,l)                  &
                                          + inp % ext % bijk (ipr,jpr,kpr) &
                                          * inp % ext % rand (ipr,j+jpr,k+kpr,l)
                         end do
                      end do
                   end do
                end do

                T = inp % ext % T (j)

                rho = P / ( T * Ws_i )

                call ha_scalar ( thd , T , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * Ya (l)
                end do

                ux = inp % ext % uin (j) + ( inp % ext % A11 (j) * randu (j,k,1) ) * factor
                vy = 0.0_dp              + ( inp % ext % A21 (j) * randu (j,k,1) + &
                                             inp % ext % A22 (j) * randu (j,k,2) ) * factor
                wz = 0.0_dp              + ( inp % ext % A31 (j) * randu (j,k,1) + &
                                             inp % ext % A32 (j) * randu (j,k,2) + &
                                             inp % ext % A33 (j) * randu (j,k,3) ) * factor

                v (sx-ng:sx,j,k,1) = rho
                v (sx-ng:sx,j,k,2) = rho * ux
                v (sx-ng:sx,j,k,3) = rho * vy
                v (sx-ng:sx,j,k,4) = rho * wz
                v (sx-ng:sx,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
                do l = 1 , nrv+npv+nvv
                   v (sx-ng:sx,j,k,niv+l) = rho * Ya (l)
                end do

             end do
          end do

          deallocate( randu )

       else

          call abort_mpi ('bc_syntheticturb ''' // trim (stg_type) // ''' is not defined')

       end if

    else

       write (*,*) 'face ' , face , ' has not been declared in bc_syntheticturb'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if


  end subroutine bc_syntheticturb


!> \brief Shift and fill the random field for Digital filter BC.
!!
!! After each iteration shift the hole random field of the digital filter
!! from rand (-nxfw:nxfw-1,:,:) <- rand (-nxfw+1:nxfw,:,:) and fill the 
!! last (first ?!) plane rand (nxfw,:,:) with a new field of random number.
!! Correspond to the step (g) of the Klein et al. algorithme in their article.
!!
!! References :
!! -# M. Klein, A. Sadiki and J. Janicka
!!    "A digital filter based generation of inflow data for spatially developing
!!     direct numerical or large eddy simulation", 
!!    Journal of computational physics, 2003
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_shift_fill_rfield ( inp )


    type (inp_type) , intent (inout)   :: inp  !< input derived type

    logical          :: loop
    integer (ip)     :: j , k , l
    real (dp)        :: epsv


    ! Test if the neigh (i) == MPI_PROC_NULL and if BC for STG is needed
    ! This to avoid allocating useless arrays (inner process)
    !===================================================================
    loop = .false.
    call apply_stg ( loop )
    if ( .not. loop ) return ! leave the subroutine

    inp % ext % rand ( -nxfw : nxfw-1 ,:,:,:) = &
          inp % ext % rand ( -nxfw+1 : nxfw ,:,:,:)

    do l = 1 , ndim
       do k = -nzfw+sz , nzfw+ez
          do j = -nyfw+sy , nyfw+ey
             call random_number (epsv)
             inp % ext % rand (nxfw,j,k,l) = epsv + epsv - 1.0_dp
          end do
       end do
    end do


  end subroutine bc_shift_fill_rfield


!> \brief Shock BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_shock ( face , inp , adi , thd , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    type (inp_type) , intent (in)                                  :: inp  !< input derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip)                                                   :: i , j , k , l
    real (dp)                                                      :: v1 , T1 , P1 , rho1 , h1 , W_i
    real (dp)                                                      :: beta , theta
    real (dp)                                                      :: v2 , T2 , P2 , rho2 , h2 , ret2
    real (dp)                                                      :: ru2 , rv2 , rw2
    real (dp) , dimension (nrv)                                    :: Y1 , has


    if ( face == S ) then


       beta = 33.0_dp * acos (-1.0_dp) / 180.0_dp

       v1 = 1634.00_dp  / adi % u_ref
       T1 = 1475.0_dp   / adi % T_ref
       P1 = 94232.25_dp / adi % P_ref

       Y1 (1:nrv) = inp % Y0o (1:nrv)
       call Wmix_i_scalar ( thd , Y1 , W_i )

       call RH_variables ( thd , Y1 , P1 , T1 , rho1 , v1 , h1 , beta , &
                           P2 , T2 , rho2 , v2 , h2 ,theta )

       rho2 = P2 / ( T2 * W_i )
       ru2  = rho2 * v2 * cos (theta)
       rv2  = rho2 * v2 * sin (theta)
       rw2  = rho2 * v (sx,sy,sz,4) / v (sx,sy,sz,1)

       call ha_scalar ( thd , T2 , has )
       h2 = 0.0_dp
       do l = 1 , nrv
          h2 = h2 + has (l) * Y1 (l)
       end do

       ret2 = rho2 * h2 - P2 + 0.5_dp * ( ru2*ru2 + rv2*rv2 + rw2*rw2 ) / rho2

       Y1 (1:nrv) = rho2 * Y1 (1:nrv)

       do k = sz , ez
          do j = sy-ng , sy
             do i = sx , ex
                v (i,j,k,1) = rho2
                v (i,j,k,2) = ru2
                v (i,j,k,3) = rv2
                v (i,j,k,4) = rw2
                v (i,j,k,5) = ret2
                do l = niv+1 , nv
                   v (i,j,k,l) = Y1 (l-niv)
                end do
             end do
          end do
       end do


    else

       write (*,*) 'face ' , face , ' has not been declared in bc_shock'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if


  end subroutine bc_shock


!> \brief Post-shock + Slip wall-BC.
!!
!! Boundary condition to use with the DMR (Double Mach Reflection)
!! initialization.  It separates the bottom domain boundary (\f$ y
!! \f$-boundary) in 2 parts, the first boundary part is applied the
!! post-shock conditions between \f$ [x_{\min}, x_0] \f$ with \f$ x_0
!! \f$ is the abscisse of the shock impact and the second boundary
!! part is applied the slipwall conditions between \f$ ]x_0, x_{\max}]
!! \f$.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_postshock_slipwall ( face , adi , thd , x , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip)    :: i , j , k , l

    real (dp)       :: x0 ! abscisse of shock impact
    real (dp)       :: thetaw, theta ! angle of the angled wedges & angle(wall, v2)
    real (dp)       :: Ms, Us ! shock Mach & velocity
    real (dp)       :: v1 , T1 , P1 , rho1 , hm1 , cp1, gamma1
    real (dp)       :: v2 , T2 , P2 , rho2 , hm2
    real (dp)       :: Ya (nrv+npv+nvv) , Ws_i

    !====== Parameters of the DMR problem ======
        ! Species : O2   N2 (in chem.inp)
        Ya (1)  = 0.00_dp ! 02 (0.23_dp for air)
        Ya (2)  = 1.00_dp ! N2 (0.77_dp for air)

        ! Angled wedge
        x0      = 0.0_dp / adi % L_ref
        thetaw  = 36.0_dp * pi/180.0_dp

        Ms = 4.5_dp

        ! Pre-shock conditions
        P1 = 101325.0_dp / adi % P_ref
        T1 = 300.0_dp    / adi % T_ref
        v1 = 0.0_dp      / adi % u_ref ! velocity magnitude
    !===========================================

    call Wmix_i_scalar ( thd , Ya , Ws_i )

    ! Calculate of pre-shock conditions (RIGHT)
        rho1 = P1 / ( T1 * Ws_i )

        call cp_scalar ( thd , T1 , Ya , cp1 )
        gamma1 = thd % gam2 * Ws_i / cp1
        gamma1 = 1.0_dp / ( 1.0_dp - gamma1 )

    Us = Ms * sqrt ( gamma1 * p1 / rho1 )
         !state '1' because of shock velocity is measured from the speed of sound upstream of shock

    ! Calculate of post-shock conditions (LEFT)
        ! Using RH subroutine, by assuming the shock is STEADY, in shock point of view
            v1 = Us ! in shock point of view (STEADY shock)

            call RH_variables ( thd, Ya, p1, T1, rho1, v1, hm1, pi/2.0_dp, &
                                p2, T2, rho2, v2, hm2, theta )
            ! where vnorme2 is now the velocity relative to the reference related to the STEADY shock
            v2 = Us - v2 ! conversion to the reference related to the laboratory (UNSTEADY shock)
            v1 = 0.0_dp / adi % u_ref ! put back the rest state upstream of the shock (UNSTEADY shock)
        ! END of using RH subroutine

    ! End of calculation of pre- and post-shock conditions

    if ( face == S ) then

        do k = sz , ez ! z-direction
            do i = sx , ex ! x-direction

                if ( x(i) <= x0 ) then
                    ! Post-shock conditions on [Xmin, X0]

                    do j = sy-ng , sy ! y-direction (ghost points)
                        v (i,j,k,1) =   rho2
                        v (i,j,k,2) =   rho2 * v2 * cos(thetaw) ! velocity relative to cartesian coord.
                        v (i,j,k,3) = - rho2 * v2 * sin(thetaw) ! < 0 because flow direction downwards
                        v (i,j,k,4) =   0.0_dp
                        v (i,j,k,5) =   rho2 * hm2 - P2 + 0.5_dp * rho2 * v2*v2
                        do l = niv+1 , nv
                           v (i,j,k,l) = rho2 * Ya (l-niv)
                        end do
                    end do

                else
                    ! Slipwall conditions on ]X0, Xmax]

                    do j = 1 , ng ! y-direction (ghost points & its symmetry)
                        v (i,sy-j,k,1) =   v (i,sy+j,k,1) ! rho    is even
                        v (i,sy-j,k,2) =   v (i,sy+j,k,2) ! rho*u  is even
                        v (i,sy-j,k,3) = - v (i,sy+j,k,3) ! rho*v  is odd
                        v (i,sy-j,k,4) =   v (i,sy+j,k,4) ! rho*w  is even
                        v (i,sy-j,k,5) =   v (i,sy+j,k,5) ! rho*et is even
                        do l = niv+1 , nv
                            v (i,sy-j,k,l) = v (i,sy+j,k,l) ! rho*chi is even
                        end do
                    end do


                end if

            end do
        end do

    else

       write (*,*) 'face ' , face , ' has not been declared in bc_postshock_slipwall'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if


  end subroutine bc_postshock_slipwall


!> \brief "Cheng inflow" BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_inflowcheng ( face , adi , thd , y , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip)                               :: i , j , k , l
    real (dp)                                  :: u1_ , u2_ , cs1 , cs2 , uc , t1_ , t2_ , p1_ , p2_
    real (dp)                                  :: Deltau , Sigmau , d2O0_i , dw0_i
    real (dp)                                  :: P , T , rho , Ws_i , hm , y0 , ux , vy , wz
    real (dp)                                  :: Ya (nrv+npv+nvv) , has (nrv)


    if ( face == W ) then


       u1_ = 2462.74_dp / adi % u_ref ; cs1 = 2238.85_dp / adi % u_ref
       u2_ = 1213.96_dp / adi % u_ref ; cs2 = 838.11_dp  / adi % u_ref

       t1_ = 875.0_dp  / adi % T_ref
       t2_ = 1750.0_dp / adi % T_ref

       p1_ = 109.5e3_dp / adi % P_ref
       p2_ = 109.5e3_dp / adi % P_ref

       Deltau = u1_ - u2_
       Sigmau = u1_ + u2_

       uc = ( cs1 * u1_ + cs2 * u2_ ) / ( cs1 + cs2 )

       dw0_i  = 1.839e-4_dp / adi % L_ref
       d2O0_i = 0.25_dp * dw0_i

       ! invert
       d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
       dw0_i  = 1.0_dp / dw0_i


       y0 = 0.0_dp

       ux = 0.0_dp
       vy = 0.0_dp
       wz = 0.0_dp


       do k = sz , ez
          do j = sy , ey
             do i = sx-ng , sx


                ux = 0.5_dp * ( Sigmau    + Deltau    * tanh ( y(j) * d2O0_i ) )
                T  = 0.5_dp * ( (t1_+t2_) + (t1_-t2_) * tanh ( y(j) * d2O0_i ) )
                P  = 0.5_dp * ( (p1_+p2_) + (p1_-p2_) * tanh ( y(j) * d2O0_i ) )

                Ya (:)  = 0.0_dp
!                Ya (2)  = 0.5_dp * ( (0.071_dp) + (0.071_dp) * tanh ( y(j) * d2O0_i ) )
                Ya (2)  = 0.5_dp * ( (1.000_dp) + (1.000_dp) * tanh ( y(j) * d2O0_i ) )
                Ya (4)  = 0.5_dp * ( (0.245_dp) - (0.245_dp) * tanh ( y(j) * d2O0_i ) )
                Ya (6)  = 0.5_dp * ( (0.175_dp) - (0.175_dp) * tanh ( y(j) * d2O0_i ) )
                Ya (9)  = 0.5_dp * ( (0.580_dp) - (0.580_dp) * tanh ( y(j) * d2O0_i ) )
!                Ya (9)  = 0.5_dp * ( (1.509_dp) + (0.349_dp) * tanh ( y(j) * d2O0_i ) )
!                call X_to_Y ( thd , Ya )

                call Wmix_i_scalar ( thd , Ya , Ws_i )
                rho = P / ( T * Ws_i )

                call ha_scalar ( thd , T , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * Ya (l)
                end do


                v (i,j,k,1) = rho
                v (i,j,k,2) = rho * ux
                v (i,j,k,3) = rho * vy
                v (i,j,k,4) = rho * wz
                v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
                do l = 1 , nrv+npv+nvv
                   v (i,j,k,niv+l) = rho * Ya (l)
                end do


             end do
          end do
       end do


    else

       write (*,*) 'face ' , face , ' has not been declared in bc_inflow'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if


  end subroutine bc_inflowcheng


!> \brief "Miller inflow" BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_inflowmiller ( face , adi , thd , y , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip)                               :: i , j , k , l
    real (dp)                                  :: u1_ , u2_ , cs1 , cs2 , uc , t1_ , t2_ , p1_ , p2_
    real (dp)                                  :: Deltau , Sigmau , d2O0_i , dw0_i
    real (dp)                                  :: P , T , rho , Ws_i , hm , y0 , ux , vy , wz
    real (dp)                                  :: Ya (nrv+npv+nvv) , has (nrv)


    if ( face == W ) then


       u1_ = 168.00_dp / adi % u_ref ; cs1 = 351.32_dp / adi % u_ref
       u2_ = 953.11_dp / adi % u_ref ; cs2 = 770.27_dp / adi % u_ref

       t1_ = 269.0_dp  / adi % T_ref
       t2_ = 1475.0_dp / adi % T_ref

       p1_ = 94232.25_dp / adi % P_ref
       p2_ = 94232.25_dp / adi % P_ref

       Deltau = u1_ - u2_
       Sigmau = u1_ + u2_

       uc = ( cs1 * u1_ + cs2 * u2_ ) / ( cs1 + cs2 )

       dw0_i  = 4.71e-5_dp / adi % L_ref
       d2O0_i = 0.25_dp * dw0_i

       ! invert
       d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
       dw0_i  = 1.0_dp / dw0_i


       y0 = 0.0_dp

       ux = 0.0_dp
       vy = 0.0_dp
       wz = 0.0_dp


       do k = sz , ez
          do j = sy , ey
             do i = sx-ng , sx


                ux = 0.5_dp * ( Sigmau    + Deltau    * tanh ( y(j) * d2O0_i ) )
                T  = 0.5_dp * ( (t1_+t2_) + (t1_-t2_) * tanh ( y(j) * d2O0_i ) )
                P  = 0.5_dp * ( (p1_+p2_) + (p1_-p2_) * tanh ( y(j) * d2O0_i ) )

                Ya (:)  = 0.0_dp
                Ya (1)  = 0.5_dp * ( (0.100_dp) + (0.100_dp) * tanh ( y(j) * d2O0_i ) )
                Ya (2)  = 0.5_dp * ( (0.230_dp) - (0.230_dp) * tanh ( y(j) * d2O0_i ) )
                Ya (3)  = 0.5_dp * ( (0.250_dp) - (0.250_dp) * tanh ( y(j) * d2O0_i ) )
                Ya (4)  = 0.5_dp * ( (1.420_dp) + (0.380_dp) * tanh ( y(j) * d2O0_i ) )
                call X_to_Y ( thd , Ya )

                call Wmix_i_scalar ( thd , Ya , Ws_i )
                rho = P / ( T * Ws_i )

                call ha_scalar ( thd , T , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * Ya (l)
                end do


                v (i,j,k,1) = rho
                v (i,j,k,2) = rho * ux
                v (i,j,k,3) = rho * vy
                v (i,j,k,4) = rho * wz
                v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
                do l = 1 , nrv+npv+nvv
                   v (i,j,k,niv+l) = rho * Ya (l)
                end do


             end do
          end do
       end do


    else

       write (*,*) 'face ' , face , ' has not been declared in bc_inflow'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if


  end subroutine bc_inflowmiller


!> \brief "Cheng-Miller inflow" BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_inflowmixchengmiller ( face , inp , adi , thd , y , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    type (inp_type) , intent (in)                                  :: inp  !< input derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip)                               :: i , j , k , l
    real (dp)                                  :: u1_ , u2_ , t1_ , t2_ , p1_ , p2_
    real (dp)                                  :: Deltau , Sigmau , d2O0_i , dw0_i
    real (dp)                                  :: P , T , rho , Ws_i , hm , y0 , ux , vy , wz
    real (dp)                                  :: Ya (nrv+npv+nvv) , has (nrv)


    if ( face == W ) then


       u1_ = 669.14_dp / adi % u_ref
       u2_ = 2185.35_dp / adi % u_ref

       t1_ = 545.0_dp  / adi % T_ref
       t2_ = 1475.0_dp / adi % T_ref

       p1_ = 94232.25_dp / adi % P_ref
       p2_ = 94232.25_dp / adi % P_ref

       Deltau = u1_ - u2_
       Sigmau = u1_ + u2_

       dw0_i  = 6.28e-5_dp / adi % L_ref
       d2O0_i = 0.25_dp * dw0_i

       ! invert
       d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )
       dw0_i  = 1.0_dp / dw0_i


       y0 = 0.0_dp

       ux = 0.0_dp
       vy = 0.0_dp
       wz = 0.0_dp


       do k = sz , ez
          do j = sy , ey
             do i = sx-ng , sx


                ux = 0.5_dp * ( Sigmau    + Deltau    * tanh ( y(j) * d2O0_i ) )
                T  = 0.5_dp * ( (t1_+t2_) + (t1_-t2_) * tanh ( y(j) * d2O0_i ) )
                P  = 0.5_dp * ( (p1_+p2_) + (p1_-p2_) * tanh ( y(j) * d2O0_i ) )

                do l = 1,nrv-1
                   Ya (l) = 0.5_dp * ( inp % Y1f(l)+inp % Y0o(l) +                           &
                                     (inp % Y1f(l)-inp % Y0o(l)) * tanh ( y(j) * d2O0_i ) )
                end do
                Ya (nrv) = 1.0_dp - sum ( Ya (1:nrv-1) )
                if ( npv > 0 ) then
                   do l = nrv+1,nrv+npv
                      Ya (l) = 0.5_dp * ( inp % Y1f(l)+inp % Y0o(l) +                        &
                                        (inp % Y1f(l)-inp % Y0o(l)) * tanh ( y(j) * d2O0_i ) )
                   end do
                end if
!                call X_to_Y ( thd , Ya )

                call Wmix_i_scalar ( thd , Ya , Ws_i )
                rho = P / ( T * Ws_i )

                call ha_scalar ( thd , T , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * Ya (l)
                end do


                v (i,j,k,1) = rho
                v (i,j,k,2) = rho * ux
                v (i,j,k,3) = rho * vy
                v (i,j,k,4) = rho * wz
                v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
                do l = 1 , nrv+npv+nvv
                   v (i,j,k,niv+l) = rho * Ya (l)
                end do


             end do
          end do
       end do


    else

       write (*,*) 'face ' , face , ' has not been declared in bc_inflow'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if


  end subroutine bc_inflowmixchengmiller


!> \brief "Bogey inflow" BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_inflowbogey ( face , adi , thd , y , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip)                               :: i , j , k , l
    real (dp)                                  :: u1_ , u2_ , cs1 , cs2 , uc , p1_ , p2_ , r1_ , r2_
    real (dp)                                  :: Deltau , Sigmau , d2O0_i , dw0_i
    real (dp)                                  :: P , T , rho , Ws_i , hm , y0 , ux , vy , wz
    real (dp)                                  :: Ya (nrv+npv+nvv) , has (nrv)


    if ( face == W ) then


       u1_ = 100.0_dp / adi % u_ref ; cs1 = 338.767_dp / adi % u_ref
       u2_ = 50.0_dp  / adi % u_ref ; cs2 = 338.767_dp / adi % u_ref

       p1_ = 1.0e5_dp / adi % P_ref
       p2_ = 1.0e5_dp / adi % P_ref

       r1_ = 1.22_dp / adi % rho_ref
       r2_ = 1.22_dp / adi % rho_ref

       Deltau = u1_ - u2_
       Sigmau = u1_ + u2_

       uc = ( cs1 * u1_ + cs2 * u2_ ) / ( cs1 + cs2 )

       dw0_i  = 1.6e-3_dp / adi % L_ref
       d2O0_i = 0.25_dp * dw0_i

       ! invert
       dw0_i  = 1.0_dp / dw0_i
       d2O0_i = 1.0_dp / ( d2O0_i + d2O0_i )

       y0 = 0.0_dp

       ux = 0.0_dp
       vy = 0.0_dp
       wz = 0.0_dp


       do k = sz , ez
          do j = sy , ey
             do i = sx-ng , sx


                ux  = 0.5_dp * ( Sigmau    + Deltau    * tanh ( y(j) * d2O0_i ) )
                rho = 0.5_dp * ( (r1_+r2_) + (r1_-r2_) * tanh ( y(j) * d2O0_i ) )
                P   = 0.5_dp * ( (p1_+p2_) + (p1_-p2_) * tanh ( y(j) * d2O0_i ) )

                Ya (1) = 0.21_dp
                Ya (2) = 0.79_dp
                call X_to_Y ( thd , Ya )

                call Wmix_i_scalar ( thd , Ya , Ws_i )
                T = P / ( rho * Ws_i )

                call ha_scalar ( thd , T , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * Ya (l)
                end do


                v (i,j,k,1) = rho
                v (i,j,k,2) = rho * ux
                v (i,j,k,3) = rho * vy
                v (i,j,k,4) = rho * wz
                v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
                do l = 1 , nrv+npv+nvv
                   v (i,j,k,niv+l) = rho * Ya (l)
                end do


             end do
          end do
       end do


    else

       write (*,*) 'face ' , face , ' has not been declared in bc_inflow'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if


  end subroutine bc_inflowbogey


!> \brief Jet inflow BC.(modifié par Romain)
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_inflowjet ( face , adi , thd , y , z , v , time )


    integer (ip)                                  , intent (in)    :: face !< face domain
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    type (thd_type)                               , intent (in)    :: thd  !< thermodynamic derived type
    real (dp)                                     , intent (in)    :: time !< time
    real (dp) , allocatable , dimension (:)       , intent (in)    :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip)                    :: i , j , k , l
    real (dp) , parameter           :: NPR = 15.0_dp
    real (dp)                       :: hbuse , tmaxy
    real (dp)                       :: u_1 , u_2 , v_1 , v_2 , &
                                       T_1 , T_2 , P_1 , P_2 , &
                                       Ya_1, Ya_2
    real (dp)                       :: delt , y0 , rrr , cs , gamma , cps
    real (dp)                       :: P , T , rho , hm , Ws_i , ux , vy , wz
    real (dp) , dimension (nrv)     :: has
    real (dp) , dimension (nrv+npv+nvv) :: Ya

    hbuse = 0.0010_dp / adi % L_ref
    tmaxy = 0.0001_dp / adi % time_ref

    if ( face == W ) then

       v_1 =  0.0_dp / adi % u_ref

       u_2 = 20.0_dp / adi % u_ref
       v_2 =  0.0_dp / adi % u_ref

       T_1 = 1000.0_dp  / adi % t_ref
       T_2 =  300.0_dp  / adi % t_ref

       P_1 = 101325.0_dp * NPR / adi % p_ref
       P_2 = 101325.0_dp       / adi % p_ref

       Ya_1 = 1.0_dp ! Y targetted in injection
       Ya_2 = 0.0_dp ! Y at initial time

       delt = 1.0_dp / (0.04_dp * hbuse)
       y0   = hbuse  / 2.0_dp

       ux = 0.0_dp
       vy = 0.0_dp
       wz = 0.0_dp

       if ( time > tmaxy) then

       else if ( time <= tmaxy ) then

          P_1  = ( P_2  - ( P_2  - P_1  ) * ( time / tmaxy) * ( time / tmaxy) * ( time / tmaxy))
          Ya_1 = ( Ya_2 - ( Ya_2 - Ya_1 ) * ( time / tmaxy) * ( time / tmaxy) * ( time / tmaxy))
          T_1  = ( T_2  - ( T_2  - T_1  ) * ( time / tmaxy) * ( time / tmaxy) * ( time / tmaxy))

       else if ( time == 0 ) then

          P_1  = P_2
          Ya_1 = Ya_2
          T_1  = T_2

       end if

       if      (thd % Nspc == 2) then
          Ya (1) = 0.233_dp ; Ya (2) = 0.767_dp
       else if (thd % Nspc == 3) then
          Ya (1) = Ya_1 ; Ya (2) = 0.0_dp ; Ya (3) = 1.0_dp - Ya_1
       else if (thd % Nspc == 9) then
          Ya (2) = Ya_1 ; Ya (4) = 0.0_dp ; Ya (9) = 1.0_dp - Ya_1
       else
          write (*,*) 'module_BCs : Problem with species initialisation'
       end if

       call Wmix_i_scalar ( thd , Ya , Ws_i )
       call cp_scalar     ( thd , T_1 , Ya , cps )

       gamma = thd % gam2 * Ws_i / cps
       gamma = 1.0_dp / ( 1.0_dp - gamma )
       cs    = sqrt ( gamma *  Ws_i * T_1 )
       u_1   = 1.01_dp * cs

       do k = sz , ez
          do j = sy , ey
             do i = sx-ng , sx

                rrr = sqrt ( y(j) * y(j) + z(k) * z(k) ) ! 3D axi
!                rrr = sqrt ( y(j) * y(j) )               ! 3D plan

                ux = 0.5_dp * ( (u_1 + u_2) - (u_1 - u_2) * tanh ( ( rrr - y0 ) * delt ) )
                T  = 0.5_dp * ( (T_1 + T_2) - (T_1 - T_2) * tanh ( ( rrr - y0 ) * delt ) )
                P  = 0.5_dp * ( (P_1 + P_2) - (P_1 - P_2) * tanh ( ( rrr - y0 ) * delt ) )

                if      (thd % Nspc == 2) then ! Boundary composition for air/air + xi calculation
                   Ya(:)= 0.0_dp
                   Ya(1)= 0.233_dp
                   Ya(2)= 1.0_dp - Ya (1)
                   Ya(3)= 0.5_dp * ( (1.0_dp + 0.0_dp) - (1.0_dp - 0.0_dp) * tanh ( ( rrr - y0 ) * delt ) )
                else if (thd % Nspc == 3) then ! Boundary composition for H2/air + xi calculation
                   Ya(:)= 0.0_dp
                   Ya(1)= 0.5_dp * ( (Ya_1          + 0.000_dp) - (Ya_1          - 0.000_dp) * tanh ( ( rrr - y0 ) * delt ) )
                   Ya(2)= 0.5_dp * ( (0.0_dp        + 0.233_dp) - (0.0_dp        - 0.233_dp) * tanh ( ( rrr - y0 ) * delt ) )
                   Ya(3)= 0.5_dp * ( (1.0_dp - Ya_1 + 0.767_dp) - (1.0_dp - Ya_1 - 0.767_dp) * tanh ( ( rrr - y0 ) * delt ) )
                   Ya(4)= 0.5_dp * ( (1.0_dp        + 0.000_dp) - (1.0_dp        - 0.000_dp) * tanh ( ( rrr - y0 ) * delt ) )
                else if (thd % Nspc == 9) then ! Boundary composition for H2/air reactive + xi calculation
                   Ya(:) = 0.0_dp
                   Ya(2) = 0.5_dp * ( (Ya_1          + 0.000_dp) - (Ya_1          - 0.000_dp) * tanh ( ( rrr - y0 ) * delt ) )
                   Ya(4) = 0.5_dp * ( (0.0_dp        + 0.233_dp) - (0.0_dp        - 0.233_dp) * tanh ( ( rrr - y0 ) * delt ) )
                   Ya(9) = 0.5_dp * ( (1.0_dp - Ya_1 + 0.767_dp) - (1.0_dp - Ya_1 - 0.767_dp) * tanh ( ( rrr - y0 ) * delt ) )
                   Ya(10)= 0.5_dp * ( (1.0_dp        + 0.000_dp) - (1.0_dp        - 0.000_dp) * tanh ( ( rrr - y0 ) * delt ) )
                else
                  write (*,*) 'module_BCs : Problem with species initialisation'
                end if


                call Wmix_i_scalar ( thd , Ya , Ws_i )
                rho = P / ( T * Ws_i )

                call ha_scalar ( thd , T , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * Ya (l)
                end do

                v (i,j,k,1) = rho
                v (i,j,k,2) = rho * ux
                v (i,j,k,3) = rho * vy
                v (i,j,k,4) = rho * wz
                v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
                do l = 1 , nrv+npv+nvv
                   v (i,j,k,niv+l) = rho * Ya (l)
                end do

             end do
          end do
       end do

    else

       write (*,*) 'face ' , face , ' has not been declared in bcjet'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if


  end subroutine bc_inflowjet


!> \brief "Supersonic inflow" BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_supersonicflow ( face , inp , adi , thd , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    type (inp_type) , intent (in)                                  :: inp  !< input derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip)                               :: j , k , l
    real (dp)                                  :: P , T , rho , ux , vy , wz
    real (dp)                                  :: Ws_i , hm , Ma
    real (dp)                                  :: Ya (nrv+npv+nvv) , has (nrv)

    !======== Parameters of the problem ========
!       X  =  X/Xtot   *   Xtot   / adi % X_ref ! where X/Xtot is function of gamma and Mach number
        ! Values of P and T for a
        ! Mach number = 2 and T0 = 1695 K, P0 = 0.41e6 Pa and gamma = 1.336 (air)
        Ma = 2.0_dp
!        P  = 0.1295322834_dp * 0.41e6_dp / adi % P_ref
!        T  = 0.5981507118_dp * 1695_dp   / adi % T_ref
        ! Mach number = 2 and T0 = 1695 K, P0 = 0.41e6 Pa and gamma = 1.272 (vitiated air)
        P  = 0.056e6_dp / adi % P_ref
        T  = 1108.0_dp  / adi % T_ref

        vy = 0.0_dp
        wz = 0.0_dp

        Ya(:) = 0.0_dp
        ! Species in chem.inp
        ! air only
!        Ya (inp % index_O2) = 0.233_dp ! O2 (0.233_dp for air)
!        Ya (inp % index_N2) = 0.767_dp ! N2 (0.767_dp for air)

        ! vitiated air for H2/air combustion
        Ya (inp % index_O2)  = 0.2527_dp ! O2
        Ya (inp % index_H2O) = 0.1631_dp ! H2O
        Ya (inp % index_N2)  = 0.5842_dp ! N2


        ! Coefficient for the tanh function for velocity profil, 
        ! to have 2mm of boundary layer thickness. It's proportional to the BLT.
!        coeff_delta99 = 0.042693309_dp
    !===========================================

!    delty = 1.0_dp / (coeff_delta99 * Ly)
!    deltz = 1.0_dp / (coeff_delta99 * Lz)

    call Wmix_i_scalar ( thd , Ya , Ws_i )

!    rho = P / ( T * Ws_i )

!    call cp_scalar ( thd , T , Ya , cp )
!    gamma0 = thd % gam2 * Ws_i / cp
!    gamma0 = 1.0_dp / ( 1.0_dp - gamma0 )
!    ux_inf = Ma * sqrt ( gamma0 * P / rho )

!    call ha_scalar ( thd , T , has )
!    hm = 0.0_dp
!    do l = 1 , nrv
!       hm = hm + has (l) * Ya (l)
!    end do

    if ( face == W ) then

       do k = sz , ez
          do j = sy , ey

!             ux = sqrt (     ux_inf * tanh ( y(j) * delty )                    & ! (-y,+y) = (adiabaticwall,symmetry)
!                         * - ux_inf * tanh ( (abs(z(k))-0.5_dp*Lz) * deltz ) )   ! (-z,+z) = adiabaticwall

!             ux = ux_inf * tanh ( y(j) * delty )   ! (-y,+y) = (adiabaticwall,symmetry)
!                                                   ! (-z,+z) = extrapolation

!             ldim = y(j) * adi % L_ref
!             if ( 0.0_dp <= ldim .and. ldim < 2.024E-04 ) then
!                ux = -4.565E+17 * ldim**4 + 3.255E+14 * ldim**3 - 8.711E+10 * ldim**2 + 1.167E+07 * ldim
!             else if ( ldim < 9.657E-04 ) then
!                ux = -2.240E+08 * ldim**2 + 6.902E+05 * ldim + 6.191E+02
!             else if ( ldim < 3.104E-03 ) then
!                ux = 3.980E+10 * ldim**3 - 3.310E+08 * ldim**2 + 9.166E+05 * ldim + 4.521E+02
!             else
!                ux = 7.117E+02 * ldim + 1.297E+03
!             end if
!             ux = ux / adi % u_ref

             T = inp % ext % T (j)
             ux = inp % ext % uin (j)

             rho = P / ( T * Ws_i )

             call ha_scalar ( thd , T , has )
             hm = 0.0_dp
             do l = 1 , nrv
                hm = hm + has (l) * Ya (l)
             end do

             v (sx-ng:sx,j,k,1) = rho
             v (sx-ng:sx,j,k,2) = rho * ux
             v (sx-ng:sx,j,k,3) = rho * vy
             v (sx-ng:sx,j,k,4) = rho * wz
             v (sx-ng:sx,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
             do l = 1 , nrv+npv+nvv
                v (sx-ng:sx,j,k,niv+l) = rho * Ya (l)
             end do

          end do
       end do

    else

       write (*,*) 'face ' , face , ' has not been declared in bc_supersonicflow'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if


  end subroutine bc_supersonicflow


!> \brief Wall injection BC.
!!
!! Jet injection across a (Slip) wall in +/- y-boundary ( face == N or S )
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_wallinjection ( face , inp , adi , thd , x , z , time , v )


    integer (ip)                                  , intent (in)    :: face !< face domain
    type (inp_type)                               , intent (in)    :: inp  !< input derived type
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    type (thd_type)                               , intent (in)    :: thd  !< thermodynamic derived type
    real (dp)                                     , intent (in)    :: time !< time
    real (dp) , allocatable , dimension (:)       , intent (in)    :: x    !< y-coordinate array
    real (dp) , allocatable , dimension (:)       , intent (in)    :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip)                    :: i , j , k , l , &
                                       lH2 , lxi , lsgsvar ! H2, passive scalar, SGS scalar variance index
    integer (ip)                    :: sj , ej , wj ! index of N or S side (start, end and wall index)
    real (dp)                       :: signe ! injection direction
    real (dp)                       :: delt , rinj , dinj, tmaxy_i , x0 , z0 , beta
    integer (ip)                    :: ntmaxy
    real (dp)                       :: vinj , Tinj , Pinj , Yainj , Minj ! steady values
    real (dp)                       :: vext , Text , Pext ! values for the injection with tanh function
    real (dp)                       :: rrr , wrk , pwr
    real (dp)                       :: cp0 , gamma0
    real (dp)                       :: rho , hm , Ws_i , ux , vy , wz
    real (dp) , dimension (nrv)     :: has
    real (dp) , dimension (nrv+npv+nvv) :: Ya

    if      ( face == N ) then
       signe = - 1.0_dp
       sj = ey
       ej = ey+ng
       wj = ey
    else if ( face == S ) then
       signe =   1.0_dp
       sj = sy-ng
       ej = sy
       wj = sy
    else
       write (*,*) 'face ' , face , ' has not been declared in bc_wallinjection'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if


    !======== Parameters of the problem ========
        ! (Slip)wall / co-flow
        vext = 0.0_dp / adi % u_ref

!        Pext = 101325.0_dp / adi % p_ref ! Ambient
!        Text = 300.0_dp    / adi % t_ref

        ! Values of Pext and Text for a
        ! Mach number = 2 and T0 = 1695 K, P0 = 0.41e6 Pa and gamma = 1.336 (air)
!        Pext  = 0.1295322834_dp * 0.41e6_dp / adi % P_ref
!        Text  = 0.5981507118_dp * 1695_dp   / adi % T_ref
        ! Mach number = 2 and T0 = 1695 K, P0 = 0.41e6 Pa and gamma = 1.272 (vitiated air)
        Pext  = 0.056e6_dp   / adi % P_ref
        Text  = 1203.0_dp    / adi % T_ref
        ! Stop conditions
!        Pext  = 347335.410774691_dp / adi % P_ref
!        Text  = 1743.0953656247_dp  / adi % T_ref

        ! Injector geometry
        x0   = 0.2000_dp / adi % L_ref ! if it modified, changed also in module_viscflux.f90 !!!
        z0   = 0.0000_dp / adi % L_ref ! if it modified, changed also in module_viscflux.f90 !!!
        dinj = 0.0020_dp / adi % L_ref ! if it modified, changed also in module_viscflux.f90 !!!

        ! Injection conditions
        Minj = 1.0_dp
        Pinj = 502918.273101923_dp / adi % p_ref ! NPR = 18 in Pamb = 0.13*0.41e6 Pa (H2, gamma = 1.425)
        Tinj = 248.3767836263_dp   / adi % t_ref
        Yainj = 1.0_dp

        ! Perturbation properties
!        xa        = 0.122_dp / adi % L_ref
!        xb        = 0.122_dp / adi % L_ref + dinj
!        Aperturb  = 0.04_dp
!        uinf      = 1300.0_dp / adi % u_ref
!        betafreq  = 75000_dp * adi % time_ref ! in rad/s
!        Pptb      = 0.0586e6_dp / adi % P_ref
!        Tptb      = 1203_dp   / adi % T_ref

!        Pinj = 16e5_dp             / adi % p_ref ! NPR = 30 (air)
!        Tinj = 1000_dp             / adi % t_ref
        tmaxy_i = 2.0e-7_dp  / adi % time_ref ! time injection etablishement
        ntmaxy = 5

        Ya (:) = 0.0_dp
        if      (thd % Nspc == 2) then ! air only
           Ya (inp % index_O2) = 0.233_dp ! O2
           Ya (inp % index_N2) = 0.767_dp ! N2
        else
           lH2 = inp % index_H2
           Ya (lH2) = 1.0_dp  ! H2
        end if

        if ( npv == 1 ) then
           lxi = thd % Nspc + 1
        else if ( npv > 1 ) then
           call abort_mpi ('ERROR: BCs wallinjection npv > 1')
        end if

        if ( nvv == 1 ) then
           lsgsvar = nrv + npv + nvv
        else if ( nvv > 1 ) then
           call abort_mpi ('ERROR: BCs wallinjection nvv > 1')
        end if


        pwr = 0.13_dp ! for r/rinj injection function, these coef depends on the pressure
    !===========================================

    tmaxy_i = float(ntmaxy) / tmaxy_i
    rinj = 0.5_dp * dinj
    delt = 1.0_dp / (0.035_dp * dinj) ! proportional to the gradient of the tanh/erf function
    beta = 1.0_dp                     ! coefficient which indicate the attachment point with the wall

    call Wmix_i_scalar ( thd , Ya , Ws_i )
    call cp_scalar ( thd , Tinj , Ya , cp0 )

    gamma0 = thd % gam2 * Ws_i / cp0
    gamma0 = 1.0_dp / ( 1.0_dp - gamma0 )
    vinj = signe * Minj * sqrt ( gamma0 * Ws_i * Tinj )

!    vinj = 1435.0_dp / adi % u_ref ! to conserve Qm = 1.85 g/s

    ! Initialization
    ux = 0.0_dp
    vy = 0.0_dp
    wz = 0.0_dp

    ! gradual increase of Tinj, Pinj and Yainj
    wrk   = time * tmaxy_i
    if ( wrk < float(ntmaxy) ) then
       wrk   = 1.0_dp - exp ( - wrk )
       Tinj  = ( Tinj - Text ) * wrk + Text
       Pinj  = ( Pinj - Pext ) * wrk + Pext
       vinj  =   vinj          * wrk
    end if


    do k = sz , ez
       do i = sx , ex

          rrr = sqrt ( (x(i)-x0)*(x(i)-x0) + (z(k)-z0)*(z(k)-z0) )

          if ( rrr > rinj ) then

!             if ( x(i) < xa .or. xb < x(i) ) then 

                ! wall face S/N  v1.6.1
                 do l = 1 , nv
                    do j = 1 , ng
                       v (i,wj-int(signe)*j,k,l) = v (i,wj+int(signe)*j,k,l)
                    end do
                 end do
                 v (i,wj,k,2:4) = 0.0_dp

!             else if ( xa <= x(i) .and. x(i) <= xb ) then

!                ! vitiated air for H2/air combustion
!                Ya (:) = 0.0_dp
!                Ya (inp % index_O2)  = 0.2527_dp ! O2
!                Ya (inp % index_H2O) = 0.1631_dp ! H2O
!                Ya (inp % index_N2)  = 0.5842_dp ! N2

!                ! perturbation
!                vy = vperturb ( inp , xa , xb , Aperturb , uinf , betafreq , x(i) , z(k) , time )
!                vy = abs (vy)
!!                vy = Aperturb * uinf

!                call Wmix_i_scalar ( thd , Ya , Ws_i )
!                rho = Pptb / ( Tptb * Ws_i )

!                call ha_scalar ( thd , Tptb , has )
!                hm = 0.0_dp
!                do l = 1 , nrv
!                   hm = hm + has (l) * Ya (l)
!                end do

!                do j = sj , ej

!                   v (i,j,k,1) = rho
!                   v (i,j,k,2) = rho * ux
!                   v (i,j,k,3) = rho * vy
!                   v (i,j,k,4) = rho * wz
!                   v (i,j,k,5) = rho * hm - Pptb + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
!                   do l = 1 , nrv+npv+nvv
!                      v (i,j,k,niv+l) = rho * Ya (l)
!                   end do

!                end do ! y-direction

!             end if

          else 

             ! uniform injection
             Ya (:) = 0.0_dp
             Ya (lH2) = Yainj
             if ( npv == 1 ) Ya (lxi) = Yainj        ! \tilde{xi}

!             Ya (6) = 0.0_dp            ! xi_variance
!             Ya (7) = Yainj * Yainj     ! \tilde{xi^2}
!             Ya (8) = Yainj * Yainj     ! \tilde{xi}^2
!             Ya (9) = 0.0_dp            ! departure from maximal variance

             ! tanh/erf injection
!             Ya (lH2)  = 0.5_dp * ( (Yainj   + 0.0_dp   ) - (Yainj   - 0.0_dp   ) * erf ( ( rrr - rinj * beta ) * delt ) )
!             Ya (2)    = 0.5_dp * ( (0.0_dp  + 0.2527_dp) - (0.0_dp  - 0.2527_dp) * erf ( ( rrr - rinj * beta ) * delt ) ) ! O2
!             Ya (3)    = 0.5_dp * ( (0.0_dp  + 0.1631_dp) - (0.0_dp  - 0.1631_dp) * erf ( ( rrr - rinj * beta ) * delt ) ) ! H2O
!             Ya (4)    = 0.5_dp * ( (0.0_dp  + 0.5842_dp) - (0.0_dp  - 0.5842_dp) * erf ( ( rrr - rinj * beta ) * delt ) ) ! N2
!             Ya (lxi)  = 0.5_dp * ( (Yainj   + 0.0_dp   ) - (Yainj   - 0.0_dp   ) * erf ( ( rrr - rinj * beta ) * delt ) )

!             Ya (lH2)  = 0.0_dp    - (Yainj   - 0.0_dp   ) * erf ( ( rrr - rinj * beta ) * delt )
!             Ya (2)    = 0.2527_dp - (0.0_dp  - 0.2527_dp) * erf ( ( rrr - rinj * beta ) * delt ) ! O2
!             Ya (3)    = 0.1631_dp - (0.0_dp  - 0.1631_dp) * erf ( ( rrr - rinj * beta ) * delt ) ! H2O
!             Ya (4)    = 0.5842_dp - (0.0_dp  - 0.5842_dp) * erf ( ( rrr - rinj * beta ) * delt ) ! N2
!             Ya (lxi)  = 0.0_dp    - (Yainj   - 0.0_dp   ) * erf ( ( rrr - rinj * beta ) * delt )

!             if ( nvv > 0 ) then
!                Ya (lsgsvar) = Ya (lxi) * ( 1.0_dp - Ya (lxi) )
!                Ya (7)      = Ya (lxi) * Ya (lxi) + Ya (6)
!                Ya (8)      = Ya (lxi) * Ya (lxi)
!                Ya (9)      = Ya (lxi) * ( 1.0_dp - Ya (lxi) ) - Ya (6)
!             end if

             vy = 0.0_dp - (vinj   - 0.0_dp ) * erf ( ( rrr - rinj * beta ) * delt )


             call Wmix_i_scalar ( thd , Ya , Ws_i )
             rho = Pinj / ( Tinj * Ws_i )

             call ha_scalar ( thd , Tinj , has )
             hm = 0.0_dp
             do l = 1 , nrv
                hm = hm + has (l) * Ya (l)
             end do

             v (i,sj:ej,k,1) = rho
             v (i,sj:ej,k,2) = rho * ux
             v (i,sj:ej,k,3) = rho * vy
             v (i,sj:ej,k,4) = rho * wz
             v (i,sj:ej,k,5) = rho * hm - Pinj + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
             do l = 1 , nrv+npv+nvv
                v (i,sj:ej,k,niv+l) = rho * Ya (l)
             end do


          end if !( rrr > rinj )

       end do ! x-direction
    end do ! z-direction


  end subroutine bc_wallinjection


!> \brief Wall perturbation BC.
!!
!! Jet injection across a (Slip) wall in +/- y-boundary ( face == N or S )
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_wallperturbation ( face , inp , adi , thd , x , v , time )


    integer (ip)                                  , intent (in)    :: face !< face domain
    type (inp_type)                               , intent (in)    :: inp  !< input derived type
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    type (thd_type)                               , intent (in)    :: thd  !< thermodynamic derived type
    real (dp)                                     , intent (in)    :: time !< time
    real (dp) , allocatable , dimension (:)       , intent (in)    :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip)                    :: i , j , k , l , &
                                       lxi ! passive scalar index
    integer (ip)                    :: sj , ej , wj ! index of N or S side (start, end and wall index)
    real (dp)                       :: signe ! injection direction
    real (dp)                       :: tstop , x0 , epsx , epsv , noise
    real (dp)                       :: vinj , Tinj , Pinj
    real (dp)                       :: P , T , rho , hm , Ws_i , ux , vy , wz
    real (dp) , dimension (nrv)     :: has
    real (dp) , dimension (nrv+npv+nvv) :: Ya


    if      ( face == N ) then
       signe = - 1.0_dp
       sj = ey
       ej = ey+ng
       wj = ey
    else if ( face == S ) then
       signe =   1.0_dp
       sj = sy-ng
       ej = sy
       wj = sy
    else
       write (*,*) 'face ' , face , ' has not been declared in bc_wallperturbation'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if


    !======== Parameters of the problem ========
        epsx = 1.0e-5_dp

        ! jet perturbation features
        x0   = 0.25_dp * Lx
        Pinj = 1.0e5_dp    / adi % p_ref
        Tinj = 300_dp      / adi % t_ref
        vinj = 1.0_dp      / adi % u_ref
        noise = 1.0_dp     ! Noise amplitude in percent
        tstop = 1.0e-4_dp  / adi % time_ref ! duration of the perturbation

        Ya(:) = 0.0_dp
        if      (thd % Nspc == 2) then ! air only
            Ya (inp % index_O2) = 0.233_dp ! O2
            Ya (inp % index_N2) = 0.767_dp ! N2
        else
            write (*,*) 'bc_wallperturbation : Problem with species initialisation'
            call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        if ( npv == 1 ) then
            lxi = thd % Nspc + 1
        else
            write (*,*) 'bc_wallperturbation : wallperturbation can not support more than one passive scalar'
            call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if
    !===========================================


    ! Initialization
    ux = 0.0_dp
    vy = 0.0_dp
    wz = 0.0_dp


    if ( time > tstop ) then ! ordinary adiabaticwall


       do k = sz , ez
          do i = sx , ex

             ! wall face S  v1.6.1
              do l = 1 , nv
                 do j = 1 , ng
                    v (i,sy-j,k,l) = v (i,sy+j,k,l)
                 end do
              end do
              do l = 2 , 4
                 do j = sy , sy
                    v (i,j,k,l) = 0.0_dp
                 end do
              end do

          end do
       end do


    else ! perturbation wall


       do k = sz , ez
          do i = sx , ex

             if ( abs ( x0 - x(i) ) > epsx ) then ! wall

                  ! wall face S  v1.6.1
                   do l = 1 , nv
                      do j = 1 , ng
                         v (i,sy-j,k,l) = v (i,sy+j,k,l)
                      end do
                   end do
                   do l = 2 , 4
                      do j = sy , sy
                         v (i,j,k,l) = 0.0_dp
                      end do
                   end do

             else ! perturbation line

                vy = signe * vinj
                T  = Tinj
                P  = Pinj

                ! add white noise
                call random_number (epsv)
                epsv = noise * ( epsv + epsv - 1.0_dp )
                vy = vy * ( 1.0_dp + epsv ) 

                call Wmix_i_scalar ( thd , Ya , Ws_i )
                rho = P / ( T * Ws_i )

                call ha_scalar ( thd , T , has )
                hm = 0.0_dp
                do l = 1 , nrv
                   hm = hm + has (l) * Ya (l)
                end do

                do j = sj , ej

                   v (i,j,k,1) = rho
                   v (i,j,k,2) = rho * ux
                   v (i,j,k,3) = rho * vy
                   v (i,j,k,4) = rho * wz
                   v (i,j,k,5) = rho * hm - P + 0.5_dp * rho * ( ux*ux + vy*vy + wz*wz )
                   do l = 1 , nrv+npv+nvv
                      v (i,j,k,niv+l) = rho * Ya (l)
                   end do

                end do ! y-direction

             end if !( abs ( x0 - x(i) ) > epsx )

          end do ! x-direction
       end do ! z-direction

    end if ! (time > tstop)


  end subroutine bc_wallperturbation


!> \brief Calculate wall-normal velocity perturbation.
!!
!! Blowing and suction function to help transition from laminar to turbulent flow
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function vperturb ( inp , xa , xb , Aperturb , uinf , beta , x , z , t )


    type (inp_type) , intent (in) :: inp  !< input derived type
    real (dp) , intent (in) :: xa , xb    !< perturbation boundaries
    real (dp) , intent (in) :: Aperturb   !< amplification factor
    real (dp) , intent (in) :: uinf       !< longitudinal velocity magnitude
    real (dp) , intent (in) :: beta       !< fundamental temporal frequency of disturbance
    real (dp) , intent (in) :: x          !< x-coordinate array
    real (dp) , intent (in) :: z          !< z-coordinate array
    real (dp) , intent (in) :: t          !< time
    real (dp)               :: vperturb   !< wall-normal velocity perturbation

    real (dp)               :: theta , eps
    real (dp)               :: fx , gz , ht ! spatial and temporal elementary function
    integer (ip)            :: i

    theta = pi2 * ( x - xa ) / ( xb - xa )
    fx = 4.0_dp * sin(theta) * ( 1.0_dp - cos(theta) ) * sqrt27_i

    gz = 0.0_dp
    do i = 1 , inp % ptb % zlmax
       call random_number (eps)
       gz = gz + inp % ptb % zl (i) * sin ( pi2 * i * (z/Lz + eps) )
    end do

    ht = 0.0_dp
    do i = 1 , inp % ptb % mmax
       call random_number (eps)
       ht = ht + inp % ptb % tm (i) * sin ( i * (beta*t  + pi2*eps) )
    end do

    vperturb = Aperturb * uinf * fx * gz * ht


  end function vperturb


!> \brief Calculate streamwise velocity perturbation.
!!
!! Use a specific random value profil for the streamwise velocity in the
!! whole plan of the flow, here BC = W
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function uperturb ( inp , y , z )


    type (inp_type) , intent (in) :: inp  !< input derived type
    real (dp) , intent (in) :: y          !< y-coordinate array
    real (dp) , intent (in) :: z          !< z-coordinate array
    real (dp)               :: uperturb   !< wall-normal velocity perturbation

    real (dp)               :: eps
    real (dp)               :: fy , gz ! spatial elementary function
    integer (ip)            :: i

    fy = 0.0_dp
    do i = 1 , inp % ptb % ylmax
       call random_number (eps)
       fy = fy + inp % ptb % yl (i) * sin ( pi2 * i * (y/Ly + eps) )
    end do

    gz = 0.0_dp
    do i = 1 , inp % ptb % zlmax
       call random_number (eps)
       gz = gz + inp % ptb % zl (i) * sin ( pi2 * i * (z/Lz + eps) )
    end do

    uperturb = fy * gz


  end function uperturb


!> \brief PSR BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine bc_premix ( face , v )


    integer (ip) , intent (in)                                     :: face !< face domain
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip) :: i , j , k


    if ( face == W ) then


       do k = sz , ez
          do j = sy , ey
             do i = 0 , ng
                v (sx-i,j,k,1)  = 19.1479207567555
                v (sx-i,j,k,2)  = 7.898997438620545E-002
                v (sx-i,j,k,3)  = 0.0_dp
                v (sx-i,j,k,4)  = 0.0_dp
                v (sx-i,j,k,5)  = -41.6382366337754
                v (sx-i,j,k,6)  = 0.0_dp
                v (sx-i,j,k,7)  = 0.331462115522950
                v (sx-i,j,k,8)  = 0.0_dp
                v (sx-i,j,k,9)  = 4.38438573417748
                v (sx-i,j,k,10) = 0.0_dp
                v (sx-i,j,k,11) = 0.0_dp
                v (sx-i,j,k,12) = 0.0_dp
                v (sx-i,j,k,13) = 0.0_dp
                v (sx-i,j,k,14) = 14.4320729070551
             end do
          end do
       end do


    else if ( face == E ) then


       do k = sz , ez
          do j = sy , ey
             do i = 0 , ng
                v (ex+i,j,k,1)  = 3.44960033237647
                v (ex+i,j,k,2)  = 7.898997438622181E-002
                v (ex+i,j,k,3)  = 0.0_dp
                v (ex+i,j,k,4)  = 0.0_dp
                v (ex+i,j,k,5)  = -41.9798174104654
                v (ex+i,j,k,6)  = 5.488991459764549E-007
                v (ex+i,j,k,7)  = 1.673624482286092E-005
                v (ex+i,j,k,8)  = 1.070214895214190E-004
                v (ex+i,j,k,9)  = 0.314913367487909
                v (ex+i,j,k,10) = 2.606437663612932E-003
                v (ex+i,j,k,11) = 0.531820062372202
                v (ex+i,j,k,12) = 4.256906039223869E-006
                v (ex+i,j,k,13) = 2.907000739338948E-007
                v (ex+i,j,k,14) = 2.60013161061314
             end do
          end do
       end do


    end if


  end subroutine bc_premix


!> \brief N-S non reflecting BC.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
 subroutine bc_nscbc_out ( face , adi , thd , dt , grid , T , W_i , cp , v)


    integer (ip)                                  , intent (in)    :: face !< face domain
    type (adi_type)                               , intent (in)    :: adi  !< non-dimensional derived type
    type (thd_type)                               , intent (in)    :: thd  !< thermodynamic derived type
    real (dp)                                     , intent (in)    :: dt   !< time step
    type (inp_grid) , intent (in)                                  :: grid !< grid derived type
    real (dp) , allocatable , dimension (:,:,:)   , intent (in)    :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:)   , intent (in)    :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:)   , intent (in)    :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array

    integer (ip)                    :: i , j , k , l , m
    real (dp)                       :: rho , rho_i , ux , vy , wz , Pres , Temp
    real (dp)                       :: gam , hm , machm , mach
    real (dp)                       :: cs , pinf
    real (dp)                       :: has (nrv)
    real (dp)                       :: drho , dux , dvy , dwz , dPres , dtime
    real (dp) , dimension (nrv+npv+nvv) :: Ya , dYa
    real (dp) , dimension (nrv+npv+nvv) :: L5_Ya , d5_Ya
    real (dp)                       :: L1 , L2 , L3 , L4 , L5
    real (dp)                       :: d1 , d2 , d3 , d4 , d5

    real (dp) , dimension (-4:4) , parameter :: &

       ! fourth-order one-side finite difference (LEFT)
         alph4ol = (/  0.0e0_dp               ,   0.0e0_dp              ,   0.0e0_dp                 , &
                       0.0e0_dp               , -25.0e0_dp / 12.0e0_dp  ,  48.0e0_dp / 12.0e0_dp     , &
                       -36.0e0_dp / 12.0e0_dp ,  16.0e0_dp / 12.0e0_dp  ,  -3.0e0_dp / 12.0e0_dp  /) , &

       ! fourth-order one-side finite difference (RIGHT)
         alph4or = (/  -3.0e0_dp / 12.0e0_dp   ,  16.0e0_dp / 12.0e0_dp  , -36.0e0_dp / 12.0e0_dp  , &
                       48.0e0_dp / 12.0e0_dp   , -25.0e0_dp / 12.0e0_dp  ,   0.0e0_dp              , &
                       0.0e0_dp                ,   0.0e0_dp              ,   0.0e0_dp               /)


    real (dp) , parameter :: sigma1 = 0.28_dp ! Tester influence paramètre

    Pinf =       101325.0_dp / adi % p_ref

    dtime = dt

    if (face == E) then

       do k = sz , ez
          do j = sy , ey


             Machm=0
             rho   = v (ex,j,k,1)
             rho_i = 1.0_dp / rho
             ux    = v (ex,j,k,2) * rho_i
             vy    = v (ex,j,k,3) * rho_i
             wz    = v (ex,j,k,4) * rho_i
             Pres  = v (ex,j,k,1) * T (ex,j,k) * W_i (ex,j,k)
             gam   = thd % gam2 * W_i (ex,j,k) / cp (ex,j,k)
             gam   = 1.0_dp / ( 1.0_dp - gam )

             cs=sqrt ( gam * Pres * rho_i )
             Mach=sqrt(ux*ux+vy*vy+wz*wz)/cs

             Machm=max(Machm,Mach) ! use of the local value of the Mach number



             drho = 0.0_dp
             dux  = 0.0_dp
             dvy  = 0.0_dp
             dwz  = 0.0_dp
             dPres= 0.0_dp
             dYa  = 0.0_dp

             do m = -4 , 0

                rho   = v (ex+m,j,k,1)
                rho_i = 1.0_dp / rho
                ux    = v (ex+m,j,k,2) * rho_i
                vy    = v (ex+m,j,k,3) * rho_i
                wz    = v (ex+m,j,k,4) * rho_i
                Pres  = v (ex+m,j,k,1) * T (ex+m,j,k) * W_i (ex+m,j,k)
                do l = 1 , nrv+npv+nvv
                   Ya (l) = v (ex+m,j,k,niv+l) * rho_i
                end do

                drho = drho + alph4or(m) * rho
                dux  = dux  + alph4or(m) * ux
                dvy  = dvy  + alph4or(m) * vy
                dwz  = dwz  + alph4or(m) * wz
                dPres= dPres+ alph4or(m) * Pres
                do l = 1 , nrv+npv+nvv
                   dYa (l) = dYa (l) + alph4or(m) * Ya (l)
                end do

             end do

             drho = - drho * grid % dx_i (ex)
             dux  = - dux  * grid % dx_i (ex)
             dvy  = - dvy  * grid % dx_i (ex)
             dwz  = - dwz  * grid % dx_i (ex)
             dPres= - dPres* grid % dx_i (ex)
             do l = 1 , nrv+npv+nvv
                dYa (l) = - dYa (l) * grid % dx_i (ex)
             end do

             L1    = 0.0_dp
             L2    = 0.0_dp
             L3    = 0.0_dp
             L4    = 0.0_dp
             L5    = 0.0_dp
             L5_Ya = 0.0_dp

             gam   = thd % gam2 * W_i (ex,j,k) / cp (ex,j,k)
             gam   = 1.0_dp / ( 1.0_dp - gam )
             cs    = sqrt ( gam * Pres * rho_i )

             if ( ux - cs < 0.0_dp ) then
                L1 = sigma1 * (1 - Machm * Machm) * cs / Lx * (Pres - Pinf)
                !write (*,*) j , k , 'L1 = ' , L1
             else
                L1 = (ux - cs) * (dPres - rho * cs * dux)
             end if

             if ( ux + cs > 0.0_dp ) then
                L5 = (ux + cs) * (dPres + rho * cs * dux)
             else
                L5 = (ux + cs) * (2.0_dp * dPres - L1 / (ux - cs))
             end if

             if ( ux > 0.0_dp ) then
                L2 = ux * (cs * cs * drho - dPres)
                L3 = ux * dvy
                L4 = ux * dwz
                do l = 1 , nrv+npv+nvv
                   L5_Ya (l) = ux * dYa (l)
                end do

             else
                L2 = ux * cs * cs * drho - &
                     0.5_dp * ux * (L5 / (ux + cs) + L1 / (ux - cs))
             end if

             d2 = 0.5_dp * (L5+L1)
             d1 = (L2 + d2) / (cs * cs)
             d3 = 0.5_dp * (L5-L1) / (cs * rho)
             d4 = L3
             d5 = L4
             do l = 1 , nrv+npv+nvv
                d5_Ya (l) = L5_Ya (l)
             end do

             Temp = T (ex,j,k)

             call ha_scalar ( thd , Temp , has )
             hm = 0.0_dp
             do l = 1 , nrv
                hm = hm + has (l) * Ya (l)
             end do

             do i = ex , ex+ng
                v (i,j,k,1) = rho - d1 * dtime
                v (i,j,k,2) = (rho * ux) - (ux * d1 + rho * d3) * dtime
                v (i,j,k,3) = (rho * vy) - (vy * d1 + rho * d4) * dtime
                v (i,j,k,4) = (rho * wz) - (wz * d1 + rho * d5) * dtime
                v (i,j,k,5) = (rho * hm - Pres + 0.5_dp * rho * (ux*ux+vy*vy+wz*wz)) - &
                              (0.5_dp*(ux*ux+vy*vy+wz*wz)*d1+rho*(ux*d3+vy*d4+wz*d5) + &
                              d2/(gam-1.0_dp))*dtime
                do l = 1 , nrv+npv+nvv
                   v (i,j,k,niv+l) = (rho * Ya (l)) - (Ya(l) * d1 + rho * d5_Ya (l)) * dtime
                end do
             end do

          end do
       end do


    else

       write (*,*) 'face ' , face , ' has not been declared in bc_subso_outflow'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

    end if

  end subroutine bc_nscbc_out


end module BCs
