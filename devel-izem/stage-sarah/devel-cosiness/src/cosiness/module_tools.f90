!------------------------------------------------------------------------------
! MODULE: tools
!------------------------------------------------------------------------------
!> \brief General tools and utilities.
!!
!! This module contains all the tools and utilites for the solver and
!! post-treatment.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module tools

  use parameters
  use deriv
  use input
  use type_thd
  use adim
  use parallel , only : rank , mpicode


  implicit none


  real (dp) , parameter , private :: epsi     = 1.0e-10_dp     !< to avoid divisions by zero
  real (dp) , parameter , private :: eps30    = 1.0e-30_dp


contains


!> \brief Create the necessary repertories for COSINESS.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine createrep (inp)


    type (inp_type) , intent (in) :: inp !< input derived type
    logical                       :: existed


#ifdef __INTEL_COMPILER
#define _DIR_ directory
#else
#define _DIR_ file    
#endif

    inquire ( _DIR_ = trim (dir_restart) , exist = existed )
    if (.not.existed) call system ( 'mkdir ' // trim (dir_restart) )

    if ( sum ( sum ( inp % probe_coord , dim=1 ) ) > 0.0 ) then
       inquire ( _DIR_ = trim (dir_probes) , exist = existed )
       if (.not.existed) call system ( 'mkdir ' // trim (dir_probes) )
    end if

    if ( inp % nvolume > 0 ) then
       inquire ( _DIR_ = trim (dir_statvol) , exist = existed )
       if (.not.existed) call system ( 'mkdir ' // trim (dir_statvol) )
    end if

    if ( inp % nxystat > 0 ) then
       inquire ( _DIR_ = trim (dir_statXY) , exist = existed )
       if (.not.existed) call system ( 'mkdir ' // trim (dir_statXY) )
    end if

    if ( inp % nxzstat > 0 ) then
       inquire ( _DIR_ = trim (dir_statXZ) , exist = existed )
       if (.not.existed) call system ( 'mkdir ' // trim (dir_statXZ) )
    end if

    if ( inp % nyzstat > 0 ) then
       inquire ( _DIR_ = trim (dir_statYZ) , exist = existed )
       if (.not.existed) call system ( 'mkdir ' // trim (dir_statYZ) )
    end if


  end subroutine createrep


!> \brief Provide metrics for the derivatives.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine metrics ( xt , yt , zt , x , y , z , dx_i , dy_i , dz_i , delta )


    real (dp) , allocatable , dimension (:) , intent (in)         :: xt    !< absolute x-coordinate array (m/m)
    real (dp) , allocatable , dimension (:) , intent (in)         :: yt    !< absolute y-coordinate array (m/m)
    real (dp) , allocatable , dimension (:) , intent (in)         :: zt    !< absolute z-coordinate array (m/m)
    real (dp) , allocatable , dimension (:) , intent (inout)      :: x     !< x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (inout)      :: y     !< y-coordinate array
    real (dp) , allocatable , dimension (:) , intent (inout)      :: z     !< z-coordinate array
    real (dp) , allocatable , dimension (:) , intent (inout)      :: dx_i  !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (inout)      :: dy_i  !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (inout)      :: dz_i  !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)  :: delta !< filter-width


    integer (ip) , parameter                                   :: order = 2*ng+1
    integer (ip)                                               :: ok
    integer (ip)                                               :: l
    integer (ip)                                               :: i  , j  , k
    real (dp)                                                  :: ndim_i
    real (dp) , allocatable , dimension (:)                    :: ddummy_i


    ndim_i = 1.0_dp / float (ndim)

    allocate ( x (sx-ng:ex+ng) , dx_i (sx-ng:ex+ng)               , &
               y (sy-ng:ey+ng) , dy_i (sy-ng:ey+ng)               , &
               z (sz-ng:ez+ng) , dz_i (sz-ng:ez+ng)               , &
               delta (sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )   , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate metrics')

    x = 0.0_dp
    y = 0.0_dp
    z = 0.0_dp

    ! Transfer the good indices
    do i = sx-ng , ex+ng
       l = max ( 1 , min ( i , ntx ) )
       x (l) = xt (l)
    end do
    if ( neigh (W) == MPI_PROC_NULL .or. bc (W) == periodic ) then
       do l = 1,ng
          x (sx-l) = x (sx) + x (sx) - x (sx+l)
       end do
    end if
    if ( neigh (E) == MPI_PROC_NULL .or. bc (E) == periodic ) then
       do l = 1,ng
          x (ex+l) = x (ex) + x (ex) - x (ex-l)
       end do
    end if


    do j = sy-ng , ey+ng
       l = max ( 1 , min ( j , nty ) )
       y (l) = yt (l)
    end do
    if ( neigh (S) == MPI_PROC_NULL .or. bc (S) == periodic ) then
       do l = 1,ng
          y (sy-l) = y (sy) + y (sy) - y (sy+l)
       end do
    end if
    if ( neigh (N) == MPI_PROC_NULL .or. bc (N) == periodic ) then
       do l = 1,ng
          y (ey+l) = y (ey) + y (ey) - y (ey-l)
       end do
    end if


    do k = sz-ng , ez+ng
       l = max ( 1 , min ( k , ntz ) )
       z (l) = zt (l)
    end do
    if ( neigh (B) == MPI_PROC_NULL .or. bc (B) == periodic ) then
       do l = 1,ng
          z (sz-l) = z (sz) + z (sz) - z (sz+l)
       end do
    end if
    if ( neigh (F) == MPI_PROC_NULL .or. bc (F) == periodic ) then
       do l = 1,ng
          z (ez+l) = z (ez) + z (ez) - z (ez-l)
       end do
    end if


    ! Calculate the inverse of the spatial discretization
    ! 1st : Calculate dx, dy, dz
    dx_i (sx-ng:ex+ng) = 1.0_dp
    dy_i (sy-ng:ey+ng) = 1.0_dp
    dz_i (sz-ng:ez+ng) = 1.0_dp


    if ( ndim >= 1 ) then


       if ( ex-sx+1 < order ) call abort_mpi ('error: not enough points in x-direction')

       allocate ( ddummy_i (sx-ng:ex+ng) , stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate metrics 2')

       ddummy_i (sx-ng:ex+ng) = dx_i (sx-ng:ex+ng)
       call dscalar ( sx-ng , ex+ng , ddummy_i , x , dx_i )

       deallocate (ddummy_i)


    end if


    if ( ndim >= 2 ) then


       if ( ey-sy+1 < order ) call abort_mpi ('not enough points in y-direction')

       allocate ( ddummy_i (sy-ng:ey+ng) , stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate metrics 3')

       ddummy_i (sy-ng:ey+ng) = dy_i (sy-ng:ey+ng)
       call dscalar ( sy-ng , ey+ng , ddummy_i , y , dy_i )

       deallocate (ddummy_i)


    end if


    if ( ndim >= 3 ) then


       if ( ez-sz+1 < order ) call abort_mpi ('not enough points in z-direction')

       allocate ( ddummy_i (sz-ng:ez+ng) , stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate metrics 4')

       ddummy_i (sz-ng:ez+ng) = dz_i (sz-ng:ez+ng)
       call dscalar ( sz-ng , ez+ng , ddummy_i , z , dz_i )

       deallocate (ddummy_i)


    end if


    ! Calculate the filter width with dx, dy, dz
    do k = sz-ng , ez+ng
       do j = sy-ng , ey+ng
          do i = sx-ng , ex+ng
             ! For dim = 1 : dy_i = dz_i = 1.0_dp
             ! For dim = 2 : dz_i = 1.0_dp
             delta (i,j,k) = ( dx_i (i) * dy_i (j) * dz_i (k) ) ** ndim_i
          end do
       end do
    end do


    ! 2nd : Calculate 1/dx, 1/dy, 1/dz
    dx_i (:) = 1.0_dp / dx_i (:)
    dy_i (:) = 1.0_dp / dy_i (:)
    dz_i (:) = 1.0_dp / dz_i (:)


  end subroutine metrics


!> \brief Calculate constants used in the simulation.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine constants ( adi )


    type (adi_type) , intent (inout) :: adi !< non-dimensional derived type


    ! Gamma definitions
    adi % gamma_inf = 1.40_dp
    adi % gamma     = adi % gamma_inf
    adi % gm1       = adi % gamma - 1.0_dp
    adi % gm1_i     = 1.0_dp / adi % gm1


    ! Non dimensional numbers
    adi % re     = 65000_dp
    adi % pr     = 0.72_dp
    adi % le     = 1.0_dp
    adi % ma     = 2.25_dp
    adi % ttrd   = 2.0_dp / 3.0_dp
    adi % ggmopr = adi % gamma * adi % gm1_i / adi % Pr
    adi % sqgmr  = sqrt (adi % gamma) * adi % ma / adi % re
    adi % Sc     = adi % Le * adi % Pr


    ! Infinite and reference paremeters


    ! Reference temperature

    adi % T_inf   = 169.44_dp
    adi % T_ref   = adi % T_inf

    ! Reference lentgh
    adi % L_ref   = 25.4e-3_dp

    ! Reference Universal Gas Constant
    adi % R_ref   = 8.31451_dp
    adi % R_inf   = adi % R_ref
    adi % r_m_ref = 287.15_dp
    adi % r_m_inf = adi % r_m_ref

    ! Reference molecular weight
    adi % W_ref   = adi % R_ref / adi % r_m_ref
    adi % W_inf   = adi % W_ref

    ! Reference specific heat
    adi % cp_ref  = adi % gamma_inf * &
                  ( adi % R_ref / adi % W_ref ) / &
                  ( adi % gamma_inf - 1.0_dp )
    adi % cp_inf  = adi % cp_ref

    ! Reference viscosity ( Shutherland's law )
    if ( adi % T_ref < 110.4_dp ) then
       adi % mu_inf = vand (adi % T_ref)
       adi % mu_ref = adi % mu_inf
    else
       adi % mu_inf = suth (adi % T_ref)
       adi % mu_ref = adi % mu_inf
    end if

    ! Infinity speed of sound
    adi % c_inf = sqrt ( adi % gamma_inf *             &
                       ( adi % R_inf / adi % W_inf ) * &
                         adi % T_inf )

    ! Reference velocity
    adi % u_ref = adi % c_inf / sqrt ( adi % gamma_inf )
    adi % u_inf = adi % c_inf * adi % ma

    ! Reference time
    adi % time_ref = adi % L_ref / adi % u_ref

    ! Reference density
    adi % rho_ref = adi % re * adi % mu_ref / &
                  ( adi % L_ref * adi % u_inf )
    adi % rho_inf = adi % rho_ref

    ! Reference pressure
    ! here p_ref matches p_inf although the way each one is calculated
    ! is different and therefore there are small round-off errors
    adi % p_ref = adi % rho_ref * ( adi % u_ref * adi % u_ref )
    adi % p_inf = adi % rho_inf * ( adi % R_inf / adi % W_inf ) * &
                  adi % t_inf

    ! Reference lambda
    adi % lbda_ref = adi % cp_ref * adi % mu_ref / &
                     adi % pr
    adi % lbda_inf = adi % lbda_ref

    ! Reference diffusion coefficient
    adi % D_ref = adi % lbda_ref / ( adi % rho_ref * &
                  adi % cp_ref * adi % Le )
    adi % D_inf = adi % D_ref

    ! ! Display adi type
    ! if ( rank == rank_default )  then
    !    call adi_display (adi)
    !    call adi_write (adi)
    ! end if


  end subroutine constants


!> \brief Calculate viscosity using Sutherland's law.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function suth (T)


    real (dp) , intent (in) :: T !< temperature
    real (dp) :: suth !< viscosity


    suth = 1.4580e-6_dp * T ** (1.5_dp)  / ( T + 110.4_dp ) ! PM fix


  end function suth


!> \brief Calculate viscosity using Van's law.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  function vand (T)


    real (dp) , intent (in) :: T !< temperature
    real (dp) :: vand !< viscosity


    vand = 0.6938730e-6_dp * T


  end function vand


!> \brief Estimate times for the simulation.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine counter ( total_progress , sub_progress , sub_time_invested , comm_time , phys_time )


    real (dp) , intent (in)                            :: total_progress    !< total progress of the simulation
    real (dp) , intent (in)                            :: sub_progress      !< current progress of the simulation
    real (dp) , intent (in)                            :: sub_time_invested !< time invested in the current progress
    real (dp) , intent (in)                            :: comm_time         !< time invisted in communication in the current progress
    real (dp) , intent (in)                            :: phys_time         !< total time of the simulation


    character ( len = 50 ) , parameter                 :: format_time = '(I8,A,1X,F9.2,1X,A,1X,F9.2,1X,A,1X,F9.2,A)'

    real (dp) , dimension (3) , parameter              :: time_cuts = (/ 86400.0_dp , 3600.0_dp , 60.0_dp /)

    character ( len = 50 ) , dimension (4)             :: phrase_time

    character ( len = 50 )                             :: phrase_phys_time

    real (dp)                                          :: percent_progress , remaining_time , phys_time_dim , parallel_perf


    phrase_time (1) = 'days running,'
    phrase_time (2) = 'hours running,'
    phrase_time (3) = 'minutes running,'
    phrase_time (4) = 'seconds running,'

    remaining_time   = ( 1.0_dp - total_progress ) * sub_time_invested / sub_progress
    percent_progress = 100.0_dp * total_progress

    parallel_perf    = 100.0_dp * ( 1.0_dp - comm_time / sub_time_invested )

    if ( phys_time > time_cuts (1) ) then
       phys_time_dim    = phys_time / time_cuts (1)
       phrase_phys_time = phrase_time (1)
    else if ( phys_time > time_cuts (2) ) then
       phys_time_dim    = phys_time / time_cuts (2)
       phrase_phys_time = phrase_time (2)
    else if ( phys_time > time_cuts (3) ) then
       phys_time_dim    = phys_time / time_cuts (3)
       phrase_phys_time = phrase_time (3)
    else
       phys_time_dim    = phys_time
       phrase_phys_time = phrase_time (4)
    end if


    if ( remaining_time > time_cuts (1) ) then

       remaining_time = remaining_time / time_cuts (1)
       write ( * , format_time ) int (percent_progress) , '% completed,' ,  &
                                 phys_time_dim  , trim (phrase_phys_time) , &
                                 remaining_time , 'days left,' ,         &
                                 parallel_perf , '% parallel efficiency'

    else if ( remaining_time > time_cuts (2) ) then

       remaining_time = remaining_time / time_cuts (2)
       write ( * , format_time ) int (percent_progress) , '% completed,' ,  &
                                 phys_time_dim  , trim (phrase_phys_time) , &
                                 remaining_time , 'hours left,' ,        &
                                 parallel_perf , '% parallel efficiency'

    else if ( remaining_time > time_cuts (3) ) then

       remaining_time = remaining_time / time_cuts (3)
       write ( * , format_time ) int (percent_progress) , '% completed,' ,  &
                                 phys_time_dim  , trim (phrase_phys_time) , &
                                 remaining_time , 'minutes left,' ,      &
                                 parallel_perf , '% parallel efficiency'

    else

       write ( * , format_time ) int (percent_progress) , '% completed,' ,  &
                                 phys_time_dim  , trim (phrase_phys_time) , &
                                 remaining_time , 'seconds left,' ,      &
                                 parallel_perf , '% parallel efficiency'

    end if


  end subroutine counter


!   subroutine saving_ascii ( time , dt , adi , x , y , z , T , W_i , v )


!     real (dp) , intent (in)                                     :: time
!     real (dp) , intent (in)                                     :: dt
!     type (adi_type) , intent (in)                               :: adi
!     real (dp) , allocatable , dimension (:) , intent (in)       :: x , y , z
!     real (dp) , allocatable , dimension (:,:,:) , intent (in)   :: T , W_i
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (in) :: v


!     real (dp) :: P , rho_i


!     P     = v (sx,sy,sz,1) * T (sx,sy,sz) * W_i (sx,sy,sz)
!     rho_i = 1.0_dp / v (sx,sy,sz,1)


!     write ( unit_saveascii , * ) time * adi % time_ref , &
!                                  dt * adi % time_ref , &
!                                  v (sx,sy,sz,1) * adi % rho_ref, &
!                                  v (sx,sy,sz,2) * rho_i * adi % u_ref , &
!                                  v (sx,sy,sz,3) * rho_i * adi % u_ref , &
!                                  v (sx,sy,sz,4) * rho_i * adi % u_ref , &
!                                  v (sx,sy,sz,5) * rho_i * adi % u_ref * adi % u_ref , &
!                                  T (sx,sy,sz) * adi % T_ref , &
!                                  P * adi % p_ref , &
!                                  v (sx,sy,sz,6:nv) * rho_i


!   end subroutine saving_ascii


!   subroutine saving_ascii2 ( time , dt , thd , adi , x , y , z , T , W_i , cp , v )


!     real (dp) , intent (in)                                     :: time
!     real (dp) , intent (in)                                     :: dt
!     type (thd_type) , intent (in)                               :: thd
!     type (adi_type) , intent (in)                               :: adi
!     real (dp) , allocatable , dimension (:) , intent (in)       :: x , y , z
!     real (dp) , allocatable , dimension (:,:,:) , intent (in)   :: T , W_i , cp
!     real (dp) , allocatable , dimension (:,:,:,:) , intent (in) :: v


!     integer (ip)                                                :: ok , i , j , k , l

!     integer (ip)                                                :: i0 , i1 , &
!                                                                    j0 , j1 , &
!                                                                    k0 , k1

!     real (dp) , allocatable , dimension (:,:,:)                 :: rho_i , P , gam , cs


!     allocate ( rho_i (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
!                cs    (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
!                P     (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
!                gam   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
!                stat = ok )
!     if ( ok > 0 ) call abort_mpi ('error allocate saving_ascii2')


!     do l = -ndim , ndim

!        call domain_select ( l , i0 , i1 , j0 , j1 , k0 , k1 )

!        do k = k0 , k1
!           do j = j0 , j1
!              do i = i0 , i1

!                 rho_i (i,j,k) = 1.0_dp / v (i,j,k,1)
!                 P (i,j,k)     = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
!                 gam (i,j,k)   = thd % gam2 * W_i (i,j,k) / cp (i,j,k)
!                 gam (i,j,k)   = 1.0_dp / ( 1.0_dp - gam (i,j,k) )
!                 cs  (i,j,k)   = sqrt ( gam (i,j,k) * P (i,j,k) * rho_i (i,j,k) )

!              end do
!           end do
!        end do

!     end do


!     ! do i = sx , ex
!     !    do j = sy , ey

!     !       write ( unit_saveascii , '(50(1X,1PE20.12))' ) x (i) * adi % L_ref , &
!     !                                                      y (j) * adi % L_ref , &
!     !                                                      v (i,j,sz,1) * adi % rho_ref , &
!     !                                                      v (i,j,sz,2) * rho_i (i,j,sz) * adi % u_ref , &
!     !                                                      v (i,j,sz,3) * rho_i (i,j,sz) * adi % u_ref , &
!     !                                                      v (i,j,sz,4) * rho_i (i,j,sz) * adi % u_ref , &
!     !                                                      v (i,j,sz,5) * rho_i (i,j,sz) * adi % u_ref * adi % u_ref , &
!     !                                                      T (i,j,sz) * adi % T_ref , &
!     !                                                      P (i,j,sz) * adi % p_ref , &
!     !                                                      cp(i,j,sz) * adi % cp_ref , &
!     !                                                      gam(i,j,sz) , &
!     !                                                      cs(i,j,sz) * adi % u_ref , &
!     !                                                      W_i (i,j,sz) / adi % W_ref , &
!     !                                                      v (i,j,sz,6:nv) * rho_i (i,j,sz)

!     !    end do

!     !    write ( unit_saveascii , * )

!     ! end do


!     do i = sx , ex

!        write ( unit_saveascii , '(50(1X,1PE20.12))' ) time * adi % time_ref , &
!                                                       x (i) * adi % L_ref , &
!                                                       v (i,sy,sz,1) * adi % rho_ref , &
!                                                       v (i,sy,sz,2) * rho_i (i,sy,sz) * adi % u_ref , &
!                                                       v (i,sy,sz,3) * rho_i (i,sy,sz) * adi % u_ref , &
!                                                       v (i,sy,sz,4) * rho_i (i,sy,sz) * adi % u_ref , &
!                                                       v (i,sy,sz,5) * rho_i (i,sy,sz) * adi % u_ref * adi % u_ref , &
!                                                       T (i,sy,sz) * adi % T_ref , &
!                                                       P (i,sy,sz) * adi % p_ref , &
!                                                       cp(i,sy,sz) * adi % cp_ref , &
!                                                       gam(i,sy,sz) , &
!                                                       cs(i,sy,sz) * adi % u_ref , &
!                                                       W_i (i,sy,sz) / adi % W_ref , &
!                                                       v (i,sy,sz,6:nv) * rho_i (i,sy,sz)

!        end do


!        deallocate ( rho_i , cs , P , gam )


!   end subroutine saving_ascii2


!> \brief Save restart file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine save_restart ( number , time , dtime , dt_4 , adi , v )


    integer (ip) , intent (in)                                  :: number !< number of restart file
    real (dp) , intent (in)                                     :: time   !< time
    real (dp) , intent (in)                                     :: dtime  !< time step
    real (dp) , dimension (:) , intent (in)                     :: dt_4   !< array of different time steps
    type (adi_type) , intent (in)                               :: adi    !< non-dimensional derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in) :: v      !< conserved variables array


    character (len_default) , parameter                         :: format_exit = '(I8,1X,10(1X,1PE18.8))'
    character (len_default)                                     :: restart_ascii , rank_ascii , adr_file
    integer (ip)                                                :: ok , i , j , k , l
    integer (kind=8)                                            :: reclmax , rec_index


    ! timing file

    if ( rank == rank_default ) then
       adr_file = trim (dir_parent) // trim (file_time_rest)
       open ( unit_time_rest , file = adr_file , status = 'old' , position = 'append' , iostat = ok )
       if ( ok /= 0 ) call abort_mpi ('ERROR opening ' // trim (adr_file))
       write ( unit_time_rest , format_exit ) number , time * adi % time_ref , dtime * adi % time_ref , dt_4 / dtime
       close ( unit_time_rest )
    end if


    write ( restart_ascii , format_restart ) number
    adr_file = trim (dir_parent) // trim (dir_restart) // trim (file_restart) // '_' // trim (restart_ascii)
    if ( rank == rank_default ) write (*,*) 'writing file ' , trim(adr_file)


    if (ind_files) then ! writing individual binary files per process

       reclmax = int ( (ex-sx+1) * (ey-sy+1) * (ez-sz+1) * nv * nbit , kind = 8 )

       write ( rank_ascii , format_restart ) rank
       adr_file = trim (adr_file) // '_' // trim (rank_ascii)

       rec_index = 1

    else ! writing one single binary file

       reclmax = int ( nxmax * nymax * nzmax * nv * nbit , kind = 8 )

       rec_index = rank + 1

    end if


    open ( unit_restart , file   = trim (adr_file) , &
                          access = 'direct'        , &
                          recl   = reclmax         , &
                          form   = 'unformatted'   , &
                          status = 'unknown'       , &
                          iostat = ok              )

    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (adr_file))

    write ( unit_restart , rec = rec_index )&
         (((( v (i,j,k,l) , i = sx , ex ) , &
                            j = sy , ey ) , &
                            k = sz , ez ) , &
                            l = 1  , nv )

    close (unit_restart)


  end subroutine save_restart


!> \brief Locate statistical planes or volumes.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine locate_stats ( inp , adi , grid )


    type (inp_type) , intent (inout)                            :: inp  !< input derived type
    type (adi_type) , intent (in)                               :: adi  !< non-dimensional derived type
    type (inp_grid) , intent (in)                               :: grid !< grid derived type


    integer (ip)                                                :: i , j , k
    integer (ip)                                                :: volume , plane , rank_plane
    integer (ip) , dimension (ndimmax)                          :: coordmin , coordmax , coordmin2 , coordmax2
    logical                                                     :: log_volx , log_voly , log_volz
    character (len_default) , parameter                         :: format_exit = '(I7,1X,I3,1X,A3,1X,3(1X,1PE18.8))'
    real (dp)                                                   :: L_ref


    L_ref = adi % L_ref


    ! volumes


    if ( rank == rank_default .and. inp % nvolume > 0 ) &
         write (*,*) 'volumes ************************************************************'

    log_volume (:) = .false.

    do volume = 1 , inp % nvolume

       log_volx = .false. ; log_voly = .false. ; log_volz = .false.

       coordmin (3) = coords (3) ; coordmax (3) = coords (3)
       coordmin (2) = coords (2) ; coordmax (2) = coords (2)
       coordmin (1) = coords (1) ; coordmax (1) = coords (1)

       ! x-direction
       if ( inp % x_volmin (volume) <= grid % x (sx) .and. &
            inp % x_volmax (volume) >= grid % x (ex) ) then

          log_volx = .true. ! entirely included

       else if ( inp % x_volmin (volume) >= grid % x (sx) .and. &
                 inp % x_volmin (volume) <= grid % x (ex) ) then

          log_volx = .true. ! WEST boundary

       else if ( inp % x_volmax (volume) >= grid % x (sx) .and. &
                 inp % x_volmax (volume) <= grid % x (ex) ) then

          log_volx = .true. ! EST boundary

       else if ( ( inp % x_volmin (volume) > grid % x (sx) .and.  &
                   inp % x_volmax (volume) > grid % x (ex) ) .or. &
                 ( inp % x_volmin (volume) < grid % x (sx) .and.  &
                   inp % x_volmax (volume) < grid % x (ex) ) ) then

          ! excluded

       else

          write (*,*) rank , 'there is a problem with volumes in x-direction'
          call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

       end if


       ! y-direction
       if ( inp % y_volmin (volume) <= grid % y (sy) .and. &
            inp % y_volmax (volume) >= grid % y (ey) ) then

          log_voly = .true. ! entirely included

       else if ( inp % y_volmin (volume) >= grid % y (sy) .and. &
                 inp % y_volmin (volume) <= grid % y (ey) ) then

          log_voly = .true. ! SOUTH boundary

       else if ( inp % y_volmax (volume) >= grid % y (sy) .and. &
                 inp % y_volmax (volume) <= grid % y (ey) ) then

          log_voly = .true. ! NORTH boundary

       else if ( ( inp % y_volmin (volume) > grid % y (sy) .and.  &
                   inp % y_volmax (volume) > grid % y (ey) ) .or. &
                 ( inp % y_volmin (volume) < grid % y (sy) .and.  &
                   inp % y_volmax (volume) < grid % y (ey) ) ) then

          ! excluded

       else

          write (*,*) rank , 'there is a problem with volumes in y-direction'
          call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

       end if


       ! z-direction
       if ( inp % z_volmin (volume) <= grid % z (sz) .and. &
            inp % z_volmax (volume) >= grid % z (ez) ) then

          log_volz = .true. ! entirely included

       else if ( inp % z_volmin (volume) >= grid % z (sz) .and. &
                 inp % z_volmin (volume) <= grid % z (ez) ) then

          log_volz = .true. ! BEHIND boundary

       else if ( inp % z_volmax (volume) >= grid % z (sz) .and. &
                 inp % z_volmax (volume) <= grid % z (ez) ) then

          log_volz = .true. ! FRONT boundary

       else if ( ( inp % z_volmin (volume) > grid % z (sz) .and.  &
                   inp % z_volmax (volume) > grid % z (ez) ) .or. &
                 ( inp % z_volmin (volume) < grid % z (sz) .and.  &
                   inp % z_volmax (volume) < grid % z (ez) ) ) then

          ! excluded

       else

          write (*,*) rank , 'there is a problem with volumes in z-direction'
          call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

       end if

       if ( log_volx .and. log_voly .and. log_volz ) log_volume (volume) = .true.

       if ( .not. log_volume (volume) ) then
          coordmin (3) = dims(3)+1 ; coordmax (3) = -1
          coordmin (2) = dims(2)+1 ; coordmax (2) = -1
          coordmin (1) = dims(1)+1 ; coordmax (1) = -1
       end if

       call mpi_allreduce ( coordmin , coordmin2 , ndimmax , MPI_INTEGER , & ! communicate the minimum
                            MPI_MIN , MPI_COMM_WORLD , mpicode )

       call mpi_allreduce ( coordmax , coordmax2 , ndimmax , MPI_INTEGER , & ! communicate the maximum
                            MPI_MAX , MPI_COMM_WORLD , mpicode )

       coord_vol_min (volume,3) = coordmin2 (3) ; coord_vol_max (volume,3) = coordmax2 (3)
       coord_vol_min (volume,2) = coordmin2 (2) ; coord_vol_max (volume,2) = coordmax2 (2)
       coord_vol_min (volume,1) = coordmin2 (1) ; coord_vol_max (volume,1) = coordmax2 (1)

       if ( coords (3) == coord_vol_min (volume,3) .and. &
            coords (2) == coord_vol_min (volume,2) .and. &
            coords (1) == coord_vol_min (volume,1) ) then ! locating the first rank on the volume
          write (*,'(I7,1X,I3,3(1X,I7))') rank , volume , coordmax2 (3) - coordmin2 (3) + 1 , &
                                                 coordmax2 (2) - coordmin2 (2) + 1 , &
                                                 coordmax2 (1) - coordmin2 (1) + 1
          write (*,'(I7,1X,I3,3(1X,1PE18.8))') &
             rank , volume , grid % x(sx) * L_ref , grid % y(sy) * L_ref , grid % z(sz) * L_ref
       end if

       if ( coords (3) == coord_vol_max (volume,3) .and. &
            coords (2) == coord_vol_max (volume,2) .and. &
            coords (1) == coord_vol_max (volume,1) ) then ! locating the last rank on the volume
          write (*,'(I7,1X,I3,3(1X,1PE18.8))') &
             rank , volume , grid % x(ex) * L_ref , grid % y(ey) * L_ref , grid % z(ez) * L_ref
       end if

       call mpi_barrier ( MPI_COMM_WORLD , mpicode )

    end do


    call mpi_barrier ( MPI_COMM_WORLD , mpicode )


    if ( rank == rank_default .and. inp % nvolume > 0 ) then
       write (*,*) 'volumes ************************************************************'
       write (*,*)
    end if


    ! XY planes


    log_planeXY (:) = .false.
    do plane = 1 , inp % nxystat

       if ( inp % z_xystat (plane) >= grid % z (sz) .and. &
            inp % z_xystat (plane) <  grid % z (ez+1) ) then ! A.Techer modif <= ez to < ez+1

          log_planeXY (plane) = .true.

          do k = ez , sz , -1
             if ( inp % z_xystat (plane) >= grid % z (k) .and. &
                  inp % z_xystat (plane) <  grid % z (k+1) ) inp % k_xystat (plane) = k ! A.Techer modif
          end do

          if ( inp % k_xystat (plane) == sz ) inp % k_xystat (plane) = sz+1
          if ( inp % k_xystat (plane) == ez ) inp % k_xystat (plane) = ez-1

       else if ( inp % z_xystat (plane) <  grid % z (sz) .or. &
                 inp % z_xystat (plane) >= grid % z (ez+1) ) then

          ! this rank does not contain any planes

       else

          write (*,*) rank , 'there is a problem with XY planes'
          call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

       end if

    end do


    if ( rank == rank_default .and. inp % nxystat > 0 ) &
         write (*,*) 'XY planes ************************************************************'


    ! reference the planes
    call mpi_barrier ( MPI_COMM_WORLD , mpicode )

    rank_plane = ( coords (3) + 1 ) * ( coords (2) + 1 )

    if ( rank_plane == 1 ) then

       do plane = 1 , inp % nxystat

          if ( inp % k_xystat (plane) - ng-1 >= sz-ng .and. &
               inp % k_xystat (plane) + ng+1 <= ez+ng ) then
             write ( * , format_exit ) rank , plane , 'Z' , grid % z ( inp % k_xystat (plane) - ng-1 ) * L_ref , &
                                                            grid % z ( inp % k_xystat (plane)        ) * L_ref , &
                                                            grid % z ( inp % k_xystat (plane) + ng+1 ) * L_ref
          end if

       end do

    end if


    call mpi_barrier ( MPI_COMM_WORLD , mpicode )


    if ( rank == rank_default .and. inp % nxystat > 0 ) then
       write (*,*) 'XY planes ************************************************************'
       write (*,*)
    end if


    call mpi_barrier ( MPI_COMM_WORLD , mpicode )


    ! XZ planes


    log_planeXZ (:) = .false.
    do plane = 1 , inp % nxzstat

       if ( inp % y_xzstat (plane) >= grid % y (sy) .and. &
            inp % y_xzstat (plane) <  grid % y (ey+1) ) then ! A.Techer modif

          log_planeXZ (plane) = .true.

          do j = ey , sy , -1
             if ( inp % y_xzstat (plane) >= grid % y (j) .and. &
                  inp % y_xzstat (plane) <  grid % y (j+1) ) inp % j_xzstat (plane) = j
          end do

          if ( inp % j_xzstat (plane) == sy ) inp % j_xzstat (plane) = sy+1
          if ( inp % j_xzstat (plane) == ey ) inp % j_xzstat (plane) = ey-1

       else if ( inp % y_xzstat (plane) <  grid % y (sy) .or. &
                 inp % y_xzstat (plane) >= grid % y (ey+1) ) then

          ! this rank does not contain any planes

       else

          write (*,*) rank , 'there is a problem with XZ planes'
          call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

       end if

    end do


    if ( rank == rank_default .and. inp % nxzstat > 0 ) &
         write (*,*) 'XZ planes ************************************************************'


    ! reference the planes
    call mpi_barrier ( MPI_COMM_WORLD , mpicode )

    rank_plane = ( coords (3) + 1 ) * ( coords (1) + 1 )

    if ( rank_plane == 1 ) then

       do plane = 1 , inp % nxzstat

          if ( inp % j_xzstat (plane) - ng-1 >= sy-ng .and. &
               inp % j_xzstat (plane) + ng+1 <= ey+ng ) then
             write ( * , format_exit ) rank , plane , 'Y' , grid % y ( inp % j_xzstat (plane) - ng-1 ) * L_ref , &
                                                            grid % y ( inp % j_xzstat (plane)        ) * L_ref , &
                                                            grid % y ( inp % j_xzstat (plane) + ng+1 ) * L_ref

          end if

       end do

    end if


    call mpi_barrier ( MPI_COMM_WORLD , mpicode )


    if ( rank == rank_default .and. inp % nxzstat > 0 ) then
       write (*,*) 'XZ planes ************************************************************'
       write (*,*)
    end if


    call mpi_barrier ( MPI_COMM_WORLD , mpicode )


    ! YZ planes


    log_planeYZ (:) = .false.
    do plane = 1 , inp % nyzstat

       if ( inp % x_yzstat (plane) >= grid % x (sx) .and. &
            inp % x_yzstat (plane) <  grid % x (ex+1) ) then

          log_planeYZ (plane) = .true.

          do i = ex , sx , -1
             if ( inp % x_yzstat (plane) >= grid % x (i) .and. &
                  inp % x_yzstat (plane) <  grid % x (i+1) ) inp % i_yzstat (plane) = i
          end do

          if ( inp % i_yzstat (plane) == sx ) inp % i_yzstat (plane) = sx+1
          if ( inp % i_yzstat (plane) == ex ) inp % i_yzstat (plane) = ex-1

       else if ( inp % x_yzstat (plane) <  grid % x (sx) .or. &
                 inp % x_yzstat (plane) >= grid % x (ex+1) ) then

          ! this rank does not contain any planes

       else

          write (*,*) rank , 'there is a problem with YZ planes'
          call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

       end if

    end do


    if ( rank == rank_default .and. inp % nyzstat > 0 ) &
         write (*,*) 'YZ planes ************************************************************'


    ! reference the planes
    call mpi_barrier ( MPI_COMM_WORLD , mpicode )

    rank_plane = ( coords (2) + 1 ) * ( coords (1) + 1 )

    if ( rank_plane == 1 ) then

       do plane = 1 , inp % nyzstat

          if ( inp % i_yzstat (plane) - ng-1 >= sx-ng .and. &
               inp % i_yzstat (plane) + ng+1 <= ex+ng ) then
             write ( * , format_exit ) rank , plane , 'X' , grid % x ( inp % i_yzstat (plane) - ng-1 ) * L_ref , &
                                                            grid % x ( inp % i_yzstat (plane)        ) * L_ref , &
                                                            grid % x ( inp % i_yzstat (plane) + ng+1 ) * L_ref
          end if

       end do

    end if


    call mpi_barrier ( MPI_COMM_WORLD , mpicode )


    if ( rank == rank_default .and. inp % nyzstat > 0 ) then
       write (*,*) 'YZ planes ************************************************************'
       write (*,*)
    end if


  end subroutine locate_stats


!> \brief Save statistical planes or volumes.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine save_stats ( number , time , dtime , dt_4 , inp , adi , v )


    integer (ip) , intent (in)                                  :: number !< number of stastical file
    real (dp) , intent (in)                                     :: time   !< time
    real (dp) , intent (in)                                     :: dtime  !< time step
    real (dp) , dimension (:) , intent (in)                     :: dt_4   !< array of different time steps
    type (inp_type) , intent (inout)                            :: inp    !< input derived type
    type (adi_type) , intent (in)                               :: adi    !< non-dimensional derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in) :: v      !< conserved variables array


    character (len_default) , parameter                         :: format_exit = '(I8,1X,10(1X,1PE18.8))'
    character (len_default)                                     :: stat_ascii , plane_ascii , rank_ascii , adr_file
    integer (ip)                                                :: ok , volume , plane , rank_volume , rank_plane
    integer (ip) , dimension (ndimmax)                          :: dimmaxvol
    integer (kind=8)                                            :: reclmax , rec_index
    integer (ip)                                                :: sx0 , ex0 , &
                                                                   sy0 , ey0 , &
                                                                   sz0 , ez0
    integer (ip)                                                :: i , j , k , l


    ! timing file

    if ( rank == rank_default ) then
       adr_file = trim (dir_parent) // trim (file_time_stat)
       open ( unit = unit_time_stat , file = adr_file , status = 'old' , position = 'append' , iostat = ok )
       if ( ok /= 0 ) call abort_mpi ('ERROR opening ' // trim (adr_file))
       write ( unit_time_stat , format_exit ) number , time * adi % time_ref , dtime * adi % time_ref , dt_4 / dtime
       close ( unit_time_stat )
    end if

    write ( stat_ascii , format_restart ) number
    write ( rank_ascii , format_restart ) rank


    ! volumes


    do volume = 1 , inp % nvolume

       write ( plane_ascii , format_nplane ) volume
       adr_file = trim (dir_parent) // trim (dir_statvol) // trim (file_statvol) &
       &          // '_' // trim (plane_ascii) // '_' // trim (stat_ascii)

       if ( coords (3) == coord_vol_min (volume,3) .and. &
            coords (2) == coord_vol_min (volume,2) .and. &
            coords (1) == coord_vol_min (volume,1) ) then ! locating the first rank on the volume
          write (*,*) 'writing file ' , trim (adr_file)
       end if


       if ( log_volume (volume) ) then


          if (ind_files) then ! writing individual binary files per process

             reclmax = int ( (ex-sx+1) * (ey-sy+1) * (ez-sz+1) * nv * nbit , kind = 8 )

             adr_file = trim (adr_file) // '_' // trim (rank_ascii)

             rec_index = 1

          else ! writing one single binary file

             reclmax = int ( nxmax * nymax * nzmax * nv * nbit , kind = 8 )

             dimmaxvol (3) = coord_vol_max (volume,3) - coord_vol_min (volume,3) + 1
             dimmaxvol (2) = coord_vol_max (volume,2) - coord_vol_min (volume,2) + 1
             dimmaxvol (1) = coord_vol_max (volume,1) - coord_vol_min (volume,1) + 1

             rank_volume = ( coords (3) - coord_vol_min (volume,3) ) +                                   &
                           ( coords (2) - coord_vol_min (volume,2) ) * dimmaxvol (3) +                   &
                           ( coords (1) - coord_vol_min (volume,1) ) * ( dimmaxvol (3) * dimmaxvol (2) )
             rec_index = rank_volume + 1

          end if


          open ( unit_statvol + volume , file   = trim (adr_file) , &
                                         access = 'direct'        , &
                                         recl   = reclmax         , &
                                         form   = 'unformatted'   , &
                                         status = 'unknown'       , &
                                         iostat = ok              )

          if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (adr_file))

          write ( unit_statvol + volume , rec = rec_index ) &
                (((( v (i,j,k,l) , i = sx , ex ) , &
                                   j = sy , ey ) , &
                                   k = sz , ez ) , &
                                   l = 1  , nv )

          close ( unit_statvol + volume )


       end if


    end do


    ! XY planes


    do plane = 1 , inp % nxystat

       write ( plane_ascii , format_nplane ) plane
       adr_file = trim (dir_parent) // trim (dir_statXY) // trim (file_statXY) &
       &          // '_' // trim (plane_ascii) // '_' // trim (stat_ascii)

       if ( rank == rank_default ) write (*,*) 'writing file ' , trim (adr_file)

       if ( log_planeXY (plane) ) then

          sz0 = inp % k_xystat (plane) - ng-1
          ez0 = inp % k_xystat (plane) + ng+1

          if (ind_files) then ! writing individual binary files per process

             reclmax = int ( (ex-sx+1) * (ey-sy+1) * ( ng+3+ng ) * nv * nbit , kind = 8 )

             adr_file = trim (adr_file) // '_' // trim (rank_ascii)

             rec_index = 1

          else ! writing one single binary file

             reclmax = int ( nxmax * nymax * ( ng+3+ng ) * nv * nbit , kind = 8 )

             rank_plane = coords (3) + coords (2) * dims (3) + 1
             rec_index = rank_plane

          end if


          open ( unit_statXY + plane , file   = trim (adr_file) , &
                                       access = 'direct'        , &
                                       recl   = reclmax         , &
                                       form   = 'unformatted'   , &
                                       status = 'unknown'       , &
                                       iostat = ok              )

          if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (adr_file))

          write ( unit_statXY + plane , rec = rec_index ) &
               (((( v (i,j,k,l) , i = sx  , ex  ) , &
                                  j = sy  , ey  ) , &
                                  k = sz0 , ez0 ) , &
                                  l = 1   , nv  )

          close ( unit_statXY + plane )


       end if


    end do


    ! XZ planes


    do plane = 1 , inp % nxzstat

       write ( plane_ascii , format_nplane ) plane
       adr_file = trim (dir_parent) // trim (dir_statXZ) // trim (file_statXZ) &
       &          // '_' // trim (plane_ascii) // '_' // trim (stat_ascii)

       if ( rank == rank_default ) write (*,*) 'writing file ' , trim (adr_file)

       if ( log_planeXZ (plane) ) then

          sy0 = inp % j_xzstat (plane) - ng-1
          ey0 = inp % j_xzstat (plane) + ng+1

          if (ind_files) then ! writing individual binary files per process

             reclmax = int ( (ex-sx+1) * ( ng+3+ng ) * (ez-sz+1) * nv * nbit , kind = 8 )

             adr_file = trim (adr_file) // '_' // trim (rank_ascii)

             rec_index = 1

          else ! writing one single binary file

             reclmax = int ( nxmax * ( ng+3+ng ) * nzmax * nv * nbit , kind = 8 )

             rank_plane = coords (3) + coords (1) * dims (3) + 1
             rec_index = rank_plane

          end if


          open ( unit_statXZ + plane , file   = trim (adr_file) , &
                                       access = 'direct'        , &
                                       recl   = reclmax         , &
                                       form   = 'unformatted'   , &
                                       status = 'unknown'       , &
                                       iostat = ok              )

          if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (adr_file))

          write ( unit_statXZ + plane , rec = rec_index ) &
               (((( v (i,j,k,l) , i = sx  , ex  ) , &
                                  j = sy0 , ey0 ) , &
                                  k = sz  , ez  ) , &
                                  l = 1   , nv  )

          close ( unit_statXZ + plane )


       end if


    end do


    ! YZ planes


    do plane = 1 , inp % nyzstat

       write ( plane_ascii , format_nplane ) plane
       adr_file = trim (dir_parent) // trim (dir_statYZ) // trim (file_statYZ) &
       &          // '_' // trim (plane_ascii) // '_' // trim (stat_ascii)

       if ( rank == rank_default ) write (*,*) 'writing file ' , trim (adr_file)

       if ( log_planeYZ (plane) ) then

          sx0 = inp % i_yzstat (plane) - ng-1
          ex0 = inp % i_yzstat (plane) + ng+1

          if (ind_files) then ! writing individual binary files per process

             reclmax = int ( ( ng+3+ng ) * (ey-sy+1) * (ez-sz+1) * nv * nbit , kind = 8 )

             adr_file = trim (adr_file) // '_' // trim (rank_ascii)

             rec_index = 1

          else ! writing one single binary file

             reclmax = int ( ( ng+3+ng ) * nymax * nzmax * nv * nbit , kind = 8 )

             rank_plane = coords (2) + coords (1) * dims (2) + 1
             rec_index = rank_plane

          end if


          open ( unit_statYZ + plane , file   = trim (adr_file) , &
                                       access = 'direct'        , &
                                       recl   = reclmax         , &
                                       form   = 'unformatted'   , &
                                       status = 'unknown'       , &
                                       iostat = ok              )

          if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (adr_file))

          write ( unit_statYZ + plane , rec = rec_index ) &
               (((( v (i,j,k,l) , i = sx0 , ex0 ) , &
                                  j = sy  , ey  ) , &
                                  k = sz  , ez  ) , &
                                  l = 1   , nv  )

          close ( unit_statYZ + plane )


       end if


    end do


  end subroutine save_stats


!> \brief Generate random numbers in parallel.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine randomize (inp)


    type (inp_type) , intent (in)                              :: inp !< input derived type


    integer (ip) , dimension (2)                               :: wrk
    integer (ip) , dimension (:) , allocatable                 :: seed

    integer (ip)                                               :: i


    if ( rank == rank_default ) then
       call random_seed  ( size  = wrk(1) )
       call system_clock ( count = wrk(2) )
!       call date_and_time ( values = time )
    end if

    call MPI_BCAST ( wrk , 2 , MPI_INTEGER , 0 , MPI_COMM_WORLD , mpicode )

    allocate ( seed ( wrk(1) ) )

    do i = 1 , wrk(1)

       seed (i) = nv +                                      &
                  nderiv +                                  &
                  ntx * nty * ntz +                         &
                  (dims(1)+1) * (dims(2)+1) * (dims(3)+1) + &
                  inp % walltime +                          &
                  inp % itshow + inp % itmax

       seed (i) = wrk(2) + seed (i) * (i-1) * (nproc+1)

    end do

    call random_seed ( put = seed )
    deallocate (seed)


  end subroutine randomize


!> \brief Locate probes.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine locate_probes ( inp , adi , grid )


    type (inp_type) , intent (in) :: inp  !< input derived type
    type (adi_type) , intent (in) :: adi  !< non-dimensional derived type
    type (inp_grid) , intent (in) :: grid !< grid derived type


    integer (ip)                        :: i , j , k , px , py , pz
    integer (ip)                        :: iprobe
    real (dp)                           :: L_ref
    real (dp) , dimension (ndimmax)     :: pt
    logical                             :: cond
    character (len_default) , parameter :: format = ' (5X,A15,3(1PE18.8)) '

    L_ref = adi % L_ref


    if ( rank == rank_default ) &
       write (*,*) 'Probes (names,x,y,z) *************************************************'

    call mpi_barrier ( MPI_COMM_WORLD , mpicode )

    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

       if ( pt (1) >= grid % x (sx) .and. pt (1) < grid % x (ex+1) .and. &
            pt (2) >= grid % y (sy) .and. pt (2) < grid % y (ey+1) .and. &
            pt (3) >= grid % z (sz) .and. pt (3) < grid % z (ez+1) ) cond = .true.

       if (cond) then

          ! calculate points
          do i = ex , sx , -1
             if ( grid % x (i) <= pt (1) .and. pt (1) < grid % x (i+1 )) px = i
          end do
          px = min(px,ex) ; px = max(px,sx)
          do j = ey , sy , -1
             if ( grid % y (j) <= pt (2) .and. pt (2) < grid % y (j+1) ) py = j
          end do
          py = min(py,ey) ; py = max(py,sy)
          if ( ndim == 3 ) then
             do k = ez , sz , -1
                if ( grid % z (k) <= pt (3) .and. pt (3) < grid % z (k+1) ) pz = k
             end do
          else
             pz = sz
          end if
          pz = min(pz,ez) ; pz = max(pz,sz)

          write ( * , format ) inp % probe_name (iprobe) , &
                               grid % x (px) * L_ref , grid % y (py) * L_ref , grid % z (pz) * L_ref

       end if

       call mpi_barrier ( MPI_COMM_WORLD , mpicode )

    end do

    if ( rank == rank_default ) then
       write (*,*) 'Probes ***************************************************************'
       write (*,*)
    end if


  end subroutine locate_probes


!> \brief Open probes files.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine open_probes ( inp , adi , grid )


    type (inp_type) , intent (in)                               :: inp !< input derived type
    type (adi_type) , intent (in)                               :: adi !< non-dimensional derived type
    type (inp_grid) , intent (in)                               :: grid !< grid derived type


    integer (ip)                        :: i , j , k , px , py , pz
    integer (ip)                        :: currunit , iprobe , ok
    real (dp)                           :: L_ref
    logical                             :: cond
    real (dp) , dimension (ndimmax)     :: pt
    character (len_default)             :: adr_file , ascii
    character (len_default) , parameter :: format = ' ( A10 , 1X , 3 ( 1PE16.8 , 1X ) ) '
    character (len_default) , parameter :: format_head  = ' ( 30 ( A16 , 1X ) ) '
    character (len_default) , parameter :: format_ascii = ' ( I2.2 ) '


    L_ref = adi % L_ref


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

       if ( pt (1) >= grid % x (sx) .and. pt (1) < grid % x (ex+1) .and. &
            pt (2) >= grid % y (sy) .and. pt (2) < grid % y (ey+1) .and. &
            pt (3) >= grid % z (sz) .and. pt (3) < grid % z (ez+1) ) cond = .true.

       if (cond) then

          currunit = unit_probes + iprobe - 1

          adr_file = trim (dir_parent) // trim (dir_probes) // trim (inp % probe_name (iprobe)) // '.out'
          open ( unit   = currunit        , &
                 file   = trim (adr_file) , &
                 form   = 'formatted'     , &
                 status = 'unknown'       , &
                 iostat = ok )
          if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (adr_file))

          ! calculate points
          do i = ex , sx , -1
             if ( grid % x (i) <= pt (1) .and. pt (1) < grid % x (i+1 )) px = i
          end do
          px = min(px,ex) ; px = max(px,sx)
          do j = ey , sy , -1
             if ( grid % y (j) <= pt (2) .and. pt (2) < grid % y (j+1) ) py = j
          end do
          py = min(py,ey) ; py = max(py,sy)
          if ( ndim == 3 ) then
             do k = ez , sz , -1
                if ( grid % z (k) <= pt (3) .and. pt (3) < grid % z (k+1) ) pz = k
             end do
          else
             pz = sz
          end if
          pz = min(pz,ez) ; pz = max(pz,sz)

          write ( currunit , format ) 'coordinates' , grid % x (px) * L_ref , grid % y (py) * L_ref , grid % z (pz) * L_ref
          write ( currunit , format_head , advance='no' ) &
          'time(s)', 'dtime(s)', 'rho(-)', 'rho*u(-)', 'rho*v(-)', 'rho*w(-)', 'rho*et(-)'
          do i = 1 , nrv+npv+nvv
             write ( ascii , format_ascii ) i
             write ( currunit , format_head , advance='no' ) 'rho*Y' // trim(ascii) // '(-)'
          end do
          write ( currunit , * )
          write ( currunit , * )

          close ( unit = currunit ) ! A.Techer

       end if

    end do


  end subroutine open_probes


!> \brief Close probes files.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine close_probes ( inp , grid )


    type (inp_type) , intent (in)                               :: inp !< input derived type
    type (inp_grid) , intent (in)                               :: grid !< grid derived type


    integer (ip)                       :: iprobe
    logical                            :: cond
    real (dp) , dimension (ndimmax)    :: pt


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

       if ( pt (1) >= grid % x (sx) .and. pt (1) < grid % x (ex+1) .and. &
            pt (2) >= grid % y (sy) .and. pt (2) < grid % y (ey+1) .and. &
            pt (3) >= grid % z (sz) .and. pt (3) < grid % z (ez+1) ) cond = .true.

       if (cond) close ( unit = unit_probes + iprobe - 1 )

    end do


  end subroutine close_probes


!> \brief Plot probes files.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine plot_probes ( time , dtime , inp , adi , grid , v )


    real (dp) , intent (in)                                              :: time  !< time
    real (dp) , intent (in)                                              :: dtime !< time step
    type (inp_type) , intent (in)                                        :: inp   !< input derived type
    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    type (inp_grid) , intent (in)                                        :: grid  !< grid derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v     !< conserved variables array


    integer (ip)                        :: i , j , k , px , py , pz
    integer (ip)                        :: iprobe , currunit , ok
    real (dp) , dimension (ndimmax)     :: pt
    logical                             :: cond
    character (len_default)             :: adr_file
    character (len_default) , parameter :: format = ' ( 30 ( 1PE16.8 , 1X ) ) '


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

       if ( pt (1) >= grid % x (sx) .and. pt (1) < grid % x (ex+1) .and. &
            pt (2) >= grid % y (sy) .and. pt (2) < grid % y (ey+1) .and. &
            pt (3) >= grid % z (sz) .and. pt (3) < grid % z (ez+1) ) cond = .true.

       if (cond) then

          ! calculate points
          do i = ex , sx , -1
             if ( grid % x (i) <= pt (1) .and. pt (1) < grid % x (i+1 )) px = i
          end do
          px = min(px,ex) ; px = max(px,sx)
          do j = ey , sy , -1
             if ( grid % y (j) <= pt (2) .and. pt (2) < grid % y (j+1) ) py = j
          end do
          py = min(py,ey) ; py = max(py,sy)
          if ( ndim == 3 ) then
             do k = ez , sz , -1
                if ( grid % z (k) <= pt (3) .and. pt (3) < grid % z (k+1) ) pz = k
             end do
          else
             pz = sz
          end if
          pz = min(pz,ez) ; pz = max(pz,sz)

          currunit = unit_probes + iprobe - 1

          adr_file = trim (dir_parent) // trim (dir_probes) // trim (inp % probe_name (iprobe)) // '.out'
          open ( unit   = currunit        , &
                 file   = trim (adr_file) , &
                 form   = 'formatted'     , &
                 status = 'old'           , &
                 position = 'append'      , &
                 iostat = ok )
          if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (adr_file))

          write ( currunit , format ) time * adi % time_ref ,  &
                                      dtime * adi % time_ref , &
                                      v (px,py,pz,:)

          close ( unit = currunit ) ! A.Techer

       end if

    end do


  end subroutine plot_probes


!> \brief Shock detector.
!!
!! This subroutine is used in the solver.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine shock_det_slv ( v , T , W_i , crit )


    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)                   :: v    !< conserved variables array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)                     :: T    !< temperature
    real (dp) , dimension (:,:,:) , allocatable , intent (in)                     :: W_i  !< inverted molar mass
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (inout)  :: crit !< criteria


    integer (ip) , parameter                  :: stencil_m1 = ng+ng
    integer (ip)                              :: st
    integer (ip)                              :: i , i_l , i_r , i_s , &
                                                 j , j_l , j_r , j_s , &
                                                 k , k_l , k_r , k_s
    real (dp)                                 :: drho , dpres , p_p , p_m


    crit (sx:ex,sy:ey,sz:ez) = 1.0_dp


    if ( ndim >= 1 ) then
       ! X-direction criteria
       do k = sz , ez ! loop in the z-direction
          do j = sy , ey ! loop in the y-direction

             do i_l = sx-1 , ex ! loop on the cell faces

                i_r = i_l + 1

                do st = 1 , stencil_m1

                   i_s = i_l + st - ng

                   ! density criteria
                   drho = abs ( v (min(i_s+1,ex),j,k,1) - v (i_s,j,k,1) ) / &
                          min ( v (min(i_s+1,ex),j,k,1) , v (i_s,j,k,1) )

                   ! pressure criteria
                   p_p   = v (min(i_s+1,ex),j,k,1) * T (min(i_s+1,ex),j,k) * W_i (min(i_s+1,ex),j,k)
                   p_m   = v (i_s,j,k,1) * T (i_s,j,k) * W_i (i_s,j,k)
                   dpres = abs ( p_p - p_m ) / min ( p_p , p_m )

                   if ( drho > percent_weight_SGS .and. dpres > percent_weight_SGS ) crit (i_l,j,k) = 0.0_dp

                end do

                ! activate 2D WENO at the boundaries (not 3D because of periodicity)
                ! if ( i_l < 1+ng )   crit (i_l,j,k) = 1.0_dp
                ! if ( i_l > ntx-ng ) crit (i_l,j,k) = 1.0_dp
                ! if ( j   < 1+ng )   crit (i_l,j,k) = 1.0_dp
                ! if ( j   > nty-ng ) crit (i_l,j,k) = 1.0_dp

             end do

          end do
       end do
    end if


    if ( ndim >= 2 ) then
       ! Y-direction criteria
       do k = sz , ez ! loop in the z-direction
          do i = sx , ex ! loop in the x-direction

             do j_l = sy-1 , ey ! loop on the cell faces

                j_r = j_l + 1

                do st = 1 , stencil_m1

                   j_s = j_l + st - ng

                   ! density criteria
                   drho = abs ( v (i,min(j_s+1,ey),k,1) - v (i,j_s,k,1) ) / &
                          min ( v (i,min(j_s+1,ey),k,1) , v (i,j_s,k,1) )

                   ! pressure criteria
                   p_p   = v (i,min(j_s+1,ey),k,1) * T (i,min(j_s+1,ey),k) * W_i (i,min(j_s+1,ey),k)
                   p_m   = v (i,j_s,k,1) * T (i,j_s,k) * W_i (i,j_s,k)
                   dpres = abs ( p_p - p_m ) / min ( p_p , p_m )

                   if ( drho > percent_weight_SGS .and. dpres > percent_weight_SGS ) crit (i,j_l,k) = 0.0_dp

                end do

                ! activate 2D WENO at the boundaries (not 3D because of periodicity)
                ! if ( i   < 1+ng )   crit (i,j_l,k) = 2.0_dp
                ! if ( i   > ntx-ng ) crit (i,j_l,k) = 2.0_dp
                ! if ( j_l < 1+ng )   crit (i,j_l,k) = 2.0_dp
                ! if ( j_l > nty-ng ) crit (i,j_l,k) = 2.0_dp

             end do

          end do
       end do
    end if


    if ( ndim == 3 ) then
       ! Z-direction criteria
       do j = sy , ey ! loop in the y-direction
          do i = sx , ex ! loop in the x-direction

             do k_l = sz-1 , ez ! loop on the cell faces

                k_r = k_l + 1

                do st = 1 , stencil_m1

                   k_s = k_l + st - ng

                   ! density criteria
                   drho = abs ( v (i,j,min(k_s+1,ez),1) - v (i,j,k_s,1) ) / &
                          min ( v (i,j,min(k_s+1,ez),1) , v (i,j,k_s,1) )

                   ! pressure criteria
                   p_p = v (i,j,min(k_s+1,ez),1) * T (i,j,min(k_s+1,ez)) * W_i (i,j,min(k_s+1,ez))
                   p_m = v (i,j,k_s,1) * T (i,j,k_s) * W_i (i,j,k_s)
                   dpres = abs ( p_p - p_m ) / min ( p_p , p_m )

                   if ( drho > percent_weight_SGS .and. dpres > percent_weight_SGS ) crit (i,j_l,k) = 0.0_dp

                end do

                ! activate 2D WENO at the boundaries (not 3D because of periodicity)
                ! if ( i < 1+ng )   crit (i,j,k_l) = 3.0_dp
                ! if ( i > ntx-ng ) crit (i,j,k_l) = 3.0_dp
                ! if ( j < 1+ng )   crit (i,j,k_l) = 3.0_dp
                ! if ( j > nty-ng ) crit (i,j,k_l) = 3.0_dp

             end do

          end do
       end do
    end if


  end subroutine shock_det_slv


!> \brief Ducros shock detector [Ducros et al., 1999]
!!
!! This subroutine is used in the solver.
!! The variable 'crit' is between 0 (for shock region) and 1 (for fully turbulent region)
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine shock_det_ducros_slv ( domain_id , dudx , crit )


    integer (ip) , intent (in)                                                    :: domain_id !< subdomain selection
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)                   :: dudx !< array of derivative variables
    real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (inout)  :: crit !< criteria

    integer (ip)                                      :: i , j , k
    integer (ip)                                      :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                         :: div , vort , wrk 


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    if ( ndim == 1 ) then ! 1D problem

       call abort_mpi ('Error: Can not used Ducros shock sensor in 1D simulation')

    else if ( ndim == 2 ) then ! 2D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                div = dudx (i,j,k,1) + dudx (i,j,k,4)
                div = div * div

                wrk = dudx (i,j,k,3) - dudx (i,j,k,2)

                vort = wrk * wrk

                crit (i,j,k) = div / ( div + vort + eps30 )

!                if ( crit (i,j,k) > percent_weight_SGS ) then
!                   crit (i,j,k) = 0.0_dp
!                else
!                   crit (i,j,k) = 1.0_dp
!                end if

!                crit (i,j,k) = ( 1.0_dp - crit (i,j,k) ) ** alpha

                crit (i,j,k) = 1.0_dp - crit (i,j,k)

             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                div = dudx (i,j,k,1) + dudx (i,j,k,5) + dudx (i,j,k,9)
                div = div * div

                wrk  = dudx (i,j,k,8) - dudx (i,j,k,6)
                vort = wrk * wrk

                wrk  = dudx (i,j,k,3) - dudx (i,j,k,7)
                wrk  = wrk * wrk
                vort = vort + wrk

                wrk  = dudx (i,j,k,4) - dudx (i,j,k,2)
                wrk  = wrk * wrk
                vort = vort + wrk

                crit (i,j,k) = div / ( div + vort + eps30 )

!                if ( crit (i,j,k) > percent_weight_SGS ) then
!                   crit (i,j,k) = 0.0_dp
!                else
!                   crit (i,j,k) = 1.0_dp
!                end if

!                crit (i,j,k) = ( 1.0_dp - crit (i,j,k) ) ** alpha

                crit (i,j,k) = 1.0_dp - crit (i,j,k)

             end do
          end do
       end do


    end if


  end subroutine shock_det_ducros_slv


!> \brief Ducros shock detector for post-processing
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine shock_det_ducros_post ( dx_i , dy_i , dz_i , ux , vy , wz , crit )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: crit !< criteria


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: div , vort
    real (dp) , dimension (:,:,:) , allocatable         :: dudx , dudy , dudz , &
                                                           dvdx , dvdy , dvdz , &
                                                           dwdx , dwdy , dwdz

    allocate ( dudx   (sx:ex,sy:ey,sz:ez) , &
               dudy   (sx:ex,sy:ey,sz:ez) , &
               dudz   (sx:ex,sy:ey,sz:ez) , &
               dvdx   (sx:ex,sy:ey,sz:ez) , &
               dvdy   (sx:ex,sy:ey,sz:ez) , &
               dvdz   (sx:ex,sy:ey,sz:ez) , &
               dwdx   (sx:ex,sy:ey,sz:ez) , &
               dwdy   (sx:ex,sy:ey,sz:ez) , &
               dwdz   (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate shock_det_ducros_post')


    if ( ndim == 1 ) then ! 1D problem

       call abort_mpi ('Error: Can not used Ducros shock sensor in 1D simulation')

    else if ( ndim == 2 ) then ! 2D problem


       call comm_one (ux) ; call comm_one (vy)

       call dx ( dx_i , ux , dudx )
       call dx ( dx_i , vy , dvdx )

       call dy ( dy_i , ux , dudy )
       call dy ( dy_i , vy , dvdy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                div = dudx (i,j,k) + dvdy (i,j,k)
                div = div * div

                vort = dvdx (i,j,k) - dudy (i,j,k)
                vort = vort * vort

                crit (i,j,k) = div / ( div + vort + eps30 )

!                if ( crit (i,j,k) > percent_weight_SGS ) then
!                   crit (i,j,k) = 0.0_dp
!                else
!                   crit (i,j,k) = 1.0_dp
!                end if

!                crit (i,j,k) = ( 1.0_dp - crit (i,j,k) ) ** alpha

                crit (i,j,k) = 1.0_dp - crit (i,j,k)

             end do
          end do
       end do

       deallocate ( dudx , dudy , &
                    dvdx , dvdy )


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

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                div = dudx (i,j,k) + dvdy (i,j,k) + dwdz (i,j,k)
                div = div * div

                vort = ( dwdy (i,j,k) - dvdz (i,j,k) ) * &
                       ( dwdy (i,j,k) - dvdz (i,j,k) )

                vort = vort + ( dudz (i,j,k) - dwdx (i,j,k) ) * &
                              ( dudz (i,j,k) - dwdx (i,j,k) )

                vort = vort + ( dvdx (i,j,k) - dudy (i,j,k) ) * &
                              ( dvdx (i,j,k) - dudy (i,j,k) )

                crit (i,j,k) = div / ( div + vort + eps30 )

!                if ( crit (i,j,k) > percent_weight_SGS ) then
!                   crit (i,j,k) = 0.0_dp
!                else
!                   crit (i,j,k) = 1.0_dp
!                end if

!                crit (i,j,k) = ( 1.0_dp - crit (i,j,k) ) ** alpha

                crit (i,j,k) = 1.0_dp - crit (i,j,k)

             end do
          end do
       end do

       deallocate ( dudx , dudy , dudz , &
                    dvdx , dvdy , dvdz , &
                    dwdx , dwdy , dwdz )


    end if




  end subroutine shock_det_ducros_post


!> \brief Stoichiometric mixture fraction
!!
!! Calculate the stoichiometric mixture fraction (zst) for the non-premixed combustion 
!! betweebb H2 and O2 by using the compositions at each inlets, which are defined 
!! in "input.dat" file.
!!
!! The stoichiometric reaction is given by : 
!!   H2 + 0.5 ( O2 + 3.76 N2 ) = H2O + 3.76/2 N2
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine stoichiometric_mf (inp)


    type (inp_type) , intent (inout) :: inp !< input derived type

    real (dp)   :: rst , phi


    rst = 8.0_dp ! ( YO2 / YH2 )_stoichio 

    phi = rst * inp % Y1f ( inp % index_H2 ) / inp % Y0o ( inp % index_O2 )

    ! OUTPUT : zst = stoichiometric mixture_fraction (-)
    zst = 1.0_dp / ( 1.0_dp + phi )

    if ( rank == rank_default ) then 
       write (*,*)
       write (*,*) "Stoichiometric mixture fraction = " , zst
       write (*,*)
    end if


  end subroutine stoichiometric_mf


!> \brief Calculate pdf-beta.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine pdfbeta ( x_ , a_ , b_ , pdf_ )


    use nr_subroutines , only : beta_s


    real (dp) , intent (in)                                        :: x_   !< independent variable
    real (dp) , intent (in)                                        :: a_   !< a parameter
    real (dp) , intent (in)                                        :: b_   !< b parameter
    real (dp) , intent (inout)                                     :: pdf_ !< pdf-beta


    pdf_ = beta_s ( a_ , b_ )
    pdf_ = max ( pdf_ , epsi )
    pdf_ = ( x_**(a_-1.0_dp) ) * ( (1.0_dp-x_)**(b_-1.0_dp) ) / pdf_


  end subroutine pdfbeta


!> \brief Calculate Gaussian distribution.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine gaussian ( x_ , mean , var , gauss )


    real (dp) , intent (in)                                        :: x_    !< independent variable
    real (dp) , intent (in)                                        :: mean  !< average value
    real (dp) , intent (in)                                        :: var   !< variance value
    real (dp) , intent (inout)                                     :: gauss !< Gaussian distribution


    real (dp) , parameter              :: two_pi = pi + pi
    real (dp)                          :: myvar


    myvar = max ( var , epsi )

    gauss = ( x_ - mean ) * (  x_ - mean ) / myvar
    gauss = exp ( - 0.5_dp * gauss ) / sqrt ( myvar * two_pi )

!    write (*,*) x_ , mean , gauss


  end subroutine gaussian


!> \brief Calculate log-normal distribution.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine log_normal ( x_ , mean , var , gauss )


    real (dp) , intent (in)                                        :: x_    !< independent variable
    real (dp) , intent (in)                                        :: mean  !< average value
    real (dp) , intent (in)                                        :: var   !< variance value
    real (dp) , intent (inout)                                     :: gauss !< log-normal distribution


    real (dp) , parameter              :: two_pi = pi + pi
    real (dp)                          :: myvar


    myvar = max ( var , epsi )

    gauss = ( log(x_) - mean ) * (  log(x_) - mean ) / myvar
    gauss = exp ( - 0.5_dp * gauss ) / ( x_ * sqrt ( myvar * two_pi ) )

!    write (*,*) x_ , mean , gauss


  end subroutine log_normal


end module tools
