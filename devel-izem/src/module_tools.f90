!------------------------------------------------------------------------------
! MODULE: tools
!------------------------------------------------------------------------------
!> \brief General tools and utilities.
!!
!! This module contains all the tools and utilites for the solver and
!! post-treatment.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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


 !  !> \brief Select the parameters of the simulation.
 !  !!
 !  !!
 !  !> \author
 !  !! Modeling, Simulation and Data Analysis (MSDA)\n
 !  !! Université Mohammed VI Polytechnique (UM6P)\n
 !  !! Benguerir, Morocco

 !  subroutine set_parameters


 !    integer (ip) :: ok


 !    open ( unit = unit_inp , file = file_inp , status = 'old' , iostat = ok )
 !    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (file_inp))

 !    read ( unit_inp , * ) ndim ! dimension of the problem
 !    read ( unit_inp , * ) nrv  ! reactive variables
 !    read ( unit_inp , * ) npv  ! passive variables
 !    read ( unit_inp , * ) nvv  ! variance variables

 !    close (unit_inp)

 !    nv = niv + nrv + npv + nvv

 !    !     ! number of bits for binary files depending on the compiler:
 !    !     ! i) __INTEL_COMPILER,
 !    !     ! ii) __GFORTRAN__
 !    !     nbit = 8 ! by default: works with gfortran and IBM specific compiler

 !    ! #ifdef __INTEL_COMPILER
 !    !     nbit = 2
 !    ! #endif


 ! end subroutine set_parameters

 !> \brief print code
 !!
 !!
 !> \author
 !! Modeling, Simulation and Data Analysis (MSDA)\n
 !! Université Mohammed VI Polytechnique (UM6P)\n
 !! Benguerir, Morocco

 subroutine print_code

    integer , dimension (8)                       :: valeurs

    if ( rank == rank_default ) then
      write (*,*)
      call print_date_and_time
      write (*,*) "=================================================== "
      write (*,*) "     ▄█   ▄███████▄     ▄████████    ▄▄▄▄███▄▄▄▄    "
      write (*,*) "    ███  ██▀     ▄██   ███    ███  ▄██▀▀▀███▀▀▀██▄  "
      write (*,*) "   ███▌       ▄███▀   ███    █▀   ███   ███   ███   "
      write (*,*) "   ███▌  ▀█▀▄███▀▄▄  ▄███▄▄▄      ███   ███   ███   "
      write (*,*) "   ███▌   ▄███▀   ▀ ▀▀███▀▀▀      ███   ███   ███   "
      write (*,*) "   ███  ▄███▀         ███    █▄   ███   ███   ███   "
      write (*,*) "   ███  ███▄     ▄█   ███    ███  ███   ███   ███   "
      write (*,*) "   █▀    ▀████████▀   ██████████   ▀█   ███   █▀    "
      write (*,*) "                                                    "
      write (*,*) "=================================================== "
   end if

end subroutine print_code

!> \brief print date and time
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco


subroutine print_date_and_time()

 integer , dimension (8) :: valeurs

 write (*,*) '==================================================='
 call date_and_time(values=valeurs)
 print '(a, i4, 5(a, i2.2))',' ', valeurs(1), '/', valeurs(2), '/', valeurs(3), ' ', &
 &                                valeurs(5), ':', valeurs(6), ':', valeurs(7)
 ! write (*,*) '==================================================='

end subroutine print_date_and_time   


!> \brief Create the necessary repertories for izem.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine createrep (inp)


  type (inp_type) , intent (in) :: inp !< input derived type
  logical                       :: existed


#ifdef __INTEL_COMPILER
#define _DIR_ directory
#else ! __GFORTRAN__ .or. XLF 
#define _DIR_ file    
#endif

  inquire ( _DIR_ = trim (dir_restart) , exist = existed )
  if (.not.existed) call system ( 'mkdir ' // trim (dir_restart) )

end subroutine createrep


!> \brief Normalize variales for the simulation.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine adimensionalize ( inp , adi , x , y , z , dx_i , dy_i , dz_i )


  type (inp_type) , intent (inout)                           :: inp  !< input derived type
  type (adi_type) , intent (in)                              :: adi  !< non-dimensional derived type
  real (dp) , allocatable , dimension (:) , intent (inout)   :: x    !< x-coordinate array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: y    !< y-coordinate array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: z    !< z-coordinate array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: dx_i !< inverted dx array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: dy_i !< inverted dy array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: dz_i !< inverted dz array


  integer (ip) :: l
  real (dp)    :: L_i


  L_i = 1.0_dp / adi % L_ref

  inp % initial_time = inp % initial_time / adi % time_ref
  inp % timing (:)   = inp % timing (:) / adi % time_ref

  x (:) = x (:) * L_i
  y (:) = y (:) * L_i
  z (:) = z (:) * L_i

  dx_i (:) = dx_i (:) * adi % L_ref
  dy_i (:) = dy_i (:) * adi % L_ref
  dz_i (:) = dz_i (:) * adi % L_ref

  do l = 1 , inp % nvolume
     inp % x_volmin (l) = inp % x_volmin (l) * inp % dim_length_coord * L_i
     inp % x_volmax (l) = inp % x_volmax (l) * inp % dim_length_coord * L_i
     inp % y_volmin (l) = inp % y_volmin (l) * inp % dim_length_coord * L_i
     inp % y_volmax (l) = inp % y_volmax (l) * inp % dim_length_coord * L_i
     inp % z_volmin (l) = inp % z_volmin (l) * inp % dim_length_coord * L_i
     inp % z_volmax (l) = inp % z_volmax (l) * inp % dim_length_coord * L_i
  end do

  do l = 1 , inp % nxystat
     inp % z_xystat (l) = inp % z_xystat (l) * inp % dim_length_coord * L_i
  end do

  do l = 1 , inp % nxzstat
     inp % y_xzstat (l) = inp % y_xzstat (l) * inp % dim_length_coord * L_i
  end do

  do l = 1 , inp % nyzstat
     inp % x_yzstat (l) = inp % x_yzstat (l) * inp % dim_length_coord * L_i
  end do


end subroutine adimensionalize


!> \brief Provide metrics for the derivatives.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine metrics ( x , y , z , xt , yt , zt , dx_i , dy_i , dz_i )


  real (dp) , allocatable , dimension (:) , intent (inout)   :: x    !< x-coordinate array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: y    !< y-coordinate array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: z    !< z-coordinate array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: xt   !< absolute x-coordinate array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: yt   !< absolute y-coordinate array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: zt   !< absolute z-coordinate array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: dx_i !< inverted dx array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: dy_i !< inverted dy array
  real (dp) , allocatable , dimension (:) , intent (inout)   :: dz_i !< inverted dz array

  character (len_default)                                    :: word
!  integer (ip) , parameter                                   :: order = 2*ng+1
! order not a parameter because ng is not a parameter anymore
  integer (ip)                                               :: order 
  integer (ip)                                               :: ok
  integer (ip)                                               :: l
  integer (ip)                                               :: i  , j  , k
  real (dp) , allocatable , dimension (:)                    :: dxt_i , dyt_i , dzt_i , &
  ddummy_i
  character (len_default) , parameter                        :: format_exit = '(2(1X,A,1PE18.8))'


  allocate ( xt (1-ng:ntx+ng) , dxt_i (1-ng:ntx+ng) , x (sx-ng:ex+ng) , dx_i (sx-ng:ex+ng) , &
   yt (1-ng:nty+ng) , dyt_i (1-ng:nty+ng) , y (sy-ng:ey+ng) , dy_i (sy-ng:ey+ng) , &
   zt (1-ng:ntz+ng) , dzt_i (1-ng:ntz+ng) , z (sz-ng:ez+ng) , dz_i (sz-ng:ez+ng) , &
   stat = ok )
  if ( ok > 0 ) call abort_mpi ('error allocate metrics')
 
  ! Initializing order variable
  order = 2*ng+1

  xt = 0.0_dp ; x = 0.0_dp
  yt = 0.0_dp ; y = 0.0_dp
  zt = 0.0_dp ; z = 0.0_dp


  open ( unit_grid , file = file_grid , status = 'old' , iostat = ok )

  if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (file_grid))


  ! read (unit_grid,*) ! number of processus in x-direction
  ! read (unit_grid,*) ! number of processus in y-direction
  ! read (unit_grid,*) ! number of processus in z-direction
  ! read (unit_grid,*) ! space
  ! do l = 1 , nneigh
  !    read (unit_grid,*) ! boundaries conditions
  ! end do
  ! read (unit_grid,*) ! space

  read (unit_grid,*) word
  if (word /= 'X-DIRECTION') call abort_mpi ('error: key word ''x-direction'' is missing')
  do i = 1 , ntx
     read (unit_grid,*) xt(i)
  end do

  read (unit_grid,*) word
  if (word /= 'Y-DIRECTION') call abort_mpi ('error: key word ''y-direction'' is missing')
  do j = 1 , nty
     read (unit_grid,*) yt(j)
  end do

  read (unit_grid,*) word
  if (word /= 'Z-DIRECTION') call abort_mpi ('error: key word ''z-direction'' is missing')
  do k = 1 , ntz
     read (unit_grid,*) zt(k)
  end do

  read (unit_grid,*) word
  if (word /= 'END') call abort_mpi ('error: key word ''END'' is missing')

  close (unit_grid)


  ! Checking the points
  if ( rank == rank_default ) then 
    write (*,*) '==================================================='
    write (*,'(A,1X,I8)')            " MPI processes          :" , nproc
    write (*,'(A,1X,I8)')            " # points in x-direction:" , ntx
    write (*,'(A,1X,I8)')            " # points in y-direction:" , nty
    write (*,'(A,1X,I8)')            " # points in z-direction:" , ntz
    write (*,'(A,1X,I8)')            " # procs in  x-direction:" , dims(3)
    write (*,'(A,1X,I8)')            " # procs in  y-direction:" , dims(2)
    write (*,'(A,1X,I8)')            " # procs in  z-direction:" , dims(1)
    write (*,'(A,1X,1(1X,1PE18.8))') " Grid xmin              :" , xt (1)
    write (*,'(A,1X,1(1X,1PE18.8))') " Grid xmax              :" , xt (ntx)
    write (*,'(A,1X,1(1X,1PE18.8))') " Grid ymin              :" , yt (1)
    write (*,'(A,1X,1(1X,1PE18.8))') " Grid ymax              :" , yt (nty)
    write (*,'(A,1X,1(1X,1PE18.8))') " Grid zmin              :" , zt (1)
    write (*,'(A,1X,1(1X,1PE18.8))') " Grid zmax              :" , zt (ntz)
    write (*,*) '==================================================='
 end if

 Lx = xt (ntx) - xt (1)
 Ly = yt (nty) - yt (1)
 Lz = zt (ntz) - zt (1)

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
dxt_i (1-ng:ntx+ng) = 1.0_dp
dyt_i (1-ng:nty+ng) = 1.0_dp
dzt_i (1-ng:ntz+ng) = 1.0_dp


if ( ndim >= 1 ) then


  if ( ex-sx+1 < order ) call abort_mpi ('error: not enough points in x-direction')

  allocate ( ddummy_i (1-ng:ntx+ng) , stat = ok )
  if ( ok > 0 ) call abort_mpi ('error allocate metrics 2')

  ddummy_i (1-ng:ntx+ng) = dxt_i (1-ng:ntx+ng)
  call dscalar ( 1 , ntx , ddummy_i , xt , dxt_i )

  deallocate (ddummy_i)


end if


if ( ndim >= 2 ) then


  if ( ey-sy+1 < order ) call abort_mpi ('not enough points in y-direction')

  allocate ( ddummy_i (1-ng:nty+ng) , stat = ok )
  if ( ok > 0 ) call abort_mpi ('error allocate metrics 3')

  ddummy_i (1-ng:nty+ng) = dyt_i (1-ng:nty+ng)
  call dscalar ( 1 , nty , ddummy_i , yt , dyt_i )

  deallocate (ddummy_i)


end if


if ( ndim >= 3 ) then


  if ( ez-sz+1 < order ) call abort_mpi ('not enough points in z-direction')

  allocate ( ddummy_i (1-ng:ntz+ng) , stat = ok )
  if ( ok > 0 ) call abort_mpi ('error allocate metrics 4')

  ddummy_i (1-ng:ntz+ng) = dzt_i (1-ng:ntz+ng)
  call dscalar ( 1 , ntz , ddummy_i , zt , dzt_i )

  deallocate (ddummy_i)


end if


dxt_i (1:ntx) = 1.0_dp / dxt_i (1:ntx)
dyt_i (1:nty) = 1.0_dp / dyt_i (1:nty)
dzt_i (1:ntz) = 1.0_dp / dzt_i (1:ntz)


dx_i = 1.0_dp ; dx_i (sx-ng:ex+ng) = dxt_i (sx-ng:ex+ng)
dy_i = 1.0_dp ; dy_i (sy-ng:ey+ng) = dyt_i (sy-ng:ey+ng)
dz_i = 1.0_dp ; dz_i (sz-ng:ez+ng) = dzt_i (sz-ng:ez+ng)


deallocate ( dxt_i , dyt_i , dzt_i )


end subroutine metrics


!> \brief Calculate constants used in the simulation.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine constants ( inp , adi )


  type (inp_type) , intent (inout) :: inp !< input derived type
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
  !    adi % L_ref   = 25.4e-3_dp
  adi % L_ref   = 0.00004839025546795050_dp
  Lx = Lx / adi % L_ref
  Ly = Ly / adi % L_ref
  Lz = Lz / adi % L_ref

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

  ! probes
  ! inp % probe_coord (:,:) = inp % probe_coord (:,:) / adi % L_ref

  ! ! Display adi type
  ! if ( rank == rank_default )  then
  !    call adi_display (adi)
  !    call adi_write (adi)
  ! end if


end subroutine constants


!> \brief Calculate viscosity using Sutherland's law.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

function suth (T)


  real (dp) , intent (in) :: T !< temperature
  real (dp) :: suth !< viscosity


  suth = 1.4580e-6_dp * T ** (1.5_dp)  / ( T + 110.4_dp ) ! PM fix


end function suth


!> \brief Calculate viscosity using Van's law.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

function vand (T)


  real (dp) , intent (in) :: T !< temperature
  real (dp) :: vand !< viscosity


  vand = 0.6938730e-6_dp * T


end function vand


!> \brief Estimate times for the simulation.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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



!> \brief Save THI  file.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine dump_thi_plan ( number , time , v )


  integer (ip) , intent (in)                                  :: number !< number of restart file
  real (dp) , intent (in)                                     :: time   !< time
  real (dp) , allocatable , dimension (:,:,:,:) , intent (in) :: v      !< conserved variables array


  character (len_default) , parameter                         :: format_exit = '(I8,1X,10(1X,1PE18.8))'
  character (len_default)                                     :: restart_ascii , rank_ascii
  integer (ip)                                                :: ok , i , j , k , l
  integer (kind=8)                                            :: reclmax,begin_write
  character (len_default) , parameter                         :: format_time = '(1PE18.8)'
  character(len_default)                                      :: filename

  write ( restart_ascii , format_restart ) number

  reclmax = int ( 1* nbit , kind = 8 )

  filename=trim (dir_restart) //'THI' // '_' //trim (restart_ascii)

  open ( unit_restart , file   = filename  , access    = 'direct'        , &
    recl   = reclmax   , form      = 'unformatted'   , &
    status = 'unknown' , iostat    = ok              )

  if ( ok /= 0 ) then
    write (*,*) 'error opening ' , filename
    call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
 end if

 write ( unit_restart , rec = rank+1 )   & ! it must start at 1
 (((( v (i,j,k,l) , i = sx , ex ) , &
    j = sy , ey ) , &
 k = sz , ez ) , &
 l = 1  , nv )

 close (unit_restart)



end subroutine dump_thi_plan


!> \brief Locate statistical planes or volumes.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine locate_stats ( inp , adi , x , y , z )


  type (inp_type) , intent (inout)                            :: inp !< input derived type
  type (adi_type) , intent (in)                               :: adi !< non-dimensional derived type
  real (dp) , allocatable , dimension (:) , intent (in)       :: x   !< x-coordinate array
  real (dp) , allocatable , dimension (:) , intent (in)       :: y   !< y-coordinate array
  real (dp) , allocatable , dimension (:) , intent (in)       :: z   !< z-coordinate array


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
     if ( inp % x_volmin (volume) <= x (sx) .and. &
      inp % x_volmax (volume) >= x (ex) ) then

        log_volx = .true. ! entirely included

     else if ( inp % x_volmin (volume) >= x (sx) .and. &
       inp % x_volmin (volume) <= x (ex) ) then

        log_volx = .true. ! WEST boundary

     else if ( inp % x_volmax (volume) >= x (sx) .and. &
       inp % x_volmax (volume) <= x (ex) ) then

        log_volx = .true. ! EST boundary

     else if ( ( inp % x_volmin (volume) > x (sx) .and.  &
        inp % x_volmax (volume) > x (ex) ) .or. &
     ( inp % x_volmin (volume) < x (sx) .and.  &
        inp % x_volmax (volume) < x (ex) ) ) then

        ! excluded

     else

        write (*,*) rank , 'there is a problem with volumes in x-direction'
        call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

     end if


     ! y-direction
     if ( inp % y_volmin (volume) <= y (sy) .and. &
      inp % y_volmax (volume) >= y (ey) ) then

        log_voly = .true. ! entirely included

     else if ( inp % y_volmin (volume) >= y (sy) .and. &
       inp % y_volmin (volume) <= y (ey) ) then

        log_voly = .true. ! SOUTH boundary

     else if ( inp % y_volmax (volume) >= y (sy) .and. &
       inp % y_volmax (volume) <= y (ey) ) then

        log_voly = .true. ! NORTH boundary

     else if ( ( inp % y_volmin (volume) > y (sy) .and.  &
        inp % y_volmax (volume) > y (ey) ) .or. &
     ( inp % y_volmin (volume) < y (sy) .and.  &
        inp % y_volmax (volume) < y (ey) ) ) then

        ! excluded

     else

        write (*,*) rank , 'there is a problem with volumes in y-direction'
        call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )

     end if


     ! z-direction
     if ( inp % z_volmin (volume) <= z (sz) .and. &
      inp % z_volmax (volume) >= z (ez) ) then

        log_volz = .true. ! entirely included

     else if ( inp % z_volmin (volume) >= z (sz) .and. &
       inp % z_volmin (volume) <= z (ez) ) then

        log_volz = .true. ! BEHIND boundary

     else if ( inp % z_volmax (volume) >= z (sz) .and. &
       inp % z_volmax (volume) <= z (ez) ) then

        log_volz = .true. ! FRONT boundary

     else if ( ( inp % z_volmin (volume) > z (sz) .and.  &
        inp % z_volmax (volume) > z (ez) ) .or. &
     ( inp % z_volmin (volume) < z (sz) .and.  &
        inp % z_volmax (volume) < z (ez) ) ) then

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
        write (*,'(I7,1X,I3,3(1X,1PE18.8))') rank , volume , x(sx) * L_ref , y(sy) * L_ref , z(sz) * L_ref
     end if

     if ( coords (3) == coord_vol_max (volume,3) .and. &
      coords (2) == coord_vol_max (volume,2) .and. &
      coords (1) == coord_vol_max (volume,1) ) then ! locating the last rank on the volume
        write (*,'(I7,1X,I3,3(1X,1PE18.8))') rank , volume , x(ex) * L_ref , y(ey) * L_ref , z(ez) * L_ref
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

     if ( inp % z_xystat (plane) >= z (sz)   .and. &
      inp % z_xystat (plane) <= z (ez+1) ) then ! A.Techer modif ez -> ez+1

        log_planeXY (plane) = .true.

        do k = ez+1 , sz , -1 ! A.Techer modif ez -> ez+1
           if ( z (k) >= inp % z_xystat (plane) ) inp % k_xystat (plane) = k
        end do

        if ( inp % k_xystat (plane) == sz ) inp % k_xystat (plane) = sz+1
        if ( inp % k_xystat (plane) == ez+1 .or. &
         inp % k_xystat (plane) == ez ) inp % k_xystat (plane) = ez-1

     else if ( inp % z_xystat (plane) < z (sz) .or. &
       inp % z_xystat (plane) > z (ez) ) then

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
           write ( * , format_exit ) rank , plane , 'Z' , z ( inp % k_xystat (plane) - ng-1 ) * L_ref , &
           z ( inp % k_xystat (plane)        ) * L_ref , &
           z ( inp % k_xystat (plane) + ng+1 ) * L_ref
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

     if ( inp % y_xzstat (plane) >= y (sy) .and. &
      inp % y_xzstat (plane) <= y (ey) ) then

        log_planeXZ (plane) = .true.

        do j = ey , sy , -1
           if ( y (j) >= inp % y_xzstat (plane) ) inp % j_xzstat (plane) = j
        end do

        if ( inp % j_xzstat (plane) == sy ) inp % j_xzstat (plane) = sy+1
        if ( inp % j_xzstat (plane) == ey ) inp % j_xzstat (plane) = ey-1

     else if ( inp % y_xzstat (plane) < y (sy) .or. &
       inp % y_xzstat (plane) > y (ey) ) then

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
           write ( * , format_exit ) rank , plane , 'Y' , y ( inp % j_xzstat (plane) - ng-1 ) * L_ref , &
           y ( inp % j_xzstat (plane)        ) * L_ref , &
           y ( inp % j_xzstat (plane) + ng+1 ) * L_ref

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

     if ( inp % x_yzstat (plane) >= x (sx) .and. &
      inp % x_yzstat (plane) <= x (ex) ) then

        log_planeYZ (plane) = .true.

        do i = ex , sx , -1
           if ( x (i) >= inp % x_yzstat (plane) ) inp % i_yzstat (plane) = i
        end do

        if ( inp % i_yzstat (plane) == sx ) inp % i_yzstat (plane) = sx+1
        if ( inp % i_yzstat (plane) == ex ) inp % i_yzstat (plane) = ex-1

     else if ( inp % x_yzstat (plane) < x (sx) .or. &
       inp % x_yzstat (plane) > x (ex) ) then

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
           write ( * , format_exit ) rank , plane , 'X' , x ( inp % i_yzstat (plane) - ng-1 ) * L_ref , &
           x ( inp % i_yzstat (plane)        ) * L_ref , &
           x ( inp % i_yzstat (plane) + ng+1 ) * L_ref
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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine save_stats ( number , time , dtime , dt_4 , inp , adi , v )


  integer (ip) , intent (in)                                  :: number !< number of stastical file
  real (dp) , intent (in)                                     :: time   !< time
  real (dp) , intent (in)                                     :: dtime  !< time step
  real (dp) , dimension (:) , intent (in)                     :: dt_4   !< array of different time steps
  type (inp_type) , intent (inout)                            :: inp    !< input derived type
  type (adi_type) , intent (in)                               :: adi    !< non-dimensional derived type
  real (dp) , allocatable , dimension (:,:,:,:) , intent (in) :: v      !< conserved variables array


  character (len_default) , parameter                         :: format_exit = '(I8,1X,10(1X,1PE18.8))'
  character (len_default)                                     :: stat_ascii , plane_ascii , rank_ascii
  integer (ip)                                                :: ok , volume , plane , rank_volume , rank_plane
  integer (ip) , dimension (ndimmax)                          :: dimmaxvol
  integer (kind=8)                                            :: reclmax
  integer (ip)                                                :: sx0 , ex0 , &
  sy0 , ey0 , &
  sz0 , ez0
  integer (ip)                                                :: i , j , k , l


  ! timing file


  if ( rank == rank_default ) then
   open ( unit = unit_time_stat , file = &
    & trim (dir_parent) // trim (file_time_stat) , status = 'old' , position = 'append' , iostat = ok )
   if ( ok /= 0 ) &
   & call abort_mpi ('ERROR opening ' // trim (dir_parent) // trim (file_time_stat))
   write ( unit_time_stat , format_exit ) number ,  &
   time * adi % time_ref , dtime * adi % time_ref , &
   dt_4 / dtime
   close ( unit_time_stat )
end if

write ( stat_ascii , format_restart ) number
write ( rank_ascii , format_restart ) rank


! volumes


do volume = 1 , inp % nvolume

  write ( plane_ascii , format_nplane ) volume

  if ( coords (3) == coord_vol_min (volume,3) .and. &
   coords (2) == coord_vol_min (volume,2) .and. &
   coords (1) == coord_vol_min (volume,1) ) then ! locating the first rank on the volume
     write (*,*) 'writing file ' , trim (dir_parent) // trim (dir_statvol) // &
     trim (file_statvol) // '_' // &
     trim (plane_ascii) // '_' //  &
     trim (stat_ascii)
  end if


  if ( log_volume (volume) ) then


     if (ind_files) then ! writing individual binary files per process


        reclmax = int ( (ex-sx+1) * (ey-sy+1) * (ez-sz+1) * nv * nbit , kind = 8 )

        open ( unit_statvol , file = trim (dir_statvol) //         &
         trim (file_statvol) // '_' // &
         trim (plane_ascii) // '_' //  &
         trim (stat_ascii)  // '_' //  &
         trim (rank_ascii)           , &
         access    = 'direct'        , &
         recl      = reclmax         , &
         form      = 'unformatted'   , &
         status    = 'unknown'       , &
         iostat    = ok              )

        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_parent) // &
           trim (dir_statvol) // trim (file_statvol) // '_' // &
           trim (plane_ascii) // '_' // trim (stat_ascii) // '_' // &
           trim (rank_ascii)
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        write ( unit_statvol , rec = 1 )         &
        (((( v (i,j,k,l) , i = sx , ex ) , &
          j = sy , ey ) , &
        k = sz , ez ) , &
        l = 1  , nv )

        close (unit_statvol)


     else ! writing one single binary file


        reclmax = int ( nxmax * nymax * nzmax * nv * nbit , kind = 8 )

        open ( unit_statvol + volume , file = trim (dir_statvol) //         &
         trim (file_statvol) // '_' // &
         trim (plane_ascii) // '_' //  &
         trim (stat_ascii)           , &
         access    = 'direct'        , &
         recl      = reclmax         , &
         form      = 'unformatted'   , &
         status    = 'unknown'       , &
         iostat    = ok              )

        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_parent) // trim (dir_statvol) // trim (file_statvol) // '_' // &
           trim (plane_ascii) // '_' // trim (stat_ascii)
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        dimmaxvol (3) = coord_vol_max (volume,3) - coord_vol_min (volume,3) + 1
        dimmaxvol (2) = coord_vol_max (volume,2) - coord_vol_min (volume,2) + 1
        dimmaxvol (1) = coord_vol_max (volume,1) - coord_vol_min (volume,1) + 1

        rank_volume = ( coords (3) - coord_vol_min (volume,3) ) +                                   &
        ( coords (2) - coord_vol_min (volume,2) ) * dimmaxvol (3) +                   &
        ( coords (1) - coord_vol_min (volume,1) ) * ( dimmaxvol (3) * dimmaxvol (2) )

        write ( unit_statvol + volume , rec = rank_volume + 1 ) &
        (((( v (i,j,k,l) , i = sx , ex ) ,                &
          j = sy , ey ) ,                &
        k = sz , ez ) ,                &
        l = 1  , nv )

        close ( unit_statvol + volume )


     end if


  end if


end do


! XY planes


do plane = 1 , inp % nxystat

  write ( plane_ascii , format_nplane ) plane

  if ( rank == rank_default ) write (*,*) 'writing file ' , trim (dir_parent) //         &
  trim (dir_statXY) //         &
  trim (file_statXY) // '_' // &
  trim (plane_ascii) // '_' // &
  trim (stat_ascii)

  if ( log_planeXY (plane) ) then


     if (ind_files) then ! writing individual binary files per process


        reclmax = int ( (ex-sx+1) * (ey-sy+1) * ( ng+3+ng ) * nv * nbit , kind = 8 )

        open ( unit_statXY + plane , file = trim (dir_statXY) //          &
           trim (file_statXY) // '_' //  &
           trim (plane_ascii) // '_' //  &
           trim (stat_ascii)  // '_' //  &
           trim (rank_ascii)           , &
           access    = 'direct'               , &
           recl      = reclmax                , &
           form      = 'unformatted'          , &
           status    = 'unknown'              , &
           iostat    = ok                     )

        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_parent) // &
           trim (dir_statXY) // trim (file_statXY) // '_' // &
           trim (plane_ascii) // '_' // trim (stat_ascii) // &
           trim (rank_ascii)
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        sz0 = inp % k_xystat (plane) - ng-1
        ez0 = inp % k_xystat (plane) + ng+1

        write ( unit_statXY + plane , rec = 1  )  &
        (((( v (i,j,k,l) , i = sx  , ex  ) , &
           j = sy  , ey  ) , &
        k = sz0 , ez0 ) , &
        l = 1   , nv  )

        close ( unit_statXY + plane )


     else ! writing one single binary file


        reclmax = int ( nxmax * nymax * ( ng+3+ng ) * nv * nbit , kind = 8 )

        open ( unit_statXY + plane , file = trim (dir_statXY) //          &
           trim (file_statXY) // '_' //  &
           trim (plane_ascii) // '_' //  &
           trim (stat_ascii)           , &
           access    = 'direct'               , &
           recl      = reclmax                , &
           form      = 'unformatted'          , &
           status    = 'unknown'              , &
           iostat    = ok                     )

        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_parent) // &
           trim (dir_statXY) // trim (file_statXY) // '_' // &
           trim (plane_ascii) // '_' // trim (stat_ascii)
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        sz0 = inp % k_xystat (plane) - ng-1
        ez0 = inp % k_xystat (plane) + ng+1
        rank_plane = coords (3) + coords (2) * dims (3) + 1

        write ( unit_statXY + plane , rec = rank_plane ) &
        (((( v (i,j,k,l) , i = sx  , ex  ) ,        &
           j = sy  , ey  ) ,        &
        k = sz0 , ez0 ) ,        &
        l = 1   , nv  )

        close ( unit_statXY + plane )


     end if


  end if


end do


! XZ planes


do plane = 1 , inp % nxzstat

  write ( plane_ascii , format_nplane ) plane

  if ( rank == rank_default ) write (*,*) 'writing file ' , trim (dir_parent) //         &
  trim (dir_statXZ) //         &
  trim (file_statXZ) // '_' // &
  trim (plane_ascii) // '_' // &
  trim (stat_ascii)

  if ( log_planeXZ (plane) ) then


     if (ind_files) then ! writing individual binary files per process


        reclmax = int ( (ex-sx+1) * ( ng+3+ng ) * (ez-sz+1) * nv * nbit , kind = 8 )

        open ( unit_statXZ + plane , file = trim (dir_statXZ) //          &
           trim (file_statXZ) // '_' //  &
           trim (plane_ascii) // '_' //  &
           trim (stat_ascii)  // '_' //  &
           trim (rank_ascii)           , &
           access    = 'direct'               , &
           recl      = reclmax                , &
           form      = 'unformatted'          , &
           status    = 'unknown'              , &
           iostat    = ok                     )

        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_parent) //                              &
           trim (dir_statXZ) // trim (file_statXZ) // '_' // &
           trim (plane_ascii) // '_' // trim (stat_ascii) // &
           trim (rank_ascii)
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        sy0 = inp % j_xzstat (plane) - ng-1
        ey0 = inp % j_xzstat (plane) + ng+1

        write ( unit_statXZ + plane , rec = 1 )   &
        (((( v (i,j,k,l) , i = sx  , ex  ) , &
           j = sy0 , ey0 ) , &
        k = sz  , ez  ) , &
        l = 1   , nv  )

        close ( unit_statXZ + plane )


     else ! writing one single binary file


        reclmax = int ( nxmax * ( ng+3+ng ) * nzmax * nv * nbit , kind = 8 )

        open ( unit_statXZ + plane , file = trim (dir_statXZ) //          &
           trim (file_statXZ) // '_' //  &
           trim (plane_ascii) // '_' //  &
           trim (stat_ascii)           , &
           access    = 'direct'               , &
           recl      = reclmax                , &
           form      = 'unformatted'          , &
           status    = 'unknown'              , &
           iostat    = ok                     )

        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_parent) //                              &
           trim (dir_statXZ) // trim (file_statXZ) // '_' // &
           trim (plane_ascii) // '_' // trim (stat_ascii)
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        sy0 = inp % j_xzstat (plane) - ng-1
        ey0 = inp % j_xzstat (plane) + ng+1
        rank_plane = coords (3) + coords (1) * dims (3) + 1

        write ( unit_statXZ + plane , rec = rank_plane ) &
        (((( v (i,j,k,l) , i = sx  , ex  ) ,        &
           j = sy0 , ey0 ) ,        &
        k = sz  , ez  ) ,        &
        l = 1   , nv  )

        close ( unit_statXZ + plane )


     end if


  end if


end do


! YZ planes


do plane = 1 , inp % nyzstat

  write ( plane_ascii , format_nplane ) plane

  if ( rank == rank_default ) write (*,*) 'writing file ' , trim (dir_parent) //         &
  trim (dir_statYZ) //         &
  trim (file_statYZ) // '_' // &
  trim (plane_ascii) // '_' // &
  trim (stat_ascii)

  if ( log_planeYZ (plane) ) then


     if (ind_files) then ! writing individual binary files per process


        reclmax = int ( ( ng+3+ng ) * (ey-sy+1) * (ez-sz+1) * nv * nbit , kind = 8 )

        open ( unit_statYZ + plane , file = trim (dir_statYZ) //          &
           trim (file_statYZ) // '_' //  &
           trim (plane_ascii) // '_' //  &
           trim (stat_ascii) // '_' //   &
           trim (rank_ascii)           , &
           access    = 'direct'               , &
           recl      = reclmax                , &
           form      = 'unformatted'          , &
           status    = 'unknown'              , &
           iostat    = ok                     )

        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_parent) //                              &
           trim (dir_statYZ) // trim (file_statYZ) // '_' // &
           trim (plane_ascii) // '_' // trim (stat_ascii) // &
           trim (rank_ascii)
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        sx0 = inp % i_yzstat (plane) - ng-1
        ex0 = inp % i_yzstat (plane) + ng+1

        write ( unit_statYZ + plane , rec = 1 )   &
        (((( v (i,j,k,l) , i = sx0 , ex0 ) , &
           j = sy  , ey  ) , &
        k = sz  , ez  ) , &
        l = 1   , nv  )

        close ( unit_statYZ + plane )


     else ! writing one single binary file


        reclmax = int ( ( ng+3+ng ) * nymax * nzmax * nv * nbit , kind = 8 )

        open ( unit_statYZ + plane , file = trim (dir_statYZ) //          &
           trim (file_statYZ) // '_' //  &
           trim (plane_ascii) // '_' //  &
           trim (stat_ascii)           , &
           access    = 'direct'               , &
           recl      = reclmax                , &
           form      = 'unformatted'          , &
           status    = 'unknown'              , &
           iostat    = ok                     )

        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_parent) //                              &
           trim (dir_statYZ) // trim (file_statYZ) // '_' // &
           trim (plane_ascii) // '_' // trim (stat_ascii)
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        sx0 = inp % i_yzstat (plane) - ng-1
        ex0 = inp % i_yzstat (plane) + ng+1
        rank_plane = coords (2) + coords (1) * dims (2) + 1

        write ( unit_statYZ + plane , rec = rank_plane ) &
        (((( v (i,j,k,l) , i = sx0 , ex0 ) ,        &
           j = sy  , ey  ) ,        &
        k = sz  , ez  ) ,        &
        l = 1   , nv  )

        close ( unit_statYZ + plane )


     end if


  end if


end do


end subroutine save_stats


!> \brief Generate random numbers in parallel.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine locate_probes ( inp , adi , x , y , z )


  type (inp_type) , intent (in)                                        :: inp   !< input derived type
  type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
  real (dp) , dimension (:) , allocatable , intent (in)                :: x     !< x-coordinate array
  real (dp) , dimension (:) , allocatable , intent (in)                :: y     !< y-coordinate array
  real (dp) , dimension (:) , allocatable , intent (in)                :: z     !< z-coordinate array


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

     if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
      pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
      pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true.

     if (cond) then

        ! calculate points
        do i = ex+1 , sx , -1
           if ( x (i) >= pt (1) ) px = i
        end do
        px = min(px,ex) ; px = max(px,sx)
        do j = ey+1 , sy , -1
           if ( y (j) >= pt (2) ) py = j
        end do
        py = min(py,ey) ; py = max(py,sy)
        if ( ndim == 3 ) then
           do k = ez+1 , sz , -1
              if ( z (k) >= pt (3) ) pz = k
           end do
        else
           pz = sz
        end if
        pz = min(pz,ez) ; pz = max(pz,sz)

        write ( * , format ) inp % probe_name (iprobe) , &
        x (px) * L_ref , y (py) * L_ref , z (pz) * L_ref

     end if

  end do

  call mpi_barrier ( MPI_COMM_WORLD , mpicode )

  if ( rank == rank_default ) then
     write (*,*) 'Probes ***************************************************************'
     write (*,*)
  end if


end subroutine locate_probes


!> \brief Open probes files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine open_probes ( inp , adi , x , y , z )


  type (inp_type) , intent (in)                                        :: inp !< input derived type
  type (adi_type) , intent (in)                                        :: adi !< non-dimensional derived type
  real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
  real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
  real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


  integer (ip)                        :: i , j , k , px , py , pz
  integer (ip)                        :: currunit , iprobe , ok
  real (dp)                           :: L_ref
  logical                             :: cond
  real (dp) , dimension (ndimmax)     :: pt
  character (len_default) , parameter :: format = ' ( A10 , 1X , 3 ( 1PE16.8 , 1X ) ) '


  L_ref = adi % L_ref


  do iprobe = 1 , inp % nprobes

     pt (1) = inp % probe_coord ( iprobe , 1 )
     pt (2) = inp % probe_coord ( iprobe , 2 )
     pt (3) = inp % probe_coord ( iprobe , 3 )

     cond = .false.

     if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
      pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
      pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true.

     if (cond) then

        currunit = unit_probes + iprobe - 1

        open ( unit   = currunit ,                                   &
          file   = trim (dir_probes) //                         &
          trim (inp % probe_name (iprobe)) // '.out' , &
          form   = 'formatted'    ,                             &
          status = 'unknown'      ,                             &
          iostat = ok )
        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_probes) // trim (inp % probe_name (iprobe))
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

        ! calculate points
        do i = ex+1 , sx , -1
           if ( x (i) >= pt (1) ) px = i
        end do
        px = min(px,ex) ; px = max(px,sx)
        do j = ey+1 , sy , -1
           if ( y (j) >= pt (2) ) py = j
        end do
        py = min(py,ey) ; py = max(py,sy)
        if ( ndim == 3 ) then
           do k = ez+1 , sz , -1
              if ( z (k) >= pt (3) ) pz = k
           end do
        else
           pz = sz
        end if
        pz = min(pz,ez) ; pz = max(pz,sz)

        write ( currunit , format ) 'coordinates' , x (px) * L_ref , y (py) * L_ref , z (pz) * L_ref
        write ( currunit , * )

        close ( unit = currunit ) ! A.Techer

     end if

  end do


end subroutine open_probes


!> \brief Close probes files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine close_probes ( inp , x , y , z )


  type (inp_type) , intent (in)                                        :: inp !< input derived type
  real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
  real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
  real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


  integer (ip)                       :: iprobe
  logical                            :: cond
  real (dp) , dimension (ndimmax)    :: pt


  do iprobe = 1 , inp % nprobes

     pt (1) = inp % probe_coord ( iprobe , 1 )
     pt (2) = inp % probe_coord ( iprobe , 2 )
     pt (3) = inp % probe_coord ( iprobe , 3 )

     cond = .false.

     if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
      pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
      pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true.

     if (cond) close ( unit = unit_probes + iprobe - 1 )

  end do


end subroutine close_probes


!> \brief Plot probes files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine plot_probes ( time , dtime , inp , adi , x , y , z , v )


  real (dp) , intent (in)                                              :: time  !< time
  real (dp) , intent (in)                                              :: dtime !< time step
  type (inp_type) , intent (in)                                        :: inp   !< input derived type
  type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
  real (dp) , dimension (:) , allocatable , intent (in)                :: x     !< x-coordinate array
  real (dp) , dimension (:) , allocatable , intent (in)                :: y     !< y-coordinate array
  real (dp) , dimension (:) , allocatable , intent (in)                :: z     !< z-coordinate array
  real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v     !< conserved variables array


  integer (ip)                        :: i , j , k , px , py , pz
  integer (ip)                        :: iprobe , currunit , ok
  real (dp) , dimension (ndimmax)     :: pt
  logical                             :: cond
  character (len_default) , parameter :: format = ' ( 30 ( 1PE16.8 , 1X ) ) '


  do iprobe = 1 , inp % nprobes

     pt (1) = inp % probe_coord ( iprobe , 1 )
     pt (2) = inp % probe_coord ( iprobe , 2 )
     pt (3) = inp % probe_coord ( iprobe , 3 )

     cond = .false.

     if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
      pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
      pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true.

     if (cond) then

        ! calculate points
        do i = ex+1 , sx , -1
           if ( x (i) >= pt (1) ) px = i
        end do
        px = min(px,ex) ; px = max(px,sx)
        do j = ey+1 , sy , -1
           if ( y (j) >= pt (2) ) py = j
        end do
        py = min(py,ey) ; py = max(py,sy)
        if ( ndim == 3 ) then
           do k = ez+1 , sz , -1
              if ( z (k) >= pt (3) ) pz = k
           end do
        else
           pz = sz
        end if
        pz = min(pz,ez) ; pz = max(pz,sz)

        currunit = unit_probes + iprobe - 1

        open ( unit   = currunit ,                                   &
          file   = trim (dir_probes) //                         &
          trim (inp % probe_name (iprobe)) // '.out' , &
          form   = 'formatted'    ,                             &
          status = 'old'          ,                             &
          position = 'append'     ,                             &
          iostat = ok )
        if ( ok /= 0 ) then
           write (*,*) 'error opening ' , trim (dir_probes) // trim (inp % probe_name (iprobe))
           call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
        end if

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine shock_det_slv ( v , T , W_i , crit )


  real (dp) , dimension (:,:,:,:) , allocatable , intent (in)                   :: v    !< conserved variables array
  real (dp) , dimension (:,:,:) , allocatable , intent (in)                     :: T    !< temperature
  real (dp) , dimension (:,:,:) , allocatable , intent (in)                     :: W_i  !< inverted molar mass
  real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (inout)  :: crit !< criteria


  !integer (ip) , parameter                  :: stencil_m1 = ng+ng
  ! stencil_m1 not a parameter because ng is not a parameter anymore
  integer (ip)                              :: stencil_m1
  integer (ip)                              :: st
  integer (ip)                              :: i , i_l , i_r , i_s , &
  j , j_l , j_r , j_s , &
  k , k_l , k_r , k_s
  real (dp)                                 :: drho , dpres , p_p , p_m

  ! Initializing stencil_m1
  stencil_m1 = ng+ng
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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine shock_det_ducros_slv ( domain_id , dudx , crit )


  integer (ip) , intent (in)                                                    :: domain_id !< subdomain selection
  real (dp) , dimension (:,:,:,:) , allocatable , intent (in)                   :: dudx !< array of derivative variables
  real (dp) , dimension (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , intent (inout)  :: crit !< criteria

  integer (ip)                                      :: ok , i , j , k
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
              if ( crit (i,j,k) > percent_weight_SGS ) then
                 crit (i,j,k) = 0.0_dp
              else
                 crit (i,j,k) = 1.0_dp
              end if

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
              if ( crit (i,j,k) > percent_weight_SGS ) then
                 crit (i,j,k) = 0.0_dp
              else
                 crit (i,j,k) = 1.0_dp
              end if

           end do
        end do
     end do


  end if


end subroutine shock_det_ducros_slv


!> \brief Ducros shock detector for post-processing
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
              if ( crit (i,j,k) > percent_weight_SGS ) then
                 crit (i,j,k) = 0.0_dp
              else
                 crit (i,j,k) = 1.0_dp
              end if

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
              if ( crit (i,j,k) > percent_weight_SGS ) then
                 crit (i,j,k) = 0.0_dp
              else
                 crit (i,j,k) = 1.0_dp
              end if

           end do
        end do
     end do

     deallocate ( dudx , dudy , dudz , &
       dvdx , dvdy , dvdz , &
       dwdx , dwdy , dwdz )


  end if




end subroutine shock_det_ducros_post


!> \brief Calculate the filter-width.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine filter_width ( domain_id , dx_i , dy_i , dz_i , delta )


  integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
  real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i      !< inverted dx array
  real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i      !< inverted dy array
  real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i      !< inverted dz array
  real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: delta     !< filter-width

  integer (ip)                                                   :: i , j , k
  integer (ip)                                                   :: i0 , i1 , j0 , j1 , k0 , k1
  real (dp)                                                      :: onethirds = 1.0_dp / 3.0_dp


  call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


  if      ( ndim == 1 ) then ! 1D problem

     do k = k0 , k1
        do j = j0 , j1
           do i = i0 , i1
              delta (i,j,k) = 1.0_dp / dx_i (i)
           end do
        end do
     end do

  else if ( ndim == 2 ) then ! 2D problem

     do k = k0 , k1
        do j = j0 , j1
           do i = i0 , i1
              delta (i,j,k) = 1.0_dp / sqrt ( dx_i (i) * dy_i (j) )
           end do
        end do
     end do

  else if ( ndim == 3 ) then ! 3D problem

     do k = k0 , k1
        do j = j0 , j1
           do i = i0 , i1
              delta (i,j,k) = 1.0_dp / ( dx_i (i) * dy_i (j) * dz_i (k) ) ** onethirds
           end do
        end do
     end do

  end if


end subroutine filter_width


!> \brief Stoichiometric mixture fraction
!!
!! Calculate the stoichiometric mixture fraction (zst) for the non-premixed combustion 
!! betweebb H2 and O2 by using the compositions at each inlets, which are defined 
!! in "input.dat" file.
!!
!! The stoichiometric reaction is given by : 
!!   H2 + 0.5 ( O2 + 3.76 N2 ) = H2O + 3.76/2 N2
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

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
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine gaussian ( x_ , mean , var , gauss )


  real (dp) , intent (in)                                        :: x_    !< independent variable
  real (dp) , intent (in)                                        :: mean  !< average value
  real (dp) , intent (in)                                        :: var   !< variance value
  real (dp) , intent (inout)                                     :: gauss !< Gaussian distribution


  real (dp) , parameter              :: two_pi = ( acos(-1.0_dp) ) + ( acos(-1.0_dp) )
  real (dp)                          :: myvar


  myvar = max ( var , epsi )

  gauss = ( x_ - mean ) * (  x_ - mean ) / myvar
  gauss = exp ( - 0.5_dp * gauss ) / sqrt ( myvar * two_pi )

  !    write (*,*) x_ , mean , gauss


end subroutine gaussian


!> \brief Calculate log-normal distribution.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine log_normal ( x_ , mean , var , gauss )


  real (dp) , intent (in)                                        :: x_    !< independent variable
  real (dp) , intent (in)                                        :: mean  !< average value
  real (dp) , intent (in)                                        :: var   !< variance value
  real (dp) , intent (inout)                                     :: gauss !< log-normal distribution


  real (dp) , parameter              :: two_pi = ( acos(-1.0_dp) ) + ( acos(-1.0_dp) )
  real (dp)                          :: myvar


  myvar = max ( var , epsi )

  gauss = ( log(x_) - mean ) * (  log(x_) - mean ) / myvar
  gauss = exp ( - 0.5_dp * gauss ) / ( x_ * sqrt ( myvar * two_pi ) )

  !    write (*,*) x_ , mean , gauss


end subroutine log_normal


!> \brief compute velocity tensor
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine compute_duixj ( dx_i , dy_i , dz_i , v , duixj)

  real (dp) , allocatable , dimension (:) , intent (in) :: dx_i !< inverted dx array
  real (dp) , allocatable , dimension (:) , intent (in) :: dy_i !< inverted dy array
  real (dp) , allocatable , dimension (:) , intent (in) :: dz_i !< inverted dz array
  real (dp) , allocatable , dimension (:,:,:,:) , intent (in) :: v    !< conserved variables array
  real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: duixj   !< velocity tensor

  ! local declaration
  integer (ip) :: ok , domain_id , i , j , k , l , i0 , i1 , j0 , j1 , k0 , k1
  real (dp) :: rho_i 
  real (dp) , dimension (:,:,:) , allocatable :: ux , vy , wz

  if ( ndim == 1 ) then

     allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , stat = ok )
     if ( ok > 0 ) stop 'error allocate ux in compute_duixj'

     do domain_id = -ndim , ndim
        call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
        do k = k0 , k1
           do j = j0 , j1
              do i = i0 , i1
                 rho_i      = 1.0_dp / v (i,j,k,1)
                 ux (i,j,k) = v (i,j,k,2) * rho_i
              end do
           end do
        end do
     end do

     call dx_fixed1 ( dx_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,1) )

     deallocate (ux)

     ! communicate the first derivative
     call comm_derivVel (duixj)

  end if 

  if ( ndim == 2 ) then

     allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , stat = ok )
     if ( ok > 0 ) stop 'error allocate ux in compute_duixj'
     allocate ( vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , stat = ok )
     if ( ok > 0 ) stop 'error allocate vy in compute_duixj'

     do domain_id = -ndim , ndim
        call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
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

     call dx_fixed1 ( dx_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,1) ) ! du/dx = duidj (:,:,:,1)
     call dy_fixed1 ( dy_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) ) ! du/dy = duidj (:,:,:,2)

     call dx_fixed1 ( dx_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,3) ) ! dv/dx = duidj (:,:,:,3)
     call dy_fixed1 ( dy_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,4) ) ! dv/dy = duidj (:,:,:,4)

     deallocate (ux , vy)

     ! communicate the first derivative
     call comm_derivVel (duixj)

  end if

  if ( ndim == 3 ) then

     allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , stat = ok )
     if ( ok > 0 ) stop 'error allocate ux in compute_duixj'
     allocate ( vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , stat = ok )
     if ( ok > 0 ) stop 'error allocate vy in compute_duixj'
     allocate ( wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , stat = ok )
     if ( ok > 0 ) stop 'error allocate wz in compute_duixj'

     do domain_id = -ndim , ndim
        call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
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

     call dx_fixed1 ( dx_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,1) ) ! du/dx = duidj (:,:,:,1)
     call dy_fixed1 ( dy_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,2) ) ! du/dy = duidj (:,:,:,2)
     call dz_fixed1 ( dz_i , ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,3) ) ! du/dz = duidj (:,:,:,3)

     call dx_fixed1 ( dx_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,4) ) ! dv/dx = duidj (:,:,:,4)
     call dy_fixed1 ( dy_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,5) ) ! dv/dy = duidj (:,:,:,5)
     call dz_fixed1 ( dz_i , vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,6) ) ! dv/dz = duidj (:,:,:,6)

     call dx_fixed1 ( dx_i , wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,7) ) ! dw/dx = duidj (:,:,:,7)
     call dy_fixed1 ( dy_i , wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,8) ) ! dw/dy = duidj (:,:,:,8)
     call dz_fixed1 ( dz_i , wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
        &             duixj (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,9) ) ! dw/dz = duidj (:,:,:,9)

     deallocate (ux , vy , wz)

     ! communicate the first derivative
     call comm_derivVel (duixj)

  end if

end subroutine compute_duixj


!> \brief Ducros shock detector [Ducros et al., 1999]
!!
!! This subroutine is used in the solver.
!! The variable 'crit' is between 0 (for shock region) and 1 (for fully turbulent region)
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine shock_sensor_ducros ( duixj , shock_sensor )

  real (dp) , allocatable , dimension (:,:,:,:)  , intent (in) :: duixj !< velocity tensor
  logical , allocatable , dimension (:,:,:) , intent (inout) :: shock_sensor !< shock criterion

  integer (ip) :: ok , i , j , k
  real (dp) :: div , vort , wrk , ducloc

  if ( ndim == 1 ) then ! 1D problem

     call abort_mpi ('Error: Can not used Ducros shock sensor in 1D simulation')

  else if ( ndim == 2 ) then ! 2D problem

     do k = sz , ez
        do j = sy , ey
           do i = sx , ex
              div = duixj (i,j,k,1) + duixj (i,j,k,4)
              div = div * div
              wrk = duixj (i,j,k,3) - duixj (i,j,k,2)
              vort = wrk * wrk
              ducloc = div / ( div + vort + eps30 )
              if ( ducloc > ducros_treduc ) then
                 shock_sensor (i,j,k) = .true.
              end if
           end do
        end do
     end do

  else if ( ndim == 3 ) then ! 3D problem

     do k = sz , ez
        do j = sy , ey
           do i = sx , ex
              div = duixj (i,j,k,1) + duixj (i,j,k,5) + duixj (i,j,k,9)
              div = div * div
              wrk  = duixj (i,j,k,8) - duixj (i,j,k,6)
              vort = wrk * wrk
              wrk  = duixj (i,j,k,3) - duixj (i,j,k,7)
              wrk  = wrk * wrk
              vort = vort + wrk
              wrk  = duixj (i,j,k,4) - duixj (i,j,k,2)
              wrk  = wrk * wrk
              vort = vort + wrk
              ducloc = div / ( div + vort + eps30 )
              if ( ducloc > ducros_treduc ) then
                 shock_sensor (i,j,k) = .true.
              end if
           end do
        end do
     end do

  end if

end subroutine shock_sensor_ducros

!> \brief modified Adams & Sharif indicator
!!
!! This subroutine is used in the solver.
!! The variable 'adams_shariff_sensor' is True (for shock region) 
!! and False (for fully turbulent region)
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine shock_sensor_modified_adams_shariff ( v , T , W_i , shock_sensor)

   real (dp) , allocatable , dimension (:,:,:,:) , intent (in)     :: v   !< conserved variables array
   real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
   real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
   logical , allocatable , dimension (:,:,:) , intent (inout)     :: shock_sensor !< shock criterion

   integer (ip) :: ok , st , i , j , k
   integer (ip) :: i_s , i_c , i_l , j_s , j_c , j_l , k_s , k_c , k_l 
   real (dp) :: dpres , drho
   real (dp) , dimension (:,:,:) , allocatable   :: P , rho_i


   allocate ( rho_i  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)    , &
      &       P      (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)    , &
      stat = ok )
   if ( ok > 0 ) stop 'error allocate ducros_sensor_3D'

      do k = sz , ez
         do j = sy , ey
            do i = sx , ex
               rho_i (i,j,k) = 1.0_dp / v (i,j,k,1)
               P (i,j,k) = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
            end do
         end do
      end do

   if ( ndim >= 1 ) then
      do k = sz , ez
         do j = sy , ey
            do i = sx , ex-1
               drho = abs ( v  (min(i+1,ex),j,k,1) - v (i,j,k,1) ) * &
               &      min ( rho_i (i,j,k) , rho_i (min(i+1,ex),j,k) )
               ! pressure criteria
               dpres = abs ( P (min(i+1,ex),j,k) - P (i,j,k) ) / &
               &       min ( P (i,j,k) , P (min(i+1,ex),j,k) )

               if ( drho > max_rel_weight .and. dpres > max_rel_weight ) then
                  shock_sensor (i,j,k) = .true.
               endif
            end do
         end do
      enddo
   endif

   if ( ndim >= 2 ) then
      do k = sz , ez
         do i = sx , ex
            do j = sy , ey-1
               drho = abs ( v  (i,min(j+1,ey),k,1) - v (i,j,k,1) ) * &
               &      min ( rho_i (i,j,k) , rho_i (i,min(j+1,ey),k) )
               ! pressure criteria
               dpres = abs ( P (i,min(j+1,ey),k) - P (i,j,k) ) / &
               &       min ( P (i,j,k) , P (i,min(j+1,ey),k) )

               if ( drho > max_rel_weight .and. dpres > max_rel_weight ) then
                  shock_sensor (i,j,k) = .true.
               endif
            end do
         end do
      enddo
   endif

   if ( ndim == 3 ) then
      do i = sx , ex
         do j = sy , ey
            do k = sz , ez-1
               drho = abs ( v  (i,j,min(k+1,ez),1) - v (i,j,k,1) ) * &
               &      min ( rho_i (i,j,k) , rho_i (i,j,min(k+1,ez)) )
               ! pressure criteria
               dpres = abs ( P (i,j,min(k+1,ez)) - P (i,j,k) ) / &
               &       min ( P (i,j,k) , P (i,j,min(k+1,ez)) )

               if ( drho > max_rel_weight .and. dpres > max_rel_weight ) then
                  shock_sensor (i,j,k) = .true.
               endif
            end do
         end do
      enddo
   endif

end subroutine shock_sensor_modified_adams_shariff

end module tools
