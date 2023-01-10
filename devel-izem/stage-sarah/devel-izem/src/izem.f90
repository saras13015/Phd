!------------------------------------------------------------------------------
! PROGRAM izem
!------------------------------------------------------------------------------
!> \brief Solver.
!!
!! Contains the main program.
!!
!!   Input files:\n
!!     * chem.inp       Chemical species, Chemical reactions and kinetic schemes\n
!!     * grid.dat       the computational domain with the mesh and boundary conditions\n
!!     * input.dat      input data for the calculation
!!                      (initial condition,problem, numeric, compute parameters, ...)\n
!!     * therm.dat      thermodynamic table data (cp, h, s) for all the species used\n
!!     * tran.dat       transport data (mu, diff) for all the species used\n
!!     * syntturb.dat   external input datas to create synthetic turbulence
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

program izem

  use parameters
  use type_thd
  use parsing
  use input
  use adim
  use parallel
  use eg_lib
  use weno
  use deriv
  use thdtools
  use tools
  use ICs
  use solver
  use solver_reaction
  use iohdf5
  use io

  implicit none

  type (inp_type)                                        :: inp
  type (adi_type)                                        :: adi
  type (thd_type)                                        :: thd
  type (cfg_type)                                        :: pars
  real (dp) , allocatable , dimension (:)                :: x , y , z
  real (dp) , allocatable , dimension (:)                :: xt , yt , zt
  real (dp) , allocatable , dimension (:)                :: dx_i , dy_i , dz_i
  real (dp) , allocatable , dimension (:,:,:,:)          :: v , ha , duixj
  real (dp) , allocatable , dimension (:,:,:)            :: T , W_i , cp
  real (dp) , allocatable , dimension (:,:,:)            :: mu_SGS
  logical , allocatable , dimension (:,:,:)              :: shock_sensor  
  integer (ip) , allocatable , dimension (:)             :: iweg
  real (dp)    , allocatable , dimension (:)             :: rweg
  integer (ip)                                           :: ok , ite , nb_restart , nb_stat
  real (dp)                                              :: t_wall , t_par , t_tmp1 , t_tmp2 , t_in_loop , t_out_loop
  real (dp)                                              :: time , dtmin , DTmax
  real (dp) , dimension (4)                              :: dt_4
  logical                                                :: loop , write_stat , existed
  character (len_default) , parameter                    :: format_exit = '(I8,1X,10(1X,1PE18.8))'


  ! Check if the file_stop exists before
  inquire ( file = trim (file_stop) , exist = existed )
  if (existed) then
    open ( unit = unit_stop , file = trim (file_stop) )
    close ( unit = unit_stop , status = "delete" )
  end if
  existed = .false.

  ! MPI (Message Passing Interface : for parallel computing) initialisation
  call init_mpi

  ! print code info
  call print_code

  ! init hdf5 interface
  call hdf_init()

  ! ! get the parameters for the problem
  ! call set_parameters

  ! read input data from file
  call readinp (inp , pars)

  ! topology creation
  call topology

  ! index for every domain
  call domain

  ! looking for neighboors
  call neighborhood


  ! initialisation of CHEMKIN interpreter
  if ( rank == rank_default ) then
    write (*,*) '==================================================='
  endif
  if ( rank == rank_default ) call tranfit ! create todo.thd and a "dummy" binary file
  if ( rank == rank_default ) then
    write (*,*) '==================================================='
  endif
  call mpi_barrier ( MPI_COMM_WORLD , mpicode )


  ! evaluation of the computational mesh and metrics
  call metrics ( x , y , z , xt , yt , zt , dx_i , dy_i , dz_i )

  ! create necessary repertories to save datas
  if ( rank == rank_default ) call createrep (inp)

  ! creation of derivated types for communication
  call type_conserved
  call type_derivVel
  if (vis) then
    call type_derivative
    if (LES) call type_derivative_SGS
  end if

  ! useful constants
  call constants ( inp , adi )

  ! adimensionalize the dimensional values given in the input.dat file
  call adimensionalize ( inp , adi , x , y , z , dx_i , dy_i , dz_i )

  ! definition of WENO parameters
  call wenopar ( inp % opt , inp % optord )

  ! evaluation of thd parameters
  call thd_reset (thd)
  call thd_init  ( thd , adi , file_tranfit ) ; nreac = thd % nreac
  if ( thd % nspc /= nrv .and. rank == rank_default ) STOP 'ERROR: thd % nspc /= nrv'

  ! initialisation of reactive variables
  if ( reaction .and. inp % reac_selec == dvode ) then ! CHEMKIN library
    if ( rank == rank_default ) call ckintp39 ! overwrite the previous chem.bin file from tranfit with 3.9 interpreter
    call mpi_barrier ( MPI_COMM_WORLD , mpicode )
    call initchemkin
  else if ( reaction .and. inp % reac_selec == MIL ) then ! Model Intermittent Lagrangien
    call stoichiometric_mf (inp)
    call load_table (inp)
  end if

  ! initialisation of EGlib
  if ( vis .and. eglib ) then
    if ( .not. reaction .and. rank == rank_default ) call ckintp39 ! overwrite the previous chem.bin file from tranfit with 3.9 interpreter
    if ( rank == rank_default ) call tran
    if ( rank == rank_default ) call egfrmc
    call mpi_barrier ( MPI_COMM_WORLD , mpicode )
    call init_eglib ( inp , rweg , iweg )
  end if

  ! ! prepare the stats
  ! if (inp % stat) then
  !   call locate_stats ( inp , adi , x , y , z )
  ! end if

  ! if ( inp % nprobes > 1 ) then
  !   call locate_probes ( inp , adi , x , y , z )
  ! end if

  ! ! opening probes files
  ! call open_probes ( inp , adi , x , y , z )

  ! opening the time files to store the iterations, time and dtmin
  if ( rank == rank_default ) then

    open ( unit = unit_time_rest , file = trim (dir_parent) // trim (file_time_rest) , status = 'new' , iostat = ok )
    if ( ok /= 0 ) call abort_mpi ('ERROR opening ' // trim (dir_parent) // trim (file_time_rest))
    close (unit_time_rest)

    open ( unit = unit_time_stat , file = trim (dir_parent) // trim (file_time_stat) , status = 'new' , iostat = ok )
    if ( ok /= 0 ) call abort_mpi ('ERROR opening ' // trim (dir_parent) // trim (file_time_stat))
    close (unit_time_stat)

  end if

  ! initializing some parameters
  if ( inp % read_restart ) then
    nb_stat    = inp % number_restart
    nb_restart = inp % number_restart
  else
    nb_stat    = 0
    nb_restart = 0
  end if
  time  = inp % initial_time
  dtmin = 1.0_dp
  dt_4  = 1.0_dp
  DTmax = 1.0_dp
  if ( inp % ini_sol ) then
    loop = .false.
    nb_restart = -1
  else
    loop = .true.
  end if
  write_stat = .false.

  ! initialisation of the solution
  if (inp % read_restart) then
    ! call read_restart ( inp , adi , thd , time , x , y , z , dx_i , dy_i , dz_i , T , W_i , cp , ha , v )
    call read_restart_hdf5 ( inp , adi , thd , time , x , y , z , dx_i , dy_i , dz_i , T , W_i , cp , ha , v )

  else
    if ( rank == rank_default ) then
      write (*,*) 'initializing the problem: ', trim (inp % init)
      write (*,*) '==================================================='
    endif
    call init_selector ( inp , adi , thd , x , y , z , T , W_i , cp , ha , v )
  end if


  if ( rank == rank_default .and. loop ) write (*,*) 'entering the main loop ...'
  if ( rank == rank_default )  write (*,*) '==================================================='

  ! if ( rank == rank_default .and. loop ) write (*,'(A8,1X,10(1X,A18))') &
  ! 'ite' , 'time' , 'dtmin' ,     &
  ! 'dt_CFL*' , 'dt_mass_diff*' , 'dt_therm_diff*' , 'dt_chem*'


  ! taking the times
  t_in_loop = MPI_WTIME()
  t_wall    = t_in_loop
  t_tmp1    = t_in_loop
  t_tmp2    = time
  tcomm     = 0.0_dp
  t_par     = tcomm

  allocate ( mu_SGS (sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng) , stat = ok )
  if ( ok > 0 ) call abort_mpi ('error allocate mu_SGS')
  allocate ( shock_sensor (sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng) , stat = ok )
  if ( ok > 0 ) call abort_mpi ('error allocate shock_sensor')
  allocate ( duixj (sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndim * ndim) , stat = ok )
  if ( ok > 0 ) call abort_mpi ('error allocate duixj')

  mu_SGS = 0.0_dp
  duixj = 0.0_dp
  shock_sensor = .false.

  ! main loop
  ite = 0
  do while (loop)

    ite = ite + 1
    if ( inp % itlim .and. ite >= inp % itmax ) then
      loop = .false.
      if ( rank == rank_default ) write (*,*) 'number of iterations has been reached'
    end if

    ! random number evaluation every 1000 time steps
    if ( ite == 1 .or. mod (ite,1000) == 0 ) call randomize (inp)

    ! calculate time step
    call timestep ( inp , adi , thd , dx_i , dy_i , dz_i , cp , W_i , T , v , mu_SGS , DTmax , dt_4 , dtmin )
    ! modify dtmin to adjust saving times
    if ( time + dtmin > inp % timing (nb_stat+1) .and. time - dtmin < inp % timing (nb_stat+1) ) then
      dtmin = inp % timing (nb_stat+1) - time
      if (inp % stat) write_stat = .true.
      if ( nb_stat+1 >= inp % nstat ) then
        loop = .false.
        if ( rank == rank_default ) write (*,*) 'number of stat files has been reached'
      end if
    end if

    if (reaction) then ! Strang splitting

      call reaction_selec ( inp , thd , adi , time , 0.5_dp * dtmin , x , y , z , &
        dx_i , dy_i , dz_i , mu_SGS , T , W_i , cp , ha , v , DTmax )
      call rk3 ( inp , thd , adi , time , dtmin , x , y , z , dx_i , dy_i , dz_i , mu_SGS , &
        rweg , iweg , T , W_i , cp , ha , shock_sensor , duixj , ite , v )
      time = time + dtmin
      call reaction_selec ( inp , thd , adi , time , 0.5_dp * dtmin , x , y , z , &
        dx_i , dy_i , dz_i , mu_SGS , T , W_i , cp , ha , v , DTmax )

    else

      call rk3 ( inp , thd , adi , time , dtmin , x , y , z , dx_i , dy_i , dz_i , mu_SGS , &
        rweg , iweg , T , W_i , cp , ha , shock_sensor , duixj , ite , v )
      time = time + dtmin

    end if

    ! storing stats when conditions satisfy
    if (write_stat) then

      if ( rank == rank_default) then
        write (*,'(A,1X,I8)')            " iteration              :" , ite
        write (*,'(A,1X,1(1X,1PE18.8))') " time                   :" , time * adi % time_ref
        write (*,'(A,1X,1(1X,1PE18.8))') " Minimum dt             :" , dtmin * adi % time_ref
        write (*,'(A,1X,1(1X,1PE18.8))') " Convective dt          :" , dt_4(1) / dtmin
        write (*,'(A,1X,1(1X,1PE18.8))') " Mass diffusivity dt    :" , dt_4(2) / dtmin
        write (*,'(A,1X,1(1X,1PE18.8))') " Thermal diffusivity dt :" , dt_4(3) / dtmin
        write (*,'(A,1X,1(1X,1PE18.8))') " Chemical criteria dt   :" , dt_4(4) / dtmin
        write (*,*)
        write (*,*) '==================================================='
      end if


      ! if ( rank == rank_default ) write ( * , format_exit ) ite , &
      ! &                           time * adi % time_ref , dtmin * adi % time_ref , dt_4 / dtmin

      ! ! write probes files
      ! call plot_probes ( time , dtmin , inp , adi , x , y , z , v )

      nb_stat = nb_stat + 1
      if ( rank == rank_default ) then
        open ( unit = unit_time_stat , file = &
          & trim (dir_parent) // trim (file_time_stat) , status = 'old' , position = 'append' , iostat = ok )
        if ( ok /= 0 ) call abort_mpi ('ERROR opening ' // trim (dir_parent) // trim (file_time_stat))
        write ( unit_time_stat , format_exit ) nb_stat ,  &
        time * adi % time_ref , dtmin * adi % time_ref , dt_4 / dtmin
        close ( unit_time_stat )
      end if

      ! call save_stats ( nb_stat , time , dtmin , dt_4 , inp , adi , v )
      nb_restart = nb_restart + 1
      call save_restart_hdf5 ( adi , thd , nb_stat , time , v )
      write_stat = .false.

      if ( rank == rank_default) then
        write (*,'(A,1X,I8)')            " iteration              :" , ite
        write (*,'(A,1X,1(1X,1PE18.8))') " time                   :" , time * adi % time_ref
        write (*,'(A,1X,1(1X,1PE18.8))') " Minimum dt             :" , dtmin * adi % time_ref
        write (*,'(A,1X,1(1X,1PE18.8))') " Convective dt          :" , dt_4(1) / dtmin
        write (*,'(A,1X,1(1X,1PE18.8))') " Mass diffusivity dt    :" , dt_4(2) / dtmin
        write (*,'(A,1X,1(1X,1PE18.8))') " Thermal diffusivity dt :" , dt_4(3) / dtmin
        write (*,'(A,1X,1(1X,1PE18.8))') " Chemical criteria dt   :" , dt_4(4) / dtmin
        write (*,*)
        write (*,*) '==================================================='
      end if



      ! if ( rank == rank_default ) write (*,'(A8,1X,10(1X,A18))') 'ite' , 'time' , 'dtmin' ,     &
      ! 'dt_CFL*' , 'dt_mass_diff*' , 'dt_therm_diff*' , 'dt_chem*'

    end if

    ! evaluate times every iteration including savings, etc.
    call mpi_barrier ( MPI_COMM_WORLD , mpicode ) ! to synchronize time for each rank
    t_wall = MPI_WTIME()
    t_par  = t_par + tcomm

    t_tmp1 = t_wall - t_tmp1
    t_tmp2 = time - t_tmp2
    ! show information at the screen
    if ( mod (ite,inp % itshow) == 0 .and. rank == rank_default ) then
      t_tmp1 = t_wall - t_tmp1
      t_tmp2 = time - t_tmp2


      if ( rank == rank_default) then
        write (*,'(A,1X,I8)')            " iteration              :" , ite
        write (*,'(A,1X,1(1X,1PE18.8))') " time                   :" , time * adi % time_ref
        write (*,'(A,1X,1(1X,1PE18.8))') " Minimum dt             :" , dtmin * adi % time_ref
        write (*,'(A,1X,1(1X,1PE18.8))') " Convective dt          :" , dt_4(1) / dtmin
        write (*,'(A,1X,1(1X,1PE18.8))') " Mass diffusivity dt    :" , dt_4(2) / dtmin
        write (*,'(A,1X,1(1X,1PE18.8))') " Thermal diffusivity dt :" , dt_4(3) / dtmin
        write (*,'(A,1X,1(1X,1PE18.8))') " Chemical criteria dt   :" , dt_4(4) / dtmin
        write (*,*)
        write (*,*) '==================================================='
      end if



      ! if ( rank == rank_default ) write ( * , format_exit ) ite , time * adi % time_ref , dtmin * adi % time_ref , &
      ! dt_4 / dtmin
      ! call counter ( time / inp % timing (inp % nstat) , t_tmp2 / inp % timing (inp % nstat) , &
      !   &  t_tmp1 , tcomm , t_wall-t_in_loop )
      t_tmp1 = t_wall
      t_tmp2 = time
      tcomm  = 0.0_dp
    end if

    ! cheking the wall time and STOP file every 20 iterations
    if ( ( t_wall - t_in_loop ) > 0.95_dp * (inp % walltime) ) then
      loop = .false.
      if ( rank == rank_default ) write (*,*) 'walltime has been reached'
      if ( rank == rank_default )  write (*,*) '==================================================='
    end if
    if ( mod (ite,20) == 0 ) then
      inquire ( file = trim ( file_stop ) , exist = existed )
      if (existed) then
        loop = .false.
        if ( rank == rank_default ) write (*,*) 'request to stop the program'
      end if
    end if

  end do ! main loop



  if ( rank == rank_default .and. .not. loop ) write (*,*) 'exit the main loop ...'
  if ( rank == rank_default )  write (*,*) '==================================================='

  ! measure the loop time
  t_out_loop = MPI_WTIME()

  nb_restart = nb_restart + 1
  call save_restart_hdf5 ( adi , thd , nb_restart , time , v )

  ! this is necessary to avoid closing temporal files
  call mpi_barrier ( MPI_COMM_WORLD , mpicode )

  if ( rank == rank_default ) then

    t_out_loop = t_out_loop - t_in_loop
    t_in_loop  = 100.0_dp * ( 1.0_dp - t_par / t_out_loop )

    ! write (*,*)
    write (*,*) '==================================================='
    write (*,'(A,ES12.3E1,A)') ' program finished in         = ' , &
    &                             t_out_loop / 60.0_dp , ' minutes'
    write (*,'(A,ES12.3E1,A)') ' communications time was     = ' , &
    &                             t_par / 60.0_dp , ' minutes'
    write (*,'(A,ES12.3E1,A)') ' average parallel efficiency = ' , &
    &                             t_in_loop , ' %'
    write (*,*) '==================================================='
    call print_code
    write (*,*) 'The program is finished, you can go home           '
  end if

  deallocate ( v , ha , T , W_i , cp , mu_SGS , duixj , shock_sensor) 

  ! detroy hdf5
  call hdf_destroy

  ! MPI exit
  call end_mpi

end program izem

