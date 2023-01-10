!------------------------------------------------------------------------------
! PROGRAM postizem
!------------------------------------------------------------------------------
!> \brief Post-treatment.
!!
!! Contains the main program.
!!
!!   Input files:\n
!!     * chem.inp    chemical species, chemical reactions and kinetic schemes\n
!!     * grid.dat    the computational domain with the mesh and boundary conditions\n
!!     * input.dat   input data for the calculation
!!                   (problem, numeric, compute parameters, ...)\n
!!     * post.dat    input data for the post-treatment of
!!                   instantaneous fields\n
!!     * therm.dat   thermodynamic table data (cp, h, s) for all the species used\n
!!     * tran.dat    transport data (mu, diff) for all the species used
!!     * tabchem.dat tabulated chemical time scale created by the solver IZEM for reaction MIL
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

program postizem


  use parameters
  use parallel
  use type_thd
  use parsing
  use input
  use adim
  use eg_lib
  use deriv
  use thdtools
  use tools
  use tools_post
  use variables
  use statistics
  use paraview
  use iohdf5
  use io_post


  implicit none

  type (inp_type)                                        :: inp
  type (adi_type)                                        :: adi
  type (thd_type)                                        :: thd
  type (cfg_type)                                        :: pars
  type (sim_type)                                        :: sim
  type (inst_type)                                       :: inst
  type (reyavg_type)                                     :: reyavg
  type (favavg_type)                                     :: favavg
  type (reyfluc_type)                                    :: reyfluc
  type (favfluc_type)                                    :: favfluc
  type (stat_type)                                       :: stat

  real (dp) , allocatable , dimension (:)                :: x  , y  , z
  real (dp) , allocatable , dimension (:)                :: xt , yt , zt

  real (dp) , allocatable , dimension (:)                :: dx_i , dy_i , dz_i

  real (dp) , allocatable , dimension (:)                :: t , dt

  ! EGlib
  integer (ip) , allocatable , dimension (:)             :: iweg
  real (dp)    , allocatable , dimension (:)             :: rweg

  integer (ip)                                           :: ifile , i , j , k , ok , nb_stats

  ! time measurement
  integer   , dimension (8)                              :: valeurs
  real (dp)                                              :: t_par , t_in_loop , t_out_loop

  logical                                                :: existed

  ! MPI (Message Passing Interface : for parallel computing) initialisation
  call init_mpi

  ! print code info
  call print_code



#ifdef __INTEL_COMPILER
#define _DIR_ directory
#else ! __GFORTRAN__ .or. XLF 
#define _DIR_ file    
#endif

  ! create necessary repertories to plot
  if ( rank == rank_default ) then
     inquire ( _DIR_ = trim (dir_plot) , exist = existed )
     if (.not.existed) call system ( 'mkdir ' // trim (dir_plot) )
  end if

  ! init hdf5 interface
  call hdf_init()

  ! read input data from file
  call readinp (inp , pars)

  ! ! get the parameters for the problem
  ! call set_parameters

  ! topology creation
  call topology
  call degenerate_topology

  ! index for every domain
  call domain

  ! looking for neighboors
  call neighborhood


  ! initialisation of CHEMKIN interpreter
  if ( rank == rank_default ) call tranfit ! create todo.thd and a "dummy" binary file
  call mpi_barrier ( MPI_COMM_WORLD , mpicode )

  
  ! read post data from file
  call readpost (inp)

  ! evaluation of the computational mesh and metrics
  call metrics ( x , y , z , xt , yt , zt , dx_i , dy_i , dz_i )

  
 ! creation of derivated types for communication
  call type_conserved
  call type_one
  call type_edges
  if (vis) call type_derivative


  ! useful constants
  call constants ( inp , adi )

  ! adimensionalize the dimensional values given in the input.dat file
  call adimensionalize_post ( inp , adi , x , y , z , xt , yt , zt , dx_i , dy_i , dz_i )

  ! evaluation of thd parameters
  call thd_reset (thd)
  call thd_init  ( thd , adi , file_tranfit )
  if ( thd % nspc /= inp % nrv ) call abort_mpi ('error: thd % nspc /= inp % nrv')
  if ( thd % nreac /= inp % nreac ) call abort_mpi ('error: thd % nreac /= inp % nreac')

  ! initialisation of reactive variables
  if ( reaction .and. ( inp % reac_selec == dvode ) .and. .not. inp % read_stat ) then ! CHEMKIN library
     if ( rank == rank_default ) call ckintp39 ! overwrite the previous chem.bin file from tranfit with 3.9 interpreter
     call mpi_barrier ( MPI_COMM_WORLD , mpicode )
     call initchemkin
  else if ( reaction .and. ( inp % reac_selec == MIL ) .and. .not. inp % read_stat ) then ! Model Intermittent Lagrangien
     call stoichiometric_mf (inp)
     call load_table (inp)
  end if

  ! initialisation of EGlib
  if ( vis .and. eglib ) then
     if ( .not. reaction .and. rank == rank_default ) call ckintp39 ! overwrite the previous chem.bin file from tranfit with 3.9 interpreter
     if ( rank == rank_default ) call tran   ! create tran.bin file
     if ( rank == rank_default ) call egfrmc ! convert tran.bin into Linkeg
     call mpi_barrier ( MPI_COMM_WORLD , mpicode )
     call init_eglib ( inp , rweg , iweg )
  end if

  ! transform restart fields (sequential algorithm: not parallel)
  if ( inp % trafo ) then
     t_in_loop = MPI_WTIME()
     do ifile = inp % start_file , inp % end_file , inp % skip_file
        call transform ( ifile , inp , thd )
        t_out_loop = MPI_WTIME()
     end do
     write (*,*) 'program finished in ' , (t_out_loop - t_in_loop) / 60.0_dp , ' minutes'
     call end_mpi
     call abort_mpi (' ')
  end if

  ! taking the times
  t_in_loop = MPI_WTIME()
  tcomm     = 0.0_dp
  t_par     = tcomm

  ! reading statistical times
  call readtimes ( inp , adi , sim , t , dt )


  ! main loop
  if ( .not. inp % temp_avg .and. .not. inp % spat_avg ) then ! instantaneous plots

     ! allocate variables
     call alloc_inst_type (inst)

     if ( inp % read_stat ) then
        if ( rank == rank_default ) write (*,*) 'reading previous statistics ...'
        call read_inst_type ( inp , inst )
        call vtr_paraview ( sim % correct_endfile + 1 , inp , thd , adi , sim , &
                            x , y , z , xt , yt , zt , dx_i , dy_i , dz_i ,     &
                            inst , reyavg , favavg , stat )
        call pvtr_paraview ( sim % correct_endfile + 1 , inp )
     else
        if ( rank == rank_default ) write (*,*) 'reading files ', trim( inp % plane_type ) , '...'
        do ifile = inp % start_file , inp % end_file , inp % skip_file

           if ( rank == rank_default .and. mod (ifile,inp % itshow) == 0 ) write (*,*) 'reading file' , ifile

           t_out_loop = MPI_WTIME()
           t_par = t_par + tcomm

           ! call read_stat ( ifile , inp , thd , adi , dx_i , dy_i , dz_i , rweg , iweg , inst )
           call read_restart_hdf5 ( ifile , inp , thd , adi , dx_i , dy_i , dz_i , rweg , iweg , inst )

          ! call data_turbulent ( ifile , t(ifile) ,  adi , thd , inst , dx_i , dy_i , dz_i )
            nb_stats = inp % end_file
!            call layers_turbulent ( inp % end_file , ifile ,t(ifile) , adi , thd , inst , xt , yt , zt , dx_i , dy_i , dz_i , x , y , z )

            call vtr_paraview ( ifile , inp , thd , adi , sim ,                &
                               x , y , z , xt , yt , zt , dx_i , dy_i , dz_i , &
                               inst , reyavg , favavg , stat )
            call pvtr_paraview ( ifile , inp )

        end do
        !call write_inst_type (inst)
     end if

     if ( rank == rank_default ) write (*,*) 'writing paraview files ...'
     call pvd_paraview ( inp , adi , sim , t )


  else ! statistics


     ! allocate variables
     call alloc_inst_type (inst)
     call alloc_reyavg_type (reyavg)
     call alloc_favavg_type (favavg)
     call alloc_reyfluc_type (reyfluc)
     call alloc_favfluc_type (favfluc)
     call alloc_stat_type (stat)


     ! open ASCII files
     call open_convergence ( inp , adi , x , y , z )
     if ( inp % pdfs )     call open_pdfs ( inp , adi , x , y , z )
     if ( inp % spectras ) call open_spectras ( inp , adi , x , y , z )
     if ( inp % cond_avg ) call open_cond_avgs ( inp , adi , x , y , z )


     ! reading previous stats
     if ( inp % read_stat ) then


        if ( rank == rank_default ) write (*,*) 'reading previous statistics ...'
        call read_reyavg_type ( inp , reyavg )
        call read_favavg_type ( inp , favavg )
        call read_stat_type ( inp , stat )


     else



        ! first reading: calculate only the reynolds and favre averages
        if ( rank == rank_default ) write (*,*) 'entering first loop ...'
        do ifile = inp % start_file , inp % end_file , inp % skip_file

           t_out_loop = MPI_WTIME()
           t_par = t_par + tcomm

           if ( mod (ifile,inp % itshow) == 0 ) then
              if ( rank == rank_default ) write (*,*) 'first loop, file ' , ifile , &
                 ( t_out_loop - t_in_loop ) / 60.0_dp , ' minutes running,' ,   &
                 100.0_dp * ( 1.0_dp - tcomm / ( t_out_loop - t_in_loop ) ) ,       &
                 '% parallel efficiency'
              tcomm = 0.0_dp
           end if

           ! call read_stat ( ifile , inp , thd , adi , dx_i , dy_i , dz_i , rweg , iweg , inst )
           call read_restart_hdf5 ( ifile , inp , thd , adi , dx_i , dy_i , dz_i , rweg , iweg , inst )
           call averages ( ifile , inp , thd , adi , sim , dx_i , dy_i , dz_i , dt , &
                           inst , reyavg , favavg )
           if ( inp % pdfs )     call pdfs ( ifile , inp , thd , adi , sim , dx_i , dy_i , dz_i , dt , inst , favavg , stat )
           if ( inp % spectras ) call spectras ( ifile , inp , thd , sim , x , y , z , dx_i , dy_i , dz_i , inst , stat )

        end do


        ! second reading: calculate fluctuations and statistical variables from averages
        if ( inp % second_loop ) then
           if ( rank == rank_default ) write (*,*) 'entering second loop ...'
           do ifile = inp % start_file , inp % end_file , inp % skip_file

              t_out_loop = MPI_WTIME()
              t_par = t_par + tcomm

              if ( mod (ifile,inp % itshow) == 0 ) then
                 if ( rank == rank_default ) write (*,*) 'second loop, file ' , ifile , &
                    ( t_out_loop - t_in_loop ) / 60.0_dp , ' minutes running,' ,    &
                    100.0_dp * ( 1.0_dp - tcomm / ( t_out_loop - t_in_loop ) ) ,        &
                    '% parallel efficiency'
                 tcomm = 0.0_dp
              end if

              ! call read_stat ( ifile , inp , thd , adi , dx_i , dy_i , dz_i , rweg , iweg , inst )
               call read_restart_hdf5 ( ifile , inp , thd , adi , dx_i , dy_i , dz_i , rweg , iweg , inst )
              call fluctuations ( inp , thd , adi , dx_i , dy_i , dz_i ,       &
                                  inst , reyavg , favavg , reyfluc , favfluc )
              call stats ( ifile , inp , thd , adi , sim , xt , yt , zt , dx_i , dy_i , dz_i , dt , &
                           inst , reyavg , favavg , reyfluc , favfluc , stat )
              call plot_convergence ( ifile , inp , adi , sim , x , y , z , t , dt , stat )

           end do
        end if


        ! writing binary files to recover the stats
        if ( rank == rank_default ) write (*,*) 'writing statistical files ...'
!        call write_reyavg_type (reyavg)
!        call write_favavg_type (favavg)
!        call write_stat_type (stat)

     end if


     ! writing paraview files
     if ( rank == rank_default ) write (*,*) 'writing paraview files ...'
     call vtr_paraview ( sim % correct_endfile , inp , thd , adi , sim , &
                         x , y , z , xt , yt , zt , dx_i , dy_i , dz_i , &
                         inst , reyavg , favavg , stat )
     call pvtr_paraview ( sim % correct_endfile , inp )
     call pvd_paraview ( inp , adi , sim , t )


     ! writing ASCII files
     if ( ( inp % pdfs .or. inp % spectras .or. inp % cond_avg ) .and. rank == rank_default ) &
        write (*,*) 'writing ASCII files ...'
     if ( inp % pdfs )     call plot_pdfs ( inp , adi , x , y , z , favavg , stat )
     if ( inp % spectras ) call plot_spectras ( inp , adi , sim , x , y , z , t , stat )
     if ( inp % cond_avg ) call plot_cond_avgs ( inp , adi , x , y , z , favavg , stat )


     ! close ASCII files
     if ( .not. inp % read_stat ) then
        call close_convergence ( inp , x , y , z )
        if ( inp % pdfs )     call close_pdfs ( inp , x , y , z )
        if ( inp % spectras ) call close_spectras ( inp , x , y , z )
        if ( inp % cond_avg ) call close_cond_avgs ( inp , x , y , z )
     end if


  end if


  ! this is necessary to avoid closing temporal files
  call mpi_barrier ( MPI_COMM_WORLD , mpicode )

  ! total time
  t_out_loop = MPI_WTIME()

  if ( rank == rank_default ) then
     t_out_loop = t_out_loop - t_in_loop
     t_in_loop  = 100.0_dp * ( 1.0_dp - t_par / t_out_loop )
     write (*,*) 'program finished in ' , t_out_loop / 60.0_dp , ' minutes'
     write (*,*) 'communications time was ' , t_par / 60.0_dp , ' minutes'
     write (*,*) 'average parallel efficiency of ' , t_in_loop , '%'
  end if

  ! detroy hdf5
  call hdf_destroy

  ! MPI exit
  call end_mpi

end program postizem
