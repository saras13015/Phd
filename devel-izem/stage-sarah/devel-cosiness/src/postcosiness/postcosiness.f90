!------------------------------------------------------------------------------
! PROGRAM POSTCOSINESS
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
!!     * tabchem.dat tabulated chemical time scale created by the solver for reaction MIL
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
program postcosiness


  use parameters
  use parallel
  use type_thd
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

  implicit none

  type (inp_type)                                        :: inp
  type (adi_type)                                        :: adi
  type (thd_type)                                        :: thd
  type (sim_type)                                        :: sim
  type (inst_type)                                       :: inst
  type (reyavg_type)                                     :: reyavg
  type (favavg_type)                                     :: favavg
  type (reyfluc_type)                                    :: reyfluc
  type (favfluc_type)                                    :: favfluc
  type (stat_type)                                       :: stat
  type (inp_grid)                                        :: grid


  real (dp) , allocatable , dimension (:)                :: t , dt

  ! EGlib
  integer (ip) , allocatable , dimension (:)             :: iweg
  real (dp)    , allocatable , dimension (:)             :: rweg

  integer (ip)                                           :: ifile

  ! time measurement
  integer   , dimension (8)                              :: valeurs
  real (dp)                                              :: t_par , t_in_loop , t_out_loop

  real (dp)                                              :: wrk

  character (len_default)                                :: adr_file
  character (len_default) , parameter                    :: format_header = '(A30,A14,A16)' , &
                                                            format_exit   = '(A30,F14.3,F16.3)'

  logical                                                :: verbose

  ! MPI (Message Passing Interface : for parallel computing) initialisation
  call init_mpi

  ! POSTCOSINESS version and date
  if ( rank == rank_default ) then
     write (*,*)
     write (*,*) '============================'
     write (*,*) 'POSTCOSINESS'
     call date_and_time ( values = valeurs )
     print " (a ,2( i2 .2 , a ) , i4 ,a ,3( i2 .2 , a ) , a ) " , &
              " " , valeurs (3) , "/" , valeurs (2) , "/" , valeurs (1) , & ! date
              " " , valeurs (5) , ":" , valeurs (6) , ":" , valeurs (7)     ! time
     write (*,*) '============================'
     write (*,*)
  end if


  ! useful constants
  call constants ( adi )

  ! initialisation of CHEMKIN interpreter
  if ( rank == rank_default ) call tranfit ! create todo.thd and a "dummy" binary file
  call mpi_barrier ( MPI_COMM_WORLD , mpicode )

  ! evaluation of thd parameters
  call thd_reset (thd)
  call thd_init  ( thd , adi , file_tranfit ) ; nreac = thd % nreac

  ! ! index of the main species (if H2 or O2 are not respectively the fuel and oxydizer, you must change them)
  ! call thd_species_index ( thd , 'H2'  , inp % index_H2  ) ! Fuel
  ! call thd_species_index ( thd , 'O2'  , inp % index_O2  ) ! Oxydizer
  ! call thd_species_index ( thd , 'H2O' , inp % index_H2O ) ! Product
  ! call thd_species_index ( thd , 'N2'  , inp % index_N2  ) ! Diluent
  ! inp % index_fuel = inp % index_H2
  ! inp % index_oxidizer = inp % index_O2

  ! read grid.dat
  call readgrid ( adi , grid )

  ! read files input.dat AND post.dat
  call readinp ( adi , thd , inp )
  call readpost ( adi , inp )

  ! reading statistical times
  sim % nfiles = ( inp % end_file - inp % start_file ) / inp % skip_file + 1
!  if ( .not. inp % trafo ) call readtimes ( inp , adi , sim , t , dt )
  call readtimes ( inp , adi , sim , t , dt )


  ! transform restart fields (can be placed here if don't need topology !)
  if ( inp % trafo ) then

     if ( rank == rank_default ) write (*,*) 'transforming files...'
     t_in_loop = MPI_WTIME()

     if ( nproc > sim % nfiles .and. rank == rank_default ) then
        write (*,*) trim(str(nproc)) , " nproc > " , trim(str(sim % nfiles)) , " sim % nfiles" , &
                    " ==> ranks " , trim(str(sim % nfiles)) , " to " , trim(str(nproc)) , " are useless"
        write (*,*)
     end if

     verbose = .true.
     if ( rank <= sim % nfiles - 1 ) then
        wrk = sim % nfiles / float( min( nproc , sim % nfiles ) )
        sim % sfile = inp % start_file + floor (  float(rank)            * wrk          ) * inp % skip_file
        sim % efile = inp % start_file + floor ( (float(rank) + 1.0_dp ) * wrk - 1.0_dp ) * inp % skip_file
        do ifile = sim % sfile , sim % efile , inp % skip_file
           call transform ( ifile , inp , thd , grid , verbose )
           verbose = .false.
        end do
     end if

     ! this is necessary to avoid closing temporal files
     call mpi_barrier ( MPI_COMM_WORLD , mpicode )
     ! total time
     t_out_loop = MPI_WTIME()
     if ( rank == rank_default ) then
        write (*,*) 'program finished in ' , (t_out_loop - t_in_loop) / 60.0_dp , ' minutes'
        write (*,*)
     end if
     stop

   end if


  ! create MPI topology
  call topology
  call degenerate_topology
  ! index for every subdomain
  call subdomain
  ! looking for neighboors
  call neighborhood
  ! evaluation of the computational mesh and metrics
  call metrics ( grid % xt   , grid % yt   , grid % zt   , &
                 grid % x    , grid % y    , grid % z    , &
                 grid % dx_i , grid % dy_i , grid % dz_i , &
                 grid % delta                            )

  ! create necessary repertories
  if ( rank == rank_default ) call createrep_post (inp)

  ! creation of derivated types for communication
  call type_conserved
  call type_one
  call type_edges
  if (vis) call type_derivative

  ! initialisation of reactive variables
  if ( reaction ) then 
     if (( inp % reac_selec == PSR .or. inp % reac_selec == hybrid_PSR_MIL ) ) then ! CHEMKIN library
        if ( rank == rank_default ) call ckintp39 ! overwrite the previous chem.bin file from tranfit with 3.9 interpreter
        call mpi_barrier ( MPI_COMM_WORLD , mpicode )
        call initchemkin
     end if
     if ( inp % reac_selec == MIL .or. inp % reac_selec == hybrid_PSR_MIL ) then ! Model Intermittent Lagrangien
        call stoichiometric_mf (inp)
!        call load_table (inp)
     end if
  end if

  ! initialisation of EGlib
  if ( vis .and. eglib ) then
     if ( .not. reaction .and. rank == rank_default ) call ckintp39 ! overwrite the previous chem.bin file from tranfit with 3.9 interpreter
     if ( rank == rank_default ) call tran   ! create tran.bin file
     if ( rank == rank_default ) call egfrmc ! convert tran.bin into Linkeg
     call mpi_barrier ( MPI_COMM_WORLD , mpicode )
     call init_eglib ( inp , rweg , iweg )
  end if


  ! transform restart fields (need to be placed here if need topology !)


  ! taking the times
  t_in_loop = MPI_WTIME()
  tcomm     = 0.0_dp
  t_par     = tcomm


  ! main loop
  if ( .not. inp % temp_avg .and. .not. inp % spat_avg ) then ! instantaneous plots


     ! allocate variables
     call alloc_inst_type (inst)

     if ( inp % read_stat ) then
        if ( rank == rank_default ) write (*,*) 'reading previous statistics ...'
        call read_inst_type ( inp , inst )
        call vtr_paraview ( inp % end_file + 1 , inp , thd , adi , sim , grid , inst , reyavg , favavg , stat )
        call pvtr_paraview ( inp % end_file + 1 , inp )
     else
        if ( rank == rank_default ) write (*,*) 'reading files'
        do ifile = inp % start_file , inp % end_file , inp % skip_file

           if ( rank == rank_default .and. mod (ifile,inp % itshow) == 0 ) then
              call address_file ( ifile , inp , adr_file )
              write (*,format_exit) trim(adr_file)
           end if

           t_out_loop = MPI_WTIME()
           t_par = t_par + tcomm

           call read_stat ( ifile , inp , thd , adi , sim , grid , rweg , iweg , inst )
           call vtr_paraview ( ifile , inp , thd , adi , sim , grid , inst , reyavg , favavg , stat )
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
     if ( inp % second_loop ) then
        call alloc_reyfluc_type (reyfluc)
        call alloc_favfluc_type (favfluc)
     end if
     call alloc_stat_type (stat)


     ! open ASCII files
     call open_convergence ( inp , adi , grid )
     if ( inp % pdfs )     call open_pdfs ( inp , adi , grid )
     if ( inp % spectras ) call open_spectras ( inp , adi , grid )
     if ( inp % cond_avg ) call open_cond_avgs ( inp , adi , grid )


     ! reading previous stats
     if ( inp % read_stat ) then


        if ( rank == rank_default ) write (*,*) 'reading previous statistics ...'
        call read_reyavg_type ( inp , reyavg )
        call read_favavg_type ( inp , favavg )
        call read_stat_type ( inp , stat )


     else

        if ( rank == rank_default ) write (*,format_header) 'reading file' , 't_run (min)' , 'para.eff. (%)'

        ! first reading: calculate only the reynolds and favre averages
        if ( rank == rank_default ) write (*,*) 'entering first loop ...'
        do ifile = inp % start_file , inp % end_file , inp % skip_file

           t_out_loop = MPI_WTIME()
           t_par = t_par + tcomm

           if ( mod (ifile,inp % itshow) == 0 ) then
              if ( rank == rank_default ) then
                 call address_file ( ifile , inp , adr_file )
                 write (*,format_exit) trim(adr_file) , ( t_out_loop - t_in_loop ) / 60.0_dp , &
                                       100.0_dp * ( 1.0_dp - tcomm / ( t_out_loop - t_in_loop ) )
              end if
              tcomm = 0.0_dp
           end if

           call read_stat ( ifile , inp , thd , adi , sim , grid , rweg , iweg , inst )
           call averages ( ifile , inp , thd , adi , sim , grid , dt , inst , reyavg , favavg )
           if ( inp % pdfs )     call pdfs ( ifile , inp , thd , &
            & adi , sim , grid % dx_i , grid % dy_i , grid % dz_i , dt , inst , favavg , stat )
           if ( inp % spectras ) call spectras ( ifile , inp , &
            & thd , sim , grid , inst , stat )

        end do

        call wall_delta_plus ( adi , reyavg % rho , reyavg % mu , favavg % ux , grid )

        ! second reading: calculate fluctuations and statistical variables from averages
        if ( inp % second_loop ) then
           if ( rank == rank_default ) write (*,*) 'entering second loop ...'
           do ifile = inp % start_file , inp % end_file , inp % skip_file

              t_out_loop = MPI_WTIME()
              t_par = t_par + tcomm

              if ( mod (ifile,inp % itshow) == 0 ) then
                 if ( rank == rank_default ) then
                    call address_file ( ifile , inp , adr_file )
                    write (*,format_exit) trim(adr_file) , ( t_out_loop - t_in_loop ) / 60.0_dp , &
                                          100.0_dp * ( 1.0_dp - tcomm / ( t_out_loop - t_in_loop ) )
                 end if
                 tcomm = 0.0_dp
              end if

              call read_stat ( ifile , inp , thd , adi , sim , grid , rweg , iweg , inst )
              call fluctuations ( inp , thd , adi , grid % dx_i , grid % dy_i , grid % dz_i ,       &
                                  inst , reyavg , favavg , reyfluc , favfluc )
              call stats ( ifile , inp , thd , adi , sim , grid , dt , &
                           inst , reyavg , favavg , reyfluc , favfluc , stat )
              call plot_convergence ( ifile , inp , adi , sim , grid , t , dt , stat )

           end do
        end if ! second loop


        ! writing binary files to recover the stats
        if ( rank == rank_default ) write (*,*) 'writing statistical files ...'
!        call write_reyavg_type (reyavg)
!        call write_favavg_type (favavg)
!        call write_stat_type (stat)


     end if ! reading previous stats


     ! writing paraview files
     if ( rank == rank_default ) write (*,*) 'writing paraview files ...'
     call vtr_paraview ( inp % end_file , inp , thd , adi , sim , grid , inst , reyavg , favavg , stat )
     call pvtr_paraview ( inp % end_file , inp )
     call pvd_paraview ( inp , adi , sim , t )


     ! writing ASCII files
     if ( ( inp % pdfs .or. inp % spectras .or. inp % cond_avg ) .and. rank == rank_default ) &
        write (*,*) 'writing ASCII files ...'
     if ( inp % pdfs )     call plot_pdfs ( inp , adi , grid , favavg , stat )
     if ( inp % spectras ) call plot_spectras ( inp , adi , sim , grid , t , stat )
     if ( inp % cond_avg ) call plot_cond_avgs ( inp , adi , grid , favavg , stat )


     ! close ASCII files
     if ( .not. inp % read_stat ) then
        call close_convergence ( inp , grid )
        if ( inp % pdfs )     call close_pdfs ( inp , grid )
        if ( inp % spectras ) call close_spectras ( inp , grid )
        if ( inp % cond_avg ) call close_cond_avgs ( inp , grid )
     end if


  end if ! instantaneous/statistics plots


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
     write (*,*)
  end if

  ! MPI exit
  call end_mpi

end program postcosiness
