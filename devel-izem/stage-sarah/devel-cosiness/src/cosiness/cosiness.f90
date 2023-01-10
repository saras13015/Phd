!------------------------------------------------------------------------------
! PROGRAM COSINESS
!------------------------------------------------------------------------------
!> \brief Solver.
!!
!! Contains the main program.
!!
!!   Input files:\n
!!     * chem.inp chemical species, chemical reactions and kinetic schemes\n
!!     * grid.dat the computational domain with the mesh and boundary conditions\n
!!     * input.datinput data for the calculation
!! (initial condition, problem, numeric, compute parameters, ...)\n
!!     * therm.datthermodynamic table data(cp, h, s) for all the species used\n
!!     * tran.dat transport data(mu, diff) for all the species used\n
!!     * syntturb.dat   external input datas to create synthetic turbulence
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
program cosiness


    use parameters
    use type_thd
    use input
    use adim
    use parallel
    use hdf5
    use h5lt
    use IOHDF5
    use eg_lib
    use weno
    use deriv
    use thdtools
    use tools
    use ICs
    use solver
    use solver_reaction

    implicit none

    type(inp_type) :: inp
    type(adi_type) :: adi
    type(thd_type) :: thd
    type(inp_grid) :: grid

    real(dp),allocatable,dimension(:,:,:,:) :: v,ha
    real(dp),allocatable,dimension(:,:,:) :: T,W_i,cp

    ! LES
    real(dp),allocatable,dimension(:,:,:) :: mu_SGS

    ! EGlib
    integer(ip),allocatable,dimension(:) :: iweg
    real(dp),allocatable,dimension(:) :: rweg

    ! number of iterations in time
    integer(ip) :: ok,ite,nb_restart,nb_stat

    ! time measurement
    integer,dimension(8) :: val
    real(dp) :: t_wall,t_par,t_tmp1,t_tmp2,t_in_loop,t_out_loop
    real(dp) :: time,dtmin,DTmax
    real(dp),dimension(4) :: dt_4
    logical :: loop,write_stat,existed
    character(len_default) :: adr_file
    character(len_default),parameter :: format_times= '(I8,1X,10(1X,1PE18.8))'
    character(len_default),parameter :: format_times_head='(A8,1X,10(1X,A18))'


    ! Check if the file_stop exists before
    inquire(file=trim(file_stop),exist=existed)
    if (existed) then
        open(unit=unit_stop,file=trim(file_stop))
        close(unit=unit_stop,status="delete")
    endif
    existed=.false.

! MPI(Message Passing Interface : for parallel computing) initialisation
    call init_mpi

    ! COSINESS version and date
    if(rank==rank_default) then
        write(*,*)
        write(*,*) '============================'
        write(*,*) 'COSINESS'
        call date_and_time(values=val)
        print "(a,2(i2 .2,a),i4,a,3(i2 .2,a),a) ",&
        &     " ",val(3),"/",val(2),"/",val(1)," ",val(5),":",val(6),":",val(7)
        write(*,*) '============================'
        write(*,*)
    endif

    ! useful constants
    call constants(adi)

    ! initialisation of CHEMKIN interpreter
    if(rank==rank_default) call tranfit ! create todo.thd and a "dummy" binary file
    call mpi_barrier(MPI_COMM_WORLD,mpicode)

    ! evaluation of thd parameters
    call thd_reset(thd)
    call thd_init (thd,adi,file_tranfit) ; nreac=thd%nreac
    !  if(rank==rank_default) call thd_write(thd,file_thd) ! not necessary

    ! ! index of the main species(if H2 or O2 are not respectively the fuel and oxydizer, you must change them)
    ! call thd_species_index(thd,'H2',inp%index_H2)   ! Fuel
    ! call thd_species_index(thd,'O2',inp%index_O2)   ! Oxydizer
    ! call thd_species_index(thd,'H2O',inp%index_H2O) ! Product
    ! call thd_species_index(thd,'N2',inp%index_N2)   ! Diluent

    ! read file input.dat
    call readinp(adi,thd,inp)

    ! read file grid.dat
    call readgrid(adi,grid)

    ! create MPI topology
    call topology

    ! index for every subdomain
    call subdomain

    ! looking for neighboors
    call neighborhood

    ! evaluation of the computational mesh and metrics
    call metrics(grid%xt,grid%yt,grid%zt,grid%x,grid%y,grid%z,grid%dx_i,grid%dy_i,grid%dz_i,grid%delta)

    ! ! read external input data from the file syntturb.dat
    ! call readinpext(adi,grid,inp)

    ! create necessary repertories to save datas
    if(rank==rank_default) call createrep(inp)

    ! creation of derivated types for communication
    call type_conserved
    if(vis) then
        call type_derivative
        if(LES) call type_derivative_SGS
    endif

    ! definition of WENO parameters
    call wenopar(inp%opt,inp%optord)

    ! initialisation of reactive variables
    if(reaction) then 
        if(inp%reac_selec==PSR  .or. &
            inp%reac_selec==hybrid_PSR_MIL) then ! CHEMKIN library
            if(rank==rank_default) call ckintp39 ! overwrite the previous chem.bin file from tranfit with 3.9 interpreter
            call mpi_barrier(MPI_COMM_WORLD,mpicode)
            call initchemkin
        endif
        if(inp%reac_selec==MIL .or. inp%reac_selec==hybrid_PSR_MIL) then ! Model Intermittent Lagrangien
            call stoichiometric_mf(inp)
            !  call load_table(inp)
        endif
    endif

    ! initialisation of EGlib
    if(vis .and. eglib) then
        if(.not. reaction .and. rank==rank_default) call ckintp39 ! overwrite the previous chem.bin file from tranfit with 3.9 interpreter
        if(rank==rank_default) call tran
        if(rank==rank_default) call egfrmc
        call mpi_barrier(MPI_COMM_WORLD,mpicode)
        call init_eglib(inp,rweg,iweg)
    endif

    ! prepare the stats
    if(inp%stat) then
        call locate_stats(inp,adi,grid)
    endif

    if(inp%nprobes > 1) then
        call locate_probes(inp,adi,grid)
    endif

    ! opening probes files
    call open_probes(inp,adi,grid)

    ! opening the time files to store the iterations, time and dtmin
    if(rank==rank_default) then

        adr_file=trim(dir_parent) // trim(file_time_rest)
        open(unit=unit_time_rest,file=adr_file,status='new',iostat=ok)
        if(ok /= 0) call abort_mpi('ERROR opening ' // trim(adr_file))
        close(unit_time_rest)

        adr_file=trim(dir_parent) // trim(file_time_stat)
        open(unit=unit_time_stat,file=adr_file,status='new',iostat=ok)
        if(ok /= 0) call abort_mpi('ERROR opening ' // trim(adr_file))
        close(unit_time_stat)

    endif

    ! initializing some parameters
    if(inp%read_restart) then
        nb_stat=inp%number_stat
        nb_restart=inp%number_restart
    else
        nb_stat=0
        nb_restart=0
    endif
    time=inp%initial_time
    dtmin=1.0_dp
    dt_4=1.0_dp
    DTmax=1.0_dp
    if(inp%ini_sol) then
        loop=.false.
        nb_restart=-1
    else
        loop=.true.
    endif
    write_stat=.false.


    ! initialisation of the solution
    !  call init_klein_random_field(inp)
    !  call init_perturbation(inp)

    if(inp%read_restart) then
        call read_restart(inp,adi,thd,time,grid,T,W_i,cp,ha,v)
    else
        if(rank==rank_default) write(*,*) 'initializing the problem: ', trim(inp%init)
        call init_selector(inp,adi,thd,grid,T,W_i,cp,ha,v)
    endif


    if(rank==rank_default .and. loop) write(*,*) 'entering the main loop ...'
    if(rank==rank_default .and. loop) write(*,format_times_head) &
    'ite','time','dtmin',    &
    'dt_CFL*','dt_mass_diff*','dt_therm_diff*','dt_chem*'


    ! taking the times
    t_in_loop=MPI_WTIME()
    t_wall=t_in_loop
    t_tmp1=t_in_loop
    t_tmp2=time
    tcomm =0.0_dp
    t_par =tcomm

    allocate(mu_SGS(sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)    ,& 
        stat=ok)
    if(ok > 0) call abort_mpi('error allocate COSINESS.f90')

    mu_SGS=0.0_dp


    ! main loop
    ite=0
    do while(loop)

        ite=ite + 1
        if(inp%itlim .and. ite >= inp%itmax) then
            loop=.false.
            if(rank==rank_default) write(*,*) 'number of iterations has been reached'
        endif

        ! random number evaluation every 1000 time steps
        if(ite==1 .or. mod(ite,1000)==0) call randomize(inp)

        ! calculate time step
        call timestep(inp,adi,thd,&
            grid,cp,W_i,T,&
            v,mu_SGS,&
            DTmax,dt_4,dtmin)

        ! modify dtmin to adjust saving times
        if(time + dtmin > inp%timing(nb_stat+1) .and.   &
            time - dtmin < inp%timing(nb_stat+1)) then
            dtmin=inp%timing(nb_stat+1) - time
            if(inp%stat) write_stat=.true.
            if(nb_stat+1 >= inp%nstat) then
                loop=.false.
                if(rank==rank_default) write(*,*) 'number of stat files has been reached'
            endif
        endif


        if(reaction) then ! Strang splitting

            call reaction_selec(inp,thd,adi,grid,time,0.5_dp * dtmin,&
                mu_SGS,T,W_i,cp,ha,v,DTmax)
            call rk3(inp,thd,adi,time,dtmin,grid,mu_SGS,&
                rweg,iweg,T,W_i,cp,ha,v)
            time=time + dtmin
            call reaction_selec(inp,thd,adi,grid,time,0.5_dp * dtmin,&
                mu_SGS,T,W_i,cp,ha,v,DTmax)

        else

            !  if(nvv > 0 .and. LES) & 
            ! call source_sink(inp,thd,adi,time,0.5_dp * dtmin,&
            !grid,mu_SGS,T,W_i,cp,ha,v)
            call rk3(inp,thd,adi,time,dtmin,grid,mu_SGS,&
                rweg,iweg,T,W_i,cp,ha,v)
            time=time + dtmin
            !  if(nvv > 0 .and. LES) & 
            ! call source_sink(inp,thd,adi,time,0.5_dp * dtmin,&
            !grid,mu_SGS,T,W_i,cp,ha,v)

        endif


        ! For synthetic turbulent generator, shift the bc random field
        !     call bc_shift_fill_rfield(inp)


        !     ! clean the exit
        !     if(ndim==3 .and. mod(ite,1000)==0) call clean_exit(inp,adi,thd,time,dtmin,x,y,z,&
        !dx_i,dy_i,dz_i,T,W_i,cp,ha,v)


        !     ! write probes files
        !     call plot_probes(time,dtmin,inp,adi,grid,v)


        ! storing stats when conditions satisfy
        if(write_stat) then

            if(rank==rank_default) write(*,format_times) ite,&
            time * adi%time_ref,   &
            dtmin * adi%time_ref,  &
            dt_4 / dtmin

            ! write probes files
            call plot_probes(time,dtmin,inp,adi,grid,v)

            nb_stat=nb_stat + 1
            call save_stats(nb_stat,time,dtmin,dt_4,inp,adi,v)
            write_stat=.false.

            if(ndim < 3) then ! in 1D and 2D problems store the whole domain
                nb_restart=nb_restart + 1
                call save_restart(nb_restart,time,dtmin,dt_4,adi,v)
            endif

            !  if(rank==rank_default) write(*,format_times_head) &
            !     'ite','time','dtmin',    &
            !     'dt_CFL*','dt_mass_diff*','dt_therm_diff*','dt_chem*'

        endif


        ! storing restart files at selected iterations
        if(mod(ite,inp%freq_restart)==0) then
            nb_restart=nb_restart + 1
            call save_restart(nb_restart,time,dtmin,dt_4,adi,v)

            if(rank==rank_default) write(*,format_times_head) &
            'ite','time','dtmin',    &
            'dt_CFL*','dt_mass_diff*','dt_therm_diff*','dt_chem*'
        endif


        ! evaluate times every iteration including savings, etc.
        call mpi_barrier(MPI_COMM_WORLD,mpicode) ! to synchronize time for each rank
        t_wall=MPI_WTIME()
        t_par=t_par + tcomm


        ! show information at the screen
        if(mod(ite,inp%itshow)==0 .and. rank==rank_default) then
            t_tmp1=t_wall - t_tmp1
            t_tmp2=time - t_tmp2
            if(rank==rank_default) write(*,format_times) ite,   &
            time * adi%time_ref, &
            dtmin * adi%time_ref,&
            dt_4 / dtmin
            call counter(time / inp%timing(inp%nstat),  &
                t_tmp2 / inp%timing(inp%nstat),&
                t_tmp1,tcomm,t_wall-t_in_loop)
            t_tmp1=t_wall
            t_tmp2=time
            tcomm=0.0_dp
        endif


        ! cheking the wall time and STOP file every 20 iterations
        !     if((t_wall - t_in_loop) > 0.95_dp *(inp%walltime)) then
        if((t_wall - t_in_loop) > inp%walltime) then
            loop=.false.
            if(rank==rank_default) write(*,*) 'walltime has been reached'
        endif
        if(mod(ite,20)==0) then
            inquire(file=trim(file_stop),exist=existed)
            if(existed) then
                loop=.false.
                if(rank==rank_default) write(*,*) 'request to stop the program'
            endif
        endif

        inquire(file='CFL.inp',exist=existed)
        if(existed) then
            open(unit=unit_inp,file='CFL.inp')
            read(unit_inp,*) inp%CFL
            if(rank==rank_default) write(*,*) 'WARNING: Change CFL value ',CFL,' to ',inp%CFL
            CFL=inp%CFL
            call mpi_barrier(MPI_COMM_WORLD,mpicode) ! Wait all procs before to close and delete file
            close(unit=unit_inp,status="delete")
        endif


        !     call saving_ascii2(time,dtmin,thd,adi,x,y,z,T,W_i,cp,v)

    end do ! main loop

    if(rank==rank_default .and. .not. loop) write(*,*) 'exit the main loop ...'


    ! measure the loop time
    t_out_loop=MPI_WTIME()

    !  call saving_ascii2(time,dtmin,thd,adi,x,y,z,T,W_i,cp,v)
    nb_restart=nb_restart + 1
    call save_restart(nb_restart,time,dtmin,dt_4,adi,v)

    ! this is necessary to avoid closing temporal files
    call mpi_barrier(MPI_COMM_WORLD,mpicode) ! not usefull anymore ?

    ! closing the temporal files, this is VERY IMPORTANT!
    !  call close_probes(inp,x,y,z)
    !  close(unit_saveascii) ! old variable NOT USED anymore.
    !  close(unit_time_rest)
    !  close(unit_time_stat)


    if(rank==rank_default) then
        t_out_loop=t_out_loop - t_in_loop
        t_in_loop=100.0_dp *(1.0_dp - t_par / t_out_loop)
        write(*,*) 'program finished in ',t_out_loop / 60.0_dp,' minutes'
        write(*,*) 'communications time was ',t_par / 60.0_dp,' minutes'
        write(*,*) 'average parallel efficiency of ',t_in_loop,'%'
        write(*,*)
    endif

    deallocate(v,ha,T,W_i,cp,mu_SGS) ! deallocate initialized variables


    ! MPI exit
    call end_mpi

end program cosiness