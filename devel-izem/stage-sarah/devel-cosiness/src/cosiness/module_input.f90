!------------------------------------------------------------------------------
! MODULE: input
!------------------------------------------------------------------------------
!> \brief Input module.
!!
!! This module contains all the operations related to parallel
!! communications through the Message Passing Interface (MPI) library.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module input


  use parameters
  use parallel
  use adim
  use type_thd


  implicit none


  integer (ip) , parameter , private :: nplanemax  = 20    , & !< maximum number of planes
                                        nstatmax   = 10000 , & !< maximum number of stat files
                                        nprobesmax = 100   , & !< maximum number of probes
                                        nvarsmax   = 200       !< maximum number of stat variables


!> \brief External input datas (simulations, experiments, ...) derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type ext_type


     !> Synthetic turbulence generator
     real (dp) , dimension (:)         , allocatable   :: uin , vin , win !< inlet velocity components
     real (dp) , dimension (:)         , allocatable   :: T   !< temperature
     real (dp) , dimension (:)         , allocatable   :: tke !< turbulent kinetic energy
     real (dp) , dimension (:)         , allocatable   :: R12 , R23 , R13 !< Reynolds tensor
!     real (dp) , dimension (:,:,:)     , allocatable   :: Rij !< Reynolds tensor

     real (dp) , dimension (:)         , allocatable   :: A11 , A21 , A31 , & ! < Cholesky decomposition
                                                                A22 , A32 , &
                                                                      A33


     !> Klein digital filter
     integer (ip) :: nxls , nyls , nzls !< size of the gaussian length scale for Klein digital filter
     integer (ip) :: nxfw , nyfw , nzfw !< Klein filter width (such that nxfw >= 2*nxls)
     real (dp) , dimension (:,:,:,:)   , allocatable   :: rand !< random number fields for the ndim direction
     real (dp) , dimension (:,:,:)     , allocatable   :: bijk !< convolution of the 3 one-dimensional filters


     !> Table: Reversed chemical time scale in function of the mixture fraction
     !  Tab features: rctmax = R_eversed C_hemical T_ime scale is max
     real (dp)                     :: rctmax , z_rctmax
     integer (ip)                  :: index_rctmax , nval

     real (dp) , dimension (:,:)       , allocatable :: table


  end type ext_type


!> \brief Perturbation derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type perturb_type


     !> Blowing and suction function for helping transition laminar to turbulence
     real (dp) , dimension (:)         , allocatable   :: yl , zl , tm !< ponderation factor of Fourier decomposition
     integer (ip)                                      :: ylmax , zlmax , mmax !< summation limit


  end type perturb_type



!> \brief Input derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type inp_type


     logical               :: ind_files !< write individual binary files for each process

     integer (ip)          :: dim !< dimension

     integer (ip)          :: nrv !< number of reactive variables

     integer (ip)          :: npv !< number of passive variables

     integer (ip)          :: nvv !< number of variance variables

     logical               :: ini_sol !< store the initial solution

     character (len_default) :: init !< name of the initial condition

     logical               :: perturb !< perturbation

     !> filter variables
     logical               :: filter
     real (dp)             :: fil_xini , fil_xend

     !> viscous variables
     logical               :: vis
     logical               :: vis_dt

     !> LES variables
     logical                 :: LES , tau_iso_SGS_switch
     character (len_default) :: mu_SGS_model , tau_iso_SGS_model
     real (dp)               :: mu_SGS_factor , tau_iso_SGS_factor , &
                                Pr_SGS , Sc_SGS , &
                                percent_weight_SGS

     !> reactive variables
     logical                        :: reaction
     character (len_default)        :: reac_selec
     logical                        :: filter_reaction
     integer (ip)                   :: index_H2 , index_N2 , index_O2 , index_H2O
     real (dp) , dimension (nrvmax) :: Y1f , Y0o

     !> EGLIB variables
     logical               :: eglib
     integer (ip)          :: eglib_precision

     logical               :: forcesum !< force sumY=1

     !> stability criteria
     logical               :: fixeddt
     real (dp)             :: dt
     real (dp)             :: CFL , CFLmin
     real (dp)             :: Fo
     real (dp)             :: CPR
     real (dp)             :: dtmin_allowed
     logical               :: dtlimit

     !> WENO variables
     logical               :: weights
     logical               :: opt
     integer (ip)          :: optord
     real (dp)             :: percent_weight

     !> print iterations at the screen
     integer (ip)          :: itshow
     integer (ip)          :: itmax
     logical               :: itlim

     integer (ip)          :: walltime !< wall time (in minutes)

     !> restart files variables
     logical                               :: read_restart
     character (len_default)               :: name_restart
     integer (ip)                          :: number_restart
     integer (ip)                          :: freq_restart
     real (dp)                             :: time_offset , initial_time , dtime

     !> stat files variables
     logical                               :: stat
     integer (ip)                          :: nstat , number_stat
     integer (ip)                          :: nvolume , nxystat , nxzstat , nyzstat
!     integer (ip) , dimension (nplanemax)  :: i_volmin , j_volmin , k_volmin ! USELESS ? not used...
!     integer (ip) , dimension (nplanemax)  :: i_volmax , j_volmax , k_volmax ! USELESS ? not used...
     integer (ip) , dimension (nplanemax)  :: i_yzstat , j_xzstat , k_xystat
     real (dp)                             :: dim_length_coord
     real (dp) , dimension (nplanemax)     :: x_volmin , y_volmin , z_volmin
     real (dp) , dimension (nplanemax)     :: x_volmax , y_volmax , z_volmax
     real (dp) , dimension (nplanemax)     :: x_yzstat , y_xzstat , z_xystat
     real (dp) , dimension (nstatmax)      :: timing

     !> post.dat file _additional_ variables
     integer (ip)                          :: nprocx , nprocy , nprocz
     integer (ip)                          :: ghostptx , ghostpty , ghostptz
     logical                               :: read_BE , write_BE
     logical                               :: trafo , nondim_grid
     real (dp)                             :: length_ref
     real (dp)                             :: umax
     real (dp)                             :: umin
     real (dp)                             :: uc
     real (dp)                             :: rhom
     real (dp)                             :: mum
     real (dp)                             :: Yst

     logical                               :: dis_rate_therm

     integer (ip)                          :: nss , index_fuel , index_oxidizer , index_Da , index_scalar

     integer (ip)                          :: nbins

     character (len_default)               :: plane_type
     integer (ip)                          :: plane_number

     logical                               :: temp_avg
     logical                               :: spat_avg
     logical                               :: cond_avg
     logical                               :: pdfs
     logical                               :: spectras

     logical                               :: similarity_x

     logical                               :: second_loop

     logical                               :: read_stat

     integer (ip)                          :: start_file , end_file , skip_file

     real (dp)                             :: s_x , e_x , s_y , e_y , s_z , e_z
!     integer (ip)                          :: sx , ex , nx , sy , ey , ny , sz , ez , nz ! useless ?

     integer (ip)                                     :: nvars !< number of statistical variables
     character (len_default) , dimension (nvarsmax)   :: var_name !< names of the statistical variables

     real (dp) , dimension (3)                        :: corrspec_coord !< autocorrelation

     !> probes, used for ASCII files, e.g. convergences, pdfs, spectras...
     integer (ip)                                     :: nprobes
     character (len_default) , dimension (nprobesmax) :: probe_name
     real (dp) , dimension (nprobesmax,3)             :: probe_coord


     !> input sub derived type
     type (ext_type)                       :: ext
     type (perturb_type)                   :: ptb


  end type inp_type


!> \brief Grid derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type inp_grid

     real (dp) , allocatable , dimension (:)       :: xt   , yt   , zt   !< absolute xyz-coordinate array in the hole domain
     real (dp) , allocatable , dimension (:)       :: x    , y    , z    !< xyz-coordinate array in the subdomain
     real (dp) , allocatable , dimension (:)       :: dx_i , dy_i , dz_i !< inverted dx array in the subdomain
     real (dp) , allocatable , dimension (:,:,:)   :: delta              !< filter width

     real (dp) , allocatable , dimension (:,:,:,:) :: YZ_dxplus, XZ_dxplus , XY_dxplus !< adimensionnal distance wall 

  end type inp_grid


contains


!> \brief Read the file grid.dat.
!!
!! Read the file grid.dat, define absolute coordinates (xt,yt,zt),
!! adimensionalize them and check perdiodic BCs.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine readgrid ( adi , grid )


    type (adi_type) , intent (in)                              :: adi  !< non-dimensional derived type
    type (inp_grid) , intent (inout)                           :: grid !< grid derived type

    logical , parameter          :: reorganisation = .false.
    character (len_default)      :: word
    integer (ip)                 :: i , j , k , ok , ndir
    real (dp)                    :: Lref_i

    character (len_default) , parameter :: format_exit1 = '(3(1X,A,I7))'
    character (len_default) , parameter :: format_exit2 = '(1X,A,2(EN25.10E3))'


    open ( unit_grid , file = file_grid , status = 'old' , action = 'read' , iostat = ok )
    if ( ok /= 0 ) call abort_mpi ( 'error opening ' // trim (file_grid) )

    read (unit_grid,*) ! number of processus in x-direction (0 if automatic)
    read (unit_grid,*) ! number of processus in y-direction (0 if automatic)
    read (unit_grid,*) ! number of processus in z-direction (0 if automatic)

    do i = 1 , nneigh
       read (unit_grid,*) bc(i) ! +/-(x,y,z) boundary conditions
    end do

    ! count number of point in each direction
    do
       read (unit_grid,*) word
       if (word == 'x-direction') exit
    end do

    ! x-direction
    i = 0
    ndir = 1
    do
       read (unit_grid,*) word
       if (word == 'y-direction') exit
       i = i + 1
    end do
    ntx = i ! count number of point in x-direction

    ! y-direction
    j = 0
    ndir = ndir +1
    do
       read (unit_grid,*) word
       if (word == 'z-direction') exit
       j = j + 1
    end do
    nty = j

    ! z-direction
    k = 0
    ndir = ndir +1
    do
       read (unit_grid,*) word
       if (word == 'END') exit
       k = k + 1
    end do
    ntz = k

    if (ndir /= 3) &
       call abort_mpi ( 'Missing direction in grid.dat : ndir = ' // trim(str(ndir)) )


    ! go back up to the file
    rewind ( unit_grid )


    allocate ( grid % xt (1-ng:ntx+ng) , &
               grid % yt (1-ng:nty+ng) , &
               grid % zt (1-ng:ntz+ng) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate readgrid')

    grid % xt = 0.0_dp
    grid % yt = 0.0_dp
    grid % zt = 0.0_dp

    read (unit_grid,*) dims(3) ! number of processus in x-direction
    read (unit_grid,*) dims(2) ! number of processus in y-direction
    read (unit_grid,*) dims(1) ! number of processus in z-direction

    do i = 1 , nneigh
       read (unit_grid,*) word ! boundary conditions
    end do

    read (unit_grid,*) word
    if (word /= 'x-direction') call abort_mpi ('error: key word ''x-direction'' is missing')
    do i = 1 , ntx
       read (unit_grid,*) grid % xt(i)
    end do

    read (unit_grid,*) word
    if (word /= 'y-direction') call abort_mpi ('error: key word ''y-direction'' is missing')
    do j = 1 , nty
       read (unit_grid,*) grid % yt(j)
    end do

    read (unit_grid,*) word
    if (word /= 'z-direction') call abort_mpi ('error: key word ''z-direction'' is missing')
    do k = 1 , ntz
       read (unit_grid,*) grid % zt(k)
    end do

    read (unit_grid,*) word
    if (word /= 'END') call abort_mpi ('error: key word ''END'' is missing')

    close (unit_grid)


    Lx = grid % xt (ntx) - grid % xt (1)
    Ly = grid % yt (nty) - grid % yt (1)
    Lz = grid % zt (ntz) - grid % zt (1)


    ! Display the points boundaries
    if ( rank == rank_default ) then
       write (*,*)
       write (*,format_exit1) 'MPI processes = ' , nproc
       write (*,format_exit1) 'Topology:    nx  = ' , dims(3) , ', ny  = ' , dims(2) , ', nz  = ' , dims(1)
       write (*,format_exit1) 'Domain size: ntx = ' , ntx     , ', nty = ' , nty     , ', ntz = ' , ntz

       write (*,format_exit2) 'xmin , xmax = ' , grid % xt (1) , grid % xt (ntx)
       write (*,format_exit2) 'ymin , ymax = ' , grid % yt (1) , grid % yt (nty)
       write (*,format_exit2) 'zmin , zmax = ' , grid % zt (1) , grid % zt (ntz)
       write (*,*)
    end if


    ! Checking for periodic boundary conditions
    do i = 1 , 5 , 2
       if ( bc (i) == periodic .and. bc(i+1) /= bc(i) ) &
          call abort_mpi ('boundary conditions must appear in pairs')
    end do

    period(:) = .false.
    if ( ( bc(W) == bc(E) ) .and. ( bc(W) == periodic ) ) period(3) = .true.
    if ( ( bc(S) == bc(N) ) .and. ( bc(S) == periodic ) ) period(2) = .true.
    if ( ( bc(B) == bc(F) ) .and. ( bc(B) == periodic ) ) period(1) = .true.


    ! Adimensionalize here to avoid adimensionalizing all other variables
    Lref_i = 1.0_dp / adi % L_ref
    grid % xt (:) = grid % xt (:) * Lref_i
    grid % yt (:) = grid % yt (:) * Lref_i
    grid % zt (:) = grid % zt (:) * Lref_i
    Lx = Lx / adi % L_ref
    Ly = Ly / adi % L_ref
    Lz = Lz / adi % L_ref


  end subroutine readgrid


!> \brief Read the file input.dat.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine readinp ( adi , thd , inp )


    type (adi_type) , intent (in)                                  :: adi !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd !< thermodynamic derived type
    type (inp_type) , intent (inout)                               :: inp !< input derived type


    integer (ip) :: ok , l , ind
    logical      :: loop
    real (dp)    :: Lref_i , time_ref_i


    ! number of bits for binary files depending on the compiler:
    ! i) __INTEL_COMPILER,
    ! ii) __GFORTRAN__
    nbit = 8 ! by default: works with gfortran and IBM specific compiler

! #ifdef __INTEL_COMPILER
!     nbit = 2
! #endif

    open ( unit = unit_inp , file = file_inp , status = 'old' , action = 'read' , iostat = ok )
    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (file_inp))


    read ( unit_inp , * ) inp % dim                            !dimension of the problem
    read ( unit_inp , * ) inp % nrv                            !reactive variables (= nb of SPECIES in chem.inp)
    read ( unit_inp , * ) inp % npv                            !passive variables (= 0 or 1)
    read ( unit_inp , * ) inp % nvv                            !variance variables for LES (= 0 or 1)
    read ( unit_inp , * ) inp % ind_files                      !manipulate individual files (false by default)
    read ( unit_inp , * ) inp % ini_sol                        !storing only the initial solution
    read ( unit_inp , * ) inp % init                           !initialization (Diffusion, Bogey, Fu, Miller, MIXChengMiller, Premix, Jet, DIFFChengMiller, Fedkiw, DMR, Ambient)
    read ( unit_inp , * ) inp % read_restart , inp % name_restart , inp % number_restart !reading a previous restart file, name, number
    read ( unit_inp , * ) inp % initial_time                   !initial simulation time
    read ( unit_inp , * ) inp % perturb                        !perturbation
    read ( unit_inp , * ) inp % filter                         !filter
    read ( unit_inp , * ) inp % fil_xini , inp % fil_xend      !   filter_start, filter_end
    read ( unit_inp , * ) inp % vis , inp % vis_dt             !viscous problem (true=NS/false=Euler), viscous time step (true by defaulut if viscous problem)
    read ( unit_inp , * ) inp % eglib , inp % eglib_precision  !eglib detailed transport (true/false), presision (1,2,3)
    read ( unit_inp , * ) inp % LES                            !LES simulation
    read ( unit_inp , * ) inp % mu_SGS_model , inp % mu_SGS_factor                                      !   Viscosity SGS model (Smagorinsky 0.18, Kolmogorov 1.4)
    read ( unit_inp , * ) inp % tau_iso_SGS_switch , inp % tau_iso_SGS_model , inp % tau_iso_SGS_factor !   Isotropic tensor SGS model (switch, name, factor)
    read ( unit_inp , * ) inp % percent_weight_SGS             !      weights density percent criteria for discontinuitie zones
    read ( unit_inp , * ) inp % Sc_SGS                         !   SGS Schmidt (0.7 to 1.0)
    read ( unit_inp , * ) inp % Pr_SGS                         !   SGS Prandtl (0.7 to 1.0) -> (0.3 to 0.9)
    read ( unit_inp , * ) inp % reaction , inp % reac_selec    !reactive problem (true/false), solver (PSR=dvode, MIL, hybrid_PSR_MIL)
    read ( unit_inp , * ) inp % filter_reaction                !   filter reaction zones to reduce computational costs
    read ( unit_inp , * ) inp % Y1f (1:inp % nrv + inp % npv)  !   pure fuel composition (calculate: E.R., Detector_reaction, Zm, Deq, reactivity table)
    read ( unit_inp , * ) inp % Y0o (1:inp % nrv + inp % npv)  !   pure oxydizer composition (idem)
    read ( unit_inp , * ) inp % forcesum                       !force sumY=1 ! (true by default)
    read ( unit_inp , * ) inp % fixeddt , inp % dt             !   fixed time step (true/false, value in seconde)
    read ( unit_inp , * ) inp % CFLmin , inp % CFL             !   CFL (min,max<=0.9) convective
    read ( unit_inp , * ) inp % Fo                             !   Fo (<=0.1) viscous
    read ( unit_inp , * ) inp % CPR                            !   DTpermis (<=20K) chemical
    read ( unit_inp , * ) inp % dtlimit , inp % dtmin_allowed  !   imposed time step minimum value (true/false, value) (recommand 1e-13 min)
    read ( unit_inp , * ) inp % weights , inp % percent_weight !WENO nonlinear weights (true/false), density percent criteria (value)
    read ( unit_inp , * ) inp % opt , inp % optord             !WENO optimum weights (true/false), order (1, 3, 5)
    read ( unit_inp , * ) inp % itshow                         !iterations to show information
    read ( unit_inp , * ) inp % freq_restart                   !iterations frequency to store restart files
    read ( unit_inp , * ) inp % walltime                       !walltime in minutes of calculation
    read ( unit_inp , * ) inp % itlim , inp % itmax            !iterations limit (true/false), maximum iterations (value)
    read ( unit_inp , * ) inp % time_offset                    !time to start recording stat files (seconde)
    read ( unit_inp , * ) inp % dtime                          !   dt between stat files
    read ( unit_inp , * ) inp % stat , inp % number_stat       !   store stat files, number of the last stat file
    read ( unit_inp , * ) inp % nstat                          !   number of stats
    read ( unit_inp , * ) inp % dim_length_coord               !adimensional length for coordinates

    read ( unit_inp , * ) !
    read ( unit_inp , * ) ! Storing volume files (xmin,xmax,ymin,ymax,zmin,zmax coordinates)

    inp % nvolume = 0 ; l = 1
    do
       if ( inp % nvolume > nplanemax ) call abort_mpi ('maximum number of volumes reached')
       read ( unit_inp , * , iostat = ok ) inp % x_volmin (l) , inp % x_volmax (l) , &
                                           inp % y_volmin (l) , inp % y_volmax (l) , &
                                           inp % z_volmin (l) , inp % z_volmax (l)
       if ( ok /=0 ) exit  ! Storing stat files in XY-plane (Z-coordinates)
       inp % nvolume = inp % nvolume + 1
       l = l + 1
    end do

    inp % nxystat = 0 ; l = 1
    do
       if ( inp % nxystat > nplanemax ) call abort_mpi ('maximum number of XY-planes reached')
       read ( unit_inp , * , iostat = ok ) inp % z_xystat (l)
       if ( ok /=0 ) exit  ! Storing stat files in XZ-plane (Y-coordinates)
       inp % nxystat = inp % nxystat + 1
       l = l + 1
    end do

    inp % nxzstat = 0 ; l = 1
    do
       if ( inp % nxzstat > nplanemax ) call abort_mpi ('maximum number of XZ-planes reached')
       read ( unit_inp , * , iostat = ok ) inp % y_xzstat (l)
       if ( ok /=0 ) exit  ! Storing stat files in YZ-plane (X-coordinates)
       inp % nxzstat = inp % nxzstat + 1
       l = l + 1
    end do

    inp % nyzstat = 0 ; l = 1
    do
       if ( inp % nyzstat > nplanemax ) call abort_mpi ('maximum number of YZ-planes reached')
       read ( unit_inp , * , iostat = ok ) inp % x_yzstat (l)
       if ( ok /=0 ) exit  ! Probe coordinates (name,x,y,z)
       inp % nyzstat = inp % nyzstat + 1
       l = l + 1
    end do


    ind = 1 ; loop = .true.
    inp % nprobes = 0
    inp % probe_coord = 0.0_dp
    do while (loop)

       if ( inp % nprobes > nprobesmax ) call abort_mpi ('maximum number of probes reached')
       read ( unit_inp , * ) inp % probe_name (ind)

       if ( inp % probe_name (ind) == 'END' ) then

          inp % nprobes = ind-1
          loop = .false.

       else

          backspace (unit_inp)
          read ( unit_inp , * ) inp % probe_name (ind) , &
                                inp % probe_coord (ind,:)

       end if

       ind = ind+1

    end do

    close (unit_inp)


    ! A little bit of post-treatment

    if ( inp % reac_selec == 'dvode' ) inp % reac_selec = 'PSR'

    if ( thd % nspc /= inp % nrv ) &
       call abort_mpi ('error: ' // trim(str(thd % nspc)) // ' thd % nspc /= ' // trim(str(inp % nrv)) // 'inp % nrv')

    if ( npv > 1 ) then
       write (*,*) 'the code is not prepared to support more than one passive scalar'
       write (*,*) 'please reconsider this parameter'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if

    if ( inp % stat .and. inp % dim == 3 .and. &
         ( inp % nvolume == 0 .and. inp % nxystat == 0 .and. inp % nxzstat == 0 .and. inp % nyzstat == 0 ) ) then
       write (*,*) 'error: you try to save statistic files (volumes or planes) for 3D simulation without defining volumes or planes'
       write (*,*) 'please reconsider this parameter'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if

    if ( nvv > 1 ) then
       write (*,*) 'the code is not prepared to support more than one variance scalar'
       write (*,*) 'please reconsider this parameter'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if

    if ( nvv >= 1 .and. ( nvv /= npv ) ) then
       write (*,*) 'error: nvv /= npv, you need as many npv (passive scalare) as nvv (variance)'
       write (*,*) 'please reconsider this parameter'
       ! if you want more than 1 nvv, you need to arrange the vector of conservative variable like :
       ! rho, rhou, rhov, rhow, rhoet, rhoYa, rhoPs1, rhoPs2, rhoPsn, rhoVar1, rhoVar2, ..., rhoVarn
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if

    if ( nvv >= 1 .and. ( .not. inp % LES ) ) then
       write (*,*) 'error: nvv >=1 and LES is desactivated, you need LES'
       write (*,*) 'please reconsider this parameter'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if

    if ( inp % tau_iso_SGS_switch .and. ( .not. inp % LES ) ) then
       write (*,*) 'error: input.dat -> calculate SGS isotropic tensor without activate LES'
       write (*,*) 'please reconsider this parameter'
       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
    end if

    if ( inp % eglib ) then
       nderiv = inp % dim * ( inp % dim + 2 ) +       &
                inp % dim * ( inp % nrv + inp % npv )
    else
       nderiv = inp % dim * ( inp % dim + 1 ) +       &
                inp % dim * ( inp % nrv + inp % npv )
    end if

    nderivSGS = 0
    if ( inp % LES ) then
       nderiv = nderiv +     &
                inp % dim * inp % nvv + & ! variance variable (for LES only)
                inp % dim +             & ! W_i derivatives
                1                         ! shock detector

       nderivSGS = inp % dim ! isotropic tensor derivatives
    end if

    ind_files = inp % ind_files

    ndim = inp % dim ! dimension of the problem
    nrv  = inp % nrv ! reactive variables
    npv  = inp % npv ! passive variables
    nvv  = inp % nvv ! variance variables
    nv   = niv + nrv + npv + nvv

    weno_avg        = inp % weights
    vis             = inp % vis
    vis_dt          = inp % vis_dt
    reaction        = inp % reaction
    filter_reaction = inp % filter_reaction
    eglib           = inp % eglib
    eglib_precision = inp % eglib_precision
    LES             = inp % LES
    forcesum        = inp % forcesum

    CFL    = inp % CFL
    CFLmin = inp % CFLmin
    Fo     = inp % Fo
    CPR    = inp % CPR

    inp % percent_weight = inp % percent_weight / 100.0_dp
    max_rel_weight       = inp % percent_weight

    percent_weight_SGS   = inp % percent_weight_SGS / 100.0_dp

    inp % walltime = inp % walltime * 60 ! minutes -> seconds

    if ( inp % nstat > nstatmax ) call abort_mpi ('maximum number of statistics reached')

    if ( inp % read_restart ) then
       inp % timing (1:inp % number_stat+1) = inp % time_offset
       do l = inp % number_stat+2 , inp % nstat
          inp % timing (l) = inp % timing (l-1) + inp % dtime
       end do
    else
       inp % timing (1) = inp % time_offset
       do l = 2 , inp % nstat
          inp % timing (l) = inp % timing (l-1) + inp % dtime
       end do
    end if

    ! Velocity gradient tensor & Kronecker delta Dij
    allocate ( duidxj  ( ndim , ndim )             , &
               kdelta  ( ndim , ndim )             , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate input')
    kdelta = 0.0_dp
    do l = 1 , ndim
       kdelta (l,l) = 1.0_dp
    end do

    ! Adimensionalize
    Lref_i = 1.0_dp / adi % L_ref
    time_ref_i = 1.0_dp / adi % time_ref

    inp % initial_time = inp % initial_time * time_ref_i
    inp % timing (:)   = inp % timing (:) * time_ref_i

    do l = 1 , inp % nprobes
       inp % probe_coord (l,:) = inp % probe_coord (l,:) * inp % dim_length_coord * Lref_i
    end do

    do l = 1 , inp % nvolume
       inp % x_volmin (l) = inp % x_volmin (l) * inp % dim_length_coord * Lref_i
       inp % x_volmax (l) = inp % x_volmax (l) * inp % dim_length_coord * Lref_i
       inp % y_volmin (l) = inp % y_volmin (l) * inp % dim_length_coord * Lref_i
       inp % y_volmax (l) = inp % y_volmax (l) * inp % dim_length_coord * Lref_i
       inp % z_volmin (l) = inp % z_volmin (l) * inp % dim_length_coord * Lref_i
       inp % z_volmax (l) = inp % z_volmax (l) * inp % dim_length_coord * Lref_i
    end do

    do l = 1 , inp % nxystat
       inp % z_xystat (l) = inp % z_xystat (l) * inp % dim_length_coord * Lref_i
    end do

    do l = 1 , inp % nxzstat
       inp % y_xzstat (l) = inp % y_xzstat (l) * inp % dim_length_coord * Lref_i
    end do

    do l = 1 , inp % nyzstat
       inp % x_yzstat (l) = inp % x_yzstat (l) * inp % dim_length_coord * Lref_i
    end do


  end subroutine readinp


!> \brief Define in which BCs to apply Synthetic Turbulent Generator.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine apply_stg ( apply )


    logical , intent (inout) :: apply

    integer (ip) :: i

    do i = 1 , nneigh
       if ( neigh (i) == MPI_PROC_NULL .and. &
            (    bc (W) == syntheticturb     &
            .or. bc (W) == supersonicflow    ) ) apply = .true. ; exit ! Apply, on WEST face here.
    end do


  end subroutine apply_stg


!> \brief Read the file syntturb.dat.
!!
!! This datas are come from previous simulations or experiments. 
!! And use to feed the boundary conditions to generate synthetic turbulence at the inlet(s).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine readinpext ( adi , grid , inp )


    type (adi_type)                         , intent (in)      :: adi  !< non-dimensional derived type
    type (inp_grid)                         , intent (in)      :: grid !< grid derived type
    type (inp_type)                         , intent (inout)   :: inp  !< input derived type

    integer (ip)                                   :: ok , j
    logical                                        :: loop

    ! Variables in input external file (grid1)
    integer (ip)                                   :: ntyg1
    real (dp) , allocatable , dimension (:)        :: yg1 , uing1 , Tg1 , tkeg1 , R12g1 , R23g1 , R13g1

    real (dp)                                      :: dummy , Rii , A11_i
    real (dp) , parameter                          :: twothirds = 2.0_dp / 3.0_dp   , &
                                                      eps10     = 1.0e-10_dp
    real (dp)                                      :: uref_i , uref2_i , Tref_i , Lref_i


    ! Test: existence of the synthetic turbulence condition and apply to physical BC 
    !       or IC in case you don't read restart file
    !================================================================================
    loop = .false.

    if ( .not. inp % read_restart ) then
       loop = .true. ! Apply, on the whole domain for IC
    else
       call apply_stg ( loop )
    end if

    if ( .not. loop ) return ! leave the subroutine


    ! Define variable according to the external input file (grid1)
    ! Need to be adapted, if you change or add face to generate synthetic turbulence
    !================================================================================

    open ( unit = unit_ext , file = file_syntturb , status = 'old' , action = 'read' , iostat = ok )
    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (file_syntturb))


       ! count the number of points in Y-direction
       read ( unit_ext , * ) ! Header...
       ntyg1 = 0
       do
          read ( unit_ext , * , iostat = ok ) dummy
          if ( ok /=0 ) exit
          ntyg1 = ntyg1 + 1
       end do

       ! Allocate buffer variables, for parralel
       allocate ( yg1      ( 1:ntyg1 )                                , &
                  uing1    ( 1:ntyg1 )                                , &
                  Tg1      ( 1:ntyg1 )                                , &
                  tkeg1    ( 1:ntyg1 )                                , &
                  R12g1    ( 1:ntyg1 )                                , &
                  R13g1    ( 1:ntyg1 )                                , &
                  R23g1    ( 1:ntyg1 )                                , &
!                  Rij    ( ndim , ndim , 1:ntyg1 )                  , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate readinpext')


       ! Read data from input file
       rewind ( unit_ext ) ! go back up to the file
       read ( unit_ext , * )  ! Y-direction     ux     temperature     tke     Rxy     Ryz     Rxz

       Lref_i = 1.0_dp / adi % L_ref
       do j = 1 , ntyg1
          read ( unit_ext , * ) yg1 (j) , uing1 (j) , Tg1 (j) , tkeg1 (j) , R12g1 (j) , R23g1 (j) , R13g1 (j)
          ! Adimensionalize
          yg1 (j) = yg1 (j) * Lref_i
       end do

    close (unit_ext)


    ! Allocate and extrapolate variables for grid2 (grid.dat)
    !========================================================

    allocate ( inp % ext % uin    ( sy:ey )                       , &
               inp % ext % T      ( sy:ey )                       , &
               inp % ext % tke    ( sy:ey )                       , &
               inp % ext % R12    ( sy:ey )                       , &
               inp % ext % R23    ( sy:ey )                       , &
               inp % ext % R13    ( sy:ey )                       , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate inp % ext type')


    call extrapolate_data_1d ( ntyg1 , yg1 , uing1 , .true. , .true. , grid % yt , inp % ext % uin )
    call extrapolate_data_1d ( ntyg1 , yg1 , Tg1   , .true. , .true. , grid % yt , inp % ext % T   )
    call extrapolate_data_1d ( ntyg1 , yg1 , tkeg1 , .true. , .true. , grid % yt , inp % ext % tke )
    call extrapolate_data_1d ( ntyg1 , yg1 , R12g1 , .true. , .true. , grid % yt , inp % ext % R12 )
    call extrapolate_data_1d ( ntyg1 , yg1 , R23g1 , .true. , .true. , grid % yt , inp % ext % R23 )
    call extrapolate_data_1d ( ntyg1 , yg1 , R13g1 , .true. , .true. , grid % yt , inp % ext % R13 )

    deallocate ( yg1 , uing1 , Tg1 , tkeg1 , R12g1 , R13g1 , R23g1 ) !, Rij )


    ! Adimensionalize
    !================

    uref_i  = 1.0_dp / adi % u_ref
    uref2_i = uref_i * uref_i
    Tref_i  = 1.0_dp / adi % T_ref

    inp % ext % uin (sy:ey) = inp % ext % uin (sy:ey) * uref_i
    inp % ext % T   (sy:ey) = inp % ext % T   (sy:ey) * Tref_i

    inp % ext % tke (sy:ey) = inp % ext % tke (sy:ey) * uref2_i
    inp % ext % R12 (sy:ey) = inp % ext % R12 (sy:ey) * uref2_i
    inp % ext % R23 (sy:ey) = inp % ext % R23 (sy:ey) * uref2_i
    inp % ext % R13 (sy:ey) = inp % ext % R13 (sy:ey) * uref2_i


    ! Test if the neigh (i) == MPI_PROC_NULL and if BC for STG is needed
    ! This to avoid allocating useless arrays (inner process)
    !===================================================================

    loop = .false.
    call apply_stg ( loop )
    if ( .not. loop ) return ! leave the subroutine


    allocate ( inp % ext % A11    ( sy:ey )                       , &
               inp % ext % A21    ( sy:ey )                       , &
               inp % ext % A31    ( sy:ey )                       , &
               inp % ext % A22    ( sy:ey )                       , &
               inp % ext % A32    ( sy:ey )                       , &
               inp % ext % A33    ( sy:ey )                       , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate inp % ext % Aij')


    ! Cholesky decomposition
    do j = sy , ey

       Rii = twothirds * inp % ext % tke (j)

       inp % ext % A11 (j) = sqrt ( Rii )
       A11_i = 1.0_dp / max ( inp % ext % A11 (j) , eps10 )
       inp % ext % A21 (j) = inp % ext % R12 (j) * A11_i
       inp % ext % A31 (j) = inp % ext % R13 (j) * A11_i

       inp % ext % A22 (j) = sqrt ( Rii - inp % ext % A21 (j) * inp % ext % A21 (j) )
       inp % ext % A32 (j) = ( inp % ext % R23 (j) - inp % ext % A21 (j) * inp % ext % A31 (j) ) &
                           / max ( inp % ext % A22 (j) , eps10 )

       inp % ext % A33 (j) = sqrt ( Rii - inp % ext % A31 (j) * inp % ext % A31 (j) &
                           - inp % ext % A32 (j) * inp % ext % A32 (j) )

    end do


  end subroutine readinpext


!> \brief 1D Extrapolation from the external input datas to the new grid
!!
!! 1D extrapolation data from the grid1 (n1,x1,y1) to the new grid2 (n2,x2,y2)
!!
!! In case where grid2 is larger than grid1: 's_extra' and 'e_extra' are
!! logical to define if you want to extrapolate or copy the first and last 
!! value to the rest of the domain.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine extrapolate_data_1d ( n1 , x1 , y1 , s_extra , e_extra , x2 , y2 )


    integer (ip)                               , intent (in)    :: n1         !< number of points in initial grid1
    real(dp)     , allocatable , dimension (:) , intent (in)    :: x1         !< grid1 coordinate
    real(dp)     , allocatable , dimension (:) , intent (in)    :: y1         !< grid1 values to extrapolate
    logical                                    , intent (in)    :: s_extra    !< start grid2 extrapolate data from grid1
    logical                                    , intent (in)    :: e_extra    !< end   grid2 extrapolate data from grid1
    real(dp)     , allocatable , dimension (:) , intent (in)    :: x2         !< grid2 coordinate
    real(dp)     , allocatable , dimension (:) , intent (inout) :: y2         !< grid2 extrapolated values

    integer (ip) :: i , k
    integer (ip) :: si , ei ! start/end points of the intersection grid1 and grid2 (common domain)
    real(dp)     :: t


    si = sy
    ei = ey

    ! Begining of the grid2 BC
    if ( x2(si) < x1(1) ) then

       if (s_extra) then
          ! Extrapole from the initial value
          do while ( x2(si) < x1(1) )
             t = ( x2(si) - x1(1) ) / ( x1(2) - x1(1) )
             y2(si) = ( 1.0_dp - t ) * y1(1) + t * y1(2)
             si = si + 1
          end do
       else
          ! Repeat the initial value at each point
          do while ( x2(si) < x1(1) )
             y2(si) = y1(1)
             si = si + 1
          end do
       end if

    end if


    ! End of the grid2 BC
    if ( x2(ei) > x1(n1) ) then

       if (e_extra) then
          ! Extrapole from the last value
          do while ( x2(ei) > x1(n1) )
             t = ( x2(ei) - x1(n1-1) ) / ( x1(n1) - x1(n1-1) )
             y2(ei) = ( 1.0_dp - t ) * y1(n1-1) + t * y1(n1)
             ei = ei - 1
          end do
       else
          ! Repeat the last value at each point
          do while ( x2(ei) > x1(n1) )
             y2(ei) = y1(n1)
             ei = ei - 1
          end do
       end if

    end if


    ! Intersection grid1 and grid2 (common domain)
    do i = si , ei
       do k = 2, n1

          if ( x1(k-1) <= x2(i) .and. x2(i) <= x1(k) ) then
             ! grid2 point between two grid1 points: extrapole from the previous value
             t = ( x2(i) - x1(k-1) ) / ( x1(k) - x1(k-1) )
             y2(i) = ( 1.0_dp - t ) * y1(k-1) + t * y1(k)
             exit
          end if

       end do
    end do


  end subroutine extrapolate_data_1d


!> \brief Read the file post.dat
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine readpost ( adi , inp )


    type (adi_type) , intent (in)    :: adi !< non-dimensional derived type
    type (inp_type) , intent (inout) :: inp !< input derived type


    integer (ip) :: ok , i , ind
    logical      :: loop

    character (len_default) , parameter :: format = ' ( I2.2 ) '
    character (len_default)             :: ascii 
    real (dp)                           :: Lref_i


    open ( unit = unit_post , file = file_post , status = 'old' , action = 'read' , iostat = ok )
    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (file_post))

    read ( unit_post , * ) inp % nprocx             !actual number of processes in x-direction
    read ( unit_post , * ) inp % nprocy             !actual number of processes in y-direction
    read ( unit_post , * ) inp % nprocz             !actual number of processes in z-direction
    read ( unit_post , * ) inp % ghostptx           !additional phantom points for volumes in x-direction
    read ( unit_post , * ) inp % ghostpty           !additional phantom points for volumes in y-direction
    read ( unit_post , * ) inp % ghostptz           !additional phantom points for volumes in z-direction
    read ( unit_post , * ) inp % read_BE            !read big_endian binary files
    read ( unit_post , * ) inp % write_BE           !write big_endian binary files
    read ( unit_post , * ) inp % trafo              !transform the grid
    read ( unit_post , * ) inp % nondim_grid        !nondimensional grid
    read ( unit_post , * ) inp % length_ref         !reference length
    read ( unit_post , * ) inp % umax               !umax
    read ( unit_post , * ) inp % umin               !umin
    read ( unit_post , * ) inp % uc                 !uc
    read ( unit_post , * ) inp % rhom               !rho avg
    read ( unit_post , * ) inp % mum                !mu avg
    read ( unit_post , * ) inp % Yst                !Yst
    read ( unit_post , * ) inp % dis_rate_therm     !same Schmidt and Lewis equal to unity (for the diffusion coefficient)
    read ( unit_post , * ) inp % nss                !number of steady state species (for reduced mechanism)
    read ( unit_post , * ) inp % index_Da           !index species (Damkohler)
    read ( unit_post , * ) inp % index_scalar       !species number for scalar properties (0: mixture fraction)
    read ( unit_post , * ) inp % nbins              !number of bins to create pdfs, conditional averages...
    read ( unit_post , * ) inp % plane_type         !plane type: XY, XZ, YZ, restart, volume
    read ( unit_post , * ) inp % plane_number       !plane number
    read ( unit_post , * ) inp % temp_avg           !temporal average
    read ( unit_post , * ) inp % spat_avg           !spatial average
    read ( unit_post , * ) inp % cond_avg           !conditional average
    read ( unit_post , * ) inp % pdfs               !pdfs
    read ( unit_post , * ) inp % spectras           !spectras
    read ( unit_post , * ) inp % similarity_x       !similarity solutions in x-direction
    read ( unit_post , * ) inp % second_loop        !second loop
    read ( unit_post , * ) inp % read_stat          !read previous stat
    read ( unit_post , * ) inp % start_file         !start file (file restart in the times_rest.out)
    read ( unit_post , * ) inp % end_file           !end file   (idem)
    read ( unit_post , * ) inp % skip_file          !skip file  (idem)

    read ( unit_post , * ) ! space
    read ( unit_post , * ) ! autocorrelation/spectra coordinates (x,y,z; only one auto-correlation/spectra per posttreatment)

    read ( unit_post , * ) inp % corrspec_coord (:)

    read ( unit_post , * ) ! space
    read ( unit_post , * ) ! probe coordinates (name,x,y,z)

    ind = 1 ; loop = .true.
    do while (loop)

       read ( unit_post , * ) inp % probe_name (ind)

       if ( inp % probe_name (ind) == 'END' ) then

          inp % nprobes = ind-1
          loop = .false.

       else

          backspace (unit_post)
          read ( unit_post , * ) inp % probe_name (ind) , &
                                 inp % probe_coord (ind,:)

       end if

       ind = ind+1

    end do

    read ( unit_post , * ) ! space
    read ( unit_post , * ) ! names of the variables to plot

    ind = 1 ; loop = .true.
    do while (loop)

       read ( unit_post , * ) inp % var_name (ind)

       if ( inp % var_name (ind) == 'W_scal_i' ) then

          do i = 1 , nreac
             write ( ascii , format ) i
             inp % var_name (ind) = 'W_scal_i' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'Mass_Fractions' ) then

          do i = 1 , nrv+npv+nvv
             write ( ascii , format ) i
             inp % var_name (ind) = 'Y' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'Mass_Fractions_steady' ) then

          do i = 1 , inp % nss
             write ( ascii , format ) i
             inp % var_name (ind) = 'Yss' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'ra_Mass_Frac' ) then

          do i = 1 , nrv+npv+nvv
             write ( ascii , format ) i
             inp % var_name (ind) = 'ra_Y' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'fa_Mass_Frac' ) then

          do i = 1 , nrv+npv+nvv
             write ( ascii , format ) i
             inp % var_name (ind) = 'fa_Y' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'ra_Mass_Frac_steady' ) then

          do i = 1 , inp % nss
             write ( ascii , format ) i
             inp % var_name (ind) = 'ra_Yss' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'fa_Mass_Frac_steady' ) then

          do i = 1 , inp % nss
             write ( ascii , format ) i
             inp % var_name (ind) = 'fa_Yss' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'rv_Mass_Frac' ) then

          do i = 1 , nrv+npv+nvv
             write ( ascii , format ) i
             inp % var_name (ind) = 'rv_Y' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'fv_Mass_Frac' ) then

          do i = 1 , nrv+npv+nvv
             write ( ascii , format ) i
             inp % var_name (ind) = 'fv_Y' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'fv_Mass_Frac_steady' ) then

          do i = 1 , inp % nss
             write ( ascii , format ) i
             inp % var_name (ind) = 'fv_Yss' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'Productions' ) then ! Production rates loop

          do i = 1 , nrv
             write ( ascii , format ) i
             inp % var_name (ind) = 'W' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'Productions_SGS' ) then ! Production rates loop

          do i = 1 , nrv
             write ( ascii , format ) i
             inp % var_name (ind) = 'W' // trim (ascii) // '_SGS'
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'ra_Productions' ) then ! Production rates loop

          do i = 1 , nrv
             write ( ascii , format ) i
             inp % var_name (ind) = 'ra_W' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'fa_Productions' ) then ! Production rates loop

          do i = 1 , nrv
             write ( ascii , format ) i
             inp % var_name (ind) = 'fa_W' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'rv_Productions' ) then ! Production rates loop

          do i = 1 , nrv
             write ( ascii , format ) i
             inp % var_name (ind) = 'rv_W' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'fv_Productions' ) then ! Production rates loop

          call abort_mpi ( 'variable ''fv_Productions'' does not exist in module_stat_plot.f90' )
          do i = 1 , nrv
             write ( ascii , format ) i
             inp % var_name (ind) = 'fv_W' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'SSparameter' ) then ! Steady state parameter

          do i = 1 , nrv
             write ( ascii , format ) i
             inp % var_name (ind) = 'SS' // trim (ascii)
             ind = ind+1
          end do

          ind = ind-1

       else if ( inp % var_name (ind) == 'END' ) then

          inp % nvars = ind-1
          loop = .false.

       end if

       ind = ind+1

    end do

    close (unit_post)


    ! verifying the number of variables to plot
    if ( inp % nvars > nvarsmax ) &
       call abort_mpi ( 'the number of variables to be plotted exceeds ' // trim(str(nvarsmax)) )

    ! A little bit of posttreatment
    ind = ( inp % end_file - inp % start_file ) / inp % skip_file + 1
    inp % end_file = inp % start_file + ( ind - 1 ) * inp % skip_file ! corrected inp % end_file
    nss    = inp % nss
    nbins  = inp % nbins
    ntimes = inp % end_file

    ! Checking compatibilities between MPI division and plane_types
    if ( nproc == 1 ) then ! sequential problem

       if ( inp % plane_type == 'XY' .and. inp % nprocz /= 1 ) &
          call abort_mpi ('error: post-process plan XY with more than 1 proc in z-direction')

       if ( inp % plane_type == 'XZ' .and. inp % nprocy /= 1 ) &
          call abort_mpi ('error: post-process plan XZ with more than 1 proc in y-direction')

       if ( inp % plane_type == 'YZ' .and. inp % nprocx /= 1 ) &
          call abort_mpi ('error: post-process plan YZ with more than 1 proc in x-direction')

    else

       if ( inp % plane_type == 'XY' .and. dims(1) /= 1 ) &
          call abort_mpi ('error: post-process plan XY with more than 1 proc in z-direction')

       if ( inp % plane_type == 'XZ' .and. dims(2) /= 1 ) &
          call abort_mpi ('error: post-process plan XZ with more than 1 proc in y-direction')

       if ( inp % plane_type == 'YZ' .and. dims(3) /= 1 ) &
          call abort_mpi ('error: post-process plan YZ with more than 1 proc in x-direction')

    end if


    ! Adimensionalize
    Lref_i = 1.0_dp / adi % L_ref
    inp % corrspec_coord (:) = inp % corrspec_coord (:) * inp % dim_length_coord * Lref_i
    do i = 1 , inp % nprobes
       inp % probe_coord (i,:) = inp % probe_coord (i,:) * inp % dim_length_coord * Lref_i
    end do


  end subroutine readpost


end module input
