!------------------------------------------------------------------------------
! MODULE: input
!------------------------------------------------------------------------------
!> \brief Input module.
!!
!! This module contains all the operations related to parallel
!! communications through the Message Passing Interface (MPI) library.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module input


   use parameters
   use parallel
   use adim
   use parsing


   implicit none

   real (dp) , parameter , private :: pi = acos(-1.0_dp) !< pi definition


   integer (ip) , parameter , private :: nplanemax  = 20    , & !< maximum number of planes
   nstatmax   = 10000 , & !< maximum number of stat files
   nprobesmax = 20    , & !< maximum number of probes
   nvarsmax   = 200       !< maximum number of stat variables


   !> \brief external input datas (simulations, experiments, ...) derived type declaration.
   !!
   !!
   !> \author
   !! modeling, simulation and data analysis (msda)\n
   !! université mohammed vi polytechnique (um6p)\n
   !! benguerir, morocco

   type ext_type


   !> synthetic turbulence generator
   real (dp) , dimension (:,:,:)     , allocatable   :: uin , vin , win !< inlet velocity components
   real (dp) , dimension (:,:,:)     , allocatable   :: tke !< turbulent kinetic energy


   !> table: reversed chemical time scale in function of the mixture fraction
   !  tab features: rctmax = r_eversed c_hemical t_ime scale is max
   real (dp)                     :: rctmax , z_rctmax
   integer (ip)                  :: index_rctmax , nval

   real (dp) , dimension (:,:)       , allocatable :: table

   !> klein turbulence generators
   real (dp) , dimension (:,:,:)     , allocatable   :: b_ijk
   real (dp) , dimension (:,:,:,:)   , allocatable   :: rij  !< reynolds tensor
   real (dp) , dimension (:,:)       , allocatable   :: a11 , a22 , a33 , a21 , a31 , a32
   real (dp) , dimension (:,:,:)     , allocatable   :: randx , randy , randz

   !> thi initialisation
   real (dp) , dimension (:,:,:,:)   , allocatable   :: var_init


   end type ext_type


   !> \brief input derived type declaration.
   !!
   !!
   !> \author
   !! modeling, simulation and data analysis (msda)\n
   !! université mohammed vi polytechnique (um6p)\n
   !! benguerir, morocco

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

   !> skewness angle variables
   real (dp)             :: skewangle

   !> viscous variables
   logical               :: vis
   logical               :: vis_dt

   !> les variables
   logical                 :: les , tau_iso_sgs_switch
   character (len_default) :: mu_sgs_model , tau_iso_sgs_model
   real (dp)               :: mu_sgs_factor , tau_iso_sgs_factor , &
   pr_sgs , sc_sgs , cte_darlyharlow_sgs , &
   percent_weight_sgs

   !> reactive variables
   logical                        :: reaction
   character (len_default)        :: reac_selec
   logical                        :: dvodesimp
   integer (ip)                   :: index_h2 , index_n2 , index_o2 , index_h2o
   real (dp) , allocatable , dimension (:) :: y1f , y0o

   !> eglib variables
   logical               :: eglib
   integer (ip)          :: eglib_precision

   logical               :: forcesum

   !> stability criteria
   logical               :: fixeddt
   real (dp)             :: dt
   real (dp)             :: cfl , cflmin
   real (dp)             :: fo
   real (dp)             :: cpr
   real (dp)             :: dtmin_allowed
   logical               :: dtlimit

   !> weno variables
   logical :: weights
   logical :: opt
   integer (ip) :: optord
   real (dp) :: percent_weight
   real (dp)             :: ducros_treduc
   character (len_default) :: shock_sensor_type

   !> print iterations at the screen
   integer (ip)          :: itshow
   integer (ip)          :: itmax
   logical               :: itlim
   integer (ip)          :: walltime

   !> restart files variables
   logical                               :: read_restart
   character (len_default)               :: name_restart
   integer (ip)                          :: number_restart
   integer (ip)                          :: freq_restart
   real (dp)                             :: time_offset , initial_time , dtime

   !> stat files variables
   logical                               :: stat
   integer (ip)                          :: nstat
   integer (ip)                          :: nvolume , nxystat , nxzstat , nyzstat
   integer (ip) , dimension (nplanemax)  :: i_volmin , j_volmin , k_volmin
   integer (ip) , dimension (nplanemax)  :: i_volmax , j_volmax , kvol_max
   integer (ip) , dimension (nplanemax)  :: i_yzstat , j_xzstat , k_xystat
   real (dp)                             :: dim_length_coord
   real (dp) , dimension (nplanemax)     :: x_volmin(1) , y_volmin(1) , z_volmin(1)
   real (dp) , dimension (nplanemax)     :: x_volmax(1) , y_volmax(1) , z_volmax(1)
   real (dp) , allocatable , dimension(:):: x_yzstat , y_xzstat , z_xystat
   real (dp) , dimension (nstatmax)      :: timing

   !> post.dat file _additional_ variables
   integer (ip)                          :: nprocx , nprocy , nprocz
   integer (ip)                          :: ghostptx , ghostpty , ghostptz
   logical                               :: read_be , write_be
   logical                               :: trafo , nondim_grid
   real (dp)                             :: length_ref
   real (dp)                             :: umax
   real (dp)                             :: umin
   real (dp)                             :: uc
   real (dp)                             :: rhom
   real (dp)                             :: mum
   real (dp)                             :: yst

   logical                               :: dis_rate_therm

   integer (ip)                          :: nreac , index_fuel , index_oxidizer , index_da , index_scalar

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
   integer (ip)                          :: sx , ex , nx , sy , ey , ny , sz , ez , nz

   integer (ip)                                     :: nvars !< number of statistical variables
   character (len_default) , dimension (nvarsmax)   :: var_name !< names of the statistical variables

   real (dp) , dimension (3)                        :: corrspec_coord !< autocorrelation

   !> probes, used for ascii files, e.g. convergences, pdfs, spectras...
   integer (ip)                                     :: nprobes
   character (len_default) , allocatable , dimension (:) :: probe_name
   real (dp) , allocatable , dimension (:,:)             :: probe_coord


   !> input sub derived type
   type (ext_type)                       :: ext


   end type inp_type


   contains

   !> \brief Read input file
   !!
   !> \author
   !! Modeling, Simulation and Data Analysis (MSDA)\n
   !! Université Mohammed VI Polytechnique (UM6P)\n
   !! Benguerir, Morocco

   subroutine readinp (inp , pars)

      type (inp_type) , intent (inout) :: inp !< input derived type
      type (cfg_type) , intent (inout) :: pars !< parsing derived type

      call parsing_input (inp , pars)
      call assessing_input (inp)

   end subroutine readinp


   !> \brief Read input file
   !!
   !> \author
   !! Modeling, Simulation and Data Analysis (MSDA)\n
   !! Université Mohammed VI Polytechnique (UM6P)\n
   !! Benguerir, Morocco

   subroutine parsing_input (inp , pars)

      type (inp_type) , intent (inout) :: inp !< input derived type
      type (cfg_type) , intent (inout) :: pars !< parsing derived type

      ! local
      real (dp) :: stor_volume(6) 
      integer (ip) :: nb_probs_x , nb_probs_y , nb_probs_z , y1f_size , y0o_size 

      call CFG_add (pars , "FILENAME","shearlayer.ini", "DEFAUT SIMILATION PARAMETERS")

      call CFG_add (pars , "GENERAL%DIM" , 1 , "DIMENSION OF THE PROBLEM")
      call CFG_add (pars , "GENERAL%NRV" , 2 , "NUMBER OF REACTIVE VARIABLES")
      call CFG_add (pars , "GENERAL%NPV" , 0 , "NUMBER OF PASSIVE VARIABLES (0 OR 1)")
      call CFG_add (pars , "GENERAL%NVV" , 0 , "NUMBER OF VARIANCE VARIABLES USED FOR LES (0 OR 1)")
      call CFG_add (pars , "GENERAL%EGLIB" , .FALSE. , "EGLIB DETAILED TRANSPORT (TRUE/FALSE)")
      call CFG_add (pars , "GENERAL%EGLIB_PRECISION" , 1 , "EGLIB DETAILED TRANSPORT PRESISION (1,2,3)")
      call CFG_add (pars , "GENERAL%VIS" , .TRUE. , "VISCOUS PROBLEM (TRUE=NS/FALSE=EULER)")
      call CFG_add (pars , "GENERAL%INIT" , 'AMBIENT' , "INITIALIZATION")
      call CFG_add (pars , "GENERAL%FORCESUM" , .TRUE. , "FORCE SUM OF Y TO 1")

      call CFG_add (pars , "MONITORING%IND_FILES" , .FALSE. , "MANIPULATE INDIVIDUAL FILES (FALSE BY DEFAULT)")
      call CFG_add (pars , "MONITORING%INI_SOL" , .FALSE. , "STORING ONLY THE INITIAL SOLUTION")
      call CFG_add (pars , "MONITORING%FILTER" , .FALSE. , "FILTER SOLUTIONS AS A SPONGE ZONE")
      call CFG_add (pars , "MONITORING%FIL_XINI" , 0.5_dp , "FILTER SOLUTIONS START IN X-DIRECTION")
      call CFG_add (pars , "MONITORING%FIL_XEND" , 1.0_dp , "FILTER SOLUTIONS END IN X-DIRECTION")
      call CFG_add (pars , "MONITORING%ITSHOW" , 1 , "ITERATIONS TO SHOW INFORMATION")
      call CFG_add (pars , "MONITORING%WALLTIME" , 1000 , "WALLTIME IN MINUTES OF CALCULATION")
      call CFG_add (pars , "MONITORING%ITLIM" , .FALSE. , "ITERATIONS LIMIT (TRUE/FALSE)")
      call CFG_add (pars , "MONITORING%ITMAX" , 10 , "MAXIMUM ITERATIONS (VALUE)")
      call CFG_add (pars , "MONITORING%READ_RESTART" , .FALSE. , "READING A PREVIOUS RESTART FILE")
      call CFG_add (pars , "MONITORING%NAME_RESTART" , "RESTART" , "READING A PREVIOUS RESTART NAME")
      call CFG_add (pars , "MONITORING%NUMBER_RESTART" , 1 , "READING A PREVIOUS RESTART NUMBER")
      call CFG_add (pars , "MONITORING%FREQ_RESTART" , 1 , "ITERATIONS FREQUENCY TO STORE RESTART FILES")
      call CFG_add (pars , "MONITORING%INITIAL_TIME" , 0.0_dp , "INITIAL SIMULATION TIME")
      call CFG_add (pars , "MONITORING%TIME_OFFSET" , 1.0e-10_dp , "TIME TO START RECORDING STAT FILES (SECONDE)")
      call CFG_add (pars , "MONITORING%DTIME" , 1.0e-7_dp , "DT BETWEEN STAT FILES")
      call CFG_add (pars , "MONITORING%STAT" , .FALSE. , "STORE STAT FILES")
      call CFG_add (pars , "MONITORING%NSTAT" , 1000 , "NUMBER OF STATS")
      call CFG_add (pars , "MONITORING%VOLUME" , [0.0_dp , 0.0_dp , 0.0_dp , 1.0_dp , 1.0_dp , 1.0_dp ] , "VOLUME FILES (XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX COORDINATES)" , dynamic_size=.False.)
      call CFG_add (pars , "MONITORING%Z_XYSTAT" , [1.0_dp] , "STORING STAT FILES IN XY-PLANE (Z-COORDINATES)" , dynamic_size=.TRUE.)
      call CFG_add (pars , "MONITORING%Y_XZSTAT" , [1.0_dp] , "STORING STAT FILES IN XZ-PLANE (Y-COORDINATES)" , dynamic_size=.TRUE.)
      call CFG_add (pars , "MONITORING%X_YZSTAT" , [1.0_dp] , "STORING STAT FILES IN YZ-PLANE (X-COORDINATES)" , dynamic_size=.TRUE.)
      call CFG_add (pars , "MONITORING%PROBE_NAME" , ["PROBE1","PROBE2"] , "PROBE COORDINATES (NAME)" , dynamic_size=.TRUE.)
      call CFG_add (pars , "MONITORING%X-PROBE_COORD" , [0.0_dp , 1.0_dp] , "PROBE COORDINATES (X COORDINATES)" , dynamic_size=.TRUE.)
      call CFG_add (pars , "MONITORING%Y-PROBE_COORD" , [0.0_dp , 1.0_dp] , "PROBE COORDINATES (Y COORDINATES)" , dynamic_size=.TRUE.)
      call CFG_add (pars , "MONITORING%Z-PROBE_COORD" , [0.0_dp , 1.0_dp] , "PROBE COORDINATES (Z COORDINATES)" , dynamic_size=.TRUE.)
      call CFG_add (pars , "MONITORING%NVOLUME", 1 , "STORING VOLUME")

      call CFG_add (pars , "TIMESTEPPING%VIS_DT" , .TRUE. , "CONSIDER VISCOUS TIME STEP")
      call CFG_add (pars , "TIMESTEPPING%FIXEDDT" , .FALSE. , "FIXED TIME STEP (TRUE/FALSE)")
      call CFG_add (pars , "TIMESTEPPING%DT" , 1.0e-10_dp , "FIXED TIME STEP (VALUE IN SECONDE)")
      call CFG_add (pars , "TIMESTEPPING%CFLMIN" , 0.1_dp , "CONVECTIVE CFL (MIN)")
      call CFG_add (pars , "TIMESTEPPING%CFL" , 0.5_dp , "CONVECTIVE CFL (MIN)")
      call CFG_add (pars , "TIMESTEPPING%FO" , 0.1_dp , "FO (<=0.1) VISCOUS")
      call CFG_add (pars , "TIMESTEPPING%CPR" , 20.0_dp , "DTPERMIS (<=20K) CHEMICAL")
      call CFG_add (pars , "TIMESTEPPING%DTLIMIT" , .FALSE. , "IMPOSED TIME STEP MINIMUM VALUE (TRUE/FALSE)")
      call CFG_add (pars , "TIMESTEPPING%DTMIN_ALLOWED" , 1.0e-10_dp , "IMPOSED TIME STEP MINIMUM VALUE (VALUE-RECOMMAND 1E-13 MIN)")

      call CFG_add (pars , "IFLES%LES" , .FALSE. , "LES SIMULATION (TRUE/FALSE)")
      call CFG_add (pars , "IFLES%MU_SGS_MODEL" , "SMAGORINSKY" , "VISCOSITY SGS MODEL (SMAGORINSKY, KOLMOGOROV, WALE)")
      call CFG_add (pars , "IFLES%MU_SGS_FACTOR" , 0.18_dp , "VISCOSITY SGS MODEL (SMAGORINSKY 0.18, KOLMOGOROV 1.4, WALE 0.18)")
      call CFG_add (pars , "IFLES%TAU_ISO_SGS_SWITCH" , .FALSE. , "ISOTROPIC TENSOR SGS MODEL (SWITCH)")
      call CFG_add (pars , "IFLES%TAU_ISO_SGS_MODEL" , "YOSHIZAWA" , "ISOTROPIC TENSOR SGS MODEL (NAME)")
      call CFG_add (pars , "IFLES%TAU_ISO_SGS_FACTOR" , 1.0_dp , "ISOTROPIC TENSOR SGS MODEL (FACTOR)")
      call CFG_add (pars , "IFLES%CTE_DARLYHARLOW_SGS" , 0.08_dp , "DARLY & HARLOW CONSTANT")
      call CFG_add (pars , "IFLES%PERCENT_WEIGHT_SGS" , 70.0_dp , "WEIGHTS DENSITY PERCENT CRITERIA FOR DISCONTINUITIE ZONES")
      call CFG_add (pars , "IFLES%SC_SGS" , 0.7_dp , "SGS SCHMIDT (0.7 TO 1.0)")
      call CFG_add (pars , "IFLES%PR_SGS" , 0.7_dp , "SGS PRANDTL (0.7 TO 1.0) -> (0.3 TO 0.9)")

      call CFG_add (pars , "REACTING%REACTION" , .FALSE. , "ACTIVATE REACTION")
      call CFG_add (pars , "REACTING%REAC_SELEC" , "DVODE" , "REACTING SOLVER USED (DVODE, MIL)")
      call CFG_add (pars , "REACTING%DVODESIMP" , .FALSE. , "SIMPLIFY DVODE TO REDUCE COMPUTATIONAL COSTS")
      call CFG_add (pars , "REACTING%INDEX_N2" , 1 , "N2 INDEX (DVODE, IF SIMPLIFY=TRUE OR TO CALCULATE: DETECTOR_REACTION)")
      call CFG_add (pars , "REACTING%INDEX_H2" , 1 , "H2 INDEX (DVODE, IF SIMPLIFY=TRUE OR TO CALCULATE: DETECTOR_REACTION)")
      call CFG_add (pars , "REACTING%INDEX_O2" , 1 , "O2 INDEX (DVODE, IF SIMPLIFY=TRUE OR TO CALCULATE: DETECTOR_REACTION)")
      call CFG_add (pars , "REACTING%INDEX_H2O" , 1 , "H2O INDEX (DVODE, IF SIMPLIFY=TRUE OR TO CALCULATE: DETECTOR_REACTION)")
      call CFG_add (pars , "REACTING%YST" , 0.41_dp , "YST")
      call CFG_add (pars , "REACTING%DIS_RATE_THERM" , .TRUE. , "SAME SCHMIDT AND LEWIS EQUAL TO UNITY (FOR THE DIFFUSION COEFFICIENT)")
      call CFG_add (pars , "REACTING%NREAC" , 21 , "NUMBER OF ELEMENTARY REACTION STEPS")
      call CFG_add (pars , "REACTING%INDEX_FUEL" ,  2 , "INDEX FUEL (TAKENO)")
      call CFG_add (pars , "REACTING%INDEX_OXIDIZER" , 4 , "INDEX OXIDIZER (TAKENO)")
      call CFG_add (pars , "REACTING%INDEX_DA" , 1 , "INDEX SPECIES (DAMKOHLER")
      call CFG_add (pars , "REACTING%INDEX_SCALAR" , 1 , "SPECIES NUMBER FOR SCALAR PROPERTIES (0: MIXTURE FRACTION)")

      call CFG_add (pars , "WENO%WEIGHTS" , .TRUE. , "WENO NONLINEAR WEIGHTS (TRUE/FALSE)")
      call CFG_add (pars , "WENO%PERCENT_WEIGHT" , 15.0_dp , "WENO DENSITY PERCENT CRITERIA (VALUE)")
      call CFG_add (pars , "WENO%OPT" , .FALSE. , "WENO OPTIMUM WEIGHTS (TRUE/FALSE)")
      call CFG_add (pars , "WENO%OPTORD" , 1 , "WENO OPTIMUM ORDER (1, 3, 5)")
      call CFG_add (pars , "WENO%DUCROS_TREDUC" , 0.2_dp , "DUCROS SENSOR TRESHOLD (should be between 0 and 1")
      call CFG_add (pars , "WENO%SHOCK_SENSOR_TYPE" , 'ADAMSSHARIFF' , "DUCROS SENSOR SENSOR")



      call CFG_add (pars , "GRID%DIM_LENGTH_COORD" , 1.0_dp , "ADIMENSIONAL LENGTH FOR COORDINATES")
      call CFG_add (pars , "GRID%NONDIM_GRID" , .FALSE. , "NONDIMENSIONAL GRID")
      call CFG_add (pars , "GRID%LENGTH_REF" , 1.44E-04_dp , "REFERENCE LENGTH")
      call CFG_add (pars , "GRID%BCXMIN" , "PERIODIC" , "XMIN BOUNDARY CONDITION")
      call CFG_add (pars , "GRID%BCXMAX" , "PERIODIC" , "XMAX BOUNDARY CONDITION")
      call CFG_add (pars , "GRID%BCYMIN" , "PERIODIC" , "YMIN BOUNDARY CONDITION")
      call CFG_add (pars , "GRID%BCYMAX" , "PERIODIC" , "YMAX BOUNDARY CONDITION")
      call CFG_add (pars , "GRID%BCZMIN" , "PERIODIC" , "ZMIN BOUNDARY CONDITION")
      call CFG_add (pars , "GRID%BCZMAX" , "PERIODIC" , "ZMAX BOUNDARY CONDITION")
      call CFG_add (pars , "GRID%NUMPROCX" , 1 , "NUMBER OF PROCS IN X-DIRECTION")
      call CFG_add (pars , "GRID%NUMPROCY" , 1 , "NUMBER OF PROCS IN Y-DIRECTION")
      call CFG_add (pars , "GRID%NUMPROCZ" , 1 , "NUMBER OF PROCS IN Z-DIRECTION")

      call CFG_add (pars , "POSTPARALLEL%NPROCX" , 1 , "ACTUAL NUMBER OF PROCESSES IN X-DIRECTION")
      call CFG_add (pars , "POSTPARALLEL%NPROCY" , 1 , "ACTUAL NUMBER OF PROCESSES IN Y-DIRECTION")
      call CFG_add (pars , "POSTPARALLEL%NPROCZ" , 1 , "ACTUAL NUMBER OF PROCESSES IN Z-DIRECTION")
      call CFG_add (pars , "POSTPARALLEL%GHOSTPTX" , 0 , "ADDITIONAL PHANTOM POINTS FOR VOLUMES IN X-DIRECTION")
      call CFG_add (pars , "POSTPARALLEL%GHOSTPTY" , 0 , "ADDITIONAL PHANTOM POINTS FOR VOLUMES IN Y-DIRECTION")
      call CFG_add (pars , "POSTPARALLEL%GHOSTPTZ" , 0 , "ADDITIONAL PHANTOM POINTS FOR VOLUMES IN Z-DIRECTION")

      call CFG_add (pars , "DUMP%READ_BE" , .FALSE. , "READ BIG_ENDIAN BINARY FILES")
      call CFG_add (pars , "DUMP%WRITE_BE" , .FALSE. , "WRITE BIG_ENDIAN BINARY FILES")
      call CFG_add (pars , "DUMP%TRAFO" , .FALSE. , "TRANSFORM THE GRID")

      call CFG_add (pars , "POSTPRCESSING%NBINS" , 200 , "NUMBER OF BINS TO CREATE PDFS, CONDITIONAL AVERAGES...")
      call CFG_add (pars , "POSTPRCESSING%PLANE_TYPE" , 'RESTART' , "PLANE TYPE: XY, XZ, YZ, RESTART, VOLUME")
      call CFG_add (pars , "POSTPRCESSING%PLANE_NUMBER" , 0 , "PLANE TYPE: XY, XZ, YZ, RESTART, VOLUME")
      call CFG_add (pars , "POSTPRCESSING%TEMP_AVG" , .FALSE. , "TEMPORAL AVERAGE")
      call CFG_add (pars , "POSTPRCESSING%SPAT_AVG" , .FALSE. , "SPATIAL AVERAGE")
      call CFG_add (pars , "POSTPRCESSING%COND_AVG" , .FALSE. , "CONDITIONAL AVERAGE")
      call CFG_add (pars , "POSTPRCESSING%PDFS" , .FALSE. , "PDFS")
      call CFG_add (pars , "POSTPRCESSING%SPECTRAS" , .FALSE. , "SPECTRAS")
      call CFG_add (pars , "POSTPRCESSING%SIMILARITY_X" , .FALSE. , "SIMILARITY SOLUTIONS IN X-DIRECTION")
      call CFG_add (pars , "POSTPRCESSING%SECOND_LOOP" , .FALSE. , "SECOND LOOP")
      call CFG_add (pars , "POSTPRCESSING%READ_STAT" , .FALSE. , "READ PREVIOUS STAT")
      call CFG_add (pars , "POSTPRCESSING%START_FILE" , 1 , "START FILE")
      call CFG_add (pars , "POSTPRCESSING%END_FILE" , 1 , "START FILE")
      call CFG_add (pars , "POSTPRCESSING%SKIP_FILE" , 1 , "SKIP FILE")

      call CFG_add (pars , "SHEARLAYER%Y1F" , [0.0_dp , 0.05_dp , 0.0_dp , 0.0_dp , 0.0_dp , 0.0_dp , 0.0_dp , 0.0_dp , 0.95_dp] , "PURE FUEL COMPOSITION" , dynamic_size=.TRUE.)
      call CFG_add (pars , "SHEARLAYER%Y0O" , [5.6e-7_dp , 0.0_dp , 1.55e-4_dp , 0.278_dp , 1.83e-3_dp , 0.17_dp , 5.1e-6_dp , 2.5e-7_dp , 0.55000_dp] , "PURE OXYDIZER COMPOSITION" , dynamic_size=.TRUE.)
      call CFG_add (pars , "SHEARLAYER%PERTURB" , .FALSE. , "ACTIVATE INLET PERTURBATION FOR SHEAR LAYER")
      call CFG_add (pars , "SHEARLAYER%SKEWANGLE" , 0.0_dp , "SKENWESS ANGLE FOR SHEAR LAYERS")
      call CFG_add (pars , "SHEARLAYER%UMAX" , 1634.0_dp , "UMAX FOR SHEARLAYER")
      call CFG_add (pars , "SHEARLAYER%UMIN" , 973.3_dp , "UMIN FOR SHEARLAYER")
      call CFG_add (pars , "SHEARLAYER%UC" , 1208.18_dp , "UC FOR SHEARLAYER")
      call CFG_add (pars , "SHEARLAYER%RHOM" , 0.2785_dp , "RHO AVG FOR SHEARLAYER")
      call CFG_add (pars , "SHEARLAYER%MUM" , 4.1455e-5_dp , "MU AVG FOR SHEARLAYER")
      ! call CFG_add (pars , "POSTPRCESSING%VAR_NAME" , ['DENSITY', 'VELOCITY-X' , 'PRESSURE' , 'TEMPERATURE'] , 'VARIABLES TO POSTPROCESS', dynamic_size=.TRUE.)

      if (rank == rank_default) then
         ! write (*,*)
         ! write (*,*)
         write (*,*) '==================================================='
         ! write (*,*)
      endif

      call CFG_sort ( pars )

      call CFG_update_from_arguments(pars)

      call CFG_get (pars , "GENERAL%DIM" , inp % dim)
      call CFG_get (pars , "GENERAL%NRV" , inp % nrv)
      call CFG_get (pars , "GENERAL%NPV" , inp % npv)
      call CFG_get (pars , "GENERAL%NVV" , inp % nvv)
      call CFG_get (pars , "GENERAL%VIS" , inp % vis)
      call CFG_get (pars , "GENERAL%EGLIB" , inp % eglib)
      call CFG_get (pars , "GENERAL%EGLIB_PRECISION" , inp % eglib_precision)
      call CFG_get (pars , "GENERAL%INIT" , inp % init)
      call CFG_get (pars , "GENERAL%FORCESUM", inp % forcesum)

      call CFG_get (pars , "IFLES%LES" , inp % les)
      call CFG_get (pars , "IFLES%MU_SGS_MODEL" , inp % mu_sgs_model)
      call CFG_get (pars , "IFLES%MU_SGS_FACTOR" , inp % mu_sgs_factor)
      call CFG_get (pars , "IFLES%TAU_ISO_SGS_SWITCH" , inp % tau_iso_sgs_switch)
      call CFG_get (pars , "IFLES%TAU_ISO_SGS_MODEL" , inp % tau_iso_sgs_model)
      call CFG_get (pars , "IFLES%TAU_ISO_SGS_FACTOR" , inp % tau_iso_sgs_factor)
      call CFG_get (pars , "IFLES%CTE_DARLYHARLOW_SGS" , inp % cte_darlyharlow_sgs)
      call CFG_get (pars , "IFLES%PERCENT_WEIGHT_SGS" , inp % percent_weight_sgs)
      call CFG_get (pars , "IFLES%SC_SGS" , inp % sc_sgs)
      call CFG_get (pars , "IFLES%PR_SGS" , inp % pr_sgs)

      call CFG_get (pars , "REACTING%REACTION" , inp % reaction)
      call CFG_get (pars , "REACTING%REAC_SELEC" , inp % reac_selec)
      call CFG_get (pars , "REACTING%DVODESIMP" , inp % dvodesimp)
      call CFG_get (pars , "REACTING%INDEX_H2", inp % index_h2)
      call CFG_get (pars , "REACTING%INDEX_N2", inp % index_n2)
      call CFG_get (pars , "REACTING%INDEX_O2", inp % index_o2)
      call CFG_get (pars , "REACTING%INDEX_H2O", inp % index_h2o)
      call CFG_get (pars , "REACTING%YST" , inp % yst)
      call CFG_get (pars , "REACTING%DIS_RATE_THERM" , inp % dis_rate_therm)
      call CFG_get (pars , "REACTING%NREAC" , inp % nreac)
      call CFG_get (pars , "REACTING%INDEX_FUEL" , inp % index_fuel)
      call CFG_get (pars , "REACTING%INDEX_OXIDIZER" , inp % index_oxidizer)
      call CFG_get (pars , "REACTING%INDEX_DA" , inp % index_da)
      call CFG_get (pars , "REACTING%INDEX_SCALAR" , inp % index_scalar)

      call CFG_get (pars , "TIMESTEPPING%VIS_DT" , inp % vis_dt)
      call CFG_get (pars , "TIMESTEPPING%FIXEDDT", inp % fixeddt)
      call CFG_get (pars , "TIMESTEPPING%DT", inp % dt)
      call CFG_get (pars , "TIMESTEPPING%CFLMIN", inp % cflmin)
      call CFG_get (pars , "TIMESTEPPING%CFL", inp % cfl)
      call CFG_get (pars , "TIMESTEPPING%FO", inp % fo)
      call CFG_get (pars , "TIMESTEPPING%CPR", inp % cpr)
      call CFG_get (pars , "TIMESTEPPING%DTLIMIT", inp % dtlimit)
      call CFG_get (pars , "TIMESTEPPING%DTMIN_ALLOWED", inp % dtmin_allowed)

      call CFG_get (pars , "WENO%WEIGHTS", inp % weights)
      call CFG_get (pars , "WENO%PERCENT_WEIGHT", inp % percent_weight)
      call CFG_get (pars , "WENO%OPT", inp % opt)
      call CFG_get (pars , "WENO%OPTORD", inp % optord)
      call CFG_get (pars , "WENO%DUCROS_TREDUC" , inp % ducros_treduc)
      call CFG_get (pars , "WENO%SHOCK_SENSOR_TYPE" , inp % shock_sensor_type)


      call CFG_get (pars , "MONITORING%IND_FILES" , inp % ind_files)
      call CFG_get (pars , "MONITORING%INI_SOL" , inp % ini_sol)
      call CFG_get (pars , "MONITORING%FILTER" , inp % filter)
      call CFG_get (pars , "MONITORING%FIL_XINI" , inp % fil_xini)
      call CFG_get (pars , "MONITORING%FIL_XEND" , inp % fil_xend)
      call CFG_get (pars , "MONITORING%WALLTIME", inp % walltime)
      call CFG_get (pars , "MONITORING%ITLIM", inp % itlim)
      call CFG_get (pars , "MONITORING%ITMAX", inp % itmax)
      call CFG_get (pars , "MONITORING%ITSHOW", inp % itshow)
      call CFG_get (pars , "MONITORING%READ_RESTART", inp % read_restart)
      call CFG_get (pars , "MONITORING%NAME_RESTART", inp % name_restart)
      call CFG_get (pars , "MONITORING%NUMBER_RESTART", inp % number_restart)
      call CFG_get (pars , "MONITORING%FREQ_RESTART", inp % freq_restart)
      call CFG_get (pars , "MONITORING%INITIAL_TIME", inp % initial_time)
      call CFG_get (pars , "MONITORING%TIME_OFFSET", inp % time_offset)
      call CFG_get (pars , "MONITORING%DTIME", inp % dtime)
      call CFG_get (pars , "MONITORING%STAT", inp % stat)
      call CFG_get (pars , "MONITORING%NSTAT", inp % nstat)
      call CFG_get (pars , "MONITORING%NVOLUME", inp % nvolume)
      call CFG_get (pars , "MONITORING%VOLUME", stor_volume)
      call CFG_get_size(pars, "MONITORING%Z_XYSTAT", inp % nxystat)
      call CFG_get_size(pars, "MONITORING%Y_XZSTAT", inp % nxzstat)
      call CFG_get_size(pars, "MONITORING%X_YZSTAT", inp % nyzstat)

      call CFG_get (pars , "GRID%DIM_LENGTH_COORD", inp % dim_length_coord)
      call CFG_get (pars , "GRID%NONDIM_GRID" , inp % nondim_grid)
      call CFG_get (pars , "GRID%LENGTH_REF" , inp % length_ref)
      call CFG_get (pars , "GRID%BCXMIN" , bc(1))
      call CFG_get (pars , "GRID%BCXMAX" , bc(2))
      call CFG_get (pars , "GRID%BCYMIN" , bc(3))
      call CFG_get (pars , "GRID%BCYMAX" , bc(4))
      call CFG_get (pars , "GRID%BCZMIN" , bc(5))
      call CFG_get (pars , "GRID%BCZMAX" , bc(6))
      call CFG_get (pars , "GRID%NUMPROCX" , dims(3))
      call CFG_get (pars , "GRID%NUMPROCY" , dims(2))
      call CFG_get (pars , "GRID%NUMPROCZ" , dims(1))

      if ( inp % nvolume > 1 ) then
         write (*,*) 'no need to store more than one volume'
         write (*,*) 'please reconsider this parameter'
         call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
      end if

      if ( inp % stat .and. inp % dim == 3 .and. &
         & (inp % nvolume == 0 .and. inp % nxystat == 0 .and. &
            & inp % nxzstat == 0 .and. inp % nyzstat == 0 ) ) then
         write (*,*) 'error: you try to save statistic files for 3D simulation &
         &            without defining volumes or planes'
         write (*,*) 'please reconsider this parameter'
         call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
      end if

      inp % X_VOLMIN (1) = stor_volume(1)
      inp % X_VOLMAX (1) = stor_volume(2)
      inp % Y_VOLMIN (1) = stor_volume(3)
      inp % Y_VOLMAX (1) = stor_volume(4)
      inp % Z_VOLMIN (1) = stor_volume(5)
      inp % Z_VOLMAX (1) = stor_volume(6)

      allocate (inp % z_xystat(inp % nxystat))
      allocate (inp % y_xzstat(inp % nxzstat))
      allocate (inp % x_yzstat(inp % nyzstat))
      call CFG_get (pars , "MONITORING%X_YZSTAT", inp % x_yzstat)
      call CFG_get (pars , "MONITORING%Z_XYSTAT", inp % z_xystat)
      call CFG_get (pars , "MONITORING%Y_XZSTAT", inp % y_xzstat)


      call CFG_get_size(pars, "MONITORING%PROBE_NAME", inp % nprobes)
      call CFG_get_size(pars, "MONITORING%X-PROBE_COORD", nb_probs_x)
      call CFG_get_size(pars, "MONITORING%X-PROBE_COORD", nb_probs_y)
      call CFG_get_size(pars, "MONITORING%X-PROBE_COORD", nb_probs_z)

      if (nb_probs_x /= inp % nprobes .or.  nb_probs_y /= inp % nprobes .or. nb_probs_z /= inp % nprobes) then
         write (*,*) 'error: x-y or z COORDINATES should have the same size as probes_names'
         write (*,*) 'please reconsider this parameter'
         call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode)
      end if
      allocate (inp % probe_name(inp % nprobes))
      allocate (inp % probe_coord(inp % nprobes,3))
      call CFG_get (pars , "MONITORING%PROBE_NAME", inp % probe_name)
      call CFG_get (pars , "MONITORING%X-PROBE_COORD", inp % probe_coord(:,1))
      call CFG_get (pars , "MONITORING%Y-PROBE_COORD", inp % probe_coord(:,2))
      call CFG_get (pars , "MONITORING%Z-PROBE_COORD", inp % probe_coord(:,3))

      call CFG_get_size(pars, "SHEARLAYER%Y1F", y1f_size)
      call CFG_get_size(pars, "SHEARLAYER%Y0O", y0o_size)
      ! if (y1f_size /= (inp % nrv + inp % npv) .or. y0o_size /= inp % nrv + inp % npv) then
      !    write (*,*) 'error: Y1F or Y0O should have the same size as nrv + npv'
      !    write (*,*) 'please reconsider this parameter'
      !    call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode)
      ! end if
      allocate (inp % y1f(y1f_size))
      allocate (inp % y0o(y0o_size))
      call CFG_get (pars , "SHEARLAYER%Y1F", inp % y1f)
      call CFG_get (pars , "SHEARLAYER%Y0O", inp % y0o)
      call CFG_get (pars , "SHEARLAYER%PERTURB" , inp % perturb)
      call CFG_get (pars , "SHEARLAYER%SKEWANGLE", inp % skewangle)
      call CFG_get (pars , "SHEARLAYER%UMAX" , inp % umax)
      call CFG_get (pars , "SHEARLAYER%UMIN" , inp % umin)
      call CFG_get (pars , "SHEARLAYER%UC" , inp % uc)
      call CFG_get (pars , "SHEARLAYER%RHOM" , inp % rhom)
      call CFG_get (pars , "SHEARLAYER%MUM" , inp % mum)

      call CFG_get (pars , "POSTPARALLEL%NPROCX" , inp % nprocx)
      call CFG_get (pars , "POSTPARALLEL%NPROCY" , inp % nprocy)
      call CFG_get (pars , "POSTPARALLEL%NPROCZ" , inp % nprocz)
      call CFG_get (pars , "POSTPARALLEL%GHOSTPTX" , inp % ghostptx)
      call CFG_get (pars , "POSTPARALLEL%GHOSTPTY" , inp % ghostpty)
      call CFG_get (pars , "POSTPARALLEL%GHOSTPTZ" , inp % ghostptz)

      call CFG_get (pars , "DUMP%READ_BE" , inp % read_be)
      call CFG_get (pars , "DUMP%WRITE_BE" , inp % write_be)
      call CFG_get (pars , "DUMP%TRAFO" , inp % trafo)

      call CFG_get (pars , "POSTPRCESSING%NBINS" , inp % nbins)
      call CFG_get (pars , "POSTPRCESSING%PLANE_TYPE" , inp % plane_type)
      call CFG_get (pars , "POSTPRCESSING%PLANE_NUMBER" , inp % plane_number)
      call CFG_get (pars , "POSTPRCESSING%TEMP_AVG" , inp % temp_avg)
      call CFG_get (pars , "POSTPRCESSING%SPAT_AVG" , inp % spat_avg)
      call CFG_get (pars , "POSTPRCESSING%COND_AVG" , inp % cond_avg)
      call CFG_get (pars , "POSTPRCESSING%PDFS" , inp % pdfs)
      call CFG_get (pars , "POSTPRCESSING%SPECTRAS" , inp % spectras)
      call CFG_get (pars , "POSTPRCESSING%SIMILARITY_X" , inp % similarity_x)
      call CFG_get (pars , "POSTPRCESSING%SECOND_LOOP" , inp % second_loop)
      call CFG_get (pars , "POSTPRCESSING%READ_STAT" , inp % read_stat)
      call CFG_get (pars , "POSTPRCESSING%START_FILE" , inp % start_file)
      call CFG_get (pars , "POSTPRCESSING%END_FILE" , inp % end_file)
      call CFG_get (pars , "POSTPRCESSING%SKIP_FILE" , inp % skip_file)
      ! call CFG_get_size(pars, "MONITORING%PROBE_NAME", inp % nvars)
      ! allocate (inp % var_name (inp % nvars))
      ! call CFG_get (pars , "POSTPRCESSING%VAR_NAME" , inp % var_name)




      if (rank == rank_default) then
       call CFG_write(pars, "stdout", custom_first=.TRUE.) 
       call CFG_write(pars, "output.cfg")          
       call CFG_write_markdown(pars, "output.md") 
    endif

    if (rank == rank_default) then
      ! write (*,*)
      write (*,*) '==================================================='
   endif


end subroutine parsing_input


!> \brief Assess input file
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine assessing_input (inp)

   type (inp_type) , intent (inout) :: inp !< input derived type

   integer (ip) :: l


   if ( npv > 1 ) then
      write (*,*) 'the code is not prepared to support more than one passive scalar'
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
      nderiv = inp % dim * ( inp % dim + 2 ) + inp % dim * ( inp % nrv + inp % npv )
   else
      nderiv = inp % dim * ( inp % dim + 1 ) + inp % dim * ( inp % nrv + inp % npv )
   end if

   nderivSGS = 0
   if ( inp % LES ) then
      nderiv = nderiv +     &
      &        ndim * nvv + & ! variance variable (for LES only)
      &        ndim +       & ! W_i derivatives
      &        1              ! shock detector

      nderivSGS = inp % dim ! isotropic tensor derivatives
   end if

   ndim = inp % dim ! dimension of the problem
   nrv  = inp % nrv ! reactive variables
   npv  = inp % npv ! passive variables
   nvv  = inp % nvv ! variance variables
   nv   = niv + nrv + npv + nvv

   ind_files = inp % ind_files
   weno_avg        = inp % weights
   vis             = inp % vis
   vis_dt          = inp % vis_dt
   reaction        = inp % reaction
   dvodesimp       = inp % dvodesimp
   eglib           = inp % eglib
   eglib_precision = inp % eglib_precision
   LES             = inp % LES
   forcesum        = inp % forcesum

   CFL    = inp % CFL
   CFLmin = inp % CFLmin
   Fo     = inp % Fo
   CPR    = inp % CPR

   inp % percent_weight = inp % percent_weight
   max_rel_weight = inp % percent_weight
   ducros_treduc = inp % ducros_treduc
   percent_weight_SGS   = inp % percent_weight_SGS
   shock_sensor_type = inp % shock_sensor_type
   inp % walltime = inp % walltime * 60 ! minutes -> seconds

   if ( inp % read_restart) then
      inp % timing (1:inp % number_restart+1) = inp % time_offset
      do l = inp % number_restart+2 , inp % nstat
         inp % timing (l) = inp % timing (l-1) + inp % dtime
      end do
   else
      inp % timing (1) = inp % time_offset
      do l = 2 , inp % nstat
         inp % timing (l) = inp % timing (l-1) + inp % dtime
      end do
   end if

   if ( inp % nvolume > nstatmax .or. &
      & inp % nxystat > nstatmax .or. &
      & inp % nxzstat > nstatmax .or. &
      & inp % nyzstat > nstatmax ) then
      write (*,*) 'maximum number of stats reached'
      call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
   end if

   nreac = inp % nreac

   ! ! verifying the number of variables to plot
   ! if ( inp % nvars > nvarsmax ) then
   !    write (*,*) 'the number of variables to be plotted exceeds ' , nvarsmax
   !    call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
   ! end if

   ! A little bit of posttreatment
   nreac  = inp % nreac
   nbins  = inp % nbins
   ntimes = inp % end_file

   if ( inp % plane_type == 'XY' .and. dims(1) /= 1 ) then
      write (*,*) 'error: post-process plan XY with more than 1 proc in z-direction'
      call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
   end if

   if ( inp % plane_type == 'XZ' .and. dims(2) /= 1 ) then
      write (*,*) 'error: post-process plan XZ with more than 1 proc in y-direction'
      call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
   end if

   if ( inp % plane_type == 'YZ' .and. dims(3) /= 1 ) then
      write (*,*) 'error: post-process plan YZ with more than 1 proc in x-direction'
      call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
   end if

end subroutine assessing_input


!> \brief Read the file post.da
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine readpost (inp)


   type (inp_type) , intent (inout) :: inp !< input derived type


   integer (ip) :: ok , i , ind
   logical      :: loop

   character (len_default) , parameter :: format = ' ( I2.2 ) '
   character (len_default)             :: ascii


   open ( unit = unit_post , file = file_post , status = 'old' , iostat = ok )
   if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (file_post))

   ind = 1 ; loop = .TRUE.
   do while (loop)

      read ( unit_post , * ) inp % var_name (ind)

      if ( inp % var_name (ind) == 'W_scal_i' ) then

         do i = 1 , inp % nreac
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

      else if ( inp % var_name (ind) == 'Productions' ) then ! Production rates loop

         do i = 1 , nrv
            write ( ascii , format ) i
            inp % var_name (ind) = 'W' // trim (ascii)
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

         do i = 1 , nrv
            write ( ascii , format ) i
            inp % var_name (ind) = 'fv_W' // trim (ascii)
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

end subroutine readpost



! !> \brief Read the file input.dat.
! !!
! !!
! !> \author
! !! Modeling, Simulation and Data Analysis (MSDA)\n
! !! Université Mohammed VI Polytechnique (UM6P)\n
! !! Benguerir, Morocco

!   subroutine readinp (inp)


!     type (inp_type) , intent (inout) :: inp !< input derived type


!     integer (ip) :: ok , l , ind
!     logical      :: loop


!     open ( unit = unit_inp , file = file_inp , status = 'old' , iostat = ok )
!     if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (file_inp))

!     read ( unit_inp , * ) inp % dim                            !dimension of the problem
!     read ( unit_inp , * ) inp % nrv                            !reactive variables (= nb of SPECIES in chem.inp)
!     read ( unit_inp , * ) inp % npv                            !passive variables (= 0 or 1)
!     ! if (rank == rank_default) print*, 'npv is', inp % npv  
!     read ( unit_inp , * ) inp % nvv                            !variance variables for LES (= 0 or 1)
!     read ( unit_inp , * ) inp % ind_files                      !manipulate individual files (false by default)
!     read ( unit_inp , * ) inp % ini_sol                        !storing only the initial solution
!     read ( unit_inp , * ) inp % init                           !initialization (Diffusion, Bogey, Fu, Miller, MIXChengMiller, Premix, Jet, DifFChengMiller, Fedkiw, DMR, Ambient)
!     read ( unit_inp , * ) inp % perturb                        !perturbation
!     read ( unit_inp , * ) inp % filter                         !filter
!     read ( unit_inp , * ) inp % fil_xini , inp % fil_xend      !   filter_start, filter_end
!     read ( unit_inp , * ) inp % vis                            !viscous problem
!     read ( unit_inp , * ) inp % vis_dt                         !viscous time step (true by defaulut)
!     read ( unit_inp , * ) inp % eglib , inp % eglib_precision  !eglib detailed transport (true/false), presision (1,2,3)
!     read ( unit_inp , * ) inp % LES                            !LES simulation
!     read ( unit_inp , * ) inp % mu_SGS_model , inp % mu_SGS_factor                                      !   Viscosity SGS model (Smagorinsky 0.18, Kolmogorov 1.4)
!     read ( unit_inp , * ) inp % tau_iso_SGS_switch , inp % tau_iso_SGS_model , inp % tau_iso_SGS_factor !   Isotropic tensor SGS model (switch, name, factor)
!     read ( unit_inp , * ) inp % Cte_DarlyHarlow_SGS            !      Darly & Harlow constant
!     read ( unit_inp , * ) inp % percent_weight_SGS             !      weights density percent criteria for discontinuitie zones
!     read ( unit_inp , * ) inp % Sc_SGS                         !   SGS Schmidt (0.7 to 1.0)
!     read ( unit_inp , * ) inp % Pr_SGS                         !   SGS Prandtl (0.7 to 1.0) -> (0.3 to 0.9)
!     read ( unit_inp , * ) inp % reaction , inp % reac_selec    !reactive problem (true/false), solver (dvode, MIL)
!     read ( unit_inp , * ) inp % dvodesimp                      !   simplify dvode to reduce computational costs
!     read ( unit_inp , * ) inp % Y1f (1:inp % nrv + inp % npv)  !   pure fuel composition (dvode, if simplify=true or to calculate: Detector_reaction, Zm, Deq)
!     read ( unit_inp , * ) inp % Y0o (1:inp % nrv + inp % npv)  !   pure oxydizer composition (idem)
!     read ( unit_inp , * ) inp % index_H2                       !   H2 index (dvode, if simplify=true or to calculate: Detector_reaction)
!     read ( unit_inp , * ) inp % index_N2                       !   N2 index (idem)
!     read ( unit_inp , * ) inp % index_O2                       !   O2 index (MIL)
!     read ( unit_inp , * ) inp % index_H2O                      !   H2O index (MIL)
!     read ( unit_inp , * ) inp % forcesum                       !force sumY=1 ! (true by default)
!     read ( unit_inp , * ) inp % fixeddt , inp % dt             !   fixed time step (true/false, value in seconde)
!     read ( unit_inp , * ) inp % CFLmin , inp % CFL             !   CFL (min,max<=0.9) convective
!     read ( unit_inp , * ) inp % Fo                             !   Fo (<=0.1) viscous
!     read ( unit_inp , * ) inp % CPR                            !   DTpermis (<=20K) chemical
!     read ( unit_inp , * ) inp % dtlimit , inp % dtmin_allowed  !   imposed time step minimum value (true/false, value) (recommand 1e-13 min)
!     read ( unit_inp , * ) inp % weights , inp % percent_weight !WENO nonlinear weights (true/false), density percent criteria (value)
!     read ( unit_inp , * ) inp % opt , inp % optord             !WENO optimum weights (true/false), order (1, 3, 5)
!     read ( unit_inp , * ) inp % itshow                         !iterations to show information
!     read ( unit_inp , * ) inp % walltime                       !walltime in minutes of calculation
!     read ( unit_inp , * ) inp % itlim , inp % itmax            !iterations limit (true/false), maximum iterations (value)
!     read ( unit_inp , * ) inp % read_restart , inp % name_restart , inp % number_restart !reading a previous restart file, name, number
!     read ( unit_inp , * ) inp % freq_restart                   !iterations frequency to store restart files
!     read ( unit_inp , * ) inp % initial_time                   !initial simulation time
!     read ( unit_inp , * ) inp % time_offset                    !time to start recording stat files (seconde)
!     read ( unit_inp , * ) inp % dtime                          !   dt between stat files
!     read ( unit_inp , * ) inp % stat                           !   store stat files
!     read ( unit_inp , * ) inp % nstat                          !   number of stats
!     read ( unit_inp , * ) inp % nvolume                        !number of volumes
!     read ( unit_inp , * ) inp % nxystat                        !number of XY-planes
!     read ( unit_inp , * ) inp % nxzstat                        !number of XZ-planes
!     read ( unit_inp , * ) inp % nyzstat                        !number of YZ-planes
!     read ( unit_inp , * ) inp % dim_length_coord               !adimensional length for coordinates
!     read ( unit_inp , * ) inp % skewangle                      !skewness angle of SL
!     ! if (rank == rank_default) print*, 'the angle is',inp % skewangle 

!     read ( unit_inp , * ) !
!     read ( unit_inp , * ) ! Storing volume files (xmin,xmax,ymin,ymax,zmin,zmax coordinates)

!     do l = 1 , inp % nvolume
!        read ( unit_inp , * ) inp % x_volmin (l) , inp % x_volmax (l) , &
!                              inp % y_volmin (l) , inp % y_volmax (l) , &
!                              inp % z_volmin (l) , inp % z_volmax (l)
!     end do

!     read ( unit_inp , * ) !
!     read ( unit_inp , * ) ! Storing stat files in XY-plane (Z-coordinates)

!     do l = 1 , inp % nxystat
!        read ( unit_inp , * ) inp % z_xystat (l)
!     end do

!     read ( unit_inp , * ) !
!     read ( unit_inp , * ) ! Storing stat files in XZ-plane (Y-coordinates)

!     do l = 1 , inp % nxzstat
!        read ( unit_inp , * ) inp % y_xzstat (l)
!     end do

!     read ( unit_inp , * ) !
!     read ( unit_inp , * ) ! Storing stat files in YZ-plane (X-coordinates)

!     do l = 1 , inp % nyzstat
!        read ( unit_inp , * ) inp % x_yzstat (l)
!     end do

!     read ( unit_inp , * ) !
!     read ( unit_inp , * ) ! Probe coordinates (name,x,y,z)

!     ind = 1 ; loop = .TRUE.
!     inp % probe_coord = 0.0_dp
!     do while (loop)

!        read ( unit_inp , * ) inp % probe_name (ind)

!        if ( inp % probe_name (ind) == 'END' ) then

!           inp % nprobes = ind-1
!           loop = .false.

!        else

!           backspace (unit_inp)
!           read ( unit_inp , * ) inp % probe_name (ind) , &
!                                 inp % probe_coord (ind,:)

!        end if

!        ind = ind+1

!     end do

!     close (unit_inp)


!     ! A little bit of post-treatment


!     if ( npv > 1 ) then
!        write (*,*) 'the code is not prepared to support more than one passive scalar'
!        write (*,*) 'please reconsider this parameter'
!        call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
!     end if

!     !if ( inp % stat .and. inp % dim == 3 .and. &
!     !     ( inp % nvolume == 0 .and. inp % nxystat == 0 .and. inp % nxzstat == 0 .and. inp % nyzstat == 0 ) ) then
!     !   write (*,*) 'error: you try to save statistic files (volumes or planes) for 3D simulation without defining volumes or planes'
!     !   write (*,*) 'please reconsider this parameter'
!     !   call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
!     !end if

! !    if ( nvv > 1 ) then
! !       write (*,*) 'the code is not prepared to support more than one variance scalar'
! !       write (*,*) 'please reconsider this parameter'
! !       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
! !    end if

! !    if ( nvv >= 1 .and. ( nvv /= npv ) ) then
! !       write (*,*) 'error: nvv /= npv, you need as many npv (passive scalare) as nvv (variance)'
! !       write (*,*) 'please reconsider this parameter'
! !       ! if you want more than 1 nvv, you need to arrange the vector of conservative variable like :
! !       ! rho, rhou, rhov, rhow, rhoet, rhoYa, rhoPs1, rhoPs2, rhoPsn, rhoVar1, rhoVar2, ..., rhoVarn
! !       call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
! !    end if

!     if ( nvv >= 1 .and. ( .not. inp % LES ) ) then
!        write (*,*) 'error: nvv >=1 and LES is desactivated, you need LES'
!        write (*,*) 'please reconsider this parameter'
!        call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
!     end if

!     if ( inp % tau_iso_SGS_switch .and. ( .not. inp % LES ) ) then
!        write (*,*) 'error: input.dat -> calculate SGS isotropic tensor without activate LES'
!        write (*,*) 'please reconsider this parameter'
!        call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
!     end if

!     if ( inp % eglib ) then
!        nderiv = inp % dim * ( inp % dim + 2 ) +       &
!                 inp % dim * ( inp % nrv + inp % npv )
!     else
!        nderiv = inp % dim * ( inp % dim + 1 ) +       &
!                 inp % dim * ( inp % nrv + inp % npv )
!     end if

!     nderivSGS = 0
!     if ( inp % LES ) then
!        nderiv = nderiv +     &
!                 ndim * nvv + & ! variance variable (for LES only)
!                 ndim +       & ! W_i derivatives
!                 1              ! shock detector

!        nderivSGS = inp % dim ! isotropic tensor derivatives
!     end if

!     ind_files = inp % ind_files

!     weno_avg        = inp % weights
!     vis             = inp % vis
!     vis_dt          = inp % vis_dt
!     reaction        = inp % reaction
!     dvodesimp       = inp % dvodesimp
!     eglib           = inp % eglib
!     eglib_precision = inp % eglib_precision
!     LES             = inp % LES
!     forcesum        = inp % forcesum

!     CFL    = inp % CFL
!     CFLmin = inp % CFLmin
!     Fo     = inp % Fo
!     CPR    = inp % CPR

!     inp % percent_weight = inp % percent_weight / 100.0_dp
!     max_rel_weight       = inp % percent_weight

!     percent_weight_SGS   = inp % percent_weight_SGS / 100.0_dp

!     inp % walltime = inp % walltime * 60 ! minutes -> seconds

!     if ( inp % read_restart) then
!        inp % timing (1:inp % number_restart+1) = inp % time_offset
!        do l = inp % number_restart+2 , inp % nstat
!           inp % timing (l) = inp % timing (l-1) + inp % dtime
!        end do
!     else
!        inp % timing (1) = inp % time_offset
!        do l = 2 , inp % nstat
!           inp % timing (l) = inp % timing (l-1) + inp % dtime
!        end do
!     end if

!     if ( inp % nvolume > nstatmax .or. &
!          inp % nxystat > nstatmax .or. &
!          inp % nxzstat > nstatmax .or. &
!          inp % nyzstat > nstatmax ) then
!        write (*,*) 'maximum number of stats reached'
!        call mpi_abort ( MPI_COMM_WORLD , 1 , mpicode )
!     end if


!   end subroutine readinp



!> \brief Read the file syntturb.dat.
!!
!! This datas are come from previous simulations or experiments. 
!! And use to feed the boundary conditions to generate synthetic turbulence at the inlet(s).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

! subroutine readinpext ( adi , xt , yt , zt , x , y , z , inp )


!    type (adi_type)                         , intent (in)      :: adi !< non-dimensional derived type
!    real (dp) , allocatable , dimension (:) , intent (in)      :: xt  !< absolute x-coordinate array
!    real (dp) , allocatable , dimension (:) , intent (in)      :: yt  !< absolute y-coordinate array
!    real (dp) , allocatable , dimension (:) , intent (in)      :: zt  !< absolute z-coordinate array
!    real (dp) , allocatable , dimension (:) , intent (in)      :: x   !< absolute x-coordinate array
!    real (dp) , allocatable , dimension (:) , intent (in)      :: y   !< absolute y-coordinate array
!    real (dp) , allocatable , dimension (:) , intent (in)      :: z   !< absolute z-coordinate array
!    type (inp_type)                         , intent (inout)   :: inp !< input derived type

!    integer (ip)                                   :: ok , i , j , k , l 
!    logical                                        :: loop

!    real (dp) , allocatable , dimension (:)        :: ytemp
!    real (dp) , allocatable , dimension (:,:,:)    :: uin , tke
!    real (dp)                                      :: dummy


!    integer (ip)                                   :: nb_stat, iter_stat , reclmax
!    real (dp) , allocatable , dimension (:,:,:,:)  :: var_init

!     nb_stat      = 128


!    allocate ( var_init(nb_stat,nty,ntz,5) , stat = ok)
!   if (rank == rank_default) then 
!    if (ok > 0) stop 'Memory allocation for var_init array has failed!'
!     reclmax = int ( nb_stat * nty * ntz * nv * nbit , kind = 8 )
!     open ( unit = unit_var_init , file      = file_var_ini    , &
!                                   access    = 'direct'        , &
!                                   recl      = reclmax         , &
!                                   form      = 'unformatted'   , &
!                                   status    = 'unknown'       , &
!                                   iostat    = ok              )

!     read ( unit = unit_var_init , rec = 1 ) ((((var_init (i,j,k,l),i=1,nb_stat),j=1,nty),k=1,ntz),l= 1,nv)
!      close(unit_var_init) 
!      end if
!      call MPI_BCAST ( var_init , nb_stat*nty*ntz*nv, MPI_DOUBLE_PRECISION , 0 , MPI_COMM_WORLD , mpicode )
!      allocate ( inp % ext % var_init(1:nb_stat,sy:ey,sz:ez,nv) , stat = ok)
!      if (ok > 0) stop 'Memory allocation for inp % ext % var_init array has failed!'
!      do i = 1 , nb_stat 
!        do j = sy , ey 
!          do k = sz , ez
!            do l = 1 , 5
!              inp % ext % var_init(i,j,k,l)=var_init(i,j,k,l) 
!            end do 
!          end do 
!        end do
!      end do

! !      deallocate(var_init)

! type (adi_type)                         , intent (in)      :: adi !< non-dimensional derived type
! real (dp) , allocatable , dimension (:) , intent (in)      :: xt  !< absolute x-coordinate array
! real (dp) , allocatable , dimension (:) , intent (in)      :: yt  !< absolute y-coordinate array
! real (dp) , allocatable , dimension (:) , intent (in)      :: zt  !< absolute z-coordinate array
! real (dp) , allocatable , dimension (:) , intent (in)      :: x   !< absolute x-coordinate array
! real (dp) , allocatable , dimension (:) , intent (in)      :: y   !< absolute y-coordinate array
! real (dp) , allocatable , dimension (:) , intent (in)      :: z   !< absolute z-coordinate array
! type (inp_type)                         , intent (inout)   :: inp !< input derived type

! integer (ip)                                   :: ok , i , j , k , l , n1 , n2 , n3
! real (dp) , allocatable , dimension (:,:,:,:)  :: var_init

! n1 = 82
! n2 = 82
! n3 = 82

! allocate ( inp % ext % var_init (n1,n2,n3,10) , stat = ok)
!     if ( ok /= 0 ) call abort_mpi ('error opening var_init.dat')
!     if (rank == rank_default) then
!        open ( unit = unit_var_init , file = file_var_ini , status = 'old' , iostat = ok )
!        do k = 1 , n1
!           do j = 1 , n2 
!              do i = 1 , n3 

!                 read( unit_var_init , * )  inp % ext % var_init(i,j,k,1)   , &
!                                                     inp % ext % var_init(i,j,k,2)   , &
!                                                        inp % ext % var_init(i,j,k,3)   , &
!                                                        inp % ext % var_init(i,j,k,4)   , &
!                                                        inp % ext % var_init(i,j,k,5)   , &
!                                            inp % ext % var_init(i,j,k,6)   , &
!                                                     inp % ext % var_init(i,j,k,7)   , &
!                                                        inp % ext % var_init(i,j,k,8)   , &
!                                                        inp % ext % var_init(i,j,k,9)   , &
!                                                        inp % ext % var_init(i,j,k,10) 
!              enddo
!           enddo 
!        enddo
!        close(unit_var_init)
!     endif

!     call MPI_BCAST ( inp % ext % var_init , n1*n2*n3*10 , MPI_DOUBLE_PRECISION , 0 , MPI_COMM_WORLD , mpicode )

! end subroutine readinpext



! subroutine interp_1d ( nd, xd, yd, ni, xi, yi )
!    implicit none

!    integer(ip) :: ni , i , k , nd
!    real(dp)    :: t
!    real(dp)    :: xd(nd) , yd(nd) , xi(ni) , yi(ni)


!    do i = 1, ni
!       yi(i) = 0.0D+00
!    end do

!    if ( nd .eq. 1 ) then
!       do i = 1, ni
!          yi(i) = yd(1)
!       end do
!    end if

!    do i = 1, ni
!       if ( xi(i) .le. xd(1) ) then
!          t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
!          yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)
!       else if ( xd(nd) .le. xi(i) ) then
!          t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
!          yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)
!       else
!          do k = 2, nd
!             if ( xd(k-1) .le. xi(i) .and. xi(i) .le. xd(k) ) then
!                t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
!                yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
!                exit
!             end if
!          end do
!       end if
!    end do

! end subroutine interp_1d






end module input
