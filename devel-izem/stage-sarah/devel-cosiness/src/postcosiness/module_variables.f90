!------------------------------------------------------------------------------
! MODULE: variables
!------------------------------------------------------------------------------
!> \brief Post-treatment variables.
!!
!! This module contains all the _additional_ variables that need to be
!! defined during the post-treatment.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module variables

  use parameters
  use input
  use parallel


  implicit none


  integer (ip) , parameter , private :: nstatmax = 10000 !< maximum number of statistical files


!> \brief Simulation derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type sim_type


     real (dp) , dimension (nstatmax)                     :: dtime_per_ite   !< non-dimensional time step of the simulation
     real (dp)                                            :: sumdtime        !< non-dimensional time of the simulation
     integer (ip)                                         :: nfiles          !< corrected total number of stat files to post-process
     integer (ip)                                         :: sfile , efile   !< start/end file for each MPI process (for subroutine transform)


  end type sim_type


!> \brief Instantaneous derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type inst_type


     !> special variable to reduce iterations
     real (dp) , dimension (:,:,:) , allocatable     :: T0

     !> vector of conserved variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: u

     !> timestep variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: dt

     !> reactive variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot_SGS
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot_hct
     real (dp) , dimension (:,:,:) , allocatable     :: omega
     real (dp) , dimension (:,:,:) , allocatable     :: omega_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: omega_hct
     real (dp) , dimension (:,:,:) , allocatable     :: reaction_det
     real (dp) , dimension (:,:,:) , allocatable     :: reaction_det_SGS
     real (dp) , dimension (:,:,:,:) , allocatable   :: Yass  ! steady state species
     real (dp) , dimension (:,:,:,:) , allocatable   :: SSa   ! chemical steady state parameter
     real (dp) , dimension (:,:,:) , allocatable     :: Da    ! Damkohler
     real (dp) , dimension (:,:,:) , allocatable     :: lambda ! reactivity rate

     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS
!     real (dp) , dimension (:,:,:) , allocatable     :: tau_iso_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: kpa
     real (dp) , dimension (:,:,:) , allocatable     :: ct
     real (dp) , dimension (:,:,:,:) , allocatable   :: dm
     real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
     real (dp) , dimension (:,:,:,:,:) , allocatable :: rd

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: H
     real (dp) , dimension (:,:,:) , allocatable     :: cp
     real (dp) , dimension (:,:,:) , allocatable     :: W_i

     !> Rate of progress of each elementary reaction
     real (dp) , dimension (:,:,:,:) , allocatable   :: W_scal_i

     !> Aditionnal variables (A.Techer)
     real (dp) , dimension (:,:,:) , allocatable     :: xiv2 ! Investigation of SGS variance modeling
     real (dp) , dimension (:,:,:) , allocatable     :: xiv3 !
     real (dp) , dimension (:,:,:) , allocatable     :: xixi !

     real (dp) , dimension (:,:,:) , allocatable     :: c         ! MIL
     real (dp) , dimension (:,:,:) , allocatable     :: zjm       !
     real (dp) , dimension (:,:,:) , allocatable     :: zjp       !
     real (dp) , dimension (:,:,:) , allocatable     :: yoequm    !
     real (dp) , dimension (:,:,:) , allocatable     :: yoburnt   !
     real (dp) , dimension (:,:,:) , allocatable     :: prob_mil  !
     real (dp) , dimension (:,:,:) , allocatable     :: wo2mil , wo2sauts , wo2bur !
     real (dp) , dimension (:,:,:) , allocatable     :: intzjmp   !

  end type inst_type


!> \brief Reynolds average derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type reyavg_type


     !> basic variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: Y
     real (dp) , dimension (:,:,:) , allocatable     :: rho
     real (dp) , dimension (:,:,:) , allocatable     :: ux
     real (dp) , dimension (:,:,:) , allocatable     :: vy
     real (dp) , dimension (:,:,:) , allocatable     :: wz
     real (dp) , dimension (:,:,:) , allocatable     :: et

     !> reactive variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot_hct
     real (dp) , dimension (:,:,:) , allocatable     :: omega
     real (dp) , dimension (:,:,:) , allocatable     :: omega_hct
     real (dp) , dimension (:,:,:,:) , allocatable   :: Yass
     real (dp) , dimension (:,:,:) , allocatable     :: Da
     real (dp) , dimension (:,:,:) , allocatable     :: Da_II
     real (dp) , dimension (:,:,:) , allocatable     :: sdr
     real (dp) , dimension (:,:,:) , allocatable     :: FI
     real (dp) , dimension (:,:,:) , allocatable     :: lambda ! reactivity rate
     real (dp) , dimension (:,:,:) , allocatable     :: proba_pz

     !> passive variables
     real (dp) , dimension (:,:,:) , allocatable     :: Zm

     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: tau_iso_SGS

     real (dp) , dimension (:,:,:,:,:) , allocatable :: tau_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: nu
     real (dp) , dimension (:,:,:) , allocatable     :: kpa
     real (dp) , dimension (:,:,:) , allocatable     :: ct
     ! real (dp) , dimension (:,:,:,:) , allocatable   :: dm
     ! real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
     ! real (dp) , dimension (:,:,:,:,:) , allocatable :: rd

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: P
     real (dp) , dimension (:,:,:) , allocatable     :: H
     real (dp) , dimension (:,:,:) , allocatable     :: cp
     real (dp) , dimension (:,:,:) , allocatable     :: W

     !> other basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: Tt
     real (dp) , dimension (:,:,:) , allocatable     :: Ht
     real (dp) , dimension (:,:,:) , allocatable     :: M
     real (dp) , dimension (:,:,:) , allocatable     :: cs
     real (dp) , dimension (:,:,:) , allocatable     :: gamma
     real (dp) , dimension (:,:,:) , allocatable     :: Pr

     !> custom variables
     real (dp) , dimension (:,:,:,:,:) , allocatable :: tau

     !> conditional averages
     real (dp) , dimension (:,:) , allocatable       :: ca_rho_Zm
     real (dp) , dimension (:,:) , allocatable       :: ca_rho_FI
     real (dp) , dimension (:,:) , allocatable       :: ca_rho_Da

     !> Aditionnal variables (A.Techer)
!     real (dp) , dimension (:,:,:) , allocatable     :: zjm
!     real (dp) , dimension (:,:,:) , allocatable     :: zjp
!     real (dp) , dimension (:,:,:) , allocatable     :: intzjmp
!     real (dp) , dimension (:,:,:) , allocatable     :: taumix

  end type reyavg_type


!> \brief Favre average derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type favavg_type


     !> basic variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: Y
     real (dp) , dimension (:,:,:) , allocatable     :: ux
     real (dp) , dimension (:,:,:) , allocatable     :: vy
     real (dp) , dimension (:,:,:) , allocatable     :: wz
     real (dp) , dimension (:,:,:) , allocatable     :: et

     !> reactive variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot_hct
     real (dp) , dimension (:,:,:) , allocatable     :: omega
     real (dp) , dimension (:,:,:) , allocatable     :: omega_hct
     real (dp) , dimension (:,:,:,:) , allocatable   :: Yass
     real (dp) , dimension (:,:,:) , allocatable     :: Da
     real (dp) , dimension (:,:,:) , allocatable     :: sdr
     real (dp) , dimension (:,:,:) , allocatable     :: FI
     real (dp) , dimension (:,:,:) , allocatable     :: lambda ! reactivity rate

     !> passive variables
     real (dp) , dimension (:,:,:) , allocatable     :: Zm

     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: tke_SGS
     real (dp) , dimension (:,:,:) , allocatable     :: nu
     real (dp) , dimension (:,:,:) , allocatable     :: kpa
     real (dp) , dimension (:,:,:) , allocatable     :: ct
     ! real (dp) , dimension (:,:,:,:) , allocatable   :: dm
     ! real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
     ! real (dp) , dimension (:,:,:,:,:) , allocatable :: rd

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: P
     real (dp) , dimension (:,:,:) , allocatable     :: H
     real (dp) , dimension (:,:,:) , allocatable     :: cp
     real (dp) , dimension (:,:,:) , allocatable     :: W

     !> other basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: Tt
     real (dp) , dimension (:,:,:) , allocatable     :: Ht
     real (dp) , dimension (:,:,:) , allocatable     :: M
     real (dp) , dimension (:,:,:) , allocatable     :: cs
     real (dp) , dimension (:,:,:) , allocatable     :: gamma
     real (dp) , dimension (:,:,:) , allocatable     :: Pr

     !> custom variables

     !> conditional averages
     real (dp) , dimension (:,:) , allocatable       :: ca_T_Zm
     real (dp) , dimension (:,:) , allocatable       :: ca_c_Zm
     real (dp) , dimension (:,:) , allocatable       :: ca_sdr_Zm
     real (dp) , dimension (:,:) , allocatable       :: ca_Da_Zm
     real (dp) , dimension (:,:) , allocatable       :: ca_omega_Zm
     real (dp) , dimension (:,:) , allocatable       :: ca_FI_Zm

     real (dp) , dimension (:,:) , allocatable       :: ca_omega_FI

     real (dp) , dimension (:,:) , allocatable       :: ca_omega_Da

     !> Aditionnal variables (A.Techer)
     real (dp) , dimension (:,:,:) , allocatable     :: xiv2 ! SGS variance
     real (dp) , dimension (:,:,:) , allocatable     :: xiv3
     real (dp) , dimension (:,:,:) , allocatable     :: xixi


  end type favavg_type


!> \brief Reynolds fluctuation derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type reyfluc_type


     !> basic variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: Y
     real (dp) , dimension (:,:,:) , allocatable     :: rho
     real (dp) , dimension (:,:,:) , allocatable     :: ux
     real (dp) , dimension (:,:,:) , allocatable     :: vy
     real (dp) , dimension (:,:,:) , allocatable     :: wz
     real (dp) , dimension (:,:,:) , allocatable     :: et

     !> reactive variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot_hct
     real (dp) , dimension (:,:,:) , allocatable     :: omega
     real (dp) , dimension (:,:,:) , allocatable     :: omega_hct
     real (dp) , dimension (:,:,:,:) , allocatable   :: Yass
     real (dp) , dimension (:,:,:) , allocatable     :: Da
     real (dp) , dimension (:,:,:) , allocatable     :: sdr
     real (dp) , dimension (:,:,:) , allocatable     :: FI

     !> passive variables
     real (dp) , dimension (:,:,:) , allocatable     :: Zm

     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS ! May be not need ?! (A.Techer)
     real (dp) , dimension (:,:,:) , allocatable     :: nu
     real (dp) , dimension (:,:,:) , allocatable     :: kpa
     real (dp) , dimension (:,:,:) , allocatable     :: ct
     ! real (dp) , dimension (:,:,:,:) , allocatable   :: dm
     ! real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
     ! real (dp) , dimension (:,:,:,:,:) , allocatable :: rd

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: P
     real (dp) , dimension (:,:,:) , allocatable     :: H
     real (dp) , dimension (:,:,:) , allocatable     :: cp
     real (dp) , dimension (:,:,:) , allocatable     :: W

     !> other basic variables

     real (dp) , dimension (:,:,:) , allocatable     :: Tt
     real (dp) , dimension (:,:,:) , allocatable     :: Ht
     real (dp) , dimension (:,:,:) , allocatable     :: M
     real (dp) , dimension (:,:,:) , allocatable     :: cs
     real (dp) , dimension (:,:,:) , allocatable     :: gamma
     real (dp) , dimension (:,:,:) , allocatable     :: Pr

     !> custom variables
     real (dp) , dimension (:,:,:,:,:) , allocatable :: tau



  end type reyfluc_type


!> \brief Favre fluctuation derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type favfluc_type


     !> basic variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: Y
     real (dp) , dimension (:,:,:) , allocatable     :: ux
     real (dp) , dimension (:,:,:) , allocatable     :: vy
     real (dp) , dimension (:,:,:) , allocatable     :: wz
     real (dp) , dimension (:,:,:) , allocatable     :: et

     !> reactive variables
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot
     real (dp) , dimension (:,:,:,:) , allocatable   :: Ydot_hct
     real (dp) , dimension (:,:,:) , allocatable     :: omega
     real (dp) , dimension (:,:,:) , allocatable     :: omega_hct
     real (dp) , dimension (:,:,:,:) , allocatable   :: Yass
     real (dp) , dimension (:,:,:) , allocatable     :: Da
     real (dp) , dimension (:,:,:) , allocatable     :: sdr
     real (dp) , dimension (:,:,:) , allocatable     :: FI

     !> passive variables
     real (dp) , dimension (:,:,:) , allocatable     :: Zm

     !> viscous variables
     real (dp) , dimension (:,:,:) , allocatable     :: mu
     real (dp) , dimension (:,:,:) , allocatable     :: mu_SGS ! May be not need ?! (A.Techer)
     real (dp) , dimension (:,:,:) , allocatable     :: nu
     real (dp) , dimension (:,:,:) , allocatable     :: kpa
     real (dp) , dimension (:,:,:) , allocatable     :: ct
     ! real (dp) , dimension (:,:,:,:) , allocatable   :: dm
     ! real (dp) , dimension (:,:,:,:) , allocatable   :: tdr
     ! real (dp) , dimension (:,:,:,:,:) , allocatable :: rd

     !> thermodynamic variables
     real (dp) , dimension (:,:,:) , allocatable     :: T
     real (dp) , dimension (:,:,:) , allocatable     :: P
     real (dp) , dimension (:,:,:) , allocatable     :: H
     real (dp) , dimension (:,:,:) , allocatable     :: cp
     real (dp) , dimension (:,:,:) , allocatable     :: W

     !> other basic variables
     real (dp) , dimension (:,:,:) , allocatable     :: Tt
     real (dp) , dimension (:,:,:) , allocatable     :: Ht
     real (dp) , dimension (:,:,:) , allocatable     :: M
     real (dp) , dimension (:,:,:) , allocatable     :: cs
     real (dp) , dimension (:,:,:) , allocatable     :: gamma
     real (dp) , dimension (:,:,:) , allocatable     :: Pr

     !> custom variables

     !> conditional fluctuations
     real (dp) , dimension (:,:,:) , allocatable     :: ca_T_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: ca_c_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: ca_sdr_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: ca_Da_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: ca_omega_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: ca_FI_Zm

     real (dp) , dimension (:,:,:) , allocatable     :: ca_omega_FI

     real (dp) , dimension (:,:,:) , allocatable     :: ca_omega_Da


  end type favfluc_type


!> \brief Statistical derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type stat_type


     !> _repeated_ reynolds averages for convergence plots
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_P
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_Y
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_ux
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_vy
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_wz

     !> other averages for convergence plots
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_uxux
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_vyvy
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_wzwz
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_uxvy
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_uxwz
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_vywz
     real (dp) , dimension (:,:,:) , allocatable     :: reyavg_enstrophy

     !> other averages
     real (dp) , dimension (:,:,:) , allocatable     :: favavg_ux_Y
     real (dp) , dimension (:,:,:) , allocatable     :: favavg_vy_Y
     real (dp) , dimension (:,:,:) , allocatable     :: favavg_wz_Y
     real (dp) , dimension (:,:,:) , allocatable     :: favavg_ux_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: favavg_vy_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: favavg_wz_Zm

     !> variances (typical)
     real (dp) , dimension (:,:,:,:) , allocatable   :: reyvar_Y
     real (dp) , dimension (:,:,:,:) , allocatable   :: favvar_Y
     real (dp) , dimension (:,:,:,:) , allocatable   :: reyvar_Yass
     real (dp) , dimension (:,:,:,:) , allocatable   :: favvar_Yass
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_Zm

     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_ux
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_vy
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_wz
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_uxvy
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_uxwz
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_vywz
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_uxvywz

     real (dp) , dimension (:,:,:) , allocatable     :: favvar_ux
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_vy
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_wz
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_uxvy
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_uxwz
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_vywz
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_uxvywz

     !> variances (other)
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_p
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_rho
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_rho_acu
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_rho_ent
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_T
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_T_acu
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_T_ent
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_W
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_du1_dx1
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_p_du1_dx1
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_rho_p
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_rho_W

     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_omega
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_Da
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_sdr
     real (dp) , dimension (:,:,:) , allocatable     :: reyvar_FI
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_omega
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_Da
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_sdr
     real (dp) , dimension (:,:,:) , allocatable     :: favvar_FI

     !> conditional variances
     real (dp) , dimension (:,:) , allocatable       :: favvar_ca_T_Zm
     real (dp) , dimension (:,:) , allocatable       :: favvar_ca_c_Zm
     real (dp) , dimension (:,:) , allocatable       :: favvar_ca_sdr_Zm
     real (dp) , dimension (:,:) , allocatable       :: favvar_ca_Da_Zm
     real (dp) , dimension (:,:) , allocatable       :: favvar_ca_omega_Zm
     real (dp) , dimension (:,:) , allocatable       :: favvar_ca_FI_Zm

     real (dp) , dimension (:,:) , allocatable       :: favvar_ca_omega_FI

     real (dp) , dimension (:,:) , allocatable       :: favvar_ca_omega_Da

     !> energy transport equation
     real (dp) , dimension (:,:,:) , allocatable     :: tke_dissip
     real (dp) , dimension (:,:,:) , allocatable     :: tke_press_strain
     real (dp) , dimension (:,:,:,:) , allocatable   :: tke_transp
     real (dp) , dimension (:,:,:) , allocatable     :: tke_ra_ff_ux
     real (dp) , dimension (:,:,:) , allocatable     :: tke_ra_ff_vy
     real (dp) , dimension (:,:,:) , allocatable     :: tke_ra_ff_wz

     !> vorticity transport equation
     real (dp) , dimension (:,:,:) , allocatable     :: vor_conv
     real (dp) , dimension (:,:,:) , allocatable     :: ens_conv
     real (dp) , dimension (:,:,:) , allocatable     :: vor_stret
     real (dp) , dimension (:,:,:) , allocatable     :: ens_stret
     real (dp) , dimension (:,:,:) , allocatable     :: vor_dila
     real (dp) , dimension (:,:,:) , allocatable     :: ens_dila
     real (dp) , dimension (:,:,:) , allocatable     :: vor_baro
     real (dp) , dimension (:,:,:) , allocatable     :: ens_baro
     real (dp) , dimension (:,:,:) , allocatable     :: vor_visc
     real (dp) , dimension (:,:,:) , allocatable     :: ens_visc

     !> scalar variance transport equation
     real (dp) , dimension (:,:,:,:) , allocatable   :: scal_var_turb_transp
     real (dp) , dimension (:,:,:,:) , allocatable   :: scal_var_turb_prod
     real (dp) , dimension (:,:,:,:) , allocatable   :: scal_var_mol_diff
     real (dp) , dimension (:,:,:) , allocatable     :: scal_var_dissip
     real (dp) , dimension (:,:,:) , allocatable     :: scal_var_dissip2
     real (dp) , dimension (:,:,:) , allocatable     :: scal_var_dissip3
     real (dp) , dimension (:,:,:) , allocatable     :: scal_var_dissip4
     real (dp) , dimension (:,:,:) , allocatable     :: scal_var_dissip5
     real (dp) , dimension (:,:,:) , allocatable     :: scal_var_dissip6

     !> spatial correlations
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_ux
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_vy
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_wz
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_p
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_Y
     real (dp) , dimension (:,:,:) , allocatable     :: tmp_corr_Zm

     !> Reynolds transport equation
     real (dp) , dimension (:,:,:) , allocatable     :: rey_dissip_11 , rey_dissip_22 , rey_dissip_33 , &
                                                        rey_dissip_12 , rey_dissip_13 , rey_dissip_23
     real (dp) , dimension (:,:,:) , allocatable     :: rey_press_strain_11 , rey_press_strain_22 , rey_press_strain_33 , &
                                                        rey_press_strain_12 , rey_press_strain_13 , rey_press_strain_23
     real (dp) , dimension (:,:,:,:) , allocatable   :: rey_transp_11 , rey_transp_22 , rey_transp_33 , &
                                                        rey_transp_12 , rey_transp_13 , rey_transp_23

     !> pdfs
     real (dp) , dimension (:,:,:) , allocatable     :: pdf_Y
     real (dp) , dimension (:,:,:) , allocatable     :: pdf_Zm
     real (dp) , dimension (:,:,:) , allocatable     :: pdf_c
     real (dp) , dimension (:,:,:) , allocatable     :: pdf_sdr
     real (dp) , dimension (:,:,:) , allocatable     :: pdf_Da
     real (dp) , dimension (:,:,:) , allocatable     :: pdf_FI

     !> conditional pdfs
     real (dp) , dimension (:,:) , allocatable       :: cpdf_omega_FI
     real (dp) , dimension (:,:) , allocatable       :: cpdf_omega_Da

     !> Fourier spectra (only the amplitude)
     real (dp) , dimension (:,:) , allocatable       :: spec_ux
     real (dp) , dimension (:,:) , allocatable       :: spec_vy
     real (dp) , dimension (:,:) , allocatable       :: spec_wz
     real (dp) , dimension (:,:) , allocatable       :: spec_Y
     real (dp) , dimension (:,:) , allocatable       :: spec_Zm
     real (dp) , dimension (:,:) , allocatable       :: spec_sdr


  end type stat_type


contains


!> \brief Allocate instantaneous derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine alloc_inst_type (inst)


    type (inst_type) , intent (inout)                            :: inst !< instantaneous derived type


    integer (ip) :: ok


    allocate ( inst % T0    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               inst % u     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv )            , &

               inst % dt    ( sx:ex , sy:ey , sz:ez , 1:5 )                               , & ! A.Techer

               inst % Ydot  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               inst % omega ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               inst % Da    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               inst % lambda ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               inst % reaction_det ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )            , &
               inst % SSa   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &

               inst % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               inst % kpa   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               inst % ct    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               inst % dm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               inst % tdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               inst % rd    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv , 1:nrv )   , &

               inst % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               inst % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               inst % cp    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               inst % W_i   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               inst % W_scal_i ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nreac )      , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate inst type')


    if ( nss > 0 ) then

       allocate ( inst % Yass  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nss )        , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate inst type 1')

    end if


    if ( LES ) then

       allocate ( inst % mu_SGS( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , & ! A.Techer
!                  inst % tau_iso_SGS( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) ,&
!                  inst % xiv2  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , & ! SGS variance
!                  inst % xiv3  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
!                  inst % xixi  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

                  inst % Ydot_SGS  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )       , &
                  inst % Ydot_hct  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )       , &
                  inst % omega_SGS ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
                  inst % omega_hct ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
                  inst % reaction_det_SGS ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )        , &

                  inst % c       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , & ! MIL
                  inst % zjm     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , &
                  inst % zjp     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , &
                  inst % yoequm  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , &
                  inst % yoburnt ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , &
                  inst % prob_mil( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , &
                  inst % wo2mil  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , &
                  inst % wo2sauts( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , &
                  inst % wo2bur  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , &
                  inst % intzjmp ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                 , &

                  stat = ok )

       if ( ok > 0 ) call abort_mpi ('error allocate inst type 2')

    end if ! LES


    inst % T0    = 0.0_dp

    ! inst % u     = 0.0_dp

    ! inst % Ydot  = 0.0_dp
    ! inst % omega = 0.0_dp

    ! inst % mu    = 0.0_dp
    ! inst % kpa   = 0.0_dp
    ! inst % ct    = 0.0_dp
    ! inst % dm    = 0.0_dp
    ! inst % tdr   = 0.0_dp
    ! inst % rd    = 0.0_dp

    ! inst % T     = 0.0_dp
    ! inst % H     = 0.0_dp
    ! inst % cp    = 0.0_dp
    ! inst % W_i   = 0.0_dp

    ! inst % W_scal_i = 0.0_dp


  end subroutine alloc_inst_type


!> \brief Write instantaneous derived type declaration into a file. (Subroutine not used... ?!)
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
!  subroutine write_inst_type (inst)


!    type (inst_type) , intent (in)                               :: inst !< instantaneous derived type


!    integer (ip)                 :: ok , i , j , k , l , m
!    character (len_default)      :: number_ascii
!
!
!    write ( number_ascii , format_restart ) rank
!
!
!    open ( unit_stat , file   = trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) , &
!                       form   = 'unformatted'   ,                                                      &
!                       status = 'unknown'       ,                                                      &
!                       iostat = ok              )
!
!    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) )
!
!
!    write (unit_stat) (((( inst % u (i,j,k,l) , i = sx , ex ) , &
!                                                j = sy , ey ) , &
!                                                k = sz , ez ) , &
!                                                l = 1  , nv )
!
!    write (unit_stat) (((( inst % Ydot (i,j,k,l) , i = sx , ex  ) , &
!                                                   j = sy , ey  ) , &
!                                                   k = sz , ez  ) , &
!                                                   l = 1  , nrv )
!
!    write (unit_stat) ((( inst % omega (i,j,k) , i = sx , ex ) , &
!                                                 j = sy , ey ) , &
!                                                 k = sz , ez )
!
!    write (unit_stat) ((( inst % mu (i,j,k) , i = sx , ex ) , &
!                                              j = sy , ey ) , &
!                                              k = sz , ez )
!
!
!    write (unit_stat) ((( inst % kpa (i,j,k) , i = sx , ex ) , &
!                                               j = sy , ey ) , &
!                                               k = sz , ez )
!
!    write (unit_stat) ((( inst % ct (i,j,k) , i = sx , ex ) , &
!                                              j = sy , ey ) , &
!                                              k = sz , ez )
!
!    write (unit_stat) (((( inst % dm (i,j,k,l) , i = sx , ex  ) , &
!                                                 j = sy , ey  ) , &
!                                                 k = sz , ez  ) , &
!                                                 l = 1  , nrv )
!
!    write (unit_stat) (((( inst % tdr (i,j,k,l) , i = sx , ex  ) , &
!                                                  j = sy , ey  ) , &
!                                                  k = sz , ez  ) , &
!                                                  l = 1  , nrv )
!
!    write (unit_stat) ((((( inst % rd (i,j,k,l,m) , i = sx , ex  ) , &
!                                                    j = sy , ey  ) , &
!                                                    k = sz , ez  ) , &
!                                                    l = 1  , nrv ) , &
!                                                    m = 1  , nrv )
!
!    write (unit_stat) ((( inst % T (i,j,k) , i = sx , ex ) , &
!                                             j = sy , ey ) , &
!                                             k = sz , ez )
!
!    write (unit_stat) ((( inst % H (i,j,k) , i = sx , ex ) , &
!                                             j = sy , ey ) , &
!                                             k = sz , ez )
!
!    write (unit_stat) ((( inst % cp (i,j,k) , i = sx , ex ) , &
!                                              j = sy , ey ) , &
!                                              k = sz , ez )
!
!    write (unit_stat) ((( inst % W_i (i,j,k) , i = sx , ex ) , &
!                                               j = sy , ey ) , &
!                                               k = sz , ez )
!
!    write (unit_stat) (((( inst % W_scal_i (i,j,k,l) , i = sx , ex    ) , &
!                                                       j = sy , ey    ) , &
!                                                       k = sz , ez    ) , &
!                                                       l = 1  , nreac )
!
!    close (unit_stat)
!
!
!  end subroutine write_inst_type


!> \brief Read instantaneous derived type declaration from a file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine read_inst_type ( inp , inst )

    type (inp_type)  , intent (in)                               :: inp  !< input derived type
    type (inst_type) , intent (inout)                            :: inst !< instantaneous derived type


    integer (ip)                    :: ok , irank , i , j , k , l , m

    integer (ip)                    :: px , py , pz , ix , fx , iy , fy , iz , fz

    integer (ip) , dimension (3)    :: mycoords , mydims

    character (len_default)         :: number_ascii


    mydims (1) = inp % nprocz ; mydims (2) = inp % nprocy ; mydims (3) = inp % nprocx
    mycoords (:) = 0
    irank        = 0


    if ( nproc == 1 ) then ! sequential problem


       do pz = 1 , inp % nprocz

          iz = ( mycoords(1)*ntz ) / mydims(1) + 1
          fz = ( ( mycoords(1)+1 )*ntz ) / mydims(1)

          mycoords (1) = mycoords (1) + 1
          mycoords (2) = 0
          mycoords (3) = 0

          do py = 1 , inp % nprocy

             iy = ( mycoords(2)*nty ) / mydims(2) + 1
             fy = ( ( mycoords(2)+1 )*nty ) / mydims(2)

             mycoords (2) = mycoords (2) + 1
             mycoords (3) = 0

             do px = 1 , inp % nprocx

                ix = ( mycoords(3)*ntx ) / mydims(3) + 1
                fx = ( ( mycoords(3)+1 )*ntx ) / mydims(3)
                mycoords (3) = mycoords (3) + 1


                write ( number_ascii , format_restart ) irank
                irank = irank + 1

                open ( unit_stat , file   = trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) , &
                                   form   = 'unformatted'   ,                                                      &
                                   status = 'unknown'       ,                                                      &
                                   iostat = ok              )

                if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // &
                    &                      trim (file_inst) // '_' // trim (number_ascii) )


                read (unit_stat) (((( inst % u (i,j,k,l) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz ) , &
                                                           l = 1  , nv )

                read (unit_stat) (((( inst % Ydot (i,j,k,l) , i = ix , fx  ) , &
                                                              j = iy , fy  ) , &
                                                              k = iz , fz  ) , &
                                                              l = 1  , nrv )

                read (unit_stat) (((( inst % Ydot_hct (i,j,k,l) , i = ix , fx  ) , &
                                                              j = iy , fy  ) , &
                                                              k = iz , fz  ) , &
                                                              l = 1  , nrv )

                read (unit_stat) ((( inst % omega (i,j,k) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            k = iz , fz )

                read (unit_stat) ((( inst % omega_hct (i,j,k) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            k = iz , fz )


                read (unit_stat) ((( inst % mu (i,j,k) , i = ix , fx ) , &
                                                         j = iy , fy ) , &
                                                         k = iz , fz )

                read (unit_stat) ((( inst % kpa (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( inst % ct (i,j,k) , i = ix , fx ) , &
                                                         j = iy , fy ) , &
                                                         k = iz , fz )

                read (unit_stat) (((( inst % dm (i,j,k,l) , i = ix , fx  ) , &
                                                            j = iy , fy  ) , &
                                                            k = iz , fz  ) , &
                                                            l = 1  , nrv )

                read (unit_stat) (((( inst % tdr (i,j,k,l) , i = ix , fx  ) , &
                                                             j = iy , fy  ) , &
                                                             k = iz , fz  ) , &
                                                             l = 1  , nrv )

                read (unit_stat) ((((( inst % rd (i,j,k,l,m) , i = ix , fx  ) , &
                                                               j = iy , fy  ) , &
                                                               k = iz , fz  ) , &
                                                               l = 1  , nrv ) , &
                                                               m = 1  , nrv )

                read (unit_stat) ((( inst % T (i,j,k) , i = ix , fx ) , &
                                                        j = iy , fy ) , &
                                                        k = iz , fz )

                read (unit_stat) ((( inst % H (i,j,k) , i = ix , fx ) , &
                                                        j = iy , fy ) , &
                                                        k = iz , fz )

                read (unit_stat) ((( inst % cp (i,j,k) , i = ix , fx ) , &
                                                         j = iy , fy ) , &
                                                         k = iz , fz )

                read (unit_stat) ((( inst % W_i (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) (((( inst % W_scal_i (i,j,k,l) , i = ix , fx    ) , &
                                                                  j = iy , fy    ) , &
                                                                  k = iz , fz    ) , &
                                                                  l = 1  , nreac )

             end do
          end do
       end do


    else ! parallel problem


       write ( number_ascii , format_restart ) rank

       open ( unit_stat , file   = trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) , &
                          form   = 'unformatted'   ,                                                      &
                          status = 'unknown'       ,                                                      &
                          iostat = ok              )

       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_inst) // '_' // trim (number_ascii) )


       read (unit_stat) (((( inst % u (i,j,k,l) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez ) , &
                                                  l = 1  , nv )

       read (unit_stat) (((( inst % Ydot (i,j,k,l) , i = sx , ex  ) , &
                                                     j = sy , ey  ) , &
                                                     k = sz , ez  ) , &
                                                     l = 1  , nrv )

       read (unit_stat) (((( inst % Ydot_hct (i,j,k,l) , i = sx , ex  ) , &
                                                     j = sy , ey  ) , &
                                                     k = sz , ez  ) , &
                                                     l = 1  , nrv )


       read (unit_stat) ((( inst % omega (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )

       read (unit_stat) ((( inst % omega_hct (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )


       read (unit_stat) ((( inst % mu (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

       read (unit_stat) ((( inst % kpa (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( inst % ct (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

       read (unit_stat) (((( inst % dm (i,j,k,l) , i = sx , ex  ) , &
                                                   j = sy , ey  ) , &
                                                   k = sz , ez  ) , &
                                                   l = 1  , nrv )

       read (unit_stat) (((( inst % tdr (i,j,k,l) , i = sx , ex  ) , &
                                                    j = sy , ey  ) , &
                                                    k = sz , ez  ) , &
                                                    l = 1  , nrv )

       read (unit_stat) ((((( inst % rd (i,j,k,l,m) , i = sx , ex  ) , &
                                                      j = sy , ey  ) , &
                                                      k = sz , ez  ) , &
                                                      l = 1  , nrv ) , &
                                                      m = 1  , nrv )

       read (unit_stat) ((( inst % T (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

       read (unit_stat) ((( inst % H (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

       read (unit_stat) ((( inst % cp (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

       read (unit_stat) ((( inst % W_i (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) (((( inst % W_scal_i (i,j,k,l) , i = sx , ex    ) , &
                                                         j = sy , ey    ) , &
                                                         k = sz , ez    ) , &
                                                         l = 1  , nreac )

    end if


    close (unit_stat)


  end subroutine read_inst_type


!> \brief Allocate Reynolds average derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine alloc_reyavg_type (reyavg)


    type (reyavg_type) , intent (out)                               :: reyavg !< Reynolds average derived type


    integer (ip) :: ok


    allocate ( reyavg % Y     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv+npv+nvv )   , &
               reyavg % rho   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % ux    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % vy    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % wz    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % et    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyavg % Ydot  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               reyavg % Ydot_hct  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )       , &
               reyavg % omega ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % omega_hct ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % Da    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % Da_II ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % sdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % FI    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % lambda( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % proba_pz ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyavg % Zm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyavg % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % nu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % kpa   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % ct    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               ! reyavg % dm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               ! reyavg % tdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               ! reyavg % rd    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv , 1:nrv )   , &

               reyavg % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % cp    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % W     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyavg % P     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % Tt    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % Ht    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % M     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % cs    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % gamma ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyavg % Pr    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyavg % tau ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax ) , &

               reyavg % ca_rho_Zm ( sx-ng:ex+ng , 1:nbins )                                 , &
               reyavg % ca_rho_FI ( sx-ng:ex+ng , 1:nbins )                                 , &
               reyavg % ca_rho_Da ( sx-ng:ex+ng , 1:nbins )                                 , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate reyavg type')

    if ( nss > 0 ) then

       allocate ( reyavg % Yass  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nss )         , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate reyavg type 2')

    end if

    if ( LES ) then

       allocate ( reyavg % mu_SGS ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , & ! A.Techer
                  reyavg % tau_iso_SGS( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                  reyavg % tau_SGS ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax )  , &

!                  reyavg % zjm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
!                  reyavg % zjp    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
!                  reyavg % intzjmp( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
!                  reyavg % taumix ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

                  stat = ok )

       if ( ok > 0 ) call abort_mpi ('error allocate reyavg type 3')

    end if


    ! reyavg % Y     = 0.0_dp
    ! reyavg % rho   = 0.0_dp
    ! reyavg % ux    = 0.0_dp
    ! reyavg % vy    = 0.0_dp
    ! reyavg % wz    = 0.0_dp
    ! reyavg % et    = 0.0_dp

    ! reyavg % Ydot  = 0.0_dp
    ! reyavg % omega = 0.0_dp
    ! reyavg % Da    = 0.0_dp
    ! reyavg % sdr   = 0.0_dp
    ! reyavg % FI    = 0.0_dp

    ! reyavg % Zm    = 0.0_dp

    ! reyavg % mu    = 0.0_dp
    ! reyavg % nu    = 0.0_dp
    ! reyavg % kpa   = 0.0_dp
    ! reyavg % ct    = 0.0_dp
    ! ! reyavg % dm    = 0.0_dp
    ! ! reyavg % tdr   = 0.0_dp
    ! ! reyavg % dr    = 0.0_dp

    ! reyavg % T     = 0.0_dp
    ! reyavg % H     = 0.0_dp
    ! reyavg % cp    = 0.0_dp
    ! reyavg % W     = 0.0_dp

    ! reyavg % P     = 0.0_dp
    ! reyavg % Tt    = 0.0_dp
    ! reyavg % Ht    = 0.0_dp
    ! reyavg % cs    = 0.0_dp
    ! reyavg % gamma = 0.0_dp
    ! reyavg % Pr    = 0.0_dp

    ! reyavg % tau   = 0.0_dp

    ! reyavg % ca_rho_Zm = 0.0_dp
    ! reyavg % ca_rho_FI = 0.0_dp
    ! reyavg % ca_rho_Da = 0.0_dp


  end subroutine alloc_reyavg_type


!> \brief Write Reynolds average derived type into a file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine write_reyavg_type (reyavg)


    type (reyavg_type) , intent (in)                                :: reyavg !< Reynolds average derived type


    integer (ip)                 :: ok , i , j , k , l , m
    character (len_default)      :: number_ascii


    write ( number_ascii , format_restart ) rank


    open ( unit_stat , file   = trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) , &
                       form   = 'unformatted'   ,                                                        &
                       status = 'unknown'       ,                                                        &
                       iostat = ok              )

    if ( ok /= 0 ) then
       call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) )
    end if


    write (unit_stat) (((( reyavg % Y (i,j,k,l) , i = sx , ex      ) , &
                                                  j = sy , ey      ) , &
                                                  k = sz , ez      ) , &
                                                  l = 1  , nrv+npv+nvv )

    write (unit_stat) ((( reyavg % rho (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

    write (unit_stat) ((( reyavg % ux (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % vy (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % wz (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % et (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) (((( reyavg % Ydot (i,j,k,l) , i = sx , ex  ) , &
                                                     j = sy , ey  ) , &
                                                     k = sz , ez  ) , &
                                                     l = 1  , nrv )

    write (unit_stat) (((( reyavg % Ydot_hct (i,j,k,l) , i = sx , ex  ) , &
                                                     j = sy , ey  ) , &
                                                     k = sz , ez  ) , &
                                                     l = 1  , nrv )


    write (unit_stat) ((( reyavg % omega (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )

    write (unit_stat) ((( reyavg % omega_hct (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )


    write (unit_stat) ((( reyavg % Da (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % sdr (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

    write (unit_stat) ((( reyavg % FI (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % Zm (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % mu (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % nu (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % kpa (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

    write (unit_stat) ((( reyavg % ct (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % T (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( reyavg % H (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( reyavg % cp (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( reyavg % W (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( reyavg % P (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( reyavg % Tt (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( reyavg % Ht (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( reyavg % M (i,j,k) , i = sx , ex ) , &
                                              j = sy , ey ) , &
                                              k = sz , ez )

   write (unit_stat) ((( reyavg % cs (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( reyavg % gamma (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

   write (unit_stat) ((( reyavg % Pr (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((((( reyavg % tau (i,j,k,l,m) , i = sx , ex )     , &
                                                      j = sy , ey )     , &
                                                      k = sz , ez )     , &
                                                      l = 1 , ndimmax ) , &
                                                      m = 1 , ndimmax )

   write (unit_stat) (( reyavg % ca_rho_Zm (i,l) , i = sx , ex   ) , &
                                                   l = 1 , nbins )

   write (unit_stat) (( reyavg % ca_rho_FI (i,l) , i = sx , ex   ) , &
                                                   l = 1 , nbins )

   write (unit_stat) (( reyavg % ca_rho_Da (i,l) , i = sx , ex   ) , &
                                                   l = 1 , nbins )


   close (unit_stat)


  end subroutine write_reyavg_type


!> \brief Read Reynolds average derived type from a file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine read_reyavg_type ( inp , reyavg )


    type (inp_type) , intent (in)                                 :: inp    !< input derived type
    type (reyavg_type) , intent (inout)                           :: reyavg !< Reynolds average derived type


    integer (ip)                    :: ok , irank , i , j , k , l , m

    integer (ip)                    :: px , py , pz , ix , fx , iy , fy , iz , fz

    integer (ip) , dimension (3)    :: mycoords , mydims

    character (len_default)         :: number_ascii


    mydims (1) = inp % nprocz ; mydims (2) = inp % nprocy ; mydims (3) = inp % nprocx
    mycoords (:) = 0
    irank      = 0


    if ( nproc == 1 ) then ! sequential problem


       do pz = 1 , inp % nprocz

          iz = ( mycoords(1)*ntz ) / mydims(1) + 1
          fz = ( ( mycoords(1)+1 )*ntz ) / mydims(1)

          mycoords (1) = mycoords (1) + 1
          mycoords (2) = 0
          mycoords (3) = 0

          do py = 1 , inp % nprocy

             iy = ( mycoords(2)*nty ) / mydims(2) + 1
             fy = ( ( mycoords(2)+1 )*nty ) / mydims(2)

             mycoords (2) = mycoords (2) + 1
             mycoords (3) = 0

             do px = 1 , inp % nprocx

                ix = ( mycoords(3)*ntx ) / mydims(3) + 1
                fx = ( ( mycoords(3)+1 )*ntx ) / mydims(3)
                mycoords (3) = mycoords (3) + 1


                write ( number_ascii , format_restart ) irank
                irank = irank + 1

                open ( unit_stat , file   = trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) , &
                                   form   = 'unformatted'   ,                                                        &
                                   status = 'unknown'       ,                                                        &
                                   iostat = ok              )

                if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // &
                    &                   trim (file_reyavg) // '_' // trim (number_ascii) )


                read (unit_stat) (((( reyavg % Y (i,j,k,l) , i = ix , fx      ) , &
                                                             j = iy , fy      ) , &
                                                             k = iz , fz      ) , &
                                                             l = 1  , nrv+npv+nvv )

                read (unit_stat) ((( reyavg % rho (i,j,k) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            k = iz , fz )

                read (unit_stat) ((( reyavg % ux (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % vy (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % wz (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % et (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) (((( reyavg % Ydot (i,j,k,l) , i = ix , fx  ) , &
                                                                j = iy , fy  ) , &
                                                                k = iz , fz  ) , &
                                                                l = 1  , nrv )

                read (unit_stat) (((( reyavg % Ydot_hct (i,j,k,l) , i = ix , fx  ) , &
                                                                j = iy , fy  ) , &
                                                                k = iz , fz  ) , &
                                                                l = 1  , nrv )


                read (unit_stat) ((( reyavg % omega (i,j,k) , i = ix , fx ) , &
                                                              j = iy , fy ) , &
                                                              k = iz , fz )

                read (unit_stat) ((( reyavg % omega_hct (i,j,k) , i = ix , fx ) , &
                                                              j = iy , fy ) , &
                                                              k = iz , fz )


                read (unit_stat) ((( reyavg % Da (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % sdr (i,j,k) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            k = iz , fz )

                read (unit_stat) ((( reyavg % FI (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % Zm (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % mu (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % nu (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % kpa (i,j,k) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            k = iz , fz )

                read (unit_stat) ((( reyavg % ct (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % T (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( reyavg % H (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( reyavg % cp (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % W (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( reyavg % P (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( reyavg % Tt (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % Ht (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % M (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( reyavg % cs (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( reyavg % gamma (i,j,k) , i = ix , fx ) , &
                                                              j = iy , fy ) , &
                                                              k = iz , fz )

                read (unit_stat) ((( reyavg % Pr (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((((( reyavg % tau (i,j,k,l,m) , i = ix , fx )     , &
                                                                  j = iy , fy )     , &
                                                                  k = iz , fz )     , &
                                                                  l = 1 , ndimmax ) , &
                                                                  m = 1 , ndimmax )

                read (unit_stat) (( reyavg % ca_rho_Zm (i,l) , i = ix , fx   ) , &
                                                               l = 1 , nbins )

                read (unit_stat) (( reyavg % ca_rho_FI (i,l) , i = ix , fx   ) , &
                                                               l = 1 , nbins )

                read (unit_stat) (( reyavg % ca_rho_Da (i,l) , i = ix , fx   ) , &
                                                               l = 1 , nbins )


             end do
          end do
       end do


    else ! parallel problem


       write ( number_ascii , format_restart ) rank

       open ( unit_stat , file   = trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) , &
                          form   = 'unformatted'   ,                                                        &
                          status = 'unknown'       ,                                                        &
                          iostat = ok              )

       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_reyavg) // '_' // trim (number_ascii) )


       read (unit_stat) (((( reyavg % Y (i,j,k,l) , i = sx , ex      ) , &
                                                    j = sy , ey      ) , &
                                                    k = sz , ez      ) , &
                                                    l = 1  , nrv+npv+nvv )

       read (unit_stat) ((( reyavg % rho (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )

       read (unit_stat) ((( reyavg % ux (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % vy (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % wz (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % et (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) (((( reyavg % Ydot (i,j,k,l) , i = sx , ex  ) , &
                                                       j = sy , ey  ) , &
                                                       k = sz , ez  ) , &
                                                       l = 1  , nrv )

       read (unit_stat) (((( reyavg % Ydot_hct (i,j,k,l) , i = sx , ex  ) , &
                                                       j = sy , ey  ) , &
                                                       k = sz , ez  ) , &
                                                       l = 1  , nrv )


       read (unit_stat) ((( reyavg % omega (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

       read (unit_stat) ((( reyavg % omega_hct (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )


       read (unit_stat) ((( reyavg % Da (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % sdr (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )

       read (unit_stat) ((( reyavg % FI (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % Zm (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % mu (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % nu (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % kpa (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )

       read (unit_stat) ((( reyavg % ct (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % T (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( reyavg % H (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( reyavg % cp (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % W (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( reyavg % P (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( reyavg % Tt (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % Ht (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % M (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( reyavg % cs (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( reyavg % gamma (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

       read (unit_stat) ((( reyavg % Pr (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((((( reyavg % tau (i,j,k,l,m) , i = sx , ex )     , &
                                                         j = sy , ey )     , &
                                                         k = sz , ez )     , &
                                                         l = 1 , ndimmax ) , &
                                                         m = 1 , ndimmax )

       read (unit_stat) (( reyavg % ca_rho_Zm (i,l) , i = sx , ex   ) , &
                                                      l = 1 , nbins )

       read (unit_stat) (( reyavg % ca_rho_FI (i,l) , i = sx , ex   ) , &
                                                      l = 1 , nbins )

       read (unit_stat) (( reyavg % ca_rho_Da (i,l) , i = sx , ex   ) , &
                                                      l = 1 , nbins )


    end if


    close (unit_stat)


  end subroutine read_reyavg_type


!> \brief Allocate Favre average derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine alloc_favavg_type (favavg)


    type (favavg_type) , intent (out)                               :: favavg !< Favre average derived type


    integer (ip) :: ok


    allocate ( favavg % Y     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv+npv+nvv )   , &
               favavg % ux    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % vy    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % wz    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % et    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favavg % Ydot  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               favavg % Ydot_hct  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               favavg % omega ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % omega_hct ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % Da    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % sdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % FI    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % lambda( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favavg % Zm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favavg % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % nu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % kpa   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % ct    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               ! favavg % dm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               ! favavg % tdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               ! favavg % rd    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv , 1:nrv )   , &

               favavg % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % cp    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % W     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favavg % P     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % Tt    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % Ht    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % M     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % cs    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % gamma ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favavg % Pr    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favavg % ca_T_Zm     ( sx-ng:ex+ng , 1:nbins )                               , &
               favavg % ca_c_Zm     ( sx-ng:ex+ng , 1:nbins )                               , &
               favavg % ca_sdr_Zm   ( sx-ng:ex+ng , 1:nbins )                               , &
               favavg % ca_Da_Zm    ( sx-ng:ex+ng , 1:nbins )                               , &
               favavg % ca_omega_Zm ( sx-ng:ex+ng , 1:nbins )                               , &
               favavg % ca_FI_Zm    ( sx-ng:ex+ng , 1:nbins )                               , &

               favavg % ca_omega_FI ( sx-ng:ex+ng , 1:nbins )                               , &

               favavg % ca_omega_Da ( sx-ng:ex+ng , 1:nbins )                               , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate favavg type')

    if ( nss > 0 ) then

       allocate ( favavg % Yass  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nss )         , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate reyavg type 2')

    end if

    if ( LES ) then

       allocate ( favavg % mu_SGS   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , & ! A.Techer
                  favavg % tke_SGS  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

!                  favavg % xiv2  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , & ! SGS variance
!                  favavg % xiv3  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
!                  favavg % xixi  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

                  stat = ok )

       if ( ok > 0 ) call abort_mpi ('error allocate favavg type 3')

    end if


    ! favavg % Y     = 0.0_dp
    ! favavg % ux    = 0.0_dp
    ! favavg % vy    = 0.0_dp
    ! favavg % wz    = 0.0_dp
    ! favavg % et    = 0.0_dp

    ! favavg % Ydot  = 0.0_dp
    ! favavg % omega = 0.0_dp
    ! favavg % Da    = 0.0_dp
    ! favavg % sdr   = 0.0_dp
    ! favavg % FI    = 0.0_dp

    ! favavg % Zm    = 0.0_dp

    ! favavg % mu    = 0.0_dp
    ! favavg % nu    = 0.0_dp
    ! favavg % kpa   = 0.0_dp
    ! favavg % ct    = 0.0_dp
    ! ! favavg % dm    = 0.0_dp
    ! ! favavg % tdr   = 0.0_dp
    ! ! favavg % rd    = 0.0_dp

    ! favavg % T     = 0.0_dp
    ! favavg % H     = 0.0_dp
    ! favavg % cp    = 0.0_dp
    ! favavg % W     = 0.0_dp

    ! favavg % P     = 0.0_dp
    ! favavg % Tt    = 0.0_dp
    ! favavg % Ht    = 0.0_dp
    ! favavg % cs    = 0.0_dp
    ! favavg % gamma = 0.0_dp
    ! favavg % Pr    = 0.0_dp

    ! favavg % ca_T_Zm     = 0.0_dp
    ! favavg % ca_c_Zm     = 0.0_dp
    ! favavg % ca_sdr_Zm   = 0.0_dp
    ! favavg % ca_omega_Zm = 0.0_dp
    ! favavg % ca_FI_Zm    = 0.0_dp

    ! favavg % ca_omega_FI = 0.0_dp

    ! favavg % ca_omega_Da = 0.0_dp


  end subroutine alloc_favavg_type


!> \brief Write Favre average derived type to a file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine write_favavg_type (favavg)


    type (favavg_type) , intent (in)                                 :: favavg !< Favre average derived type


    integer (ip)                 :: ok , i , j , k , l
    character (len_default)      :: number_ascii


    write ( number_ascii , format_restart ) rank


    open ( unit_stat , file   = trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) , &
                       form   = 'unformatted'   ,                                                        &
                       status = 'unknown'       ,                                                        &
                       iostat = ok              )

    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) )


    write (unit_stat) (((( favavg % Y (i,j,k,l) , i = sx , ex      ) , &
                                                  j = sy , ey      ) , &
                                                  k = sz , ez      ) , &
                                                  l = 1  , nrv+npv+nvv )

    write (unit_stat) ((( favavg % ux (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % vy (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % wz (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % et (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) (((( favavg % Ydot (i,j,k,l) , i = sx , ex  ) , &
                                                     j = sy , ey  ) , &
                                                     k = sz , ez  ) , &
                                                     l = 1  , nrv )

    write (unit_stat) (((( favavg % Ydot_hct (i,j,k,l) , i = sx , ex  ) , &
                                                     j = sy , ey  ) , &
                                                     k = sz , ez  ) , &
                                                     l = 1  , nrv )


    write (unit_stat) ((( favavg % omega (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )

    write (unit_stat) ((( favavg % omega_hct (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )


    write (unit_stat) ((( favavg % Da (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % sdr (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

    write (unit_stat) ((( favavg % FI (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % Zm (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % mu (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % nu (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % kpa (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

    write (unit_stat) ((( favavg % ct (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % T (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( favavg % H (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( favavg % cp (i,j,k) , i = sx , ex ) , &
                                                j = sy , ey ) , &
                                                k = sz , ez )

    write (unit_stat) ((( favavg % W (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

    write (unit_stat) ((( favavg % P (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( favavg % Tt (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( favavg % Ht (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( favavg % M (i,j,k) , i = sx , ex ) , &
                                              j = sy , ey ) , &
                                              k = sz , ez )

   write (unit_stat) ((( favavg % cs (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) ((( favavg % gamma (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

   write (unit_stat) ((( favavg % Pr (i,j,k) , i = sx , ex ) , &
                                               j = sy , ey ) , &
                                               k = sz , ez )

   write (unit_stat) (( favavg % ca_T_Zm (i,l) , i = sx , ex   ) , &
                                                 l = 1 , nbins )

   write (unit_stat) (( favavg % ca_c_Zm (i,l) , i = sx , ex   ) , &
                                                 l = 1 , nbins )

   write (unit_stat) (( favavg % ca_sdr_Zm (i,l) , i = sx , ex   ) , &
                                                   l = 1 , nbins )

   write (unit_stat) (( favavg % ca_Da_Zm (i,l) , i = sx , ex   ) , &
                                                  l = 1 , nbins )

   write (unit_stat) (( favavg % ca_omega_Zm (i,l) , i = sx , ex   ) , &
                                                     l = 1 , nbins )

   write (unit_stat) (( favavg % ca_FI_Zm (i,l) , i = sx , ex   ) , &
                                                  l = 1 , nbins )

   write (unit_stat) (( favavg % ca_omega_FI (i,l) , i = sx , ex   ) , &
                                                     l = 1 , nbins )

   write (unit_stat) (( favavg % ca_omega_Da (i,l) , i = sx , ex   ) , &
                                                     l = 1 , nbins )


   close (unit_stat)


  end subroutine write_favavg_type


!> \brief Read Favre average derived from a file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine read_favavg_type ( inp , favavg )


    type (inp_type) , intent (in)                                 :: inp    !< input derived type
    type (favavg_type) , intent (inout)                           :: favavg !< Favre average derived type


    integer (ip)                    :: ok , irank , i , j , k , l

    integer (ip)                    :: px , py , pz , ix , fx , iy , fy , iz , fz

    integer (ip) , dimension (3)    :: mycoords , mydims

    character (len_default)         :: number_ascii


    mydims (1) = inp % nprocz ; mydims (2) = inp % nprocy ; mydims (3) = inp % nprocx
    mycoords (:) = 0
    irank      = 0


    if ( nproc == 1 ) then


       do pz = 1 , inp % nprocz

          iz = ( mycoords(1)*ntz ) / mydims(1) + 1
          fz = ( ( mycoords(1)+1 )*ntz ) / mydims(1)

          mycoords (1) = mycoords (1) + 1
          mycoords (2) = 0
          mycoords (3) = 0

          do py = 1 , inp % nprocy

             iy = ( mycoords(2)*nty ) / mydims(2) + 1
             fy = ( ( mycoords(2)+1 )*nty ) / mydims(2)

             mycoords (2) = mycoords (2) + 1
             mycoords (3) = 0

             do px = 1 , inp % nprocx

                ix = ( mycoords(3)*ntx ) / mydims(3) + 1
                fx = ( ( mycoords(3)+1 )*ntx ) / mydims(3)
                mycoords (3) = mycoords (3) + 1


                write ( number_ascii , format_restart ) irank
                irank = irank + 1

                open ( unit_stat , file   = trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) , &
                                   form   = 'unformatted'   ,                                                        &
                                   status = 'unknown'       ,                                                        &
                                   iostat = ok              )

                if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // &
                    &                   trim (file_favavg) // '_' // trim (number_ascii) )


                read (unit_stat) (((( favavg % Y (i,j,k,l) , i = ix , fx      ) , &
                                                             j = iy , fy      ) , &
                                                             k = iz , fz      ) , &
                                                             l = 1  , nrv+npv+nvv )

                read (unit_stat) ((( favavg % ux (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % vy (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % wz (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % et (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) (((( favavg % Ydot (i,j,k,l) , i = ix , fx  ) , &
                                                                j = iy , fy  ) , &
                                                                k = iz , fz  ) , &
                                                                l = 1  , nrv )

                read (unit_stat) (((( favavg % Ydot_hct (i,j,k,l) , i = ix , fx  ) , &
                                                                j = iy , fy  ) , &
                                                                k = iz , fz  ) , &
                                                                l = 1  , nrv )

                read (unit_stat) ((( favavg % omega (i,j,k) , i = ix , fx ) , &
                                                              j = iy , fy ) , &
                                                              k = iz , fz )

                read (unit_stat) ((( favavg % omega_hct (i,j,k) , i = ix , fx ) , &
                                                              j = iy , fy ) , &
                                                              k = iz , fz )


                read (unit_stat) ((( favavg % Da (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % sdr (i,j,k) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            k = iz , fz )

                read (unit_stat) ((( favavg % FI (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % Zm (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % mu (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % nu (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % kpa (i,j,k) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            k = iz , fz )

                read (unit_stat) ((( favavg % ct (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % T (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( favavg % H (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( favavg % cp (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % W (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( favavg % P (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( favavg % Tt (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % Ht (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % M (i,j,k) , i = ix , fx ) , &
                                                          j = iy , fy ) , &
                                                          k = iz , fz )

                read (unit_stat) ((( favavg % cs (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) ((( favavg % gamma (i,j,k) , i = ix , fx ) , &
                                                              j = iy , fy ) , &
                                                              k = iz , fz )

                read (unit_stat) ((( favavg % Pr (i,j,k) , i = ix , fx ) , &
                                                           j = iy , fy ) , &
                                                           k = iz , fz )

                read (unit_stat) (( favavg % ca_T_Zm (i,l) , i = ix , fx  ) , &
                                                             l = 1 , nbins )

                read (unit_stat) (( favavg % ca_c_Zm (i,l) , i = ix , fx   ) , &
                                                             l = 1 , nbins )

                read (unit_stat) (( favavg % ca_sdr_Zm (i,l) , i = ix , fx   ) , &
                                                               l = 1 , nbins )

                read (unit_stat) (( favavg % ca_Da_Zm (i,l) , i = ix , fx   ) , &
                                                              l = 1 , nbins )

                read (unit_stat) (( favavg % ca_omega_Zm (i,l) , i = ix , fx   ) , &
                                                                 l = 1 , nbins )

                read (unit_stat) (( favavg % ca_FI_Zm (i,l) , i = ix , fx   ) , &
                                                              l = 1 , nbins )

                read (unit_stat) (( favavg % ca_omega_FI (i,l) , i = ix , fx   ) , &
                                                                 l = 1 , nbins )

                read (unit_stat) (( favavg % ca_omega_Da (i,l) , i = ix , fx   ) , &
                                                                 l = 1 , nbins )

             end do
          end do
       end do


    else


       write ( number_ascii , format_restart ) rank

       open ( unit_stat , file   = trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) , &
                          form   = 'unformatted'   ,                                                        &
                          status = 'unknown'       ,                                                        &
                          iostat = ok              )

       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_favavg) // '_' // trim (number_ascii) )


       read (unit_stat) (((( favavg % Y (i,j,k,l) , i = sx , ex      ) , &
                                                    j = sy , ey      ) , &
                                                    k = sz , ez      ) , &
                                                    l = 1  , nrv+npv+nvv )

       read (unit_stat) ((( favavg % ux (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % vy (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % wz (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % et (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) (((( favavg % Ydot (i,j,k,l) , i = sx , ex  ) , &
                                                       j = sy , ey  ) , &
                                                       k = sz , ez  ) , &
                                                       l = 1  , nrv )

       read (unit_stat) (((( favavg % Ydot_hct (i,j,k,l) , i = sx , ex  ) , &
                                                       j = sy , ey  ) , &
                                                       k = sz , ez  ) , &
                                                       l = 1  , nrv )


       read (unit_stat) ((( favavg % omega (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

       read (unit_stat) ((( favavg % omega_hct (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

       read (unit_stat) ((( favavg % Da (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % sdr (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )

       read (unit_stat) ((( favavg % FI (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % Zm (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % mu (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % nu (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % kpa (i,j,k) , i = sx , ex ) , &
                                                   j = sy , ey ) , &
                                                   k = sz , ez )

       read (unit_stat) ((( favavg % ct (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % T (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( favavg % H (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( favavg % cp (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % W (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( favavg % P (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( favavg % Tt (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % Ht (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % M (i,j,k) , i = sx , ex ) , &
                                                 j = sy , ey ) , &
                                                 k = sz , ez )

       read (unit_stat) ((( favavg % cs (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) ((( favavg % gamma (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

       read (unit_stat) ((( favavg % Pr (i,j,k) , i = sx , ex ) , &
                                                  j = sy , ey ) , &
                                                  k = sz , ez )

       read (unit_stat) (( favavg % ca_T_Zm (i,l) , i = sx , ex   ) , &
                                                    l = 1 , nbins )

       read (unit_stat) (( favavg % ca_c_Zm (i,l) , i = sx , ex   ) , &
                                                    l = 1 , nbins )

       read (unit_stat) (( favavg % ca_sdr_Zm (i,l) , i = sx , ex   ) , &
                                                      l = 1 , nbins )

       read (unit_stat) (( favavg % ca_Da_Zm (i,l) , i = sx , ex   ) , &
                                                     l = 1 , nbins )

       read (unit_stat) (( favavg % ca_omega_Zm (i,l) , i = sx , ex   ) , &
                                                        l = 1 , nbins )

       read (unit_stat) (( favavg % ca_FI_Zm (i,l) , i = sx , ex   ) , &
                                                     l = 1 , nbins )

       read (unit_stat) (( favavg % ca_omega_FI (i,l) , i = sx , ex   ) , &
                                                        l = 1 , nbins )

       read (unit_stat) (( favavg % ca_omega_Da (i,l) , i = sx , ex   ) , &
                                                        l = 1 , nbins )


    end if


    close (unit_stat)


  end subroutine read_favavg_type


!> \brief Allocate Reynolds fluctuation derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine alloc_reyfluc_type (reyfluc)


    type (reyfluc_type) , intent (out)                               :: reyfluc !< Reynolds fluctuation derived type


    integer (ip) :: ok


    allocate ( reyfluc % Y     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv+npv+nvv )   , &
               reyfluc % rho   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % ux    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % vy    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % wz    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % et    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyfluc % Ydot  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               reyfluc % Ydot_hct  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               reyfluc % omega ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % omega_hct ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % Da    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % sdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % FI    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyfluc % Zm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyfluc % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % nu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % kpa   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % ct    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               ! reyfluc % dm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               ! reyfluc % tdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               ! reyfluc % rd    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv , 1:nrv )   , &

               reyfluc % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % cp    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % W     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyfluc % P     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % Tt    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % Ht    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % M     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % cs    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % gamma ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               reyfluc % Pr    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               reyfluc % tau ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax ) , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate reyfluc type')

    if ( nss > 0 ) then

       allocate ( reyfluc % Yass  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nss )        , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate reyfluc type 2')

    end if

    ! reyfluc % Y     = 0.0_dp
    ! reyfluc % rho   = 0.0_dp
    ! reyfluc % ux    = 0.0_dp
    ! reyfluc % vy    = 0.0_dp
    ! reyfluc % wz    = 0.0_dp
    ! reyfluc % et    = 0.0_dp

    ! reyfluc % Ydot  = 0.0_dp
    ! reyfluc % omega = 0.0_dp
    ! reyfluc % Da    = 0.0_dp
    ! reyfluc % sdr   = 0.0_dp
    ! reyfluc % FI    = 0.0_dp

    ! reyfluc % Zm    = 0.0_dp

    ! reyfluc % mu    = 0.0_dp
    ! reyfluc % nu    = 0.0_dp
    ! reyfluc % kpa   = 0.0_dp
    ! reyfluc % ct    = 0.0_dp
    ! ! reyfluc % dm    = 0.0_dp
    ! ! reyfluc % tdr   = 0.0_dp
    ! ! reyfluc % rd    = 0.0_dp

    ! reyfluc % T     = 0.0_dp
    ! reyfluc % H     = 0.0_dp
    ! reyfluc % cp    = 0.0_dp
    ! reyfluc % W     = 0.0_dp

    ! reyfluc % P     = 0.0_dp
    ! reyfluc % Tt    = 0.0_dp
    ! reyfluc % Ht    = 0.0_dp
    ! reyfluc % cs    = 0.0_dp
    ! reyfluc % gamma = 0.0_dp
    ! reyfluc % Pr    = 0.0_dp

    ! reyfluc % tau   = 0.0_dp


  end subroutine alloc_reyfluc_type


!> \brief Allocate Favre fluctuation derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine alloc_favfluc_type (favfluc)


    type (favfluc_type) , intent (out)                               :: favfluc !< Favre fluctuation derived type


    integer (ip) :: ok


    allocate ( favfluc % Y     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv+npv+nvv )   , &
               favfluc % ux    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % vy    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % wz    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % et    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favfluc % Ydot  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               favfluc % Ydot_hct  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               favfluc % omega ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % omega_hct ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % Da    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % sdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % FI    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favfluc % Zm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favfluc % mu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % nu    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % kpa   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % ct    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               ! favfluc % dm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               ! favfluc % tdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv )           , &
               ! favfluc % rd    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv , 1:nrv )   , &

               favfluc % T     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % H     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % cp    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % W     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favfluc % P     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % Tt    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % Ht    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % M     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % cs    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % gamma ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &
               favfluc % Pr    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                   , &

               favfluc % ca_T_Zm     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               favfluc % ca_c_Zm     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               favfluc % ca_sdr_Zm   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               favfluc % ca_Da_Zm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               favfluc % ca_omega_Zm ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               favfluc % ca_FI_Zm    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &

               favfluc % ca_omega_FI ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &

               favfluc % ca_omega_Da ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate favfluc type')

    if ( nss > 0 ) then

       allocate ( favfluc % Yass  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nss )        , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate favavg type 2')

    end if

    ! favfluc % Y     = 0.0_dp
    ! favfluc % ux    = 0.0_dp
    ! favfluc % vy    = 0.0_dp
    ! favfluc % wz    = 0.0_dp
    ! favfluc % et    = 0.0_dp

    ! favfluc % Ydot  = 0.0_dp
    ! favfluc % omega = 0.0_dp
    ! favfluc % Da    = 0.0_dp
    ! favfluc % sdr   = 0.0_dp
    ! favfluc % FI    = 0.0_dp

    ! favfluc % Zm    = 0.0_dp

    ! favfluc % mu    = 0.0_dp
    ! favfluc % nu    = 0.0_dp
    ! favfluc % kpa   = 0.0_dp
    ! favfluc % ct    = 0.0_dp
    ! ! favfluc % dm    = 0.0_dp
    ! ! favfluc % tdr   = 0.0_dp
    ! ! favfluc % rd    = 0.0_dp

    ! favfluc % T     = 0.0_dp
    ! favfluc % H     = 0.0_dp
    ! favfluc % cp    = 0.0_dp
    ! favfluc % W     = 0.0_dp

    ! favfluc % P     = 0.0_dp
    ! favfluc % Tt    = 0.0_dp
    ! favfluc % Ht    = 0.0_dp
    ! favfluc % cs    = 0.0_dp
    ! favfluc % gamma = 0.0_dp
    ! favfluc % Pr    = 0.0_dp

    ! favfluc % ca_T_Zm     = 0.0_dp
    ! favfluc % ca_c_Zm     = 0.0_dp
    ! favfluc % ca_sdr_Zm   = 0.0_dp
    ! favfluc % ca_omega_Zm = 0.0_dp
    ! favfluc % ca_FI_Zm    = 0.0_dp

    ! favfluc % ca_omega_FI = 0.0_dp

    ! favfluc % ca_omega_Da = 0.0_dp


  end subroutine alloc_favfluc_type


!> \brief Allocate statistical derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine alloc_stat_type (stat)


    type (stat_type) , intent (out)                                   :: stat !< statistical derived type


    integer (ip) :: ok


    allocate ( stat % reyavg_P      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyavg_Y      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyavg_Zm     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyavg_ux     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyavg_vy     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyavg_wz     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &

               stat % reyavg_uxux      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_vyvy      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_wzwz      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_uxvy      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_uxwz      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_vywz      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyavg_enstrophy ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % favavg_ux_Y      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % favavg_vy_Y      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % favavg_wz_Y      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % favavg_ux_Zm     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % favavg_vy_Zm     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % favavg_wz_Zm     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % reyvar_Y      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv+npv+nvv )     , &
               stat % favvar_Y      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv+npv+nvv )     , &
               stat % reyvar_Zm     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_Zm     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &

               stat % reyvar_ux     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_vy     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_wz     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_uxvy   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_uxwz   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_vywz   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % reyvar_uxvywz ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &

               stat % favvar_ux     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_vy     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_wz     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_uxvy   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_uxwz   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_vywz   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
               stat % favvar_uxvywz ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &

               stat % reyvar_p         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_rho       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_rho_acu   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_rho_ent   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_T         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_T_acu     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_T_ent     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_W         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % favvar_du1_dx1   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_p_du1_dx1 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_rho_p     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % reyvar_rho_W     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % reyvar_omega ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                      , &
               stat % reyvar_Da    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                      , &
               stat % reyvar_sdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                      , &
               stat % reyvar_FI    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                      , &
               stat % favvar_omega ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                      , &
               stat % favvar_Da    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                      , &
               stat % favvar_sdr   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                      , &
               stat % favvar_FI    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                      , &

               stat % favvar_ca_T_Zm     ( sx-ng:ex+ng , 1:nbins )                                  , &
               stat % favvar_ca_c_Zm     ( sx-ng:ex+ng , 1:nbins )                                  , &
               stat % favvar_ca_sdr_Zm   ( sx-ng:ex+ng , 1:nbins )                                  , &
               stat % favvar_ca_Da_Zm    ( sx-ng:ex+ng , 1:nbins )                                  , &
               stat % favvar_ca_omega_Zm ( sx-ng:ex+ng , 1:nbins )                                  , &
               stat % favvar_ca_FI_Zm    ( sx-ng:ex+ng , 1:nbins )                                  , &

               stat % favvar_ca_omega_FI ( sx-ng:ex+ng , 1:nbins )                                  , &

               stat % favvar_ca_omega_Da ( sx-ng:ex+ng , 1:nbins )                                  , &

               stat % tke_dissip       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tke_press_strain ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tke_transp       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )        , &
               stat % tke_ra_ff_ux     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tke_ra_ff_vy     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tke_ra_ff_wz     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % vor_conv         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_conv         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % vor_stret        ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_stret        ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % vor_dila         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_dila         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % vor_baro         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_baro         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % vor_visc         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % ens_visc         ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % scal_var_turb_transp  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )   , &
               stat % scal_var_turb_prod    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )   , &
               stat % scal_var_mol_diff     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )   , &
               stat % scal_var_dissip       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               stat % scal_var_dissip2      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               stat % scal_var_dissip3      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               stat % scal_var_dissip4      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               stat % scal_var_dissip5      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &
               stat % scal_var_dissip6      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )             , &

               stat % tmp_corr_ux      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tmp_corr_vy      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tmp_corr_wz      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tmp_corr_p       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tmp_corr_Y       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % tmp_corr_Zm      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % rey_dissip_11    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_22    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_33    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_12    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_13    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &
               stat % rey_dissip_23    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                  , &

               stat % rey_press_strain_11 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_22 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_33 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_12 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_13 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &
               stat % rey_press_strain_23 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )               , &

               stat % rey_transp_11 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_22 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_33 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_12 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_13 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
               stat % rey_transp_23 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &

               stat % pdf_Y   ( sx-ng:ex+ng , sy-ng:ey+ng , 1:nbins )                               , &
               stat % pdf_Zm  ( sx-ng:ex+ng , sy-ng:ey+ng , 1:nbins )                               , &
               stat % pdf_c   ( sx-ng:ex+ng , sy-ng:ey+ng , 1:nbins )                               , &
               stat % pdf_sdr ( sx-ng:ex+ng , sy-ng:ey+ng , 1:nbins )                               , &
               stat % pdf_Da  ( sx-ng:ex+ng , sy-ng:ey+ng , 1:nbins )                               , &
               stat % pdf_FI  ( sx-ng:ex+ng , sy-ng:ey+ng , 1:nbins )                               , &

               stat % cpdf_omega_FI ( sx-ng:ex+ng , 1:nbins )                                       , &
               stat % cpdf_omega_Da ( sx-ng:ex+ng , 1:nbins )                                       , &

               stat % spec_ux  ( sz-ng:ez+ng , 1:ntimes )                                           , &
               stat % spec_vy  ( sz-ng:ez+ng , 1:ntimes )                                           , &
               stat % spec_wz  ( sz-ng:ez+ng , 1:ntimes )                                           , &
               stat % spec_Y   ( sz-ng:ez+ng , 1:ntimes )                                           , &
               stat % spec_Zm  ( sz-ng:ez+ng , 1:ntimes )                                           , &
               stat % spec_sdr ( sz-ng:ez+ng , 1:ntimes )                                           , &

               stat = ok )

    if ( ok > 0 ) call abort_mpi ('error allocate stat type')

    if ( nss > 0 ) then

       allocate ( stat % reyvar_Yass  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nss )             , &
                  stat % favvar_Yass  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nss )             , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate stat type 2')

    end if

    ! stat % reyavg_P         = 0.0_dp
    ! stat % reyavg_Y         = 0.0_dp
    ! stat % reyavg_Zm        = 0.0_dp
    ! stat % reyavg_ux        = 0.0_dp
    ! stat % reyavg_vy        = 0.0_dp
    ! stat % reyavg_wz        = 0.0_dp

    ! stat % reyavg_uxux      = 0.0_dp
    ! stat % reyavg_vyvy      = 0.0_dp
    ! stat % reyavg_wzwz      = 0.0_dp
    ! stat % reyavg_uxvy      = 0.0_dp
    ! stat % reyavg_uxwz      = 0.0_dp
    ! stat % reyavg_vywz      = 0.0_dp
    ! stat % reyavg_enstrophy = 0.0_dp

    ! stat % favavg_ux_Y      = 0.0_dp
    ! stat % favavg_vy_Y      = 0.0_dp
    ! stat % favavg_wz_Y      = 0.0_dp
    ! stat % favavg_ux_Zm     = 0.0_dp
    ! stat % favavg_vy_Zm     = 0.0_dp
    ! stat % favavg_wz_Zm     = 0.0_dp

    ! stat % reyvar_Y         = 0.0_dp
    ! stat % favvar_Y         = 0.0_dp
    ! stat % reyvar_Zm        = 0.0_dp
    ! stat % favvar_Zm        = 0.0_dp

    ! stat % reyvar_ux        = 0.0_dp
    ! stat % reyvar_vy        = 0.0_dp
    ! stat % reyvar_wz        = 0.0_dp
    ! stat % reyvar_uxvy      = 0.0_dp
    ! stat % reyvar_uxwz      = 0.0_dp
    ! stat % reyvar_vywz      = 0.0_dp
    ! stat % reyvar_uxvywz    = 0.0_dp

    ! stat % favvar_ux        = 0.0_dp
    ! stat % favvar_vy        = 0.0_dp
    ! stat % favvar_wz        = 0.0_dp
    ! stat % favvar_uxvy      = 0.0_dp
    ! stat % favvar_uxwz      = 0.0_dp
    ! stat % favvar_vywz      = 0.0_dp
    ! stat % favvar_uxvywz    = 0.0_dp

    ! stat % reyvar_p         = 0.0_dp
    ! stat % reyvar_rho       = 0.0_dp
    ! stat % reyvar_rho_acu   = 0.0_dp
    ! stat % reyvar_rho_ent   = 0.0_dp
    ! stat % reyvar_T         = 0.0_dp
    ! stat % reyvar_T_acu     = 0.0_dp
    ! stat % reyvar_T_ent     = 0.0_dp
    ! stat % reyvar_W         = 0.0_dp
    ! stat % favvar_du1_dx1   = 0.0_dp
    ! stat % reyvar_p_du1_dx1 = 0.0_dp
    ! stat % reyvar_rho_p     = 0.0_dp
    ! stat % reyvar_rho_W     = 0.0_dp

    ! stat % reyvar_omega     = 0.0_dp
    ! stat % reyvar_Da        = 0.0_dp
    ! stat % reyvar_sdr       = 0.0_dp
    ! stat % reyvar_FI        = 0.0_dp
    ! stat % favvar_omega     = 0.0_dp
    ! stat % favvar_Da        = 0.0_dp
    ! stat % favvar_sdr       = 0.0_dp
    ! stat % favvar_FI        = 0.0_dp

    ! stat % favvar_ca_T_Zm     = 0.0_dp
    ! stat % favvar_ca_c_Zm     = 0.0_dp
    ! stat % favvar_ca_sdr_Zm   = 0.0_dp
    ! stat % favvar_ca_Da_Zm    = 0.0_dp
    ! stat % favvar_ca_omega_Zm = 0.0_dp
    ! stat % favvar_ca_FI_Zm    = 0.0_dp

    ! stat % favvar_ca_omega_FI = 0.0_dp

    ! stat % favvar_ca_omega_Da = 0.0_dp

    ! stat % tke_dissip       = 0.0_dp
    ! stat % tke_press_strain = 0.0_dp
    ! stat % tke_transp       = 0.0_dp
    ! stat % tke_ra_ff_ux     = 0.0_dp
    ! stat % tke_ra_ff_vy     = 0.0_dp
    ! stat % tke_ra_ff_wz     = 0.0_dp

    ! stat % vor_conv         = 0.0_dp
    ! stat % ens_conv         = 0.0_dp
    ! stat % vor_stret        = 0.0_dp
    ! stat % ens_stret        = 0.0_dp
    ! stat % vor_dila         = 0.0_dp
    ! stat % ens_dila         = 0.0_dp
    ! stat % vor_baro         = 0.0_dp
    ! stat % ens_baro         = 0.0_dp
    ! stat % vor_visc         = 0.0_dp
    ! stat % ens_visc         = 0.0_dp

    ! stat % scal_var_turb_transp = 0.0_dp
    ! stat % scal_var_turb_prod   = 0.0_dp
    ! stat % scal_var_mol_diff    = 0.0_dp
    ! stat % scal_var_dissip      = 0.0_dp
    ! stat % scal_var_dissip2     = 0.0_dp
    ! stat % scal_var_dissip3     = 0.0_dp
    ! stat % scal_var_dissip4     = 0.0_dp
    ! stat % scal_var_dissip5     = 0.0_dp
    ! stat % scal_var_dissip6     = 0.0_dp

    ! stat % tmp_corr_ux      = 0.0_dp
    ! stat % tmp_corr_vy      = 0.0_dp
    ! stat % tmp_corr_wz      = 0.0_dp
    ! stat % tmp_corr_p       = 0.0_dp
    ! stat % tmp_corr_Y       = 0.0_dp
    ! stat % tmp_corr_Zm      = 0.0_dp

    ! stat % rey_dissip_11 = 0.0_dp       ; stat % rey_dissip_22 = 0.0_dp       ; stat % rey_dissip_33 = 0.0_dp
    ! stat % rey_dissip_12 = 0.0_dp       ; stat % rey_dissip_13 = 0.0_dp       ; stat % rey_dissip_23 = 0.0_dp
    ! stat % rey_press_strain_11 = 0.0_dp ; stat % rey_press_strain_22 = 0.0_dp ; stat % rey_press_strain_33 = 0.0_dp
    ! stat % rey_press_strain_12 = 0.0_dp ; stat % rey_press_strain_13 = 0.0_dp ; stat % rey_press_strain_23 = 0.0_dp
    ! stat % rey_transp_11 = 0.0_dp       ; stat % rey_transp_22 = 0.0_dp       ; stat % rey_transp_33 = 0.0_dp
    ! stat % rey_transp_12 = 0.0_dp       ; stat % rey_transp_13 = 0.0_dp       ; stat % rey_transp_23 = 0.0_dp

    ! stat % pdf_Y            = 0.0_dp
    ! stat % pdf_Zm           = 0.0_dp
    ! stat % pdf_c            = 0.0_dp
    ! stat % pdf_sdr          = 0.0_dp
    ! stat % pdf_Da           = 0.0_dp
    ! stat % pdf_FI           = 0.0_dp

    ! stat % cpdf_omega_FI    = 0.0_dp
    ! stat % cpdf_omega_Da    = 0.0_dp

    ! stat % spec_ux          = 0.0_dp
    ! stat % spec_vy          = 0.0_dp
    ! stat % spec_wz          = 0.0_dp
    ! stat % spec_Y           = 0.0_dp
    ! stat % spec_Zm          = 0.0_dp
    ! stat % spec_sdr         = 0.0_dp


  end subroutine alloc_stat_type


!> \brief Write statistical derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine write_stat_type (stat)


    type (stat_type) , intent (in)                                :: stat !< statistical derived type


    integer (ip)                 :: ok , i , j , k , l
    character (len_default)      :: number_ascii


    write ( number_ascii , format_restart ) rank


    open ( unit_stat , file   = trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) , &
                       form   = 'unformatted'   ,                                                      &
                       status = 'unknown'       ,                                                      &
                       iostat = ok              )

    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) )


    write (unit_stat) ((( stat % reyavg_P (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % reyavg_Y (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % reyavg_Zm (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyavg_ux (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyavg_vy (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyavg_wz (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyavg_uxux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_vyvy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_wzwz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_uxvy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_uxwz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_vywz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyavg_enstrophy (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) ((( stat % favavg_ux_Y (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % favavg_vy_Y (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % favavg_wz_Y (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % favavg_ux_Zm (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % favavg_vy_Zm (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % favavg_wz_Zm (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) (((( stat % reyvar_Y (i,j,k,l) , i = sx , ex     ) , &
                                                       j = sy , ey     ) , &
                                                       k = sz , ez     ) , &
                                                       l = 1 , nrv+npv+nvv )

    write (unit_stat) (((( stat % favvar_Y (i,j,k,l) , i = sx , ex     ) , &
                                                       j = sy , ey     ) , &
                                                       k = sz , ez     ) , &
                                                       l = 1 , nrv+npv+nvv )

    write (unit_stat) ((( stat % reyvar_Zm (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % favvar_Zm (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyvar_ux (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyvar_vy (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyvar_wz (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyvar_uxvy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyvar_uxwz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyvar_vywz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % reyvar_uxvywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % favvar_ux (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % favvar_vy (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % favvar_wz (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % favvar_uxvy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % favvar_uxwz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % favvar_vywz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % favvar_uxvywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % reyvar_p (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % reyvar_rho (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

    write (unit_stat) ((( stat % reyvar_rho_acu (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

    write (unit_stat) ((( stat % reyvar_rho_ent (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

    write (unit_stat) ((( stat % reyvar_T (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % reyvar_T_acu (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % reyvar_T_ent (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % reyvar_W (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % favvar_du1_dx1 (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

    write (unit_stat) ((( stat % reyvar_p_du1_dx1 (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) ((( stat % reyvar_rho_p (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % reyvar_rho_W (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % reyvar_omega (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % reyvar_Da (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % reyvar_sdr (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

    write (unit_stat) ((( stat % reyvar_FI (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % favvar_omega (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % favvar_Da (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % favvar_sdr (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

    write (unit_stat) ((( stat % favvar_FI (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) (( stat % favvar_ca_T_Zm (i,l) , i = sx , ex ) , &
                                                       l = 1 , nbins )

    write (unit_stat) (( stat % favvar_ca_c_Zm (i,l) , i = sx , ex ) , &
                                                       l = 1 , nbins )

    write (unit_stat) (( stat % favvar_ca_sdr_Zm (i,l) , i = sx , ex ) , &
                                                         l = 1 , nbins )

    write (unit_stat) (( stat % favvar_ca_Da_Zm (i,l) , i = sx , ex ) , &
                                                        l = 1 , nbins )

    write (unit_stat) (( stat % favvar_ca_omega_Zm (i,l) , i = sx , ex ) , &
                                                           l = 1 , nbins )

    write (unit_stat) (( stat % favvar_ca_FI_Zm (i,l) , i = sx , ex ) , &
                                                        l = 1 , nbins )

    write (unit_stat) (( stat % favvar_ca_omega_FI (i,l) , i = sx , ex ) , &
                                                           l = 1 , nbins )

    write (unit_stat) (( stat % favvar_ca_omega_Da (i,l) , i = sx , ex ) , &
                                                           l = 1 , nbins )

    write (unit_stat) ((( stat % tke_dissip (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

    write (unit_stat) ((( stat % tke_press_strain (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) (((( stat % tke_transp (i,j,k,l) , i = sx , ex      ) , &
                                                         j = sy , ey      ) , &
                                                         k = sz , ez      ) , &
                                                         l = 1  , ndimmax )

    write (unit_stat) ((( stat % tke_ra_ff_ux (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % tke_ra_ff_vy (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % tke_ra_ff_wz (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

    write (unit_stat) ((( stat % vor_conv (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % ens_conv (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % vor_stret (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % ens_stret (i,j,k) , i = sx , ex ) , &
                                                     j = sy , ey ) , &
                                                     k = sz , ez )

    write (unit_stat) ((( stat % vor_dila (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % ens_dila (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % vor_baro (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % ens_baro (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % vor_visc (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) ((( stat % ens_visc (i,j,k) , i = sx , ex ) , &
                                                    j = sy , ey ) , &
                                                    k = sz , ez )

    write (unit_stat) (((( stat % scal_var_turb_transp (i,j,k,l) , i = sx , ex      ) , &
                                                                   j = sy , ey      ) , &
                                                                   k = sz , ez      ) , &
                                                                   l = 1  , ndimmax )

    write (unit_stat) (((( stat % scal_var_turb_prod (i,j,k,l) , i = sx , ex      ) , &
                                                                 j = sy , ey      ) , &
                                                                 k = sz , ez      ) , &
                                                                 l = 1  , ndimmax )

    write (unit_stat) (((( stat % scal_var_mol_diff (i,j,k,l) , i = sx , ex      ) , &
                                                                j = sy , ey      ) , &
                                                                k = sz , ez      ) , &
                                                                l = 1  , ndimmax )

    write (unit_stat) ((( stat % scal_var_dissip (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

    write (unit_stat) ((( stat % scal_var_dissip2 (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) ((( stat % scal_var_dissip3 (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) ((( stat % scal_var_dissip4 (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) ((( stat % scal_var_dissip5 (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) ((( stat % scal_var_dissip6 (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_ux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_vy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_wz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_p (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_Y (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

    write (unit_stat) ((( stat % tmp_corr_Zm (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_11 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_22 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_33 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_12 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_13 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_dissip_23 (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_11 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_22 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_33 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_12 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_13 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) ((( stat % rey_press_strain_23 (i,j,k) , i = sx , ex ) , &
                                                               j = sy , ey ) , &
                                                               k = sz , ez )

    write (unit_stat) (((( stat % rey_transp_11 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_22 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_33 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_12 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_13 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    write (unit_stat) (((( stat % rey_transp_23 (i,j,k,l) , i = sx , ex      ) , &
                                                            j = sy , ey      ) , &
                                                            k = sz , ez      ) , &
                                                            l = 1  , ndimmax )

    ! pdfs and spectras do not need to be saved (not plotted in paraview)

    write (unit_stat) ((( stat % pdf_Y (i,j,l) , i = sx , ex   ) , &
                                                 j = sy , ey   ) , &
                                                 l = 1 , nbins )

    write (unit_stat) ((( stat % pdf_Zm (i,j,l) , i = sx , ex   ) , &
                                                  j = sy , ey   ) , &
                                                  l = 1 , nbins )

    write (unit_stat) ((( stat % pdf_c (i,j,l) , i = sx , ex   ) , &
                                                 j = sy , ey   ) , &
                                                 l = 1 , nbins )

    write (unit_stat) ((( stat % pdf_sdr (i,j,l) , i = sx , ex   ) , &
                                                   j = sy , ey   ) , &
                                                   l = 1 , nbins )

    write (unit_stat) ((( stat % pdf_Da (i,j,l) , i = sx , ex   ) , &
                                                  j = sy , ey   ) , &
                                                  l = 1 , nbins )

    write (unit_stat) ((( stat % pdf_FI (i,j,l) , i = sx , ex   ) , &
                                                  j = sy , ey   ) , &
                                                  l = 1 , nbins )

    write (unit_stat) (( stat % cpdf_omega_FI (i,l) , i = sx , ex   ) , &
                                                      l = 1 , nbins )

    write (unit_stat) (( stat % cpdf_omega_Da (i,l) , i = sx , ex   ) , &
                                                      l = 1 , nbins )

    write (unit_stat) (( stat % spec_ux (k,l) , k = sz , ez    ) , &
                                                l = 1 , ntimes )

    write (unit_stat) (( stat % spec_vy (k,l) , k = sz , ez    ) , &
                                                l = 1 , ntimes )

    write (unit_stat) (( stat % spec_wz (k,l) , k = sz , ez    ) , &
                                                l = 1 , ntimes )

    write (unit_stat) (( stat % spec_Y (k,l) , k = sz , ez    ) , &
                                               l = 1 , ntimes )

    write (unit_stat) (( stat % spec_Zm (k,l) , k = sz , ez    ) , &
                                                l = 1 , ntimes )

    write (unit_stat) (( stat % spec_sdr (k,l) , k = sz , ez    ) , &
                                                 l = 1 , ntimes )


    close (unit_stat)


  end subroutine write_stat_type


!> \brief Read statistical derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine read_stat_type ( inp , stat )


    type (inp_type) , intent (in)                                 :: inp  !< input derived type
    type (stat_type) , intent (inout)                             :: stat !< statistical derived type


    integer (ip)                    :: ok , irank , i , j , k , l

    integer (ip)                    :: px , py , pz , ix , fx , iy , fy , iz , fz

    integer (ip) , dimension (3)    :: mycoords , mydims

    character (len_default)         :: number_ascii


    mydims (1) = inp % nprocz ; mydims (2) = inp % nprocy ; mydims (3) = inp % nprocx
    mycoords (:) = 0
    irank        = 0


    if ( nproc == 1 ) then


       do pz = 1 , inp % nprocz

          iz = ( mycoords(1)*ntz ) / mydims(1) + 1
          fz = ( ( mycoords(1)+1 )*ntz ) / mydims(1)

          mycoords (1) = mycoords (1) + 1
          mycoords (2) = 0
          mycoords (3) = 0

          do py = 1 , inp % nprocy

             iy = ( mycoords(2)*nty ) / mydims(2) + 1
             fy = ( ( mycoords(2)+1 )*nty ) / mydims(2)

             mycoords (2) = mycoords (2) + 1
             mycoords (3) = 0

             do px = 1 , inp % nprocx

                ix = ( mycoords(3)*ntx ) / mydims(3) + 1
                fx = ( ( mycoords(3)+1 )*ntx ) / mydims(3)
                mycoords (3) = mycoords (3) + 1


                write ( number_ascii , format_restart ) irank
                irank = irank + 1

                open ( unit_stat , file   = trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) , &
                                   form   = 'unformatted'   ,                                                      &
                                   status = 'unknown'       ,                                                      &
                                   iostat = ok              )

                if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // &
                    &                     trim (file_stat) // '_' // trim (number_ascii) )


                read (unit_stat) ((( stat % reyavg_P (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % reyavg_Y (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % reyavg_Zm (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyavg_ux (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyavg_vy (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyavg_wz (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyavg_uxux (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_vyvy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_wzwz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_uxvy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_uxwz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_vywz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyavg_enstrophy (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) ((( stat % favavg_ux_Y (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % favavg_vy_Y (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % favavg_wz_Y (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % favavg_ux_Zm (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % favavg_vy_Zm (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % favavg_wz_Zm (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) (((( stat % reyvar_Y (i,j,k,l) , i = ix , fx     ) , &
                                                                  j = iy , fy     ) , &
                                                                  k = iz , fz     ) , &
                                                                  l = 1 , nrv+npv+nvv )

                read (unit_stat) (((( stat % favvar_Y (i,j,k,l) , i = ix , fx     ) , &
                                                                  j = iy , fy     ) , &
                                                                  k = iz , fz     ) , &
                                                                  l = 1 , nrv+npv+nvv )

                read (unit_stat) ((( stat % reyvar_Zm (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % favvar_Zm (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyvar_ux (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyvar_vy (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyvar_wz (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyvar_uxvy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyvar_uxwz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyvar_vywz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % reyvar_uxvywz (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % favvar_ux (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % favvar_vy (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % favvar_wz (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % favvar_uxvy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % favvar_uxwz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % favvar_vywz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % favvar_uxvywz (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % reyvar_p (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % reyvar_rho (i,j,k) , i = ix , fx ) , &
                                                                 j = iy , fy ) , &
                                                                 k = iz , fz )

                read (unit_stat) ((( stat % reyvar_rho_acu (i,j,k) , i = ix , fx ) , &
                                                                     j = iy , fy ) , &
                                                                     k = iz , fz )

                read (unit_stat) ((( stat % reyvar_rho_ent (i,j,k) , i = ix , fx ) , &
                                                                     j = iy , fy ) , &
                                                                     k = iz , fz )

                read (unit_stat) ((( stat % reyvar_T (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % reyvar_T_acu (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % reyvar_T_ent (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % reyvar_W (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % favvar_du1_dx1 (i,j,k) , i = ix , fx ) , &
                                                                     j = iy , fy ) , &
                                                                     k = iz , fz )

                read (unit_stat) ((( stat % reyvar_p_du1_dx1 (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) ((( stat % reyvar_rho_p (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % reyvar_rho_W (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % reyvar_omega (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % reyvar_Da (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % reyvar_sdr (i,j,k) , i = ix , fx ) , &
                                                                 j = iy , fy ) , &
                                                                 k = iz , fz )

                read (unit_stat) ((( stat % reyvar_FI (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % favvar_omega (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % favvar_Da (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % favvar_sdr (i,j,k) , i = ix , fx ) , &
                                                                 j = iy , fy ) , &
                                                                 k = iz , fz )

                read (unit_stat) ((( stat % favvar_FI (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) (( stat % favvar_ca_T_Zm (i,l) , i = ix , fx   ) , &
                                                                  l = 1 , nbins )

                read (unit_stat) (( stat % favvar_ca_c_Zm (i,l) , i = ix , fx   ) , &
                                                                  l = 1 , nbins )

                read (unit_stat) (( stat % favvar_ca_sdr_Zm (i,l) , i = ix , fx   ) , &
                                                                    l = 1 , nbins )

                read (unit_stat) (( stat % favvar_ca_Da_Zm (i,l) , i = ix , fx   ) , &
                                                                   l = 1 , nbins )

                read (unit_stat) (( stat % favvar_ca_omega_Zm (i,l) , i = ix , fx   ) , &
                                                                     l = 1 , nbins )

                read (unit_stat) (( stat % favvar_ca_FI_Zm (i,l) , i = ix , fx   ) , &
                                                                   l = 1 , nbins )

                read (unit_stat) (( stat % favvar_ca_omega_FI (i,l) , i = ix , fx   ) , &
                                                                      l = 1 , nbins )

                read (unit_stat) (( stat % favvar_ca_omega_Da (i,l) , i = ix , fx   ) , &
                                                                      l = 1 , nbins )

                read (unit_stat) ((( stat % tke_dissip (i,j,k) , i = ix , fx ) , &
                                                                 j = iy , fy ) , &
                                                                 k = iz , fz )

                read (unit_stat) ((( stat % tke_press_strain (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) (((( stat % tke_transp (i,j,k,l) , i = ix , fx      ) , &
                                                                    j = iy , fy      ) , &
                                                                    k = iz , fz      ) , &
                                                                    l = 1  , ndimmax )

                read (unit_stat) ((( stat % tke_ra_ff_ux (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % tke_ra_ff_vy (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % tke_ra_ff_wz (i,j,k) , i = ix , fx ) , &
                                                                   j = iy , fy ) , &
                                                                   k = iz , fz )

                read (unit_stat) ((( stat % vor_conv (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % ens_conv (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % vor_stret (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % ens_stret (i,j,k) , i = ix , fx ) , &
                                                                j = iy , fy ) , &
                                                                k = iz , fz )

                read (unit_stat) ((( stat % vor_dila (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % ens_dila (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % vor_baro (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % ens_baro (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % vor_visc (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) ((( stat % ens_visc (i,j,k) , i = ix , fx ) , &
                                                               j = iy , fy ) , &
                                                               k = iz , fz )

                read (unit_stat) (((( stat % scal_var_turb_transp (i,j,k,l) , i = ix , fx      ) , &
                                                                              j = iy , fy      ) , &
                                                                              k = iz , fz      ) , &
                                                                              l = 1  , ndimmax )

                read (unit_stat) (((( stat % scal_var_turb_prod (i,j,k,l) , i = ix , fx      ) , &
                                                                            j = iy , fy      ) , &
                                                                            k = iz , fz      ) , &
                                                                            l = 1  , ndimmax )

                read (unit_stat) (((( stat % scal_var_mol_diff (i,j,k,l) , i = ix , fx      ) , &
                                                                           j = iy , fy      ) , &
                                                                           k = iz , fz      ) , &
                                                                           l = 1  , ndimmax )

                read (unit_stat) ((( stat % scal_var_dissip (i,j,k) , i = ix , fx ) , &
                                                                      j = iy , fy ) , &
                                                                      k = iz , fz )

                read (unit_stat) ((( stat % scal_var_dissip2 (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) ((( stat % scal_var_dissip3 (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) ((( stat % scal_var_dissip4 (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) ((( stat % scal_var_dissip5 (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) ((( stat % scal_var_dissip6 (i,j,k) , i = ix , fx ) , &
                                                                       j = iy , fy ) , &
                                                                       k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_ux (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_vy (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_wz (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_p (i,j,k) , i = ix , fx ) , &
                                                                 j = iy , fy ) , &
                                                                 k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_Y (i,j,k) , i = ix , fx ) , &
                                                                 j = iy , fy ) , &
                                                                 k = iz , fz )

                read (unit_stat) ((( stat % tmp_corr_Zm (i,j,k) , i = ix , fx ) , &
                                                                  j = iy , fy ) , &
                                                                  k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_11 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_22 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_33 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_12 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_13 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_dissip_23 (i,j,k) , i = ix , fx ) , &
                                                                    j = iy , fy ) , &
                                                                    k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_11 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_22 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_33 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_12 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_13 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) ((( stat % rey_press_strain_23 (i,j,k) , i = ix , fx ) , &
                                                                          j = iy , fy ) , &
                                                                          k = iz , fz )

                read (unit_stat) (((( stat % rey_transp_11 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_22 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_33 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_12 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_13 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) (((( stat % rey_transp_23 (i,j,k,l) , i = ix , fx      ) , &
                                                                       j = iy , fy      ) , &
                                                                       k = iz , fz      ) , &
                                                                       l = 1  , ndimmax )

                read (unit_stat) ((( stat % pdf_Y (i,j,l) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            l = 1 , nbins )

                read (unit_stat) ((( stat % pdf_Zm (i,j,l) , i = ix , fx ) , &
                                                             j = iy , fy ) , &
                                                             l = 1 , nbins )

                read (unit_stat) ((( stat % pdf_c (i,j,l) , i = ix , fx ) , &
                                                            j = iy , fy ) , &
                                                            l = 1 , nbins )

                read (unit_stat) ((( stat % pdf_sdr (i,j,l) , i = ix , fx ) , &
                                                              j = iy , fy ) , &
                                                              l = 1 , nbins )

                read (unit_stat) ((( stat % pdf_Da (i,j,l) , i = ix , fx ) , &
                                                             j = iy , fy ) , &
                                                             l = 1 , nbins )

                read (unit_stat) ((( stat % pdf_FI (i,j,l) , i = ix , fx ) , &
                                                             j = iy , fy ) , &
                                                             l = 1 , nbins )

                read (unit_stat) (( stat % cpdf_omega_FI (i,l) , i = ix , fx ) , &
                                                                 l = 1 , nbins )

                read (unit_stat) (( stat % cpdf_omega_Da (i,l) , i = ix , fx ) , &
                                                                 l = 1 , nbins )

                read (unit_stat) (( stat % spec_ux (k,l) , k = iz , fz ) , &
                                                           l = 1 , ntimes )

                read (unit_stat) (( stat % spec_vy (k,l) , k = iz , fz ) , &
                                                           l = 1 , ntimes )

                read (unit_stat) (( stat % spec_wz (k,l) , k = iz , fz ) , &
                                                           l = 1 , ntimes )

                read (unit_stat) (( stat % spec_Y (k,l) , k = iz , fz ) , &
                                                          l = 1 , ntimes )

                read (unit_stat) (( stat % spec_Zm (k,l) , k = iz , fz ) , &
                                                           l = 1 , ntimes )

                read (unit_stat) (( stat % spec_sdr (k,l) , k = iz , fz ) , &
                                                            l = 1 , ntimes )


             end do
          end do
       end do


    else ! parallel problem


       write ( number_ascii , format_restart ) rank

       open ( unit_stat , file   = trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) , &
                          form   = 'unformatted'   ,                                                      &
                          status = 'unknown'       ,                                                      &
                          iostat = ok              )

       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_statvar) // trim (file_stat) // '_' // trim (number_ascii) )


       read (unit_stat) ((( stat % reyavg_P (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % reyavg_Y (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % reyavg_Zm (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyavg_ux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyavg_vy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyavg_wz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyavg_uxux (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_vyvy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_wzwz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_uxvy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_uxwz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_vywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyavg_enstrophy (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) ((( stat % favavg_ux_Y (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % favavg_vy_Y (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % favavg_wz_Y (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % favavg_ux_Zm (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % favavg_vy_Zm (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % favavg_wz_Zm (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) (((( stat % reyvar_Y (i,j,k,l) , i = sx , ex     ) , &
                                                         j = sy , ey     ) , &
                                                         k = sz , ez     ) , &
                                                         l = 1 , nrv+npv+nvv )

       read (unit_stat) (((( stat % favvar_Y (i,j,k,l) , i = sx , ex     ) , &
                                                         j = sy , ey     ) , &
                                                         k = sz , ez     ) , &
                                                         l = 1 , nrv+npv+nvv )

       read (unit_stat) ((( stat % reyvar_Zm (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % favvar_Zm (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyvar_ux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyvar_vy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyvar_wz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyvar_uxvy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyvar_uxwz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyvar_vywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % reyvar_uxvywz (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % favvar_ux (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % favvar_vy (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % favvar_wz (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % favvar_uxvy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % favvar_uxwz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % favvar_vywz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % favvar_uxvywz (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % reyvar_p (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % reyvar_rho (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

       read (unit_stat) ((( stat % reyvar_rho_acu (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

       read (unit_stat) ((( stat % reyvar_rho_ent (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

       read (unit_stat) ((( stat % reyvar_T (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % reyvar_T_acu (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % reyvar_T_ent (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % reyvar_W (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % favvar_du1_dx1 (i,j,k) , i = sx , ex ) , &
                                                            j = sy , ey ) , &
                                                            k = sz , ez )

       read (unit_stat) ((( stat % reyvar_p_du1_dx1 (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) ((( stat % reyvar_rho_p (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % reyvar_rho_W (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % reyvar_omega (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % reyvar_Da (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % reyvar_sdr (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

       read (unit_stat) ((( stat % reyvar_FI (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % favvar_omega (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % favvar_Da (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % favvar_sdr (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

       read (unit_stat) ((( stat % favvar_FI (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) (( stat % favvar_ca_T_Zm (i,l) , i = sx , ex   ) , &
                                                         l = 1 , nbins )

       read (unit_stat) (( stat % favvar_ca_c_Zm (i,l) , i = sx , ex   ) , &
                                                         l = 1 , nbins )

       read (unit_stat) (( stat % favvar_ca_sdr_Zm (i,l) , i = sx , ex   ) , &
                                                           l = 1 , nbins )

       read (unit_stat) (( stat % favvar_ca_Da_Zm (i,l) , i = sx , ex   ) , &
                                                          l = 1 , nbins )

       read (unit_stat) (( stat % favvar_ca_omega_Zm (i,l) , i = sx , ex   ) , &
                                                             l = 1 , nbins )

       read (unit_stat) (( stat % favvar_ca_FI_Zm (i,l) , i = sx , ex   ) , &
                                                          l = 1 , nbins )

       read (unit_stat) (( stat % favvar_ca_omega_FI (i,l) , i = sx , ex   ) , &
                                                             l = 1 , nbins )

       read (unit_stat) (( stat % favvar_ca_omega_Da (i,l) , i = sx , ex   ) , &
                                                             l = 1 , nbins )

       read (unit_stat) ((( stat % tke_dissip (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

       read (unit_stat) ((( stat % tke_press_strain (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) (((( stat % tke_transp (i,j,k,l) , i = sx , ex      ) , &
                                                           j = sy , ey      ) , &
                                                           k = sz , ez      ) , &
                                                           l = 1  , ndimmax )

       read (unit_stat) ((( stat % tke_ra_ff_ux (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % tke_ra_ff_vy (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % tke_ra_ff_wz (i,j,k) , i = sx , ex ) , &
                                                          j = sy , ey ) , &
                                                          k = sz , ez )

       read (unit_stat) ((( stat % vor_conv (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % ens_conv (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % vor_stret (i,j,k) , i = sx , ex ) , &
                                                       j = sy , ey ) , &
                                                       k = sz , ez )

       read (unit_stat) ((( stat % ens_stret (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % vor_dila (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % ens_dila (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % vor_baro (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % ens_baro (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % vor_visc (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) ((( stat % ens_visc (i,j,k) , i = sx , ex ) , &
                                                      j = sy , ey ) , &
                                                      k = sz , ez )

       read (unit_stat) (((( stat % scal_var_turb_transp (i,j,k,l) , i = sx , ex      ) , &
                                                                     j = sy , ey      ) , &
                                                                     k = sz , ez      ) , &
                                                                     l = 1  , ndimmax )

       read (unit_stat) (((( stat % scal_var_turb_prod (i,j,k,l) , i = sx , ex      ) , &
                                                                   j = sy , ey      ) , &
                                                                   k = sz , ez      ) , &
                                                                   l = 1  , ndimmax )

       read (unit_stat) (((( stat % scal_var_mol_diff (i,j,k,l) , i = sx , ex      ) , &
                                                                  j = sy , ey      ) , &
                                                                  k = sz , ez      ) , &
                                                                  l = 1  , ndimmax )

       read (unit_stat) ((( stat % scal_var_dissip (i,j,k) , i = sx , ex ) , &
                                                             j = sy , ey ) , &
                                                             k = sz , ez )

       read (unit_stat) ((( stat % scal_var_dissip2 (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) ((( stat % scal_var_dissip3 (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) ((( stat % scal_var_dissip4 (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) ((( stat % scal_var_dissip5 (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) ((( stat % scal_var_dissip6 (i,j,k) , i = sx , ex ) , &
                                                              j = sy , ey ) , &
                                                              k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_ux (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_vy (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_wz (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_p (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_Y (i,j,k) , i = sx , ex ) , &
                                                        j = sy , ey ) , &
                                                        k = sz , ez )

       read (unit_stat) ((( stat % tmp_corr_Zm (i,j,k) , i = sx , ex ) , &
                                                         j = sy , ey ) , &
                                                         k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_11 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_22 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_33 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_12 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_13 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_dissip_23 (i,j,k) , i = sx , ex ) , &
                                                           j = sy , ey ) , &
                                                           k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_11 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_22 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_33 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_12 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_13 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) ((( stat % rey_press_strain_23 (i,j,k) , i = sx , ex ) , &
                                                                 j = sy , ey ) , &
                                                                 k = sz , ez )

       read (unit_stat) (((( stat % rey_transp_11 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_22 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_33 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_12 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_13 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) (((( stat % rey_transp_23 (i,j,k,l) , i = sx , ex      ) , &
                                                              j = sy , ey      ) , &
                                                              k = sz , ez      ) , &
                                                              l = 1  , ndimmax )

       read (unit_stat) ((( stat % pdf_Y (i,j,l) , i = sx , ex   ) , &
                                                   j = sy , ey   ) , &
                                                   l = 1 , nbins )

       read (unit_stat) ((( stat % pdf_Zm (i,j,l) , i = sx , ex   ) , &
                                                    j = sy , ey   ) , &
                                                    l = 1 , nbins )

       read (unit_stat) ((( stat % pdf_c (i,j,l) , i = sx , ex   ) , &
                                                   j = sy , ey   ) , &
                                                   l = 1 , nbins )

       read (unit_stat) ((( stat % pdf_sdr (i,j,l) , i = sx , ex   ) , &
                                                     j = sy , ey   ) , &
                                                     l = 1 , nbins )

       read (unit_stat) ((( stat % pdf_Da (i,j,l) , i = sx , ex   ) , &
                                                    j = sy , ey   ) , &
                                                    l = 1 , nbins )

       read (unit_stat) ((( stat % pdf_FI (i,j,l) , i = sx , ex   ) , &
                                                    j = sy , ey   ) , &
                                                    l = 1 , nbins )

       read (unit_stat) (( stat % cpdf_omega_FI (i,l) , i = sx , ex   ) , &
                                                        l = 1 , nbins )

       read (unit_stat) (( stat % cpdf_omega_Da (i,l) , i = sx , ex   ) , &
                                                        l = 1 , nbins )

       read (unit_stat) (( stat % spec_ux (k,l) , k = sz , ez    ) , &
                                                  l = 1 , ntimes )

       read (unit_stat) (( stat % spec_vy (k,l) , k = sz , ez    ) , &
                                                  l = 1 , ntimes )

       read (unit_stat) (( stat % spec_wz (k,l) , k = sz , ez    ) , &
                                                  l = 1 , ntimes )

       read (unit_stat) (( stat % spec_Y (k,l) , k = sz , ez    ) , &
                                                 l = 1 , ntimes )

       read (unit_stat) (( stat % spec_Zm (k,l) , k = sz , ez    ) , &
                                                  l = 1 , ntimes )

       read (unit_stat) (( stat % spec_sdr (k,l) , k = sz , ez    ) , &
                                                   l = 1 , ntimes )


    end if


    close (unit_stat)


  end subroutine read_stat_type


end module variables
