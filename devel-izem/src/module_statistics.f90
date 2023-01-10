!------------------------------------------------------------------------------
! MODULE: statistics
!------------------------------------------------------------------------------
!> \brief Calculation of statistical quantities.
!!
!! This module calculates all the statistical quantities such as
!! averages, fluctuations, variances...
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module statistics


  use parameters
  use input
  use type_thd
  use variables
  use thermodynamics
  use tools_post
  use budget_tke
  use budget_Reynolds
  use budget_vorticity
  use budget_scal_var
  use autocorrelations


  implicit none


  real (dp) , parameter , private :: epsi = 1.0e-10_dp !< to avoid divisions by zero


contains


!> \brief Calculate averages.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine averages ( ifile , inp , thd , adi , sim , dx_i , dy_i , dz_i , dt , inst , reyavg , favavg )


    integer (ip) , intent (in)                                     :: ifile  !< file number
    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (sim_type) , intent (in)                                  :: sim    !< simulation derived type
    real (dp) , dimension (:) , allocatable , intent (in)          :: dx_i   !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dy_i   !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dz_i   !< inverted dz array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dt     !< time step
    type (inst_type) , intent (in)                                 :: inst   !< instantaneous derived type
    type (reyavg_type) , intent (inout)                            :: reyavg !< Reynolds average derived type
    type (favavg_type) , intent (inout)                            :: favavg !< Favre average derived type


    integer (ip)                                    :: ok , end_file , spc
    real (dp)                                       :: dtime , sumdtime
    real (dp)                                       :: min_Da_adim  , max_Da_adim , &
                                                       min_sdr_adim , max_sdr_adim
    ! temporal variables
    real (dp) , dimension (:) , allocatable         :: dy , dz
    real (dp) , dimension (:,:,:) , allocatable     :: rho , rho_i
    real (dp) , dimension (:,:,:) , allocatable     :: wrk1 , wrk2 , wrk3 , wrk4
    real (dp) , dimension (:,:) , allocatable       :: wrk_dim2
    real (dp) , dimension (:,:,:,:,:) , allocatable :: wrk_dim5


    min_Da_adim = min_Da / ( ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) ) )
    max_Da_adim = max_Da / ( ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) ) )
    if ( inp % dis_rate_therm ) then
       min_sdr_adim = min_sdr / ( ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / ( adi % L_ref * adi % L_ref ) )
       max_sdr_adim = max_sdr / ( ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / ( adi % L_ref * adi % L_ref ) )
    else
       min_sdr_adim = min_sdr / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       max_sdr_adim = max_sdr / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
    end if


    allocate  ( dy       ( sy:ey )                                                       , &
                dz       ( sz:ez )                                                       , &
                rho      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                rho_i    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk1     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk2     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk3     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk4     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk_dim2 ( sx-ng:ex+ng , nbins )                                         , &
                wrk_dim5 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax ) , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate averages')


    end_file = sim % correct_endfile !inp % end_file
    dtime    = dt (ifile)
    sumdtime = sim % sumdtime

    rho (sx:ex,sy:ey,sz:ez)   = inst % u (sx:ex,sy:ey,sz:ez,1)
    rho_i (sx:ex,sy:ey,sz:ez) = 1.0_dp / inst % u (sx:ex,sy:ey,sz:ez,1)
    dy (sy:ey)                = 1.0_dp / dy_i (sy:ey)
    dz (sz:ez)                = 1.0_dp / dz_i (sz:ez)


    ! basic variables


    ! rho (this has to be placed first because of rho division in the Favre average
    wrk1 (sx:ex,sy:ey,sz:ez) = rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , wrk2 , wrk1 , reyavg % rho ) ! here wrk2 is a dummy variable

    ! Y
    do spc = 1 , nrv+npv+nvv
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) * rho_i (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc)
       wrk3 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       wrk4 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call spatZ_avg ( inp % spat_avg , dz , wrk1 )
       call spatZ_avg ( inp % spat_avg , dz , wrk2 )
       call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk3 )
       call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , wrk4 )
       reyavg % Y (sx:ex,sy:ey,sz:ez,spc) = wrk3 (sx:ex,sy:ey,sz:ez)
       favavg % Y (sx:ex,sy:ey,sz:ez,spc) = wrk4 (sx:ex,sy:ey,sz:ez)
    end do

    if ( LES ) then

!> Aditionnal variables (A.Techer)
    ! xiv2 (variance SGS - method 2)
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % xiv2 (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
       call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % xiv2 )
    ! xiv3 (variance SGS - method 3)
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % xiv3 (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
       call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % xiv3 )
    ! xixi : conservation of the square of the filtred mixture fraction
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % xixi (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
       call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % xixi )

    ! zjm and zjp: jump position for ignition
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % zjm (sx:ex,sy:ey,sz:ez)
       call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , wrk2 , wrk1 , reyavg % zjm )
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % zjp (sx:ex,sy:ey,sz:ez)
       call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , wrk2 , wrk1 , reyavg % zjp )
    ! integration pdf between zjm and zjp: ignition probability
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % intzjmp (sx:ex,sy:ey,sz:ez)
       call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , wrk2 , wrk1 , reyavg % intzjmp )
    ! taudelta: SGS dissipation time scale
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % taumix (sx:ex,sy:ey,sz:ez)
       call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , wrk2 , wrk1 , reyavg % taumix )

    end if

    ! Zm
    call Zmsimp ( inp , thd , inst % u , wrk1 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % Zm )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % Zm )

    ! ux
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % ux )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % ux )




    ! vy
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % vy )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % vy )

    ! wz
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % wz )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % wz )

    ! et
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,5) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,5)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % et )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % et )

    ! H
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % H (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % H )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % H )

    ! T
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % T (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % T )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % T )

    ! Tt
    call Ttmix ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % Tt )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % Tt )

    ! Ht
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * inst % u (sx:ex,sy:ey,sz:ez,2) + &
                               inst % u (sx:ex,sy:ey,sz:ez,3) * inst % u (sx:ex,sy:ey,sz:ez,3) + &
                               inst % u (sx:ex,sy:ey,sz:ez,4) * inst % u (sx:ex,sy:ey,sz:ez,4)
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % H (sx:ex,sy:ey,sz:ez) + 0.5_dp * wrk1 (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % Ht )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % Ht )

    ! P
    wrk1 (sx:ex,sy:ey,sz:ez) = rho (sx:ex,sy:ey,sz:ez) * inst % T (sx:ex,sy:ey,sz:ez) * inst % W_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % P )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % P )

    ! W
    wrk1 (sx:ex,sy:ey,sz:ez) = 1.0_dp / inst % W_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % W )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % W )

    ! cs
    call soundmix ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % cs )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % cs )

    ! cp
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % cp (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % cp )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % cp )


    ! viscous variables


    ! mu
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % mu (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % mu )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % mu )

    if ( LES ) then ! mu_SGS
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % mu_SGS (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
       call spatZ_avg ( inp % spat_avg , dz , wrk1 )
       call spatZ_avg ( inp % spat_avg , dz , wrk2 )
       call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % mu_SGS )
       call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % mu_SGS )
    end if

    ! nu
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % mu (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % mu (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % nu )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % nu )

    ! kpa
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % kpa (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % kpa )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % kpa )

    ! ct
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % ct (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % ct )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % ct )


    ! non dimensional variables


    ! gamma
    call gammamix ( 0 , thd , inst % W_i , inst % cp , wrk1 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % gamma )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % gamma )

    ! M
    call machmix ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % M )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % M )

    ! Pr
    if (reaction) then
       call prandtlmix ( 0 , inst % cp , inst % mu , inst % ct , wrk1 )
       wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
       call spatZ_avg ( inp % spat_avg , dz , wrk1 )
       call spatZ_avg ( inp % spat_avg , dz , wrk2 )
       call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % Pr )
       call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % Pr )
    end if


    ! reactive variables


    ! Ydot
    do spc = 1 , nrv
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       wrk4 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call spatZ_avg ( inp % spat_avg , dz , wrk1 )
       call spatZ_avg ( inp % spat_avg , dz , wrk2 )
       call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk3 )
       call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , wrk4 )
       reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc) = wrk3 (sx:ex,sy:ey,sz:ez)
       favavg % Ydot (sx:ex,sy:ey,sz:ez,spc) = wrk4 (sx:ex,sy:ey,sz:ez)
    end do

    ! omega
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % omega (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % omega )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % omega )

    ! Damkohler
    call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                     inst % dm , inst % tdr , inst % rd , inst % Ydot , wrk2 )
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk2 (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % Da )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % Da )

    ! sdr
    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call D_equivalent ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                        inst % ct , inst % cp , inst % dm , inst % tdr , inst % rd , wrk1 , wrk2 )
    call dissip_rate ( dx_i , dy_i , dz_i , wrk2 , wrk1 , wrk3 )
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk3 (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk3 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % sdr )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % sdr )

    ! FI
    call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk1 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , reyavg % FI )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , favavg % FI )


    ! custom variables


    ! Reynolds shear stresses


    ! tau
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk3 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) * rho_i (sx:ex,sy:ey,sz:ez)
    call tau_ij ( dx_i , dy_i , dz_i , wrk1 , wrk2 , wrk3 , inst % mu , inst % kpa , wrk_dim5 )

    ! tau_11
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim5 (sx:ex,sy:ey,sz:ez,1,1)
    wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,1,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    reyavg % tau (sx:ex,sy:ey,sz:ez,1,1) = wrk2 (sx:ex,sy:ey,sz:ez)

    ! tau_12
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim5 (sx:ex,sy:ey,sz:ez,1,2)
    wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,1,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    reyavg % tau (sx:ex,sy:ey,sz:ez,1,2) = wrk2 (sx:ex,sy:ey,sz:ez)

    ! tau_13
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim5 (sx:ex,sy:ey,sz:ez,1,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,1,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    reyavg % tau (sx:ex,sy:ey,sz:ez,1,3) = wrk2 (sx:ex,sy:ey,sz:ez)

    ! tau_21
    reyavg % tau (sx:ex,sy:ey,sz:ez,2,1) = reyavg % tau (sx:ex,sy:ey,sz:ez,1,2)

    ! tau_22
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim5 (sx:ex,sy:ey,sz:ez,2,2)
    wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,2,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    reyavg % tau (sx:ex,sy:ey,sz:ez,2,2) = wrk2 (sx:ex,sy:ey,sz:ez)

    ! tau_23
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim5 (sx:ex,sy:ey,sz:ez,2,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,2,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    reyavg % tau (sx:ex,sy:ey,sz:ez,2,3) = wrk2 (sx:ex,sy:ey,sz:ez)

    ! tau_31
    reyavg % tau (sx:ex,sy:ey,sz:ez,3,1) = reyavg % tau (sx:ex,sy:ey,sz:ez,1,3)

    ! tau_32
    reyavg % tau (sx:ex,sy:ey,sz:ez,3,2) = reyavg % tau (sx:ex,sy:ey,sz:ez,2,3)

    ! tau_33
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim5 (sx:ex,sy:ey,sz:ez,3,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,3,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    reyavg % tau (sx:ex,sy:ey,sz:ez,3,3) = wrk2 (sx:ex,sy:ey,sz:ez)


    ! conditional averages


    if ( inp % cond_avg ) then

    ! Y Reynolds conditional average
    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call cond_avg ( rey_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , rho , wrk_dim2 , reyavg % ca_rho_Zm ) ! wrk_dim2 is a dummy variable

    ! T_Zm Favre conditional average
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % T (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk2 , reyavg % ca_rho_Zm , favavg % ca_T_Zm )

    ! c_Zm Favre conditional average


    ! sdr_Zm Favre conditional average
    call D_equivalent ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                        inst % ct , inst % cp , inst % dm , inst % tdr , inst % rd , wrk1 , wrk2 )
    call dissip_rate ( dx_i , dy_i , dz_i , wrk2 , wrk1 , wrk3 )
    wrk3 (sx:ex,sy:ey,sz:ez) = wrk3 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk3 , reyavg % ca_rho_Zm , favavg % ca_sdr_Zm )

    ! Da_Zm Favre conditional average
    call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                     inst % dm , inst % tdr , inst % rd , inst % Ydot , wrk2 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk2 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk2 , reyavg % ca_rho_Zm , favavg % ca_Da_Zm )

    ! omega_Zm Favre conditional average
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % omega (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk2 , reyavg % ca_rho_Zm , favavg % ca_omega_Zm )

    ! FI_Zm Favre conditional average
    call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk2 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk2 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk2 , reyavg % ca_rho_Zm , favavg % ca_FI_Zm )


    ! FI Reynolds conditional average
    call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk1 )
    call cond_avg ( rey_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_FI , max_FI ,                                         &
                    wrk1 , rho , wrk_dim2 , reyavg % ca_rho_FI ) ! wrk_dim2 is a dummy variable

    ! omega_FI Favre conditional average
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % omega (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_FI , max_FI ,                                         &
                    wrk1 , wrk2 , reyavg % ca_rho_FI , favavg % ca_omega_FI )


    ! Da Reynolds conditional average
    call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                     inst % dm , inst % tdr , inst % rd , inst % Ydot , wrk1 )
    wrk1 (sx:ex,sy:ey,sz:ez) = log10 ( wrk1 (sx:ex,sy:ey,sz:ez) )
    call cond_avg ( rey_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    log10 (min_Da_adim) , log10 (max_Da_adim) ,               &
                    wrk1 , rho , wrk_dim2 , reyavg % ca_rho_Da ) ! wrk_dim2 is a dummy variable

    ! omega_Da Favre conditional average
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % omega (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    log10 (min_Da_adim) , log10 (max_Da_adim) ,               &
                    wrk1 , wrk2 , reyavg % ca_rho_Da , favavg % ca_omega_Da )

    end if


    deallocate ( dy , dz , rho , rho_i )
    deallocate ( wrk1 , wrk2 , wrk3 , wrk4 )
    deallocate ( wrk_dim2 , wrk_dim5 )


  end subroutine averages


!> \brief Calculate probability density functions.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine pdfs ( ifile , inp , thd , adi , sim , dx_i , dy_i , dz_i , dt , inst , favavg , stat )


    integer (ip) , intent (in)                                     :: ifile  !< file number
    type (inp_type) , intent (in)                                  :: inp    !< input derived type
    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi    !< non-dimensional derived type
    type (sim_type) , intent (in)                                  :: sim    !< simulation derived type
    real (dp) , dimension (:) , allocatable , intent (in)          :: dx_i   !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dy_i   !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dz_i   !< inverted dz array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dt     !< time step
    type (inst_type) , intent (in)                                 :: inst   !< instantaneous derived type
    type (favavg_type) , intent (in)                               :: favavg !< Favre averaged derived type
    type (stat_type) , intent (inout)                              :: stat   !< statistical derived type


    integer (ip)                                    :: ok , end_file
    real (dp)                                       :: dtime , sumdtime
    real (dp)                                       :: min_Da_adim  , max_Da_adim , &
                                                       min_sdr_adim , max_sdr_adim
    ! temporal variables
    real (dp) , dimension (:,:,:) , allocatable     :: rho , rho_i
    real (dp) , dimension (:,:,:) , allocatable     :: wrk1 , wrk2 , wrk3


    min_Da_adim = min_Da / ( ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) ) )
    max_Da_adim = max_Da / ( ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) ) )
    if ( inp % dis_rate_therm ) then
       min_sdr_adim = min_sdr / ( ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / ( adi % L_ref * adi % L_ref ) )
       max_sdr_adim = max_sdr / ( ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / ( adi % L_ref * adi % L_ref ) )
    else
       min_sdr_adim = min_sdr / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       max_sdr_adim = max_sdr / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
    end if


    allocate  ( rho   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                rho_i ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                wrk1  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                wrk2  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                wrk3  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate pdfs')


    end_file = sim % correct_endfile !inp % end_file
    dtime    = dt (ifile)
    sumdtime = sim % sumdtime

    rho (sx:ex,sy:ey,sz:ez)   = inst % u (sx:ex,sy:ey,sz:ez,1)
    rho_i (sx:ex,sy:ey,sz:ez) = 1.0_dp / inst % u (sx:ex,sy:ey,sz:ez,1)


    ! typical pdfs


    ! pdf_Y
    if ( inp % index_scalar == 0 ) then
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv + inp % index_scalar ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call pdf ( decimal , ifile , end_file , min_Y , max_Y , wrk1 , stat % pdf_Y )

    ! pdf_Zm
    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call pdf ( decimal , ifile , end_file , min_Y , max_Y , wrk1 , stat % pdf_Zm )

    ! pdf_c


    ! pdf_sdr
    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call D_equivalent ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                        inst % ct , inst % cp , inst % dm , inst % tdr , inst % rd , wrk1 , wrk2 )
    call dissip_rate ( dx_i , dy_i , dz_i , wrk2 , wrk1 , wrk3 )
    wrk3 (sx:ex,sy:ey,sz:ez) = log10 ( wrk3 (sx:ex,sy:ey,sz:ez) )
    call pdf ( logarit , ifile , end_file , log10 (min_sdr_adim) , log10 (max_sdr_adim) , wrk3 , stat % pdf_sdr )

    ! pdf_Da
    call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                     inst % dm , inst % tdr , inst % rd , inst % Ydot , wrk1 )
    wrk1 (sx:ex,sy:ey,sz:ez) = log10 ( wrk1 (sx:ex,sy:ey,sz:ez) )
    call pdf ( logarit , ifile , end_file , log10 (min_Da_adim) , log10 (max_Da_adim) , wrk1 , stat % pdf_Da )

    ! pdf_FI
    call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk1 )
    call pdf ( decimal , ifile , end_file , min_FI , max_FI , wrk1 , stat % pdf_FI )


    ! conditional pdfs


    if (reaction) then ! this only has sense in a reactive problem

    ! axpdf_FI
    call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk1 )
    call pdf_conditional ( decimal , ifile , end_file , min_FI , max_FI , favavg % ca_omega_FI , &
                           wrk1 , stat % cpdf_omega_FI )

    ! axpdf_Da
    call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                     inst % dm , inst % tdr , inst % rd , inst % Ydot , wrk1 )
    wrk1 (sx:ex,sy:ey,sz:ez) = log10 ( wrk1 (sx:ex,sy:ey,sz:ez) )
    call pdf_conditional ( logarit , ifile , end_file , log10 (min_Da_adim) , log10 (max_Da_adim) , &
                           favavg % ca_omega_Da , wrk1 , stat % cpdf_omega_Da )

    end if


    deallocate ( rho , rho_i )
    deallocate ( wrk1 , wrk2 , wrk3 )


  end subroutine pdfs


!> \brief Calculate spectras.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine spectras ( ifile , inp , thd , sim , x , y , z , dx_i , dy_i , dz_i , inst , stat )


    integer (ip) , intent (in)                                     :: ifile !< file number
    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type
    type (sim_type) , intent (in)                                  :: sim   !< simulation derived type
    real (dp) , dimension (:) , allocatable , intent (in)          :: x     !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)          :: y     !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)          :: z     !< z-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dz_i  !< inverted dz array
    type (inst_type) , intent (in)                                 :: inst  !< instantaneous derived type
    type (stat_type) , intent (inout)                              :: stat  !< statistical derived type


    integer (ip)                                    :: ok , i , j , k , n
    integer (ip)                                    :: start_file , end_file
    integer (ip)                                    :: px , py , pz
    real (dp) , dimension (ndimmax)                 :: pt
    real (dp) , dimension (:) , allocatable         :: dy , dz
    real (dp) , dimension (:,:,:) , allocatable     :: rho_i
    real (dp) , dimension (:) , allocatable         :: wrk_dim1
    real (dp) , dimension (:,:,:) , allocatable     :: wrk1 , wrk2 , wrk3


    if ( inp % skip_file > 1 ) call abort_mpi ('skip_file must be one for DFT')


    allocate  ( dy       ( sy:ey )                                   , &
                dz       ( sz:ez )                                   , &
                rho_i    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                wrk_dim1 ( 1:ntimes )                                , &
                wrk1     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                wrk2     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                wrk3     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate spectras')


    start_file = inp % start_file
    end_file   = sim % correct_endfile !inp % end_file


    rho_i (sx:ex,sy:ey,sz:ez) = 1.0_dp / inst % u (sx:ex,sy:ey,sz:ez,1)
    dy (sy:ey)                = 1.0_dp / dy_i (sy:ey)
    dz (sz:ez)                = 1.0_dp / dz_i (sz:ez)


    ! spectra


    px = ex ; py = ey ; pz = ez ! initialize this to avoid bugs (verified)

    pt (1) = inp % corrspec_coord (1)
    pt (2) = inp % corrspec_coord (2)
    pt (3) = inp % corrspec_coord (3)

    ! calculate points
    do i = ex , sx , -1
       if ( x (i) > pt (1) ) px = i
    end do
    do j = ey , sy , -1
       if ( y (j) > pt (2) ) py = j
    end do
    if ( ndim == 3 ) then
       do k = ez , sz , -1
          if ( z (k) > pt (3) ) pz = k
       end do
    else
       pz = sz
    end if


    ! spec_ux
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * rho_i (sx:ex,sy:ey,sz:ez)
!    call spatYZ_avg ( inp % spat_avg , dy , dz , wrk1 ) ! new thing
    call similarity_x ( inp , wrk1 )
    stat % spec_ux (sz:ez,ifile) = wrk1 (px,py,sz:ez)
    if ( ifile == end_file ) then
!       if ( rank == rank_default ) write (*,*) 'entering spec_ux ...'
       do k = sz,ez
          wrk_dim1 (:) = stat % spec_ux (k,:)
          call DFT ( start_file , end_file , wrk_dim1 )
          stat % spec_ux (k,:) = wrk_dim1 (:)
       end do
!       if ( rank == rank_default ) write (*,*) 'averaging spec_ux ...'
       do n = 1,ntimes
          wrk1 (px,py,sz:ez) = stat % spec_ux (sz:ez,n)
          call spatZ_avg ( inp % spat_avg , dz , wrk1 )
          stat % spec_ux (sz:ez,n) = wrk1 (px,py,sz:ez)
       end do
    end if


    ! spec_vy
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * rho_i (sx:ex,sy:ey,sz:ez)
!    call spatYZ_avg ( inp % spat_avg , dy , dz , wrk1 ) ! new thing
    call similarity_x ( inp , wrk1 )
    stat % spec_vy (sz:ez,ifile) = wrk1 (px,py,sz:ez)
    if ( ifile == end_file ) then
!       if ( rank == rank_default ) write (*,*) 'entering spec_vy ...'
       do k = sz,ez
          wrk_dim1 (:) = stat % spec_vy (k,:)
          call DFT ( start_file , end_file , wrk_dim1 )
          stat % spec_vy (k,:) = wrk_dim1 (:)
       end do
!       if ( rank == rank_default ) write (*,*) 'averaging spec_vy ...'
       do n = 1,ntimes
          wrk1 (px,py,sz:ez) = stat % spec_vy (sz:ez,n)
          call spatZ_avg ( inp % spat_avg , dz , wrk1 )
          stat % spec_vy (sz:ez,n) = wrk1 (px,py,sz:ez)
       end do
    end if


    ! spec_wz
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) * rho_i (sx:ex,sy:ey,sz:ez)
!    call spatYZ_avg ( inp % spat_avg , dy , dz , wrk1 ) ! new thing
    call similarity_x ( inp , wrk1 )
    stat % spec_wz (sz:ez,ifile) = wrk1 (px,py,sz:ez)
    if ( ifile == end_file ) then
!       if ( rank == rank_default ) write (*,*) 'entering spec_wz ...'
       do k = sz,ez
          wrk_dim1 (:) = stat % spec_wz (k,:)
          call DFT ( start_file , end_file , wrk_dim1 )
          stat % spec_wz (k,:) = wrk_dim1 (:)
       end do
!       if ( rank == rank_default ) write (*,*) 'averaging spec_wz ...'
       do n = 1,ntimes
          wrk1 (px,py,sz:ez) = stat % spec_wz (sz:ez,n)
          call spatZ_avg ( inp % spat_avg , dz , wrk1 )
          stat % spec_wz (sz:ez,n) = wrk1 (px,py,sz:ez)
       end do
    end if


    ! spec_Y
    if ( inp % index_scalar == 0 ) then
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv + inp % index_scalar ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
!    call spatYZ_avg ( inp % spat_avg , dy , dz , wrk1 ) ! new thing
    call similarity_x ( inp , wrk1 )
    stat % spec_Y (sz:ez,ifile) = wrk1 (px,py,sz:ez)
    if ( ifile == end_file ) then
!       if ( rank == rank_default ) write (*,*) 'entering spec_Y ...'
       do k = sz,ez
          wrk_dim1 (:) = stat % spec_Y (k,:)
          call DFT ( start_file , end_file , wrk_dim1 )
          stat % spec_Y (k,:) = wrk_dim1 (:)
       end do
!       if ( rank == rank_default ) write (*,*) 'averaging spec_Y ...'
       do n = 1,ntimes
          wrk1 (px,py,sz:ez) = stat % spec_Y (sz:ez,n)
          call spatZ_avg ( inp % spat_avg , dz , wrk1 )
          stat % spec_Y (sz:ez,n) = wrk1 (px,py,sz:ez)
       end do
    end if


    ! spec_sdr
    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call D_equivalent ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                        inst % ct , inst % cp , inst % dm , inst % tdr , inst % rd , wrk1 , wrk2 )
    call dissip_rate ( dx_i , dy_i , dz_i , wrk2 , wrk1 , wrk3 )
    call similarity_x ( inp , wrk3 )
    stat % spec_sdr (sz:ez,ifile) = wrk3 (px,py,sz:ez)
    if ( ifile == end_file ) then
!       if ( rank == rank_default ) write (*,*) 'entering spec_sdr ...'
       do k = sz,ez
          wrk_dim1 (:) = stat % spec_sdr (k,:)
          call DFT ( start_file , end_file , wrk_dim1 )
          stat % spec_sdr (k,:) = wrk_dim1 (:)
       end do
!       if ( rank == rank_default ) write (*,*) 'averaging spec_sdr ...'
       do n = 1,ntimes
          wrk1 (px,py,sz:ez) = stat % spec_sdr (sz:ez,n)
          call spatZ_avg ( inp % spat_avg , dz , wrk1 )
          stat % spec_sdr (sz:ez,n) = wrk1 (px,py,sz:ez)
       end do
    end if


    deallocate ( dy , dz )
    deallocate ( rho_i )
    deallocate ( wrk_dim1 )
    deallocate ( wrk1 , wrk2 , wrk3 )


  end subroutine spectras


!> \brief Calculate fluctuations.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine fluctuations ( inp , thd , adi , dx_i , dy_i , dz_i , inst , reyavg , favavg , reyfluc , favfluc )


    type (inp_type) , intent (in)                                  :: inp     !< input derived type
    type (thd_type) , intent (in)                                  :: thd     !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi     !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)          :: dx_i    !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dy_i    !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dz_i    !< inverted dz array
    type (inst_type) , intent (in)                                 :: inst    !< instantaneous derived type
    type (reyavg_type) , intent (in)                               :: reyavg  !< Reynolds average derived type
    type (favavg_type) , intent (in)                               :: favavg  !< Favre average derived type
    type (reyfluc_type) , intent (inout)                           :: reyfluc !< Reynolds fluctuation derived type
    type (favfluc_type) , intent (inout)                           :: favfluc !< Reynolds fluctuation derived type


    integer (ip)                                    :: ok , spc
    real (dp)                                       :: min_Da_adim  , max_Da_adim , &
                                                       min_sdr_adim , max_sdr_adim
    ! temporal variables
    real (dp) , dimension (:,:,:) , allocatable     :: rho_i
    real (dp) , dimension (:,:,:) , allocatable     :: wrk1 , wrk2 , wrk3 , wrk4
    real (dp) , dimension (:,:,:,:,:) , allocatable :: wrk_dim5


    min_Da_adim = min_Da / ( ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) ) )
    max_Da_adim = max_Da / ( ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) ) )
    if ( inp % dis_rate_therm ) then
       min_sdr_adim = min_sdr / ( ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / ( adi % L_ref * adi % L_ref ) )
       max_sdr_adim = max_sdr / ( ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / ( adi % L_ref * adi % L_ref ) )
    else
       min_sdr_adim = min_sdr / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       max_sdr_adim = max_sdr / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
    end if


    allocate  ( rho_i    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk1     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk2     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk3     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk4     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk_dim5 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax ) , &
                stat  = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate fluctuations')


    rho_i (sx:ex,sy:ey,sz:ez) = 1.0_dp / inst % u (sx:ex,sy:ey,sz:ez,1)


    ! basic variables


    ! rho
    reyfluc % rho (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,1) - &
                                        reyavg % rho (sx:ex,sy:ey,sz:ez)

    ! Y
    do spc = 1 , nrv+npv+nvv
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) * rho_i (sx:ex,sy:ey,sz:ez)
       reyfluc % Y (sx:ex,sy:ey,sz:ez,spc) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       favfluc % Y (sx:ex,sy:ey,sz:ez,spc) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % Y (sx:ex,sy:ey,sz:ez,spc)
    end do

    ! Zm
    call Zmsimp ( inp , thd , inst % u , wrk1 )
    reyfluc % Zm (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % Zm (sx:ex,sy:ey,sz:ez)
    favfluc % Zm (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % Zm (sx:ex,sy:ey,sz:ez)

    ! ux
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * rho_i (sx:ex,sy:ey,sz:ez)
    reyfluc % ux (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % ux (sx:ex,sy:ey,sz:ez)
    favfluc % ux (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % ux (sx:ex,sy:ey,sz:ez)

    ! vy
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * rho_i (sx:ex,sy:ey,sz:ez)
    reyfluc % vy (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % vy (sx:ex,sy:ey,sz:ez)
    favfluc % vy (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % vy (sx:ex,sy:ey,sz:ez)

    ! wz
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) * rho_i (sx:ex,sy:ey,sz:ez)
    reyfluc % wz (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % wz (sx:ex,sy:ey,sz:ez)
    favfluc % wz (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % wz (sx:ex,sy:ey,sz:ez)

    ! et
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,5) * rho_i (sx:ex,sy:ey,sz:ez)
    reyfluc % et (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % et (sx:ex,sy:ey,sz:ez)
    favfluc % et (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % et (sx:ex,sy:ey,sz:ez)

    ! H
    reyfluc % H (sx:ex,sy:ey,sz:ez) = inst % H (sx:ex,sy:ey,sz:ez) - reyavg % H (sx:ex,sy:ey,sz:ez)
    favfluc % H (sx:ex,sy:ey,sz:ez) = inst % H (sx:ex,sy:ey,sz:ez) - favavg % H (sx:ex,sy:ey,sz:ez)

    ! T
    reyfluc % T (sx:ex,sy:ey,sz:ez) = reyavg % T (sx:ex,sy:ey,sz:ez) - inst % T (sx:ex,sy:ey,sz:ez)
    favfluc % T (sx:ex,sy:ey,sz:ez) = favavg % T (sx:ex,sy:ey,sz:ez) - inst % T (sx:ex,sy:ey,sz:ez)

    ! Tt
    call Ttmix ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 )
    reyfluc % Tt (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % Tt (sx:ex,sy:ey,sz:ez)
    favfluc % Tt (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % Tt (sx:ex,sy:ey,sz:ez)

    ! P
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,1) * inst % T (sx:ex,sy:ey,sz:ez) * inst % W_i (sx:ex,sy:ey,sz:ez)
    reyfluc % P (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % P (sx:ex,sy:ey,sz:ez)
    favfluc % P (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % P (sx:ex,sy:ey,sz:ez)

    ! W
    reyfluc % W (sx:ex,sy:ey,sz:ez) = 1.0_dp / inst % W_i (sx:ex,sy:ey,sz:ez) - reyavg % W (sx:ex,sy:ey,sz:ez)
    favfluc % W (sx:ex,sy:ey,sz:ez) = 1.0_dp / inst % W_i (sx:ex,sy:ey,sz:ez) - favavg % W (sx:ex,sy:ey,sz:ez)

    ! cs
    call soundmix ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 )
    reyfluc % cs (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % cs (sx:ex,sy:ey,sz:ez)
    favfluc % cs (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % cs (sx:ex,sy:ey,sz:ez)

    ! cp
    reyfluc % cp (sx:ex,sy:ey,sz:ez) = inst % cp (sx:ex,sy:ey,sz:ez) - reyavg % cp (sx:ex,sy:ey,sz:ez)
    favfluc % cp (sx:ex,sy:ey,sz:ez) = inst % cp (sx:ex,sy:ey,sz:ez) - favavg % cp (sx:ex,sy:ey,sz:ez)


    ! viscous variables


    ! mu
    reyfluc % mu (sx:ex,sy:ey,sz:ez) = inst % mu (sx:ex,sy:ey,sz:ez) - reyavg % mu (sx:ex,sy:ey,sz:ez)
    favfluc % mu (sx:ex,sy:ey,sz:ez) = inst % mu (sx:ex,sy:ey,sz:ez) - favavg % mu (sx:ex,sy:ey,sz:ez)

    ! nu
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % mu (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez)
    reyfluc % nu (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % nu (sx:ex,sy:ey,sz:ez)
    favfluc % nu (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % nu (sx:ex,sy:ey,sz:ez)

    ! kpa
    reyfluc % kpa (sx:ex,sy:ey,sz:ez) = inst % kpa (sx:ex,sy:ey,sz:ez) - reyavg % kpa (sx:ex,sy:ey,sz:ez)
    favfluc % kpa (sx:ex,sy:ey,sz:ez) = inst % kpa (sx:ex,sy:ey,sz:ez) - favavg % kpa (sx:ex,sy:ey,sz:ez)

    ! ct
    reyfluc % ct (sx:ex,sy:ey,sz:ez) = inst % ct (sx:ex,sy:ey,sz:ez) - reyavg % ct (sx:ex,sy:ey,sz:ez)
    favfluc % ct (sx:ex,sy:ey,sz:ez) = inst % ct (sx:ex,sy:ey,sz:ez) - favavg % ct (sx:ex,sy:ey,sz:ez)


    ! non dimensional variables


    ! gamma
    call gammamix ( 0 , thd , inst % W_i , inst % cp , wrk1 )
    reyfluc % gamma (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % gamma (sx:ex,sy:ey,sz:ez)
    favfluc % gamma (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % gamma (sx:ex,sy:ey,sz:ez)

    ! M
    call machmix ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 )
    reyfluc % M (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % M (sx:ex,sy:ey,sz:ez)
    favfluc % M (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % M (sx:ex,sy:ey,sz:ez)

    ! Pr
    if (reaction) then
       call prandtlmix ( 0 , inst % cp , inst % mu , inst % ct , wrk1 )
       reyfluc % Pr (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % Pr (sx:ex,sy:ey,sz:ez)
       favfluc % Pr (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % Pr (sx:ex,sy:ey,sz:ez)
    end if


    ! reactive variables


    ! Ydot
    do spc = 1 , nrv
       reyfluc % Ydot (sx:ex,sy:ey,sz:ez,spc) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc) - reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       favfluc % Ydot (sx:ex,sy:ey,sz:ez,spc) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc) - favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
    end do

    ! omega
    reyfluc % omega (sx:ex,sy:ey,sz:ez) = inst % omega (sx:ex,sy:ey,sz:ez) - reyavg % omega (sx:ex,sy:ey,sz:ez)
    favfluc % omega (sx:ex,sy:ey,sz:ez) = inst % omega (sx:ex,sy:ey,sz:ez) - favavg % omega (sx:ex,sy:ey,sz:ez)

    ! Da
    call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                     inst % dm , inst % tdr , inst % rd , inst % Ydot , wrk2 )
    reyfluc % Da (sx:ex,sy:ey,sz:ez) = wrk2 (sx:ex,sy:ey,sz:ez) - reyavg % Da (sx:ex,sy:ey,sz:ez)
    favfluc % Da (sx:ex,sy:ey,sz:ez) = wrk2 (sx:ex,sy:ey,sz:ez) - favavg % Da (sx:ex,sy:ey,sz:ez)

    ! sdr
    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call D_equivalent ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                        inst % ct , inst % cp , inst % dm , inst % tdr , inst % rd , wrk1 , wrk2 )
    call dissip_rate ( dx_i , dy_i , dz_i , wrk2 , wrk1 , wrk3 )
    reyfluc % sdr (sx:ex,sy:ey,sz:ez) = wrk3 (sx:ex,sy:ey,sz:ez) - reyavg % sdr (sx:ex,sy:ey,sz:ez)
    favfluc % sdr (sx:ex,sy:ey,sz:ez) = wrk3 (sx:ex,sy:ey,sz:ez) - favavg % sdr (sx:ex,sy:ey,sz:ez)

    ! FI
    call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk1 )
    reyfluc % FI (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - reyavg % FI (sx:ex,sy:ey,sz:ez)
    favfluc % FI (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - favavg % FI (sx:ex,sy:ey,sz:ez)


    ! custom variables


    ! Reynolds shear stresses


    ! tau
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk3 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) * rho_i (sx:ex,sy:ey,sz:ez)
    call tau_ij ( dx_i , dy_i , dz_i , wrk1 , wrk2 , wrk3 , inst % mu , inst % kpa , wrk_dim5 )

    ! tau_11
    reyfluc % tau (sx:ex,sy:ey,sz:ez,1,1) = wrk_dim5 (sx:ex,sy:ey,sz:ez,1,1) - reyavg % tau (sx:ex,sy:ey,sz:ez,1,1)

    ! tau_12
    reyfluc % tau (sx:ex,sy:ey,sz:ez,1,2) = wrk_dim5 (sx:ex,sy:ey,sz:ez,1,2) - reyavg % tau (sx:ex,sy:ey,sz:ez,1,2)

    ! tau_13
    reyfluc % tau (sx:ex,sy:ey,sz:ez,1,3) = wrk_dim5 (sx:ex,sy:ey,sz:ez,1,3) - reyavg % tau (sx:ex,sy:ey,sz:ez,1,3)

    ! tau_21
    reyfluc % tau (sx:ex,sy:ey,sz:ez,2,1) = reyfluc % tau (sx:ex,sy:ey,sz:ez,1,2)

    ! tau_22
    reyfluc % tau (sx:ex,sy:ey,sz:ez,2,2) = wrk_dim5 (sx:ex,sy:ey,sz:ez,2,2) - reyavg % tau (sx:ex,sy:ey,sz:ez,2,2)

    ! tau_23
    reyfluc % tau (sx:ex,sy:ey,sz:ez,2,3) = wrk_dim5 (sx:ex,sy:ey,sz:ez,2,3) - reyavg % tau (sx:ex,sy:ey,sz:ez,2,3)

    ! tau_31
    reyfluc % tau (sx:ex,sy:ey,sz:ez,3,1) = reyfluc % tau (sx:ex,sy:ey,sz:ez,1,3)

    ! tau_32
    reyfluc % tau (sx:ex,sy:ey,sz:ez,3,2) = reyfluc % tau (sx:ex,sy:ey,sz:ez,2,3)

    ! tau_33
    reyfluc % tau (sx:ex,sy:ey,sz:ez,3,3) = wrk_dim5 (sx:ex,sy:ey,sz:ez,3,3) - reyavg % tau (sx:ex,sy:ey,sz:ez,3,3)


    ! conditional fluctuations


    if ( inp % cond_avg ) then

    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if

    ! T_Zm
    call cond_fluc ( min_Y , max_Y , wrk1 , inst % T , favavg % ca_T_Zm , favfluc % ca_T_Zm )

    ! c_Zm


    ! sdr_Zm
    call D_equivalent ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                        inst % ct , inst % cp , inst % dm , inst % tdr , inst % rd , wrk1 , wrk2 )
    call dissip_rate ( dx_i , dy_i , dz_i , wrk2 , wrk1 , wrk3 )
    call cond_fluc ( min_Y , max_Y , wrk1 , wrk3 , favavg % ca_sdr_Zm , favfluc % ca_sdr_Zm )

    ! Da_Zm
    call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                     inst % dm , inst % tdr , inst % rd , inst % Ydot , wrk2 )
    call cond_fluc ( min_Y , max_Y , wrk1 , wrk2 , favavg % ca_Da_Zm , favfluc % ca_Da_Zm )

    ! omega_Zm
    call cond_fluc ( min_Y , max_Y , wrk1 , inst % omega , favavg % ca_omega_Zm , favfluc % ca_omega_Zm )

    ! FI_Zm
    call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk2 )
    call cond_fluc ( min_Y , max_Y , wrk1 , wrk2 , favavg % ca_FI_Zm , favfluc % ca_FI_Zm )


    ! omega_FI
    call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk1 )
    call cond_fluc ( min_FI , max_FI , wrk1 , inst % omega , favavg % ca_omega_FI , favfluc % ca_omega_FI )


    ! omega_Da
    call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                     inst % dm , inst % tdr , inst % rd , inst % Ydot , wrk1 )
    wrk1 (sx:ex,sy:ey,sz:ez) = log10 ( wrk1 (sx:ex,sy:ey,sz:ez) )
    call cond_fluc ( log10 (min_Da_adim) , log10 (max_Da_adim) , wrk1 , inst % omega , &
                     favavg % ca_omega_Da , favfluc % ca_omega_Da )

    end if


    deallocate (rho_i)
    deallocate ( wrk1 , wrk2 , wrk3 , wrk4 )
    deallocate (wrk_dim5)


  end subroutine fluctuations


!> \brief Calculate more statistical variables.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine stats ( ifile , inp , thd , adi , sim , xt , yt , zt , dx_i , dy_i , dz_i , dt , &
                     inst , reyavg , favavg , reyfluc , favfluc , stat )


    integer (ip) , intent (in)                                     :: ifile   !< file number
    type (inp_type) , intent (in)                                  :: inp     !< input derived type
    type (thd_type) , intent (in)                                  :: thd     !< thermodynamic derived type
    type (adi_type) , intent (in)                                  :: adi     !< non-dimensional derived type
    type (sim_type) , intent (in)                                  :: sim     !< simulation derived type
    real (dp) , dimension (:) , allocatable , intent (in)          :: xt      !< absolute x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)          :: yt      !< absolute y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)          :: zt      !< absolute z-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dx_i    !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dy_i    !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dz_i    !< inverted dz array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dt      !< time step
    type (inst_type) , intent (in)                                 :: inst    !< instantaneous derived type
    type (reyavg_type) , intent (inout)                            :: reyavg  !< Reynolds average derived type
    type (favavg_type) , intent (inout)                            :: favavg  !< Favre average derived type
    type (reyfluc_type) , intent (inout)                           :: reyfluc !< Reynolds fluctuation derived type
    type (favfluc_type) , intent (inout)                           :: favfluc !< Favre fluctuation derived type
    type (stat_type) , intent (inout)                              :: stat    !< statistical derived type


    integer (ip)                                    :: spc , ok , end_file
    real (dp)                                       :: dtime , sumdtime
    real (dp)                                       :: min_Da_adim  , max_Da_adim , &
                                                       min_sdr_adim , max_sdr_adim
    logical                                         :: loop
    character (len_default)                         :: name
    ! temporal variables
    real (dp) , dimension (:) , allocatable         :: dy , dz
    real (dp) , dimension (:,:,:) , allocatable     :: rho , rho_i
    real (dp) , dimension (:,:,:) , allocatable     :: wrk1 , wrk2 , wrk3 , wrk4 , wrk5 , wrk6
    real (dp) , dimension (:,:) , allocatable       :: wrk_dim2
    real (dp) , dimension (:,:,:,:) , allocatable   :: wrk_dim4
    real (dp) , dimension (:,:,:,:,:) , allocatable :: wrk_dim5


    min_Da_adim = min_Da / ( ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) ) )
    max_Da_adim = max_Da / ( ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) ) )
    if ( inp % dis_rate_therm ) then
       min_sdr_adim = min_sdr / ( ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / ( adi % L_ref * adi % L_ref ) )
       max_sdr_adim = max_sdr / ( ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / ( adi % L_ref * adi % L_ref ) )
    else
       min_sdr_adim = min_sdr / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       max_sdr_adim = max_sdr / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
    end if


    allocate  ( dy       ( sy:ey )                                                       , &
                dz       ( sz:ez )                                                       , &
                rho      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng                     ) , &
                rho_i    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng                     ) , &
                wrk1     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng                     ) , &
                wrk2     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng                     ) , &
                wrk3     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng                     ) , &
                wrk4     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng                     ) , &
                wrk5     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng                     ) , &
                wrk6     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng                     ) , &
                wrk_dim2 ( sz-ng:ez+ng , sz-ng:ez+ng                                   ) , &
                wrk_dim4 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax           ) , &
                wrk_dim5 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax ) , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate stats')


    end_file = sim % correct_endfile !inp % end_file
    dtime    = dt (ifile)
    sumdtime = sim % sumdtime

    rho (sx:ex,sy:ey,sz:ez)   = inst % u (sx:ex,sy:ey,sz:ez,1)
    rho_i (sx:ex,sy:ey,sz:ez) = 1.0_dp / inst % u (sx:ex,sy:ey,sz:ez,1)
    dy (sy:ey)                = 1.0_dp / dy_i (sy:ey)
    dz (sz:ez)                = 1.0_dp / dz_i (sz:ez)


    ! "repeated" Reynolds averages for convergence plots


    ! P
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,1) * &
                               inst % T (sx:ex,sy:ey,sz:ez) *   &
                               inst % W_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_P )

    ! Y
    if ( inp % index_scalar == 0 ) then
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv + inp % index_scalar ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_Y )

    ! Zm
    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_Zm )

    ! ux
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * rho_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_ux )

    ! vy
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * rho_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_vy )

    ! wz
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) * rho_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_wz )


    ! other averages for convergence plots


    ! uxux
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * inst % u (sx:ex,sy:ey,sz:ez,2) * &
                               rho_i (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_uxux )

    ! vyvy
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * inst % u (sx:ex,sy:ey,sz:ez,3) * &
                               rho_i (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_vyvy )

    ! wzwz
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) * inst % u (sx:ex,sy:ey,sz:ez,4) * &
                               rho_i (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_wzwz )

    ! uxvy
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * inst % u (sx:ex,sy:ey,sz:ez,3) * &
                               rho_i (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_uxvy )

    ! uxwz
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * inst % u (sx:ex,sy:ey,sz:ez,4) * &
                               rho_i (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_uxwz )

    ! vywz
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * inst % u (sx:ex,sy:ey,sz:ez,4) * &
                               rho_i (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyavg_vywz )

    ! enstrophy
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk3 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) * rho_i (sx:ex,sy:ey,sz:ez)
    call vorticity ( dx_i , dy_i , dz_i , wrk1 , wrk2 , wrk3 , wrk4 )
    wrk4 (sx:ex,sy:ey,sz:ez) = wrk4 (sx:ex,sy:ey,sz:ez) * wrk4 (sx:ex,sy:ey,sz:ez) ! vorticity is the root square of the enstrophy
    call spatZ_avg ( inp % spat_avg , dz , wrk4 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk4 , stat % reyavg_enstrophy )


    ! Reynolds and Favre variances


    ! YY
    do spc = 1 , nrv+npv+nvv
       wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % Y (sx:ex,sy:ey,sz:ez,spc) * reyfluc % Y (sx:ex,sy:ey,sz:ez,spc)
       wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % Y (sx:ex,sy:ey,sz:ez,spc) * favfluc % Y (sx:ex,sy:ey,sz:ez,spc) * &
                                  rho (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       wrk4 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call spatZ_avg ( inp % spat_avg , dz , wrk1 )
       call spatZ_avg ( inp % spat_avg , dz , wrk2 )
       call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk3 )
       call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , wrk4 )
       stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc) = wrk3 (sx:ex,sy:ey,sz:ez)
       stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc) = wrk4 (sx:ex,sy:ey,sz:ez)
    end do


    ! ZmZm
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % Zm (sx:ex,sy:ey,sz:ez) * reyfluc % Zm (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % Zm (sx:ex,sy:ey,sz:ez) * favfluc % Zm (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_Zm )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_Zm )

    ! uxux
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % ux (sx:ex,sy:ey,sz:ez) * reyfluc % ux (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ux (sx:ex,sy:ey,sz:ez) * favfluc % ux (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_ux )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_ux )

    ! vyvy
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % vy (sx:ex,sy:ey,sz:ez) * reyfluc % vy (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % vy (sx:ex,sy:ey,sz:ez) * favfluc % vy (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_vy )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_vy )

    ! wzwz
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % wz (sx:ex,sy:ey,sz:ez) * reyfluc % wz (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % wz (sx:ex,sy:ey,sz:ez) * favfluc % wz (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_wz )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_wz )

    ! uxvy
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % ux (sx:ex,sy:ey,sz:ez) * reyfluc % vy (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ux (sx:ex,sy:ey,sz:ez) * favfluc % vy (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_uxvy )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_uxvy )

    ! uxwz
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % ux (sx:ex,sy:ey,sz:ez) * reyfluc % wz (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ux (sx:ex,sy:ey,sz:ez) * favfluc % wz (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_uxwz )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_uxwz )

    ! vywz
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % vy (sx:ex,sy:ey,sz:ez) * reyfluc % wz (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % vy (sx:ex,sy:ey,sz:ez) * favfluc % wz (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_vywz )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_vywz )

    ! uxvywz
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % ux (sx:ex,sy:ey,sz:ez) * &
                               reyfluc % vy (sx:ex,sy:ey,sz:ez) * &
                               reyfluc % wz (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ux (sx:ex,sy:ey,sz:ez) * &
                               favfluc % vy (sx:ex,sy:ey,sz:ez) * &
                               favfluc % wz (sx:ex,sy:ey,sz:ez) * &
                               rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_uxvywz )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_uxvywz )

    ! rhorho
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % rho (sx:ex,sy:ey,sz:ez) * reyfluc % rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_rho )

    ! rho_acu rho_acu
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % p (sx:ex,sy:ey,sz:ez) /                                   &
                             ( reyavg % cs (sx:ex,sy:ey,sz:ez) * reyavg % cs (sx:ex,sy:ey,sz:ez) )
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * wrk1 (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_rho_acu )

    ! rho_ent rho_ent
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % p (sx:ex,sy:ey,sz:ez) /                                   &
                             ( reyavg % cs (sx:ex,sy:ey,sz:ez) * reyavg % cs (sx:ex,sy:ey,sz:ez) )
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % rho (sx:ex,sy:ey,sz:ez) - wrk1 (sx:ex,sy:ey,sz:ez)
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * wrk1 (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_rho_ent )

    ! pp
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % p (sx:ex,sy:ey,sz:ez) * reyfluc % p (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_p )

    ! TT
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % T (sx:ex,sy:ey,sz:ez) * reyfluc % T (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_T )

    ! T_acu T_acu
    call gammamix ( 0 , thd , inst % W_i , inst % cp , wrk1 )
    wrk1 (sx:ex,sy:ey,sz:ez) = ( wrk1 (sx:ex,sy:ey,sz:ez) - 1.0_dp ) / wrk1 (sx:ex,sy:ey,sz:ez)
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez) / inst % W_i (sx:ex,sy:ey,sz:ez)
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * reyfluc % p (sx:ex,sy:ey,sz:ez)
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * wrk1 (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_T_acu )

    ! T_ent T_ent
    call gammamix ( 0 , thd , inst % W_i , inst % cp , wrk1 )
    wrk1 (sx:ex,sy:ey,sz:ez) = ( wrk1 (sx:ex,sy:ey,sz:ez) - 1.0_dp ) / wrk1 (sx:ex,sy:ey,sz:ez)
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * rho_i (sx:ex,sy:ey,sz:ez) / inst % W_i (sx:ex,sy:ey,sz:ez)
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * reyfluc % p (sx:ex,sy:ey,sz:ez)
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % T (sx:ex,sy:ey,sz:ez) - wrk1 (sx:ex,sy:ey,sz:ez)
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * wrk1 (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_T_ent )

    ! WW
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % W (sx:ex,sy:ey,sz:ez) * reyfluc % W (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_W )

    ! du1/dx1 du1/dx1
    wrk1 (sx:ex,sy:ey,sz:ez) = favfluc % ux (sx:ex,sy:ey,sz:ez)
    call comm_one (wrk1) ; call dx ( dx_i , wrk1 , wrk2 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk2 (sx:ex,sy:ey,sz:ez) * wrk2 (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_du1_dx1 )

    ! pdu1/dx1
    wrk1 (sx:ex,sy:ey,sz:ez) = favfluc % ux (sx:ex,sy:ey,sz:ez)
    call comm_one (wrk1) ; call dx ( dx_i , wrk1 , wrk2 )
    wrk2 (sx:ex,sy:ey,sz:ez) = wrk2 (sx:ex,sy:ey,sz:ez) * reyfluc % p (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % reyvar_p_du1_dx1 )

    ! rho p
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % rho (sx:ex,sy:ey,sz:ez) * reyfluc % p (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_rho_p )

    ! rho W
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % rho (sx:ex,sy:ey,sz:ez) * reyfluc % W (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_rho_W )

    ! omega omega
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % omega (sx:ex,sy:ey,sz:ez) * reyfluc % omega (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % omega (sx:ex,sy:ey,sz:ez) * favfluc % omega (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_omega )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_omega )

    ! Da Da
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % Da (sx:ex,sy:ey,sz:ez) * reyfluc % Da (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % Da (sx:ex,sy:ey,sz:ez) * favfluc % Da (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_Da )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_Da )

    ! sdr sdr
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % sdr (sx:ex,sy:ey,sz:ez) * reyfluc % sdr (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % sdr (sx:ex,sy:ey,sz:ez) * favfluc % sdr (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_sdr )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_sdr )

    ! FI FI
    wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % FI (sx:ex,sy:ey,sz:ez) * reyfluc % FI (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % FI (sx:ex,sy:ey,sz:ez) * favfluc % FI (sx:ex,sy:ey,sz:ez) * rho (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % reyvar_FI )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % favvar_FI )


    ! scalar fluxes


    wrk1 (sx:ex,sy:ey,sz:ez) = favfluc % Zm (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % Y (sx:ex,sy:ey,sz:ez, max ( inp % index_scalar , 1 ) )


    ! ux
    wrk3 (sx:ex,sy:ey,sz:ez) = rho (sx:ex,sy:ey,sz:ez) * favfluc % ux (sx:ex,sy:ey,sz:ez) * &
                               wrk1 (sx:ex,sy:ey,sz:ez)
    wrk4 (sx:ex,sy:ey,sz:ez) = rho (sx:ex,sy:ey,sz:ez) * favfluc % ux (sx:ex,sy:ey,sz:ez) * &
                               wrk2 (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk3 )
    call spatZ_avg ( inp % spat_avg , dz , wrk4 )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk3 , stat % favavg_ux_Zm )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk4 , stat % favavg_ux_Y  )
    ! vy
    wrk3 (sx:ex,sy:ey,sz:ez) = rho (sx:ex,sy:ey,sz:ez) * favfluc % vy (sx:ex,sy:ey,sz:ez) * &
                               wrk1 (sx:ex,sy:ey,sz:ez)
    wrk4 (sx:ex,sy:ey,sz:ez) = rho (sx:ex,sy:ey,sz:ez) * favfluc % vy (sx:ex,sy:ey,sz:ez) * &
                               wrk2 (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk3 )
    call spatZ_avg ( inp % spat_avg , dz , wrk4 )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk3 , stat % favavg_vy_Zm )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk4 , stat % favavg_vy_Y  )
    ! wz
    wrk3 (sx:ex,sy:ey,sz:ez) = rho (sx:ex,sy:ey,sz:ez) * favfluc % wz (sx:ex,sy:ey,sz:ez) * &
                               wrk1 (sx:ex,sy:ey,sz:ez)
    wrk4 (sx:ex,sy:ey,sz:ez) = rho (sx:ex,sy:ey,sz:ez) * favfluc % wz (sx:ex,sy:ey,sz:ez) * &
                               wrk2 (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk3 )
    call spatZ_avg ( inp % spat_avg , dz , wrk4 )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk3 , stat % favavg_wz_Zm )
    call temp_avg ( fav_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk4 , stat % favavg_wz_Y  )


    ! conditional variances


    if ( inp % cond_avg ) then

    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if

    ! T_Zm T_Zm
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ca_T_Zm (sx:ex,sy:ey,sz:ez) * &
                               favfluc % ca_T_Zm (sx:ex,sy:ey,sz:ez) * &
                               rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk2 , reyavg % ca_rho_Zm , stat % favvar_ca_T_Zm )

    ! c_Zm c_Zm


    ! sdr_Zm sdr_Zm
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ca_sdr_Zm (sx:ex,sy:ey,sz:ez) * &
                               favfluc % ca_sdr_Zm (sx:ex,sy:ey,sz:ez) * &
                               rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk2 , reyavg % ca_rho_Zm , stat % favvar_ca_sdr_Zm )

    ! Da_Zm Da_Zm
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ca_Da_Zm (sx:ex,sy:ey,sz:ez) * &
                               favfluc % ca_Da_Zm (sx:ex,sy:ey,sz:ez) * &
                               rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk2 , reyavg % ca_rho_Zm , stat % favvar_ca_Da_Zm )

    ! omega_Zm omega_Zm
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ca_omega_Zm (sx:ex,sy:ey,sz:ez) * &
                               favfluc % ca_omega_Zm (sx:ex,sy:ey,sz:ez) * &
                               rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk2 , reyavg % ca_rho_Zm , stat % favvar_ca_omega_Zm )

    ! FI_Zm FI_Zm
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ca_FI_Zm (sx:ex,sy:ey,sz:ez) * &
                               favfluc % ca_FI_Zm (sx:ex,sy:ey,sz:ez) * &
                               rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_Y , max_Y ,                                           &
                    wrk1 , wrk2 , reyavg % ca_rho_Zm , stat % favvar_ca_FI_Zm )


    ! FI
    call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk1 )

    ! omega_FI omega_FI
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ca_omega_FI (sx:ex,sy:ey,sz:ez) * &
                               favfluc % ca_omega_FI (sx:ex,sy:ey,sz:ez) * &
                               rho (sx:ex,sy:ey,sz:ez)
    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    min_FI , max_FI ,                                         &
                    wrk1 , wrk2 , reyavg % ca_rho_FI , stat % favvar_ca_omega_FI )


    ! Da
    call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u , inst % T , inst % W_i , &
                     inst % dm , inst % tdr , inst % rd , inst % Ydot , wrk1 )
    wrk1 (sx:ex,sy:ey,sz:ez) = log10 ( wrk1 (sx:ex,sy:ey,sz:ez) )

    ! omega_Da omega_Da
    wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % ca_omega_Da (sx:ex,sy:ey,sz:ez) * &
                               favfluc % ca_omega_Da (sx:ex,sy:ey,sz:ez) * &
                               rho (sx:ex,sy:ey,sz:ez)

    call cond_avg ( fav_avg , ifile , end_file , dtime , sumdtime , dy , dz , &
                    log10 (min_Da_adim) , log10 (max_Da_adim) ,               &
                    wrk1 , wrk2 , reyavg % ca_rho_Da , stat % favvar_ca_omega_Da )

    end if



    ! turbulent kinetic energy and Reynolds transport equations


    name = 'TKE_Convection' ; call detect_name ( inp , name , loop )
    if (loop) then

    ! dissipation
    call tke_dissipation ( dx_i , dy_i , dz_i , favfluc % ux , favfluc % vy , favfluc % wz , &
                           reyfluc % tau , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % tke_dissip )
    call rey_dissipation_11 ( dx_i , dy_i , dz_i , favfluc % ux , &
                              reyfluc % tau , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_dissip_11 )
    call rey_dissipation_22 ( dx_i , dy_i , dz_i , favfluc % vy , &
                              reyfluc % tau , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_dissip_22 )
    call rey_dissipation_33 ( dx_i , dy_i , dz_i , favfluc % wz , &
                              reyfluc % tau , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_dissip_33 )
    call rey_dissipation_12 ( dx_i , dy_i , dz_i , favfluc % ux , favfluc % vy , &
                              reyfluc % tau , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_dissip_12 )
    call rey_dissipation_13 ( dx_i , dy_i , dz_i , favfluc % ux , favfluc % wz , &
                              reyfluc % tau , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_dissip_13 )
    call rey_dissipation_23 ( dx_i , dy_i , dz_i , favfluc % vy , favfluc % wz , &
                              reyfluc % tau , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_dissip_23 )


    ! dissipation last iteration
    if ( ifile == end_file ) then
       stat % tke_dissip (sx:ex,sy:ey,sz:ez) = stat % tke_dissip (sx:ex,sy:ey,sz:ez) / &
                                               reyavg % rho (sx:ex,sy:ey,sz:ez)
       stat % rey_dissip_11 (sx:ex,sy:ey,sz:ez) = stat % rey_dissip_11 (sx:ex,sy:ey,sz:ez) / &
                                                  reyavg % rho (sx:ex,sy:ey,sz:ez)
       stat % rey_dissip_22 (sx:ex,sy:ey,sz:ez) = stat % rey_dissip_22 (sx:ex,sy:ey,sz:ez) / &
                                                  reyavg % rho (sx:ex,sy:ey,sz:ez)
       stat % rey_dissip_33 (sx:ex,sy:ey,sz:ez) = stat % rey_dissip_33 (sx:ex,sy:ey,sz:ez) / &
                                                  reyavg % rho (sx:ex,sy:ey,sz:ez)
       stat % rey_dissip_12 (sx:ex,sy:ey,sz:ez) = stat % rey_dissip_12 (sx:ex,sy:ey,sz:ez) / &
                                                  reyavg % rho (sx:ex,sy:ey,sz:ez)
       stat % rey_dissip_13 (sx:ex,sy:ey,sz:ez) = stat % rey_dissip_13 (sx:ex,sy:ey,sz:ez) / &
                                                  reyavg % rho (sx:ex,sy:ey,sz:ez)
       stat % rey_dissip_23 (sx:ex,sy:ey,sz:ez) = stat % rey_dissip_23 (sx:ex,sy:ey,sz:ez) / &
                                                  reyavg % rho (sx:ex,sy:ey,sz:ez)
    end if


    ! pressure-strain
    call tke_pressure_strain ( dx_i , dy_i , dz_i , favfluc % ux , favfluc % vy , favfluc % wz , &
                               reyfluc % P , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % tke_press_strain )
    call rey_pressure_strain_11 ( dx_i , favfluc % ux , &
                                  reyfluc % P , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_press_strain_11 )
    call rey_pressure_strain_22 ( dy_i , favfluc % vy , &
                                  reyfluc % P , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_press_strain_22 )
    call rey_pressure_strain_33 ( dz_i , favfluc % wz , &
                                  reyfluc % P , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_press_strain_33 )
    call rey_pressure_strain_12 ( dx_i , dy_i , favfluc % ux , favfluc % vy , &
                                  reyfluc % P , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_press_strain_12 )
    call rey_pressure_strain_13 ( dx_i , dz_i , favfluc % ux , favfluc % wz , &
                                  reyfluc % P , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_press_strain_13 )
    call rey_pressure_strain_23 ( dy_i , dz_i , favfluc % vy , favfluc % wz , &
                                  reyfluc % P , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % rey_press_strain_23 )


    ! transport
    call tke_transport1 ( adi , reyfluc % ux , reyfluc % vy , reyfluc % wz , &
                          favfluc % ux , favfluc % vy , favfluc % wz ,       &
                          rho , reyfluc % P , reyfluc % tau , wrk_dim4 )
    ! 1
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % tke_transp (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % tke_transp (sx:ex,sy:ey,sz:ez,1) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 2
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % tke_transp (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % tke_transp (sx:ex,sy:ey,sz:ez,2) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 3
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % tke_transp (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % tke_transp (sx:ex,sy:ey,sz:ez,3) = wrk2 (sx:ex,sy:ey,sz:ez)

    call rey_transport1_11 ( adi , reyfluc % ux ,                           &
                             favfluc % ux , favfluc % vy , favfluc % wz ,   &
                             rho , reyfluc % P , reyfluc % tau , wrk_dim4 )
    ! 1
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_11 (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_11 (sx:ex,sy:ey,sz:ez,1) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 2
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_11 (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_11 (sx:ex,sy:ey,sz:ez,2) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 3
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_11 (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_11 (sx:ex,sy:ey,sz:ez,3) = wrk2 (sx:ex,sy:ey,sz:ez)

    call rey_transport1_22 ( adi , reyfluc % vy ,                           &
                             favfluc % ux , favfluc % vy , favfluc % wz ,   &
                             rho , reyfluc % P , reyfluc % tau , wrk_dim4 )
    ! 1
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_22 (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_22 (sx:ex,sy:ey,sz:ez,1) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 2
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_22 (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_22 (sx:ex,sy:ey,sz:ez,2) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 3
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_22 (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_22 (sx:ex,sy:ey,sz:ez,3) = wrk2 (sx:ex,sy:ey,sz:ez)

    call rey_transport1_33 ( adi , reyfluc % wz ,                           &
                             favfluc % ux , favfluc % vy , favfluc % wz ,   &
                             rho , reyfluc % P , reyfluc % tau , wrk_dim4 )
    ! 1
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_33 (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_33 (sx:ex,sy:ey,sz:ez,1) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 2
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_33 (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_33 (sx:ex,sy:ey,sz:ez,2) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 3
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_33 (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_33 (sx:ex,sy:ey,sz:ez,3) = wrk2 (sx:ex,sy:ey,sz:ez)

    call rey_transport1_12 ( adi , reyfluc % ux , reyfluc % vy ,            &
                             favfluc % ux , favfluc % vy , favfluc % wz ,   &
                             rho , reyfluc % P , reyfluc % tau , wrk_dim4 )
    ! 1
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_12 (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_12 (sx:ex,sy:ey,sz:ez,1) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 2
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_12 (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_12 (sx:ex,sy:ey,sz:ez,2) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 3
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_12 (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_12 (sx:ex,sy:ey,sz:ez,3) = wrk2 (sx:ex,sy:ey,sz:ez)

    call rey_transport1_13 ( adi , reyfluc % ux , reyfluc % wz ,            &
                             favfluc % ux , favfluc % vy , favfluc % wz ,   &
                             rho , reyfluc % P , reyfluc % tau , wrk_dim4 )
    ! 1
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_13 (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_13 (sx:ex,sy:ey,sz:ez,1) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 2
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_13 (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_13 (sx:ex,sy:ey,sz:ez,2) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 3
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_13 (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_13 (sx:ex,sy:ey,sz:ez,3) = wrk2 (sx:ex,sy:ey,sz:ez)

    call rey_transport1_23 ( adi , reyfluc % vy , reyfluc % wz ,            &
                             favfluc % ux , favfluc % vy , favfluc % wz ,   &
                             rho , reyfluc % P , reyfluc % tau , wrk_dim4 )
    ! 1
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_23 (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_23 (sx:ex,sy:ey,sz:ez,1) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 2
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_23 (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_23 (sx:ex,sy:ey,sz:ez,2) = wrk2 (sx:ex,sy:ey,sz:ez)
    ! 3
    wrk1 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk2 (sx:ex,sy:ey,sz:ez) = stat % rey_transp_23 (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , wrk2 )
    stat % rey_transp_23 (sx:ex,sy:ey,sz:ez,3) = wrk2 (sx:ex,sy:ey,sz:ez)


    ! mass flux coupling
    ! ux
    wrk1 (sx:ex,sy:ey,sz:ez) = favfluc % ux (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % tke_ra_ff_ux )
    ! vy
    wrk1 (sx:ex,sy:ey,sz:ez) = favfluc % vy (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % tke_ra_ff_vy )
    ! wz
    wrk1 (sx:ex,sy:ey,sz:ez) = favfluc % wz (sx:ex,sy:ey,sz:ez)
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % tke_ra_ff_wz )

    end if


    ! vorticity transport equation
    name = 'ENS_Convection' ; call detect_name ( inp , name , loop )
    if (loop) then

    ! vorticity vector
    call vorticity_i ( dx_i , dy_i , dz_i , inst % u , wrk_dim4 )

    ! convection
    call vorticity_convection ( dx_i , dy_i , dz_i , inst % u , wrk_dim4 , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % vor_conv )

    call enstrophy_convection ( dx_i , dy_i , dz_i , inst % u , wrk_dim4 , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % ens_conv )

    ! stretching
    call vorticity_stretching ( dx_i , dy_i , dz_i , inst % u , wrk_dim4 , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % vor_stret )

    call enstrophy_stretching ( dx_i , dy_i , dz_i , inst % u , wrk_dim4 , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % ens_stret )

    ! dilatation
    call vorticity_dilatation ( dx_i , dy_i , dz_i , inst % u , wrk_dim4 , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % vor_dila )

    call enstrophy_dilatation ( dx_i , dy_i , dz_i , inst % u , wrk_dim4 , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % ens_dila )

    ! baroclinic
    call vorticity_baroclinic ( dx_i , dy_i , inst % W_i , inst % T , inst % u , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % vor_baro )

    call enstrophy_baroclinic ( dx_i , dy_i , dz , inst % W_i , inst % T , inst % u , wrk_dim4 , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % ens_baro )

    ! viscous

    ! tau
    wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) * rho_i (sx:ex,sy:ey,sz:ez)
    wrk3 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) * rho_i (sx:ex,sy:ey,sz:ez)
    call tau_ij ( dx_i , dy_i , dz_i , wrk1 , wrk2 , wrk3 , inst % mu , inst % kpa , wrk_dim5 )

    call vorticity_viscous ( dx_i , dy_i , dz_i , inst % u , wrk_dim5 , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % vor_visc )

    call enstrophy_viscous ( dx_i , dy_i , dz_i , inst % u , wrk_dim4 , wrk_dim5 , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % ens_visc )
    
    ! volume viscosity production
    call vorticity_volume_viscosity_production ( dx_i , dy_i , dz_i , inst % u , inst % kpa , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % vor_volume )    

    ! shear viscosity production
    call vorticity_shear_viscosity_production ( dx_i , dy_i , dz_i , inst % u , inst % mu , wrk1 )
    call spatZ_avg ( inp % spat_avg , dz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % vor_shear )      

    end if


    ! scalar variance transport equation
    name = 'SCAL_VAR_Mean_Convection' ; call detect_name ( inp , name , loop )
    if (loop) then

    if ( inp % index_scalar == 0 ) then
       call Zmsimp ( inp , thd , inst % u , wrk1 )                 ! wrk1 is the instantaneous value
       wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % Zm (sx:ex,sy:ey,sz:ez) ! wrk2 is the fluctuating value
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % Zm (sx:ex,sy:ey,sz:ez)  ! wrk3 is the averaged value
    else
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv + inp % index_scalar ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = favfluc % Y (sx:ex,sy:ey,sz:ez, inp % index_scalar )
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez, inp % index_scalar )
    end if

    ! equivalent diffusion coefficient: Deq
    call Zmsimp ( inp , thd , inst % u , wrk4 )
    call D_equivalent ( inp , thd                           , &
                        dx_i      , dy_i       , dz_i       , &
                        inst % u  , inst % T   , inst % W_i , &
                        inst % ct , inst % cp  ,              &
                        inst % dm , inst % tdr , inst % rd  , &
                        wrk4      , wrk5 )
    wrk4 = wrk5 ! Deq is wrk4


    ! turbulent transport
    call scal_var_turb_transp ( rho , favfluc % ux , favfluc % vy , favfluc % wz , wrk2 , wrk_dim4 )

    ! turbulent transport 1
    wrk5 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk6 (sx:ex,sy:ey,sz:ez) = stat % scal_var_turb_transp (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , wrk6 )
    stat % scal_var_turb_transp (sx:ex,sy:ey,sz:ez,1) = wrk6 (sx:ex,sy:ey,sz:ez)

    ! turbulent transport 2
    wrk5 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk6 (sx:ex,sy:ey,sz:ez) = stat % scal_var_turb_transp (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , wrk6 )
    stat % scal_var_turb_transp (sx:ex,sy:ey,sz:ez,2) = wrk6 (sx:ex,sy:ey,sz:ez)

    ! turbulent transport 3
    wrk5 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk6 (sx:ex,sy:ey,sz:ez) = stat % scal_var_turb_transp (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , wrk6 )
    stat % scal_var_turb_transp (sx:ex,sy:ey,sz:ez,3) = wrk6 (sx:ex,sy:ey,sz:ez)


    ! turbulent production
    call scal_var_turb_prod ( rho , favfluc % ux , favfluc % vy , favfluc % wz , wrk2 , wrk_dim4 )

    ! turbulent production 1
    wrk5 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk6 (sx:ex,sy:ey,sz:ez) = stat % scal_var_turb_prod (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , wrk6 )
    stat % scal_var_turb_prod (sx:ex,sy:ey,sz:ez,1) = wrk6 (sx:ex,sy:ey,sz:ez)

    ! turbulent production 2
    wrk5 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk6 (sx:ex,sy:ey,sz:ez) = stat % scal_var_turb_prod (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , wrk6 )
    stat % scal_var_turb_prod (sx:ex,sy:ey,sz:ez,2) = wrk6 (sx:ex,sy:ey,sz:ez)

    ! turbulent production 3
    wrk5 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk6 (sx:ex,sy:ey,sz:ez) = stat % scal_var_turb_prod (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , wrk6 )
    stat % scal_var_turb_prod (sx:ex,sy:ey,sz:ez,3) = wrk6 (sx:ex,sy:ey,sz:ez)


    ! molecular diffusion
    call scal_var_mol_diff ( inp % index_scalar , thd , dx_i , dy_i , dz_i , inst % T , inst % W_i , &
                             inst % dm , inst % tdr , inst % rd , inst % u , wrk4 , wrk2 , wrk1 , wrk_dim4 )

    ! molecular diffusion 1
    wrk5 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,1)
    wrk6 (sx:ex,sy:ey,sz:ez) = stat % scal_var_mol_diff (sx:ex,sy:ey,sz:ez,1)
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , wrk6 )
    stat % scal_var_mol_diff (sx:ex,sy:ey,sz:ez,1) = wrk6 (sx:ex,sy:ey,sz:ez)

    ! molecular diffusion 2
    wrk5 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,2)
    wrk6 (sx:ex,sy:ey,sz:ez) = stat % scal_var_mol_diff (sx:ex,sy:ey,sz:ez,2)
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , wrk6 )
    stat % scal_var_mol_diff (sx:ex,sy:ey,sz:ez,2) = wrk6 (sx:ex,sy:ey,sz:ez)

    ! molecular diffusion 3
    wrk5 (sx:ex,sy:ey,sz:ez) = wrk_dim4 (sx:ex,sy:ey,sz:ez,3)
    wrk6 (sx:ex,sy:ey,sz:ez) = stat % scal_var_mol_diff (sx:ex,sy:ey,sz:ez,3)
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , wrk6 )
    stat % scal_var_mol_diff (sx:ex,sy:ey,sz:ez,3) = wrk6 (sx:ex,sy:ey,sz:ez)


    ! dissipation (actual term)
    call scal_var_dissip ( inp % index_scalar , thd , dx_i , dy_i , dz_i , inst % T , inst % W_i , &
                           inst % dm , inst % tdr , inst % rd , inst % u , wrk4 , wrk2 , wrk1 , wrk5 )
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , stat % scal_var_dissip )

    ! dissipation (instantaneous-average term)
    call scal_var_dissip ( inp % index_scalar , thd , dx_i , dy_i , dz_i , inst % T , inst % W_i , &
                           inst % dm , inst % tdr , inst % rd , inst % u , wrk4 , wrk3 , wrk1 , wrk5 )
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , stat % scal_var_dissip2 )

    ! dissipation (instantaneous-instantaneous term)
    call scal_var_dissip ( inp % index_scalar , thd , dx_i , dy_i , dz_i , inst % T , inst % W_i , &
                           inst % dm , inst % tdr , inst % rd , inst % u , wrk4 , wrk1 , wrk1 , wrk5 )
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , stat % scal_var_dissip3 )

    ! dissipation (average-average term)
    call scal_var_dissip ( inp % index_scalar , thd , dx_i , dy_i , dz_i , inst % T , inst % W_i , &
                           inst % dm , inst % tdr , inst % rd , inst % u , wrk4 , wrk3 , wrk3 , wrk5 )
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , stat % scal_var_dissip4 )

    ! dissipation (fluctuating-average term)
    call scal_var_dissip ( inp % index_scalar , thd , dx_i , dy_i , dz_i , inst % T , inst % W_i , &
                           inst % dm , inst % tdr , inst % rd , inst % u , wrk4 , wrk2 , wrk3 , wrk5 )
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , stat % scal_var_dissip5 )

    ! dissipation (fluctuating-fluctuating term)
    call scal_var_dissip ( inp % index_scalar , thd , dx_i , dy_i , dz_i , inst % T , inst % W_i , &
                           inst % dm , inst % tdr , inst % rd , inst % u , wrk4 , wrk2 , wrk2 , wrk5 )
    call spatZ_avg ( inp % spat_avg , dz , wrk5 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk5 , stat % scal_var_dissip6 )


    end if


    ! correlation functions


    name = 'CORR_ux' ; call detect_name ( inp , name , loop )
    if (loop) then

    ! ux
    call correlation_Z ( inp , xt , yt , zt , reyfluc % ux , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % tmp_corr_ux )
    if ( ifile == end_file ) call last_correlation_Z ( inp , xt , yt , zt , stat % tmp_corr_ux )

    ! vy
    call correlation_Z ( inp , xt , yt , zt , reyfluc % vy , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % tmp_corr_vy )
    if ( ifile == end_file ) call last_correlation_Z ( inp , xt , yt , zt , stat % tmp_corr_vy )

    ! wz
    call correlation_Z ( inp , xt , yt , zt , reyfluc % wz , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % tmp_corr_wz )
    if ( ifile == end_file ) call last_correlation_Z ( inp , xt , yt , zt , stat % tmp_corr_wz )

    ! p
    call correlation_Z ( inp , xt , yt , zt , reyfluc % p , wrk1 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk1 , stat % tmp_corr_p )
    if ( ifile == end_file ) call last_correlation_Z ( inp , xt , yt , zt , stat % tmp_corr_p )

    ! Y
    if ( inp % index_scalar == 0 ) then
       wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % Zm (sx:ex,sy:ey,sz:ez)
    else
       wrk1 (sx:ex,sy:ey,sz:ez) = reyfluc % Y (sx:ex,sy:ey,sz:ez, inp % index_scalar )
    end if
    call correlation_Z ( inp , xt , yt , zt , wrk1 , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % tmp_corr_Y )
    if ( ifile == end_file ) call last_correlation_Z ( inp , xt , yt , zt , stat % tmp_corr_Y )

    ! Zm
    if ( npv == 0 ) then ! calculated passive scalar
       call Zmsimp ( inp , thd , inst % u , wrk1 )
    else                 ! real passive scalar
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) * &
                                  rho_i (sx:ex,sy:ey,sz:ez)
    end if
    call correlation_Z ( inp , xt , yt , zt , wrk1 , wrk2 )
    call temp_avg ( rey_avg , ifile , end_file , dtime , sumdtime , reyavg % rho , wrk2 , stat % tmp_corr_Zm )
    if ( ifile == end_file ) call last_correlation_Z ( inp , xt , yt , zt , stat % tmp_corr_Zm )

    end if


    deallocate ( dy , dz , rho , rho_i )
    deallocate ( wrk1 , wrk2 , wrk3 , wrk4 )
    deallocate ( wrk_dim2 , wrk_dim4 , wrk_dim5 )


  end subroutine stats


!> \brief Calculate the temporal average.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine temp_avg ( type_avg , ifile , end_file , dt , sumdt , avg_rho , inst_var , avg_var )


    character (len_default) , intent (in)                              :: type_avg !< average type
    integer (ip) , intent (in)                                         :: ifile    !< file number
    integer (ip) , intent (in)                                         :: end_file !< last file number
    real (dp) , intent (in)                                            :: dt       !< time step
    real (dp) , intent (in)                                            :: sumdt    !< sum of time steps
    real (dp) , dimension (:,:,:) , allocatable , intent (in)          :: avg_rho  !< averaged density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)          :: inst_var !< instantaneous variable
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)       :: avg_var  !< averaged variable


    integer (ip) :: i , j , k
    real (dp)    :: sumdt_i


    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             avg_var (i,j,k) = avg_var (i,j,k) + inst_var (i,j,k) * dt
          end do
       end do
    end do


    if ( ifile == end_file ) then


       if ( type_avg == rey_avg ) then

          sumdt_i = 1.0_dp / sumdt
          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   avg_var (i,j,k) = avg_var (i,j,k) * sumdt_i
                end do
             end do
          end do

       else if ( type_avg == fav_avg ) then

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   avg_var (i,j,k) = avg_var (i,j,k) / ( sumdt * avg_rho (i,j,k) )
                end do
             end do
          end do

       end if


    end if


  end subroutine temp_avg


!> \brief Calculate the spatial average in z-direction.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine spatZ_avg ( spat_avg , dz , v )


    logical , intent (in)                                              :: spat_avg !< activate-deactivate spatial average
    real (dp) , dimension (:) , allocatable , intent (in)              :: dz       !< dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)       :: v        !< variable


    integer (ip)                                    :: ok , nxy , index , i , j , k
    real (dp)                                       :: sumdz_i
    real (dp) , dimension (:) , allocatable         :: wrk1 , wrk2


    if ( .not. spat_avg ) return ! leave the subroutine


    nxy = (ex-sx+1) * (ey-sy+1)

    allocate ( wrk1 ( nxy+1 ) , &
               wrk2 ( nxy+1 ) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate spatZ_avg')

    index    = 1
    wrk1 (:) = 0.0_dp

    do j = sy,ey
       do i = sx,ex
          do k = sz,ez
             wrk1 (index) = wrk1 (index) + v (i,j,k) * dz (k)
          end do
          index = index + 1
       end do
    end do
    wrk1 (index) = sum ( dz(sz:ez) )

    if ( comm1dz /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk1 , wrk2 , nxy+1 , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dz , mpicode )
    else
       wrk2 (:) = wrk1 (:)
    end if
    sumdz_i = 1.0_dp / wrk2 (nxy+1)

    index = 1
    do j = sy,ey
       do i = sx,ex
          v (i,j,sz:ez) = wrk2 (index) * sumdz_i
          index         = index + 1
       end do
    end do

    deallocate ( wrk1 , wrk2 )


  end subroutine spatZ_avg


!> \brief Calculate the spatial average in yz-plane.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine spatYZ_avg ( spat_avg , dy , dz , v )


    logical , intent (in)                                              :: spat_avg !< activate-deactivate spatial average
    real (dp) , dimension (:) , allocatable , intent (in)              :: dy       !< dy array
    real (dp) , dimension (:) , allocatable , intent (in)              :: dz       !< dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)       :: v        !< variable


    integer (ip)                                    :: ok , nxz , nx , index1 , index2 , i , j , k
    real (dp)                                       :: sumdy_i , sumdz_i
    real (dp) , dimension (:) , allocatable         :: wrk1 , wrk2 , wrk3 , wrk4


    if ( .not. spat_avg ) return ! leave the subroutine


    nx  = (ex-sx+1)
    nxz = (ex-sx+1) * (ez-sz+1)
    allocate ( wrk1 ( nxz+1 ) , &
               wrk2 ( nxz+1 ) , &
               wrk3 ( nx+1 )  , &
               wrk4 ( nx+1 )  , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate spatYZ_avg')


    ! averaging in y-direction
    index1   = 1
    wrk1 (:) = 0.0_dp
    do i = sx,ex
       do k = sz,ez
          do j = sy,ey
             wrk1 (index1) = wrk1 (index1) + v (i,j,k) * dy (j)
          end do
          index1 = index1 + 1
       end do
    end do
    wrk1 (nxz+1) = sum ( dy(sy:ey) )

    if ( comm1dy /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk1 , wrk2 , nxz+1 , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dy , mpicode )
    else
       wrk2 (:) = wrk1 (:)
    end if
    sumdy_i = 1.0_dp / wrk2 (nxz+1)

    index1 = 1
    do i = sx,ex
       do k = sz,ez
          wrk2 (index1) = wrk2 (index1) * sumdy_i
          index1        = index1 + 1
       end do
    end do


    ! averaging in z-direction
    index1 = 1 ; index2 = 1
    wrk3 (:) = 0.0_dp
    do i = sx,ex
       do k = sz,ez
          wrk3 (index2) = wrk3 (index2) + wrk2 (index1) * dz (k)
          index1 = index1 + 1
       end do
       index2 = index2 + 1
    end do
    wrk3 (nx+1) = sum ( dz(sz:ez) )

    if ( comm1dz /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk3 , wrk4 , nx+1 , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dz , mpicode )
    else
       wrk4 (:) = wrk3 (:)
    end if
    sumdz_i = 1.0_dp / wrk4 (nx+1)

    index2 = 1
    do i = sx,ex
       v (i,sy:ey,sz:ez) = wrk4 (index2) * sumdz_i
       index2            = index2 + 1
    end do


    deallocate ( wrk1 , wrk2 , wrk3 , wrk4 )


  end subroutine spatYZ_avg


!> \brief Calculate the conditional average.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine cond_avg ( type_avg , ifile , end_file , dt , sumdt , dy , dz , &
                        min_value , max_value ,                              &
                        cond_var , inst_var , avg_rho , avg_var )


    character (len_default) , intent (in)                              :: type_avg  !< average type
    integer (ip) , intent (in)                                         :: ifile     !< file number
    integer (ip) , intent (in)                                         :: end_file  !< last file number
    real (dp) , intent (in)                                            :: dt        !< time step
    real (dp) , intent (in)                                            :: sumdt     !< sum of time steps
    real (dp) , dimension (:) , allocatable , intent (in)              :: dy        !< dy array
    real (dp) , dimension (:) , allocatable , intent (in)              :: dz        !< dz array
    real (dp) , intent (in)                                            :: min_value !< minimum value of the conditional variable
    real (dp) , intent (in)                                            :: max_value !< maximum value of the conditional variable
    real (dp) , dimension (:,:,:) , allocatable , intent (in)          :: cond_var  !< conditional variable
    real (dp) , dimension (:,:,:) , allocatable , intent (in)          :: inst_var  !< instantaneous variable
    real (dp) , dimension (:,:) , allocatable , intent (in)            :: avg_rho   !< averaged density
    real (dp) , dimension (:,:) , allocatable , intent (inout)         :: avg_var   !< conditional averaged variable


    integer (ip)                                    :: ok , index1 , index2 , i , j , k , n
    integer (ip)                                    :: nxzbins , nxbins
    real (dp)                                       :: sumdt_i , wrk
    real (dp)                                       :: delta_cond , half_interval , mytarget , error
    real (dp) , dimension (:) , allocatable         :: wrk1 , wrk2 , wrk3 , wrk4


    delta_cond    = ( max_value - min_value ) / ( nbins - 1 )
    half_interval = 0.5_dp * delta_cond


    ! y-direction spatial homogeneization
    nxzbins = (ex-sx+1) * nbins * (ez-sz+1)
    allocate ( wrk1 (nxzbins+nxzbins) , &
               wrk2 (nxzbins+nxzbins) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate cond_avg')

    index1 = 1
    do i = sx,ex ! loop for each point in x-direction

       do n = 1,nbins ! loop for each bin in the cond_var space

          mytarget = min_value + (n-1) * delta_cond

          do k = sz,ez ! loop for each point in z-direction

             wrk1 (index1)         = 0.0_dp
             wrk1 (index1+nxzbins) = 0.0_dp

             do j = sy,ey ! loop for each point in y-direction
                error = abs ( cond_var (i,j,k) - mytarget )
                if ( error < half_interval ) then
                   wrk1 (index1)         = wrk1 (index1) + inst_var (i,j,k) * dy (j)
                   wrk1 (index1+nxzbins) = wrk1 (index1+nxzbins) + dy (j)
                end if
             end do

             index1 = index1 + 1 ! index of nxzbins array

          end do

       end do

    end do


    if ( comm1dy /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk1 , wrk2 , nxzbins+nxzbins , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dy , mpicode )
    else
       wrk2 (:) = wrk1 (:)
    end if


    ! avoid divisions by zero
    do n = 1,nxzbins
       if ( wrk2 (nxzbins+n) /= 0.0_dp ) wrk2 (n) = wrk2 (n) / wrk2 (nxzbins+n)
    end do


    ! sumdt_i = 0.0_dp
    ! do n = 1,nxzbins
    !    sumdt_i = sumdt_i + wrk2 (n)
    ! end do
    ! write (*,*) rank , sumdt_i
    ! stop


    ! z-direction spatial homogeneization
    nxbins = (ex-sx+1) * nbins
    allocate ( wrk3 (nxbins+nxbins) , &
               wrk4 (nxbins+nxbins) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate cond_avg_2')

    index1 = 1
    index2 = 1
    do i = sx,ex ! loop for each point in x-direction
       do n = 1,nbins ! loop for each bin in the cond_var space

          wrk3 (index2)        = 0.0_dp
          wrk3 (index2+nxbins) = 0.0_dp

          do k = sz,ez ! loop for each point in z-direction

             if ( wrk2 (index1) /= 0.0_dp ) then
                wrk3 (index2)        = wrk3 (index2) + wrk2 (index1) * dz (k)
                wrk3 (index2+nxbins) = wrk3 (index2+nxbins) + dz (k)
             end if

             index1 = index1 + 1 ! index of nxzbins array

          end do

          index2 = index2 + 1 ! index of nxbins array

       end do
    end do


    if ( comm1dz /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk3 , wrk4 , nxbins+nxbins , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dz , mpicode )
    else
       wrk4 (:) = wrk3 (:)
    end if


    ! avoid divisions by zero
    do n = 1,nxbins
       if ( wrk4 (n) /= 0.0_dp ) wrk4 (n) = wrk4 (n) / wrk4 (nxbins+n)
    end do


    ! sumdt_i = 0.0_dp
    ! do n = 1,nxbins
    !    sumdt_i = sumdt_i + wrk4 (n)
    ! end do
    ! write (*,*) rank , sumdt_i
    ! call mpi_barrier ( MPI_COMM_WORLD , mpicode )
    ! stop


    ! temporal average
    index2 = 1
    do i = sx,ex
       do n = 1,nbins
          avg_var (i,n) = avg_var (i,n) + wrk4 (index2) * dt
          index2        = index2 + 1
       end do
    end do


    if ( ifile == end_file ) then

       if ( type_avg == rey_avg ) then

          sumdt_i = 1.0_dp / sumdt
          do n = 1,nbins
             do i = sx,ex
                avg_var (i,n) = avg_var (i,n) * sumdt_i
             end do
          end do

       else if ( type_avg == fav_avg ) then

          do n = 1,nbins
             do i = sx,ex
                wrk = 1.0_dp / max ( epsi , sumdt * avg_rho (i,n) )
                avg_var (i,n) = avg_var (i,n) * wrk
             end do
          end do

       end if

    end if


    deallocate ( wrk1 , wrk2 , wrk3 , wrk4 )


  end subroutine cond_avg


!> \brief Calculate the conditional fluctuation.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine cond_fluc ( min_value , max_value , cond_var , inst_var , avg_var , fluc_var )


    real (dp) , intent (in)                                            :: min_value !< minimum value of the conditional variable
    real (dp) , intent (in)                                            :: max_value !< maximum value of the conditional variable
    real (dp) , dimension (:,:,:) , allocatable , intent (in)          :: cond_var  !< conditional variable
    real (dp) , dimension (:,:,:) , allocatable , intent (in)          :: inst_var  !< instantaneous variable
    real (dp) , dimension (:,:) , allocatable , intent (in)            :: avg_var   !< conditional averaged variable
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)       :: fluc_var  !< conditional fluctuation variable


    logical                                         :: condition
    integer (ip)                                    :: i , j , k , n
    real (dp)                                       :: match_value , half_interval , delta_cond , mytarget , error


    delta_cond    = ( max_value - min_value ) / ( nbins - 1 )
    half_interval = 0.5_dp * delta_cond


    do k = sz,ez
       do j = sy,ey
          do i = sx,ex

             n           = 1
             condition   = .true.
             match_value = 0.0_dp

             do while (condition)

                if ( n == nbins ) condition = .false.

                mytarget = min_value + (n-1) * delta_cond
                error    = abs ( cond_var (i,j,k) - mytarget )

                if ( error < half_interval ) then
                   condition   = .false.
                   match_value = avg_var (i,n)
                end if

                n = n+1

             end do

             fluc_var (i,j,k) = inst_var (i,j,k) - match_value

          end do
       end do
    end do


  end subroutine cond_fluc


!> \brief Calculate the probability density function.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine pdf ( type_scale , ifile , end_file , min_value , max_value , var , pdf_var )


    character (len_default) , intent (in)                              :: type_scale !< pdf abscise scale
    integer (ip) , intent (in)                                         :: ifile      !< file number
    integer (ip) , intent (in)                                         :: end_file   !< last file number
    real (dp) , intent (in)                                            :: min_value  !< minimum value of the pdf interval
    real (dp) , intent (in)                                            :: max_value  !< maximum value of the pdf interval
    real (dp) , dimension (:,:,:) , allocatable , intent (in)          :: var        !< variable
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)       :: pdf_var    !< pdf variable


    integer (ip)                                    :: ok , index , i , j , k , n
    integer (ip)                                    :: nxybins
    real (dp)                                       :: wrk
    real (dp)                                       :: delta_cond , half_interval , mytarget , error
    real (dp) , dimension (:) , allocatable         :: wrk1 , wrk2


    delta_cond    = ( max_value - min_value ) / ( nbins - 1 )
    half_interval = 0.5_dp * delta_cond


    ! z-direction spatial homogeneization
    nxybins = (ex-sx+1) * (ey-sy+1) * nbins
    allocate ( wrk1 (nxybins) , &
               wrk2 (nxybins) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate pdf')

    index = 1
    do j = sy,ey ! loop for each point in y-direction
       do i = sx,ex ! loop for each point in x-direction
          do n = 1,nbins ! loop for each bin in the var space

             mytarget     = min_value + (n-1) * delta_cond
             wrk1 (index) = 0.0_dp

             do k = sz,ez ! loop for each point in z-direction

                error = abs ( var (i,j,k) - mytarget )
                if ( error < half_interval ) then
                   wrk1 (index) = wrk1 (index) + 1.0_dp
                end if

             end do

             index = index + 1

          end do
       end do
    end do

    if ( comm1dz /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk1 , wrk2 , nxybins , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dz , mpicode )
    else
       wrk2 (:) = wrk1 (:)
    end if

    index = 1
    do j = sy,ey
       do i = sx,ex
          do n = 1,nbins
             pdf_var (i,j,n) = pdf_var (i,j,n) + wrk2 (index)
             index           = index + 1
          end do
       end do
    end do


    if ( ifile == end_file ) then

       if ( type_scale == logarit ) then ! this has to be physical dimensions

          do j = sy,ey ! loop for each point in y-direction
             do i = sx,ex ! loop for each point in x-direction
                wrk = 0.0_dp

                do n = 1,nbins-1 ! loop for each bin in the var space
                   mytarget = 10.0_dp ** ( min_value + delta_cond * (n) )
                   mytarget = mytarget - 10.0_dp ** ( min_value + delta_cond * (n-1) )
                   wrk = wrk + pdf_var (i,j,n) * mytarget
                end do
                wrk = wrk + pdf_var (i,j,nbins) * mytarget

                wrk             = 1.0_dp / max ( wrk , epsi )
                pdf_var (i,j,:) = pdf_var (i,j,:) * wrk
             end do
          end do

       else

          do j = sy,ey ! loop for each point in y-direction
             do i = sx,ex ! loop for each point in x-direction
                wrk = 0.0_dp
                do n = 1,nbins ! loop for each bin in the var space
                   wrk = wrk + pdf_var (i,j,n)
                end do
                wrk             = wrk * delta_cond ! verified
                wrk             = 1.0_dp / max ( wrk , epsi )
                pdf_var (i,j,:) = pdf_var (i,j,:) * wrk
             end do
          end do

       end if

    end if


    deallocate ( wrk1 , wrk2 )


  end subroutine pdf


!> \brief Calculate the conditional average.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine pdf_conditional ( type_scale , ifile , end_file , min_value , max_value , cond_avg , var , pdf_var )


    character (len_default) , intent (in)                              :: type_scale !< pdf abscise scale
    integer (ip) , intent (in)                                         :: ifile      !< file number
    integer (ip) , intent (in)                                         :: end_file   !< last file number
    real (dp) , intent (in)                                            :: min_value  !< minimum value of the pdf interval
    real (dp) , intent (in)                                            :: max_value  !< minimum value of the pdf interval
    real (dp) , allocatable , dimension (:,:) , intent (in)            :: cond_avg   !< conditional averaged value
    real (dp) , dimension (:,:,:) , allocatable , intent (in)          :: var        !< variable
    real (dp) , dimension (:,:) , allocatable , intent (inout)         :: pdf_var    !< conditional pdf variable


    integer (ip)                                    :: ok , index , i , j , k , n
    integer (ip)                                    :: nxbins
    real (dp)                                       :: wrk
    real (dp)                                       :: delta_cond , half_interval , mytarget , error
    real (dp) , dimension (:) , allocatable         :: wrk1 , wrk2


    delta_cond    = ( max_value - min_value ) / ( nbins - 1 )
    half_interval = 0.5_dp * delta_cond


    ! yz-direction spatial homogeneization
    nxbins = (ex-sx+1) * nbins
    allocate ( wrk1 (nxbins) , &
               wrk2 (nxbins) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate pdf_axial')

    index = 1
    do i = sx,ex ! loop for each point in x-direction
       do n = 1,nbins ! loop for each bin in the var space

          mytarget     = min_value + (n-1) * delta_cond
          wrk1 (index) = 0.0_dp

          do k = sz,ez ! loop for each point in z-direction
             do j = sy,ey ! loop for each point in y-direction

                error = abs ( var (i,j,k) - mytarget )
                if ( error < half_interval ) then
                   wrk1 (index) = wrk1 (index) + 1.0_dp
                end if

             end do
          end do

          index = index + 1

       end do
    end do

    if ( comm2dyz /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk1 , wrk2 , nxbins , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm2dyz , mpicode )
    else
       wrk2 (:) = wrk1 (:)
    end if

    index = 1
    do i = sx,ex
       do n = 1,nbins
          pdf_var (i,n) = pdf_var (i,n) + wrk2 (index)
          index         = index + 1
       end do
    end do


    if ( ifile == end_file ) then

       do i = sx,ex
          do n = 1,nbins
             pdf_var (i,n) = pdf_var (i,n) * cond_avg (i,n)
          end do
       end do

       if ( type_scale == logarit ) then ! this has to be physical dimensions

          do i = sx,ex ! loop for each point in x-direction
             wrk = 0.0_dp

             do n = 1,nbins-1 ! loop for each bin in the var space
                mytarget = 10.0_dp ** ( min_value + delta_cond * (n) )
                mytarget = mytarget - 10.0_dp ** ( min_value + delta_cond * (n-1) )
                wrk = wrk + pdf_var (i,n) * mytarget
             end do
             wrk = wrk + pdf_var (i,nbins) * mytarget

             wrk           = 1.0_dp / max ( wrk , epsi )
             pdf_var (i,:) = pdf_var (i,:) * wrk
          end do

       else

          do i = sx,ex ! loop for each point in x-direction
             wrk = 0.0_dp
             do n = 1,nbins ! loop for each bin in the var space
                wrk = wrk + pdf_var (i,n)
             end do
             wrk           = wrk * delta_cond ! verified
             wrk           = 1.0_dp / max ( wrk , epsi )
             pdf_var (i,:) = pdf_var (i,:) * wrk
          end do

       end if

    end if


    deallocate ( wrk1 , wrk2 )


  end subroutine pdf_conditional


end module statistics
