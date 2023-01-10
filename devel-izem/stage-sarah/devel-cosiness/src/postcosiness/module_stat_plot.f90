!------------------------------------------------------------------------------
! MODULE: stat_plot
!------------------------------------------------------------------------------
!> \brief ParaView plot module for _statistical_ variables.
!!
!! This module allows to plot all types of statistical variables in a
!! ParaView file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module stat_plot


  use parameters
  use parallel
  use input
  use type_thd ! NOT USED ?!
  use adim
  use variables
  use deriv
!  use tools     ! NOT USED
  use tools_post
  use budget_tke
  use budget_Reynolds
  use budget_scal_var
  use autocorrelations

  implicit none


contains


!> \brief Plot a variable in a VTR file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine plot_stat_vtr ( vref , var )


    real (dp) , intent (in)                                            :: vref !< reference (physical) value to multiply with
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)       :: var  !< non-dimensional(in)/dimensional(out) variable to plot

    integer (ip) , parameter :: nbitsimple = 4

    integer (ip)             :: ix , fx , nx , &
                                iy , fy , ny , &
                                iz , fz , nz


    call domain_paraview ( ix , fx , iy , fy , iz , fz )
    nx = fx - ix + 1
    ny = fy - iy + 1
    nz = fz - iz + 1


    call comm_one (var)
    call comm_edges (var)


    write (unit_plot) nx * ny * nz * nbitsimple , &
                      real ( var (ix:fx,iy:fy,iz:fz) * vref , sp )


  end subroutine plot_stat_vtr


!> \brief Selector of instantaneous variables to plot.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine plot_stat_var ( idnumber , inp , thd , adi , grid , reyavg , favavg , stat )


    integer (ip) , intent (in)                                       :: idnumber !< number of the statistical variable
    type (inp_type) , intent (in)                                    :: inp      !< input derived type
    type (thd_type) , intent (in)                                    :: thd      !< thermodynamic derived type              ! NOT USED ?!
    type (adi_type) , intent (in)                                    :: adi      !< non-dimensional derived type
    type (inp_grid) , intent (in)                                    :: grid     !< grid derived type
    type (reyavg_type) , intent (inout)                              :: reyavg   !< Reynolds averge derived type
    type (favavg_type) , intent (inout)                              :: favavg   !< Favre average derived type
    type (stat_type) , intent (inout)                                :: stat     !< statistical derived type


    character (len_default) , parameter :: format = ' ( A11 , 1X , I2 , 1X , A10 , 1X , A ) '

    integer (ip)                        :: ok , spc , m , n , i ,j ,k

    real (dp)                           :: wrk

    character (len_default)             :: name

    ! temporal variables
    real (dp) , dimension (:,:,:) , allocatable     :: wrk1 , wrk2 , wrk3 , &
                                                       wrk4 , wrk5 , wrk6 , &
                                                       wrk7 , wrk8 , wrk9 , &
                                                       wrk10

    real (dp) , dimension (:,:,:,:) , allocatable   :: wrk_dim4

    real (dp) , dimension (:,:,:,:,:) , allocatable :: wrk_dim5


    name = inp % var_name (idnumber)

    if ( rank == rank_default ) &
       write ( * , format ) 'variable #' , idnumber , ', plotting' , trim (name)

    allocate  ( wrk1     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk2     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk3     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk4     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk5     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk6     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk7     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk8     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk9     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk10    ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )                     , &
                wrk_dim4 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax )           , &
                wrk_dim5 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , ndimmax , ndimmax ) , &
                stat = ok )
    if ( ok > 0 ) stop 'error allocate plot_stat_vars'


    ! basic variables


    if ( name == 'ra_Density' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % rho_ref , wrk1 )

    else if ( name == 'ra_Velocity-X' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % ux (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'rv_Velocity-X' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_ux (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'fa_Velocity-X' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'fv_Velocity-X' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_ux (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'ra_Velocity-Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % vy (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'rv_Velocity-Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_vy (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'fa_Velocity-Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'fv_Velocity-Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_vy (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'ra_Velocity-Z' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % wz (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'rv_Velocity-Z' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_wz (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'fa_Velocity-Z' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'fv_Velocity-Z' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_wz (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'ra_Total_Energy' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % et (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'fa_Total_Energy' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % et (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'ra_Enthalpy' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % H (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'fa_Enthalpy' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % H (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'ra_Total_Enthalpy' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ht (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'fa_Total_Enthalpy' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ht (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'ra_Temperature' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % T (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % T_ref , wrk1 )

    else if ( name == 'ra_max_Temperature' ) then

       call max_val_Y ( reyavg % T , wrk1 )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % T_ref , wrk1 )

    else if ( name == 'fa_Temperature' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % T (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % T_ref , wrk1 )

    else if ( name == 'ra_Total_Temp' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Tt (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % T_ref , wrk1 )

    else if ( name == 'fa_Total_Temp' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Tt (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % T_ref , wrk1 )

    else if ( name == 'ra_Pressure' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % p (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref , wrk1 )

    else if ( name == 'fa_Pressure' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % p (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref , wrk1  )

    else if ( name == 'ra_Molar_Mass' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % W (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % W_ref , wrk1 )

    else if ( name == 'fa_Molar_Mass' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % W (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % W_ref , wrk1 )

    else if ( name == 'ra_Sound_Speed' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % cs (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'fa_Sound_Speed' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % cs (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'ra_Cp' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % cp (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % cp_ref , wrk1 )

    else if ( name == 'fa_Cp' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % cp (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % cp_ref , wrk1 )


    ! viscous variables


    else if ( name == 'ra_Dyn_Viscosity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % mu (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % mu_ref , wrk1 )

    else if ( name == 'fa_Dyn_Viscosity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % mu (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % mu_ref , wrk1 )

    else if ( name == 'ra_SGS_Viscosity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % mu_SGS (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % mu_ref , wrk1 )

    else if ( name == 'fa_SGS_Viscosity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % mu_SGS (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % mu_ref , wrk1 )

    else if ( name == 'ra_Kin_Viscosity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % nu (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % mu_ref / adi % rho_ref , wrk1 )

    else if ( name == 'fa_Kin_Viscosity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % nu (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % mu_ref / adi % rho_ref , wrk1 )

    else if ( name == 'ra_Bulk' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % kpa (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % mu_ref , wrk1 )

    else if ( name == 'fa_Bulk' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % kpa (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % mu_ref , wrk1 )

    else if ( name == 'ra_Therm_Cond' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % ct (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % lbda_ref , wrk1 )

    else if ( name == 'fa_Therm_Cond' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % ct (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % lbda_ref , wrk1 )


    ! non dimensional variables


    else if ( name == 'ra_Gamma' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % gamma (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Gamma' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % gamma (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Mach' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % M (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Mach' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % M (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Prandtl' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Pr (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % Pr , wrk1 )

    else if ( name == 'fa_Prandtl' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Pr (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % Pr , wrk1 )



    ! reactive variables


    else if ( name == 'ra_Heat_Rel' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % omega (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref / adi % time_ref , wrk1 )

    else if ( name == 'ra_Heat_Rel_hct' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % omega_hct (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref / adi % time_ref , wrk1 )


    else if ( name == 'ra_max_Heat_Rel' ) then

       call max_val_Y ( reyavg % omega , wrk1 )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref / adi % time_ref , wrk1 )

    else if ( name == 'fa_Heat_Rel' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % omega (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref / adi % time_ref , wrk1 )

    else if ( name == 'fa_Heat_Rel_hct' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % omega_hct (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref / adi % time_ref , wrk1 )


    else if ( name == 'ra_Da' ) then

       wrk = ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Da (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'ra_Da_II' ) then

       wrk = ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Da_II (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'fa_Da' ) then

       wrk = ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Da (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'ra_sdr' ) then

       if ( inp % dis_rate_therm ) then
          wrk = ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / &
                ( adi % L_ref * adi % L_ref )
       else
          wrk = adi % D_ref / ( adi % L_ref * adi % L_ref )
       end if
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % sdr (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'fa_sdr' ) then

       if ( inp % dis_rate_therm ) then
          wrk = ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / &
                ( adi % L_ref * adi % L_ref )
       else
          wrk = adi % D_ref / ( adi % L_ref * adi % L_ref )
       end if
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % sdr (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'ra_FI' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % FI (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_FI' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % FI (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_proba_pz' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % proba_pz (sx:ex,sy:ey,sz:ez)
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Reactivity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % lambda (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_Reactivity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % lambda (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )


    ! statistical variables


    else if ( name == 'domega' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       call f_domega ( inp , adi , grid % dy_i , wrk1 , wrk2 )
       call plot_stat_vtr ( adi % L_ref , wrk2 )

    else if ( name == 'dtheta' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       call f_dtheta ( inp , adi , grid % y , wrk1 , wrk2 , wrk3 )
       call plot_stat_vtr ( adi % L_ref , wrk3 )

    else if ( name == 'dmix' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       if ( inp % index_scalar == 0 ) then
          wrk2 (sx:ex,sy:ey,sz:ez) = favavg % Zm (sx:ex,sy:ey,sz:ez)
       else
          wrk2 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,inp % index_scalar)
       end if
       call f_dmix ( inp , adi , grid % y , wrk1 , wrk2 , wrk3 )
       call plot_stat_vtr ( adi % L_ref , wrk3 )

    else if ( name == 'dvis99' ) then

       if ( inp % index_scalar == 0 ) then
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Zm (sx:ex,sy:ey,sz:ez)
       else
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,inp % index_scalar)
       end if
       call dvis_99 ( inp , adi , grid % y , wrk1 , wrk2 )
       call plot_stat_vtr ( adi % L_ref , wrk2 )

    else if ( name == 'TKE' .or. name == 'K' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'fa_TKE_SGS' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % tke_SGS (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'ra_tau_iso_SGS' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % tau_iso_SGS (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref , wrk1 )

    else if ( name == 'R_11' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_ux (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'R_22' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_vy (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'R_33' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_wz (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'R_12' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxvy (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'R_13' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxwz (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'R_23' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_vywz (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'ra_tau_SGS_11' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % tau_SGS (sx:ex,sy:ey,sz:ez,1,1)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref , wrk1 )

    else if ( name == 'ra_tau_SGS_22' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % tau_SGS (sx:ex,sy:ey,sz:ez,2,2)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref , wrk1 )

    else if ( name == 'ra_tau_SGS_33' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % tau_SGS (sx:ex,sy:ey,sz:ez,3,3)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref , wrk1 )

    else if ( name == 'ra_tau_SGS_12' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % tau_SGS (sx:ex,sy:ey,sz:ez,1,2)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref , wrk1 )

    else if ( name == 'ra_tau_SGS_13' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % tau_SGS (sx:ex,sy:ey,sz:ez,1,3)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref , wrk1 )

    else if ( name == 'ra_tau_SGS_23' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % tau_SGS (sx:ex,sy:ey,sz:ez,2,3)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref , wrk1 )

    else if ( name == 'b_11' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * stat % favvar_ux (sx:ex,sy:ey,sz:ez) / &
                                  wrk1 (sx:ex,sy:ey,sz:ez) - 1.0_dp / 3.0_dp
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'b_22' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * stat % favvar_vy (sx:ex,sy:ey,sz:ez) / &
                                  wrk1 (sx:ex,sy:ey,sz:ez) - 1.0_dp / 3.0_dp
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'b_33' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * stat % favvar_wz (sx:ex,sy:ey,sz:ez) / &
                                  wrk1 (sx:ex,sy:ey,sz:ez) - 1.0_dp / 3.0_dp
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'b_12' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * stat % favvar_uxvy (sx:ex,sy:ey,sz:ez) / &
                                  wrk1 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'b_13' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * stat % favvar_uxwz (sx:ex,sy:ey,sz:ez) / &
                                  wrk1 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'b_23' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * stat % favvar_vywz (sx:ex,sy:ey,sz:ez) / &
                                  wrk1 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Lenght_Taylor' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % nu (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = stat % tke_dissip (sx:ex,sy:ey,sz:ez)
       call length_taylor ( wrk1 , wrk2 , wrk3 , wrk4 )
       call similarity_x ( inp , wrk4 )
       call plot_stat_vtr ( adi % L_ref , wrk4 ) ! dimensions verified

    else if ( name == 'Reynolds_Taylor' ) then

       wrk = adi % rho_ref * adi % L_ref * adi % u_ref / adi % mu_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % nu (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = stat % tke_dissip (sx:ex,sy:ey,sz:ez)
       call reynolds_taylor ( wrk1 , wrk2 , wrk3 , wrk4 )
       call similarity_x ( inp , wrk4 )
       call plot_stat_vtr ( wrk , wrk4 )

    else if ( name == 'Lenght_Kolmogorov' ) then

!       wrk = sqrt ( adi % mu_ref * adi % L_ref / ( adi % rho_ref * adi % u_ref ) ) ! dimensions verified
       wrk = adi % L_ref * sqrt ( adi % sqgmr ) ! A.Techer
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % nu (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = stat % tke_dissip (sx:ex,sy:ey,sz:ez)
       call length_kolmogorov ( wrk1 , wrk2 , wrk3 )
       call similarity_x ( inp , wrk3 )
       call plot_stat_vtr ( wrk , wrk3 )

    else if ( name == 'Reynolds' ) then

       wrk = adi % rho_ref * adi % L_ref * adi % u_ref / adi % mu_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       call reynolds ( inp , adi , grid % dy_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( wrk , wrk2 )

    else if ( name == 'Index_quality_Celik' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % nu (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = stat % tke_dissip (sx:ex,sy:ey,sz:ez)
       call length_kolmogorov ( wrk1 , wrk2 , wrk3 )
       call mesh_quality_celik ( adi , grid % delta , wrk3 , wrk5 )
       call plot_stat_vtr ( 1.0_dp , wrk5 )

    else if ( name == 'Index_quality_Pope' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk2 (sx:ex,sy:ey,sz:ez) = favavg % tke_SGS (sx:ex,sy:ey,sz:ez)
       call mesh_quality_pope ( adi , wrk1 , wrk2 , wrk3 )
       call plot_stat_vtr ( 1.0_dp , wrk3 )


    ! other variables


    else if ( name == 'ra_Vorticity' ) then

       call vorticity ( grid % dx_i , grid % dy_i , grid % dz_i , reyavg % ux , reyavg % vy , reyavg % wz , wrk1 )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref / adi % L_ref , wrk1 )

    else if ( name == 'fa_Vorticity' ) then

       call vorticity ( grid % dx_i , grid % dy_i , grid % dz_i , favavg % ux , favavg % vy , favavg % wz , wrk1 )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref / adi % L_ref , wrk1 )

    else if ( name == 'ra_Enstrophy' ) then

       wrk = adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref )
       call vorticity ( grid % dx_i , grid % dy_i , grid % dz_i , reyavg % ux , reyavg % vy , reyavg % wz , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * wrk1 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'fa_Enstrophy' ) then

       wrk = adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref )
       call vorticity ( grid % dx_i , grid % dy_i , grid % dz_i , favavg % ux , favavg % vy , favavg % wz , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * wrk1 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'du1_x2' ) then

       call dy ( grid % dy_i , favavg % ux , wrk1 )
       call comm_one (wrk1)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'rho_du1_x2' ) then

       call dy ( grid % dy_i , favavg % ux , wrk1 )
       call comm_one (wrk1)
       call similarity_x ( inp , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) * reyavg % rho (sx:ex,sy:ey,sz:ez)
       call plot_stat_vtr ( adi % rho_ref / adi % time_ref , wrk1 )

    else if ( name == 'Mt' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                    stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                    stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk1 (sx:ex,sy:ey,sz:ez) = sqrt ( wrk1 (sx:ex,sy:ey,sz:ez) ) / &
                                  reyavg % cs (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Mg' ) then

       call comm_one (favavg % ux) ; call dy ( grid % dy_i , favavg % ux , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) / &
                                  reyavg % cs (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % L_ref , wrk1 ) ! integral length scale must be choosen

    else if ( name == 'Filter_width' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = grid % delta (sx:ex,sy:ey,sz:ez)
       call plot_stat_vtr ( adi % L_ref , wrk1 )

    else if ( name == 'dx_plus' ) then

       do k = sz , ez
          do i = sx , ex
             wrk1 (i,sy,k) = grid % XZ_dxplus (i,k,1,1)
             wrk1 (i,ey,k) = grid % XZ_dxplus (i,k,2,1)
          end do
       end do
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'dy_plus' ) then

       do k = sz , ez
          do i = sx , ex
             wrk1 (i,sy,k) = grid % XZ_dxplus (i,k,1,2)
             wrk1 (i,ey,k) = grid % XZ_dxplus (i,k,2,2)
          end do
       end do
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'dz_plus' ) then

       do k = sz , ez
          do i = sx , ex
             wrk1 (i,sy,k) = grid % XZ_dxplus (i,k,1,3)
             wrk1 (i,ey,k) = grid % XZ_dxplus (i,k,2,3)
          end do
       end do
       call plot_stat_vtr ( 1.0_dp , wrk1 )


    ! variances


    else if ( name == 'rv_rho' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_rho (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % rho_ref * adi % rho_ref , wrk1 )

    else if ( name == 'rv_rho_acu' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_rho_acu (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % rho_ref * adi % rho_ref , wrk1 )

    else if ( name == 'rv_rho_ent' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_rho_ent (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % rho_ref * adi % rho_ref , wrk1 )

    else if ( name == 'rv_p' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_p (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % p_ref * adi % p_ref , wrk1 )

    else if ( name == 'rv_T' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_T (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % T_ref * adi % T_ref , wrk1 )

    else if ( name == 'rv_T_acu' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_T_acu (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % T_ref * adi % T_ref , wrk1 )

    else if ( name == 'rv_T_ent' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_T_ent (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % T_ref * adi % T_ref , wrk1 )

    else if ( name == 'rv_W' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_W (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % W_ref * adi % W_ref , wrk1 )

    else if ( name == 'fv_du1_dx1' ) then

       wrk = adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_du1_dx1 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'rv_p_du1_dx1' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_p_du1_dx1 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'rv_rho_p' ) then

       wrk = adi % rho_ref * adi % p_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_rho_p (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'rv_rho_W' ) then

       wrk = adi % rho_ref * adi % W_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_rho_W (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'CORR_rho_W' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_rho_W (sx:ex,sy:ey,sz:ez) / &
                           sqrt ( stat % reyvar_rho (sx:ex,sy:ey,sz:ez) *   &
                                  stat % reyvar_W (sx:ex,sy:ey,sz:ez) )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'CORR_rho_p' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_rho_p (sx:ex,sy:ey,sz:ez) / &
                           sqrt ( stat % reyvar_rho (sx:ex,sy:ey,sz:ez) *   &
                                  stat % reyvar_p (sx:ex,sy:ey,sz:ez) )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'CORR_p_du1_dx1' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_p_du1_dx1 (sx:ex,sy:ey,sz:ez) / &
                           sqrt ( stat % reyvar_p (sx:ex,sy:ey,sz:ez) *         &
                                  stat % favvar_du1_dx1 (sx:ex,sy:ey,sz:ez) )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Heat_Rel' ) then

       wrk = adi % u_ref * adi % u_ref / adi % time_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_omega (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk * wrk , wrk1 )

    else if ( name == 'fv_Heat_Rel' ) then

       wrk = adi % u_ref * adi % u_ref / adi % time_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_omega (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk * wrk , wrk1 )

    else if ( name == 'rv_Da' ) then

       wrk = ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Da (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk * wrk , wrk1 )

    else if ( name == 'fv_Da' ) then

       wrk = ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Da (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk * wrk , wrk1 )

    else if ( name == 'rv_sdr' ) then

       if ( inp % dis_rate_therm ) then
          wrk = ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / &
                ( adi % L_ref * adi % L_ref )
       else
          wrk = adi % D_ref / ( adi % L_ref * adi % L_ref )
       end if
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_sdr (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk * wrk , wrk1 )

    else if ( name == 'fv_sdr' ) then

       if ( inp % dis_rate_therm ) then
          wrk = ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / &
                ( adi % L_ref * adi % L_ref )
       else
          wrk = adi % D_ref / ( adi % L_ref * adi % L_ref )
       end if
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_sdr (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk * wrk , wrk1 )

    else if ( name == 'rv_FI' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_FI (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_FI' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_FI (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )


    ! scalar variables


    else if ( name == 'Scalar_Flux_ux_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favavg_ux_Y (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Scalar_Flux_ux_Zm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favavg_ux_Zm (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Scalar_Flux_vy_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favavg_vy_Y (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Scalar_Flux_vy_Zm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favavg_vy_Zm (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Scalar_Flux_wz_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favavg_wz_Y (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Scalar_Flux_wz_Zm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favavg_wz_Zm (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Grad_x_Zm' ) then

       if ( npv == 0 ) then ! calculated passive scalar
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Zm (sx:ex,sy:ey,sz:ez)
       else                 ! real passive scalar
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez, nrv+npv )
       end if
       call comm_one (wrk1)
       call dx ( grid % dx_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp / adi % L_ref , wrk2 )

    else if ( name == 'Grad_x_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez, max ( inp % index_scalar , 1 ) )
       call comm_one (wrk1)
       call dx ( grid % dx_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp / adi % L_ref , wrk2 )

    else if ( name == 'Grad_y_Zm' ) then

       if ( npv == 0 ) then ! calculated passive scalar
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Zm (sx:ex,sy:ey,sz:ez)
       else                 ! real passive scalar
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez, nrv+npv )
       end if
       call comm_one (wrk1)
       call dy ( grid % dy_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp / adi % L_ref , wrk2 )

    else if ( name == 'Grad_y_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez, max ( inp % index_scalar , 1 ) )
       call comm_one (wrk1)
       call dy ( grid % dy_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp / adi % L_ref , wrk2 )

    else if ( name == 'Grad_z_Zm' ) then

       if ( npv == 0 ) then ! calculated passive scalar
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Zm (sx:ex,sy:ey,sz:ez)
       else                 ! real passive scalar
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez, nrv+npv )
       end if
       call comm_one (wrk1)
       call dz ( grid % dz_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp / adi % L_ref , wrk2 )

    else if ( name == 'Grad_z_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez, max ( inp % index_scalar , 1 ) )
       call comm_one (wrk1)
       call dz ( grid % dz_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp / adi % L_ref , wrk2 )


    ! budget equations


    else if ( name == 'MASS_1' ) then

       wrk = adi % rho_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                  favavg % ux (sx:ex,sy:ey,sz:ez)
       call comm_one (wrk1)
       call dx ( grid % dx_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( wrk , wrk2 )

    else if ( name == 'MASS_2' ) then

       wrk = adi % rho_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                  favavg % vy (sx:ex,sy:ey,sz:ez)
       call comm_one (wrk1)
       call dy ( grid % dy_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( wrk , wrk2 )

    else if ( name == 'MASS_3' ) then

       wrk = adi % rho_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                  favavg % wz (sx:ex,sy:ey,sz:ez)
       call comm_one (wrk1)
       call dz ( grid % dz_i , wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( wrk , wrk2 )



    else if ( name == 'TKE_Dissipation' ) then

       wrk = adi % mu_ref * adi % u_ref * adi % u_ref / &
           ( adi % rho_ref * adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tke_dissip (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'minus_rho_TKE_Dissipation' ) then

       wrk = adi % mu_ref * adi % u_ref * adi % u_ref / &
           ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                    stat % tke_dissip (sx:ex,sy:ey,sz:ez) ! dissipation is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Dissipation_11' ) then

       wrk = adi % mu_ref * adi % u_ref * adi % u_ref / &
           ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                    stat % rey_dissip_11 (sx:ex,sy:ey,sz:ez) ! dissipation is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Dissipation_22' ) then

       wrk = adi % mu_ref * adi % u_ref * adi % u_ref / &
           ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                    stat % rey_dissip_22 (sx:ex,sy:ey,sz:ez) ! dissipation is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Dissipation_33' ) then

       wrk = adi % mu_ref * adi % u_ref * adi % u_ref / &
           ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                    stat % rey_dissip_33 (sx:ex,sy:ey,sz:ez) ! dissipation is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Dissipation_12' ) then

       wrk = adi % mu_ref * adi % u_ref * adi % u_ref / &
           ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                    stat % rey_dissip_12 (sx:ex,sy:ey,sz:ez) ! dissipation is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Dissipation_13' ) then

       wrk = adi % mu_ref * adi % u_ref * adi % u_ref / &
           ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                    stat % rey_dissip_13 (sx:ex,sy:ey,sz:ez) ! dissipation is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Dissipation_23' ) then

       wrk = adi % mu_ref * adi % u_ref * adi % u_ref / &
           ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - reyavg % rho (sx:ex,sy:ey,sz:ez) * &
                                    stat % rey_dissip_23 (sx:ex,sy:ey,sz:ez) ! dissipation is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'TKE_Pressure_Strain' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tke_press_strain (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Pressure_Strain_11' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % rey_press_strain_11 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Pressure_Strain_22' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % rey_press_strain_22 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Pressure_Strain_33' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % rey_press_strain_33 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Pressure_Strain_12' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % rey_press_strain_12 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Pressure_Strain_13' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % rey_press_strain_13 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Pressure_Strain_23' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % rey_press_strain_23 (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'TKE_Transport' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       call tke_transport2 ( grid % dx_i , grid % dy_i , grid % dz_i , stat % tke_transp , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = - wrk1 (sx:ex,sy:ey,sz:ez) ! transport is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Transport_11' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       call rey_transport2 ( grid % dx_i , grid % dy_i , grid % dz_i , stat % rey_transp_11 , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = - wrk1 (sx:ex,sy:ey,sz:ez) ! transport is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Transport_22' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       call rey_transport2 ( grid % dx_i , grid % dy_i , grid % dz_i , stat % rey_transp_22 , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = - wrk1 (sx:ex,sy:ey,sz:ez) ! transport is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Transport_33' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       call rey_transport2 ( grid % dx_i , grid % dy_i , grid % dz_i , stat % rey_transp_33 , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = - wrk1 (sx:ex,sy:ey,sz:ez) ! transport is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Transport_12' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       call rey_transport2 ( grid % dx_i , grid % dy_i , grid % dz_i , stat % rey_transp_12 , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = - wrk1 (sx:ex,sy:ey,sz:ez) ! transport is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Transport_13' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       call rey_transport2 ( grid % dx_i , grid % dy_i , grid % dz_i , stat % rey_transp_13 , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = - wrk1 (sx:ex,sy:ey,sz:ez) ! transport is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'REY_Transport_23' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       call rey_transport2 ( grid % dx_i , grid % dy_i , grid % dz_i , stat % rey_transp_23 , wrk1)
       wrk1 (sx:ex,sy:ey,sz:ez) = - wrk1 (sx:ex,sy:ey,sz:ez) ! transport is a negative term
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'TKE_Mass_Flux_Coupl' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_ux (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_vy (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_wz (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = reyavg % P (sx:ex,sy:ey,sz:ez)
       do m = 1 , ndimmax
          do n = 1 , ndimmax
             wrk5 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,m,n)
             wrk_dim5 (sx:ex,sy:ey,sz:ez,m,n) = wrk5 (sx:ex,sy:ey,sz:ez)
          end do
       end do
       call tke_mass_flux_coupl ( adi , grid % dx_i , grid % dy_i , grid % dz_i ,  &
                                  wrk1 , wrk2 , wrk3 , wrk4 , &
                                  wrk_dim5 , wrk6 )
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Mass_Flux_Coupl_11' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_ux (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = reyavg % P (sx:ex,sy:ey,sz:ez)
       do m = 1 , ndimmax
          do n = 1 , ndimmax
             wrk5 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,m,n)
             wrk_dim5 (sx:ex,sy:ey,sz:ez,m,n) = wrk5 (sx:ex,sy:ey,sz:ez)
          end do
       end do
       call rey_mass_flux_coupl_11 ( adi , grid % dx_i , grid % dy_i , grid % dz_i , &
                                     wrk1 , wrk4 ,              &
                                     wrk_dim5 , wrk6 )
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Mass_Flux_Coupl_22' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_vy (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = reyavg % P (sx:ex,sy:ey,sz:ez)
       do m = 1 , ndimmax
          do n = 1 , ndimmax
             wrk5 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,m,n)
             wrk_dim5 (sx:ex,sy:ey,sz:ez,m,n) = wrk5 (sx:ex,sy:ey,sz:ez)
          end do
       end do
       call rey_mass_flux_coupl_22 ( adi , grid % dx_i , grid % dy_i , grid % dz_i , &
                                     wrk1 , wrk4 ,              &
                                     wrk_dim5 , wrk6 )
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Mass_Flux_Coupl_33' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_wz (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = reyavg % P (sx:ex,sy:ey,sz:ez)
       do m = 1 , ndimmax
          do n = 1 , ndimmax
             wrk5 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,m,n)
             wrk_dim5 (sx:ex,sy:ey,sz:ez,m,n) = wrk5 (sx:ex,sy:ey,sz:ez)
          end do
       end do
       call rey_mass_flux_coupl_33 ( adi , grid % dx_i , grid % dy_i , grid % dz_i , &
                                     wrk1 , wrk4 ,              &
                                     wrk_dim5 , wrk6 )
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Mass_Flux_Coupl_12' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_ux (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_vy (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = reyavg % P (sx:ex,sy:ey,sz:ez)
       do m = 1 , ndimmax
          do n = 1 , ndimmax
             wrk5 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,m,n)
             wrk_dim5 (sx:ex,sy:ey,sz:ez,m,n) = wrk5 (sx:ex,sy:ey,sz:ez)
          end do
       end do
       call rey_mass_flux_coupl_12 ( adi , grid % dx_i , grid % dy_i , grid % dz_i , &
                                     wrk1 , wrk2, wrk4 ,        &
                                     wrk_dim5 , wrk6 )
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Mass_Flux_Coupl_13' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_ux (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_wz (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = reyavg % P (sx:ex,sy:ey,sz:ez)
       do m = 1 , ndimmax
          do n = 1 , ndimmax
             wrk5 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,m,n)
             wrk_dim5 (sx:ex,sy:ey,sz:ez,m,n) = wrk5 (sx:ex,sy:ey,sz:ez)
          end do
       end do
       call rey_mass_flux_coupl_12 ( adi , grid % dx_i , grid % dy_i , grid % dz_i , &
                                     wrk1 , wrk2, wrk4 ,        &
                                     wrk_dim5 , wrk6 )
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Mass_Flux_Coupl_23' ) then

       wrk = adi % p_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_vy (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = stat % tke_ra_ff_wz (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = reyavg % P (sx:ex,sy:ey,sz:ez)
       do m = 1 , ndimmax
          do n = 1 , ndimmax
             wrk5 (sx:ex,sy:ey,sz:ez) = reyavg % tau (sx:ex,sy:ey,sz:ez,m,n)
             wrk_dim5 (sx:ex,sy:ey,sz:ez,m,n) = wrk5 (sx:ex,sy:ey,sz:ez)
          end do
       end do
       call rey_mass_flux_coupl_12 ( adi , grid % dx_i , grid % dy_i , grid % dz_i , &
                                     wrk1 , wrk2, wrk4 ,        &
                                     wrk_dim5 , wrk6 )
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'TKE_Production' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = stat % favvar_ux (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = stat % favvar_vy (sx:ex,sy:ey,sz:ez)
       wrk6 (sx:ex,sy:ey,sz:ez) = stat % favvar_wz (sx:ex,sy:ey,sz:ez)
       wrk7 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxvy (sx:ex,sy:ey,sz:ez)
       wrk8 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxwz (sx:ex,sy:ey,sz:ez)
       wrk9 (sx:ex,sy:ey,sz:ez) = stat % favvar_vywz (sx:ex,sy:ey,sz:ez)
       call tke_production ( grid % dx_i , grid % dy_i , grid % dz_i , &
                             wrk1 , wrk2 , wrk3 , &
                             wrk4 , wrk5 , wrk6 , &
                             wrk7 , wrk8 , wrk9 , &
                             wrk10 )
       wrk10 (sx:ex,sy:ey,sz:ez) = wrk10 (sx:ex,sy:ey,sz:ez) * &
                                   reyavg % rho (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk10 )
       call plot_stat_vtr ( wrk , wrk10 )

    else if ( name == 'REY_Production_11' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = stat % favvar_ux (sx:ex,sy:ey,sz:ez)
       wrk7 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxvy (sx:ex,sy:ey,sz:ez)
       wrk8 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxwz (sx:ex,sy:ey,sz:ez)
       call rey_production_11 ( grid % dx_i , grid % dy_i , grid % dz_i , &
                                wrk1 ,               &
                                wrk4 ,               &
                                wrk7 , wrk8 ,        &
                                wrk10 )
       wrk10 (sx:ex,sy:ey,sz:ez) = wrk10 (sx:ex,sy:ey,sz:ez) * &
                                   reyavg % rho (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk10 )
       call plot_stat_vtr ( wrk , wrk10 )

    else if ( name == 'REY_Production_22' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = stat % favvar_vy (sx:ex,sy:ey,sz:ez)
       wrk7 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxvy (sx:ex,sy:ey,sz:ez)
       wrk8 (sx:ex,sy:ey,sz:ez) = stat % favvar_vywz (sx:ex,sy:ey,sz:ez)
       call rey_production_11 ( grid % dx_i , grid % dy_i , grid % dz_i , &
                                wrk1 ,               &
                                wrk4 ,               &
                                wrk7 , wrk8 ,        &
                                wrk10 )
       wrk10 (sx:ex,sy:ey,sz:ez) = wrk10 (sx:ex,sy:ey,sz:ez) * &
                                   reyavg % rho (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk10 )
       call plot_stat_vtr ( wrk , wrk10 )

    else if ( name == 'REY_Production_33' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = stat % favvar_wz (sx:ex,sy:ey,sz:ez)
       wrk7 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxwz (sx:ex,sy:ey,sz:ez)
       wrk8 (sx:ex,sy:ey,sz:ez) = stat % favvar_vywz (sx:ex,sy:ey,sz:ez)
       call rey_production_11 ( grid % dx_i , grid % dy_i , grid % dz_i , &
                                wrk1 ,               &
                                wrk4 ,               &
                                wrk7 , wrk8 ,        &
                                wrk10 )
       wrk10 (sx:ex,sy:ey,sz:ez) = wrk10 (sx:ex,sy:ey,sz:ez) * &
                                   reyavg % rho (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk10 )
       call plot_stat_vtr ( wrk , wrk10 )

    else if ( name == 'REY_Production_12' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = stat % favvar_ux (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = stat % favvar_vy (sx:ex,sy:ey,sz:ez)
       wrk7 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxvy (sx:ex,sy:ey,sz:ez)
       wrk8 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxwz (sx:ex,sy:ey,sz:ez)
       wrk9 (sx:ex,sy:ey,sz:ez) = stat % favvar_vywz (sx:ex,sy:ey,sz:ez)
       call rey_production_12 ( grid % dx_i , grid % dy_i , grid % dz_i , &
                                wrk1 , wrk2 ,        &
                                wrk4 , wrk5 ,        &
                                wrk7 , wrk8 , wrk9 , &
                                wrk10 )
       wrk10 (sx:ex,sy:ey,sz:ez) = wrk10 (sx:ex,sy:ey,sz:ez) * &
                                   reyavg % rho (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk10 )
       call plot_stat_vtr ( wrk , wrk10 )

    else if ( name == 'REY_Production_13' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = stat % favvar_ux (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = stat % favvar_wz (sx:ex,sy:ey,sz:ez)
       wrk7 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxvy (sx:ex,sy:ey,sz:ez)
       wrk8 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxwz (sx:ex,sy:ey,sz:ez)
       wrk9 (sx:ex,sy:ey,sz:ez) = stat % favvar_vywz (sx:ex,sy:ey,sz:ez)
       call rey_production_13 ( grid % dx_i , grid % dy_i , grid % dz_i , &
                                wrk1 , wrk2 ,        &
                                wrk4 , wrk5 ,        &
                                wrk7 , wrk8 , wrk9 , &
                                wrk10 )
       wrk10 (sx:ex,sy:ey,sz:ez) = wrk10 (sx:ex,sy:ey,sz:ez) * &
                                   reyavg % rho (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk10 )
       call plot_stat_vtr ( wrk , wrk10 )

    else if ( name == 'REY_Production_23' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = stat % favvar_vy (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = stat % favvar_wz (sx:ex,sy:ey,sz:ez)
       wrk7 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxvy (sx:ex,sy:ey,sz:ez)
       wrk8 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxwz (sx:ex,sy:ey,sz:ez)
       wrk9 (sx:ex,sy:ey,sz:ez) = stat % favvar_vywz (sx:ex,sy:ey,sz:ez)
       call rey_production_23 ( grid % dx_i , grid % dy_i , grid % dz_i , &
                                wrk1 , wrk2 ,        &
                                wrk4 , wrk5 ,        &
                                wrk7 , wrk8 , wrk9 , &
                                wrk10 )
       wrk10 (sx:ex,sy:ey,sz:ez) = wrk10 (sx:ex,sy:ey,sz:ez) * &
                                   reyavg % rho (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk10 )
       call plot_stat_vtr ( wrk , wrk10 )

    else if ( name == 'TKE_Convection' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = 0.5_dp * ( stat % favvar_ux (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_vy (sx:ex,sy:ey,sz:ez) + &
                                             stat % favvar_wz (sx:ex,sy:ey,sz:ez) )
       wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       call tke_convection ( grid % dx_i , grid % dy_i , grid % dz_i , &
                             wrk2 , wrk3 , wrk4 , &
                             wrk5 , wrk1 , wrk6 )
       wrk6 (sx:ex,sy:ey,sz:ez) = - wrk6 (sx:ex,sy:ey,sz:ez) ! LHS to RHS
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Convection_11' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_ux (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       call rey_convection ( grid % dx_i , grid % dy_i , grid % dz_i , &
                             wrk2 , wrk3 , wrk4 , &
                             wrk5 , wrk1 , wrk6 )
       wrk6 (sx:ex,sy:ey,sz:ez) = - wrk6 (sx:ex,sy:ey,sz:ez) ! LHS to RHS
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Convection_22' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_vy (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       call rey_convection ( grid % dx_i , grid % dy_i , grid % dz_i , &
                             wrk2 , wrk3 , wrk4 , &
                             wrk5 , wrk1 , wrk6 )
       wrk6 (sx:ex,sy:ey,sz:ez) = - wrk6 (sx:ex,sy:ey,sz:ez) ! LHS to RHS
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Convection_33' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_wz (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       call rey_convection ( grid % dx_i , grid % dy_i , grid % dz_i , &
                             wrk2 , wrk3 , wrk4 , &
                             wrk5 , wrk1 , wrk6 )
       wrk6 (sx:ex,sy:ey,sz:ez) = - wrk6 (sx:ex,sy:ey,sz:ez) ! LHS to RHS
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Convection_12' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxvy (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       call rey_convection ( grid % dx_i , grid % dy_i , grid % dz_i , &
                             wrk2 , wrk3 , wrk4 , &
                             wrk5 , wrk1 , wrk6 )
       wrk6 (sx:ex,sy:ey,sz:ez) = - wrk6 (sx:ex,sy:ey,sz:ez) ! LHS to RHS
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Convection_13' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_uxwz (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       call rey_convection ( grid % dx_i , grid % dy_i , grid % dz_i , &
                             wrk2 , wrk3 , wrk4 , &
                             wrk5 , wrk1 , wrk6 )
       wrk6 (sx:ex,sy:ey,sz:ez) = - wrk6 (sx:ex,sy:ey,sz:ez) ! LHS to RHS
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )

    else if ( name == 'REY_Convection_23' ) then

       wrk = adi % rho_ref * adi % u_ref * adi % u_ref * adi % u_ref / adi % L_ref
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_vywz (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = reyavg % rho (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = favavg % ux (sx:ex,sy:ey,sz:ez)
       wrk4 (sx:ex,sy:ey,sz:ez) = favavg % vy (sx:ex,sy:ey,sz:ez)
       wrk5 (sx:ex,sy:ey,sz:ez) = favavg % wz (sx:ex,sy:ey,sz:ez)
       call rey_convection ( grid % dx_i , grid % dy_i , grid % dz_i , &
                             wrk2 , wrk3 , wrk4 , &
                             wrk5 , wrk1 , wrk6 )
       wrk6 (sx:ex,sy:ey,sz:ez) = - wrk6 (sx:ex,sy:ey,sz:ez) ! LHS to RHS
       call similarity_x ( inp , wrk6 )
       call plot_stat_vtr ( wrk , wrk6 )


    else if ( name == 'VOR_Convection' ) then

       wrk = adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - stat % vor_conv (sx:ex,sy:ey,sz:ez) ! LHS to RHS
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'ENS_Convection' ) then

       wrk = adi % u_ref * adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - stat % ens_conv (sx:ex,sy:ey,sz:ez) ! LHS to RHS
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'VOR_Stretching' ) then

       wrk = adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % vor_stret (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'ENS_Stretching' ) then

       wrk = adi % u_ref * adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % ens_stret (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'VOR_Dilatation' ) then

       wrk = adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - stat % vor_dila (sx:ex,sy:ey,sz:ez) ! this term is negative
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'ENS_Dilatation' ) then

       wrk = adi % u_ref * adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = - stat % ens_dila (sx:ex,sy:ey,sz:ez) ! this term is negative
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'VOR_Baroclinic' ) then

       wrk = adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % vor_baro (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'ENS_Baroclinic' ) then

       wrk = adi % u_ref * adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % ens_baro (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'VOR_Viscous' ) then

       wrk = adi % u_ref * adi % mu_ref / ( adi % rho_ref * adi % L_ref * adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % vor_visc (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'ENS_Viscous' ) then

       wrk = adi % u_ref * adi % u_ref * adi % mu_ref / &
           ( adi % rho_ref * adi % L_ref * adi % L_ref * adi % L_ref * adi % L_ref )
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % ens_visc (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )



    else if ( name == 'SCAL_VAR_Mean_Convection' ) then

       wrk = adi % rho_ref * adi % u_ref / adi % L_ref

       if ( inp % index_scalar == 0 ) then
          wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Zm (sx:ex,sy:ey,sz:ez) * &
                                     reyavg % rho (sx:ex,sy:ey,sz:ez)
       else
          wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,inp % index_scalar) * &
                                     reyavg % rho (sx:ex,sy:ey,sz:ez)
       end if

       call scal_var_mean_conv ( grid % dx_i , grid % dy_i ,grid %  dz_i , &
                                 favavg % ux , favavg % vy , favavg % wz , &
                                 wrk1 , wrk2 )
       wrk2 (sx:ex,sy:ey,sz:ez) = - wrk2 (sx:ex,sy:ey,sz:ez) ! LHS to RHS

       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( wrk , wrk2 )

    else if ( name == 'SCAL_VAR_Turb_Transport' ) then

       wrk = adi % rho_ref * adi % u_ref / adi % L_ref
       call scal_var_turb_transp_last ( grid % dx_i , grid % dy_i , grid % dz_i , stat % scal_var_turb_transp , wrk1 )
       wrk1 (sx:ex,sy:ey,sz:ez) = - wrk1 (sx:ex,sy:ey,sz:ez) ! negative term

       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'SCAL_VAR_Turb_Production' ) then

       wrk = adi % rho_ref * adi % u_ref / adi % L_ref

       if ( inp % index_scalar == 0 ) then
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Zm (sx:ex,sy:ey,sz:ez)
       else
          wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,inp % index_scalar)
       end if

       call scal_var_turb_prod_last ( grid % dx_i , grid % dy_i , grid % dz_i , stat % scal_var_turb_prod , wrk1 , wrk2 )
       wrk2 (sx:ex,sy:ey,sz:ez) = - wrk2 (sx:ex,sy:ey,sz:ez) ! negative term

       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( wrk , wrk2 )

    else if ( name == 'SCAL_VAR_Mol_Diffusion' ) then

       if ( inp % dis_rate_therm ) then
          wrk = adi % lbda_ref / ( adi % rho_ref * adi % cp_ref )
       else
          wrk = adi % D_ref
       end if
       wrk = adi % rho_ref * wrk / ( adi % L_ref * adi % L_ref )

       call scal_var_mol_diff_last ( grid % dx_i , grid % dy_i , grid % dz_i , stat % scal_var_mol_diff , wrk1 ) ! this variable has implicit sign

       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'SCAL_VAR_Dissipation' ) then

       if ( inp % dis_rate_therm ) then
          wrk = adi % lbda_ref / ( adi % rho_ref * adi % cp_ref )
       else
          wrk = adi % D_ref
       end if
       wrk = adi % rho_ref * wrk / ( adi % L_ref * adi % L_ref )

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % scal_var_dissip (sx:ex,sy:ey,sz:ez) ! this variable has implicit sign

       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'SCAL_VAR_Dissipation2' ) then

       if ( inp % dis_rate_therm ) then
          wrk = adi % lbda_ref / ( adi % rho_ref * adi % cp_ref )
       else
          wrk = adi % D_ref
       end if
       wrk = adi % rho_ref * wrk / ( adi % L_ref * adi % L_ref )

       wrk1 (sx:ex,sy:ey,sz:ez) = - stat % scal_var_dissip2 (sx:ex,sy:ey,sz:ez) ! this variable has opposite sign

       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'SCAL_VAR_Dissipation3' ) then

       if ( inp % dis_rate_therm ) then
          wrk = adi % lbda_ref / ( adi % rho_ref * adi % cp_ref )
       else
          wrk = adi % D_ref
       end if
       wrk = adi % rho_ref * wrk / ( adi % L_ref * adi % L_ref )

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % scal_var_dissip3 (sx:ex,sy:ey,sz:ez) ! this variable has implicit sign

       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'SCAL_VAR_Dissipation4' ) then

       if ( inp % dis_rate_therm ) then
          wrk = adi % lbda_ref / ( adi % rho_ref * adi % cp_ref )
       else
          wrk = adi % D_ref
       end if
       wrk = adi % rho_ref * wrk / ( adi % L_ref * adi % L_ref )

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % scal_var_dissip4 (sx:ex,sy:ey,sz:ez) ! this variable has implicit sign

       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'SCAL_VAR_Dissipation5' ) then

       if ( inp % dis_rate_therm ) then
          wrk = adi % lbda_ref / ( adi % rho_ref * adi % cp_ref )
       else
          wrk = adi % D_ref
       end if
       wrk = adi % rho_ref * wrk / ( adi % L_ref * adi % L_ref )

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % scal_var_dissip5 (sx:ex,sy:ey,sz:ez) + & ! this variable has implicit sign
                                  stat % scal_var_dissip5 (sx:ex,sy:ey,sz:ez)

       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )

    else if ( name == 'SCAL_VAR_Dissipation6' ) then

       if ( inp % dis_rate_therm ) then
          wrk = adi % lbda_ref / ( adi % rho_ref * adi % cp_ref )
       else
          wrk = adi % D_ref
       end if
       wrk = adi % rho_ref * wrk / ( adi % L_ref * adi % L_ref )

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % scal_var_dissip6 (sx:ex,sy:ey,sz:ez) ! this variable has implicit sign

       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( wrk , wrk1 )


    ! correlations


    else if ( name == 'CORR_ux' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tmp_corr_ux (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'CORR_vy' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tmp_corr_vy (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'CORR_wz' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tmp_corr_wz (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'CORR_p' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tmp_corr_p (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'CORR_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tmp_corr_Y (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'CORR_Zm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % tmp_corr_Zm (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Lz' ) then

       call integral_length_Z ( inp , grid , stat % tmp_corr_ux , wrk1 )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( adi % L_ref , wrk1 )


    ! other correlations



!> Aditionnal variables (A.Techer)

    else if ( name == 'fa_xiv2' ) then ! (variance SGS - method 2)

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % xiv2 (sx:ex,sy:ey,sz:ez)
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_xiv3' ) then ! (variance SGS - method 3)

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % xiv3 (sx:ex,sy:ey,sz:ez)
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_xixi' ) then ! conservation of the square of the filtred mixture fraction

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % xixi (sx:ex,sy:ey,sz:ez)
       call plot_stat_vtr ( 1.0_dp , wrk1 )

!    else if ( name == 'ra_zjm' ) then ! minus jump position for ignition

!       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % zjm (sx:ex,sy:ey,sz:ez)
!       call plot_stat_vtr ( 1.0_dp , wrk1 )

!    else if ( name == 'ra_zjp' ) then ! plus jump position for ignition

!       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % zjp (sx:ex,sy:ey,sz:ez)
!       call plot_stat_vtr ( 1.0_dp , wrk1 )

!    else if ( name == 'ra_intzjmp' ) then ! ignition probability

!       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % intzjmp (sx:ex,sy:ey,sz:ez)
!       call plot_stat_vtr ( 1.0_dp , wrk1 )

!    else if ( name == 'ra_taumix' ) then ! SGS dissipation time scale

!       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % taumix (sx:ex,sy:ey,sz:ez)
!       call plot_stat_vtr ( adi % time_ref , wrk1 )


    ! species stuff - Reynolds and Favre average


    else if ( name == 'ra_Zm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Zm (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez, max ( inp % index_scalar , 1 ) )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Zm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Zm (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez, max ( inp % index_scalar , 1 ) )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Zm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Zm (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez, max ( inp % index_scalar , 1 ) )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Zm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Zm (sx:ex,sy:ey,sz:ez)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez, max ( inp % index_scalar , 1 ) )
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )


    else if ( name == 'ra_Y01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y05' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y05' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y06' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y06' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y07' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y07' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y08' ) then

       spc = 8
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y08' ) then

       spc = 8
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y09' ) then

       spc = 9
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y09' ) then

       spc = 9
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y11' ) then

       spc = 11
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y11' ) then

       spc = 11
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y12' ) then

       spc = 12
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y12' ) then

       spc = 12
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y13' ) then

       spc = 13
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y13' ) then

       spc = 13
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y14' ) then

       spc = 14
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y14' ) then

       spc = 14
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y15' ) then

       spc = 15
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y15' ) then

       spc = 15
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y16' ) then

       spc = 16
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y16' ) then

       spc = 16
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y17' ) then

       spc = 17
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y17' ) then

       spc = 17
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y18' ) then

       spc = 18
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y18' ) then

       spc = 18
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y19' ) then

       spc = 19
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y19' ) then

       spc = 19
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Y20' ) then

       spc = 20
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Y20' ) then

       spc = 20
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )


    ! species stuff - Reynolds and Favre average - Mass fraction of steady state species (for reduced chemistry)

    else if ( name == 'ra_Yss01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Yss01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Yss02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Yss02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Yss03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Yss03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'ra_Yss04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fa_Yss04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )


    ! species stuff - Reynolds and Favre variance

    else if ( name == 'rv_Y01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y05' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y05' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y06' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y06' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y07' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y07' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y08' ) then

       spc = 8
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y08' ) then

       spc = 8
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y09' ) then

       spc = 9
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y09' ) then

       spc = 9
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y11' ) then

       spc = 11
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y11' ) then

       spc = 11
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y12' ) then

       spc = 12
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y12' ) then

       spc = 12
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y13' ) then

       spc = 13
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y13' ) then

       spc = 13
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y14' ) then

       spc = 14
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y14' ) then

       spc = 14
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y15' ) then

       spc = 15
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y15' ) then

       spc = 15
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y16' ) then

       spc = 16
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y16' ) then

       spc = 16
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y17' ) then

       spc = 17
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y17' ) then

       spc = 17
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y18' ) then

       spc = 18
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y18' ) then

       spc = 18
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y19' ) then

       spc = 19
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y19' ) then

       spc = 19
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Y20' ) then

       spc = 20
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Y20' ) then

       spc = 20
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )


    ! species stuff - Reynolds and Favre variance - Mass fraction of steady state species (for reduced chemistry)

    else if ( name == 'rv_Yss01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Yss01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Yss02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Yss02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Yss03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Yss03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rv_Yss04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % reyvar_Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'fv_Yss04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Yass (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp , wrk1 )


    ! species stuff - Reynolds and Favre average - Production rate

    else if ( name == 'ra_W01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W05' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W05' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W06' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W06' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W07' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W07' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W08' ) then

       spc = 8
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W08' ) then

       spc = 8
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W09' ) then

       spc = 9
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W09' ) then

       spc = 9
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W11' ) then

       spc = 11
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W11' ) then

       spc = 11
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W12' ) then

       spc = 12
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W12' ) then

       spc = 12
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W13' ) then

       spc = 13
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W13' ) then

       spc = 13
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W14' ) then

       spc = 14
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W14' ) then

       spc = 14
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W15' ) then

       spc = 15
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W15' ) then

       spc = 15
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W16' ) then

       spc = 16
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W16' ) then

       spc = 16
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W17' ) then

       spc = 17
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W17' ) then

       spc = 17
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W18' ) then

       spc = 18
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W18' ) then

       spc = 18
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W19' ) then

       spc = 19
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W19' ) then

       spc = 19
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W20' ) then

       spc = 20
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W20' ) then

       spc = 20
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )



   ! radicals


   else if ( name == 'fa_max_OH' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call max_val_Y ( wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp , wrk2 )

   else if ( name == 'fa_max_HO2' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Y (sx:ex,sy:ey,sz:ez,spc)
       call max_val_Y ( wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp , wrk2 )

   else if ( name == 'fv_max_OH' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call max_val_Y ( wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp , wrk2 )

   else if ( name == 'fv_max_HO2' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = stat % favvar_Y (sx:ex,sy:ey,sz:ez,spc)
       call max_val_Y ( wrk1 , wrk2 )
       call similarity_x ( inp , wrk2 )
       call plot_stat_vtr ( 1.0_dp , wrk2 )


    ! HCT MODEL
    else if ( name == 'ra_W01_hct' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W01_hct' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W02_hct' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W02_hct' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W03_hct' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W03_hct' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W04_hct' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W04_hct' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W05_hct' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W05_hct' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'ra_W06_hct' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = reyavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'fa_W06_hct' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = favavg % Ydot_hct (sx:ex,sy:ey,sz:ez,spc)
       call similarity_x ( inp , wrk1 )
       call plot_stat_vtr ( 1.0_dp / adi % time_ref , wrk1 )



    else

       if ( rank == rank_default ) write (*,*) 'error: variable ' , trim (name) , &
                                               ' not found in the database ...'
       wrk1 = 0.0_dp
       call plot_stat_vtr ( 1.0_dp , wrk1 ) ! plot it, just not to have a corrupt file

    end if


    deallocate ( wrk1 , wrk2 , wrk3 , &
                 wrk4 , wrk5 , wrk6 , &
                 wrk7 , wrk8 , wrk9 , &
                 wrk10 )

    deallocate ( wrk_dim4 , wrk_dim5 )


  end subroutine plot_stat_var


end module stat_plot
