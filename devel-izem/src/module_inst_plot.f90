!------------------------------------------------------------------------------
! MODULE: inst_plot
!------------------------------------------------------------------------------
!> \brief ParaView plot module for _instantaneous_ variables.
!!
!! This module allows to plot all types of instantaneous variables in
!! a ParaView file.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module inst_plot


  use parameters
  use parallel
  use input
  use type_thd
  use adim
  use variables
  use tools
  use tools_post


  implicit none


contains


!> \brief Plot a variable in a VTR file.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine plot_inst_vtr ( vref , var )


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


    ! ! Plane XY
    ! if ( rank == 13 ) write (*,*) var (ex+1,ey+1,sz-1:ez+1)
    ! if ( rank == 13 ) write (*,*) var (sx-1,sy-1,sz-1:ez+1)
    ! if ( rank == 13 ) write (*,*) var (ex+1,sy-1,sz-1:ez+1)
    ! if ( rank == 13 ) write (*,*) var (sx-1,ey+1,sz-1:ez+1)
    ! ! Plane XZ
    ! if ( rank == 13 ) write (*,*) var (ex+1,sy-1:ey+1,ez+1)
    ! if ( rank == 13 ) write (*,*) var (sx-1,sy-1:ey+1,sz-1)
    ! if ( rank == 13 ) write (*,*) var (ex+1,sy-1:ey+1,sz-1)
    ! if ( rank == 13 ) write (*,*) var (sx-1,sy-1:ey+1,ez+1)
    ! ! Plane YZ
    ! if ( rank == 13 ) write (*,*) var (sx-1:ex+1,ey+1,ez+1)
    ! if ( rank == 13 ) write (*,*) var (sx-1:ex+1,sy-1,sz-1)
    ! if ( rank == 13 ) write (*,*) var (sx-1:ex+1,ey+1,sz-1)
    ! if ( rank == 13 ) write (*,*) var (sx-1:ex+1,sy-1,ez+1)
    ! call mpi_barrier ( MPI_COMM_WORLD , mpicode )
    ! stop

    write (unit_plot) nx * ny * nz * nbitsimple , &
                      real ( var (ix:fx,iy:fy,iz:fz) * vref , sp )


  end subroutine plot_inst_vtr



!> \brief Selector of instantaneous variables to plot.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine plot_inst_var ( ifile , idnumber , inp , thd , adi , sim , dx_i , dy_i , dz_i , inst )


    integer (ip) , intent (in)                                       :: ifile    !< file number
    integer (ip) , intent (in)                                       :: idnumber !< number of the statistical variable
    type (inp_type) , intent (in)                                    :: inp      !< input derived type
    type (thd_type) , intent (in)                                    :: thd      !< thermodynamic derived type
    type (adi_type) , intent (in)                                    :: adi      !< non-dimensional derived type
    type (sim_type) , intent (in)                                    :: sim      !< simulation derived type
    real (dp) , dimension (:) , allocatable , intent (in)            :: dx_i     !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)            :: dy_i     !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)            :: dz_i     !< inverted dz array
    type (inst_type) , intent (inout)                                :: inst     !< instantaneous derived type


    character (len_default) , parameter :: format = ' ( A11 , 1X , I2 , 1X , A10 , 1X , A20 ) '
    character (len_default)             :: name
    integer (ip)                        :: ok , spc , reac
    integer (ip)                        :: i , j , k
    real (dp)                           :: wrk

    real (dp)                           :: DTmax , dtmin ! for timesteps (A.Techer)


    ! temporal variables
    real (dp) , dimension (:,:,:) , allocatable :: wrk1 , wrk2 , wrk3 , wrk4


    name = inp % var_name (idnumber)
!    if ( rank == rank_default ) write ( * , format ) 'variable #' , idnumber , ', plotting' , trim (name)

    allocate  ( wrk1 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                wrk2 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                wrk3 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                wrk4 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate plot_inst_vars')


    ! basic variables
    if ( name == 'Density' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( adi % rho_ref , wrk1 )

    else if ( name == 'Velocity-X' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'rel_Velocity-X' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = ( inst % u (sx:ex,sy:ey,sz:ez,2) /   &
                                    inst % u (sx:ex,sy:ey,sz:ez,1) ) - &
                                    inp % uc / adi % u_ref
       call plot_inst_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Velocity-Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Velocity-Z' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Total_Energy' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,5) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'Enthalpy' ) then

       call plot_inst_vtr ( adi % u_ref * adi % u_ref , inst % H )

    else if ( name == 'Total_Enthalpy' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) * inst % u (sx:ex,sy:ey,sz:ez,2) + &
                                  inst % u (sx:ex,sy:ey,sz:ez,3) * inst % u (sx:ex,sy:ey,sz:ez,3) + &
                                  inst % u (sx:ex,sy:ey,sz:ez,4) * inst % u (sx:ex,sy:ey,sz:ez,4)
       wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) /                                        &
                                ( inst % u (sx:ex,sy:ey,sz:ez,1) * inst % u (sx:ex,sy:ey,sz:ez,1) )
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % H (sx:ex,sy:ey,sz:ez) + 0.5_dp * wrk1 (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( adi % u_ref * adi % u_ref , wrk1 )

    else if ( name == 'Temperature' ) then

       call plot_inst_vtr ( adi % T_ref , inst % T )

    else if ( name == 'Total_Temp' ) then

       call Ttmix ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 )
       call plot_inst_vtr ( adi % T_ref , wrk1 )

    else if ( name == 'Pressure' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,1) * &
                                  inst % T (sx:ex,sy:ey,sz:ez)   * &
                                  inst % W_i (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( adi % P_ref , wrk1 )

    else if ( name == 'Mean_Molar_Mass' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = 1.0_dp / inst % W_i (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( adi % W_ref , wrk1 )

    else if ( name == 'Sound_Speed' ) then

       call soundmix ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 )
       call plot_inst_vtr ( adi % u_ref , wrk1 )

    else if ( name == 'Cp' ) then

       call plot_inst_vtr ( adi % cp_ref , inst % cp )


    ! timestep variables
    else if ( name == 'dt_min' ) then

       DTmax = 50.0_dp
       dtmin = 1.0_dp
       call timestep ( inp , adi , thd , dx_i , dy_i , dz_i , inst % cp , inst % W_i , inst % T , & 
                       inst % u , DTmax , dtmin , inst % dt )
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % dt (sx:ex,sy:ey,sz:ez,5)
       call plot_inst_vtr ( adi % time_ref , wrk1 )

    else if ( name == 'dt_CFL' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % dt (sx:ex,sy:ey,sz:ez,1) !/ inst % dt (sx:ex,sy:ey,sz:ez,5)
       call plot_inst_vtr ( adi % time_ref , wrk1 )

    else if ( name == 'dt_mass' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % dt (sx:ex,sy:ey,sz:ez,2) !/ inst % dt (sx:ex,sy:ey,sz:ez,5)
       call plot_inst_vtr ( adi % time_ref , wrk1 )

    else if ( name == 'dt_therm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % dt (sx:ex,sy:ey,sz:ez,3) !/ inst % dt (sx:ex,sy:ey,sz:ez,5)
       call plot_inst_vtr ( adi % time_ref , wrk1 )

    else if ( name == 'dt_chem' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % dt (sx:ex,sy:ey,sz:ez,4) !/ inst % dt (sx:ex,sy:ey,sz:ez,5)
       call plot_inst_vtr ( adi % time_ref , wrk1 )


    ! viscous variables
    else if ( name == 'Dyn_Viscosity' ) then

       call plot_inst_vtr ( adi % mu_ref , inst % mu )

    else if ( name == 'SGS_Viscosity' ) then

       if ( .not. inp % LES ) call abort_mpi ('inp % LES == .false. can not plot SGS_Viscosity...')
       call plot_inst_vtr ( adi % mu_ref , inst % mu_SGS )

    else if ( name == 'Kin_Viscosity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % mu (sx:ex,sy:ey,sz:ez) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( adi % mu_ref / adi % rho_ref , wrk1 )

    else if ( name == 'Bulk' ) then

       call plot_inst_vtr ( adi % mu_ref , inst % kpa )

    else if ( name == 'Therm_Cond' ) then

       call plot_inst_vtr ( adi % lbda_ref , inst % ct )

    else if ( name == 'Alphaz' ) then

       if ( inp % dis_rate_therm ) then
          wrk = adi % lbda_ref / ( adi % rho_ref * adi % cp_ref )
       else
          wrk = adi % D_ref
       end if

       call Zmsimp ( inp , thd , inst % u , wrk1 )
       call D_equivalent ( inp , thd                           , &
                           dx_i      , dy_i       , dz_i       , &
                           inst % u  , inst % T   , inst % W_i , &
                           inst % ct , inst % cp  ,              &
                           inst % dm , inst % tdr , inst % rd  , &
                           wrk1      , wrk2 )
       call plot_inst_vtr ( wrk , wrk2 )

    else if ( name == 'Dissip_rate' ) then

       if ( inp % dis_rate_therm ) then
          wrk = ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / &
                ( adi % L_ref * adi % L_ref )
       else
          wrk = adi % D_ref / ( adi % L_ref * adi % L_ref )
       end if

       if ( npv == 0 ) then ! calculated passive scalar
          call Zmsimp ( inp , thd , inst % u , wrk1 )
       else                 ! real passive scalar
          wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez, niv+nrv+npv ) / &
                                     inst % u (sx:ex,sy:ey,sz:ez, 1 )
       end if

       call D_equivalent ( inp , thd                           , &
                           dx_i      , dy_i       , dz_i       , &
                           inst % u  , inst % T   , inst % W_i , &
                           inst % ct , inst % cp  ,              &
                           inst % dm , inst % tdr , inst % rd  , &
                           wrk1      , wrk2 )
       call dissip_rate ( dx_i , dy_i , dz_i , wrk2 , wrk1 , wrk3 )
       call plot_inst_vtr ( wrk , wrk3 )


    ! non dimensional variables
    else if ( name == 'Gamma' ) then

       call gammamix ( 0 , thd , inst % W_i , inst % cp , wrk1 )
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Mach' ) then

       call machmix ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 )
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'rel_Mach' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = ( inst % u (sx:ex,sy:ey,sz:ez,2) /   &
                                    inst % u (sx:ex,sy:ey,sz:ez,1) ) - &
                                    inp % uc / adi % u_ref
       wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3)
       wrk3 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4)
       call machmix_rel ( 0 , thd , inst % W_i , inst % cp , inst % T , inst % u , wrk1 , wrk2 , wrk3 , wrk4 )
       call plot_inst_vtr ( 1.0_dp , wrk4 )

    else if ( name == 'Prandtl' ) then

       if ( .not. vis ) write (*,*) 'calculating Prandtl number in an inviscid problem ...'
       call prandtlmix ( 0 , inst % cp , inst % mu , inst % ct , wrk1 )
       call plot_inst_vtr ( adi % Pr , wrk1 )


    ! other variables
    else if ( name == 'Divergence' ) then

       call divergence ( dx_i , dy_i , dz_i , inst % u , wrk1 )
       call plot_inst_vtr ( adi % u_ref / adi % L_ref , wrk1 )

    else if ( name == 'Vorticity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       wrk3 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)

       call vorticity ( dx_i , dy_i , dz_i , wrk1 , wrk2 , wrk3 , wrk4 )
       call plot_inst_vtr ( adi % u_ref / adi % L_ref , wrk4 )

    else if ( name == 'Vorticity-X' ) then

       wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       wrk3 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)

       call vorticity_x ( dy_i , dz_i , wrk2 , wrk3 , wrk4 )
       call plot_inst_vtr ( adi % u_ref / adi % L_ref , wrk4 )

    else if ( name == 'Vorticity-Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       wrk3 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)

       call vorticity_y ( dx_i , dz_i , wrk1 , wrk3 , wrk4 )
       call plot_inst_vtr ( adi % u_ref / adi % L_ref , wrk4 )

    else if ( name == 'Vorticity-Z' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)

       call vorticity_z ( dx_i , dy_i , wrk1 , wrk2 , wrk4 )
       call plot_inst_vtr ( adi % u_ref / adi % L_ref , wrk4 )

    else if ( name == 'Heat_Rel' ) then

       call plot_inst_vtr ( adi % u_ref * adi % u_ref / adi % time_ref , inst % omega )

    else if ( name == 'Q' ) then

       wrk = adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref )
       call Qcriterion ( dx_i , dy_i , dz_i , inst % u , wrk1 )
       call plot_inst_vtr ( wrk , wrk1 )

    else if ( name == 'Lambda2' ) then

       wrk = adi % u_ref * adi % u_ref / ( adi % L_ref * adi % L_ref )
       call lbda2criterion ( dx_i , dy_i , dz_i , inst % u , wrk1 )
       call plot_inst_vtr ( wrk , wrk1 )

    else if ( name == 'Takeno' ) then

       call takeno ( inp , dx_i , dy_i , dz_i , inst % omega , inst % u , wrk1 )
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Detector_shock' ) then

       call shock_det_post ( inst % u , inst % T , inst % W_i , wrk1 )
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Detector_shock_Ducros' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       wrk3 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)

       call shock_det_ducros_post ( dx_i , dy_i , dz_i , wrk1 , wrk2 , wrk3 , wrk4 )
       call plot_inst_vtr ( 1.0_dp , wrk4 )

    else if ( name == 'Detector_reaction' ) then

       call reaction_det ( ifile , inp , adi , sim , inst % u , inst % Ydot , inst % T , wrk1 )
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Damkohler' ) then

       wrk = ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )
       call Damkohler ( inp , thd , dx_i , dy_i , dz_i , inst % u ,                  &
                        inst % T , inst % W_i , inst % dm , inst % tdr , inst % rd , &
                        inst % Ydot , wrk1 )
       call plot_inst_vtr ( wrk , wrk1 )

    else if ( name == 'Filter_width' ) then

       call filter_width ( 0 , dx_i , dy_i , dz_i , wrk1 )
       call plot_inst_vtr ( adi % L_ref , wrk1 )

    else if ( name == 'Rank' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = rank
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'zjm' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % zjm (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'zjp' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % zjp (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'intzjmp' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % intzjmp (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'taumix' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % taumix (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( adi % time_ref , wrk1 )



    ! schieleren variables
    else if ( name == 'sDensity' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,1)
       call schlieren ( dx_i , dy_i , dz_i , wrk1 , wrk2 )

       call plot_inst_vtr ( 1.0_dp , wrk2 )

    else if ( name == 'sVelocity-X' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,2) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call schlieren ( dx_i , dy_i , dz_i , wrk1 , wrk2 )

       call plot_inst_vtr ( 1.0_dp , wrk2 )

    else if ( name == 'sVelocity-Y' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,3) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call schlieren ( dx_i , dy_i , dz_i , wrk1 , wrk2 )

       call plot_inst_vtr ( 1.0_dp , wrk2 )

    else if ( name == 'sVelocity-Z' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,4) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call schlieren ( dx_i , dy_i , dz_i , wrk1 , wrk2 )

       call plot_inst_vtr ( 1.0_dp , wrk2 )

    else if ( name == 'sTemperature' ) then

       call schlieren ( dx_i , dy_i , dz_i , inst % T , wrk1 )

       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'sPressure' ) then

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,1) * &
                                  inst % T (sx:ex,sy:ey,sz:ez)   * &
                                  inst % W_i (sx:ex,sy:ey,sz:ez)
       call schlieren ( dx_i , dy_i , dz_i , wrk1 , wrk2 )

       call plot_inst_vtr ( 1.0_dp , wrk2 )






    ! species stuff
    else if ( name == 'Zm' ) then

       call Zmsimp ( inp , thd , inst % u , wrk1 )
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'xiv2' ) then ! variance SGS (methode 2)

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % xiv2 (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'xiv3' ) then ! variance SGS (methode 3)

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % xiv3 (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'xixi' ) then ! conservation of the square of the filtred mixture fraction

       wrk1 (sx:ex,sy:ey,sz:ez) = inst % xixi (sx:ex,sy:ey,sz:ez)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'S_Zi' ) then ! Segregation for variance equation (methode 1)

       spc = 6 ! \xi_{v,i} = Variance_I
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       spc = 5 ! \tilde{\xi}
       wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)

       wrk3 (sx:ex,sy:ey,sz:ez) = 0.0_dp

       do k = sz , ez
          do j = sy , ey
             do i = sx , ex
                if ( wrk2 (i,j,k) > 1e-3_dp .and. wrk2 (i,j,k) < (1.0_dp - 1e-3_dp) ) then
                   wrk = wrk2 (i,j,k) * ( 1.0_dp - wrk2 (i,j,k) )
                   wrk1 (i,j,k) = max ( min ( wrk1 (i,j,k) , wrk ) , 0.0_dp )

                   wrk3 (i,j,k) = wrk1 (i,j,k) / wrk

                   wrk3 (i,j,k) = max ( min ( wrk3 (i,j,k) , 1.0_dp ) , 0.0_dp )
                   if ( wrk3 (i,j,k) /= wrk3 (i,j,k) ) wrk3 (i,j,k) = 0.0_dp
                end if
             end do
          end do
       end do

       call plot_inst_vtr ( 1.0_dp , wrk3 )

    else if ( name == 'S_Zii' ) then ! Segregation for rebuilt variance (methode 2)

       spc = 7 ! \tilde{\xi^2} = ps2
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       spc = 5 ! \tilde{\xi}
       wrk2 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)

       ! Rebuilt variance (overwrite wrk1)
       wrk1 (sx:ex,sy:ey,sz:ez) = wrk1 (sx:ex,sy:ey,sz:ez) - wrk2 (sx:ex,sy:ey,sz:ez) &
                                                           * wrk2 (sx:ex,sy:ey,sz:ez)

       wrk3 (sx:ex,sy:ey,sz:ez) = 0.0_dp

       do k = sz , ez
          do j = sy , ey
             do i = sx , ex
                if ( wrk2 (i,j,k) > 1e-3_dp .and. wrk2 (i,j,k) < (1.0_dp - 1e-3_dp) ) then
                   wrk = wrk2 (i,j,k) * ( 1.0_dp - wrk2 (i,j,k) )
                   wrk1 (i,j,k) = max ( min ( wrk1 (i,j,k) , wrk ) , 0.0_dp )

                   wrk3 (i,j,k) = wrk1 (i,j,k) / wrk

                   wrk3 (i,j,k) = max ( min ( wrk3 (i,j,k) , 1.0_dp ) , 0.0_dp )
                   if ( wrk3 (i,j,k) /= wrk3 (i,j,k) ) wrk3 (i,j,k) = 0.0_dp
                end if
             end do
          end do
       end do

       call plot_inst_vtr ( 1.0_dp , wrk3 )

    else if ( name == 'Y01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y05' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y06' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y07' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y08' ) then

       spc = 8
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y09' ) then

       spc = 9
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y11' ) then

       spc = 11
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y12' ) then

       spc = 12
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y13' ) then

       spc = 13
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y14' ) then

       spc = 14
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y15' ) then

       spc = 15
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y16' ) then

       spc = 16
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y17' ) then

       spc = 17
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y18' ) then

       spc = 18
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y19' ) then

       spc = 19
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'Y20' ) then

       spc = 20
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % u (sx:ex,sy:ey,sz:ez,niv+spc) / &
                                  inst % u (sx:ex,sy:ey,sz:ez,1)
       call plot_inst_vtr ( 1.0_dp , wrk1 )

    else if ( name == 'W01' ) then

       spc = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W02' ) then

       spc = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W03' ) then

       spc = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W04' ) then

       spc = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W05' ) then

       spc = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W06' ) then

       spc = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W07' ) then

       spc = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W08' ) then

       spc = 8
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W09' ) then

       spc = 9
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W10' ) then

       spc = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W11' ) then

       spc = 11
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W12' ) then

       spc = 12
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W13' ) then

       spc = 13
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W14' ) then

       spc = 14
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W15' ) then

       spc = 15
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W16' ) then

       spc = 16
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W17' ) then

       spc = 17
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W18' ) then

       spc = 18
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W19' ) then

       spc = 19
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W20' ) then

       spc = 20
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % Ydot (sx:ex,sy:ey,sz:ez,spc)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    ! rate of progress of each elementary step
    else if ( name == 'W_scal_i01' ) then

       reac = 1
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i02' ) then

       reac = 2
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i03' ) then

       reac = 3
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i04' ) then

       reac = 4
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i05' ) then

       reac = 5
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i06' ) then

       reac = 6
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i07' ) then

       reac = 7
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i08' ) then

       reac = 8
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i09' ) then

       reac = 9
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i10' ) then

       reac = 10
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i11' ) then

       reac = 11
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i12' ) then

       reac = 12
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i13' ) then

       reac = 13
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i14' ) then

       reac = 14
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i15' ) then

       reac = 15
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i16' ) then

       reac = 16
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i17' ) then

       reac = 17
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i18' ) then

       reac = 18
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i19' ) then

       reac = 19
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i20' ) then

       reac = 20
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i21' ) then

       reac = 21
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i22' ) then

       reac = 22
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i23' ) then

       reac = 23
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i24' ) then

       reac = 24
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i25' ) then

       reac = 25
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i26' ) then

       reac = 26
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i27' ) then

       reac = 27
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i28' ) then

       reac = 28
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i29' ) then

       reac = 29
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else if ( name == 'W_scal_i30' ) then

       reac = 30
       wrk1 (sx:ex,sy:ey,sz:ez) = inst % W_scal_i (sx:ex,sy:ey,sz:ez,reac)
       call plot_inst_vtr ( 1.0_dp / adi % time_ref , wrk1 )

    else

       if ( rank == rank_default ) write (*,*) 'error: variable ' , trim (name) , &
                                               ' not found in the database ...'
       wrk1 = 0.0_dp
       call plot_inst_vtr ( 1.0_dp , wrk1 ) ! plot it anyway to avoid a corrupt file

    end if


    deallocate ( wrk1 , wrk2 , wrk3 , wrk4 )


  end subroutine plot_inst_var


end module inst_plot
