!------------------------------------------------------------------------------
! MODULE: tools_post
!------------------------------------------------------------------------------
!> \brief Post-treatment tools and utilities.
!!
!! Contrarily to module_tools.f90, this module provides all the tools
!! and utilites used _exclusively_ during the post-treatment.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module tools_post

  use parameters
  use parallel
  use type_thd
  use input
  use adim
  use variables
  use deriv
  use eg_lib
  use tools , only : pdfbeta , gaussian , log_normal
  use SGS_models
  use BCs , only : bc_extrapolation
  use solver_reaction


  implicit none


  real (dp) , parameter , private :: epsi     = 1.0e-10_dp                !< to avoid divisions by zero
  real (dp) , parameter , private :: minYf    = 1.0e-4_dp                 !< minimum mass fraction
  real (dp) , parameter , private :: minZm    = 0.01_dp , maxZm = 0.99_dp !< mixture fraction bounds
  real (dp) , parameter , private :: minomega = 1.0e+0_dp                 !< minimum heat release


contains


!> \brief Normalize variales for the post-treatment.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine adimensionalize_post ( inp , adi , x , y , z , xt , yt , zt , dx_i , dy_i , dz_i )


    type (inp_type) , intent (inout)                           :: inp  !< input derived type
    type (adi_type) , intent (in)                              :: adi  !< non-dimensional derived type
    real (dp) , allocatable , dimension (:) , intent (inout)   :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (inout)   :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:) , intent (inout)   :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:) , intent (inout)   :: xt   !< absolute x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (inout)   :: yt   !< absolute y-coordinate array
    real (dp) , allocatable , dimension (:) , intent (inout)   :: zt   !< absolute z-coordinate array
    real (dp) , allocatable , dimension (:) , intent (inout)   :: dx_i !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (inout)   :: dy_i !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (inout)   :: dz_i !< inverted dz array


    integer (ip) :: l
    real (dp)    :: L_i


    L_i = 1.0_dp / adi % L_ref

    inp % initial_time = inp % initial_time / adi % time_ref
    inp % timing (:)   = inp % timing (:) / adi % time_ref

    x (:) = x (:) * L_i ; xt (:) = xt (:) * L_i
    y (:) = y (:) * L_i ; yt (:) = yt (:) * L_i
    z (:) = z (:) * L_i ; zt (:) = zt (:) * L_i

    dx_i (:) = dx_i (:) * adi % L_ref
    dy_i (:) = dy_i (:) * adi % L_ref
    dz_i (:) = dz_i (:) * adi % L_ref

    do l = 1 , inp % nvolume
       inp % x_volmin (l) = inp % x_volmin (l) * inp % dim_length_coord * L_i
       inp % x_volmax (l) = inp % x_volmax (l) * inp % dim_length_coord * L_i
       inp % y_volmin (l) = inp % y_volmin (l) * inp % dim_length_coord * L_i
       inp % y_volmax (l) = inp % y_volmax (l) * inp % dim_length_coord * L_i
       inp % z_volmin (l) = inp % z_volmin (l) * inp % dim_length_coord * L_i
       inp % z_volmax (l) = inp % z_volmax (l) * inp % dim_length_coord * L_i
    end do

    do l = 1 , inp % nxystat
       inp % z_xystat (l) = inp % z_xystat (l) * inp % dim_length_coord * L_i
    end do

    do l = 1 , inp % nxzstat
       inp % y_xzstat (l) = inp % y_xzstat (l) * inp % dim_length_coord * L_i
    end do

    do l = 1 , inp % nyzstat
       inp % x_yzstat (l) = inp % x_yzstat (l) * inp % dim_length_coord * L_i
    end do


    ! post.data file
    inp % s_x = inp % s_x * L_i ; inp % e_x = inp % e_x * L_i
    inp % s_y = inp % s_y * L_i ; inp % e_y = inp % e_y * L_i
    inp % s_z = inp % s_z * L_i ; inp % e_z = inp % e_z * L_i

    ! inp % probe_coord (:,:)  = inp % probe_coord (:,:) * inp % dim_length_coord !* L_i
    ! inp % corrspec_coord (:) = inp % corrspec_coord (:) * inp % dim_length_coord * L_i


  end subroutine adimensionalize_post


!> \brief Read times for statistics.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine readtimes ( inp , adi , sim , time , dtime )


    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    type (sim_type) , intent (inout)                               :: sim   !< simulation derived type
    real (dp) , allocatable , dimension (:) , intent (inout)       :: time  !< time
    real (dp) , allocatable , dimension (:) , intent (inout)       :: dtime !< time step


    integer (ip)                        :: ok , ifile , idummy , indicator
    character (len_default) , parameter :: format  = ' ( A8 , 1X , I3 ) '


    allocate ( time  (0:ntimes) , &
               dtime (0:ntimes) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate readtimes')
    time (:)  = 0.0_dp
    dtime (:) = 0.0_dp


    open ( unit = unit_time_stat   , file = trim (dir_parent) // trim (file_time_stat) , &
           form = 'formatted' ,  status = 'old' , action = 'read' , iostat = ok )

    if ( ok > 0 ) call abort_mpi ('error opening ' // trim (dir_parent) // trim (file_time_stat))

    indicator = inp % start_file
    do ifile = inp % start_file , inp % end_file , inp % skip_file
       indicator = indicator + inp % skip_file
    end do
    indicator = indicator - inp % skip_file
    sim % correct_endfile = indicator

    do ifile = 1 , inp % end_file ! read from 1 always and the whole interval
       read ( unit_time_stat , * ) idummy , time (ifile) , sim % dtime_per_ite (ifile)
    end do

    do ifile = inp % start_file , inp % end_file - inp % skip_file , inp % skip_file ! skip the last file
       dtime (ifile) =  time ( ifile + inp % skip_file ) - time (ifile)
    end do
    dtime (indicator) = dtime ( max ( 0 , indicator - inp % skip_file ) ) ! last file (A.Techer modif 1 -> 0 for initial solution)

    sim % sumdtime = 0.0_dp
    do ifile = inp % start_file , inp % end_file , inp % skip_file
       sim % sumdtime = sim % sumdtime + dtime (ifile)
    end do

    ! adimensionalize
    time (:)        = time (:) / adi % time_ref
    dtime (:)       = dtime (:) / adi % time_ref
    sim % sumdtime  = sim % sumdtime / adi % time_ref

    close (unit_time_stat)


  end subroutine readtimes


!> \brief Read statistical file.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine read_stat ( ifile , inp , thd , adi , dx_i , dy_i , dz_i , rweg , iweg , inst )


    integer (ip) , intent (in)                                     :: ifile !< file number
    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamical derived type
    type (adi_type) , intent (in)                                  :: adi   !< non-dimensional derived type
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dx_i      !< inverted dx array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dy_i      !< inverted dy array
    real (dp) , dimension (:)       , allocatable , intent (in)    :: dz_i      !< inverted dz array
    real (dp) , allocatable , dimension (:) , intent (inout)       :: rweg  !< real EGLIB work array
    integer (ip) , allocatable , dimension (:) , intent (inout)    :: iweg  !< integer EGLIB work array
    type (inst_type) , intent (inout)                              :: inst  !< instantaneous derived type


    integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip)                                                   :: reclmax , ok , i , j , k , l
    integer (ip)                                                   :: myrec , px , py , pz , ix , fx , iy , fy , iz , fz
    integer (ip) , dimension (ndimmax)                             :: mycoords , mydims
    character (len_default)                                        :: number_ascii , rank_ascii , plane_ascii , dir_stat , file_type
    real (dp)                                                      :: rho_i , rho_cgs , rho_phys , T_phys , W_cgs, &
    & W_phys , wrk
    real (dp)                                                      :: yo2 , zm , zp2 , zp2max , xb , dm_SGS , &
    & taumix , Sc_SGS_i , rvd0
    real (dp)                                                      :: beta1 , beta2 , deltazj
    real (dp) , dimension (nrv+npv+nvv)                            :: Ya
    real (dp) , dimension (nreac)                                  :: kb , kf , qi , W_scal_i
    real (dp) , dimension (:,:,:,:) , allocatable                  :: Xa , ha
    real (dp) , dimension (:,:,:)   , allocatable                  :: delta_2
    real (dp) , parameter                                          :: eps2 = 1.0e-2_dp , eps3 = 1.0e-3_dp , eps4 = 1.0e-4_dp


    write ( number_ascii , format_restart ) ifile
    write ( rank_ascii , format_restart ) rank
    write ( plane_ascii , format_nplane ) inp % plane_number


    if ( inp % plane_type == 'volume' ) then
       dir_stat  = dir_statvol
       file_type = trim (file_statvol) // '_' // trim (plane_ascii)
    else if ( inp % plane_type == 'XY' ) then
       dir_stat  = dir_statXY
       file_type = trim (file_statXY) // '_' // trim (plane_ascii)
    else if ( inp % plane_type == 'XZ' ) then
       dir_stat  = dir_statXZ
       file_type = trim (file_statXZ) // '_' // trim (plane_ascii)
    else if ( inp % plane_type == 'YZ' ) then
       dir_stat  = dir_statYZ
       file_type = trim (file_statYZ) // '_' // trim (plane_ascii)
    else if ( inp % plane_type == 'restart' ) then
       dir_stat  = dir_restart
       file_type = file_restart
    else
       call abort_mpi (trim(inp % plane_type) // ' is a not defined type of plane')
    end if


    if ( nproc == 1 ) then ! sequential problem


       if (ind_files) call abort_mpi ('cannot read individual binary files with one process!')

       if ( inp % plane_type == 'XY' .and. inp % nprocz /= 1 ) &
            call abort_mpi ('error: post-process plan XY with more than 1 proc in z-direction')

       if ( inp % plane_type == 'XZ' .and. inp % nprocy /= 1 ) &
            call abort_mpi ('error: post-process plan XZ with more than 1 proc in y-direction')

       if ( inp % plane_type == 'YZ' .and. inp % nprocx /= 1 ) &
            call abort_mpi ('error: post-process plan YZ with more than 1 proc in x-direction')


       reclmax = ceiling ( 1.0_dp * ntx / inp % nprocx + inp % ghostptx ) * &
                 ceiling ( 1.0_dp * nty / inp % nprocy + inp % ghostpty ) * &
                 ceiling ( 1.0_dp * ntz / inp % nprocz + inp % ghostptz ) * &
                 nv * nbit

       if ( inp % read_BE ) then

          open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' // trim (number_ascii) , &
                                access    = 'direct'        ,                                                   &
                                recl      = reclmax         ,                                                   &
                                form      = 'unformatted'   ,                                                   &
                                status    = 'old'           ,                                                   &
                                convert   = 'big_endian'    ,                                                   &
                                iostat    =  ok             )

       else

          open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' // trim (number_ascii) , &
                                access    = 'direct'        ,                                                   &
                                recl      = reclmax         ,                                                   &
                                form      = 'unformatted'   ,                                                   &
                                status    = 'old'           ,                                                   &
                                iostat    =  ok             )

       end if
       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_stat) // trim (file_type) // '_' // trim (number_ascii))

       mydims (3)   = inp % nprocx
       mydims (2)   = inp % nprocy
       mydims (1)   = inp % nprocz
       mycoords (:) = 0
       myrec        = 1


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

!               write (*,'(9(1X,I7))') myrec , mycoords (3) , mycoords (2) , mycoords (1) , fx-ix+1 , fy-iy+1 , fz-iz+1

               read ( unit_restart , rec = myrec )            & ! it must start at 1
                    (((( inst % u (i,j,k,l) , i = ix , fx ) , &
                                              j = iy , fy ) , &
                                              k = iz , fz ) , &
                                              l = 1  , nv )

               myrec = myrec + 1

            end do
         end do
      end do

      close (unit_restart)


    else ! parallel problem


       if (ind_files) then ! reading individual binary files per process


          reclmax = ( (ex-sx+1) + inp % ghostptx ) * &
                    ( (ey-sy+1) + inp % ghostpty ) * &
                    ( (ez-sz+1) + inp % ghostptz ) * &
                    nv * nbit

          if ( inp % read_BE ) then

             open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' //     &
                                               trim (number_ascii) // '_' // trim (rank_ascii) , &
                                   access    = 'direct'        ,                                 &
                                   recl      = reclmax         ,                                 &
                                   form      = 'unformatted'   ,                                 &
                                   status    = 'old'           ,                                 &
                                   convert   = 'big_endian'    ,                                 &
                                   iostat    =  ok             )

          else

             open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' //     &
                                               trim (number_ascii) // '_' // trim (rank_ascii) , &
                                   access    = 'direct'        ,                                 &
                                   recl      = reclmax         ,                                 &
                                   form      = 'unformatted'   ,                                 &
                                   status    = 'old'           ,                                 &
                                   iostat    =  ok             )

          end if
          if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_stat) // trim (file_type) // '_' //    &
                                                             trim (number_ascii) // '_' // trim (rank_ascii))

          read ( unit_restart , rec = 1 )                &
               (((( inst % u (i,j,k,l) , i = sx , ex ) , &
                                         j = sy , ey ) , &
                                         k = sz , ez ) , &
                                         l = 1  , nv )

          close (unit_restart)


       else ! one binary file


          reclmax = ( nxmax + inp % ghostptx ) * ( nymax + inp % ghostpty ) * ( nzmax + inp % ghostptz ) * nv * nbit

          if ( inp % read_BE ) then

             open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' // trim (number_ascii) , &
                                   access    = 'direct'        ,                                                   &
                                   recl      = reclmax         ,                                                   &
                                   form      = 'unformatted'   ,                                                   &
                                   status    = 'old'           ,                                                   &
                                   convert   = 'big_endian'    ,                                                   &
                                   iostat    =  ok             )

          else

             open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' // trim (number_ascii) , &
                                   access    = 'direct'        ,                                                   &
                                   recl      = reclmax         ,                                                   &
                                   form      = 'unformatted'   ,                                                   &
                                   status    = 'old'           ,                                                   &
                                   iostat    =  ok             )

          end if
          if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_stat) // trim (file_type) // '_' // trim (number_ascii))

          read ( unit_restart , rec = rank+1 )           & ! it must start at 1
               (((( inst % u (i,j,k,l) , i = sx , ex ) , &
                                         j = sy , ey ) , &
                                         k = sz , ez ) , &
                                         l = 1  , nv )

          close (unit_restart)


       end if ! ind_files


    end if ! (nproc == 1)

   


    ! verifying the fields
    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             if ( inst % u (i,j,k,1) <= 0.0_dp ) then
                write (*,'(A,3I7,1PE16.8)') &
                   'error: density is either zero or negative (i,j,k,inst%rho) = ' , &
                   i,j,k , inst % u (i,j,k,1)
                call abort_mpi (' ')
             end if
          end do
       end do
    end do


    allocate ( Xa   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
               ha   ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate read_stat')


    ! approximation of the temperature
    inst % T (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) = inst % T0 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)


    ! filling the ghost points to avoid errors
    call comm_cons ( inst % u )
    call prim_inv_var ( 0 , thd , inst % u , inst % W_i , inst % T , inst % cp , ha )

    ! upd_boundaries: only extrapolation
    do l = 1 , ndim+ndim
       if ( neigh (l) == MPI_PROC_NULL .or. bc (l) == periodic ) then
          call bc_extrapolation ( l , inst % u )
       end if
    end do
    do l = 1 , ndim+ndim
       call prim_inv_var ( face_domain (l) , thd , inst % u , inst % W_i , inst % T , inst % cp , ha )
    end do


    ! approximation to the temperature
    inst % T0 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) = inst % T (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)


    ! ha -> H
    do k = sz,ez
       do j = sy,ey
          do i = sx,ex

             rho_i = 1.0_dp / inst % u (i,j,k,1)

             inst % H (i,j,k) = 0.0_dp
             do l = 1 , nrv
                inst % H (i,j,k) = inst % H (i,j,k) + &
                                   inst % u (i,j,k,niv+l) * rho_i * ha (i,j,k,l)
             end do

          end do
       end do
    end do


    ! viscous variables
    if ( vis .and. .not. eglib ) then
       call molefrac ( 0 , thd , inst % W_i , inst % u , Xa )
       call prim_vis_vars ( 0 , thd , inst % u , inst % W_i , inst % T , &
                            Xa , inst % dm , inst % mu , inst % ct )
    else if ( vis .and. eglib ) then
       call egvars ( thd , inst % T , inst % W_i , inst % u , iweg , rweg , &
                     inst % rd , inst % mu , inst % kpa , inst % ct , inst % tdr )
    end if
    if ( vis .and. LES ) then
       allocate  ( delta_2 ( sx:ex , sy:ey , sz:ez )                   , &
                   stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate read_stat LES')
       call filter_width ( 0 , dx_i , dy_i , dz_i , delta_2 )
       delta_2 = delta_2 * delta_2

       call mu_SGS_selector_post ( inp , adi , dx_i , dy_i , dz_i , inst % u , inst % mu_SGS )
    end if


    ! reactive variables
    if ( reaction .and. inp % reac_selec == dvode ) then

       rho_cgs = adi % rho_ref * 1.0e-3_dp
       W_cgs   = adi % W_ref * 1.0e3_dp

       ! domain to apply combustion (not in every boundary)
       ! if ( ndim >=2 ) then
       !    ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 )
       !    iy = max ( ng+1 , sy ) ; fy = min ( ey , nty-ng-1 )
       !    iz = sz                ; fz = ez
       ! else
       !    ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 )
       !    iy = sy                ; fy = ey
       !    iz = sz                ; fz = ez
       ! end if
       ix = sx ; fx = ex
       iy = sy ; fy = ey
       iz = sz ; fz = ez

       inst % Ydot (sx:ex,sy:ey,sz:ez,:) = 0.0_dp
       inst % omega (sx:ex,sy:ey,sz:ez)  = 0.0_dp

       do k = iz , fz
          do j = iy , fy
             do i = ix , fx

                rho_i = 1.0_dp / inst % u (i,j,k,1)
                do l = 1 , nrv
                   Ya (l) = inst % u (i,j,k, niv+l ) * rho_i
                end do

                rho_phys = inst % u (i,j,k,1) * rho_cgs
                T_phys   = inst % T (i,j,k) * adi % T_ref
                W_phys   = ( 1.0_dp / inst % W_i (i,j,k) ) * W_cgs

                ! here Ya is the non dimensional production rate dY/dt
                call psr_simp_chemkin_qi ( T_phys , rho_phys , W_phys , Ya , adi % time_ref , &
                                           inp % index_scalar , kb , kf , qi , W_scal_i )
                ! call psr_simp_chemkin ( T_phys , rho_phys , Ya , adi % time_ref ) ! here adi % time_ref

                do l = 1 , nrv-1
                   inst % Ydot (i,j,k,l) = Ya (l)
                end do
                inst % Ydot (i,j,k,nrv) = 0.0_dp - sum ( inst % Ydot (i,j,k,1:nrv-1) )

                inst % omega (i,j,k) = 0.0_dp
                do l = 1 , nrv
                   inst % omega (i,j,k) = inst % omega (i,j,k) - &
                                          thd % h0fc (l) * inst % Ydot (i,j,k,l)
                end do
                inst % omega (i,j,k) = max ( 0.0_dp , inst % omega (i,j,k) )

                do l = 1 , nreac
                   inst % W_scal_i (i,j,k,l) = W_scal_i (l)
                end do

             end do
          end do
       end do

    else if ( reaction .and. inp % reac_selec == MIL ) then

       ix = sx ; fx = ex
       iy = sy ; fy = ey
       iz = sz ; fz = ez

       inst % Ydot (sx:ex,sy:ey,sz:ez,:) = 0.0_dp
       inst % omega (sx:ex,sy:ey,sz:ez)  = 0.0_dp
       Sc_sgs_i = 1.0_dp / inp % Sc_sgs
       rvd0     = adi % sqgmr / adi % Sc ! = Ma sqrt(gamma) / ( Re Sc )

       do k = iz , fz
          do j = iy , fy
             do i = ix , fx

                rho_i = 1.0_dp / inst % u (i,j,k,1)
                do l = 1 , nrv+npv+nvv
                   Ya (l) = inst % u (i,j,k, niv+l ) * rho_i
                end do

                yo2      = Ya (inp % index_O2)              ! mean oxygen mass fraction
                zm       = Ya (nrv+npv)                     ! mean mixture fractionn value
                zp2      = Ya (nrv+npv+1)                   ! mixture fraction variance value

                zp2max   = zm * ( max_Z - zm )
                zp2      = min ( max ( zp2 , epsi ) , zp2max )

                ! Relaxation time ( + eps to avoid divisions by zero )
                dm_SGS   = adi % Sc * inst % mu_SGS (i,j,k) * rho_i * Sc_SGS_i
                taumix   = delta_2 (i,j,k) / ( rvd0 * dm_SGS + epsi )
                inst % taumix (i,j,k) = taumix

                taumix   = taumix * adi % time_ref ! dimensionalize in second

                ! Jump position
                inst % zjm (i,j,k) = zm
                inst % zjp (i,j,k) = zm
                call jump_mp ( inp , taumix , 1.0_dp , inst % zjm (i,j,k) , inst % zjp (i,j,k) )

                inst % zjm (i,j,k) = min ( inst % zjm (i,j,k) , zm )
                inst % zjp (i,j,k) = max ( inst % zjp (i,j,k) , zm )
                deltazj            = inst % zjp (i,j,k) - inst % zjm (i,j,k)

                ! Ignition probability

!                if ( zm <= eps2 .or. zm >= (1.0_dp - eps2) ) then
                if ( zp2 <= eps4 * zp2max .or. zp2 >= (1.0_dp - eps4) * zp2max ) then
!                if ( deltazj <= epsi ) then

                   inst % intzjmp (i,j,k) = 0.0_dp

                else

                   xb = zp2max / zp2 - 1.0_dp
                   if ( xb <= 0.0_dp ) xb = 1.0e-8_dp
                   beta1 = zm * xb
                   beta2 = ( 1.0_dp - zm ) * xb

                   inst % intzjmp (i,j,k) = intmg ( inst % zjm (i,j,k) , inst % zjp (i,j,k) , beta1 , beta2 )

                end if

!                call wreac_MIL ( inp , thd , taumix , yo2 , zm , zp2 , inst % Ydot )

                inst % Ydot (i,j,k,nrv) = 0.0_dp - sum ( inst % Ydot (i,j,k,1:nrv-1) )

                inst % omega (i,j,k) = 0.0_dp
                do l = 1 , nrv
                   inst % omega (i,j,k) = inst % omega (i,j,k) - &
                                          thd % h0fc (l) * inst % Ydot (i,j,k,l)
                end do
                inst % omega (i,j,k) = max ( 0.0_dp , inst % omega (i,j,k) )

                do l = 1 , nreac
                   inst % W_scal_i (i,j,k,l) = W_scal_i (l)
                end do

             end do
          end do
       end do

    end if ! reaction


    deallocate ( Xa , ha )
    if ( allocated ( delta_2 ) ) deallocate ( delta_2 )


    ! Define the aditionnal variables (A.Techer)
    if ( LES .and. nvv == 4 ) then ! 3 methodes to calculate SGS variance

       do k = sz , ez
          do j = sy , ey
             do i = sx , ex

                rho_i = 1.0_dp / inst % u (i,j,k,1)
                wrk   = inst % u (i,j,k,niv+5) * rho_i

                inst % xiv2 (i,j,k) = inst % u (i,j,k,niv+7) * rho_i - ( wrk * wrk )
                inst % xiv3 (i,j,k) = wrk * ( 1.0_dp - wrk ) - inst % u (i,j,k,niv+9) * rho_i
                inst % xixi (i,j,k) = wrk * wrk

             end do
          end do
       end do

!    else

!       call abort_mpi ('error: module_tools_post.f90 no admited this configuration')

    end if


  end subroutine read_stat


!> \brief Transform statistical files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine transform ( ifile , inp , thd )


    integer (ip) , intent (in)                                     :: ifile !< file number
    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (thd_type) , intent (in)                                  :: thd   !< thermodynamic derived type


    integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip)                                                   :: reclmax , myrec , ok , i , j , k , l
    integer (ip)                                                   :: px , py , pz
    integer (ip)                                                   :: ix  , fx  , iy  , fy  , iz  , fz  , &
                                                                      ix0 , fx0 , iy0 , fy0 , iz0 , fz0
    integer (ip)                                                   :: myntx , mynty , myntz
    integer (ip) , dimension (3)                                   :: mycoords , mydims
    character (len_default)                                        :: number_ascii , rank_ascii , plane_ascii , dir_stat , file_type
    integer (ip)                                                   :: mynv
    real (dp) , dimension (:,:,:) , allocatable                    :: wrk1
    real (dp) , dimension (:,:,:,:) , allocatable                  :: ori_u , trafo_u


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! reading the original fields !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if ( rank == rank_default ) write (*,*) 'transforming the grid...'


    ! allocate the auxiliary vector of conserved variables
    allocate ( ori_u ( 1-ng:ntx+ng , 1-ng:nty+ng , 1-ng:ntz+ng , 1:nv ) , &
               wrk1  ( 1-ng:ntx+ng , 1-ng:nty+ng , 1-ng:ntz+ng )        , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate transform 1')


    write ( number_ascii , format_restart ) ifile
    write ( plane_ascii , format_nplane ) inp % plane_number


    if ( inp % plane_type == 'volume' ) then
       dir_stat  = dir_statvol
       file_type = trim (file_statvol) // '_' // trim (plane_ascii)
    else if ( inp % plane_type == 'XY' ) then
       dir_stat  = dir_statXY
       file_type = trim (file_statXY) // '_' // trim (plane_ascii)
    else if ( inp % plane_type == 'XZ' ) then
       dir_stat  = dir_statXZ
       file_type = trim (file_statXZ) // '_' // trim (plane_ascii)
    else if ( inp % plane_type == 'YZ' ) then
       dir_stat  = dir_statYZ
       file_type = trim (file_statYZ) // '_' // trim (plane_ascii)
    else if ( inp % plane_type == 'restart' ) then
       dir_stat  = dir_restart
       file_type = file_restart
    else
       call abort_mpi ('type of plane not defined')
    end if


    if ( nproc == 1 ) then ! sequential problem


       if (ind_files) then ! reading individual binary files per process


          dims (1) = inp % nprocz ; dims (2) = inp % nprocy ; dims (3) = inp % nprocx
          mycoords (:) = 0
          myrec        = 1

          do pz = 1 , inp % nprocz

             iz = ( mycoords(1)*ntz ) / dims(1) + 1
             fz = ( ( mycoords(1)+1 )*ntz ) / dims(1)

             mycoords (1) = mycoords (1) + 1
             mycoords (2) = 0
             mycoords (3) = 0

             do py = 1 , inp % nprocy

                iy = ( mycoords(2)*nty ) / dims(2) + 1
                fy = ( ( mycoords(2)+1 )*nty ) / dims(2)

                mycoords (2) = mycoords (2) + 1
                mycoords (3) = 0

                do px = 1 , inp % nprocx

                   ix = ( mycoords(3)*ntx ) / dims(3) + 1
                   fx = ( ( mycoords(3)+1 )*ntx ) / dims(3)
                   mycoords (3) = mycoords (3) + 1


                   write ( rank_ascii , format_restart ) myrec - 1
                   reclmax = (fx-ix+1) * (fy-iy+1) * (fz-iz+1) * nv * nbit

                   if ( inp % read_BE ) then

                      open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' //     &
                                                        trim (number_ascii) // '_' // trim (rank_ascii) , &
                                            access    = 'direct'        ,                                 &
                                            recl      = reclmax         ,                                 &
                                            form      = 'unformatted'   ,                                 &
                                            status    = 'old'           ,                                 &
                                            convert   = 'big_endian'    ,                                 &
                                            iostat    =  ok             )

                   else

                      open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' //     &
                                                        trim (number_ascii) // '_' // trim (rank_ascii) , &
                                                        access    = 'direct'        ,                     &
                                                        recl      = reclmax         ,                     &
                                                        form      = 'unformatted'   ,                     &
                                                        status    = 'old'           ,                     &
                                                        iostat    =  ok             )

                   end if
                   if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_stat) // trim (file_type) // '_' //    &
                                                                      trim (number_ascii) // '_' // trim (rank_ascii))


                   read ( unit_restart , rec = 1 )             &
                        (((( ori_u (i,j,k,l) , i = ix , fx ) , &
                                               j = iy , fy ) , &
                                               k = iz , fz ) , &
                                               l = 1  , nv )

                   close (unit_restart)

                   myrec = myrec + 1

                end do
             end do
          end do


       else ! one binary file


          reclmax = ceiling ( 1.0_dp * ntx / inp % nprocx + inp % ghostptx ) * &
                    ceiling ( 1.0_dp * nty / inp % nprocy + inp % ghostpty ) * &
                    ceiling ( 1.0_dp * ntz / inp % nprocz + inp % ghostptz ) * &
                    nv * nbit

          if ( inp % read_BE ) then ! read Big_Endian format files

             open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' // trim (number_ascii) , &
                                   access    = 'direct'        ,                                                   &
                                   recl      = reclmax         ,                                                   &
                                   form      = 'unformatted'   ,                                                   &
                                   status    = 'old'           ,                                                   &
                                   convert   = 'big_endian'    ,                                                   &
                                   iostat    =  ok             )

          else

             open ( unit_restart , file      = trim (dir_stat) // trim (file_type) // '_' // trim (number_ascii) , &
                                   access    = 'direct'        ,                                                   &
                                   recl      = reclmax         ,                                                   &
                                   form      = 'unformatted'   ,                                                   &
                                   status    = 'old'           ,                                                   &
                                   iostat    =  ok             )

          end if
          if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_stat) // trim (file_type) // '_' // trim (number_ascii))

          dims (1) = inp % nprocz ; dims (2) = inp % nprocy ; dims (3) = inp % nprocx
          mycoords (:) = 0
          myrec        = 1

          do pz = 1 , inp % nprocz

             iz = ( mycoords(1)*ntz ) / dims(1) + 1
             fz = ( ( mycoords(1)+1 )*ntz ) / dims(1)

             mycoords (1) = mycoords (1) + 1
             mycoords (2) = 0
             mycoords (3) = 0

             do py = 1 , inp % nprocy

                iy = ( mycoords(2)*nty ) / dims(2) + 1
                fy = ( ( mycoords(2)+1 )*nty ) / dims(2)

                mycoords (2) = mycoords (2) + 1
                mycoords (3) = 0

                do px = 1 , inp % nprocx

                   ix = ( mycoords(3)*ntx ) / dims(3) + 1
                   fx = ( ( mycoords(3)+1 )*ntx ) / dims(3)
                   mycoords (3) = mycoords (3) + 1

                   read ( unit_restart , rec = myrec )         & ! it must start at 1
                        (((( ori_u (i,j,k,l) , i = ix , fx ) , &
                                               j = iy , fy ) , &
                                               k = iz , fz ) , &
                                               l = 1  , nv )

                   myrec = myrec + 1

                end do
             end do
          end do

          close (unit_restart)


       end if


    else ! parallel problem


       call abort_mpi ('no support to transform parallel problems!')


    end if


    ! verifying the fields
    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             if ( ori_u (i,j,k,1) <= 0.0_dp ) then
                write (*,'(A,3I7,1PE16.8)') &
                   'error: density is either zero or negative (i,j,k,ori_rho) = ' , &
                   i,j,k , ori_u (i,j,k,1)
                call abort_mpi (' ')
             end if
          end do
       end do
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! modifying the original fields !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! new variables
    mynv = nrv + npv + nvv
    mynv = mynv + niv


    ! trafo_dims in each direction
    mydims (3) = 6   !inp % nprocx
    mydims (2) = 4   !inp % nprocy
    mydims (1) = 1   !inp % nprocz


    ! total number of points in each direction
    ix0 = 1     ; fx0 = ntx     ; myntx = fx0-ix0+1
!    ix0 = 203   ; fx0 = 1300    ; myntx = fx0-ix0+1 ! A.Techer: JISCF simulation
    iy0 = 1     ; fy0 = nty     ; mynty = fy0-iy0+1
    iz0 = 1     ; fz0 = ntz     ; myntz = fz0-iz0+1
!    iz0 = 86    ; fz0 = 96      ; myntz = fz0-iz0+1


    ! check the choice of mydims
    
    if ( ndim >= 1 .and. mydims (3) > (fx0-ix0+1) / (2*ng+1) ) &
    & call abort_mpi ('error: not enough points OR too many MPI process in x-direction')
    if ( ndim >= 2 .and. mydims (2) > (fy0-iy0+1) / (2*ng+1) ) &
    & call abort_mpi ('error: not enough points OR too many MPI process in y-direction')
    if ( ndim >= 3 .and. mydims (1) > (fz0-iz0+1) / (2*ng+1) ) &
    & call abort_mpi ('error: not enough points OR too many MPI process in z-direction')



    ! modify ori_u

    ! problem at the outlet
    ! do l = 1 , nv
    !    do k = sz , ez
    !       do j = sy , ey
    !          do i = ex-11 , ex
    !             ori_u (i,j,k,l) = ori_u (ex-12,j,k,l)
    !          end do
    !       end do
    !    end do
    ! end do

    ! problem at the bottom
    ! do l = 1 , nv
    !    do k = sz , ez
    !       do j = sy , sy+8
    !          do i = sx , ex
    !             ori_u (i,j,k,l) = ori_u (i,sy+9,k,l)
    !          end do
    !       end do
    !    end do
    ! end do

    ! mixture fraction
    ! call Zmsimp ( inp , thd , ori_u , wrk1 )


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! writing the modified fields: one single binary file !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    reclmax = ceiling ( 1.0_dp * myntx / mydims (3) ) * &
              ceiling ( 1.0_dp * mynty / mydims (2) ) * &
              ceiling ( 1.0_dp * myntz / mydims (1) ) * &
              mynv * nbit


    ! create the new restart file
    if ( inp % write_BE ) then ! write Big_Endian format files

       open ( unit_restart , file    = trim (file_type) // '_' // &
                                       trim (number_ascii)      , &
                             access  = 'direct'                 , &
                             recl    = reclmax                  , &
                             form    = 'unformatted'            , &
                             status  = 'unknown'                , &
                             convert = 'big_endian'             , &
                             iostat  = ok                       )

    else

       open ( unit_restart , file    = trim (file_type) // '_' // &
                                       trim (number_ascii)      , &
                             access  = 'direct'                 , &
                             recl    = reclmax                  , &
                             form    = 'unformatted'            , &
                             status  = 'unknown'                , &
                             iostat  = ok                       )

    end if
    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (file_type) // '_' // trim (number_ascii))


    ! allocate the auxiliary vector of conserved variables
    allocate ( trafo_u ( 1-ng:myntx+ng , 1-ng:mynty+ng , 1-ng:myntz+ng , 1:mynv ) , stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate transform 2')

    ! trafo_u (1:1638,1:716,1:179,1:nv) = ori_u (1:1638,1:716,1:179,1:nv)
    ! do i = 1639,1656
    !    trafo_u (i,1:716,1:179,1:nv) = trafo_u (1638,1:716,1:179,1:nv)
    ! end do
    ! trafo_u (1:1656,1:716,180,1:nv) = trafo_u (1:1656,1:716,179,1:nv)

!     ! A.Techer: JISCF simulation
!     trafo_u (1:1300,1:124,1:264,1:nv) = ori_u (1:1300,1:124,1:264,1:nv)
!     do i = 1301,1330
!        trafo_u (i,1:124,1:264,1:nv) = trafo_u (1300,1:124,1:264,1:nv)
!     end do

    mycoords (:) = 0
    myrec        = 1


    ! loop
    do pz = 1 , mydims (1)

       iz = ( mycoords(1)*myntz ) / mydims(1) + iz0
       fz = ( ( mycoords(1)+1 )*myntz ) / mydims(1) + iz0-1

       mycoords (1) = mycoords (1) + 1
       mycoords (2) = 0
       mycoords (3) = 0

       do py = 1 , mydims (2)

          iy = ( mycoords(2)*mynty ) / mydims(2) + iy0
          fy = ( ( mycoords(2)+1 )*mynty ) / mydims(2) + iy0-1

          mycoords (2) = mycoords (2) + 1
          mycoords (3) = 0

          do px = 1 , mydims (3)

             ix = ( mycoords(3)*myntx ) / mydims(3) + ix0
             fx = ( ( mycoords(3)+1 )*myntx ) / mydims(3) + ix0-1
             mycoords (3) = mycoords (3) + 1


!             ! initialization
!             trafo_u (ix:fx,iy:fy,iz:fz,1:mynv) = 0.0_dp

!             ! 5 first conserved variables
!             trafo_u (ix:fx,iy:fy,iz:fz,1:niv) = ori_u (ix:fx,iy:fy,iz:fz,1:niv)

!             ! Mass fractions variables
!              trafo_u (ix:fx,iy:fy,iz:fz,niv+2) = ori_u (ix:fx,iy:fy,iz:fz,niv+1)
!              trafo_u (ix:fx,iy:fy,iz:fz,niv+4) = ori_u (ix:fx,iy:fy,iz:fz,niv+2)
!              trafo_u (ix:fx,iy:fy,iz:fz,niv+6) = ori_u (ix:fx,iy:fy,iz:fz,niv+3)
!              trafo_u (ix:fx,iy:fy,iz:fz,niv+9) = ori_u (ix:fx,iy:fy,iz:fz,niv+4)
!              trafo_u (ix:fx,iy:fy,iz:fz,niv+10) = ori_u (ix:fx,iy:fy,iz:fz,niv+5)

!              trafo_u (ix:fx,iy:fy,iz:fz,niv+1:niv+6) = ori_u (ix:fx,iy:fy,iz:fz,niv+1:niv+6)

!              trafo_u (ix:fx,iy:fy,iz:fz,niv+7) = ori_u (ix:fx,iy:fy,iz:fz,niv+5) * ori_u (ix:fx,iy:fy,iz:fz,niv+5) &
!                                                / ori_u (ix:fx,iy:fy,iz:fz,1) &
!                                                + ori_u (ix:fx,iy:fy,iz:fz,niv+6)

!              trafo_u (ix:fx,iy:fy,iz:fz,niv+8) = ori_u (ix:fx,iy:fy,iz:fz,niv+5) * ori_u (ix:fx,iy:fy,iz:fz,niv+5) &
!                                                / ori_u (ix:fx,iy:fy,iz:fz,1)

!              trafo_u (ix:fx,iy:fy,iz:fz,niv+9) = ori_u (ix:fx,iy:fy,iz:fz,niv+5) &
!                                                - ori_u (ix:fx,iy:fy,iz:fz,niv+5) * ori_u (ix:fx,iy:fy,iz:fz,niv+5) &
!                                                  / ori_u (ix:fx,iy:fy,iz:fz,1) &
!                                                - ori_u (ix:fx,iy:fy,iz:fz,niv+6)


             ! mixture fraction
             ! trafo_u (ix:fx,iy:fy,iz:fz,mynv) = ori_u (ix:fx,iy:fy,iz:fz,1) * &
             !                                    wrk1 (ix:fx,iy:fy,iz:fz)

             ! write (*,*) myrec , mynv
             ! write (*,*) mydims (:)
             ! write (*,*) mycoords (:)
             ! write (*,*) ix , fx
             ! write (*,*) iy , fy
             ! write (*,*) iz , fz


             write (*,'(10(I5))') myrec , px , py , pz ! , ix,fx,iy,fy,iz,fz

             ! read (*,*)

             write ( unit_restart , rec = myrec )         & ! it must start at 1
                   (((( ori_u (i,j,k,l) , i = ix , fx ) , &
                                          j = iy , fy ) , &
                                          k = iz , fz ) , &
                                          l = 1  , mynv )


!              write ( unit_restart , rec = myrec )           &  it must start at 1
!                    (((( trafo_u (i,j,k,l) , i = ix , fx ) , &
!                                             j = iy , fy ) , &
!                                             k = iz , fz ) , &
!                                             l = 1  , mynv )

             myrec = myrec + 1

          end do
       end do
    end do


    close (unit_restart)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! writing the modified fields: individual binary file per process !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!     ! allocate the auxiliary vector of conserved variables
!     allocate ( trafo_u ( 1-ng:myntx+ng , 1-ng:mynty+ng , 1-ng:myntz+ng , 1:mynv ) , stat = ok )
!     if ( ok > 0 ) call abort_mpi ('error allocate transform 2')

!     ! trafo_u (1:1638,1:716,1:179,1:nv) = ori_u (1:1638,1:716,1:179,1:nv)
!     ! do i = 1639,1656
!     !    trafo_u (i,1:716,1:179,1:nv) = trafo_u (1638,1:716,1:179,1:nv)
!     ! end do
!     ! trafo_u (1:1656,1:716,180,1:nv) = trafo_u (1:1656,1:716,179,1:nv)

!     mycoords (:) = 0
!     myrec        = 1


!     ! loop
!     do pz = 1 , mydims (1)

!        iz = ( mycoords(1)*myntz ) / mydims(1) + iz0
!        fz = ( ( mycoords(1)+1 )*myntz ) / mydims(1) + iz0-1

!        mycoords (1) = mycoords (1) + 1
!        mycoords (2) = 0
!        mycoords (3) = 0

!        do py = 1 , mydims (2)

!           iy = ( mycoords(2)*mynty ) / mydims(2) + iy0
!           fy = ( ( mycoords(2)+1 )*mynty ) / mydims(2) + iy0-1

!           mycoords (2) = mycoords (2) + 1
!           mycoords (3) = 0

!           do px = 1 , mydims (3)

!              ix = ( mycoords(3)*myntx ) / mydims(3) + ix0
!              fx = ( ( mycoords(3)+1 )*myntx ) / mydims(3) + ix0-1
!              mycoords (3) = mycoords (3) + 1


!              ! initialization
!              !trafo_u (ix:fx,iy:fy,iz:fz,1:mynv) = 0.0_dp

!              ! 5 first conserved variables
! !            ! trafo_u (ix:fx,iy:fy,iz:fz,1:nv) = ori_u (ix:fx,iy:fy,iz:fz,1:nv)

!              ! Mass fractions variables
!              ! trafo_u (ix:fx,iy:fy,iz:fz,niv+2) = ori_u (ix:fx,iy:fy,iz:fz,niv+1)
!              ! trafo_u (ix:fx,iy:fy,iz:fz,niv+4) = ori_u (ix:fx,iy:fy,iz:fz,niv+2)
!              ! trafo_u (ix:fx,iy:fy,iz:fz,niv+6) = ori_u (ix:fx,iy:fy,iz:fz,niv+3)
!              ! trafo_u (ix:fx,iy:fy,iz:fz,niv+9) = ori_u (ix:fx,iy:fy,iz:fz,niv+4)

!              ! mixture fraction
!              ! trafo_u (ix:fx,iy:fy,iz:fz,mynv) = ori_u (ix:fx,iy:fy,iz:fz,1) * &
!              !                                    wrk1 (ix:fx,iy:fy,iz:fz)

!              ! write (*,*) myrec , mynv
!              ! write (*,*) mydims (:)
!              ! write (*,*) mycoords (:)
!              ! write (*,*) ix , fx
!              ! write (*,*) iy , fy
!              ! write (*,*) iz , fz


!              write (*,'(10(I5))') myrec , px , py , pz ! , ix,fx,iy,fy,iz,fz


!              write ( rank_ascii , format_restart ) myrec - 1
!              reclmax = (fx-ix+1) * (fy-iy+1) * (fz-iz+1) * mynv * nbit


!              ! create the new restart file
!              if ( inp % write_BE ) then

!                 open ( unit_restart , file    = trim (file_type)    // '_' // &
!                                                 trim (number_ascii) // '_' // &
!                                                 trim (rank_ascii)           , &
!                                       access  = 'direct'                    , &
!                                       recl    = reclmax                     , &
!                                       form    = 'unformatted'               , &
!                                       status  = 'unknown'                   , &
!                                       convert = 'big_endian'                , &
!                                       iostat  = ok                          )

!              else

!                 open ( unit_restart , file    = trim (file_type)    // '_' // &
!                                                 trim (number_ascii) // '_' // &
!                                                 trim (rank_ascii)           , &
!                                       access  = 'direct'                    , &
!                                       recl    = reclmax                     , &
!                                       form    = 'unformatted'               , &
!                                       status  = 'unknown'                   , &
!                                       iostat  = ok                          )

!              end if
!              if ( ok /= 0 ) then
!                 write (*,*) 'error opening ' , trim (file_type) // '_' // trim (number_ascii) // '_' // trim (rank_ascii)
!                 stop
!              end if


!              write ( unit_restart , rec = 1 )             &
!                    (((( ori_u (i,j,k,l) , i = ix , fx ) , &
!                                           j = iy , fy ) , &
!                                           k = iz , fz ) , &
!                                           l = 1  , mynv )


!              ! write ( unit_restart , rec = 1 )               &
!              !       (((( trafo_u (i,j,k,l) , i = ix , fx ) , &
!              !                                j = iy , fy ) , &
!              !                                k = iz , fz ) , &
!              !                                l = 1  , mynv )

!              close (unit_restart)

!              myrec = myrec + 1

!           end do
!        end do
!     end do


    !!!!!!!!!!!!!!!!!!
    ! deleting stuff !
    !!!!!!!!!!!!!!!!!!


    deallocate ( ori_u , wrk1 , trafo_u )


    write (*,*) 'file ' , trim (file_type) // '_' // trim (number_ascii) , ' written'


  end subroutine transform


!> \brief Calculate Schlieren.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine schlieren ( dx_i , dy_i , dz_i , var , svar )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: var  !< variable
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: svar !< Schieleren


    real (dp) , parameter                               :: alpha = 5.0_dp , &
                                                           beta  = 1.0_dp
    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: norm_max , norm_max_tmp
    real (dp) , dimension (:,:,:) , allocatable         :: dvarx , dvary , dvarz


    allocate ( dvarx (sx:ex,sy:ey,sz:ez) , &
               dvary (sx:ex,sy:ey,sz:ez) , &
               dvarz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate schieleren')


    ! communicate always before calculate a derivate
    call comm_one (var)

    if ( ndim == 1 ) then
       call dx ( dx_i , var , dvarx )
       dvary = 0.0_dp
       dvarz = 0.0_dp
    else if ( ndim == 2 ) then
       call dx ( dx_i , var , dvarx )
       call dy ( dy_i , var , dvary )
       dvarz = 0.0_dp
    else if ( ndim == 3 ) then
       call dx ( dx_i , var , dvarx )
       call dy ( dy_i , var , dvary )
       call dz ( dz_i , var , dvarz )
    end if


    norm_max = 0.0_dp
    do k = sz,ez
       do j = sy,ey
          do i = sx,ex

             svar (i,j,k) = sqrt ( dvarx (i,j,k) * dvarx (i,j,k) + &
                                   dvary (i,j,k) * dvary (i,j,k) + &
                                   dvarz (i,j,k) * dvarz (i,j,k) )

             norm_max     = max ( norm_max , svar (i,j,k) )

          end do
       end do
    end do


    ! MPI global norm max
    call mpi_allreduce ( norm_max , norm_max_tmp , 1 , MPI_DOUBLE_PRECISION , &
                         MPI_MAX , MPI_COMM_WORLD , mpicode )
    norm_max = norm_max_tmp


    if ( norm_max /= 0.0_dp ) then
       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                svar (i,j,k) = beta * exp ( - alpha * svar (i,j,k) / norm_max )
             end do
          end do
       end do
    else
!       if ( rank == rank_default ) write (*,*) 'schieleren problem'
       svar (sx:ex,sy:ey,sz:ez) = 1.0_dp
    end if


    deallocate ( dvarx , dvary , dvarz )


  end subroutine schlieren


!> \brief Calculate vorticity magnitude.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine vorticity ( dx_i , dy_i , dz_i , ux , vy , wz , vort )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vort !< vorticity magnitude


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dxvary , dxvarz , &
                                                           dyvarx , dyvarz , &
                                                           dzvarx , dzvary


    allocate ( dxvary (sx:ex,sy:ey,sz:ez) , &
               dxvarz (sx:ex,sy:ey,sz:ez) , &
               dyvarx (sx:ex,sy:ey,sz:ez) , &
               dyvarz (sx:ex,sy:ey,sz:ez) , &
               dzvarx (sx:ex,sy:ey,sz:ez) , &
               dzvary (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate vorticity')


    if ( ndim == 2 ) then ! 2D problem


       call comm_one (ux) ; call comm_one (vy)

       call dx ( dx_i , vy , dxvary )
       call dy ( dy_i , ux , dyvarx )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                vort (i,j,k) = dxvary (i,j,k) - dyvarx (i,j,k)
                vort (i,j,k) = abs ( vort (i,j,k) )
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( dx_i , vy , dxvary )
       call dx ( dx_i , wz , dxvarz )
       call dy ( dy_i , ux , dyvarx )
       call dy ( dy_i , wz , dyvarz )
       call dz ( dz_i , ux , dzvarx )
       call dz ( dz_i , vy , dzvary )


       ! this is the norm of the vorticity (root of the enstrophy)
       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                vort (i,j,k) = ( dyvarz (i,j,k) - dzvary (i,j,k) ) * &
                               ( dyvarz (i,j,k) - dzvary (i,j,k) )

                vort (i,j,k) = vort (i,j,k)                      + &
                             ( dzvarx (i,j,k) - dxvarz (i,j,k) ) * &
                             ( dzvarx (i,j,k) - dxvarz (i,j,k) )

                vort (i,j,k) = vort (i,j,k)                      + &
                             ( dxvary (i,j,k) - dyvarx (i,j,k) ) * &
                             ( dxvary (i,j,k) - dyvarx (i,j,k) )

                vort (i,j,k) = sqrt ( vort (i,j,k) )

             end do
          end do
       end do


    end if


    deallocate ( dxvary , dxvarz , dyvarx , dyvarz , dzvarx , dzvary )


  end subroutine vorticity


!> \brief Calculate vorticity-x.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine vorticity_x ( dy_i , dz_i , vy , wz , vort )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vort !< vorticity magnitude


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dyvarz , dzvary


    allocate ( dyvarz (sx:ex,sy:ey,sz:ez) , &
               dzvary (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate vorticity_x')


    call comm_one (vy) ; call comm_one (wz)

    call dy ( dy_i , wz , dyvarz )
    call dz ( dz_i , vy , dzvary )


    do k = sz,ez
       do j = sy,ey
          do i = sx,ex

             vort (i,j,k) = ( dyvarz (i,j,k) - dzvary (i,j,k) )

          end do
       end do
    end do


    deallocate ( dyvarz , dzvary )


  end subroutine vorticity_x


!> \brief Calculate vorticity-Y.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine vorticity_y ( dx_i , dz_i , ux , wz , vort )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vort !< vorticity magnitude


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dxvarz , dzvarx


    allocate ( dxvarz (sx:ex,sy:ey,sz:ez) , &
               dzvarx (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate vorticity_y')


    call comm_one (ux) ; call comm_one (wz)

    call dx ( dx_i , wz , dxvarz )
    call dz ( dz_i , ux , dzvarx )


    do k = sz,ez
       do j = sy,ey
          do i = sx,ex

             vort (i,j,k) = ( dzvarx (i,j,k) - dxvarz (i,j,k) )

          end do
       end do
    end do


    deallocate ( dxvarz , dzvarx )


  end subroutine vorticity_y


!> \brief Calculate vorticity-Z.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine vorticity_z ( dx_i , dy_i , ux , vy , vort )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vort !< vorticity magnitude


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: dxvary , dyvarx


    allocate ( dxvary (sx:ex,sy:ey,sz:ez) , &
               dyvarx (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate vorticity_z')


    call comm_one (ux) ; call comm_one (vy)

    call dx ( dx_i , vy , dxvary )
    call dy ( dy_i , ux , dyvarx )


    do k = sz,ez
       do j = sy,ey
          do i = sx,ex

             vort (i,j,k) = ( dxvary (i,j,k) - dyvarx (i,j,k) )

          end do
       end do
    end do


    deallocate ( dxvary , dyvarx )


  end subroutine vorticity_z


!> \brief Calculate divergence.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine divergence ( dx_i , dy_i , dz_i , v , div )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dx array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v    !< variable
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: div  !< divergence


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: rho_i
    real (dp) , dimension (:,:,:) , allocatable         :: ux , vy , wz
    real (dp) , dimension (:,:,:) , allocatable         :: dux_varx , dvy_vary , dwz_varz


    allocate ( ux       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz       (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dux_varx (sx:ex,sy:ey,sz:ez)                   , &
               dvy_vary (sx:ex,sy:ey,sz:ez)                   , &
               dwz_varz (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate divergence')


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
             end do
          end do
       end do

       call comm_one (ux) ; call comm_one (vy)

       call dx ( dx_i , ux , dux_varx )
       call dy ( dy_i , vy , dvy_vary )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                div (i,j,k) = dux_varx (i,j,k) + dvy_vary (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                rho_i      = 1.0_dp / v (i,j,k,1)
                ux (i,j,k) = v (i,j,k,2) * rho_i
                vy (i,j,k) = v (i,j,k,3) * rho_i
                wz (i,j,k) = v (i,j,k,4) * rho_i
             end do
          end do
       end do

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)

       call dx ( dx_i , ux , dux_varx )
       call dy ( dy_i , vy , dvy_vary )
       call dz ( dz_i , wz , dwz_varz )


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                div (i,j,k) = dux_varx (i,j,k) + dvy_vary (i,j,k) + dwz_varz (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( ux , vy , wz , dux_varx , dvy_vary , dwz_varz )


  end subroutine divergence


!> \brief Calculate vorticity thickness of the shear layer.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine f_domega ( inp , adi , dy_i , ux , domega )


    type (inp_type) , intent (in)                                        :: inp    !< input derived type
    type (adi_type) , intent (in)                                        :: adi    !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i   !< inverted dy array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux     !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: domega !< thickness

    integer (ip)                                    :: ok , nxz , index , i , k
    real (dp)                                       :: umax , umin , der_max_i , deltaU
    real (dp) , dimension (:) , allocatable         :: wrk_min1 , wrk_min2 , wrk_max1 , wrk_max2


    nxz = (ex-sx+1) * (ez-sz+1)

    allocate ( wrk_min1 (nxz) , &
               wrk_min2 (nxz) , &
               wrk_max1 (nxz) , &
               wrk_max2 (nxz) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate f_domega')


    umax   = inp % umax / adi % u_ref
    umin   = inp % umin / adi % u_ref
    deltaU = ( umax - umin )


    ! communicate always before calculate a derivate
    call comm_one (ux)
    call dy ( dy_i , ux , domega )


    index = 1

    do i = sx,ex
       do k = sz,ez
          wrk_min1 (index) = minval ( domega (i,sy:ey,k) )
          wrk_max1 (index) = maxval ( domega (i,sy:ey,k) )
          index            = index + 1
       end do
    end do


    if ( comm1dy /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk_min1 , wrk_min2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_MIN , comm1dy , mpicode )
       call mpi_allreduce ( wrk_max1 , wrk_max2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_MAX , comm1dy , mpicode )
    else
       wrk_min2 (:) = wrk_min1 (:)
       wrk_max2 (:) = wrk_max1 (:)
    end if


    index = 1

    do i = sx,ex
       do k = sz,ez
          der_max_i          = 1.0_dp / max ( abs (wrk_min2(index)) , abs (wrk_max2(index)) )
          domega (i,sy:ey,k) = deltaU * der_max_i
          index              = index + 1
       end do
    end do


    deallocate ( wrk_min1 , wrk_min2 , wrk_max1 , wrk_max2 )


  end subroutine f_domega


!> \brief Calculate momentum thickness of the shear layer.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine f_dtheta ( inp , adi , y , rho , ux , domega )


    type (inp_type) , intent (in)                                        :: inp    !< input derived type
    type (adi_type) , intent (in)                                        :: adi    !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: y      !< inverted dy array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho    !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux     !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: domega !< thickness

    integer (ip)                                    :: ok , nxz , index , i , j , k , iy , fy
    real (dp)                                       :: deltaU , umax , umin , rhom , denom_i , u
    real (dp) , dimension (:) , allocatable         :: wrk1 , wrk2


    nxz = (ex-sx+1) * (ez-sz+1)

    allocate ( wrk1 (nxz) , &
               wrk2 (nxz) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate f_dtheta')


    iy = max ( ng+1 , sy ) ; fy = min ( nty-ng-1 , ey ) ! is not worth to consider points close to the boundaries


    umax = inp % umax / adi % u_ref
    umin = inp % umin / adi % u_ref
    rhom = inp % rhom / adi % rho_ref


    deltaU  = umax - umin
    denom_i = 1.0_dp / ( rhom * deltaU * deltaU )


    index    = 1
    wrk1 (:) = 0.0_dp

    do i = sx,ex
       do k = sz,ez
          do j = iy,fy

             u = ux (i,j,k) ! u can be lower or higher than the maxima at the end of the domain
                            ! being said that, if u<umin or u>umax at the end of the domain this
                            ! expression will not work as it has been seen in reactive problems

             wrk1 (index) = wrk1 (index) + rho (i,j,k) * &
                          ( umax - u )                 * &
                          ( u - umin )                 * &
                            abs ( y(j+1)-y(j) )

          end do
          index = index + 1
       end do
    end do


    if ( comm1dy /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk1 , wrk2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dy , mpicode )
    else
       wrk2 (:) = wrk1 (:)
    end if


    index = 1

    do i = sx,ex
       do k = sz,ez
          domega (i,sy:ey,k) = wrk2 (index) * denom_i
          index              = index + 1
       end do
    end do


    deallocate ( wrk1 , wrk2 )


  end subroutine f_dtheta


!> \brief Calculate mixture thickness of the shear layer.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine f_dmix ( inp , adi , y , rho , Ya , dmix )


    type (inp_type) , intent (in)                                        :: inp  !< input derived type
    type (adi_type) , intent (in)                                        :: adi  !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: y    !< y-coordinate
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho  !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: Ya   !< species array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: dmix !< thickness


    integer (ip)                                    :: ok , nxz , index , i , j , k , iy , fy
    real (dp)                                       :: rhom , denom_i
    real (dp) , dimension (:) , allocatable         :: wrk1 , wrk2


    rhom    = inp % rhom / adi % rho_ref
    denom_i = 1.0_dp / ( inp % rhom * inp % Yst )


    nxz = (ex-sx+1) * (ez-sz+1)
    allocate ( wrk1 (nxz) , &
               wrk2 (nxz) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate f_dmix')

    iy = max ( ng+1 , sy ) ; fy = min ( nty-ng-1 , ey ) ! is not worth to consider points close to the boundaries

    index    = 1
    wrk1 (:) = 0.0_dp

    do i = sx,ex
       do k = sz,ez
          do j = iy,fy
             wrk1 (index) = wrk1 (index) +             &
                            rho (i,j,k) * Ya (i,j,k) * &
                            abs ( y(j+1)-y(j) )
          end do
          index = index + 1
       end do
    end do


    if ( comm1dy /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk1 , wrk2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dy , mpicode )
    else
       wrk2 (:) = wrk1 (:)
    end if


    index = 1

    do i = sx,ex
       do k = sz,ez
          dmix (i,sy:ey,k) = wrk2 (index) * denom_i
          index            = index + 1
       end do
    end do


    deallocate ( wrk1 , wrk2 )


  end subroutine f_dmix


!> \brief Calculate 99% thickness of the shear layer.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dvis_99 ( inp , adi , y , Ya , d99 )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    type (adi_type) , intent (in)                                        :: adi !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: Ya  !< species array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: d99 !< thickness


    integer (ip)                                    :: ok , nxz , index , i , j , k , iy , fy
    real (dp)                                       :: rhom , denom_i , wrkmin , wrkmax
    real (dp) , dimension (:) , allocatable         :: wrk1 , wrkmin1 , wrkmax1 , &
                                                       wrk2 , wrkmin2 , wrkmax2

    rhom    = inp % rhom / adi % rho_ref
    denom_i = 1.0_dp / ( inp % rhom * inp % Yst )


    nxz = (ex-sx+1) * (ez-sz+1)
    allocate ( wrk1    (nxz) , &
               wrk2    (nxz) , &
               wrkmin1 (nxz) , &
               wrkmin2 (nxz) , &
               wrkmax1 (nxz) , &
               wrkmax2 (nxz) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate d99')

    iy = max ( ng+1 , sy ) ; fy = min ( nty-ng-1 , ey ) ! is not worth to consider points close to the boundaries

    index = 1
    do i = sx,ex
       do k = sz,ez
          wrkmin1 (index) = minval ( Ya (i,iy:fy,k) )
          wrkmax1 (index) = maxval ( Ya (i,iy:fy,k) )
          index = index + 1
       end do
    end do


    if ( comm1dy /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrkmin1 , wrkmin2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_MIN , comm1dy , mpicode )
       call mpi_allreduce ( wrkmax1 , wrkmax2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_MAX , comm1dy , mpicode )
    else
       wrkmin2 (:) = wrkmin1 (:)
       wrkmax2 (:) = wrkmax1 (:)
    end if


    index = 1
    wrk1 (:) = 0.0_dp
    do i = sx,ex
       do k = sz,ez
          wrkmin = max ( 0.01_dp * wrkmin2 (index) , 0.01_dp ) ! for zero values
          wrkmax = 0.99_dp * wrkmax2 (index)
          do j = iy,fy
             if ( Ya (i,j,k) > wrkmin .and. Ya (i,j,k) < wrkmax ) then
                wrk1 (index) = wrk1 (index) + abs ( y(j+1)-y(j) )
             end if
          end do
          index = index + 1
       end do
    end do


    if ( comm1dy /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk1 , wrk2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dy , mpicode )
    else
       wrk2 (:) = wrk1 (:)
    end if


    index = 1
    do i = sx,ex
       do k = sz,ez
          d99 (i,sy:ey,k) = wrk2 (index)
          index           = index + 1
       end do
    end do


    deallocate ( wrk1 , wrkmin1 , wrkmax1 , &
                 wrk2 , wrkmin2 , wrkmax2 )


  end subroutine dvis_99


!> \brief Calculate maximum value of a quantity in y-direction in parallel.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine max_val_Y ( var , maxvar )


    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: var    !< variable
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: maxvar !< maximum variable

    integer (ip)                                    :: ok , nxz , index , i , k
    real (dp) , dimension (:) , allocatable         :: wrk_max1 , wrk_max2


    nxz = (ex-sx+1) * (ez-sz+1)

    allocate ( wrk_max1 (nxz) , &
               wrk_max2 (nxz) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate max_Y')


    index = 1
    do i = sx,ex
       do k = sz,ez
          wrk_max1 (index) = maxval ( var (i,sy:ey,k) )
          index            = index + 1
       end do
    end do


    if ( comm1dy /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk_max1 , wrk_max2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_MAX , comm1dy , mpicode )
    else
       wrk_max2 (:) = wrk_max1 (:)
    end if


    index = 1
    do i = sx,ex
       do k = sz,ez
          maxvar (i,sy:ey,k) = wrk_max2(index)
          index              = index + 1
       end do
    end do


    deallocate ( wrk_max1 , wrk_max2 )


  end subroutine max_val_Y


!> \brief Calculate viscous stress tensor for Newtonian fluids.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine tau_ij ( dx_i , dy_i , dz_i , ux , vy , wz , mu , kpa , tau )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux   !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: vy   !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: wz   !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: mu   !< dynamic viscosity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: kpa  !< bulk viscosity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (inout)     :: tau  !< Reynolds stress tensor


    integer (ip)                                        :: ok , i , j , k
    real (dp) , parameter                               :: twothirds = 2.0_dp / 3.0_dp
    real (dp)                                           :: wrk
    real (dp) , dimension (:,:,:) , allocatable         :: dx_ux , dy_ux , dz_ux , &
                                                           dx_vy , dy_vy , dz_vy , &
                                                           dx_wz , dy_wz , dz_wz


    allocate ( dx_ux (sx:ex,sy:ey,sz:ez) , &
               dy_ux (sx:ex,sy:ey,sz:ez) , &
               dz_ux (sx:ex,sy:ey,sz:ez) , &
               dx_vy (sx:ex,sy:ey,sz:ez) , &
               dy_vy (sx:ex,sy:ey,sz:ez) , &
               dz_vy (sx:ex,sy:ey,sz:ez) , &
               dx_wz (sx:ex,sy:ey,sz:ez) , &
               dy_wz (sx:ex,sy:ey,sz:ez) , &
               dz_wz (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate tau_ij')


    if ( ndim == 2 ) then ! 2D problem


       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                wrk = ( kpa (i,j,k) - twothirds * mu (i,j,k) ) * &
                      ( dx_ux (i,j,k) + dy_vy (i,j,k) )

                tau (i,j,k,1,1) = wrk + mu (i,j,k) * ( dx_ux (i,j,k) + dx_ux (i,j,k) )
                tau (i,j,k,1,2) = mu (i,j,k) * ( dy_ux (i,j,k) + dx_vy (i,j,k) )
                tau (i,j,k,1,3) = 0.0_dp

                tau (i,j,k,2,1) = tau (i,j,k,1,2)
                tau (i,j,k,2,2) = wrk + mu (i,j,k) * ( dy_vy (i,j,k) + dy_vy (i,j,k) )
                tau (i,j,k,2,3) = 0.0_dp

                tau (i,j,k,3,:) = 0.0_dp

             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                wrk = ( kpa (i,j,k) - twothirds * mu (i,j,k) ) * &
                      ( dx_ux (i,j,k) + dy_vy (i,j,k) + dz_wz (i,j,k) )

                tau (i,j,k,1,1) = wrk + mu (i,j,k) * ( dx_ux (i,j,k) + dx_ux (i,j,k) )
                tau (i,j,k,1,2) = mu (i,j,k) * ( dy_ux (i,j,k) + dx_vy (i,j,k) )
                tau (i,j,k,1,3) = mu (i,j,k) * ( dz_ux (i,j,k) + dx_wz (i,j,k) )

                tau (i,j,k,2,1) = tau (i,j,k,1,2)
                tau (i,j,k,2,2) = wrk + mu (i,j,k) * ( dy_vy (i,j,k) + dy_vy (i,j,k) )
                tau (i,j,k,2,3) = mu (i,j,k) * ( dz_vy (i,j,k) + dy_wz (i,j,k) )

                tau (i,j,k,3,1) = tau (i,j,k,1,3)
                tau (i,j,k,3,2) = tau (i,j,k,2,3)
                tau (i,j,k,3,3) = wrk + mu (i,j,k) * ( dz_wz (i,j,k) + dz_wz (i,j,k) )

             end do
          end do
       end do


    end if


    deallocate ( dx_ux , dy_ux , dz_ux , &
                 dx_vy , dy_vy , dz_vy , &
                 dx_wz , dy_wz , dz_wz )


  end subroutine tau_ij


!> \brief Calculate Taylor micro-scale.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine length_taylor ( K_ , nu , eps , lt )


    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: K_  !< square of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: nu  !< kinematic viscosity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: eps !< TKE dissipation
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: lt  !< length


    integer (ip) :: i , j , k


    lt (sx:ex,sy:ey,sz:ez) = 0.0_dp
    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             if ( eps (i,j,k) /= 0.0_dp ) &
                  lt (i,j,k) = sqrt ( abs ( 10.0_dp * nu (i,j,k) * K_ (i,j,k) / eps (i,j,k) ) )
          end do
       end do
    end do


  end subroutine length_taylor


!> \brief Calculate Kolmogorov integral length.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine length_kolmogorov ( nu , eps , lk )


    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: nu  !< kinematic viscosity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: eps !< TKE dissipation
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: lk  !< length


    integer (ip) :: i , j , k


    lk (sx:ex,sy:ey,sz:ez) = 0.0_dp
    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             if ( eps (i,j,k) /= 0.0_dp ) &
                  lk (i,j,k) = ( nu (i,j,k) * nu (i,j,k) * nu (i,j,k) / eps (i,j,k) ) ** 0.25_dp
          end do
       end do
    end do


  end subroutine length_kolmogorov


!> \brief Calculate Reynolds based on Taylor micro-scale.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine reynolds_taylor ( K_ , nu , eps , rt )


    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: K_  !< square of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: nu  !< kinematic viscosity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: eps !< TKE dissipation
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: rt  !< Reynolds number


    integer (ip) :: i , j , k


    rt (sx:ex,sy:ey,sz:ez) = 0.0_dp
    do k = sz,ez
       do j = sy,ey
          do i = sx,ex
             if ( eps (i,j,k) /= 0.0_dp ) &
                  rt (i,j,k) = ( K_ (i,j,k) + K_ (i,j,k) ) * sqrt ( 5.0_dp / ( nu (i,j,k) * eps (i,j,k) ) )
          end do
       end do
    end do


  end subroutine reynolds_taylor


!> \brief Calculate Taylor based on shear layer thickness.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine reynolds ( inp , adi , dy_i , ux , reynolds_ )


    type (inp_type) , intent (in)                                        :: inp       !< input derived type
    type (adi_type) , intent (in)                                        :: adi       !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: ux        !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: reynolds_ !< Reynolds number

    integer (ip)                                    :: ok , nxz , index , i , k
    real (dp)                                       :: der_max_i , deltaU2 , umax , umin , rhom , mum_i
    real (dp) , dimension (:) , allocatable         :: wrk_min1 , wrk_min2 , wrk_max1 , wrk_max2


    nxz = (ex-sx+1) * (ez-sz+1)

    allocate ( wrk_min1 (nxz) , &
               wrk_min2 (nxz) , &
               wrk_max1 (nxz) , &
               wrk_max2 (nxz) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate reynolds')


    umax    = inp % umax / adi % u_ref
    umin    = inp % umin / adi % u_ref
    rhom    = inp % rhom / adi % rho_ref
    mum_i   = 1.0_dp / ( inp % mum / adi % mu_ref )
    deltaU2 = ( umax - umin ) * ( umax - umin )


    ! communicate always before calculate a derivate
    call comm_one (ux)
    call dy ( dy_i , ux , reynolds_ )


    index = 1

    do i = sx,ex
       do k = sz,ez
          wrk_min1 (index) = minval ( reynolds_ (i,sy:ey,k) )
          wrk_max1 (index) = maxval ( reynolds_ (i,sy:ey,k) )
          index            = index + 1
       end do
    end do


    if ( comm1dy /= MPI_COMM_NULL ) then
       call mpi_allreduce ( wrk_min1 , wrk_min2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_MIN , comm1dy , mpicode )
       call mpi_allreduce ( wrk_max1 , wrk_max2 , nxz , MPI_DOUBLE_PRECISION , &
                            MPI_MAX , comm1dy , mpicode )
    else
       wrk_min2 (:) = wrk_min1 (:)
       wrk_max2 (:) = wrk_max1 (:)
    end if


    index = 1

    do i = sx,ex
       do k = sz,ez
          der_max_i             = 1.0_dp / max ( abs (wrk_min2(index)) , abs (wrk_max2(index)) )
          index                 = index + 1
          reynolds_ (i,sy:ey,k) = deltaU2 * der_max_i * rhom * mum_i
       end do
    end do


    deallocate ( wrk_min1 , wrk_min2 , wrk_max1 , wrk_max2 )


  end subroutine reynolds


!> \brief Calculate Q criterion.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine Qcriterion ( dx_i , dy_i , dz_i , v , Q )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v    !< conserved variables array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: Q    !< Q criterion


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: ux    , vy    , wz    , &
                                                           dx_ux , dy_ux , dz_ux , &
                                                           dx_vy , dy_vy , dz_vy , &
                                                           dx_wz , dy_wz , dz_wz


    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_ux (sx:ex,sy:ey,sz:ez) ,                &
               dy_ux (sx:ex,sy:ey,sz:ez) ,                &
               dz_ux (sx:ex,sy:ey,sz:ez) ,                &
               dx_vy (sx:ex,sy:ey,sz:ez) ,                &
               dy_vy (sx:ex,sy:ey,sz:ez) ,                &
               dz_vy (sx:ex,sy:ey,sz:ez) ,                &
               dx_wz (sx:ex,sy:ey,sz:ez) ,                &
               dy_wz (sx:ex,sy:ey,sz:ez) ,                &
               dz_wz (sx:ex,sy:ey,sz:ez) ,                &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate Qcriterion')


    if ( ndim == 2 ) then ! 2D problem

       ux (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,2) / v (sx:ex,sy:ey,sz:ez,1)
       vy (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,3) / v (sx:ex,sy:ey,sz:ez,1)

       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                Q (i,j,k) = dx_ux (i,j,k) * dy_vy (i,j,k) - dx_vy (i,j,k) * dy_ux (i,j,k)

             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       ux (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,2) / v (sx:ex,sy:ey,sz:ez,1)
       vy (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,3) / v (sx:ex,sy:ey,sz:ez,1)
       wz (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,4) / v (sx:ex,sy:ey,sz:ez,1)

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

!                Q (i,j,k) = dx_ux (i,j,k) * dy_vy (i,j,k) - dx_vy (i,j,k) * dy_ux (i,j,k) + &
!                            dx_ux (i,j,k) * dz_wz (i,j,k) - dx_wz (i,j,k) * dz_ux (i,j,k) + &
!                            dy_vy (i,j,k) * dz_wz (i,j,k) - dy_wz (i,j,k) * dz_vy (i,j,k)

                ! Error in previous definition of the Q-criterion. It was defined for incompressible case.
                Q (i,j,k) = - 0.5_dp * ( dx_ux (i,j,k) * dx_ux (i,j,k) + &
                                         dy_vy (i,j,k) * dy_vy (i,j,k) + &
                                         dz_wz (i,j,k) * dz_wz (i,j,k) ) &
                            - dy_ux (i,j,k) * dx_vy (i,j,k)              &
                            - dz_ux (i,j,k) * dx_wz (i,j,k)              &
                            - dz_vy (i,j,k) * dy_wz (i,j,k)

             end do
          end do
       end do

    end if


    deallocate ( ux    , vy    , wz    , &
                 dx_ux , dy_ux , dz_ux , &
                 dx_vy , dy_vy , dz_vy , &
                 dx_wz , dy_wz , dz_wz )


  end subroutine Qcriterion


!> \brief Calculate Lambda-2 criterion.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine lbda2criterion ( dx_i , dy_i , dz_i , v , lbda2 )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i  !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v     !< conserved variables array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: lbda2 !< lambda-2 criterion


    integer (ip)                                        :: ok , i , j , k
    real (dp) , parameter                               :: onethird = 1.0_dp / 3.0_dp , pi = ( acos(-1.0_dp) )
    real (dp)                                           :: a_ , b_ , c_ , d_ , p_ , q_ , sign_ , phi_ , Delta_ , tmp
    real (dp)                                           :: S_ , O_
    real (dp)                                           :: M11 , M12 , M13 , &
                                                           M21 , M22 , M23 , &
                                                           M31 , M32 , M33
    real (dp) , dimension (3)                           :: x , z
    real (dp) , dimension (:,:,:) , allocatable         :: ux    , vy    , wz    , &
                                                           dx_ux , dy_ux , dz_ux , &
                                                           dx_vy , dy_vy , dz_vy , &
                                                           dx_wz , dy_wz , dz_wz


    allocate ( ux (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               vy (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wz (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_ux (sx:ex,sy:ey,sz:ez) ,                &
               dy_ux (sx:ex,sy:ey,sz:ez) ,                &
               dz_ux (sx:ex,sy:ey,sz:ez) ,                &
               dx_vy (sx:ex,sy:ey,sz:ez) ,                &
               dy_vy (sx:ex,sy:ey,sz:ez) ,                &
               dz_vy (sx:ex,sy:ey,sz:ez) ,                &
               dx_wz (sx:ex,sy:ey,sz:ez) ,                &
               dy_wz (sx:ex,sy:ey,sz:ez) ,                &
               dz_wz (sx:ex,sy:ey,sz:ez) ,                &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate lbda2criterion')


    if ( ndim == 2 ) then ! 2D problem

       ux (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,2) / v (sx:ex,sy:ey,sz:ez,1)
       vy (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,3) / v (sx:ex,sy:ey,sz:ez,1)

       call comm_one (ux) ; call comm_one (vy)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                ! I think there is no sense in calculating lbda2-criterion in 2D systems
                ! because there are not 3 roots in the characteristic polynomial
                lbda2 (i,j,k) = 0.0_dp


                ! ! M11
                ! ! S_
                ! S_ = ( dx_ux (i,j,k) + dx_ux (i,j,k) ) * ( dx_ux (i,j,k) + dx_ux (i,j,k) ) + &
                !      ( dy_ux (i,j,k) + dx_vy (i,j,k) ) * ( dx_vy (i,j,k) + dy_ux (i,j,k) )
                ! ! O_
                ! O_ = ( dx_ux (i,j,k) - dx_ux (i,j,k) ) * ( dx_ux (i,j,k) - dx_ux (i,j,k) ) + &
                !      ( dy_ux (i,j,k) - dx_vy (i,j,k) ) * ( dx_vy (i,j,k) - dy_ux (i,j,k) )
                ! M11 = 0.25_dp * ( S_ + O_ )

                ! ! M12
                ! ! S_
                ! S_ = ( dx_ux (i,j,k) + dx_ux (i,j,k) ) * ( dy_ux (i,j,k) + dx_vy (i,j,k) ) + &
                !      ( dy_ux (i,j,k) + dx_vy (i,j,k) ) * ( dy_vy (i,j,k) + dy_vy (i,j,k) )
                ! ! O_
                ! O_ = ( dx_ux (i,j,k) - dx_ux (i,j,k) ) * ( dy_ux (i,j,k) - dx_vy (i,j,k) ) + &
                !      ( dy_ux (i,j,k) - dx_vy (i,j,k) ) * ( dy_vy (i,j,k) - dy_vy (i,j,k) )
                ! M12 = 0.25_dp * ( S_ + O_ )

                ! ! M21
                ! ! S_
                ! S_ = ( dx_vy (i,j,k) + dy_ux (i,j,k) ) * ( dx_ux (i,j,k) + dx_ux (i,j,k) ) + &
                !      ( dy_vy (i,j,k) + dy_vy (i,j,k) ) * ( dx_vy (i,j,k) + dy_ux (i,j,k) )
                ! ! O_
                ! O_ = ( dx_vy (i,j,k) - dy_ux (i,j,k) ) * ( dx_ux (i,j,k) - dx_ux (i,j,k) ) + &
                !      ( dy_vy (i,j,k) - dy_vy (i,j,k) ) * ( dx_vy (i,j,k) - dy_ux (i,j,k) )
                ! M21 = 0.25_dp * ( S_ + O_ )

                ! ! M22
                ! ! S_
                ! S_ = ( dx_vy (i,j,k) + dy_ux (i,j,k) ) * ( dy_ux (i,j,k) + dx_vy (i,j,k) ) + &
                !      ( dy_vy (i,j,k) + dy_vy (i,j,k) ) * ( dy_vy (i,j,k) + dy_vy (i,j,k) )
                ! ! O_
                ! O_ = ( dx_vy (i,j,k) - dy_ux (i,j,k) ) * ( dy_ux (i,j,k) - dx_vy (i,j,k) ) + &
                !      ( dy_vy (i,j,k) - dy_vy (i,j,k) ) * ( dy_vy (i,j,k) - dy_vy (i,j,k) )
                ! M22 = 0.25_dp * ( S_ + O_ )

                ! ! second degree polynomial coefficients: ax^3+bx^2+cx+d
                ! a_ = 1.0_dp
                ! b_ = - ( M11 + M22 )
                ! c_ = M11 * M22 - M12 * M21
                ! d_ =

                ! ! discriminant
                ! Delta_ = 18.0_dp * a_ * b_ * c_ * d_ - &
                !          4.0_dp * b_ * b_ * b_ * d_ +  &
                !          b_ * b_ * c_ * c_ -           &
                !          4.0_dp * a_ * c_ * c_ * c_ -  &
                !          27.0_dp * a_ * a_ * d_ * d_

                ! ! roots
                ! if ( Delta_ > 0.0_dp ) then ! three distinct real roots

                !    p_ = ( 3.0_dp * a_ * c_ - b_ * b_ ) / ( 3.0_dp * a_ * a_ )
                !    q_ = ( 2.0_dp * b_ * b_ * b_ - 9.0_dp * a_ * b_ * c_ + &
                !         27.0_dp * a_ * a_ * d_ ) /                        &
                !       ( 27.0_dp * a_ * a_ * a_ )

                !    if ( q_ > 0.0_dp ) then
                !       sign_ = -1.0_dp
                !    else
                !       sign_ = 1.0_dp
                !    end if

                !    phi_ = sqrt ( ( 27.0_dp * q_ * q_ ) / &
                !         ( - 4.0_dp * p_ * p_ * p_ ) )
                !    phi_ = acos (phi_)

                !    tmp = sign_ * 2.0_dp * sqrt ( -onethird * p_ )

                !    z (1) = tmp * cos ( onethird * ( phi_           ) )
                !    z (2) = tmp * cos ( onethird * ( phi_ + pi      ) )
                !    z (3) = tmp * cos ( onethird * ( phi_ + pi + pi ) )

                !    tmp = - b_ / ( 3.0_dp * a_ )

                !    x (1) = z (1) + tmp
                !    x (2) = z (2) + tmp
                !    x (3) = z (3) + tmp

                !    lbda2 (i,j,k) = x (2)

                ! else

                !    lbda2 (i,j,k) = 0.0_dp

                ! end if

             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       ux (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,2) / v (sx:ex,sy:ey,sz:ez,1)
       vy (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,3) / v (sx:ex,sy:ey,sz:ez,1)
       wz (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,4) / v (sx:ex,sy:ey,sz:ez,1)

       call comm_one (ux) ; call comm_one (vy) ; call comm_one (wz)
       call dx ( dx_i , ux , dx_ux ) ; call dy ( dy_i , ux , dy_ux ) ; call dz ( dz_i , ux , dz_ux )
       call dx ( dx_i , vy , dx_vy ) ; call dy ( dy_i , vy , dy_vy ) ; call dz ( dz_i , vy , dz_vy )
       call dx ( dx_i , wz , dx_wz ) ; call dy ( dy_i , wz , dy_wz ) ; call dz ( dz_i , wz , dz_wz )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                ! M11
                ! S_
                S_ = ( dx_ux (i,j,k) + dx_ux (i,j,k) ) * ( dx_ux (i,j,k) + dx_ux (i,j,k) ) + &
                     ( dy_ux (i,j,k) + dx_vy (i,j,k) ) * ( dx_vy (i,j,k) + dy_ux (i,j,k) ) + &
                     ( dz_ux (i,j,k) + dx_wz (i,j,k) ) * ( dz_ux (i,j,k) + dx_wz (i,j,k) )
                ! O_
                O_ = ( dx_ux (i,j,k) - dx_ux (i,j,k) ) * ( dx_ux (i,j,k) - dx_ux (i,j,k) ) + &
                     ( dy_ux (i,j,k) - dx_vy (i,j,k) ) * ( dx_vy (i,j,k) - dy_ux (i,j,k) ) + &
                     ( dz_ux (i,j,k) - dx_wz (i,j,k) ) * ( dz_ux (i,j,k) - dx_wz (i,j,k) )
                M11 = 0.25_dp * ( S_ + O_ )

                ! M12
                ! S_
                S_ = ( dx_ux (i,j,k) + dx_ux (i,j,k) ) * ( dy_ux (i,j,k) + dx_vy (i,j,k) ) + &
                     ( dy_ux (i,j,k) + dx_vy (i,j,k) ) * ( dy_vy (i,j,k) + dy_vy (i,j,k) ) + &
                     ( dz_ux (i,j,k) + dx_wz (i,j,k) ) * ( dy_wz (i,j,k) + dz_vy (i,j,k) )
                ! O_
                O_ = ( dx_ux (i,j,k) - dx_ux (i,j,k) ) * ( dy_ux (i,j,k) - dx_vy (i,j,k) ) + &
                     ( dy_ux (i,j,k) - dx_vy (i,j,k) ) * ( dy_vy (i,j,k) - dy_vy (i,j,k) ) + &
                     ( dz_ux (i,j,k) - dx_wz (i,j,k) ) * ( dy_wz (i,j,k) - dz_vy (i,j,k) )
                M12 = 0.25_dp * ( S_ + O_ )

                ! M13
                ! S_
                S_ = ( dx_ux (i,j,k) + dx_ux (i,j,k) ) * ( dz_ux (i,j,k) + dx_wz (i,j,k) ) + &
                     ( dy_ux (i,j,k) + dx_vy (i,j,k) ) * ( dz_vy (i,j,k) + dy_wz (i,j,k) ) + &
                     ( dz_ux (i,j,k) + dx_wz (i,j,k) ) * ( dz_wz (i,j,k) + dz_wz (i,j,k) )
                ! O_
                O_ = ( dx_ux (i,j,k) - dx_ux (i,j,k) ) * ( dz_ux (i,j,k) - dx_wz (i,j,k) ) + &
                     ( dy_ux (i,j,k) - dx_vy (i,j,k) ) * ( dz_vy (i,j,k) - dy_wz (i,j,k) ) + &
                     ( dz_ux (i,j,k) - dx_wz (i,j,k) ) * ( dz_wz (i,j,k) - dz_wz (i,j,k) )
                M13 = 0.25_dp * ( S_ + O_ )

                ! M21
                ! S_
                S_ = ( dx_vy (i,j,k) + dy_ux (i,j,k) ) * ( dx_ux (i,j,k) + dx_ux (i,j,k) ) + &
                     ( dy_vy (i,j,k) + dy_vy (i,j,k) ) * ( dx_vy (i,j,k) + dy_ux (i,j,k) ) + &
                     ( dz_vy (i,j,k) + dy_wz (i,j,k) ) * ( dx_wz (i,j,k) + dz_ux (i,j,k) )
                ! O_
                O_ = ( dx_vy (i,j,k) - dy_ux (i,j,k) ) * ( dx_ux (i,j,k) - dx_ux (i,j,k) ) + &
                     ( dy_vy (i,j,k) - dy_vy (i,j,k) ) * ( dx_vy (i,j,k) - dy_ux (i,j,k) ) + &
                     ( dz_vy (i,j,k) - dy_wz (i,j,k) ) * ( dx_wz (i,j,k) - dz_ux (i,j,k) )
                M21 = 0.25_dp * ( S_ + O_ )

                ! M22
                ! S_
                S_ = ( dx_vy (i,j,k) + dy_ux (i,j,k) ) * ( dy_ux (i,j,k) + dx_vy (i,j,k) ) + &
                     ( dy_vy (i,j,k) + dy_vy (i,j,k) ) * ( dy_vy (i,j,k) + dy_vy (i,j,k) ) + &
                     ( dz_vy (i,j,k) + dy_wz (i,j,k) ) * ( dy_wz (i,j,k) + dz_vy (i,j,k) )
                ! O_
                O_ = ( dx_vy (i,j,k) - dy_ux (i,j,k) ) * ( dy_ux (i,j,k) - dx_vy (i,j,k) ) + &
                     ( dy_vy (i,j,k) - dy_vy (i,j,k) ) * ( dy_vy (i,j,k) - dy_vy (i,j,k) ) + &
                     ( dz_vy (i,j,k) - dy_wz (i,j,k) ) * ( dy_wz (i,j,k) - dz_vy (i,j,k) )
                M22 = 0.25_dp * ( S_ + O_ )

                ! M23
                ! S_
                S_ = ( dx_vy (i,j,k) + dy_ux (i,j,k) ) * ( dz_ux (i,j,k) + dx_wz (i,j,k) ) + &
                     ( dy_vy (i,j,k) + dy_vy (i,j,k) ) * ( dz_vy (i,j,k) + dy_wz (i,j,k) ) + &
                     ( dz_vy (i,j,k) + dy_wz (i,j,k) ) * ( dz_wz (i,j,k) + dz_wz (i,j,k) )
                ! O_
                O_ = ( dx_vy (i,j,k) - dy_ux (i,j,k) ) * ( dz_ux (i,j,k) - dx_wz (i,j,k) ) + &
                     ( dy_vy (i,j,k) - dy_vy (i,j,k) ) * ( dz_vy (i,j,k) - dy_wz (i,j,k) ) + &
                     ( dz_vy (i,j,k) - dy_wz (i,j,k) ) * ( dz_wz (i,j,k) - dz_wz (i,j,k) )
                M23 = 0.25_dp * ( S_ + O_ )

                ! M31
                ! S_
                S_ = ( dx_wz (i,j,k) + dz_ux (i,j,k) ) * ( dx_ux (i,j,k) + dx_ux (i,j,k) ) + &
                     ( dy_wz (i,j,k) + dz_vy (i,j,k) ) * ( dx_vy (i,j,k) + dy_ux (i,j,k) ) + &
                     ( dz_wz (i,j,k) + dz_wz (i,j,k) ) * ( dx_wz (i,j,k) + dz_ux (i,j,k) )
                ! O_
                O_ = ( dx_wz (i,j,k) - dz_ux (i,j,k) ) * ( dx_ux (i,j,k) - dx_ux (i,j,k) ) + &
                     ( dy_wz (i,j,k) - dz_vy (i,j,k) ) * ( dx_vy (i,j,k) - dy_ux (i,j,k) ) + &
                     ( dz_wz (i,j,k) - dz_wz (i,j,k) ) * ( dx_wz (i,j,k) - dz_ux (i,j,k) )
                M31 = 0.25_dp * ( S_ + O_ )

                ! M32
                ! S_
                S_ = ( dx_wz (i,j,k) + dz_ux (i,j,k) ) * ( dy_ux (i,j,k) + dx_vy (i,j,k) ) + &
                     ( dy_wz (i,j,k) + dz_vy (i,j,k) ) * ( dy_vy (i,j,k) + dy_vy (i,j,k) ) + &
                     ( dz_wz (i,j,k) + dz_wz (i,j,k) ) * ( dy_wz (i,j,k) + dz_vy (i,j,k) )
                ! O_
                O_ = ( dx_wz (i,j,k) - dz_ux (i,j,k) ) * ( dy_ux (i,j,k) - dx_vy (i,j,k) ) + &
                     ( dy_wz (i,j,k) - dz_vy (i,j,k) ) * ( dy_vy (i,j,k) - dy_vy (i,j,k) ) + &
                     ( dz_wz (i,j,k) - dz_wz (i,j,k) ) * ( dy_wz (i,j,k) - dz_vy (i,j,k) )
                M32 = 0.25_dp * ( S_ + O_ )

                ! M33
                ! S_
                S_ = ( dx_wz (i,j,k) + dz_ux (i,j,k) ) * ( dz_ux (i,j,k) + dx_wz (i,j,k) ) + &
                     ( dy_wz (i,j,k) + dz_vy (i,j,k) ) * ( dz_vy (i,j,k) + dy_wz (i,j,k) ) + &
                     ( dz_wz (i,j,k) + dz_wz (i,j,k) ) * ( dz_wz (i,j,k) + dz_wz (i,j,k) )
                ! O_
                O_ = ( dx_wz (i,j,k) - dz_ux (i,j,k) ) * ( dz_ux (i,j,k) - dx_wz (i,j,k) ) + &
                     ( dy_wz (i,j,k) - dz_vy (i,j,k) ) * ( dz_vy (i,j,k) - dy_wz (i,j,k) ) + &
                     ( dz_wz (i,j,k) - dz_wz (i,j,k) ) * ( dz_wz (i,j,k) - dz_wz (i,j,k) )
                M33 = 0.25_dp * ( S_ + O_ )


                ! third degree polynomial coefficients: ax^3+bx^2+cx+d
                a_ = 1.0_dp
                b_ = - ( M11 + M22 + M33 )
                c_ = M11 * M22 + M11 * M33 + M22 * M33 - &
                   ( M12 * M21 + M13 * M31 + M23 * M32 )
                d_ = M13 * ( M22 * M31 - M21 * M32 ) + &
                     M23 * ( M11 * M32 - M12 * M31 ) + &
                     M33 * ( M12 * M21 - M11 * M22 )

                ! discriminant
                Delta_ = 18.0_dp * a_ * b_ * c_ * d_ - &
                         4.0_dp * b_ * b_ * b_ * d_ +  &
                         b_ * b_ * c_ * c_ -           &
                         4.0_dp * a_ * c_ * c_ * c_ -  &
                         27.0_dp * a_ * a_ * d_ * d_

                ! roots
                if ( Delta_ > 0.0_dp ) then ! three distinct real roots

                   p_ = ( 3.0_dp * a_ * c_ - b_ * b_ ) / ( 3.0_dp * a_ * a_ )
                   q_ = ( 2.0_dp * b_ * b_ * b_ - 9.0_dp * a_ * b_ * c_ + &
                        27.0_dp * a_ * a_ * d_ ) /                        &
                      ( 27.0_dp * a_ * a_ * a_ )

                   if ( q_ > 0.0_dp ) then
                      sign_ = -1.0_dp
                   else
                      sign_ = 1.0_dp
                   end if

                   phi_ = sqrt ( ( 27.0_dp * q_ * q_ ) / &
                        ( - 4.0_dp * p_ * p_ * p_ ) )
                   phi_ = acos (phi_)

                   tmp = sign_ * 2.0_dp * sqrt ( -onethird * p_ )

                   z (1) = tmp * cos ( onethird * ( phi_           ) )
                   z (2) = tmp * cos ( onethird * ( phi_ + pi      ) )
                   z (3) = tmp * cos ( onethird * ( phi_ + pi + pi ) )

                   tmp = - b_ / ( 3.0_dp * a_ )

                   x (1) = z (1) + tmp
                   x (2) = z (2) + tmp
                   x (3) = z (3) + tmp

                   lbda2 (i,j,k) = x (2)

                else

                   lbda2 (i,j,k) = 0.0_dp

                end if

             end do
          end do
       end do

    end if


    deallocate ( ux    , vy    , wz    , &
                 dx_ux , dy_ux , dz_ux , &
                 dx_vy , dy_vy , dz_vy , &
                 dx_wz , dy_wz , dz_wz )


  end subroutine lbda2criterion


!> \brief Calculate Takeno index.
!!
!! * Diffusion, null, takeno=0.
!! * Premixed, unity, takeno=1.
!! * Inert, half, takeno=0.5.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine takeno ( inp , dx_i , dy_i , dz_i , omega , v , tk )


    type (inp_type) , intent (in)                                       :: inp   !< input derived type
    real (dp) , dimension (:) , allocatable , intent (in)               :: dx_i  !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)               :: dy_i  !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)               :: dz_i  !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)           :: omega !< chemical production rate
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)         :: v     !< conserved variables array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)        :: tk    !< Takeno index


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: norm_Yf , norm_Yo , wrk
    real (dp) , dimension (:,:,:) , allocatable         :: Yf    , Yo            , &
                                                           dx_Yf , dy_Yf , dz_Yf , &
                                                           dx_Yo , dy_Yo , dz_Yo


    allocate ( Yf (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               Yo (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               dx_Yf (sx:ex,sy:ey,sz:ez) ,                &
               dy_Yf (sx:ex,sy:ey,sz:ez) ,                &
               dz_Yf (sx:ex,sy:ey,sz:ez) ,                &
               dx_Yo (sx:ex,sy:ey,sz:ez) ,                &
               dy_Yo (sx:ex,sy:ey,sz:ez) ,                &
               dz_Yo (sx:ex,sy:ey,sz:ez) ,                &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate takeno')


    if (reaction) then ! reactive problems


       if ( ndim == 2 ) then ! 2D problem

          Yf (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,niv+inp % index_fuel) / &
                                   v (sx:ex,sy:ey,sz:ez,1)
          Yo (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,niv+inp % index_oxidizer) / &
                                   v (sx:ex,sy:ey,sz:ez,1)

          call comm_one (Yf) ; call comm_one (Yo)
          call dx ( dx_i , Yf , dx_Yf ) ; call dy ( dy_i , Yf , dy_Yf )
          call dx ( dx_i , Yo , dx_Yo ) ; call dy ( dy_i , Yo , dy_Yo )

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex

                   if ( Yf (i,j,k) > minYf .and. Yo (i,j,k) > minYf .and. omega (i,j,k) > minomega ) then

                      norm_Yf = sqrt ( dx_Yf (i,j,k) * dx_Yf (i,j,k) + dy_Yf (i,j,k) * dy_Yf (i,j,k) )
                      norm_Yo = sqrt ( dx_Yo (i,j,k) * dx_Yo (i,j,k) + dy_Yo (i,j,k) * dy_Yo (i,j,k) )

                      wrk = max ( norm_Yf * norm_Yo , epsi )

                      tk (i,j,k) = 0.5_dp * ( 1.0_dp + ( dx_Yf (i,j,k) * dx_Yo (i,j,k) +         &
                                                         dy_Yf (i,j,k) * dy_Yo (i,j,k) ) / wrk )

                   else

                      tk (i,j,k) = 0.5_dp

                   end if

                end do
             end do
          end do

       else if ( ndim == 3 ) then ! 3D problem

          Yf (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,niv+inp % index_fuel) / &
                                   v (sx:ex,sy:ey,sz:ez,1)
          Yo (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,niv+inp % index_oxidizer) / &
                                   v (sx:ex,sy:ey,sz:ez,1)

          call comm_one (Yf) ; call comm_one (Yo)
          call dx ( dx_i , Yf , dx_Yf ) ; call dy ( dy_i , Yf , dy_Yf ) ; call dz ( dz_i , Yf , dz_Yf )
          call dx ( dx_i , Yo , dx_Yo ) ; call dy ( dy_i , Yo , dy_Yo ) ; call dz ( dz_i , Yo , dz_Yo )

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex

                   if ( Yf (i,j,k) > minYf .and. Yo (i,j,k) > minYf .and. omega (i,j,k) > minomega ) then

                      norm_Yf = sqrt ( dx_Yf (i,j,k) * dx_Yf (i,j,k) + dy_Yf (i,j,k) * &
                                       dy_Yf (i,j,k) + dz_Yf (i,j,k) * dz_Yf (i,j,k) )
                      norm_Yo = sqrt ( dx_Yo (i,j,k) * dx_Yo (i,j,k) + dy_Yo (i,j,k) * &
                                       dy_Yo (i,j,k) + dz_Yo (i,j,k) * dz_Yo (i,j,k) )

                      wrk = max ( norm_Yf * norm_Yo , epsi )

                      tk (i,j,k) = 0.5_dp * ( 1.0_dp + ( dx_Yf (i,j,k) * dx_Yo (i,j,k) +         &
                                                         dy_Yf (i,j,k) * dy_Yo (i,j,k) +         &
                                                         dz_Yf (i,j,k) * dz_Yo (i,j,k) ) / wrk )

                   else

                      tk (i,j,k) = 0.5_dp

                   end if

                end do
             end do
          end do

       end if


    else ! inert problems


       if ( ndim == 2 ) then ! 2D problem

          Yf (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,niv+inp % index_fuel) / &
                                   v (sx:ex,sy:ey,sz:ez,1)
          Yo (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,niv+inp % index_oxidizer) / &
                                   v (sx:ex,sy:ey,sz:ez,1)

          call comm_one (Yf) ; call comm_one (Yo)
          call dx ( dx_i , Yf , dx_Yf ) ; call dy ( dy_i , Yf , dy_Yf )
          call dx ( dx_i , Yo , dx_Yo ) ; call dy ( dy_i , Yo , dy_Yo )

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex

                   if ( Yf (i,j,k) > minYf .and. Yo (i,j,k) > minYf ) then

                      norm_Yf = sqrt ( dx_Yf (i,j,k) * dx_Yf (i,j,k) + dy_Yf (i,j,k) * dy_Yf (i,j,k) )
                      norm_Yo = sqrt ( dx_Yo (i,j,k) * dx_Yo (i,j,k) + dy_Yo (i,j,k) * dy_Yo (i,j,k) )

                      wrk = max ( norm_Yf * norm_Yo , epsi )

                      tk (i,j,k) = 0.5_dp * ( 1.0_dp + ( dx_Yf (i,j,k) * dx_Yo (i,j,k) +         &
                                                         dy_Yf (i,j,k) * dy_Yo (i,j,k) ) / wrk )

                   else

                      tk (i,j,k) = 0.5_dp

                   end if

                end do
             end do
          end do

       else if ( ndim == 3 ) then ! 3D problem

          Yf (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,niv+inp % index_fuel) / &
                                   v (sx:ex,sy:ey,sz:ez,1)
          Yo (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,niv+inp % index_oxidizer) / &
                                   v (sx:ex,sy:ey,sz:ez,1)

          call comm_one (Yf) ; call comm_one (Yo)
          call dx ( dx_i , Yf , dx_Yf ) ; call dy ( dy_i , Yf , dy_Yf ) ; call dz ( dz_i , Yf , dz_Yf )
          call dx ( dx_i , Yo , dx_Yo ) ; call dy ( dy_i , Yo , dy_Yo ) ; call dz ( dz_i , Yo , dz_Yo )

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex

                   if ( Yf (i,j,k) > minYf .and. Yo (i,j,k) > minYf ) then

                      norm_Yf = sqrt ( dx_Yf (i,j,k) * dx_Yf (i,j,k) + dy_Yf (i,j,k) * &
                                       dy_Yf (i,j,k) + dz_Yf (i,j,k) * dz_Yf (i,j,k) )
                      norm_Yo = sqrt ( dx_Yo (i,j,k) * dx_Yo (i,j,k) + dy_Yo (i,j,k) * &
                                       dy_Yo (i,j,k) + dz_Yo (i,j,k) * dz_Yo (i,j,k) )

                      wrk = max ( norm_Yf * norm_Yo , epsi )

                      tk (i,j,k) = 0.5_dp * ( 1.0_dp + ( dx_Yf (i,j,k) * dx_Yo (i,j,k) +         &
                                                         dy_Yf (i,j,k) * dy_Yo (i,j,k) +         &
                                                         dz_Yf (i,j,k) * dz_Yo (i,j,k) ) / wrk )

                   else

                      tk (i,j,k) = 0.5_dp

                   end if

                end do
             end do
          end do

       end if


    end if


    deallocate ( Yf    , Yo            , &
                 dx_Yf , dy_Yf , dz_Yf , &
                 dx_Yo , dy_Yo , dz_Yo )


  end subroutine takeno


!> \brief Shock detector.
!!
!! This subroutine is used only in post-treatment.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine shock_det_post ( v , T , W_i , crit )


    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v    !< conserved variables array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: T    !< temperature
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: W_i  !< inverted molar mass
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: crit !< criteria


    integer (ip) , parameter                  :: stencil_m1 = ng+ng
    integer (ip)                              :: st
    integer (ip)                              :: i , i_l , i_r , i_s , &
                                                 j , j_l , j_r , j_s , &
                                                 k , k_l , k_r , k_s
    real (dp)                                 :: drho , dpres , p_p , p_m


    crit (sx:ex,sy:ey,sz:ez) = 0.0_dp


    if ( ndim >= 1 ) then
       ! X-direction criteria
       do k = sz , ez ! loop in the z-direction
          do j = sy , ey ! loop in the y-direction

             do i_l = sx-1 , ex ! loop on the cell faces

                i_r = i_l + 1

                do st = 1 , stencil_m1

                   i_s = i_l + st - ng

                   ! density criteria
                   drho = abs ( v (min(i_s+1,ex),j,k,1) - v (i_s,j,k,1) ) / &
                          min ( v (min(i_s+1,ex),j,k,1) , v (i_s,j,k,1) )

                   ! pressure criteria
                   p_p   = v (min(i_s+1,ex),j,k,1) * T (min(i_s+1,ex),j,k) * W_i (min(i_s+1,ex),j,k)
                   p_m   = v (i_s,j,k,1) * T (i_s,j,k) * W_i (i_s,j,k)
                   dpres = abs ( p_p - p_m ) / min ( p_p , p_m )

                   if ( drho > max_rel_weight .and. dpres > max_rel_weight ) crit (i_l,j,k) = 1.0_dp

                end do

                ! activate 2D WENO at the boundaries (not 3D because of periodicity)
                if ( i_l < 1+ng )   crit (i_l,j,k) = 1.0_dp
                if ( i_l > ntx-ng ) crit (i_l,j,k) = 1.0_dp
                if ( j   < 1+ng )   crit (i_l,j,k) = 1.0_dp
                if ( j   > nty-ng ) crit (i_l,j,k) = 1.0_dp

             end do

          end do
       end do
    end if


    if ( ndim >= 2 ) then
       ! Y-direction criteria
       do k = sz , ez ! loop in the z-direction
          do i = sx , ex ! loop in the x-direction

             do j_l = sy-1 , ey ! loop on the cell faces

                j_r = j_l + 1

                do st = 1 , stencil_m1

                   j_s = j_l + st - ng

                   ! density criteria
                   drho = abs ( v (i,min(j_s+1,ey),k,1) - v (i,j_s,k,1) ) / &
                          min ( v (i,min(j_s+1,ey),k,1) , v (i,j_s,k,1) )

                   ! pressure criteria
                   p_p   = v (i,min(j_s+1,ey),k,1) * T (i,min(j_s+1,ey),k) * W_i (i,min(j_s+1,ey),k)
                   p_m   = v (i,j_s,k,1) * T (i,j_s,k) * W_i (i,j_s,k)
                   dpres = abs ( p_p - p_m ) / min ( p_p , p_m )

                   if ( drho > max_rel_weight .and. dpres > max_rel_weight ) crit (i,j_l,k) = 2.0_dp

                end do

                ! activate 2D WENO at the boundaries (not 3D because of periodicity)
                if ( i   < 1+ng )   crit (i,j_l,k) = 2.0_dp
                if ( i   > ntx-ng ) crit (i,j_l,k) = 2.0_dp
                if ( j_l < 1+ng )   crit (i,j_l,k) = 2.0_dp
                if ( j_l > nty-ng ) crit (i,j_l,k) = 2.0_dp

             end do

          end do
       end do
    end if


    if ( ndim == 3 ) then
       ! Z-direction criteria
       do j = sy , ey ! loop in the y-direction
          do i = sx , ex ! loop in the x-direction

             do k_l = sz-1 , ez ! loop on the cell faces

                k_r = k_l + 1

                do st = 1 , stencil_m1

                   k_s = k_l + st - ng

                   ! density criteria
                   drho = abs ( v (i,j,min(k_s+1,ez),1) - v (i,j,k_s,1) ) / &
                          min ( v (i,j,min(k_s+1,ez),1) , v (i,j,k_s,1) )

                   ! pressure criteria
                   p_p = v (i,j,min(k_s+1,ez),1) * T (i,j,min(k_s+1,ez)) * W_i (i,j,min(k_s+1,ez))
                   p_m = v (i,j,k_s,1) * T (i,j,k_s) * W_i (i,j,k_s)
                   dpres = abs ( p_p - p_m ) / min ( p_p , p_m )

                   if ( drho > max_rel_weight .and. dpres > max_rel_weight ) crit (i,j_l,k) = 3.0_dp

                end do

                ! activate 2D WENO at the boundaries (not 3D because of periodicity)
!                if ( i < 1+ng )   crit (i,j,k_l) = 3.0_dp
!                if ( i > ntx-ng ) crit (i,j,k_l) = 3.0_dp
!                if ( j < 1+ng )   crit (i,j,k_l) = 3.0_dp
!                if ( j > nty-ng ) crit (i,j,k_l) = 3.0_dp

             end do

          end do
       end do
    end if


  end subroutine shock_det_post


!> \brief Reaction detector.
!!
!! New criteria that seems to work well with the presence or not of
!! radicals:
!!
!! #- The fuel stream is composed by hydrogen and nitrogen. If the
!!    inlet quantities are changed we are not in the fuel stream but
!!    in a mixture region.
!! #- The oxidizer stream is always composed by a certain amount of
!!    nitrogen and zero hydrogen.
!!
!! If these two quantities change, then we are not in the oxidizer
!! stream.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine reaction_det ( ifile , inp , adi , sim , v , Ydot , T , crit )


    integer (ip) , intent (in)                                           :: ifile !< file number
    type (inp_type) , intent (in)                                        :: inp   !< input derived type
    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    type (sim_type) , intent (in)                                        :: sim   !< simulation derived type
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v     !< variable
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: Ydot  !< production rates array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: T     !< temperature
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: crit  !< criteria


    real (dp) , parameter         :: eps_w = 1.0e-12_dp , eps_f = 1.0e-6_dp , eps_o = 1.0e-6_dp
    integer (ip)                  :: ix , fx , &
                                     iy , fy , &
                                     iz , fz
    integer (ip)                  :: i , j , k , l
    real (dp)                     :: dt_phys , rho_phys , T_phys , rho_i , rho0 , T0_phys , Tref_i , rho_cgs , wrk
    logical                       :: call_dvode , call_dvode_w , call_dvode_f , call_dvode_o
    real (dp) , dimension (1:nrv) :: Ya , Yf , Yo


    crit (sx:ex,sy:ey,sz:ez) = 0.0_dp


    Tref_i  = 1.0_dp / adi % T_ref
    rho_cgs = adi % rho_ref * 1.0e-3_dp


    ! specify the pure fuel
    Yf (1:nrv) = inp % Y1f (1:nrv)
!    call X_to_Y ( thd , Yf )

    ! specify the pure oxidizer
    Yo (1:nrv) = inp % Y0o (1:nrv)
!    call X_to_Y ( thd , Yo )


    ! domain to apply combustion (not in every boundary)
    ! if ( ndim >= 2 ) then
    !    ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 ) ! added ng
    !    iy = max ( ng+1 , sy ) ; fy = min ( ey , nty-ng-1 ) ! added ng
    !    iz = sz                ; fz = ez
    ! else
    !    ix = max ( ng+1 , sx ) ; fx = min ( ex , ntx-ng-1 ) ! added ng
    !    iy = sy                ; fy = ey
    !    iz = sz                ; fz = ez
    ! end if
    ix = sx ; fx = ex
    iy = sy ; fy = ey
    iz = sz ; fz = ez


    if ( sim % dtime_per_ite (ifile) >= 1.0_dp ) then
       dt_phys = 1.0e-9_dp ! S.I. units
    else
       dt_phys = sim % dtime_per_ite (ifile) ! S.I. units
    end if
    if ( dt_phys == 0.0_dp ) dt_phys = 1.0e-9_dp ! S.I. units
    dt_phys = dt_phys / adi % time_ref ! because Ydot is the "non dimensional" production rate dY/dt


    if (reaction) then

       do k = iz , fz
          do j = iy , fy
             do i = ix , fx


                call_dvode   = .false.
                call_dvode_w = .false.
                call_dvode_f = .false.
                call_dvode_o = .false.


                rho_i = 1.0_dp / v (i,j,k,1)
                do l = 1 , nrv
                   Ya (l)  = v (i,j,k, niv+l ) * rho_i
                end do


                rho0     = v (i,j,k,1)
                rho_phys = rho0 * rho_cgs
                T0_phys  = T (i,j,k) * adi % T_ref
                T_phys   = T0_phys


                ! verify the production rate of each species
                do l = 1 , nrv
                   wrk = Ydot (i,j,k,l) * dt_phys ! wrk is the mass fraction variation: Yf-Yi
                   if ( abs (wrk) > eps_w ) call_dvode_w = .true.
                end do


                ! verify the pure fuel
                ! do l = 1 , nrv
                !    if ( abs ( Ya(l)-Yf(l) ) > eps_f ) call_dvode_f = .true.
                ! end do

                if ( abs ( Ya (inp % index_H2) - Yf (inp % index_H2) ) > eps_f .and. &
                     abs ( Ya (inp % index_N2) - Yf (inp % index_N2) ) > eps_f ) call_dvode_f = .true.


                ! verify the pure oxidizer
                ! do l = 1 , nrv
                !    if ( abs ( Ya(l)-Yo(l) ) > eps2 ) call_dvode_o = .true.
                ! end do
                if ( abs ( Ya (inp % index_H2) - Yo (inp % index_H2) ) > eps_o .and. &
                     abs ( Ya (inp % index_N2) - Yo (inp % index_N2) ) > eps_o ) call_dvode_o = .true.


                if ( call_dvode_w )       call_dvode = .true.
                if ( .not. call_dvode_f ) call_dvode = .false.
                if ( .not. call_dvode_o ) call_dvode = .false.


                if (call_dvode) crit (i,j,k) = 1.0_dp


             end do
          end do
       end do

    end if


  end subroutine reaction_det


!> \brief x-direction similarity average.
!!
!! This subroutine only works if nx is small, i.e. YZ planes.
!! Delta x is supposed to be constant.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine similarity_x ( inp , v )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)         :: v   !< conserved variables array


    integer (ip)                                :: i , j , k
    real (dp)                                   :: wrk1 , wrk2


    if ( nxmax > 11 .or. .not. inp % similarity_x ) return


    wrk1 = real ( ex-sx+1 , dp )
    wrk1 = 1.0_dp / wrk1

    do k = sz , ez
       do j = sy , ey
          wrk2 = 0.0_dp
          do i = sx , ex
             wrk2 = wrk2 + v (i,j,k)
          end do
          wrk2 = wrk2 * wrk1
          v (sx:ex,j,k) = wrk2
       end do
    end do


  end subroutine similarity_x


!> \brief x-direction for 2D-bins variables.
!!
!! This subroutine only works if nx is small, i.e. YZ planes.
!! Delta x is supposed to be constant.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine similarity_x_2Dbins ( inp , v )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)         :: v   !< conserved variables array


    integer (ip)                                :: i , j , n
    real (dp)                                   :: wrk1 , wrk2


    if ( nxmax > 11 .or. .not. inp % similarity_x ) return


    wrk1 = real ( ex-sx+1 , dp )
    wrk1 = 1.0_dp / wrk1

    do n = 1 , nbins
       do j = sy , ey
          wrk2 = 0.0_dp
          do i = sx , ex
             wrk2 = wrk2 + v (i,j,n)
          end do
          wrk2 = wrk2 * wrk1
          v (sx:ex,j,n) = wrk2
       end do
    end do


  end subroutine similarity_x_2Dbins


!> \brief x-direction for 1D-bins variables.
!!
!! This subroutine only works if nx is small, i.e. YZ planes.
!! Delta x is supposed to be constant.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine similarity_x_1Dbins ( inp , v )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    real (dp) , allocatable , dimension (:,:) , intent (inout)           :: v   !< conserved variables array


    integer (ip)                                :: i , n
    real (dp)                                   :: wrk1 , wrk2


    if ( nxmax > 11 .or. .not. inp % similarity_x ) return


    wrk1 = real ( ex-sx+1 , dp )
    wrk1 = 1.0_dp / wrk1

    do n = 1 , nbins
       wrk2 = 0.0_dp
       do i = sx , ex
          wrk2 = wrk2 + v (i,n)
       end do
       wrk2 = wrk2 * wrk1
       v (sx:ex,n) = wrk2
    end do


  end subroutine similarity_x_1Dbins


!> \brief Open convergence files.
!!
!! For using the probes, to estimate the convergence of the simulation
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine open_convergence ( inp , adi , x , y , z )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    type (adi_type) , intent (in)                                        :: adi !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


    integer (ip)                        :: i , j , k , px , py , pz
    integer (ip)                        :: currunit , iprobe , ok
    real (dp)                           :: L_ref , L_ref_nodim
    logical                             :: cond
    real (dp) , dimension (ndimmax)     :: pt
    character (len_default) , parameter :: format = ' ( A30 , 1X , 3 ( 1X , F16.8 ) ) '


    px = ex ; py = ey ; pz = ez ! initialize this to avoid bugs (verified)
    L_ref = adi % L_ref ; L_ref_nodim = adi % L_ref / inp % length_ref


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true.

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true.  ! A.Techer pour etre coherent avec "locate_probes"

       if (cond) then

          currunit = unit_conv + iprobe - 1

          open ( unit   = currunit ,                                    &
                 file   = trim (dir_plot) // trim (file_conv) // '_' // &
                          trim (inp % probe_name (iprobe)) // '.out' ,  &
                 form   = 'formatted'    ,                              &
                 status = 'unknown'      ,                              &
                 iostat = ok )
          if ( ok /= 0 ) call abort_mpi &
               ('error opening ' // trim (dir_plot) // trim (file_conv) // trim (inp % probe_name (iprobe)))

          ! calculate points
          do i = ex+1 , sx , -1
             if ( x (i) > pt (1) ) px = i
          end do
          px = min(px,ex) ; px = max(px,sx)
          do j = ey+1 , sy , -1
             if ( y (j) > pt (2) ) py = j
          end do
          py = min(py,ey) ; py = max(py,sy)
          if ( ndim == 3 ) then
             do k = ez+1 , sz , -1
                if ( z (k) > pt (3) ) pz = k
             end do
          else
             pz = sz
          end if
          pz = min(pz,ez) ; pz = max(pz,sz) ! A.Techer

          write ( currunit , * )      '# probe:  ' , trim (inp % probe_name (iprobe))
          write ( currunit , * )      '# points:  ' , px , py , pz
          write ( currunit , format ) '# coordinates (dim):  ' , x (px) * L_ref , y (py) * L_ref , z (pz) * L_ref
          write ( currunit , format ) '# coordinates (nodim):  ' , x (px) * L_ref_nodim , y (py) * L_ref_nodim , &
                                                                   z (pz) * L_ref_nodim
          write ( currunit , '(A)' ) '# variables:  file_number,time,ra_p,ra_Y,ra_Zm,ra_ux,ra_vy,ra_wz,' // &
                                     'var_uxux,var_vyvy,var_wzwz,var_uxvy,var_uxwz,var_vywz,'   // &
                                     'var_uxvywz,f(uxux),f(vyvy),f(wzwz),f(uxvy),f(uxwz),f(vywz),enstrophy'

       end if

    end do


  end subroutine open_convergence


!> \brief Close convergence files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine close_convergence ( inp , x , y , z )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


    integer (ip)                       :: iprobe
    logical                            :: cond
    real (dp) , dimension (ndimmax)    :: pt


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true.

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

       if (cond) close ( unit = unit_conv + iprobe - 1 )

    end do


  end subroutine close_convergence


!> \brief Plot convergence files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine plot_convergence ( ifile , inp , adi , sim , x , y , z , time , dtime , stat )


    integer (ip) , intent (in)                                           :: ifile !< file number
    type (inp_type) , intent (in)                                        :: inp   !< input derived type
    type (adi_type) , intent (in)                                        :: adi   !< non-dimensional derived type
    type (sim_type) , intent (in)                                        :: sim   !< simulation derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x     !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y     !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z     !< z-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: time  !< time
    real (dp) , dimension (:) , allocatable , intent (in)                :: dtime !< time step
    type (stat_type) , intent (in)                                       :: stat  !< statistical derived type


    integer (ip)                        :: i , j , k , px , py , pz
    integer (ip)                        :: iprobe , currunit
    real (dp)                           :: sumdt_i , P_ref , u_ref , u2_ref , u3_ref , ens_ref
    real (dp) , dimension (ndimmax)     :: pt
    real (dp) , dimension (6)           :: wrk
    logical                             :: cond
    character (len_default) , parameter :: format = ' ( I4 , 2X , 30 ( 1X , E16.8 ) ) '


    px = ex ; py = ey ; pz = ez ! initialize this to avoid bugs (verified)
    p_ref   = adi % P_ref
    u_ref   = adi % u_ref
    u2_ref  = u_ref * u_ref
    u3_ref  = u_ref * u2_ref
    ens_ref = u_ref * u_ref / ( adi % L_ref * adi % L_ref )


    if ( ifile == sim % correct_endfile ) then
       sumdt_i = 1.0_dp
    else
       sumdt_i = 0.0_dp
       do i = inp % start_file , ifile , inp % skip_file
          sumdt_i = sumdt_i + dtime (i)
       end do
       sumdt_i = 1.0_dp / sumdt_i
    end if


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true.

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

       if (cond) then

          ! calculate points
          do i = ex+1 , sx , -1
             if ( x (i) > pt (1) ) px = i
          end do
          do j = ey+1 , sy , -1
             if ( y (j) > pt (2) ) py = j
          end do
          if ( ndim == 3 ) then
             do k = ez+1 , sz , -1
                if ( z (k) > pt (3) ) pz = k
             end do
          else
             pz = sz
          end if

          ! wrk variables must tend to zero

          wrk (1) = stat % reyavg_uxux (px,py,pz) * sumdt_i - &
                    stat % reyavg_ux (px,py,pz) * sumdt_i *   &
                    stat % reyavg_ux (px,py,pz) * sumdt_i -   &
                    stat % reyvar_ux (px,py,pz) * sumdt_i

          wrk (2) = stat % reyavg_vyvy (px,py,pz) * sumdt_i - &
                    stat % reyavg_vy (px,py,pz) * sumdt_i *   &
                    stat % reyavg_vy (px,py,pz) * sumdt_i -   &
                    stat % reyvar_vy (px,py,pz) * sumdt_i

          wrk (3) = stat % reyavg_wzwz (px,py,pz) * sumdt_i - &
                    stat % reyavg_wz (px,py,pz) * sumdt_i *   &
                    stat % reyavg_wz (px,py,pz) * sumdt_i -   &
                    stat % reyvar_wz (px,py,pz) * sumdt_i

          wrk (4) = stat % reyavg_uxvy (px,py,pz) * sumdt_i - &
                    stat % reyavg_ux (px,py,pz) * sumdt_i *   &
                    stat % reyavg_vy (px,py,pz) * sumdt_i -   &
                    stat % reyvar_uxvy (px,py,pz) * sumdt_i

          wrk (5) = stat % reyavg_uxwz (px,py,pz) * sumdt_i - &
                    stat % reyavg_ux (px,py,pz) * sumdt_i *   &
                    stat % reyavg_wz (px,py,pz) * sumdt_i -   &
                    stat % reyvar_uxwz (px,py,pz) * sumdt_i

          wrk (6) = stat % reyavg_vywz (px,py,pz) * sumdt_i - &
                    stat % reyavg_vy (px,py,pz) * sumdt_i *   &
                    stat % reyavg_wz (px,py,pz) * sumdt_i -   &
                    stat % reyvar_vywz (px,py,pz) * sumdt_i

          currunit = unit_conv + iprobe - 1

          write ( currunit , format ) ifile ,                                                   &
                                      time (ifile) * adi % time_ref ,                           &
                                      stat % reyavg_P (px,py,pz)      * sumdt_i * P_ref       , &
                                      stat % reyavg_Y (px,py,pz)      * sumdt_i               , &
                                      stat % reyavg_Zm (px,py,pz)     * sumdt_i               , &
                                      stat % reyavg_ux (px,py,pz)     * sumdt_i * u_ref       , &
                                      stat % reyavg_vy (px,py,pz)     * sumdt_i * u_ref       , &
                                      stat % reyavg_wz (px,py,pz)     * sumdt_i * u_ref       , &
                                      stat % reyvar_ux (px,py,pz)     * sumdt_i * u2_ref      , &
                                      stat % reyvar_vy (px,py,pz)     * sumdt_i * u2_ref      , &
                                      stat % reyvar_wz (px,py,pz)     * sumdt_i * u2_ref      , &
                                      stat % reyvar_uxvy (px,py,pz)   * sumdt_i * u2_ref      , &
                                      stat % reyvar_uxwz (px,py,pz)   * sumdt_i * u2_ref      , &
                                      stat % reyvar_vywz (px,py,pz)   * sumdt_i * u2_ref      , &
                                      stat % reyvar_uxvywz (px,py,pz) * sumdt_i * u3_ref      , &
                                      wrk (1) * u2_ref                                        , &
                                      wrk (2) * u2_ref                                        , &
                                      wrk (3) * u2_ref                                        , &
                                      wrk (4) * u2_ref                                        , &
                                      wrk (5) * u2_ref                                        , &
                                      wrk (6) * u2_ref                                        , &
                                      stat % reyavg_enstrophy (px,py,pz) * sumdt_i * ens_ref
       end if

    end do


  end subroutine plot_convergence


!> \brief Open probability density function files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine open_pdfs ( inp , adi , x , y , z )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    type (adi_type) , intent (in)                                        :: adi !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


    integer (ip)                        :: i , j , k , px , py , pz
    integer (ip)                        :: currunit , iprobe , ok
    real (dp)                           :: L_ref , L_ref_nodim
    logical                             :: cond
    real (dp) , dimension (ndimmax)     :: pt
    character (len_default) , parameter :: format = ' ( A30 , 1X , 3 ( 1X , F16.8 ) ) '


    px = ex ; py = ey ; pz = ez ! initialize this to avoid bugs (verified)
    L_ref = adi % L_ref ; L_ref_nodim = adi % L_ref / inp % length_ref


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true.

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

       if (cond) then

          currunit = unit_pdf + iprobe - 1

          open ( unit   = currunit ,                                   &
                 file   = trim (dir_plot) // trim (file_pdf) // '_' // &
                          trim (inp % probe_name (iprobe)) // '.out' , &
                 form   = 'formatted'    ,                             &
                 status = 'unknown'      ,                             &
                 iostat = ok )
          if ( ok /= 0 ) call abort_mpi &
               ('error opening ' // trim (dir_plot) // trim (file_pdf) // trim (inp % probe_name (iprobe)))

          ! calculate points
          do i = ex+1 , sx , -1
             if ( x (i) > pt (1) ) px = i
          end do
          px = min(px,ex) ; px = max(px,sx)
          do j = ey+1 , sy , -1
             if ( y (j) > pt (2) ) py = j
          end do
          py = min(py,ey) ; py = max(py,sy)
          if ( ndim == 3 ) then
             do k = ez+1 , sz , -1
                if ( z (k) > pt (3) ) pz = k
             end do
          else
             pz = sz
          end if
          pz = min(pz,ez) ; pz = max(pz,sz) ! A.Techer

          write ( currunit , * )      '# probe:  ' , trim (inp % probe_name (iprobe))
          write ( currunit , * )      '# points:  ' , px , py! , pz
          write ( currunit , format ) '# coordinates (dim):  '   , x (px) * L_ref , y (py) * L_ref! , z (pz) * L_ref
          write ( currunit , format ) '# coordinates (nodim):  ' , x (px) * L_ref_nodim , y (py) * L_ref_nodim! , z (pz) * L_ref_nodim
          write ( currunit , '(A)' ) '# variables:  Y,pdf(Y),Zm,pdf(Zm),sdr,pdf(sdr),Da,pdf(Da),FI,pdf(FI),' // &
                                     'FI,cpdf(omega,FI),Da,cpdf(omega,Da),Y,pdfbeta(Y),sdr,gauss(sdr),Da,gauss(Da)'

       end if

    end do


  end subroutine open_pdfs


!> \brief Close probability density function files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine close_pdfs ( inp , x , y , z )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


    integer (ip)                       :: iprobe
    logical                            :: cond
    real (dp) , dimension (ndimmax)    :: pt


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true.

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

       if (cond) close ( unit = unit_pdf + iprobe - 1 )

    end do


  end subroutine close_pdfs


!> \brief Plot probability density function files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine plot_pdfs ( inp , adi , x , y , z , favavg , stat )


    type (inp_type) , intent (in)                                        :: inp    !< input derived type
    type (adi_type) , intent (in)                                        :: adi    !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x      !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y      !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z      !< z-coordinate array
    type (favavg_type) , intent (inout)                                  :: favavg !< Favre averaged derived type
    type (stat_type) , intent (inout)                                    :: stat   !< statistical derived type


    logical                             :: cond
    integer (ip)                        :: n , i , j , k , px , py , pz
    integer (ip)                        :: iprobe , currunit
    real (dp)                           :: L_ref , sdr_ref , Da_ref
    real (dp)                           :: varbin , mean , variance , wrk , a_ , b_ , pdfbeta_Y , gauss_sdr , gauss_Da
    real (dp) , dimension (ndimmax)     :: pt

    character (len_default) , parameter :: format = ' ( 30 ( 1X , E16.8 ) ) '


    L_ref = adi % L_ref
    if ( inp % dis_rate_therm ) then
       sdr_ref = ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / &
                 ( adi % L_ref * adi % L_ref )
    else
       sdr_ref = adi % D_ref / ( adi % L_ref * adi % L_ref )
    end if
    Da_ref = ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )


    call similarity_x_2Dbins ( inp , stat % pdf_Y )
    call similarity_x_2Dbins ( inp , stat % pdf_Zm )
    call similarity_x_2Dbins ( inp , stat % pdf_sdr )
    call similarity_x_2Dbins ( inp , stat % pdf_Da )
    call similarity_x_2Dbins ( inp , stat % pdf_FI )

    call similarity_x_1Dbins ( inp , stat % cpdf_omega_FI )
    call similarity_x_1Dbins ( inp , stat % cpdf_omega_Da )

    call similarity_x ( inp , favavg % Da )
    call similarity_x ( inp , favavg % sdr )
    call similarity_x ( inp , stat % favvar_Da )
    call similarity_x ( inp , stat % favvar_sdr )


    px = ex ; py = ey ; pz = ez ! initialize this to avoid bugs (verified)
    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true. ! not writing from ghost points because of communications problems

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

       if (cond) then

          if ( iprobe == 1 ) write (*,*) 'writing pdfs files ...'

          ! calculate points
          do i = ex+1 , sx , -1
             if ( x (i) > pt (1) ) px = i
          end do
          do j = ey+1 , sy , -1
             if ( y (j) > pt (2) ) py = j
          end do
          if ( ndim == 3 ) then
             do k = ez+1 , sz , -1
                if ( z (k) > pt (3) ) pz = k
             end do
          else
             pz = sz
          end if

          currunit = unit_pdf + iprobe - 1

          do n = 1,nbins


             ! pdfbeta for Y (mixture fraction)
             varbin = min_Y + (n-1) * (max_Y-min_Y) / (nbins-1)
             if ( inp % index_scalar == 0 ) then
                mean     = favavg % Zm (px,py,pz)
                variance = stat % favvar_Zm  (px,py,pz)
             else
                mean     = favavg % Y (px,py,pz, inp % index_scalar )
                variance = stat % favvar_Y  (px,py,pz, inp % index_scalar )
             end if
             wrk = mean * ( 1.0_dp - mean ) / variance - 1.0_dp ! A.Techer: error, add "- 1.0_dp" in wrk variable..
             a_  = mean * wrk
             b_  = ( 1.0_dp - mean ) * wrk
             call pdfbeta ( varbin , a_ , b_ , pdfbeta_Y )

             ! Gaussian distribution for sdr
             varbin   = log10 (min_sdr) + (n-1) * (log10(max_sdr)-log10(min_sdr) ) / (nbins-1)
             varbin   = 10.0_dp ** varbin
             mean     = favavg % sdr (px,py,pz) * sdr_ref
             variance = stat % favvar_sdr (px,py,pz) * sdr_ref * sdr_ref
!             write (*,*) 'sdr', mean, variance
!             call gaussian ( varbin , mean , variance , gauss_sdr )
             call log_normal ( varbin , mean , variance , gauss_sdr )

             ! Gaussian distribution for Da
             varbin   = log10 (min_Da) + (n-1) * (log10(max_Da)-log10(min_Da) ) / (nbins-1)
             varbin   = 10.0_dp ** varbin
             mean     = favavg % Da (px,py,pz) * Da_ref
             variance = stat % favvar_Da (px,py,pz) * Da_ref * Da_ref
!             write (*,*) 'Da', mean, variance
!             call gaussian ( varbin , mean , variance , gauss_Da )
             call log_normal ( varbin , mean , variance , gauss_Da )


             write ( currunit , format ) min_Y + (n-1) * (max_Y-min_Y) / (nbins-1)                   , & ! 01
                                         stat % pdf_Y (px,py,n)                                      , & ! 02

                                         min_Y + (n-1) * (max_Y-min_Y) / (nbins-1)                   , & ! 03
                                         stat % pdf_Zm (px,py,n)                                     , & ! 04

             ( log10 (min_sdr) + (n-1) * (log10(max_sdr)-log10(min_sdr) ) / (nbins-1) )              , & ! 05
                                         stat % pdf_sdr (px,py,n) / sdr_ref                          , & ! 06

                ( log10 (min_Da) + (n-1) * (log10(max_Da)-log10(min_Da) ) / (nbins-1) )              , & ! 07
                                         stat % pdf_Da (px,py,n) / Da_ref                            , & ! 08

                                         min_FI + (n-1) * (max_FI-min_FI) / (nbins-1)                , & ! 09
                                         stat % pdf_FI (px,py,n)                                     , & ! 10

                                         min_FI + (n-1) * (max_FI-min_FI) / (nbins-1)                , & ! 11
                                         stat % cpdf_omega_FI (px,n)                                 , & ! 12

                ( log10 (min_Da) + (n-1) * (log10(max_Da)-log10(min_Da) ) / (nbins-1) )              , & ! 13
                                         stat % cpdf_omega_Da (px,n) / Da_ref                        , & ! 14

                                         min_Y + (n-1) * (max_Y-min_Y) / (nbins-1)                   , & ! 15
                                         pdfbeta_Y                                                   , & ! 16

!                                         min_sdr + (n-1) * (max_sdr-min_sdr) / (nbins-1)             , & ! 17
             ( log10 (min_sdr) + (n-1) * (log10(max_sdr)-log10(min_sdr) ) / (nbins-1) )              , & ! 17
                                         gauss_sdr                                                   , & ! 18

!                                         min_Da + (n-1) * (max_Da-min_Da) / (nbins-1)                , & ! 19
                ( log10 (min_Da) + (n-1) * (log10(max_Da)-log10(min_Da) ) / (nbins-1) )              , & ! 19
                                         gauss_Da                                                        ! 20


          end do

       end if

    end do


  end subroutine plot_pdfs


!> \brief Open spectra files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine open_spectras ( inp , adi , x , y , z )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    type (adi_type) , intent (in)                                        :: adi !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


    integer (ip)                        :: i , j , k , px , py , pz
    integer (ip)                        :: currunit , ok
    real (dp)                           :: L_ref , L_ref_nodim
    logical                             :: cond
    real (dp) , dimension (ndimmax)     :: pt
    character (len_default) , parameter :: format = ' ( A30 , 1X , 3 ( 1X , F16.8 ) ) '


    px = ex ; py = ey ; pz = ez ! initialize this to avoid bugs (verified)
    L_ref = adi % L_ref ; L_ref_nodim = adi % L_ref / inp % length_ref


    pt (1) = inp % corrspec_coord (1)
    pt (2) = inp % corrspec_coord (2)
    pt (3) = inp % corrspec_coord (3)

    cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true.

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

    if (cond) then

       currunit = unit_spec

       open ( unit   = currunit ,                                      &
              file   = trim (dir_plot) // trim (file_spec) // '.out' , &
              form   = 'formatted'    ,                                &
              status = 'unknown'      ,                                &
              iostat = ok )
       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_plot) // trim (file_spec))

       ! calculate points
       do i = ex+1 , sx , -1
          if ( x (i) > pt (1) ) px = i
       end do
       px = min(px,ex) ; px = max(px,sx)
       do j = ey+1 , sy , -1
          if ( y (j) > pt (2) ) py = j
       end do
       py = min(py,ey) ; py = max(py,sy)
       if ( ndim == 3 ) then
          do k = ez+1 , sz , -1
             if ( z (k) > pt (3) ) pz = k
          end do
       else
          pz = sz
       end if
       pz = min(pz,ez) ; pz = max(pz,sz) ! A.Techer


       write ( currunit , * )      '# points:  ' , px , py! , pz
       write ( currunit , format ) '# coordinates (dim):  ' , x (px) * L_ref , y (py) * L_ref! , z (pz) * L_ref
       write ( currunit , format ) '# coordinates (nodim):  ' , x (px) * L_ref_nodim , y (py) * L_ref_nodim! , z (pz) * L_ref_nodim
       write ( currunit , '(A)' )  '# variables:  t,f,spec_ux,spec_vy,spec_wz,spec_Y,spec_Zm,spec_sdr'


    end if


  end subroutine open_spectras


!> \brief Close spectra files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine close_spectras ( inp , x , y , z )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


    logical                            :: cond
    real (dp) , dimension (ndimmax)    :: pt


    pt (1) = inp % corrspec_coord (1)
    pt (2) = inp % corrspec_coord (2)
    pt (3) = inp % corrspec_coord (3)

    cond = .false.

!    if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!         pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!         pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true. ! taking ghost points because the point could be between two blocks

    if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
         pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
         pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

    if (cond) close ( unit = unit_spec )


  end subroutine close_spectras


!> \brief Plot spectra files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine plot_spectras ( inp , adi , sim , x , y , z , t , stat )


    type (inp_type) , intent (in)                                        :: inp  !< input derived type
    type (adi_type) , intent (in)                                        :: adi  !< non-dimensional derived type
    type (sim_type) , intent (in)                                        :: sim  !< simulation derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x    !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y    !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z    !< z-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: t    !< time array
    type (stat_type) , intent (inout)                                    :: stat !< statistical derived type


    logical                                        :: cond
    character (len_default) , parameter            :: format = ' ( 30 ( 1X , E16.8 ) ) '
    integer (ip)                                   :: ok , n , k , px , py , pz
    integer (ip)                                   :: start_file , end_file
    integer (ip)                                   :: currunit
    real (dp)                                      :: Tperiod
    real (dp)                                      :: v_ref , sdr_ref
    real (dp) , dimension (ndimmax)                :: pt
    real (dp) , dimension (:,:,:) , allocatable    :: wrk1


    allocate  ( wrk1 ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng ) , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate plot_spectras')


    start_file = inp % start_file
    end_file   = sim % correct_endfile !inp % end_file


    Tperiod = sim % sumdtime * adi % time_ref
    v_ref   = adi % u_ref ; v_ref = v_ref * v_ref
    if ( inp % dis_rate_therm ) then
       sdr_ref = ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / &
                 ( adi % L_ref * adi % L_ref )
    else
       sdr_ref = adi % D_ref / ( adi % L_ref * adi % L_ref )
    end if
    sdr_ref = sdr_ref * sdr_ref


    px = ex ; py = ey ; pz = ez ! initialize this to avoid bugs (verified)

    pt (1) = inp % corrspec_coord (1)
    pt (2) = inp % corrspec_coord (2)
    pt (3) = inp % corrspec_coord (3)

    cond = .false.

!    if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!         pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!         pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true.

    if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
         pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
         pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

    if (cond) then

       write (*,*) 'writing spectras file ...'

       currunit = unit_spec

       k = 0
       do n = start_file , ( end_file + start_file ) / 2

!          write (*,*) n , k , 1.0_dp / Tperiod , real ( k , dp ) / Tperiod

          write ( currunit , format ) t (n) * adi % time_ref              , &
                                      real ( k , dp ) / Tperiod           , & ! hertzs
                                      stat % spec_ux  (sz,n) * v_ref      , &
                                      stat % spec_vy  (sz,n) * v_ref      , &
                                      stat % spec_wz  (sz,n) * v_ref      , &
                                      stat % spec_Y   (sz,n)              , &
                                      stat % spec_Zm  (sz,n)              , &
                                      stat % spec_sdr (sz,n) * sdr_ref
          k = k+1
       end do

    end if


  end subroutine plot_spectras


!> \brief Open conditional average files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine open_cond_avgs ( inp , adi , x , y , z )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    type (adi_type) , intent (in)                                        :: adi !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


    integer (ip)                        :: i , j , k , px , py , pz
    integer (ip)                        :: currunit , iprobe , ok
    real (dp)                           :: L_ref , L_ref_nodim
    logical                             :: cond
    real (dp) , dimension (ndimmax)     :: pt
    character (len_default) , parameter :: format = ' ( A10 , 1X , 3 ( 1X , F16.8 ) ) '


    px = ex ; py = ey ; pz = ez ! initialize this to avoid bugs (verified)
    L_ref = adi % L_ref ; L_ref_nodim = adi % L_ref / inp % length_ref


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true.

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

       if (cond) then

          currunit = unit_cavg + iprobe - 1

          open ( unit   = currunit ,                                    &
                 file   = trim (dir_plot) // trim (file_cavg) // '_' // &
                          trim (inp % probe_name (iprobe)) // '.out' ,  &
                 form   = 'formatted'    ,                              &
                 status = 'unknown'      ,                              &
                 iostat = ok )
          if ( ok /= 0 ) call abort_mpi &
               ('error opening ' // trim (dir_plot) // trim (file_conv) // trim (inp % probe_name (iprobe)))

          ! calculate points
          do i = ex+1 , sx , -1
             if ( x (i) > pt (1) ) px = i
          end do
          px = min(px,ex) ; px = max(px,sx)
          do j = ey+1 , sy , -1
             if ( y (j) > pt (2) ) py = j
          end do
          py = min(py,ey) ; py = max(py,sy)
          if ( ndim == 3 ) then
             do k = ez+1 , sz , -1
                if ( z (k) > pt (3) ) pz = k
             end do
          else
             pz = sz
          end if
          pz = min(pz,ez) ; pz = max(pz,sz) ! A.Techer

          write ( currunit , * )      'probe' , trim (inp % probe_name (iprobe))
          write ( currunit , * )      'points' , px! , py , pz
          write ( currunit , format ) 'coordinates (dim)' , x (px) * L_ref! , y (py) * L_ref , z (pz) * L_ref
          write ( currunit , format ) 'coordinates (nodim)' , x (px) * L_ref_nodim! , y (py) * L_ref_nodim , z (pz) * L_ref_nodim
          write ( currunit , '(A)' )  '# variables:  Zm,ca_T_Zm,cv_T_Zm,Zm,ca_sdr_Zm,cv_sdr_Zm,Zm,ca_Da_Zm,' // &
                                      'cv_Da_Zm,Zm,ca_omega_Zm,cv_omega_Zm,Zm,ca_FI_Zm,cv_FI_Zm,FI,ca_omega_FI,' // &
                                      'cv_omega_FI,Da,ca_omega_Da,cv_omega_Da'

       end if

    end do


  end subroutine open_cond_avgs


!> \brief Close conditional average files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine close_cond_avgs ( inp , x , y , z )


    type (inp_type) , intent (in)                                        :: inp !< input derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x   !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y   !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z   !< z-coordinate array


    integer (ip)                       :: iprobe
    logical                            :: cond
    real (dp) , dimension (ndimmax)    :: pt


    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true. ! taking ghost points because the point could be between two blocks

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

       if (cond) close ( unit = unit_cavg + iprobe - 1 )

    end do


  end subroutine close_cond_avgs


!> \brief Plot conditional average files.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine plot_cond_avgs ( inp , adi , x , y , z , favavg , stat )


    type (inp_type) , intent (in)                                        :: inp    !< input derived type
    type (adi_type) , intent (in)                                        :: adi    !< non-dimensional derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: x      !< x-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: y      !< y-coordinate array
    real (dp) , dimension (:) , allocatable , intent (in)                :: z      !< z-coordinate array
    type (favavg_type) , intent (inout)                                  :: favavg !< Favre averaged derived type
    type (stat_type) , intent (inout)                                    :: stat   !< statistical derived type


    logical                             :: cond
    integer (ip)                        :: n , i , j , k , px , py , pz
    integer (ip)                        :: iprobe , currunit
    real (dp)                           :: sdr_ref , omega_ref , Da_ref
    real (dp) , dimension (ndimmax)     :: pt

    character (len_default) , parameter :: format = ' ( 30 ( 1X , E16.8 ) ) '


    if ( inp % dis_rate_therm ) then
       sdr_ref = ( adi % lbda_ref / ( adi % rho_ref * adi % cp_ref ) ) / &
                 ( adi % L_ref * adi % L_ref )
    else
       sdr_ref = adi % D_ref / ( adi % L_ref * adi % L_ref )
    end if
    omega_ref = adi % u_ref * adi % u_ref / adi % time_ref
    Da_ref = ( 1.0_dp / adi % time_ref ) / ( adi % D_ref / ( adi % L_ref * adi % L_ref ) )


    call similarity_x_1Dbins ( inp , favavg % ca_T_Zm )
    call similarity_x_1Dbins ( inp , stat % favvar_ca_T_Zm )
    call similarity_x_1Dbins ( inp , favavg % ca_sdr_Zm )
    call similarity_x_1Dbins ( inp , stat % favvar_ca_sdr_Zm )
    call similarity_x_1Dbins ( inp , favavg % ca_Da_Zm )
    call similarity_x_1Dbins ( inp , stat % favvar_ca_Da_Zm )
    call similarity_x_1Dbins ( inp , favavg % ca_omega_Zm )
    call similarity_x_1Dbins ( inp , stat % favvar_ca_omega_Zm )
    call similarity_x_1Dbins ( inp , favavg % ca_FI_Zm )
    call similarity_x_1Dbins ( inp , stat % favvar_ca_FI_Zm )
    call similarity_x_1Dbins ( inp , favavg % ca_omega_FI )
    call similarity_x_1Dbins ( inp , stat % favvar_ca_omega_FI )
    call similarity_x_1Dbins ( inp , favavg % ca_omega_Da )
    call similarity_x_1Dbins ( inp , stat % favvar_ca_omega_Da )


    px = ex ; py = ey ; pz = ez ! initialize this to avoid bugs (verified)
    do iprobe = 1 , inp % nprobes

       pt (1) = inp % probe_coord ( iprobe , 1 )
       pt (2) = inp % probe_coord ( iprobe , 2 )
       pt (3) = inp % probe_coord ( iprobe , 3 )

       cond = .false.

!       if ( pt (1) >= x (sx-1) .and. pt (1) <= x (ex) .and. &
!            pt (2) >= y (sy-1) .and. pt (2) <= y (ey) .and. &
!            pt (3) >= z (sz-1) .and. pt (3) <= z (ez) ) cond = .true.

       if ( pt (1) >= x (sx) .and. pt (1) <= x (ex+1) .and. &
            pt (2) >= y (sy) .and. pt (2) <= y (ey+1) .and. &
            pt (3) >= z (sz) .and. pt (3) <= z (ez+1) ) cond = .true. ! A.Techer

       if (cond) then

          if ( iprobe == 1 ) write (*,*) 'writing conditional averages files ...'

          ! calculate points
          do i = ex+1 , sx , -1
             if ( x (i) > pt (1) ) px = i
          end do
          do j = ey+1 , sy , -1
             if ( y (j) > pt (2) ) py = j
          end do
          if ( ndim == 3 ) then
             do k = ez+1 , sz , -1
                if ( z (k) > pt (3) ) pz = k
             end do
          else
             pz = sz
          end if

          currunit = unit_cavg + iprobe - 1

          do n = 1,nbins

             write ( currunit , format ) min_Y + (n-1) * (max_Y-min_Y) / (nbins-1)                 , & ! 01
                                         favavg % ca_T_Zm (px,n) * adi % T_ref                     , & ! 02
                                         stat % favvar_ca_T_Zm (px,n) * adi % T_ref * adi % T_ref  , & ! 03

                                         min_Y + (n-1) * (max_Y-min_Y) / (nbins-1)                 , & ! 04
                                         favavg % ca_sdr_Zm (px,n) * sdr_ref                       , & ! 05
                                         stat % favvar_ca_sdr_Zm (px,n) * sdr_ref * sdr_ref        , & ! 06

                                         min_Y + (n-1) * (max_Y-min_Y) / (nbins-1)                 , & ! 07
                                         favavg % ca_Da_Zm (px,n) * Da_ref                         , & ! 08
                                         stat % favvar_ca_Da_Zm (px,n) * Da_ref * Da_ref           , & ! 09

                                         min_Y + (n-1) * (max_Y-min_Y) / (nbins-1)                 , & ! 10
                                         favavg % ca_omega_Zm (px,n) * omega_ref                   , & ! 11
                                         stat % favvar_ca_omega_Zm (px,n) * omega_ref * omega_ref  , & ! 12

                                         min_Y + (n-1) * (max_Y-min_Y) / (nbins-1)                 , & ! 13
                                         favavg % ca_FI_Zm (px,n)                                  , & ! 14
                                         stat % favvar_ca_FI_Zm (px,n)                             , & ! 15

                                         min_FI + (n-1) * (max_FI-min_FI) / (nbins-1)              , & ! 16
                                         favavg % ca_omega_FI (px,n) * omega_ref                   , & ! 17
                                         stat % favvar_ca_omega_FI (px,n) * omega_ref * omega_ref  , & ! 18

                ( log10 (min_Da) + (n-1) * (log10(max_Da)-log10(min_Da) ) / (nbins-1) )            , & ! 19
                                         favavg % ca_omega_Da (px,n) * omega_ref                   , & ! 20
                                         stat % favvar_ca_omega_Da (px,n) * omega_ref * omega_ref      ! 21

          end do

       end if

    end do


  end subroutine plot_cond_avgs


!> \brief Name detection utility.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine detect_name ( inp , name , condition )


    type (inp_type) , intent (in)                         :: inp       !< input derived type
    character (len_default) , intent (in)                 :: name      !< name string
    logical , intent (inout)                              :: condition !< match


    integer (ip) :: i


    condition = .false.


    do i = 1 , inp % nvars
       if ( inp % var_name (i) == trim (name) ) condition = .true.
    end do


  end subroutine detect_name


!> \brief Calculate Damkohler number.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine Damkohler ( inp , thd , dx_i , dy_i , dz_i , v , T , W_i , dm , tdr , rd , omegak , Da )


    type (inp_type) , intent (in)                                        :: inp    !< input derived type
    type (thd_type) , intent (in)                                        :: thd    !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)                :: dx_i   !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (in)                :: dy_i   !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (in)                :: dz_i   !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)          :: v      !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)            :: T      !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)            :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)          :: dm     !< diffusion matrix
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)          :: tdr    !< thermal diffusivity
    real (dp) , allocatable , dimension (:,:,:,:,:) , intent (in)        :: rd     !< density times diffusion
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)          :: omegak !< production rates array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)         :: Da     !< Damkohler number


    integer (ip)                                                         :: ok , i , j , k
    integer (ip)                                                         :: myspecies
    real (dp) , dimension (:,:,:) , allocatable                          :: Zm
    real (dp) , dimension (:,:,:) , allocatable                          :: wrk1 , wrk2 , wrk3 , &
                                                                            wrk4 , wrk5 , wrk6
    real (dp) , dimension (:,:,:,:) , allocatable                        :: rYVx , rYVy , rYVz


    if ( .not. reaction ) return


    myspecies = inp % index_Da


    allocate ( Zm   (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk1 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk2 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk3 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk4 (sx:ex,sy:ey,sz:ez)                   , &
               wrk5 (sx:ex,sy:ey,sz:ez)                   , &
               wrk6 (sx:ex,sy:ey,sz:ez)                   , &
               rYVx (sx:ex,sy:ey,sz:ez,nrv)               , &
               rYVy (sx:ex,sy:ey,sz:ez,nrv)               , &
               rYVz (sx:ex,sy:ey,sz:ez,nrv)               , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate Damkohler')


    if ( ndim == 2 ) then


       if ( vis .and. .not. eglib ) then
          call rhoYaVa_mix ( thd , dx_i , dy_i , dx_i , v , W_i , dm , rYVx , rYVy , rYVx )
       else if ( vis .and. eglib ) then
          call rhoYaVa_eglib ( thd , dx_i , dy_i , dx_i , T , W_i , tdr , rd , v , rYVx , rYVy , rYVx )
       end if

       wrk1 (sx:ex,sy:ey,sz:ez) = rYVx (sx:ex,sy:ey,sz:ez,myspecies)
       call comm_one (wrk1)
       call dx ( dx_i , wrk1 , wrk4 )
       wrk2 (sx:ex,sy:ey,sz:ez) = rYVy (sx:ex,sy:ey,sz:ez,myspecies)
       call comm_one (wrk2)
       call dy ( dy_i , wrk2 , wrk5 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk1 (i,j,k) = sqrt ( wrk4 (i,j,k) * wrk4 (i,j,k) + &
                                      wrk5 (i,j,k) * wrk5 (i,j,k) )
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       if ( vis .and. .not. eglib ) then
          call rhoYaVa_mix ( thd , dx_i , dy_i , dz_i , v , W_i , dm , rYVx , rYVy , rYVz )
       else if ( vis .and. eglib ) then
          call rhoYaVa_eglib ( thd , dx_i , dy_i , dz_i , T , W_i , tdr , rd , v , rYVx , rYVy , rYVz )
       end if

       wrk1 (sx:ex,sy:ey,sz:ez) = rYVx (sx:ex,sy:ey,sz:ez,myspecies)
       call comm_one (wrk1)
       call dx ( dx_i , wrk1 , wrk4 )
       wrk2 (sx:ex,sy:ey,sz:ez) = rYVy (sx:ex,sy:ey,sz:ez,myspecies)
       call comm_one (wrk2)
       call dy ( dy_i , wrk2 , wrk5 )
       wrk3 (sx:ex,sy:ey,sz:ez) = rYVz (sx:ex,sy:ey,sz:ez,myspecies)
       call comm_one (wrk3)
       call dz ( dz_i , wrk3 , wrk6 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk1 (i,j,k) = sqrt ( wrk4 (i,j,k) * wrk4 (i,j,k) + &
                                      wrk5 (i,j,k) * wrk5 (i,j,k) + &
                                      wrk6 (i,j,k) * wrk6 (i,j,k) )
             end do
          end do
       end do


    end if


    if ( npv == 0 ) then
       call Zmsimp ( inp , thd , v , Zm )
    else
       Zm (sx:ex,sy:ey,sz:ez) = v (sx:ex,sy:ey,sz:ez,niv+nrv+npv) / &
                                v (sx:ex,sy:ey,sz:ez,1)
    end if


    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             if ( Zm (i,j,k) > minZm .and. Zm (i,j,k) < maxZm ) then ! Da=0 in pure fuel and pure oxydant
                Da (i,j,k) = max ( epsi , omegak (i,j,k,myspecies) ) ! Da always positive
                Da (i,j,k) = v (i,j,k,1) * Da (i,j,k) / wrk1 (i,j,k)
             else
                Da (i,j,k) = epsi
             end if
          end do
       end do
    end do


    deallocate (Zm)
    deallocate ( wrk1 , wrk2 , wrk3 , wrk4 , wrk5 , wrk6 )
    deallocate ( rYVx , rYVy , rYVz )


  end subroutine Damkohler


!> \brief Calculate _numerical_ mixture fraction.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine Zmsimp ( inp , thd , v , Zm )


    type (inp_type) , intent (in)                                  :: inp !< input derived type
    type (thd_type) , intent (in)                                  :: thd !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v   !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: Zm  !< mixture fraction


    integer (ip)                :: i , j , k , l , m , ne
    real (dp)                   :: wrk , rho_i , num , denom , aj , aj0 , aj1
    real (dp) , dimension (nrv) :: a0 , a1


    ne = thd % nele


    ! specify the pure fuel
    a1 (1:nrv) = inp % Y1f (1:nrv)
!    call X_to_Y ( thd , a1 )

    ! specify the pure oxidizer
    a0 (1:nrv) = inp % Y0o (1:nrv)
!    call X_to_Y ( thd , a0 )


    do k = sz , ez
       do j = sy , ey
          do i = sx , ex

             rho_i = 1.0_dp / v (i,j,k,1)
             num   = 0.0_dp
             denom = 0.0_dp

             do l = 1 , ne

                aj0 = 0.0_dp ; aj1 = 0.0_dp ; aj = 0.0_dp
                do m = 1 , nrv
                   wrk = thd % KNCF (l,m) * thd % AWTc (l) * thd % Wc_i (m)
                   aj0 = aj0 + wrk * a0 (m)
                   aj1 = aj1 + wrk * a1 (m)
                   aj  = aj  + wrk * v (i,j,k,niv+m) * rho_i
                end do

                num   = num   + abs ( aj  - aj0 )
                denom = denom + abs ( aj1 - aj0 )

             end do

             Zm (i,j,k) = num / denom

          end do
       end do
    end do


  end subroutine Zmsimp


!> \brief Calculate _equivalent_ diffusion matrix.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine D_equivalent ( inp , thd , dx_i , dy_i , dz_i , v , T , W_i , ct , cp , dm , tdr , rd , Zm , Deq )


    type (inp_type) , intent (in)                                  :: inp  !< input derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , dimension (:) , allocatable , intent (in)          :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dy_i !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dz_i !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v    !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: ct   !< thermal conductivity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: dm   !< diffusion matrix
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: tdr  !< thermal diffusivity
    real (dp) , allocatable , dimension (:,:,:,:,:) , intent (in)  :: rd   !< density times diffusion
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: Zm   !< numerical mixture fraction
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: Deq  !< equivalent diffusion matrix


    integer (ip)                                        :: ok , i , j , k , l , m , ne
    real (dp)                                           :: wrk , rho_i , num , denom , aj0 , aj1 , ajx , ajy , ajz
    real (dp) , dimension (nrv)                         :: a0 , a1
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3
    real (dp) , dimension (:,:,:,:) , allocatable       :: rYVx , rYVy , rYVz


    if ( inp % dis_rate_therm ) then ! hypothesis: same Schmidt for all species and Lewis equal to unity


       Deq (sx:ex,sy:ey,sz:ez) = ct (sx:ex,sy:ey,sz:ez) / &
                               ( v (sx:ex,sy:ey,sz:ez,1) * cp (sx:ex,sy:ey,sz:ez) )


    else


       allocate ( wrk1 (sx:ex,sy:ey,sz:ez)     , &
                  wrk2 (sx:ex,sy:ey,sz:ez)     , &
                  wrk3 (sx:ex,sy:ey,sz:ez)     , &
                  rYVx (sx:ex,sy:ey,sz:ez,nrv) , &
                  rYVy (sx:ex,sy:ey,sz:ez,nrv) , &
                  rYVz (sx:ex,sy:ey,sz:ez,nrv) , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate D_equivalent')


       ne = thd % nele


       ! specify the pure fuel
       a1 (1:nrv) = inp % Y1f (1:nrv)
       !    call X_to_Y ( thd , Yf )

       ! specify the pure oxidizer
       a0 (1:nrv) = inp % Y0o (1:nrv)
       !    call X_to_Y ( thd , Yo )


       if ( ndim == 2 ) then ! 2D problem


          ! Zm gradient
          call comm_one (Zm)
          call dx ( dx_i , Zm , wrk1 )
          call dy ( dy_i , Zm , wrk2 )

          ! Zm norm
          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   wrk1 (i,j,k) = sqrt ( wrk1 (i,j,k) * wrk1 (i,j,k) + &
                                         wrk2 (i,j,k) * wrk2 (i,j,k) )
                end do
             end do
          end do

          rYVz (sx:ex,sy:ey,sz:ez,:) = 0.0_dp
          if ( vis .and. .not. eglib ) then
             call rhoYaVa_mix ( thd , dx_i , dy_i , dy_i , v , W_i , dm , rYVx , rYVy , rYVx )
          else if ( vis .and. eglib ) then
             call rhoYaVa_eglib ( thd , dx_i , dy_i , dx_i , T , W_i , tdr , rd , v , rYVx , rYVy , rYVx )
          end if


       else if ( ndim == 3 ) then ! 3D problem


          ! Zm gradient
          call comm_one (Zm)
          call dx ( dx_i , Zm , wrk1 )
          call dy ( dy_i , Zm , wrk2 )
          call dz ( dz_i , Zm , wrk3 )

          ! Zm gradient norm
          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   wrk1 (i,j,k) = sqrt ( wrk1 (i,j,k) * wrk1 (i,j,k) + &
                                         wrk2 (i,j,k) * wrk2 (i,j,k) + &
                                         wrk3 (i,j,k) * wrk3 (i,j,k) )
                end do
             end do
          end do

          if ( vis .and. .not. eglib ) then
             call rhoYaVa_mix ( thd , dx_i , dy_i , dz_i , v , W_i , dm , rYVx , rYVy , rYVz )
          else if ( vis .and. eglib ) then
             call rhoYaVa_eglib ( thd , dx_i , dy_i , dz_i , T , W_i , tdr , rd , v , rYVx , rYVy , rYVz )
          end if


       end if


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                rho_i = 1.0_dp / v (i,j,k,1)
                num   = 0.0_dp
                denom = 0.0_dp

                do l = 1 , ne

                   aj0 = 0.0_dp ; aj1 = 0.0_dp
                   ajx = 0.0_dp ; ajy = 0.0_dp ; ajz = 0.0_dp

                   do m = 1 , nrv

                      wrk = thd % KNCF (l,m) * thd % AWTc (l) * thd % Wc_i (m)

                      aj0 = aj0 + wrk * a0 (m)
                      aj1 = aj1 + wrk * a1 (m)

                      ajx = ajx + wrk * rYVx (i,j,k,m)
                      ajy = ajy + wrk * rYVy (i,j,k,m)
                      ajz = ajz + wrk * rYVz (i,j,k,m)

                   end do

                   num   = num   + sqrt ( ajx*ajx + ajy*ajy + ajz*ajz )
                   denom = denom + abs ( aj1 - aj0 )

                end do

                denom = denom * wrk1 (i,j,k) * v (i,j,k,1)

                if ( Zm (i,j,k) >= 0.001_dp .and. Zm (i,j,k) <= 0.999 ) then
                   Deq (i,j,k) = num / denom
                else
                   Deq (i,j,k) = 0.0_dp
                end if

             end do
          end do
       end do

       deallocate ( wrk1 , wrk2 , wrk3 )
       deallocate ( rYVx , rYVy , rYVz )


    end if


  end subroutine D_equivalent


!> \brief Calculate scalar disspation rate.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine dissip_rate ( dx_i , dy_i , dz_i , Deq , Zm , dr )


    real (dp) , dimension (:) , allocatable , intent (in)          :: dx_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dy_i !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)          :: dz_i !< inverted dx array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: Deq  !< equivalent diffusion matrix
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: Zm   !< numerical mixture fraction
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: dr   !< scalar dissipation rate


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3


    allocate ( wrk1 (sx:ex,sy:ey,sz:ez) , &
               wrk2 (sx:ex,sy:ey,sz:ez) , &
               wrk3 (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate dissip_rate')


    if ( ndim == 2 ) then ! 2D problem

       ! Zm gradient
       call comm_one (Zm)
       call dx ( dx_i , Zm , wrk1 )
       call dy ( dy_i , Zm , wrk2 )

       ! Zm gradient norm
       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dr (i,j,k) = wrk1 (i,j,k) * wrk1 (i,j,k) + &
                             wrk2 (i,j,k) * wrk2 (i,j,k)
                dr (i,j,k) = dr (i,j,k) * Deq (i,j,k)
                dr (i,j,k) = max ( dr (i,j,k) , epsi )
             end do
          end do
       end do

    else if ( ndim == 3 ) then ! 3D problem

       ! Zm gradient
       call comm_one (Zm)
       call dx ( dx_i , Zm , wrk1 )
       call dy ( dy_i , Zm , wrk2 )
       call dz ( dz_i , Zm , wrk3 )

       ! Zm gradient norm
       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                dr (i,j,k) = wrk1 (i,j,k) * wrk1 (i,j,k) + &
                             wrk2 (i,j,k) * wrk2 (i,j,k) + &
                             wrk3 (i,j,k) * wrk3 (i,j,k)
                dr (i,j,k) = dr (i,j,k) * Deq (i,j,k)
                dr (i,j,k) = max ( dr (i,j,k) , epsi )
             end do
          end do
       end do

    end if


    deallocate ( wrk1 , wrk2 , wrk3 )


  end subroutine dissip_rate


!> \brief Calculate discrete Fourier transform.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine DFT ( start , end , dft_var )


    integer (ip) , intent (in)                                     :: start   !< start interval
    integer (ip) , intent (in)                                     :: end     !< end interval
    real (dp) , allocatable , dimension (:) , intent (inout)       :: dft_var !< discrete Fourier transform


    real (dp) , parameter                           :: two_pi = ( acos(-1.0_dp) ) + ( acos(-1.0_dp) )
    integer (ip)                                    :: ok , l , m , k , n , Ntot
    real (dp)                                       :: wrk_r , wrk_i , wrk1 , wrk2
    real (dp) , dimension (:,:) , allocatable       :: wrk


    allocate ( wrk ( ntimes , 2 ) , stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate DFT')


    Ntot = end - start + 1


    k = 0
    do m = start , end

       wrk1 = two_pi * real ( k , dp ) / Ntot

       wrk_r = 0.0_dp
       wrk_i = 0.0_dp

       n = 0
       do l = start , end

          wrk2 = wrk1 * real ( n , dp )

          wrk_r = wrk_r + dft_var (l) * cos (wrk2)
          wrk_i = wrk_i + dft_var (l) * sin (wrk2)

          n = n+1

       end do

       wrk (m,1) =   wrk_r / Ntot
       wrk (m,2) = - wrk_i / Ntot

       k = k+1

    end do


    ! the spectra is the square of the amplitude
    do m = start , end
       dft_var (m) = wrk (m,1) * wrk (m,1) + &
                     wrk (m,2) * wrk (m,2)
    end do


    deallocate (wrk)


  end subroutine DFT


!> \brief Calculate _simplified_ species diffusion flux.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rhoYaVa_mix ( thd , dx_i , dy_i , dz_i , v , W_i , dm , rhoYVx , rhoYVy , rhoYVz )


    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i   !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i   !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i   !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v      !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: dm     !< diffusion matrix
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: rhoYVx !< x-component of species diffusion flux
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: rhoYVy !< y-component of species diffusion flux
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: rhoYVz !< z-component of species diffusion flux


    integer (ip)                                  :: nfluid
    integer (ip)                                  :: ok , domain_id , i , j , k , l
    real (dp)                                     :: wrk , Ya , rho_i
    real (dp)                                     :: vdcx , vdtx , &
                                                     vdcy , vdty , &
                                                     vdcz , vdtz
    real (dp) , dimension (:,:,:,:) , allocatable :: Xa
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable :: fd


    nfluid = ndim * ( ndim + 1 )


    if ( ndim == 2 ) then ! 2D problem


       ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
       allocate ( Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate rhoYaVa_mix_2D')


       do domain_id = -ndim , ndim
          call molefrac ( domain_id , thd , W_i , v , Xa )
       end do


       fd (:,:,:,:) = 0.0_dp
       do l = 1 , nrv+npv
          call dx_fixed1 ( dx_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) )
          call dy_fixed1 ( dy_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) )
       end do


       ! communicate the first derivative
       call comm_deriv (fd)


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                rho_i = 1.0_dp / v (i,j,k,1)

                ! Correction velocity
                vdcx = 0.0_dp ; vdcy = 0.0_dp
                do l = 1 , nrv
                   wrk  = dm (i,j,k,l) * thd % Wc (l) * W_i (i,j,k)
                   vdcx = vdcx + wrk * fd (i,j,k,nfluid+ndim*l-1)
                   vdcy = vdcy + wrk * fd (i,j,k,nfluid+ndim*l  )
                end do

                do l = 1 , nrv

                   wrk = - dm (i,j,k,l) * thd % Wc (l) * W_i (i,j,k)

                   Ya = v (i,j,k,niv+l) * rho_i

                   vdtx = vdcx * Ya ; vdty = vdcy * Ya

                   ! plus diffussion velocity
                   vdtx = vdtx + wrk * fd (i,j,k,nfluid+ndim*l-1)
                   vdty = vdty + wrk * fd (i,j,k,nfluid+ndim*l  )

                   ! finally ... chi j = rho * Yalpha * V alpha j
                   rhoYVx (i,j,k,l) = v (i,j,k,1) * vdtx
                   rhoYVy (i,j,k,l) = v (i,j,k,1) * vdty

                end do

             end do
          end do
       end do


       deallocate ( Xa , fd )


    else if ( ndim == 3 ) then ! 3D problem


       ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
       allocate ( Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate rhoYaVa_mix_3D')


       do domain_id = -ndim , ndim
          call molefrac ( domain_id , thd , W_i , v , Xa )
       end do


       fd (:,:,:,:) = 0.0_dp
       do l = 1 , nrv+npv
          call dx_fixed1 ( dx_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-2) )
          call dy_fixed1 ( dy_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) )
          call dz_fixed1 ( dz_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) )
       end do


       ! communicate the first derivative
       call comm_deriv (fd)


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                rho_i = 1.0_dp / v (i,j,k,1)

                ! Correction velocity
                vdcx = 0.0_dp ; vdcy = 0.0_dp ; vdcz = 0.0_dp
                do l = 1 , nrv
                   wrk  = dm (i,j,k,l) * thd % Wc (l) * W_i (i,j,k)
                   vdcx = vdcx + wrk * fd (i,j,k,nfluid+ndim*l-2)
                   vdcy = vdcy + wrk * fd (i,j,k,nfluid+ndim*l-1)
                   vdcz = vdcz + wrk * fd (i,j,k,nfluid+ndim*l  )
                end do

                do l = 1 , nrv

                   wrk = - dm (i,j,k,l) * thd % Wc (l) * W_i (i,j,k)

                   Ya = v (i,j,k,niv+l) * rho_i

                   vdtx = vdcx * Ya ; vdty = vdcy * Ya ; vdtz = vdcz * Ya

                   ! plus diffussion velocity
                   vdtx = vdtx + wrk * fd (i,j,k,nfluid+ndim*l-2)
                   vdty = vdty + wrk * fd (i,j,k,nfluid+ndim*l-1)
                   vdtz = vdtz + wrk * fd (i,j,k,nfluid+ndim*l  )

                   ! finally ... chi j = rho * Yalpha * V alpha j
                   rhoYVx (i,j,k,l) = v (i,j,k,1) * vdtx
                   rhoYVy (i,j,k,l) = v (i,j,k,1) * vdty
                   rhoYVz (i,j,k,l) = v (i,j,k,1) * vdtz

                end do

             end do
          end do
       end do


       deallocate ( Xa , fd )


    end if


  end subroutine rhoYaVa_mix


!> \brief Calculate _detailed_ species diffusion flux.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine rhoYaVa_eglib ( thd , dx_i , dy_i , dz_i , T , W_i , tdr , rd , v , rhoYVx , rhoYVy , rhoYVz )


    type (thd_type) , intent (in)                                  :: thd    !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i   !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i   !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i   !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T      !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i    !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: tdr    !< thermal diffusivity
    real (dp) , allocatable , dimension (:,:,:,:,:) , intent (in)  :: rd     !< density times diffusion
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v      !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: rhoYVx !< x-component of species diffusion flux
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: rhoYVy !< y-component of species diffusion flux
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: rhoYVz !< z-component of species diffusion flux


    integer (ip)                                    :: nfluid
    integer (ip)                                    :: ok , domain_id , i , j , k , l , m
    integer (ip)                                    :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                       :: rho_i , Ts , Ps
    real (dp) , dimension (nrv)                     :: Xa_ , d_x , d_y , d_z
    real (dp)                                       :: wrk1 , wrk2 , wrk3 , wrk4 , wrk5 , wrk6
    ! additional variables
    real (dp) , dimension (:,:,:) , allocatable     :: P
    real (dp) , dimension (:,:,:,:) , allocatable   :: Xa
    ! vector containing the first derivative to be communicated
    real (dp) , dimension (:,:,:,:) , allocatable   :: fd


    nfluid = ndim * ( ndim + 2 )


    if ( ndim == 2 ) then ! 2D problem


       ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
       allocate ( P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
                  Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate rhoYaVa_eglib_2D')


       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          call molefrac ( domain_id , thd , W_i , v , Xa )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   P (i,j,k) = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
                end do
             end do
          end do

       end do


       fd (:,:,:,:) = 0.0_dp
       call dx_fixed1 ( dx_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,5) )
       call dy_fixed1 ( dy_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,6) )

       call dx_fixed1 ( dx_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,7) )
       call dy_fixed1 ( dy_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)     , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,8) ) ! nfluid = 8 in 2D

       do l = 1 , nrv+npv
          call dx_fixed1 ( dx_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) )
          call dy_fixed1 ( dy_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) )
       end do

       deallocate ( P , Xa )


       ! communicate the first derivative
       call comm_deriv (fd)


       ! rhoYVx , rhoYVy
       rhoYVx (sx:ex,sy:ey,sz:ez,:) = 0.0_dp
       rhoYVy (sx:ex,sy:ey,sz:ez,:) = 0.0_dp


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                rho_i = 1.0_dp / v (i,j,k,1)
                Ps    = 1.0_dp / ( v (i,j,k,1) * T (i,j,k) * W_i (i,j,k) )
                Ts    = 1.0_dp / T (i,j,k)

                wrk1 = fd (i,j,k,5) * Ts
                wrk2 = fd (i,j,k,6) * Ts

                wrk3 = fd (i,j,k,7) * Ps
                wrk4 = fd (i,j,k,8) * Ps

                do l = 1 , nrv
                   ! Xa_
                   Xa_ (l) = v (i,j,k, l+niv ) /                        &
                           ( v (i,j,k,1) * W_i (i,j,k) * thd % Wc (l) )
                   ! d_x
                   d_x (l) = fd (i,j,k,nfluid+ndim*l-1) +          &
                           ( Xa_ (l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk3
                   ! d_y
                   d_y (l) = fd (i,j,k,nfluid+ndim*l  ) +          &
                           ( Xa_ (l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk4
                end do

                do l = 1 , nrv

                   fd (i,j,k,nfluid+ndim*l-1) = 0.0_dp ! rYV_x
                   fd (i,j,k,nfluid+ndim*l  ) = 0.0_dp ! rYV_y

                   do m = 1 , nrv
                      fd (i,j,k,nfluid+ndim*l-1) = fd (i,j,k,nfluid+ndim*l-1) +               &
                                                   rd (i,j,k,l,m) *                           &
                                                 ( d_x (m) + Xa_ (m) * tdr (i,j,k,m) * wrk1 )
                      fd (i,j,k,nfluid+ndim*l  ) = fd (i,j,k,nfluid+ndim*l  ) +               &
                                                   rd (i,j,k,l,m) *                           &
                                                 ( d_y (m) + Xa_ (m) * tdr (i,j,k,m) * wrk2 )
                   end do

                   rhoYVx (i,j,k,l) = - fd (i,j,k,nfluid+ndim*l-1)
                   rhoYVy (i,j,k,l) = - fd (i,j,k,nfluid+ndim*l  )

                end do

             end do
          end do
       end do


       deallocate (fd)


    else if ( ndim == 3 ) then ! 3D problem


       ! ALLOCATION MUST REMAIN AS IN 3D BECAUSE OF THE DX_FIXED AND DX_FIXED_FIXED
       allocate ( P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)         , &
                  Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nrv+npv) , &
                  fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nderiv)  , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate rhoYaVa_eglib_3D')


       do domain_id = -ndim , ndim

          call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )

          call molefrac ( domain_id , thd , W_i , v , Xa )

          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   P (i,j,k) = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
                end do
             end do
          end do

       end do


       fd (:,:,:,:) = 0.0_dp
       call dx_fixed1 ( dx_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,10) )
       call dy_fixed1 ( dy_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,11) )
       call dz_fixed1 ( dz_i , T  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,12) )

       call dx_fixed1 ( dx_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,13) )
       call dy_fixed1 ( dy_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,14) )
       call dz_fixed1 ( dz_i , P  (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng)      , &
                               fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,15) ) ! nfluid = 15 in 3D

       do l = 1 , nrv+npv
          call dx_fixed1 ( dx_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-2) )
          call dy_fixed1 ( dy_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l-1) )
          call dz_fixed1 ( dz_i                                                       , &
                           Xa (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,l)                 , &
                           fd (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng,nfluid+ndim*l  ) )
       end do

       deallocate ( P , Xa )


       ! communicate the first derivative
       call comm_deriv (fd)


       ! rhoYVx , rhoYVy , rhoYVz
       rhoYVx (sx:ex,sy:ey,sz:ez,:) = 0.0_dp
       rhoYVy (sx:ex,sy:ey,sz:ez,:) = 0.0_dp
       rhoYVz (sx:ex,sy:ey,sz:ez,:) = 0.0_dp


       ! rYV_x , rYV_y , rYV_z
       do k = sz,ez
          do j = sy,ey
             do i = sx,ex

                rho_i = 1.0_dp / v (i,j,k,1)
                Ps    = 1.0_dp / ( v (i,j,k,1) * T (i,j,k) * W_i (i,j,k) )
                Ts    = 1.0_dp / T (i,j,k)

                wrk1 = fd (i,j,k,10) * Ts
                wrk2 = fd (i,j,k,11) * Ts
                wrk3 = fd (i,j,k,12) * Ts

                wrk4 = fd (i,j,k,13) * Ps
                wrk5 = fd (i,j,k,14) * Ps
                wrk6 = fd (i,j,k,15) * Ps

                do l = 1 , nrv
                   ! Xa_
                   Xa_ (l) = v (i,j,k, l+niv ) /                        &
                           ( v (i,j,k,1) * W_i (i,j,k) * thd % Wc (l) )
                   ! d_x
                   d_x (l) = fd (i,j,k,nfluid+ndim*l-2) +          &
                           ( Xa_ (l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk4
                   ! d_y
                   d_y (l) = fd (i,j,k,nfluid+ndim*l-1) +          &
                           ( Xa_ (l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk5
                   ! d_z
                   d_z (l) = fd (i,j,k,nfluid+ndim*l  ) +          &
                           ( Xa_ (l) - v (i,j,k,niv+l) * rho_i ) * &
                             wrk6
                end do

                do l = 1 , nrv

                   fd (i,j,k,nfluid+ndim*l-2) = 0.0_dp ! rYV_x
                   fd (i,j,k,nfluid+ndim*l-1) = 0.0_dp ! rYV_y
                   fd (i,j,k,nfluid+ndim*l  ) = 0.0_dp ! rYV_z

                   do m = 1 , nrv

                      fd (i,j,k,nfluid+ndim*l-2) = fd (i,j,k,nfluid+ndim*l-2) +               &
                                                   rd (i,j,k,l,m) *                           &
                                                 ( d_x (m) + Xa_ (m) * tdr (i,j,k,m) * wrk1 )
                      fd (i,j,k,nfluid+ndim*l-1) = fd (i,j,k,nfluid+ndim*l-1) +               &
                                                   rd (i,j,k,l,m) *                           &
                                                 ( d_y (m) + Xa_ (m) * tdr (i,j,k,m) * wrk2 )
                      fd (i,j,k,nfluid+ndim*l  ) = fd (i,j,k,nfluid+ndim*l  ) +               &
                                                   rd (i,j,k,l,m) *                           &
                                                 ( d_z (m) + Xa_ (m) * tdr (i,j,k,m) * wrk3 )

                   end do

                   rhoYVx (i,j,k,l) = - fd (i,j,k,nfluid+ndim*l-2)
                   rhoYVy (i,j,k,l) = - fd (i,j,k,nfluid+ndim*l-1)
                   rhoYVz (i,j,k,l) = - fd (i,j,k,nfluid+ndim*l  )

                end do

             end do
          end do
       end do


       deallocate (fd)


    end if


  end subroutine rhoYaVa_eglib


!> \brief Simulate extrapolation BCs for post-treatment.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

!  subroutine bc_extrapolation ( face , v )                                                           ! Not need anymore because of use solver_reaction which uses BCs


!    integer (ip) , intent (in)                                     :: face !< face domain
!    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


!    integer (ip) :: i , j , k , l


!    if ( face == W ) then


!       do l = 1 , nv
!          do k = sz , ez
!             do j = sy , ey
!                do i = 1 , ng
!                   v (sx-i,j,k,l) = v (sx,j,k,l)
!                end do
!             end do
!          end do
!       end do


!    else if ( face == E ) then


!       do l = 1 , nv
!          do k = sz , ez
!             do j = sy , ey
!                do i = 1 , ng
!                   v (ex+i,j,k,l) = v (ex,j,k,l)
!                end do
!             end do
!          end do
!       end do


!    else if ( face == S ) then


!       do l = 1 , nv
!          do k = sz , ez
!             do j = 1 , ng
!                do i = sx , ex
!                   v (i,sy-j,k,l) = v (i,sy,k,l)
!                end do
!             end do
!          end do
!       end do


!    else if ( face == N ) then


!       do l = 1 , nv
!          do k = sz , ez
!             do j = 1 , ng
!                do i = sx , ex
!                   v (i,ey+j,k,l) = v (i,ey,k,l)
!                end do
!             end do
!          end do
!       end do


!    else if ( face == B ) then


!       do l = 1 , nv
!          do k = 1 , ng
!             do j = sy , ey
!                do i = sx , ex
!                   v (i,j,sz-k,l) = v (i,j,sz,l)
!                end do
!             end do
!          end do
!       end do


!    else if ( face == F ) then


!       do l = 1 , nv
!          do k = 1 , ng
!             do j = sy , ey
!                do i = sx , ex
!                   v (i,j,ez+k,l) = v (i,j,ez,l)
!                end do
!             end do
!          end do
!       end do


!    end if


!  end subroutine bc_extrapolation


!> \brief Calculate the time step.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine timestep ( inp , adi , thd , dx_i , dy_i , dz_i , cp , W_i , T , v , DTmax , dtmin , dt )


    type (inp_type) , intent (in)                                  :: inp   !< input derived type
    type (adi_type)                                                :: adi   !< non-dimensional derived type
    type (thd_type)                                                :: thd   !< thermodynamic derived type
    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i  !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i  !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i  !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp    !< heat capacity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i   !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T     !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v     !< conserved variables array
    real (dp) , intent (in)                                        :: DTmax !< maximum temperature increase
    real (dp) , intent (in)                                        :: dtmin !< minimum previous time step
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: dt    !< array of different time steps


    integer (ip)                                          :: ok , i , j , k
    real (dp)                                             :: rho_i                                  , &
                                                             dtcx_i  , dtcy_i   , dtcz_i            , &
                                                             dtc , gamma  , P                       , &
                                                             csound
    real (dp) , allocatable , dimension (:,:,:)           :: ux , vy , wz , hmin
    real (dp) , allocatable , dimension (:,:,:,:)         :: Xa , dm
    real (dp) , allocatable , dimension (:,:,:)           :: ct

    ! SGS variables
    real (dp) , allocatable , dimension (:,:,:)           :: mut
    real (dp)                                             :: dm_sgs , ct_sgs

    ! vector containing the first derivative
    real (dp) , allocatable , dimension (:,:,:,:)         :: fd


    allocate  ( ux   ( sx:ex , sy:ey , sz:ez ) , &
                vy   ( sx:ex , sy:ey , sz:ez ) , &
                wz   ( sx:ex , sy:ey , sz:ez ) , &
                hmin ( sx:ex , sy:ey , sz:ez )                   , &
!                fd   ( sx:ex , sy:ey , sz:ez , ndim*ndim)           , &
                stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate timestep')


    hmin = 1.0_dp
    dtc  = 1.0_dp
    dt   = 1.0_dp


    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             rho_i      = 1.0_dp / v (i,j,k,1)
             ux (i,j,k) = v (i,j,k,2) * rho_i
             vy (i,j,k) = v (i,j,k,3) * rho_i
             wz (i,j,k) = v (i,j,k,4) * rho_i
          end do
       end do
    end do


    ! Convective criteria
    do k = sz , ez
       do j = sy , ey
          do i = sx , ex

             hmin (i,j,k) = max ( dx_i (i) , dy_i (j) , dz_i (k) )
             hmin (i,j,k) = 1.0_dp / hmin (i,j,k)

             gamma  = thd % gam2 * W_i (i,j,k) / cp (i,j,k)
             gamma  = 1.0_dp / ( 1.0_dp - gamma )
             P      = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             csound = sqrt ( gamma * P * rho_i )

             dtcx_i = ( abs ( ux (i,j,k) ) + csound ) * dx_i (i)
             dtcy_i = ( abs ( vy (i,j,k) ) + csound ) * dy_i (j)
             dtcz_i = ( abs ( wz (i,j,k) ) + csound ) * dz_i (k)

             dtc = max ( dtcx_i , dtcy_i , dtcz_i )

             dt (i,j,k,1) = CFL / dtc

          end do
       end do
    end do


    ! To estimate the time step, the viscous criteria will be
    ! the one from simple transport even if EGLIB is being used
    if ( vis .and. vis_dt ) then

       allocate ( ct  ( sx:ex , sy:ey , sz:ez )         , &
                  Xa  ( sx:ex , sy:ey , sz:ez , 1:nrv ) , &
                  dm  ( sx:ex , sy:ey , sz:ez , 1:nrv ) , &
!                  mut ( sx:ex , sy:ey , sz:ez )         , &
                  stat = ok )
       if ( ok > 0 ) call abort_mpi ('error allocate timestep 1')


       call molefrac            ( 0 , thd , W_i , v , Xa )
       call prim_vis_vars_wo_mu ( 0 , thd , v , W_i , T , Xa , dm , ct )
       deallocate (Xa)


       if ( inp % LES ) then

          if ( rank == rank_default ) then 
             write (*,*) 'WARNING: call timestep in postreatment with LES calculation'
          end if

!          calculate fd... by using dx, dx_fixed1 or dx_fixed2 ...
!
!          call SGS_selector ( inp , adi , 0 , dx_i , dy_i , dz_i , v , fd , shk , mut )
!
!          dm_sgs = maxval ( adi % Sc * mut (sx:ex,sy:ey,sz:ez) &
!                            / ( v (sx:ex,sy:ey,sz:ez,1) * inp % Sc_sgs ) )
!          dm (sx:ex,sy:ey,sz:ez,1:nrv) = dm (sx:ex,sy:ey,sz:ez,1:nrv) + dm_sgs
!
!          ct_sgs = maxval ( adi % Pr * mut (sx:ex,sy:ey,sz:ez) &
!                            * cp (sx:ex,sy:ey,sz:ez) / inp % Pr_sgs )
!          ct (sx:ex,sy:ey,sz:ez) = ct (sx:ex,sy:ey,sz:ez) + ct_sgs

       end if


       do k = sz , ez
          do j = sy , ey
             do i = sx , ex

                ! Mass diffusivity
                dt (i,j,k,2) = Fo * hmin (i,j,k) * hmin (i,j,k) / maxval( dm (i,j,k,1:nrv) )
                dt (i,j,k,2) = dt (i,j,k,2) * adi % L_ref * adi % L_ref / ( adi % time_ref * adi % D_ref )

                ! Thermal diffusivity
                dt (i,j,k,3) = Fo * hmin (i,j,k) * hmin (i,j,k) * v (i,j,k,1) * cp (i,j,k) / ct (i,j,k)
                dt (i,j,k,3) = dt (i,j,k,3) * adi % L_ref * adi % L_ref * adi % rho_ref * adi % cp_ref / &
                               ( adi % lbda_ref * adi % time_ref )

             end do
          end do
       end do

       deallocate ( dm , ct )

    end if


    ! chemical criteria
    if (reaction) then

       call abort_mpi ('error: calculate dt_chem, DTmax and dtmin are unkonwn...')

       ! if 2 convective and 1 reactive: dtmin and dTmax is somehow overstimated
       ! if 1 convective and 2 reactive: 0.5 * dtmin and dTmax is exact

       dt (i,j,k,4) = dtmin * CPR / DTmax       ! maximum CPR degrees Kelvin of temperature
       dt (i,j,k,4) = min ( dt (i,j,k,4) , 1.0_dp )

       ! The minimum chemical time step is CFLmin
       dt (i,j,k,4) = max ( dt (i,j,k,4) , dt (i,j,k,1) * CFLmin / CFL )

    end if


    do k = sz , ez
       do j = sy , ey
          do i = sx , ex

             dt (i,j,k,5) = minval ( dt (i,j,k,1:4) )

          end do
       end do
    end do

    write (*,'(1X,A,1X,1PE10.2,4(1X,I10))') 'Localization dtmin, i, j, k, dt_index =' ,&
                                       minval ( dt (sx:ex,sy:ey,sz:ez,1:4) ) * adi % time_ref , &
                                       minloc ( dt (sx:ex,sy:ey,sz:ez,1:4) )

    deallocate ( ux , vy , wz , hmin )!, fd , mut )


  end subroutine timestep


end module tools_post
