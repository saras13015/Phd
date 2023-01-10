module IO

    use parameters
    use input
    use adim  
    use parallel
    use solver
    use type_thd
    use tools
    use BCs    
    use iohdf5

    implicit none

    contains

    !     !> \brief Create the necessary repertories for izem.
    !     subroutine createrep 

    !         logical                       :: existed

    ! #ifdef __INTEL_COMPILER
    ! #define _DIR_ directory
    ! #else
    ! #define _DIR_ file    
    ! #endif

    !         inquire ( _DIR_ = trim (dir_restart) , exist = existed )
    !         if (.not.existed) call system ( 'mkdir ' // trim (dir_restart) )

    !         inquire ( _DIR_ = trim (dir_plot) , exist = existed )
    !         if (.not.existed) call system ( 'mkdir ' // trim (dir_plot) )

    !     end subroutine createrep



    subroutine get_variable_names ( var_names )

        character (len_default) , intent(inout) :: var_names (nv)

        var_names(1)   = "rho"
        var_names(2:4) = (/ "rhoU" , "rhoV" , "rhoW"/)
        var_names(5)   = "rhoE"

    end subroutine get_variable_names


    !> \brief Save restart file.
    !!
    subroutine save_restart_hdf5 ( adi , thd , number , time , v )

        type (adi_type)                               , intent (in)   :: adi  
        type (thd_type)                               , intent (in)   :: thd
        integer (ip)                                  , intent (in)   :: number !< number of restart file
        real (dp)                                     , intent (in)   :: time   !< time
        real (dp) , allocatable , dimension (:,:,:,:) , intent (in)   :: v      !< conserved variables array


        character (len_default) :: restart_ascii  , adr_file
        character (len_default) :: var_names (nv)
        integer (ip)            :: i

        write ( restart_ascii , format_restart ) number
        adr_file = trim (dir_restart) //   trim (file_restart) // '_' // trim (restart_ascii)

        adr_file = trim (adr_file) // ".h5"

        if ( rank == rank_default ) write (*,*) '....writing file ' , trim ( adr_file )
        if ( rank == rank_default )  write (*,*) '==================================================='

        call get_variable_names ( var_names )

        call hdf_create_file ( adr_file )

        call hdf_open_or_create_group ( "/consvar" )

        do i = 1 , niv

            call hdf_write_r3d_array ( v(sx:ex,sy:ey,sz:ez,i) , var_names(i) , 3 , &
                &                 (/ ntx , nty , ntz /) , subdom_size , subdom_offset)

        enddo


        do i = 1 , nrv

            call hdf_write_r3d_array ( v(sx:ex,sy:ey,sz:ez,niv+i) , thd % NameSpecy(i) , 3 , &
                &                 (/ ntx , nty , ntz /) , subdom_size , subdom_offset)

        enddo


        call hdf_close_group()

        call hdf_open_or_create_group("/simulation")

        call hdf_write_r1d_array ( (/time /) , "time" , 1 , (/1/) , (/1/) )

        call hdf_close_group() 

        call hdf_close_file()

    end subroutine save_restart_hdf5


    subroutine read_restart_hdf5 ( inp , adi , thd , time , x , y , z , &
        dx_i , dy_i , dz_i , T , W_i , cp , ha , v )


    type (inp_type) , intent (inout)                               :: inp  !< input derived type
    type (adi_type) , intent (in)                                  :: adi  !< non-dimensional derived type
    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , intent (in)                                        :: time !< time
    real (dp) , allocatable , dimension (:) , intent (in)          :: x    !< x-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: y    !< y-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: z    !< z-coordinate array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dx_i !< inverted dx array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dy_i !< inverted dy array
    real (dp) , allocatable , dimension (:) , intent (in)          :: dz_i !< inverted dz array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T    !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i  !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp   !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha   !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v    !< conserved variables array


    integer (ip) , dimension (ndimmax+ndimmax) , parameter         :: face_domain = (/ -1 , 1 , -2 , 2 , -3 , 3 /)
    integer (ip)                                                   :: ok , i , j , k , l , ite
    integer (kind=8)                                               :: reclmax
    character (len_default)                                        :: number_ascii , rank_ascii , adr_file
    real    (dp)                                                   :: dtmin
    real    (dp) , dimension (4)                                   :: dt_4
    real (dp) , allocatable , dimension (:,:,:)                    :: mu_SGS !< turbulent viscosity (for LES)
    character (len_default) :: var_names (nv)
    character (len_default) :: restart_ascii


    dt_4 = 1.0_dp

    allocate ( v       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nv  ) , &
        ha      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng , 1:nrv ) , &
        T       ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
        W_i     ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
        cp      ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
        mu_SGS  ( sx-ng:ex+ng , sy-ng:ey+ng , sz-ng:ez+ng )         , &
        stat = ok )


    if ( ok > 0 ) call abort_mpi ('error allocate read_restart')

    v   = 0.0_dp
    ha  = 0.0_dp
    T   = 0.0_dp
    W_i = 0.0_dp
    cp  = 0.0_dp
    mu_SGS = 0.0_dp

    write ( number_ascii , format_restart ) inp % number_restart
    adr_file = trim (dir_parent) //  trim (inp % name_restart) // '_' //  trim (number_ascii) // '.h5'

    if ( rank == rank_default ) write (*,*) 'reading file ' , trim ( adr_file )

    call get_variable_names ( var_names )

    call hdf_open_file ( adr_file, readonly=.true.)

    call hdf_open_group ("/consvar")

    do i = 1, niv
        call hdf_read_r3d_array ( var_names(i) , v (sx:ex,sy:ey,sz:ez,i), subdom_size , subdom_offset)
    enddo

    do i = 1 , nrv
        call hdf_read_r3d_array ( thd % NameSpecy(i) , v (sx:ex,sy:ey,sz:ez,niv+i), subdom_size , subdom_offset)
    enddo

    call hdf_close_group()

    call hdf_close_file()

    call comm_cons (v)
    call prim_inv_var_err_mach ( 0 , thd , v , W_i , T , cp , ha )
    call upd_extrapolation (v)
    call timestep ( inp , adi , thd , dx_i , dy_i , dz_i , cp , W_i , T , v , mu_SGS , 1.0_dp , dt_4 , dtmin )
    call upd_boundaries ( inp , adi , thd , time , dtmin , x , y , z , &
        dx_i , dy_i , dz_i , T , W_i , cp , ha , ite , v )
    do l = 1 , ndim+ndim
        call prim_inv_var_err_mach ( face_domain (l) , thd , v , W_i , T , cp , ha )
    end do

    deallocate ( mu_SGS )


end subroutine read_restart_hdf5



end module IO
