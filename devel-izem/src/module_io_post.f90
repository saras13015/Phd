module IO_post

    use parameters
    use input
    use adim  
    use parallel
    use variables
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



    subroutine post_get_variable_names ( var_names )

        character (len_default) , intent(inout) :: var_names (nv)

        var_names(1)   = "rho"
        var_names(2:4) = (/ "rhoU" , "rhoV" , "rhoW"/)
        var_names(5)   = "rhoE"

    end subroutine post_get_variable_names


    subroutine read_restart_hdf5 ( ifile , inp , thd , adi , dx_i , dy_i , dz_i , rweg , iweg , inst )


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
    integer (ip)                                                   :: ok , i , j , k , l
    integer (ip)                                                   :: myrec , px , py , pz , ix , &
    & fx , iy , fy , iz , fz

    character (len_default)                                        :: number_ascii , rank_ascii , &
    & plane_ascii , dir_stat , file_type , adr_file
    real (dp)                                                      :: rho_i , rho_cgs , rho_phys , T_phys , &
    & W_cgs , W_phys , wrk
    real (dp)                                                      :: yo2 , zm , zp2 , zp2max , xb , dm_SGS , &
    & taumix , Sc_SGS_i , rvd0
    real (dp)                                                      :: beta1 , beta2 , deltazj
    real (dp) , dimension (nrv+npv+nvv)                            :: Ya
    real (dp) , dimension (nreac)                                  :: kb , kf , qi , W_scal_i
    real (dp) , dimension (:,:,:,:) , allocatable                  :: Xa , ha
    real (dp) , dimension (:,:,:)   , allocatable                  :: delta_2
    real (dp) , parameter                                          :: eps2 = 1.0e-2_dp , eps3 = 1.0e-3_dp , eps4 = 1.0


        character (len_default)  :: var_names (nv)



        write ( number_ascii , format_restart ) ifile


        adr_file = trim (dir_restart) // trim (file_restart) // '_' // trim (number_ascii) // '.h5'

        ! adr_file = trim (dir_parent) // trim ( dir_restart ) //  trim ( name)



        if ( rank == rank_default ) write (*,*) '    reading previous restart file: ' , trim ( adr_file )

        call post_get_variable_names (var_names)

        call hdf_open_file ( adr_file, readonly=.true.)

        call hdf_open_group ("/consvar")

        do i = 1, nv
            call hdf_read_r3d_array ( var_names(i) , inst % u (sx:ex,sy:ey,sz:ez,i) , subdom_size , subdom_offset)
        enddo

        do i = 1 , nrv
            call hdf_read_r3d_array ( thd % NameSpecy(i) , inst % u  (sx:ex,sy:ey,sz:ez,niv+i) , subdom_size , subdom_offset)
        enddo



        call hdf_close_group()

        ! call hdf_open_group ("/simulation")

        ! call hdf_read_r1d_array ( "Time", time , (/1/) )

        ! call hdf_close_group() 

        call hdf_close_file()


    ! verifying the fields
    do k = sz , ez
       do j = sy , ey
          do i = sx , ex
             if ( inst % u (i,j,k,1) <= 0.0_dp ) then
                write (*,'(A,3I7,1PE16.8)') &
                   'error: density is either zero or negative (i,j,k,rho) = ' , &
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
   

       call molefrac ( 0 , thd , inst % W_i , inst % u , Xa )
       call prim_vis_vars ( 0 , thd , inst % u , inst % W_i , inst % T , &
                            Xa , inst % dm , inst % mu , inst % ct )


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
    endif

    deallocate ( Xa , ha )


    end subroutine read_restart_hdf5

end module IO_post