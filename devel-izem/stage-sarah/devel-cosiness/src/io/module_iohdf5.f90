!------------------------------------------------------------------------------
! module: iohdf5
!------------------------------------------------------------------------------
!> \brief 
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)

module IOHDF5

    use parameters
    use input
    use parallel
    use hdf5
    use h5lt
    implicit none

    ! percision common parameters
    integer,parameter::rper=dp
    integer,parameter::iper=ip
    integer,parameter::longlen=1024
    integer,parameter::shortlen=64

    ! for configurations
    character(len=50),private::keyword='#hdf5 io#'

    ! hdf5 file ids and sizes
    integer(hid_t)::file_id       ! file identifier
    integer(hid_t)::dset_id       ! dataset identifier
    integer(hid_t)::filespace     ! dataspace identifier in file
    integer(hid_t)::memspace      ! dataspace identifier in memory
    integer(hid_t)::plist_id      ! property list identifier
    integer(hid_t)::dist_id       ! property list identifiler for dataset

    ! 3d datarank
    integer,private,parameter::rank3=3            ! dataset rank
    integer(hsize_t),dimension(rank3)::dimsf      ! dataset dimensions

    ! in the file.
    integer(hsize_t),dimension(rank3)::lcl_dims          ! chunks dimensions
    integer(hssize_t),dimension(rank3)::offst_blck
    integer(hsize_t),dimension(rank3)::count
    integer(hsize_t),dimension(rank3)::blockh5
    integer(hsize_t),dimension(rank3)::stride
    integer(hsize_t),dimension(1)::t_dims,l_dims,l_count,l_block,l_stride,l_offst
    integer(hsize_t),dimension(1)::dimsx,dimsy,dimsz

    ! mpi commutor and info
    integer::h5comm,h5info,mpi_rank

    ! folder path 
    character(len=longlen)::path='./post01/'

    logical,private::l_write_init=.true.
    logical::l_h5_io_init=.true.
    logical::l_io_h5=.true.

    ! location of this cpu
    integer(kind=iper)::nx_proc=0,ny_proc=0,nz_proc=0

    ! compact treatments
    logical::l_h5_compact=.false.

    ! 4d data rank...
    integer,private,parameter::rank4=4          ! dataset rank
    integer(hsize_t),dimension(rank4)::dimsf4             ! 4d dataset dimensions
    integer(hssize_t),dimension(rank4)::offst_blck4        ! 4d offsets for 
    integer(hsize_t),dimension(rank4)::count4

    ! compact file dimensions
    integer::n_fields, &      ! number of 3d fields
    &        n_real,   &      ! number of reals
    &        n_int            ! number of ints 
    real(kind=ra),allocatable::data_fields(:,:,:,:),& ! data for 3d fields
    &                          data_real(:)            ! data for real   
    integer(kind=ia),allocatable::data_int(:)          ! data for int
    real(kind=ra)::time_tps(1)

    contains

    !> \brief Conf for hdf5
    !!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
    subroutine h5_io_conf(file_data)
        character(len=25),intent(in):: file_data
        integer::iunit
        character(len=40)::messs
        if(rank == rank_default) then
            call file_mark(file_data,trim(keyword),iunit)
            if(iunit>0) then
                read(iunit,*) l_io_h5
                read(iunit,*) l_h5_compact
                read(iunit,'(a)') messs
                path=trim(adjustl(messs))
                close(iunit)       
            end if
        else
            l_io_h5=.false.
            l_h5_compact=.false.       
        end if

#ifdef mpi
        call mpi_bcast(l_io_h5,1,mpi_logical,0,mpi_comm_world,ierr)
        call mpi_bcast(l_h5_compact,1,mpi_logical,0,mpi_comm_world,ierr)
        call mpi_bcast(path,longlen,mpi_character,0,mpi_comm_world,ierr)
#endif
        ! write(6,*) l_io_h5,l_h5_compact,my_rank,path

    end subroutine h5_io_conf


    !> \brief Init compact reading
    !!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
    subroutine h5_compact_init(nd,ns)
        integer,intent(in)::nd,ns
        logical,save::l_allo=.false.

        if(l_allo) return
        l_allo=.true. 

        call system("mkdir -p ./"//trim(adjustl(path))//'/compact')
        ! fields
        n_fields=3+4+nd+2*ns+1+ns+nd+1+1+2 ! rhomem(3)+ rho+u(nd)+,t,p,h,+ys,wdots(2*ns)+hr+dydt(ns)+dudt(nd)+dtdt+lvf+(pscalar+dpassivdt)
        allocate(data_fields(n_fields,count(1),count(2),count(3)))
        data_fields=0.d00
        ! reals
        n_real=1+3+9! p0lmn,rhomemtime,adis
        allocate(data_real(n_real));data_real=-1000.
        ! ints,
        n_int=3
        allocate(data_int(n_int))

        dimsf4(1)=n_fields
        dimsf4(2:4)=dimsf(1:3)
        offst_blck4(:)=0
        offst_blck4(2:4)=offst_blck(:)
        count4(1)=n_fields
        count4(2:4)=count(:)

    end subroutine h5_compact_init

    !> \brief Writing the 4d compact file
    !!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
    subroutine h5_write_4d_file(filename,tps)
        implicit none
        character(len=*),                    intent(in)::filename ! file name
        real(rper),                          intent(in)::tps    !
        real(rper)      ::data(1)
        integer:: rank1d
        integer(hsize_t)::dim1d(1)
        integer(hid_t)::file_id         ! file identifier 
        integer(hid_t)::plist_id        ! property list identifier 
        integer(size_t):: attr_dim      ! attribute shape 
        integer::error                  ! error flags

        !
        character(len=longlen)::pathfile,dsetname
        real(rper)::tic,toc


        tic=mpi_wtime()
        call system("mkdir -p ./"//trim(adjustl(path))//'/compact')
        ! initialize fortran predefined datatypes
        call h5open_f(error)
        ! setup file access property list with parallel i/o access.
        call h5pcreate_f(h5p_file_access_f,plist_id,error)
        call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)

        pathfile=trim(adjustl(path))//'/compact/'//trim(adjustl(filename))

        call h5fcreate_f(trim(pathfile),h5f_acc_trunc_f,file_id,error,access_prp=plist_id)
        call h5pclose_f(plist_id,error)

        ! attribute,time
        data(1)=tps
        attr_dim=1
        call h5ltset_attribute_double_f(file_id,'/','time',data,attr_dim,error)
        !call h5fcreate_f(trim(pathfile),h5f_acc_trunc_f,file_id,error)
        ! write float fields
        rank1d=1
        dim1d=n_real
        call h5screate_simple_f(rank1d,dim1d,filespace,error)

        ! create the dataset with default properties.
        dsetname="reals"
        call h5dcreate_f(file_id,dsetname,h5t_native_double,filespace,dset_id,error)

        ! create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 

        !call h5pset_dxpl_mpio_f(dlist_id,h5fd_mpio_independent_f,error)
        call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)

        ! write the dataset with first cpu
        call h5dwrite_f(dset_id,h5t_native_double,data_real,dim1d,error,xfer_prp=plist_id)

        ! close dataspaces.
        call h5sclose_f(filespace,error)

        ! close the dataset and property list.
        call h5dclose_f(dset_id,error)

        !-------------------------------------------------------------------------
        ! write int fields

        rank1d=1
        dim1d=n_int
        call h5screate_simple_f(rank1d,dim1d,filespace,error)

        ! create the dataset with default properties.
        dsetname="ints"
        call h5dcreate_f(file_id,dsetname,h5t_native_integer,filespace,dset_id,error)
        call h5dwrite_f(dset_id,h5t_native_integer,data_int,dim1d,error,xfer_prp=plist_id)

        ! close dataspaces.
        call h5sclose_f(filespace,error)
        ! close the dataset and property list.
        call h5dclose_f(dset_id,error)


        !-------------------------------------------------------------------------
        ! write 4d fields
        ! create the data space for the  dataset. 
        call h5screate_simple_f(rank4,dimsf4,filespace,error)
        !write (6,*) rank,dimsf,count
        ! create the dataset with default properties.
        dsetname="fields"
        call h5dcreate_f(file_id,dsetname,h5t_native_double,filespace,dset_id,error)
        call h5sclose_f(filespace,error)

        call h5screate_simple_f(rank4,count4,memspace,error) 
        ! select hyperslab in the file.
        call h5dget_space_f(dset_id,filespace,error)
        call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offst_blck4,count4,error)
        ! create property list for collective dataset write
        !call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 
        !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_independent_f,error)
        !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)
        ! write the dataset collectively.

        call h5dwrite_f(dset_id,h5t_native_double,data_fields,dimsf4,error,&
            &          file_space_id=filespace,mem_space_id=memspace,&
            &          xfer_prp=plist_id)
        ! write the dataset independently. 
        !call h5dwrite_f(dset_id,h5t_native_double,data,dimsf,error,&
        !     &          file_space_id=filespace,mem_space_id=memspace)
        !
        !!$    t_dims=1; tps(:)=localtps; attr_dim=1
        !!$    call h5ltset_attribute_double_f(file_id,dsetname,'time',tps,attr_dim,error)    
        ! close dataspaces.
        call h5sclose_f(filespace,error)
        call h5sclose_f(memspace,error)
        ! close the dataset and property list.
        call h5dclose_f(dset_id,error)
        call h5pclose_f(plist_id,error)
        call h5fclose_f(file_id,error)
        call h5close_f(error)
        toc=mpi_wtime()
        call frline("write compact,t",toc-tic)
    end subroutine h5_write_4d_file


    !> \brief Reading the 4d compact file
    !!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
    subroutine h5_read_4d_file(filename,tps)
        character(len=*),                    intent(in)::filename ! file name
        real(rper),                          intent(out):: tps    !
        real(rper)      ::data(1)
        integer:: rank1d
        integer(hsize_t)::dim1d(1),count4_read(4)
        integer(hid_t)::file_id       ! file identifier 
        !integer(hid_t)::dset_id       ! dataset identifier 

        integer(hid_t)::plist_id        ! property list identifier 
        integer(size_t):: attr_dim      ! attribute shape 
        integer::error                  ! error flags
     
        ! mpi definitions and calls.
        character(len=longlen)::pathfile,dsetname
        real(rper)::tic,toc

        tic=mpi_wtime()
        call system("mkdir -p ./"//trim(adjustl(path))//'/compact')
        ! initialize fortran predefined datatypes
        call h5open_f(error)
        ! setup file access property list with parallel i/o access.
        call h5pcreate_f(h5p_file_access_f,plist_id,error)
        call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)

        pathfile=trim(adjustl(path))//'/compact/'//trim(adjustl(filename))
        call h5fopen_f(trim(pathfile),&
            &         h5f_acc_rdonly_f,file_id,error)
        !call h5fcreate_f(trim(pathfile),h5f_acc_trunc_f,file_id,error,access_prp=plist_id)
        call h5pclose_f(plist_id,error)

        ! attribute,time

        attr_dim=1
        data=-100.d00
        call h5ltget_attribute_double_f(file_id,'/','time',data,error)
        tps=data(1)
        ! create property list for collective dataset write
        call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 
        !call h5pset_dxpl_mpio_f(dlist_id,h5fd_mpio_independent_f,error)
        call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)

        !call h5fcreate_f(trim(pathfile),h5f_acc_trunc_f,file_id,error)
        ! -----------------------------------------------------------------------
        ! read float fields
        rank1d=1
        dim1d=n_real

        dsetname="reals"


        call h5dopen_f(file_id,dsetname,dset_id,error)

        call h5dread_f(dset_id,h5t_native_double,data_real,dim1d,error)!,xfer_prp=plist_id)

        ! close the dataset and property list.
        call h5dclose_f(dset_id,error)

        !-------------------------------------------------------------------------
        ! read int fields

        rank1d=1
        dim1d=n_int
        ! create the dataset with default properties.
        dsetname="ints"
        call h5dopen_f(file_id,dsetname,dset_id,error)

        call h5dread_f(dset_id,h5t_native_integer,data_int,dim1d,error)!,xfer_prp=plist_id)

        ! close dataspaces.
        call h5dclose_f(dset_id,error)


        !-------------------------------------------------------------------------
        ! read 4d fields
        ! create the data space for the  dataset. 


        dsetname="fields"
        call h5dopen_f(file_id,dsetname,dset_id,error)
        call h5dget_space_f(dset_id,filespace,error)
        count4_read=count4
        call h5screate_simple_f(rank4,count4_read,memspace,error) 
        ! select hyperslab in the file.

        call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offst_blck4,count4_read,error)
        ! create property list for collective dataset write
        !call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 
        !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_independent_f,error)
        !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)
        ! write the dataset collectively.

        call h5dread_f(dset_id,h5t_native_double,data_fields(1:count4_read(1),:,:,:),count4_read,error,&
            &          file_space_id=filespace,mem_space_id=memspace)!,&
        !         &          xfer_prp=plist_id)
        ! write the dataset independently. 
        !call h5dwrite_f(dset_id,h5t_native_double,data,dimsf,error,&
        !     &          file_space_id=filespace,mem_space_id=memspace)
        !
        !!$    t_dims=1; tps(:)=localtps; attr_dim=1
        !!$    call h5ltset_attribute_double_f(file_id,dsetname,'time',tps,attr_dim,error)    
        ! close dataspaces.
        call h5sclose_f(filespace,error)
        call h5sclose_f(memspace,error)
        ! close the dataset and property list.
        call h5dclose_f(dset_id,error)
        call h5pclose_f(plist_id,error)
        call h5fclose_f(file_id,error)
        call h5close_f(error)
        !write(6,*) maxval(data_fields(38,:,:,:))
        toc=mpi_wtime()
        call frline("read hdf comp,t",toc-tic)

    end subroutine h5_read_4d_file



    !> \brief Conf root folder for hdf5
    !! inputs:
    !! folder: root folder for hdf
    !!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
    subroutine h5_conf_folder(folder)
        character(len=*),intent(in)::folder
        path=folder
    end subroutine h5_conf_folder

    !> \brief Init field writing
    !! inputs: global nx,ny,nz,local nx,ny,nz
    !!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)


    subroutine h5_initio_field( nx_h5,  ny_h5, nz_h5,    &
        &                       nx_h5_l,ny_h5_l,nz_h5_l,nd,ns)
    implicit none
    !
    integer(kind=iper),intent(in)::nx_h5_l  ,ny_h5_l  ,nz_h5_l
    integer(kind=iper),intent(in)::nx_h5    ,ny_h5    ,nz_h5,nd,ns
    !
    integer::ierr
    !
    dimsf   =(/nx_h5,ny_h5,nz_h5/)
    lcl_dims=(/nx_h5_l,ny_h5_l,nz_h5_l/)
    count   =(/nx_h5_l,ny_h5_l,nz_h5_l/)
    blockh5 =(/1,1,1/)
    stride  =(/1,1,1/)

    !
    nx_proc=coords(3)
    if(ndim>1) then
        ny_proc=coords(2)
        if(ndim>2) then
            nz_proc=coords(1)
        end if
    end if
    !
    offst_blck(1)=nx_proc*(nx_h5_l)
    offst_blck(2)=ny_proc*(ny_h5_l)
    offst_blck(3)=nz_proc*(nz_h5_l)

    !
    !mpi_rank=rang

    h5comm=mpi_comm_world
    !
    ! h5info
    h5info=mpi_info_null
    call mpi_info_create(h5info,ierr)
    !write(6,*) ierr
    call mpi_info_set(h5info,"romio_cb_write","automatic",ierr)
    !write(6,*) ierr
    call mpi_info_set(h5info,"romio_ds_write","disable",ierr)
    !write(6,*) ierr
    call mpi_info_set(h5info,"romio_cb_read","automatic",ierr)
    !write(6,*) ierr
    call mpi_info_set(h5info,"romio_ds_read","disable",ierr)
    !write(6,*) ierr
    call mpi_info_set(h5info,"striping_unit","67108864",ierr)
    call mpi_info_set(h5info,"striping_factor","60",ierr)
    call mpi_comm_rank(mpi_comm_world,mpi_rank,ierr)
    call system("mkdir -p ./"//trim(adjustl(path)))
    if(l_h5_compact) then
        call h5_compact_init(nd,ns)    
    end if


    !write(6,*) lcl_dims
end subroutine h5_initio_field

! -----------------------------------------------------------------------
! create a file
! inputs
!       filename: file name of h5
!       foldername: foldername inside root
!       tps:        physical time,to writen as attri
subroutine h5_create_file(filename,foldername,tps)

    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::foldername
    real(rper),                          intent(in)::tps    !
    real(rper)      ::data(1)
    integer(hid_t)::file_id         ! file identifier 
    integer(hid_t)::plist_id        ! property list identifier 
    integer(size_t):: attr_dim      ! attribute shape 
    integer::error                  ! error flags
    
    character(len=longlen)::pathfile
    !return
    call system("mkdir -p ./"//trim(adjustl(path))//'/'//trim(adjustl(foldername)))

    if(l_write_init) then
        l_write_init=.false.
    end if
    !h5comm=mpi_comm_world
    !h5info=mpi_info_null
    !
    ! initialize fortran predefined datatypes
    call h5open_f(error)
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    call h5fcreate_f(trim(pathfile),h5f_acc_trunc_f,file_id,error,access_prp=plist_id)
    !call h5fcreate_f(trim(pathfile),h5f_acc_trunc_f,file_id,error)

    data(1)=tps
    attr_dim=1
    call h5ltset_attribute_double_f(file_id,'/','time',data,attr_dim,error)
    call h5pclose_f(plist_id,error)
    call h5fclose_f(file_id,error)
    call h5close_f(error)


end subroutine h5_create_file


! ------------------------------------------------------------------------
! write a 3d field to hdf files
! input:
!        filename
!        foldername
!        dsetname
!        nxl,nyl,nzl : local shape
!        data: data block to write,shape(nxl,nyl,nzl)
! ------------------------------------------------------------------------
subroutine h5_write_field(filename,foldername,dsetname,nxl,nyl,nzl,data)
    implicit none
    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::dsetname
    character(len=*),                    intent(in)::foldername
    integer(kind=iper),                  intent(in)::nxl,nyl,nzl ! local block size
    real(kind=rper),dimension(nxl,nyl,nzl),intent(in)::data! data to write allocatable
    !
    integer(hid_t)::file_id       ! file identifier 
    integer(hid_t)::dset_id       ! dataset identifier 
    integer(hid_t)::filespace     ! dataspace identifier in file 
    integer(hid_t)::memspace      ! dataspace identifier in memory
    integer(hid_t)::plist_id      ! property list identifier 
    integer::error                ! error flags
    character(len=longlen)::pathfile


    ! initialize fortran predefined datatypes
    call h5open_f(error) 
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    ! create the file collectively.
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    ! try if exist a file with acc_rdwr
    call h5fopen_f(pathfile,h5f_acc_rdwr_f,file_id,error,access_prp=plist_id)
    !    write(6,*) error

    call h5pclose_f(plist_id,error)
    ! create the data space for the  dataset. 
    call h5screate_simple_f(rank,dimsf,filespace,error)
    !write (6,*) rank,dimsf,count
    ! create the dataset with default properties.
    call h5dcreate_f(file_id,dsetname,h5t_native_double,filespace,dset_id,error)
    call h5sclose_f(filespace,error)

    call h5screate_simple_f(rank,count,memspace,error) 
    ! select hyperslab in the file.
    call h5dget_space_f(dset_id,filespace,error)
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offst_blck,count,error)
    ! create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 
    !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_independent_f,error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)
    ! write the dataset collectively. 
    call h5dwrite_f(dset_id,h5t_native_double,data,dimsf,error,&
        file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
    ! write the dataset independently. 
    !call h5dwrite_f(dset_id,h5t_native_double,data,dimsf,error,&
    !     &          file_space_id=filespace,mem_space_id=memspace)
    !
    !!$    t_dims=1; tps(:)=localtps; attr_dim=1
    !!$    call h5ltset_attribute_double_f(file_id,dsetname,'time',tps,attr_dim,error)    
    ! close dataspaces.
    call h5sclose_f(filespace,error)
    call h5sclose_f(memspace,error)
    ! close the dataset and property list.
    call h5dclose_f(dset_id,error)
    call h5pclose_f(plist_id,error)
    ! close the file
    call h5fclose_f(file_id,error)
    call h5close_f(error)
end subroutine h5_write_field


! ------------------------------------------------------------------------
! write a 1d array to hdf files
! input:
!        filename
!        foldername
!        dsetname
!        n: data size
!        data: data block to write,shape(n)
! ------------------------------------------------------------------------
subroutine h5_write_real(filename,foldername,dsetname,n,data)
    implicit none
    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::dsetname
    character(len=*),                    intent(in)::foldername
    integer(kind=iper),                  intent(in)::n ! local block size
    real(kind=rper),dimension(n),       intent(in)::data! data to write allocatable
    !
    integer(hid_t)::file_id           ! file identifier 
    integer(hid_t)::dset_id           ! dataset identifier 
    integer(hid_t)::filespace         ! dataspace identifier in file 
    integer(hid_t)::plist_id          ! property list identifier 
    integer::error                    ! error flags
    character(len=longlen)::pathfile
    integer:: rank1d
    integer(hsize_t)::dim1d(1)
    !return
    !
    !
    !call get_coords(rang,ix,iy,iz)
    !
    ! initialize fortran predefined datatypes
    call h5open_f(error) 
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    ! create the file collectively.
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    ! try if exist a file with acc_rdwr
    call h5fopen_f(pathfile,h5f_acc_rdwr_f,file_id,error,access_prp=plist_id)


    call h5pclose_f(plist_id,error)
    ! create the data space for the  dataset.
    rank1d=1
    dim1d=n
    call h5screate_simple_f(rank1d,dim1d,filespace,error)
    ! create the dataset with default properties.
    call h5dcreate_f(file_id,dsetname,h5t_native_double,filespace,dset_id,error)
    ! create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 
    !call h5pset_dxpl_mpio_f(dlist_id,h5fd_mpio_independent_f,error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)
    ! write the dataset with first cpu

    call h5dwrite_f(dset_id,h5t_native_double,data,dim1d,error,xfer_prp=plist_id)
    ! write the dataset independently. 
    !call h5dwrite_f(dset_id,h5t_native_double,data,dimsf,error,&
    !     &          file_space_id=filespace,mem_space_id=memspace)
    !
    !!$    t_dims=1; tps(:)=localtps; attr_dim=1
    !!$    call h5ltset_attribute_double_f(file_id,dsetname,'time',tps,attr_dim,error)    
    ! close dataspaces.
    call h5sclose_f(filespace,error)
    ! close the dataset and property list.
    call h5dclose_f(dset_id,error)
    call h5pclose_f(plist_id,error)
    ! close the file
    call h5fclose_f(file_id,error)
    call h5close_f(error)
end subroutine h5_write_real





! ------------------------------------------------------------------------
! write  integer array to hdf files
! input:
!        filename
!        foldername
!        dsetname
!        n: data size
!        data: data block to write,shape(n)
! ------------------------------------------------------------------------
subroutine h5_write_int(filename,foldername,dsetname,n,data)
    implicit none
    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::dsetname
    character(len=*),                    intent(in)::foldername
    integer(kind=iper),                  intent(in)::n ! local block size
    integer(kind=iper),dimension(n),    intent(in)::data! data to write allocatable
    !
    integer(hid_t)::file_id       ! file identifier 
    integer(hid_t)::dset_id       ! dataset identifier 
    integer(hid_t)::filespace     ! dataspace identifier in file 
    integer(hid_t)::plist_id      ! property list identifier 
    integer::error                ! error flags

    character(len=longlen)::pathfile
    integer:: rank1d
    integer(hsize_t)::dim1d(1)
    !return
    !
    !
    !call get_coords(rang,ix,iy,iz)
    !
    ! initialize fortran predefined datatypes
    call h5open_f(error) 
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    ! create the file collectively.
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    ! try if exist a file with acc_rdwr
    call h5fopen_f(pathfile,h5f_acc_rdwr_f,file_id,error,access_prp=plist_id)
    !    write(6,*) error

    call h5pclose_f(plist_id,error)
    ! create the data space for the  dataset.
    rank1d=1
    dim1d=n
    call h5screate_simple_f(rank1d,dim1d,filespace,error)
    ! create the dataset with default properties.
    call h5dcreate_f(file_id,dsetname,h5t_native_integer,filespace,dset_id,error)
    ! create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 
    !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_independent_f,error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)
    ! write the dataset collectively.

    call h5dwrite_f(dset_id,h5t_native_integer,data,dim1d,error,&
        &          xfer_prp=plist_id)

    ! write the dataset independently. 
    !call h5dwrite_f(dset_id,h5t_native_double,data,dimsf,error,&
    !     &          file_space_id=filespace,mem_space_id=memspace)
    !
    !!$    t_dims=1; tps(:)=localtps; attr_dim=1
    !!$    call h5ltset_attribute_double_f(file_id,dsetname,'time',tps,attr_dim,error)    
    ! close dataspaces.
    call h5sclose_f(filespace,error)
    ! close the dataset and property list.
    call h5dclose_f(dset_id,error)
    call h5pclose_f(plist_id,error)
    ! close the file
    call h5fclose_f(file_id,error)
    call h5close_f(error)
end subroutine h5_write_int



! ------------------------------------------------------------------------
! read a 1d array to hdf files
! input:
!        filename
!        foldername
!        dsetname
!        n: data size
!        data: data block to write,shape(n)
! ------------------------------------------------------------------------
subroutine h5_read_real(filename,foldername,dsetname,n,data)
    implicit none
    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::dsetname
    character(len=*),                    intent(in)::foldername
    integer(kind=iper),                  intent(in)::n ! local block size
    real(kind=rper),dimension(n),       intent(out)::data! data to write allocatable
    !
    integer(hid_t)::file_id       ! file identifier 
    integer(hid_t)::dset_id       ! dataset identifier 
    integer(hid_t)::plist_id      ! property list identifier 
    integer::error                ! error flags

    character(len=longlen)::pathfile
    integer:: rank1d
    integer(hsize_t)::dim1d(1)
    !
    !
    !call get_coords(rang,ix,iy,iz)
    !
    ! initialize fortran predefined datatypes
    call h5open_f(error) 
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    ! create the file collectively.
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    ! try if exist a file with acc_rdwr
    call h5fopen_f(pathfile,h5f_acc_rdonly_f,file_id,error)!,access_prp=plist_id)
    !    write(6,*) error

    call h5pclose_f(plist_id,error)
    ! create the data space for the  dataset.
    rank1d=1
    dim1d=n
    ! create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 
    !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_independent_f,error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)
    ! write the dataset collectively.
    call h5dopen_f(file_id,dsetname,dset_id,error)

    call h5dread_f(dset_id,h5t_native_double,data,dim1d,error)!,&
    !         &          xfer_prp=plist_id)
    ! write the dataset independently. 
    !call h5dwrite_f(dset_id,h5t_native_double,data,dimsf,error,&
    !     &          file_space_id=filespace,mem_space_id=memspace)
    !
    !!$    t_dims=1; tps(:)=localtps; attr_dim=1
    !!$    call h5ltset_attribute_double_f(file_id,dsetname,'time',tps,attr_dim,error)    

    ! close the dataset and property list.
    call h5dclose_f(dset_id,error)
    call h5pclose_f(plist_id,error)
    ! close the file
    call h5fclose_f(file_id,error)
    call h5close_f(error)
end subroutine h5_read_real

! ------------------------------------------------------------------------
! read a 1d integer array to hdf files
! input:
!        filename
!        foldername
!        dsetname
!        n: data size
!        data: data block to write,shape(n)
! ------------------------------------------------------------------------
subroutine h5_read_int(filename,foldername,dsetname,n,data)
    implicit none
    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::dsetname
    character(len=*),                    intent(in)::foldername
    integer(kind=iper),                  intent(in)::n ! local block size
    integer(kind=iper),dimension(n),    intent(out)::data! data to write allocatable
    !
    integer(hid_t)::file_id       ! file identifier 
    integer(hid_t)::dset_id       ! dataset identifier 
    integer(hid_t)::plist_id      ! property list identifier 
    integer::error                ! error flags
    character(len=longlen)::pathfile
    integer:: rank1d
    integer(hsize_t)::dim1d(1)
    !
    !
    !call get_coords(rang,ix,iy,iz)
    !
    ! initialize fortran predefined datatypes
    call h5open_f(error) 
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    ! create the file collectively.
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    ! try if exist a file with acc_rdwr
    call h5fopen_f(pathfile,h5f_acc_rdonly_f,file_id,error)!,access_prp=plist_id)
    !    write(6,*) error

    call h5pclose_f(plist_id,error)
    ! create the data space for the  dataset.
    rank1d=1
    dim1d=n
    ! create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 
    !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_independent_f,error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)
    ! write the dataset collectively.
    call h5dopen_f(file_id,dsetname,dset_id,error)

    call h5dread_f(dset_id,h5t_native_integer,data,dim1d,error)!,&
    !         &          xfer_prp=plist_id)
    ! write the dataset independently. 
    !call h5dwrite_f(dset_id,h5t_native_double,data,dimsf,error,&
    !     &          file_space_id=filespace,mem_space_id=memspace)
    !
    !!$    t_dims=1; tps(:)=localtps; attr_dim=1
    !!$    call h5ltset_attribute_double_f(file_id,dsetname,'time',tps,attr_dim,error)    

    ! close the dataset and property list.
    call h5dclose_f(dset_id,error)
    call h5pclose_f(plist_id,error)
    ! close the file
    call h5fclose_f(file_id,error)
    call h5close_f(error)
end subroutine h5_read_int

! ------------------------------------------------------------
! write a real attribute
! inputs:
!        filename
!        foldername: foldername in root
!        dsetname:   name of the attribute
!        n:          size of the attrib
!        data:       real array to write size=(n)
subroutine h5_write_attri_r(filename,foldername,dsetname,n,data)
    implicit none

    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::dsetname
    character(len=*),                    intent(in)::foldername
    integer(kind=iper),                  intent(in)::n! data size
    real(kind=rper),dimension(n),intent(in)::data! attrib data to write 
    !
    integer(hid_t)::file_id       ! file identifier 
    integer(hid_t)::plist_id      ! property list identifier 
    integer(size_t):: attr_dim    ! attribute shape 
    integer::error                ! error flags
    character(len=longlen)::pathfile

    ! initialize fortran predefined datatypes
    call h5open_f(error) 
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    ! create the file collectively.
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    ! try if exist a file with acc_rdwr
    call h5fopen_f(pathfile,h5f_acc_rdwr_f,file_id,error,access_prp=plist_id)
    !    write(6,*) error
    call h5pclose_f(plist_id,error)

    !     &          file_space_id=filespace,mem_space_id=memspace)
    !
    t_dims=1; attr_dim=n

    call h5ltset_attribute_double_f(file_id,'/',trim(adjustl(dsetname)),data,attr_dim,error)    
    ! close dataspaces.
    ! close the file
    call h5fclose_f(file_id,error)
    call h5close_f(error)
end subroutine h5_write_attri_r

! ------------------------------------------------------------
! read a real attribute
! inputs:
!        filename
!        foldername: foldername in root
!        dsetname:   name of the attribute
!        n:          size of the attrib
! outputs:
!        data:       real array to write size=(n)
subroutine h5_read_attri_r(filename,foldername,dsetname,n,data)

    use hdf5
    use h5lt
    !use mpi_mod,only : cmp
    implicit none

    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::dsetname
    character(len=*),                    intent(in)::foldername
    integer(kind=iper),                  intent(in)::n! data size
    real(kind=rper),dimension(n),intent(out)::data! attrib data to write 
    !
    integer(hid_t)::file_id       ! file identifier 
    integer(hid_t)::plist_id      ! property list identifier 
    integer(size_t):: attr_dim    ! attribute shape 
    integer::error                ! error flags
    character(len=longlen)::pathfile

    ! initialize fortran predefined datatypes
    call h5open_f(error) 
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    ! create the file collectively.
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    ! try if exist a file with acc_rdwr
    call h5fopen_f(pathfile,h5f_acc_rdonly_f,file_id,error)!,access_prp=plist_id)
    !     &          file_space_id=filespace,mem_space_id=memspace)
    !
    t_dims=1; attr_dim=n
    call h5ltget_attribute_double_f(file_id,'/',trim(adjustl(dsetname)),data,error)    
    ! close dataspaces.

    ! close the file
    call h5pclose_f(plist_id,error)
    call h5fclose_f(file_id,error)
    call h5close_f(error)
end subroutine h5_read_attri_r


! ------------------------------------------------------------
! read a int attribute
! inputs:
!        filename
!        foldername: foldername in root
!        dsetname:   name of the attribute
!        n:          size of the attrib
! outputs:
!        data:       real array to write size=(n)
subroutine h5_read_attri_i(filename,foldername,dsetname,n,data)
    implicit none

    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::dsetname
    character(len=*),                    intent(in)::foldername
    integer(kind=iper),                  intent(in)::n! data size
    integer(kind=iper),dimension(n),intent(out)::data! attrib data to write 
    !
    integer(hid_t)::file_id         ! file identifier 
    integer(hid_t)::plist_id        ! property list identifier 
    integer(size_t):: attr_dim      ! attribute shape 
    integer::error                  ! error flags
    character(len=longlen)::pathfile

    ! initialize fortran predefined datatypes
    call h5open_f(error) 
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    ! create the file collectively.
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    ! try if exist a file with acc_rdwr
    call h5fopen_f(pathfile,h5f_acc_rdonly_f,file_id,error)!,access_prp=plist_id)
    !    write(6,*) error
    !     &          file_space_id=filespace,mem_space_id=memspace)
    !
    t_dims=1; attr_dim=n

    call h5ltget_attribute_int_f(file_id,'/',trim(adjustl(dsetname)),data,error)
    ! close dataspaces.

    ! close the file
    call h5pclose_f(plist_id,error)
    call h5fclose_f(file_id,error)
    call h5close_f(error)
end subroutine h5_read_attri_i


! -----------------------------------------------------------------
! read hdf5
! input:
!       filename: file name from root path
!       dsetname: dataset name
!       nxl,nyl,nzl: local data size
!       data: output data
subroutine h5_read_field(filename,foldername,dsetname,nxl,nyl,nzl,data)
    use filetools
    !
    implicit none
    !
    character(len=*),intent(in)::filename,foldername
    character(len=*),intent(in)::dsetname
    integer(kind=iper),intent(in)::nxl,nyl,nzl
    !
    ! in the file.
    real(kind=rper),intent(out)::data (nxl,nyl,nzl)  ! data to read
    integer::error                                   ! error flags
    real::tic
    !
    !----------------------------------------------------------------------
    ! initialize hdf5 library and fortran interfaces.
    ! - here we open the fortran interface...
    !----------------------------------------------------------------------

    call h5open_f(error) 

    !----------------------------------------------------------------------
    ! file access (existing file-read only access)
    ! - but first i generate the parallel file access property
    !----------------------------------------------------------------------
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    call h5fopen_f(trim(adjustl(path))//'/'//trim(adjustl(foldername))//&
        &         trim(adjustl(filename)),&
        &         h5f_acc_rdonly_f,file_id,error)!,access_prp=plist_id)
    call h5pclose_f(plist_id,error)

    !----------------------------------------------------------------------
    !   hdf5 data access
    !----------------------------------------------------------------------

    !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_independent_f,error)  
    call h5dopen_f(file_id,dsetname,dset_id,error)


! ----------------------------------------------------------------------
    ! hdf5 dataset access
    ! ---------------------------------------------------------------------

    call h5dget_space_f(dset_id,filespace,error)
    call h5screate_simple_f(3,lcl_dims,memspace,error)

    !----------------------------------------------------------------------
    ! select hyperslab in the file.
    ! select this subset of the variable's space in the file 
    !----------------------------------------------------------------------

    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offst_blck,count,error)

    !----------------------------------------------------------------------
    ! read the data.  we're writing it from memory,where it is saved 
    ! in native_integer format 
    !----------------------------------------------------------------------
    call cpu_time(tic)
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_independent_f,error)  
    call h5dread_f(dset_id,h5t_native_double,data,lcl_dims,error, &
        &         file_space_id=filespace,mem_space_id=memspace)!,&
    !         &         xfer_prp=plist_id)
    if (error /= 0) then
        write(6,*) 'hdfreaderror'
        stop
    endif
    !----------------------------------------------------------------------
    ! close the opened material
    !----------------------------------------------------------------------
    call h5sclose_f(memspace,error)
    call h5dclose_f(dset_id,error) ! close the dataset.
    call h5pclose_f(plist_id,error) ! close dataset property list
    call h5sclose_f(filespace,error)  ! close dataspaces.
    call h5fclose_f(file_id,error) ! close the file.
    call h5close_f(error) ! close fortran interfaces and hdf5 library.
    !
end subroutine h5_read_field

! ------------------------------------------------------------
! get a hdf file name from simulation name,var name and inumfile
!
!
! -------------------------------------------------------------
subroutine h5_file_name(simulname,varname,inumfile,flname)
    implicit none
    character(len=*),intent(in)::simulname,varname !max len(6)
    integer(iper),intent(in)::inumfile
    character(len=longlen),intent(out)::flname  
    character(len=6)::chinum
    !
    write(chinum,'(i6.6)') inumfile

    flname=trim(adjustl(simulname))//"."//trim(adjustl(varname))//'.'//&
    &   trim(adjustl(chinum))//".h5"
    !write(6,*) trim(flname)
    !
end subroutine h5_file_name

! ------------------------------------------------------------------------
! write a 4d field to hdf files
! input:
!        filename
!        foldername
!        dsetname
!        nxl,nyl,nzl : local shape
!        data: data block to write,shape(nd,nxl,nyl,nzl)
! ------------------------------------------------------------------------
subroutine h5_write_field_4d(filename,foldername,dsetname,nd,nxl,nyl,nzl,data)
    use hdf5
    use h5lt
    use mpi
    !use mpi_mod,only : cmp
    implicit none
    character(len=*),                    intent(in)::filename ! file name
    character(len=*),                    intent(in)::dsetname
    character(len=*),                    intent(in)::foldername
    integer(kind=iper),                  intent(in)::nd,nxl,nyl,nzl ! local block size
    real(kind=rper),dimension(nd,nxl,nyl,nzl),intent(in)::data! data to write allocatable
    !
    integer(hid_t)::file_id       ! file identifier 
    integer(hid_t)::dset_id       ! dataset identifier 
    integer(hid_t)::filespace     ! dataspace identifier in file 
    integer(hid_t)::memspace      ! dataspace identifier in memory
    integer(hid_t)::plist_id      ! property list identifier 
    integer::error                ! error flags

    real(dp)::tic,toc
    character(len=longlen)::pathfile
    !
    tic=mpi_wtime()
    !call get_coords(rang,ix,iy,iz)
    !
    ! initialize fortran predefined datatypes
    call h5open_f(error) 
    ! setup file access property list with parallel i/o access.
    call h5pcreate_f(h5p_file_access_f,plist_id,error)
    call h5pset_fapl_mpio_f(plist_id,h5comm,h5info,error)
    ! create the file collectively.
    pathfile=trim(adjustl(path))//'/'//trim(adjustl(foldername))//'/'//trim(adjustl(filename))

    ! try if exist a file with acc_rdwr
    call h5fopen_f(pathfile,h5f_acc_rdwr_f,file_id,error,access_prp=plist_id)
    !    write(6,*) error

    call h5pclose_f(plist_id,error)
    ! create the data space for the  dataset. 
    call h5screate_simple_f(rank4,dimsf4,filespace,error)
    !write (6,*) rank,dimsf,count
    ! create the dataset with default properties.
    call h5dcreate_f(file_id,dsetname,h5t_native_double,filespace,dset_id,error)
    call h5sclose_f(filespace,error)

    call h5screate_simple_f(rank4,count4,memspace,error) 
    ! select hyperslab in the file.
    call h5dget_space_f(dset_id,filespace,error)
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offst_blck4,count4,error)
    ! create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,error) 
    !call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_independent_f,error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,error)
    ! write the dataset collectively.

    call h5dwrite_f(dset_id,h5t_native_double,data,dimsf4,error,&
        &          file_space_id=filespace,mem_space_id=memspace,&
        &          xfer_prp=plist_id)
    ! write the dataset independently. 
    !call h5dwrite_f(dset_id,h5t_native_double,data,dimsf,error,&
    !     &          file_space_id=filespace,mem_space_id=memspace)
    !
    !!$    t_dims=1; tps(:)=localtps; attr_dim=1
    !!$    call h5ltset_attribute_double_f(file_id,dsetname,'time',tps,attr_dim,error)    
    ! close dataspaces.
    call h5sclose_f(filespace,error)
    call h5sclose_f(memspace,error)
    ! close the dataset and property list.
    call h5dclose_f(dset_id,error)
    call h5pclose_f(plist_id,error)
    ! close the file
    call h5fclose_f(file_id,error)
    call h5close_f(error)
    toc=mpi_wtime()
    !call frline("write 4d",real(toc-tic,kind=rper))
end subroutine h5_write_field_4d

!-------------------------------------------------------|
!-------------------------------------------------------|
!- Mark detection in a data file
!-------------------------------------------------------|
!-------------------------------------------------------|
! iflag = -1 : file not found
! iflag = -2 : marker not found
! iflag > 0  : unit number where to read the data
!-------------------------------------------------------|
subroutine file_mark(filename,mark,flag)
    implicit none
    character (len=*), intent(in)     :: filename,mark
    integer(IA), intent(out)          :: flag
    !-
    character (len=50) :: readmark
    integer(IA)        :: err
    close(unit_data_file)
    call exist_file(filename)
    open(unit=unit_data_file,file=filename, form='formatted', status='old',iostat=err)
    if(err /= 0) then
        flag = -1
    else
        readmark = ''
        do while((readmark /= '# END FILE #').and.(readmark /= mark))
            read (unit_data_file,'(a50)') readmark
        enddo
        if (readmark == mark) then
            read (unit_data_file,'(a50)')
            flag    = unit_data_file
        elseif (readmark == '# END FILE #') then
            if (rank==0) write(iu_out,*) 'Unabled to find flag', '"' // trim(adjustl(mark)) // '" in ' // filename
            flag    = -2
            close(unit_data_file)
        else
            write(iu_out , *) 'SECURITY : '// filename // ' must contain # END FILE #'
            stop
            ! call stopline('SECURITY : '// filename // ' must contain # END FILE #')
        endif
    endif
    !-
end subroutine file_mark

!-------------------------------------------------------|
subroutine exist_file(filename)
    character (len=*), intent(in)     :: filename
    logical                           :: ex
    inquire(file = fileName, exist=ex)
    if (ex.eqv..FALSE.) then
        write(iu_out,*) char(27)//'[5m'//char(27)//'[1m'//               &
        'no such file or directory : ', filename//char(27)//'[m'
        stop
    endif
end subroutine exist_file

subroutine frline(ch,f)
    character(len=*), intent(in) :: ch
    real(dp), intent(in) :: f
    !-
    if (rank == rank_default) then
        130    format('x...',a25,' : ',e15.6)
        write(iu_out,130) ch,f
    endif
end subroutine frline


!> \brief Write the global grid to hdf file for posttreatment
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)


! subroutine Grid_write_hdf5(Grid)
!     type (inp_grid), intent(in)              :: grid

!     character(len=Shortlen)::foldername,dsetname
!     character(len=Longlen) :: HDF5FileName

!     ! should be put in other places
!     real(dp) ::Tps,Attri1(1)
!     integer(ip) :: N
!     Attri1=0.d00

!     Tps=0.d00
!     HDF5filename="mesh.h5"
!     foldername="../"
!     call H5_create_file(HDF5filename,foldername,tps)

!     dsetname="xcoords"
!     call H5_write_real(HDF5filename,foldername,dsetname,ntx,grid%xt)

!     dsetname="ycoords"
!     if(ndim>1) then
!         call H5_write_real(HDF5filename,foldername,dsetname,nty,grid%yt)

!     else
!         N=1
!         call H5_write_real(HDF5filename,foldername,dsetname,1,Attri1)
!     end if

!     dsetname="zcoords"
!     if(ndim>2) then
!         call H5_write_real(HDF5filename,foldername,dsetname,ntz,grid%zt)
!     else
!         call H5_write_real(HDF5filename,foldername,dsetname,1,Attri1)
!     end if
! end subroutine Grid_write_hdf5



end module iohdf5
