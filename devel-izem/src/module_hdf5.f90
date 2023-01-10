module iohdf5

!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

   !------------------------------------------------------------------------------------   
   use hdf5
   use mpi, only:MPI_INFO_NULL  
   use parameters,   only: dp, ip
   use parallel,     only: print_warning, error_and_exit,MPI_COMM_WORLD

   implicit none

   integer, parameter :: HDF_INT = 1
   integer, parameter :: HDF_REAL = 2
   integer, parameter :: HDF_CHAR = 3
   integer, parameter :: HDF_CURRENT = 4
   integer, parameter :: HDF_MAX_DIMS = 5
   integer, parameter :: HDF_LEN = 500

   integer, private :: hdf_ierror

   integer(HID_T), private :: hdf_file_id
   integer(HID_T), private :: hdf_group_id
   integer(HID_T), private :: hdf_dataset_id
   character(len=HDF_LEN), private :: hdf_dataset_name
   integer(HID_T), private :: hdf_dataset_properties_id
   integer, private :: hdf_compression
   integer(HID_T), private :: hdf_attribute_id
   character(len=HDF_LEN), private :: hdf_attribute_name
   integer(HID_T), private :: hdf_attribute_type_id
   integer(HID_T), private :: hdf_dataspace_id
   integer(HID_T), private :: hdf_dataspace_slab_id
   integer(HID_T), private :: hdf_datatype_id
   integer(HSIZE_T), private :: hdf_dataset_dim(HDF_MAX_DIMS)
   integer(HSIZE_T), private :: hdf_dataset_dim_max(HDF_MAX_DIMS)
   integer(HSIZE_T), private :: hdf_dataset_offset(HDF_MAX_DIMS)

   integer(HID_T), private :: hdf_parallel_access
   integer(HID_T), private :: hdf_parallel_io
   integer(HID_T), private :: hdf_filespace_id

   !===========================================================
   type hdf_group_t
   integer(HID_T) :: id
   character(len=HDF_LEN) :: name
   type(hdf_group_t), pointer :: parent
   type(hdf_group_t), pointer :: child
   end type hdf_group_t
   !===========================================================

   type(hdf_group_t), pointer :: hdf_first_group => NULL()

   type(hdf_group_t), pointer :: hdf_active_group => NULL()


   contains

   !================================================================
   !                        HDF INTERFACE
   !================================================================

   !=================================
   subroutine hdf_init()
      implicit none
      call h5open_f (hdf_ierror)
      ! call MPI_Info_create(hdf_mpi_info, hdf_ierror)
   end subroutine hdf_init
   !=================================

   !=================================
   subroutine hdf_destroy()
      implicit none
      call h5close_f (hdf_ierror)
      call hdf_error(hdf_ierror,"Unable to destroy the HDF session")
   end subroutine hdf_destroy
   !=================================

   !=================================
   subroutine hdf_error(error,themessage,is_fatal)
      implicit none
      ! ----------------------
      integer, intent(in) :: error
      character(len=*)   :: themessage
      logical, optional :: is_fatal
      ! ----------------------
      logical :: fatal
      ! ----------------------
      ! parameters
      fatal = .false.
      if (present(is_fatal)) then
         fatal = is_fatal
      end if
      ! check if there is an error
      if (error < 0) then
         call h5eprint_f (hdf_ierror)
         if (fatal) then
            call error_and_exit(themessage)
         else
            call print_warning(themessage)
         end if
      end if
   end subroutine hdf_error
   !=================================


   !================================================================
   !                        GROUP OPERATIONS
   !================================================================

   !=================================
   subroutine hdf_create_group(groupname)

      implicit none

      ! ----------------------
      character(len=*) :: groupname
      ! ----------------------
      integer :: ires
      character(HDF_LEN) :: message
      ! ----------------------
      ! create the group
      call h5gcreate_f (hdf_active_group%id,groupname,hdf_group_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to create the HDF group "//trim(groupname),is_fatal=.true.)

      ! allocate the group id
      allocate(hdf_active_group%child, stat=ires)
      if (ires /= 0) then
         write(message,'(a)')"[MEMORY error] Unable to allocate memory in hdf_create_group"
         call error_and_exit(message)
      endif
      hdf_active_group%child%id = hdf_group_id
      hdf_active_group%child%name = trim(groupname)
      NULLIFY(hdf_active_group%child%parent)
      NULLIFY(hdf_active_group%child%child)
      hdf_active_group%child%parent => hdf_active_group
      hdf_active_group => hdf_active_group%child
   end subroutine hdf_create_group
   !=================================

   !=================================
   subroutine hdf_open_group(groupname)

      implicit none

      ! ----------------------
      character(len=*) :: groupname
      ! ----------------------
      integer :: ires
      character(HDF_LEN) :: message
      ! ----------------------
      ! open the group
      call h5gopen_f (hdf_active_group%id,groupname,hdf_group_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to open the HDF group "//trim(groupname),is_fatal=.true.)

      ! allocate the group id
      allocate(hdf_active_group%child, stat=ires)
      if (ires /= 0) then
         write(message,'(a)')"[MEMORY error] Unable to allocate memory in hdf_open_group"
         call error_and_exit(message)
      endif
      hdf_active_group%child%id = hdf_group_id
      hdf_active_group%child%name = trim(groupname)
      NULLIFY(hdf_active_group%child%parent)
      NULLIFY(hdf_active_group%child%child)
      hdf_active_group%child%parent => hdf_active_group
      hdf_active_group => hdf_active_group%child
   end subroutine hdf_open_group
   !=================================   

   !=================================   
   subroutine hdf_open_or_create_group(groupname)

      implicit none

      ! ----------------------
      character(len=*) :: groupname
      ! ----------------------
      integer :: ires
      character(HDF_LEN) :: message
      ! ----------------------

      ! create the group
      if (hdf_find_group(groupname)) then
         call h5gopen_f (hdf_active_group%id,groupname,hdf_group_id,hdf_ierror)
      else
         call h5gcreate_f (hdf_active_group%id,groupname,hdf_group_id,hdf_ierror)
      endif

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to create or open the HDF group "//trim(groupname),is_fatal=.true.)

      ! allocate the group id
      allocate(hdf_active_group%child, stat=ires)
      if (ires /= 0) then
         write(message,'(a)')"[MEMORY error] Unable to allocate memory in hdf_open_or_create_group"
         call error_and_exit(message)
      endif

      hdf_active_group%child%id = hdf_group_id
      hdf_active_group%child%name = trim(groupname)
      NULLIFY(hdf_active_group%child%parent)
      NULLIFY(hdf_active_group%child%child)
      hdf_active_group%child%parent => hdf_active_group
      hdf_active_group => hdf_active_group%child
   end subroutine hdf_open_or_create_group
   !=================================   

   !=================================   
   function hdf_find_group(groupname) result(found)

      implicit none

      ! ----------------------
      character(len=*), intent(in) :: groupname
      logical :: found
      ! ----------------------
      character(len=HDF_LEN) :: resname
      integer :: nd,objtype
      ! ----------------------


      ! search the dataset in the group
      found = .false.
      do nd=1,hdf_count_elements()
         call h5gget_obj_info_idx_f (hdf_active_group%parent%id,trim(hdf_active_group%name),nd-1, &
            resname,objtype,hdf_ierror)
         if ((hdf_ierror==0).and.(trim(resname)==trim(groupname)).and.(objtype==H5G_GROUP_F)) then
            found = .true.
         end if
      end do
   end function hdf_find_group
   !=================================   

   !=================================   
   function hdf_count_elements() result(nd)

      implicit none

      ! ----------------------
      integer :: nd
      ! ----------------------
      ! ----------------------

      ! count the datasets
      call h5gn_members_f (hdf_active_group%parent%id,trim(hdf_active_group%name),nd,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to count the HDF elements")
   end function hdf_count_elements
   !=================================   

   !=================================    
   function hdf_find_element_from_nb(nb,elementname) result(found)

      implicit none

      ! ----------------------
      integer :: nb
      character(len=*) :: elementname
      logical :: found
      ! ----------------------
      character(len=HDF_LEN) :: resname
      integer :: objtype
      ! ----------------------


      ! get the element in the group
      found = .false.
      resname = ""
      call h5gget_obj_info_idx_f (hdf_active_group%parent%id,trim(hdf_active_group%name), &
         nb-1,resname,objtype,hdf_ierror)
      if (hdf_ierror==0) then
         found = .true.
         elementname = trim(resname)
      else
         write(*,*)hdf_active_group%parent%id,trim(hdf_active_group%name)," ",trim(resname)
      end if
   end function hdf_find_element_from_nb

   !=================================   

   !=================================   
   subroutine hdf_close_group()

      implicit none
      ! close the group
      call h5gclose_f (hdf_active_group%id,hdf_ierror)

      ! deallocate the group id
      hdf_active_group => hdf_active_group%parent
      deallocate(hdf_active_group%child)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to close the HDF group",is_fatal=.true.)
   end subroutine hdf_close_group
   !================================= 


   !================================================================
   !                        FILE OPERATIONS  
   !================================================================

   !=================================
   subroutine hdf_create_file(filename,compression)

      implicit none

      ! ----------------------
      character(len=*) :: filename
      integer, intent(in), optional :: compression
      ! ----------------------
      integer :: ires
      character(len=HDF_LEN) :: message
      ! ----------------------

      ! parameters
      hdf_compression = 0
      if (present(compression)) then
         hdf_compression = compression
      end if
      if ((hdf_compression>9).or.(hdf_compression<0)) then
         call error_and_exit("Issue with the compression level")
      end if

      ! create the file
      call H5Pcreate_f (H5P_FILE_ACCESS_F, hdf_parallel_access, hdf_ierror)
      call H5Pset_fapl_mpio_f (hdf_parallel_access, MPI_COMM_WORLD, MPI_INFO_NULL, hdf_ierror)
      call H5Fcreate_f (filename, H5F_ACC_TRUNC_F, hdf_file_id, hdf_ierror, access_prp=hdf_parallel_access)

      ! allocate the group id
      allocate(hdf_first_group, stat=ires)
      if (ires /= 0) then
         write(message,'(a)')"[MEMORY error] Unable to allocate memory in hdf_create_file"
         call error_and_exit(message)
      endif
      hdf_first_group%id = hdf_file_id
      NULLIFY(hdf_first_group%parent)
      NULLIFY(hdf_first_group%child)
      hdf_active_group => hdf_first_group
      hdf_active_group%parent => hdf_first_group
      hdf_active_group%name = "/"

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to create the HDF file "//trim(filename),is_fatal=.true.)
   end subroutine hdf_create_file
   !=================================

   !=================================
   subroutine hdf_open_file(filename,compression,readonly)

      implicit none

      ! ----------------------
      character(len=*) :: filename
      integer, intent(in), optional :: compression
      logical, intent(in), optional :: readonly
      ! ----------------------
      logical :: rdnly
      integer :: ires
      character(len=HDF_LEN) :: message

      ! ----------------------

      ! parameters
      hdf_compression = 0
      if (present(compression)) then
         hdf_compression = compression
      end if
      if ((hdf_compression>9).or.(hdf_compression<0)) then
         call error_and_exit("Issue with the compression level")
      end if
      rdnly = .false.
      if (present(readonly)) then
         rdnly = readonly
      end if

      ! open the file
      if (rdnly) then
         call h5fopen_f (filename,H5F_ACC_RDONLY_F,hdf_file_id,hdf_ierror)
      else
         call h5fopen_f (filename,H5F_ACC_RDWR_F,hdf_file_id,hdf_ierror)
      end if

      ! allocate the group id
      allocate(hdf_first_group, stat=ires)
      if (ires /= 0) then
         write(message,'(a)')"[MEMORY error] Unable to allocate memory in hdf_open_file"
         call error_and_exit(message)
      endif
      hdf_first_group%id = hdf_file_id
      NULLIFY(hdf_first_group%parent)
      NULLIFY(hdf_first_group%child)
      hdf_active_group => hdf_first_group
      hdf_active_group%parent => hdf_first_group
      hdf_active_group%name = "/"

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to open the HDF file "//trim(filename),is_fatal=.true.)
   end subroutine hdf_open_file
   !=================================

   !=================================
   subroutine hdf_close_file()
      implicit none
      ! close the file
      call h5fclose_f (hdf_file_id,hdf_ierror)
      ! deallocate the group id
      deallocate(hdf_first_group)
      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to close the HDF file",is_fatal=.true.)
   end subroutine hdf_close_file
   !=================================

   !================================================================
   !                        DATASET OPERATIONS  
   !================================================================

   !==================================
   subroutine hdf_create_dataset(datasetname,datatype)

      implicit none

      ! ----------------------
      character(len=*) :: datasetname
      integer, intent(in) :: datatype
      ! ----------------------

      ! get the data type
      select case(datatype)
      case (HDF_INT)
         hdf_datatype_id = H5T_NATIVE_INTEGER
      case (HDF_REAL)
         hdf_datatype_id = H5T_NATIVE_DOUBLE
      case (HDF_CHAR)
         hdf_datatype_id = H5T_NATIVE_CHARACTER
      case (HDF_CURRENT)
         ! use the current data type
      end select

      ! create the dataset
      call h5dcreate_f (hdf_active_group%id,datasetname,hdf_datatype_id, &
         hdf_dataspace_id,hdf_dataset_id,hdf_ierror,hdf_dataset_properties_id)

      call h5pcreate_f (H5P_DATASET_XFER_F, hdf_parallel_io, hdf_ierror)
      call h5pset_dxpl_mpio_f (hdf_parallel_io, H5FD_MPIO_COLLECTIVE_F, hdf_ierror)


      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to create the HDF dataset "//trim(datasetname),is_fatal=.true.)

      ! store the name
      hdf_dataset_name = datasetname
   end subroutine hdf_create_dataset
   !=================================

   !=================================
   subroutine hdf_open_dataset(datasetname,ierror)

      implicit none

      ! ----------------------
      character(len=*) :: datasetname
      integer, optional :: ierror
      ! ----------------------
      ! open the dataset
      call h5dopen_f (hdf_active_group%id,datasetname,hdf_dataset_id,hdf_ierror)

      ! check if there is an error
      if (present(ierror)) then
         ierror = hdf_ierror
      else
         call hdf_error(hdf_ierror,"Unable to open the HDF dataset "//trim(datasetname))
      end if

      ! store the name
      hdf_dataset_name = datasetname
   end subroutine hdf_open_dataset
   !=================================
   !=================================
   subroutine hdf_close_dataset()

      implicit none

      ! close the dataset
      call h5dclose_f (hdf_dataset_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to close the HDF dataset",is_fatal=.true.)
   end subroutine hdf_close_dataset
   !=================================

   !=================================
   function hdf_find_dataset(datasetname) result(found)

      implicit none

      ! ----------------------
      character(len=*), intent(in) :: datasetname
      logical :: found
      ! ----------------------
      character(len=HDF_LEN) :: resname
      integer :: nd,objtype
      ! ----------------------

      ! search the dataset in the group
      found = .false.
      do nd=1,hdf_count_elements()
         call h5gget_obj_info_idx_f (hdf_active_group%parent%id,trim(hdf_active_group%name),nd-1, &
            resname,objtype,hdf_ierror)
         if ((hdf_ierror==0).and.(trim(resname)==trim(datasetname)).and.(objtype==H5G_DATASET_F)) then
            found = .true.
         end if
      end do
   end function hdf_find_dataset
   !=================================


   !================================================================
   !                        DATASPACE OPERATIONS  
   !================================================================

   !=================================
   subroutine hdf_set_dims(ndim,dims)

      implicit none
      ! ----------------------
      integer,intent(in) :: ndim
      integer,intent(in) :: dims(ndim)
      ! ----------------------
      integer :: i
      ! ----------------------

      if (ndim > HDF_MAX_DIMS) then
         call error_and_exit("In hdf_set_dims, dimensions exceed MAX_DIMS")
      elseif (size(dims, dim=1) /= ndim) then
         call error_and_exit("In hdf_set_dims, dimensions number is non conforming to array size")
      endif

      hdf_dataset_dim(1:HDF_MAX_DIMS) = 0
      do i=1, ndim
         hdf_dataset_dim(i) = dims(i)
      enddo
   end subroutine hdf_set_dims
   !=================================

   !=================================
   subroutine hdf_set_offsets(ndim,offsets)

      implicit none
      ! ----------------------
      integer,intent(in) :: ndim
      integer,intent(in) :: offsets(ndim)
      ! ----------------------
      integer :: i
      ! ----------------------

      if (ndim > HDF_MAX_DIMS) then
         call error_and_exit("In hdf_set_offsets, dimensions exceed MAX_DIMS")
      elseif (size(offsets, dim=1) /= ndim) then
         call error_and_exit("In hdf_set_offsets, dimensions number is non conforming to array size")
      endif

      hdf_dataset_offset(1:HDF_MAX_DIMS) = 0
      do i=1, ndim
         hdf_dataset_offset(i) = offsets(i)
      enddo
   end subroutine hdf_set_offsets
   !=================================
   !=================================
   subroutine hdf_create_space(ndim,space_kind,dims)
      implicit none
      ! ----------------------
      integer,intent(in) :: ndim
      integer,intent(in) :: space_kind
      integer,intent(in) :: dims(ndim)
      ! ----------------------

      call hdf_set_dims(ndim, dims)
      select case (space_kind)
      case (1) ! dataspace
         call h5screate_simple_f (ndim,hdf_dataset_dim,hdf_dataspace_id,hdf_ierror)
         call hdf_error(hdf_ierror,"Unable to create the HDF dataspace")
      case (2) ! filespace
         call h5screate_simple_f (ndim,hdf_dataset_dim,hdf_filespace_id,hdf_ierror)
         call hdf_error(hdf_ierror,"Unable to create the HDF filespace")
      case (3) ! dataspace_slab
         call h5screate_simple_f (ndim,hdf_dataset_dim,hdf_dataspace_slab_id,hdf_ierror)
         call hdf_error(hdf_ierror,"Unable to create the HDF dataspace slab")
      case default
         call error_and_exit("Unknown HDF space")
      end select
   end subroutine hdf_create_space
   !=================================

   !=================================
   subroutine hdf_create_dataspace(ndim,dims)

      implicit none
      ! ----------------------
      integer,intent(in) :: ndim
      integer,intent(in) :: dims(ndim)
      ! ----------------------

      call hdf_create_space(ndim,1,dims)
   end subroutine hdf_create_dataspace
   !=================================

   !=================================
   subroutine hdf_create_filespace(ndim,dims)

      implicit none
      ! ----------------------
      integer,intent(in) :: ndim
      integer,intent(in) :: dims(ndim)
      ! ----------------------

      call hdf_create_space(ndim,2,dims)

   end subroutine hdf_create_filespace
   !=================================

   !=================================
   subroutine hdf_create_dataspace_slab(ndim,dims)

      implicit none
      ! ----------------------
      integer,intent(in) :: ndim
      integer,intent(in) :: dims(ndim)
      ! ----------------------

      call hdf_create_space(ndim,3,dims)

   end subroutine hdf_create_dataspace_slab
   !=================================

   !=================================
   subroutine hdf_select_hyperslab(slab_kind, ndim,dims,offsets)
      implicit none

      ! ----------------------
      integer,intent(in)           :: slab_kind
      integer,intent(in)           :: ndim
      integer,intent(in)           :: dims(ndim)
      integer,intent(in), optional :: offsets(ndim)

      ! ----------------------
      ! ----------------------
      ! get dimensions
      call hdf_set_dims(ndim, dims)
      if (present(offsets)) then
         call hdf_set_offsets(ndim, offsets)
      else
         hdf_dataset_offset(1:HDF_MAX_DIMS) = 0
      endif

      select case (slab_kind)
      case (1) ! dataspace
         call h5sselect_hyperslab_f (hdf_dataspace_id,H5S_SELECT_SET_F, &
            hdf_dataset_offset,hdf_dataset_dim,hdf_ierror)
         call hdf_error(hdf_ierror,"Unable to select the HDF dataspace slab")
      case (2) ! filespace
         call h5sselect_hyperslab_f (hdf_filespace_id,H5S_SELECT_SET_F, &
            hdf_dataset_offset,hdf_dataset_dim,hdf_ierror)
         call hdf_error(hdf_ierror,"Unable to select the HDF filespace slab")
      case (3) ! dataspace_slab
         call h5sselect_hyperslab_f (hdf_dataspace_slab_id,H5S_SELECT_SET_F, &
            hdf_dataset_offset,hdf_dataset_dim,hdf_ierror)
         call hdf_error(hdf_ierror,"Unable to select the HDF slab")
      end select

   end subroutine hdf_select_hyperslab
   !=================================

   !=================================
   subroutine hdf_select_dataspace_slab(ndim, dims, offsets)
      implicit none

      ! ----------------------
      integer,intent(in)           :: ndim
      integer,intent(in)           :: dims(ndim)
      integer,intent(in), optional :: offsets(ndim)
      ! ----------------------

      if (present(offsets)) then
         call hdf_select_hyperslab(1, ndim,dims,offsets)
      else
         call hdf_select_hyperslab(1, ndim,dims)
      endif
   end subroutine hdf_select_dataspace_slab
   !=================================

   !=================================
   subroutine hdf_select_filespace(ndim, dims, offsets)
      implicit none

      ! ----------------------
      integer,intent(in)           :: ndim
      integer,intent(in)           :: dims(ndim)
      integer,intent(in), optional :: offsets(ndim)
      ! ----------------------

      if (present(offsets)) then
         call hdf_select_hyperslab(2, ndim,dims,offsets)
      else
         call hdf_select_hyperslab(2, ndim,dims)
      endif
   end subroutine hdf_select_filespace
   !=================================

   !=================================
   subroutine hdf_select_slab(ndim, dims)

      implicit none

      ! ----------------------
      integer,intent(in)           :: ndim
      integer,intent(in)           :: dims(ndim)
      ! ----------------------

      call hdf_select_hyperslab(3, ndim,dims)
   end subroutine hdf_select_slab
   !=================================

   !=================================
   subroutine hdf_get_dataspace()

      implicit none
      ! get the dataspace
      call h5dget_space_f (hdf_dataset_id,hdf_dataspace_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to get the HDF dataspace")
   end subroutine hdf_get_dataspace
   !=================================

   !=================================
   subroutine hdf_close_dataspace()

      implicit none

      ! close the dataspace
      call h5sclose_f (hdf_dataspace_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to close the HDF dataspace")
   end subroutine hdf_close_dataspace
   !=================================

   !=================================
   subroutine hdf_close_filespace()

      implicit none

      ! close the filespace
      call h5sclose_f (hdf_filespace_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to close the HDF filespace")
   end subroutine hdf_close_filespace
   !=================================

   !=================================
   subroutine hdf_close_dataspace_slab()

      implicit none

      ! close the dataspace_slab
      call h5sclose_f (hdf_dataspace_slab_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to close the HDF dataspace_slab")
   end subroutine hdf_close_dataspace_slab
   !=================================



   !================================================================
   !                        WRITERS OPERATIONS  
   !================================================================

   !=================================  
   subroutine check_array_dims(ndims, arr_rank, arr_dims,g_dims,l_dims,name, offs)
      implicit none
      ! ----------------------
      integer, intent(in) :: ndims
      integer, intent(in) :: arr_rank
      integer, intent(in) :: arr_dims(ndims)
      integer, intent(in) :: g_dims(ndims)
      integer, intent(in) :: l_dims(ndims)
      character(len=*),intent(in) :: name
      integer, intent(in), optional :: offs(ndims)
      ! ----------------------
      integer :: i
      character(len=2) :: char_dim
      ! ----------------------
      if (arr_rank /= ndims) then
         call error_and_exit("in hdf_check_array ("//trim(name)//"): array rank mismatch with ndims")
      endif
      do i=1, ndims
         if (arr_dims(i) /= l_dims(i)) then
            write(char_dim, '(I2)') i
            call error_and_exit("in hdf_check_array ("//trim(name)//"): array rank for dimension "//&
               trim(char_dim))
         endif 
         if (present(offs)) then
            if ((arr_dims(i) + offs(i)) > g_dims(i)) then
               write(char_dim, '(I2)') i
               call error_and_exit("in hdf_check_array ("//trim(name)//"): array local size + offset "//&
                  "exeed global size for dimension"// trim(char_dim))           
            endif
         endif
      enddo
   end subroutine check_array_dims
   !=================================  

   !=================================  
   subroutine hdf_create_writing_spaces(ndims,name, glob_dims, loc_dims, datatype, offsets)

      implicit none

      ! ----------------------
      integer, intent(in) :: ndims
      character(len=*),intent(in) :: name
      integer, intent(in) :: glob_dims(ndims)
      integer, intent(in) :: loc_dims(ndims)
      integer, intent(in) :: datatype
      integer, intent(in), optional :: offsets(ndims)
      ! ----------------------

      call hdf_create_dataspace(ndims, glob_dims)
      call hdf_create_dataset(name,datatype)
      call hdf_create_filespace(ndims, glob_dims)
      if (present(offsets)) then
         call hdf_select_filespace(ndims, loc_dims, offsets)
      else
         call hdf_select_filespace(ndims, loc_dims, offsets)
      endif
      call hdf_create_dataspace_slab(ndims,loc_dims)
      call hdf_select_dataspace_slab(ndims, loc_dims)

   end subroutine hdf_create_writing_spaces
   !=================================  
   subroutine close_spaces()

      implicit none

      call hdf_close_dataspace()
      call hdf_close_filespace()
      call hdf_close_dataspace_slab()

   end subroutine close_spaces
   !================================

   !=================================  
   subroutine hdf_write_r3d_array(array, datasetname,ndims,glob_dims,loc_dims,offsets)

      implicit none

      ! ----------------------
      real(dp) :: array(:,:,:)
      character(len=*) :: datasetname
      integer, intent(in) :: ndims
      integer, intent(in) :: glob_dims(ndims)
      integer, intent(in) :: loc_dims(ndims)
      integer, intent(in), optional :: offsets(ndims)
      ! ----------------------
      integer :: i, arr_rank, arr_dims(ndims), offs(ndims)
      ! ----------------------

      ! check
      offs(1:ndims) = 0
      arr_rank = rank(array)
      do i=1, ndims
         arr_dims(i) = size(array,dim=i)
      enddo
      if (present(offsets)) then
         do i=1, ndims
            offs(i) = offsets(i)
         enddo
      endif

      call check_array_dims(ndims, arr_rank, arr_dims,glob_dims,loc_dims,datasetname, offs)
      call hdf_create_writing_spaces(ndims,datasetname,glob_dims, loc_dims, HDF_REAL, offs)
      call h5dwrite_f (hdf_dataset_id,H5T_NATIVE_DOUBLE,array,hdf_dataset_dim,hdf_ierror, &
         hdf_dataspace_slab_id,file_space_id=hdf_filespace_id,xfer_prp = hdf_parallel_io)

      call hdf_close_dataset()
      call close_spaces()
      call hdf_error(hdf_ierror,"Unable to write the real array")
   end subroutine hdf_write_r3d_array
   !=================================  

   !=================================  
   subroutine hdf_write_r2d_array(array, datasetname,ndims,glob_dims,loc_dims,offsets)

      implicit none

      ! ----------------------
      real(dp) :: array(:,:)
      character(len=*) :: datasetname
      integer, intent(in) :: ndims
      integer, intent(in) :: glob_dims(ndims)
      integer, intent(in) :: loc_dims(ndims)
      integer, intent(in), optional :: offsets(ndims)
      ! ----------------------
      integer :: i, arr_rank, arr_dims(ndims), offs(ndims)
      ! ----------------------

      ! check
      offs(1:ndims) = 0
      arr_rank = rank(array)
      do i=1, ndims
         arr_dims(i) = size(array,dim=i)
      enddo
      if (present(offsets)) then
         do i=1, ndims
            offs(i) = offsets(i)
         enddo
      endif

      call check_array_dims(ndims, arr_rank, arr_dims,glob_dims,loc_dims,datasetname, offs)
      call hdf_create_writing_spaces(ndims,datasetname,glob_dims, loc_dims, HDF_REAL, offs)
      call h5dwrite_f (hdf_dataset_id,H5T_NATIVE_DOUBLE,array,hdf_dataset_dim,hdf_ierror, &
         hdf_dataspace_slab_id,file_space_id=hdf_filespace_id,xfer_prp = hdf_parallel_io)

      call hdf_close_dataset()
      call close_spaces()
      call hdf_error(hdf_ierror,"Unable to write the real array")
   end subroutine hdf_write_r2d_array
   !=================================  

   !=================================  
   subroutine hdf_write_r1d_array(array, datasetname,ndims,glob_dims,loc_dims,offsets)

      implicit none

      ! ----------------------
      real(dp) :: array(:)
      character(len=*) :: datasetname
      integer, intent(in) :: ndims
      integer, intent(in) :: glob_dims(ndims)
      integer, intent(in) :: loc_dims(ndims)
      integer, intent(in), optional :: offsets(ndims)
      ! ----------------------
      integer :: i, arr_rank, arr_dims(ndims), offs(ndims)
      ! ----------------------

      ! check
      offs(1:ndims) = 0
      arr_rank = rank(array)
      do i=1, ndims
         arr_dims(i) = size(array,dim=i)
      enddo
      if (present(offsets)) then
         do i=1, ndims
            offs(i) = offsets(i)
         enddo
      endif

      call check_array_dims(ndims, arr_rank, arr_dims,glob_dims,loc_dims,datasetname, offs)
      call hdf_create_writing_spaces(ndims,datasetname,glob_dims, loc_dims, HDF_REAL, offs)
      call h5dwrite_f (hdf_dataset_id,H5T_NATIVE_DOUBLE,array,hdf_dataset_dim,hdf_ierror, &
         hdf_dataspace_slab_id,file_space_id=hdf_filespace_id,xfer_prp = hdf_parallel_io)

      call hdf_close_dataset()
      call close_spaces()
      call hdf_error(hdf_ierror,"Unable to write the real array")
   end subroutine hdf_write_r1d_array
   !=================================  

   !=================================  
   subroutine hdf_write_r4d_array(array, datasetname,ndims,glob_dims,loc_dims,offsets)

      implicit none

      ! ----------------------
      real(dp) :: array(:,:,:,:)
      character(len=*) :: datasetname
      integer, intent(in) :: ndims
      integer, intent(in) :: glob_dims(ndims)
      integer, intent(in) :: loc_dims(ndims)
      integer, intent(in), optional :: offsets(ndims)
      ! ----------------------
      integer :: i, arr_rank, arr_dims(ndims), offs(ndims)
      ! ----------------------

      ! check
      offs(1:ndims) = 0
      arr_rank = rank(array)
      do i=1, ndims
         arr_dims(i) = size(array,dim=i)
      enddo
      if (present(offsets)) then
         do i=1, ndims
            offs(i) = offsets(i)
         enddo
      endif

      call check_array_dims(ndims, arr_rank, arr_dims,glob_dims,loc_dims,datasetname, offs)
      call hdf_create_writing_spaces(ndims,datasetname,glob_dims, loc_dims, HDF_REAL, offs)
      call h5dwrite_f (hdf_dataset_id,H5T_NATIVE_DOUBLE,array,hdf_dataset_dim,hdf_ierror, &
         hdf_dataspace_slab_id , file_space_id = hdf_filespace_id , xfer_prp = hdf_parallel_io)

      call hdf_close_dataset()
      call close_spaces()
      call hdf_error(hdf_ierror,"Unable to write the real array")
   end subroutine hdf_write_r4d_array
   !=================================  

   !=================================  
   subroutine hdf_write_r5d_array(array, datasetname,ndims,glob_dims,loc_dims,offsets)

      implicit none

      ! ----------------------
      real(dp) :: array(:,:,:,:,:)
      character(len=*) :: datasetname
      integer, intent(in) :: ndims
      integer, intent(in) :: glob_dims(ndims)
      integer, intent(in) :: loc_dims(ndims)
      integer, intent(in), optional :: offsets(ndims)
      ! ----------------------
      integer :: i, arr_rank, arr_dims(ndims), offs(ndims)
      ! ----------------------

      ! check
      offs(1:ndims) = 0
      arr_rank = rank(array)
      do i=1, ndims
         arr_dims(i) = size(array,dim=i)
      enddo
      if (present(offsets)) then
         do i=1, ndims
            offs(i) = offsets(i)
         enddo
      endif

      call check_array_dims(ndims, arr_rank, arr_dims,glob_dims,loc_dims,datasetname, offs)
      call hdf_create_writing_spaces(ndims,datasetname,glob_dims, loc_dims, HDF_REAL, offs)
      call h5dwrite_f (hdf_dataset_id,H5T_NATIVE_DOUBLE,array,hdf_dataset_dim,hdf_ierror, &
         hdf_dataspace_slab_id,file_space_id=hdf_filespace_id,xfer_prp = hdf_parallel_io)

      call hdf_close_dataset()
      call close_spaces()
      call hdf_error(hdf_ierror,"Unable to write the real array")
   end subroutine hdf_write_r5d_array
   !=================================  

   !=================================  
   subroutine hdf_write_i1d_array(array, datasetname,ndims,glob_dims,loc_dims,offsets)

      implicit none

      ! ----------------------
      integer(ip) :: array(:)
      character(len=*) :: datasetname
      integer, intent(in) :: ndims
      integer, intent(in) :: glob_dims(ndims)
      integer, intent(in) :: loc_dims(ndims)
      integer, intent(in), optional :: offsets(ndims)
      ! ----------------------
      integer :: i, arr_rank, arr_dims(ndims), offs(ndims)
      ! ----------------------

      ! check
      offs(1:ndims) = 0
      arr_rank = rank(array)
      do i=1, ndims
         arr_dims(i) = size(array,dim=i)
      enddo
      if (present(offsets)) then
         do i=1, ndims
            offs(i) = offsets(i)
         enddo
      endif

      call check_array_dims(ndims, arr_rank, arr_dims,glob_dims,loc_dims,datasetname, offs)
      call hdf_create_writing_spaces(ndims,datasetname,glob_dims, loc_dims, HDF_INT, offs)
      call h5dwrite_f (hdf_dataset_id,H5T_NATIVE_INTEGER,array,hdf_dataset_dim,hdf_ierror, &
         hdf_dataspace_slab_id,file_space_id=hdf_filespace_id,xfer_prp = hdf_parallel_io)

      call hdf_close_dataset()
      call close_spaces()
      call hdf_error(hdf_ierror,"Unable to write the real array")
   end subroutine hdf_write_i1d_array
   !=================================  
   subroutine hdf_write_i2d_array(array, datasetname,ndims,glob_dims,loc_dims,offsets)

      implicit none

      ! ----------------------
      integer(ip) :: array(:,:)
      character(len=*) :: datasetname
      integer, intent(in) :: ndims
      integer, intent(in) :: glob_dims(ndims)
      integer, intent(in) :: loc_dims(ndims)
      integer, intent(in), optional :: offsets(ndims)
      ! ----------------------
      integer :: i, arr_rank, arr_dims(ndims), offs(ndims)
      ! ----------------------

      ! check
      offs(1:ndims) = 0
      arr_rank = rank(array)
      do i=1, ndims
         arr_dims(i) = size(array,dim=i)
      enddo
      if (present(offsets)) then
         do i=1, ndims
            offs(i) = offsets(i)
         enddo
      endif

      call check_array_dims(ndims, arr_rank, arr_dims,glob_dims,loc_dims,datasetname, offs)
      call hdf_create_writing_spaces(ndims,datasetname,glob_dims, loc_dims, HDF_INT, offs)
      call h5dwrite_f (hdf_dataset_id,H5T_NATIVE_INTEGER,array,hdf_dataset_dim,hdf_ierror, &
         hdf_dataspace_slab_id,file_space_id=hdf_filespace_id,xfer_prp = hdf_parallel_io)

      call hdf_close_dataset()
      call close_spaces()
      call hdf_error(hdf_ierror,"Unable to write the real array")
   end subroutine hdf_write_i2d_array
   !=================================  
   !=================================  
   subroutine hdf_write_i3d_array(array, datasetname,ndims,glob_dims,loc_dims,offsets)

      implicit none

      ! ----------------------
      integer(ip) :: array(:,:,:)
      character(len=*) :: datasetname
      integer, intent(in) :: ndims
      integer, intent(in) :: glob_dims(ndims)
      integer, intent(in) :: loc_dims(ndims)
      integer, intent(in), optional :: offsets(ndims)
      ! ----------------------
      integer :: i, arr_rank, arr_dims(ndims), offs(ndims)
      ! ----------------------

      ! check
      offs(1:ndims) = 0
      arr_rank = rank(array)
      do i=1, ndims
         arr_dims(i) = size(array,dim=i)
      enddo
      if (present(offsets)) then
         do i=1, ndims
            offs(i) = offsets(i)
         enddo
      endif

      call check_array_dims(ndims, arr_rank, arr_dims,glob_dims,loc_dims,datasetname, offs)
      call hdf_create_writing_spaces(ndims,datasetname,glob_dims, loc_dims, HDF_INT, offs)
      call h5dwrite_f (hdf_dataset_id,H5T_NATIVE_INTEGER,array,hdf_dataset_dim,hdf_ierror, &
         hdf_dataspace_slab_id,file_space_id=hdf_filespace_id,xfer_prp = hdf_parallel_io)

      call hdf_close_dataset()
      call close_spaces()
      call hdf_error(hdf_ierror,"Unable to write the real array")
   end subroutine hdf_write_i3d_array
   !=================================  

   ! !=================================
   subroutine hdf_read_r3d_array(datasetname,r3d_array,dims,offsets)

      implicit none

      ! ----------------------
      character(len=*)    :: datasetname
      real(dp)            :: r3d_array(:,:,:)
      integer, intent(in) :: dims(:)
      integer, intent(in) :: offsets(3)
      ! ----------------------
      ! ----------------------

      ! check if the dataset exists or not
      if (hdf_find_dataset(datasetname)) then
         ! open the dataset
         call hdf_open_dataset(datasetname,hdf_ierror)
         call hdf_get_dataspace()
         call hdf_select_dataspace_slab(3, dims, offsets)
         call hdf_create_dataspace_slab(3,dims)
         call hdf_select_slab(3, dims)
         call h5dread_f (hdf_dataset_id,H5T_NATIVE_DOUBLE,r3d_array,hdf_dataset_dim,hdf_ierror, &
            hdf_dataspace_slab_id, hdf_dataspace_id)

         call hdf_error(hdf_ierror,"Unable to read r3d_array")
         call hdf_close_dataspace()
         call hdf_close_dataset()
      else
         call hdf_error(hdf_ierror,"Unable to find dataset "//trim(datasetname))
      end if

   end subroutine hdf_read_r3d_array
   ! !=================================

   !================================================================
   !                        ATTRIBUITE OPERATIONS
   !================================================================
   !====================================
   subroutine hdf_create_attribute(attributename,datatype)

      implicit none

      ! ----------------------
      character(len=*) :: attributename
      integer, intent(in), optional :: datatype
      ! ----------------------
      ! ----------------------

      ! get the data type
      if (present(datatype)) then
         select case(datatype)
         case (HDF_INT)
            hdf_attribute_type_id = H5T_NATIVE_INTEGER
         case (HDF_REAL)
            hdf_attribute_type_id = H5T_NATIVE_DOUBLE
         case (HDF_CHAR)
            hdf_attribute_type_id = H5T_NATIVE_CHARACTER
         case (HDF_CURRENT)
            ! use the current data type
         end select
      end if

      ! create the attribute
      call h5acreate_f (hdf_dataset_id,attributename,hdf_attribute_type_id, &
         hdf_dataspace_id,hdf_attribute_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to create the HDF attribute " //&
         trim(attributename),is_fatal=.true.)

      ! store the name
      hdf_attribute_name = attributename
   end subroutine hdf_create_attribute

   !====================================

   !====================================
   function hdf_count_attributes() result(nd)

      implicit none

      ! ----------------------
      integer :: nd
      ! ----------------------
      ! ----------------------

      ! count the attributes
      call h5aget_num_attrs_f (hdf_dataset_id,nd,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to count the HDF attributes")
   end function hdf_count_attributes

   !====================================

   !====================================
   function hdf_attribute_exists(attributename) result(attribute_exists)

      implicit none

      ! ----------------------
      character(len=*) :: attributename
      logical :: attribute_exists
      ! ----------------------
      ! ----------------------

      ! check if the attribute exists
      call h5aexists_f (hdf_dataset_id,attributename,attribute_exists,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to check the HDF attribute "//trim(attributename))
   end function hdf_attribute_exists

   !====================================

   !====================================
   subroutine hdf_open_attribute(attributename,ierror)

      implicit none

      ! ----------------------
      character(len=*) :: attributename
      integer, optional :: ierror
      ! ----------------------
      ! ----------------------

      ! open the attribute
      call h5aopen_f (hdf_dataset_id,attributename,hdf_attribute_id,hdf_ierror)

      ! check if there is an error
      if (present(ierror)) then
         ierror = hdf_ierror
      else
         call hdf_error(hdf_ierror,"Unable to open the HDF attribute "//trim(attributename))
      end if

      ! store the name
      hdf_attribute_name = attributename
   end subroutine hdf_open_attribute

   !====================================

   !====================================
   function hdf_find_attribute_from_nb(nb,attrname) result(found)

      implicit none

      ! ----------------------
      integer :: nb
      character(len=*) :: attrname
      logical :: found
      ! ----------------------
      character(len=HDF_LEN) :: resname
      integer(HSIZE_T) :: n
      ! ----------------------

      ! get the element in the group
      found = .false.
      resname = ""
      n = int(nb-1,kind(n))
      call h5aget_name_by_idx_f (hdf_active_group%id,trim(hdf_dataset_name),H5_INDEX_NAME_F, &
         H5_ITER_INC_F,n,resname,hdf_ierror)

      if (hdf_ierror==0) then
         found = .true.
         attrname = trim(resname)
      else
         write(*,*)hdf_active_group%id,trim(hdf_dataset_name)," ",trim(resname)
      end if
   end function hdf_find_attribute_from_nb

   !====================================

   !====================================
   subroutine hdf_get_attribute_type(attrtype)

      implicit none

      ! ----------------------
      integer :: attrtype
      ! ----------------------
      integer :: classtype
      character(HDF_LEN) :: message
      ! ----------------------

      ! get the datatype id
      call h5aget_type_f (hdf_attribute_id,hdf_attribute_type_id,hdf_ierror)

      ! get the native datatype
      call h5tget_class_f (hdf_attribute_type_id,classtype,hdf_ierror)

      ! translate the type
      if (classtype==H5T_INTEGER_F) then
         attrtype = HDF_INT
      else if (classtype==H5T_FLOAT_F) then
         attrtype = HDF_REAL
      else if (classtype==H5T_STRING_F) then
         attrtype = HDF_CHAR
      else
         write(message,*)"Unknown attribute type ",classtype
         call error_and_exit(message)
      end if

      ! close the data type
      call h5tclose_f (hdf_attribute_type_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to get the type of the attribute")
   end subroutine hdf_get_attribute_type

   !====================================

   !====================================
   subroutine hdf_get_attribute_size(dim1)

      implicit none

      ! ----------------------
      integer :: dim1
      ! ----------------------
      integer :: ndim
      ! ----------------------

      ! get the dataspace
      call h5dget_space_f (hdf_attribute_id,hdf_dataspace_id,hdf_ierror)

      ! get the number of dimensions
      call h5sget_simple_extent_ndims_f (hdf_dataspace_id,ndim,hdf_ierror)

      ! get the size of the dimensions
      call h5sget_simple_extent_dims_f (hdf_dataspace_id,hdf_dataset_dim,hdf_dataset_dim_max,hdf_ierror)

      ! compute the length
      dim1 = int(sum(hdf_dataset_dim(1:ndim)))

      ! close dataspace
      call hdf_close_dataspace()

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to get the size of the attribute")
   end subroutine hdf_get_attribute_size

   !====================================

   !====================================
   subroutine hdf_get_attribute_size_by_name(attributename,dim1)

      implicit none

      ! ----------------------
      character(len=*) :: attributename
      integer :: dim1
      ! ----------------------
      ! ----------------------

      ! open the attribute
      call hdf_open_attribute(attributename,hdf_ierror)

      ! get the number of dimensions
      call hdf_get_attribute_size(dim1)

      ! close attribute
      call hdf_close_attribute()

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to get the size of the attribute from its name")
   end subroutine hdf_get_attribute_size_by_name

   !====================================

   !====================================
   subroutine hdf_delete_attribute(attributename)

      implicit none

      ! ----------------------
      character(len=*) :: attributename
      ! ----------------------
      ! ----------------------

      ! delete the attribute
      call h5adelete_f (hdf_dataset_id,attributename,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to delete the HDF attribute")
   end subroutine hdf_delete_attribute

   !====================================

   !====================================
   subroutine hdf_close_attribute()

      implicit none

      ! ----------------------
      ! ----------------------
      ! ----------------------

      ! close the attribute
      call h5aclose_f (hdf_attribute_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to close the HDF attribute",is_fatal=.true.)
   end subroutine hdf_close_attribute

   subroutine hdf_read_i0_attribute(attributename,i0)

      implicit none

      ! ----------------------
      character(len=*) :: attributename
      integer :: i0
      ! ----------------------
      ! ----------------------
      ! open the attribute
      call hdf_open_attribute(attributename,hdf_ierror)

      ! get the i0
      call hdf_read_current_i0_attribute(i0)

      ! close attribute
      call hdf_close_attribute()

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to read the i0 attribute")
   end subroutine hdf_read_i0_attribute

   !=================================================

   !=================================================
   subroutine hdf_write_i0_attribute(attributename,i0)

      implicit none

      ! ----------------------
      character(len=*) :: attributename
      integer:: i0
      ! ----------------------
      ! ----------------------

      ! create the dataspace for the attribute
      call hdf_create_dataspace(1,(/1/))

      ! create the attribute
      if (hdf_attribute_exists(attributename)) then
         call hdf_delete_attribute(attributename)
      end if
      call hdf_create_attribute(attributename,datatype=HDF_INT)

      ! write the attribute
      hdf_dataset_dim(1:HDF_MAX_DIMS) = 0
      hdf_dataset_dim(1) = 1
      call h5awrite_f (hdf_attribute_id,H5T_NATIVE_INTEGER,i0,hdf_dataset_dim,hdf_ierror)

      ! close attribute
      call hdf_close_attribute()

      ! close dataspace
      call hdf_close_dataspace()

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to write the i0 attribute")
   end subroutine hdf_write_i0_attribute
   !=================================================

   !=================================================
   subroutine hdf_read_r0_attribute(attributename,r0)

      implicit none

      ! ----------------------
      character(len=*) :: attributename
      real(dp) :: r0
      ! ----------------------
      ! ----------------------
      ! open the attribute
      call hdf_open_attribute(attributename,hdf_ierror)

      ! get the r0
      call hdf_read_current_r0_attribute(r0)

      ! close attribute
      call hdf_close_attribute()

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to read the r0 attribute")
   end subroutine hdf_read_r0_attribute
   !=================================================

   !=================================================
   subroutine hdf_write_r0_attribute(attributename,r0)

      implicit none

      ! ----------------------
      character(len=*) :: attributename
      real(dp) :: r0
      ! ----------------------
      ! ----------------------
      ! create the dataspace for the attribute
      call hdf_create_dataspace(1,(/1/))

      ! create the attribute
      if (hdf_attribute_exists(attributename)) then
         call hdf_delete_attribute(attributename)
      end if
      call hdf_create_attribute(attributename,datatype=HDF_REAL)

      ! write the attribute
      hdf_dataset_dim(1:HDF_MAX_DIMS) = 0
      hdf_dataset_dim(1) = 1
      call h5awrite_f (hdf_attribute_id,H5T_NATIVE_DOUBLE,r0,hdf_dataset_dim,hdf_ierror)

      ! close attribute
      call hdf_close_attribute()

      ! close dataspace
      call hdf_close_dataspace()

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to write the r0 attribute")
   end subroutine hdf_write_r0_attribute
   !=================================================

   !=================================================
   subroutine hdf_read_current_r0_attribute(r0)

      implicit none

      ! ----------------------
      real(dp) :: r0
      ! ----------------------
      ! ----------------------
      ! init
      r0 = 0.0_dp
      ! get the datatype
      call h5aget_type_f (hdf_attribute_id,hdf_attribute_type_id,hdf_ierror)

      ! read the string
      hdf_dataset_dim(1:HDF_MAX_DIMS) = 0
      hdf_dataset_dim(1) = 1
      call h5aread_f (hdf_attribute_id,H5T_NATIVE_DOUBLE,r0,hdf_dataset_dim,hdf_ierror)

      ! close the data type
      call h5tclose_f (hdf_attribute_type_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to read the r0 attribute")

   end subroutine hdf_read_current_r0_attribute
   !=================================================

   !=================================================
   subroutine hdf_read_current_i0_attribute(i0)

      implicit none

      integer :: i0
      ! ----------------------
      ! init
      i0 = 0
      ! get the datatype
      call h5aget_type_f (hdf_attribute_id,hdf_attribute_type_id,hdf_ierror)

      ! read the string
      hdf_dataset_dim(1:HDF_MAX_DIMS) = 0
      hdf_dataset_dim(1) = 1
      call h5aread_f (hdf_attribute_id,H5T_NATIVE_INTEGER,i0,hdf_dataset_dim,hdf_ierror)

      ! close the data type
      call h5tclose_f (hdf_attribute_type_id,hdf_ierror)

      ! check if there is an error
      call hdf_error(hdf_ierror,"Unable to read the i0 attribute")
   end subroutine hdf_read_current_i0_attribute
   !=================================================

   subroutine hdf_read_r1d_array(datasetname,r1d_array,dims)

      implicit none

      ! ----------------------
      character(len=*)    :: datasetname
      real(dp)            :: r1d_array
      integer, intent(in) :: dims(:)
      ! ----------------------
      ! ----------------------

      ! check if the dataset exists or not
      if (hdf_find_dataset(datasetname)) then
         ! open the dataset
         call hdf_open_dataset(datasetname,hdf_ierror)
         call hdf_get_dataspace()
         call hdf_select_dataspace_slab(1,dims)
         call hdf_create_dataspace_slab(1,dims)
         call hdf_select_slab(1, dims)
         call h5dread_f (hdf_dataset_id,H5T_NATIVE_DOUBLE,r1d_array,hdf_dataset_dim,hdf_ierror, &
            hdf_dataspace_slab_id, hdf_dataspace_id)

         call hdf_error(hdf_ierror,"Unable to read r1d_array")
         call hdf_close_dataspace()
         call hdf_close_dataset()
      else
         call hdf_error(hdf_ierror,"Unable to find dataset "//trim(datasetname))
      end if

   end subroutine hdf_read_r1d_array
   !=================================================

   !=================================================
   subroutine hdf_read_i1d_array(datasetname,i1d_array,dims)

      implicit none

      ! ----------------------
      character(len=*)    :: datasetname
      integer(ip)         :: i1d_array
      integer, intent(in) :: dims(:)
      ! ----------------------
      ! ----------------------

      ! check if the dataset exists or not
      if (hdf_find_dataset(datasetname)) then
         ! open the dataset
         call hdf_open_dataset(datasetname,hdf_ierror)
         call hdf_get_dataspace()
         call hdf_select_dataspace_slab(1,dims)
         call hdf_create_dataspace_slab(1,dims)
         call hdf_select_slab(1, dims)
         call h5dread_f (hdf_dataset_id,H5T_NATIVE_INTEGER,i1d_array,hdf_dataset_dim,hdf_ierror, &
            hdf_dataspace_slab_id, hdf_dataspace_id)

         call hdf_error(hdf_ierror,"Unable to read i1d_array")
         call hdf_close_dataspace()
         call hdf_close_dataset()
      else
         call hdf_error(hdf_ierror,"Unable to find dataset "//trim(datasetname))
      end if

   end subroutine hdf_read_i1d_array
   !=================================================

   !------------------------------------------------------------------------------------
end module iohdf5