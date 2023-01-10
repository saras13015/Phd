!------------------------------------------------------------------------------
! MODULE: paraview
!------------------------------------------------------------------------------
!> \brief Create ParaView plots
!!
!! This module contains subroutines to write ParaView plots.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module paraview


  use parameters
  use parallel
  use input
  use type_thd
  use adim
  use variables
  use inst_plot
  use stat_plot


  implicit none


contains


!> \brief Create a VTR file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine vtr_paraview ( ifile , inp , thd , adi , sim , grid , inst , reyavg , favavg , stat )


    integer (ip)     , intent (in)                                   :: ifile  !< file number
    type (inp_type)  , intent (in)                                   :: inp    !< input derived type
    type (thd_type)  , intent (in)                                   :: thd    !< thermodynamic derived type
    type (adi_type)  , intent (in)                                   :: adi    !< non-dimensional derived type
    type (sim_type)  , intent (in)                                   :: sim    !< simulation derived type
    type (inp_grid)  , intent (in)                                   :: grid   !< grid derived type
    type (inst_type) , intent (inout)                                :: inst   !< instantaneous derived type
    type (reyavg_type) , intent (inout)                              :: reyavg !< Reynolds average derived type
    type (favavg_type) , intent (inout)                              :: favavg !< Favre average derived type
    type (stat_type) , intent (inout)                                :: stat   !< statistical derived type


    integer (kind=8)                    :: offset ! for huge 3D plots
    integer (ip) , parameter            :: nbitsimple = 4
    character (len_default) , parameter :: format_int = ' ( I50 ) ' , &
                                           format_100 = ' ( 100A ) '
    integer (ip)                        :: ix , fx , nx , &
                                           iy , fy , ny , &
                                           iz , fz , nz
    integer (ip)                        :: ok , iv
    character ( len = 150 )             :: buf
    character (len_default)             :: buf1 , buf2 , buf3 , buf4 , buf5 , buf6
    character (len_default)             :: rank_ascii , file_ascii


    write ( rank_ascii , format_restart ) rank
    write ( file_ascii , format_restart ) ifile


    call domain_paraview ( ix , fx , iy , fy , iz , fz )
    nx = fx - ix + 1
    ny = fy - iy + 1
    nz = fz - iz + 1


    ! writing the VTR file
    open ( unit = unit_plot ,                              &
           file = trim (dir_plot)  // trim (file_plot) //  &
                  trim (file_ascii) // '_' //              &
                  trim (rank_ascii) // '.vtr' ,            &
           form    = 'unformatted' ,                       &
           access  = 'stream'      ,                       &
           action  = 'write'       ,                       &
           status  = 'unknown'     ,                       &
           iostat  = ok            ,                       &
           convert = 'big_endian'  )
    if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_plot) // trim (file_plot) // trim (file_ascii))

    write ( buf , format_100 ) '<?xml version="1.0"?>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '<VTKFile type="RectilinearGrid" ' ,    &
                               'version="0.1" byte_order="BigEndian">'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf1 , format_int ) ix
    write ( buf2 , format_int ) fx
    write ( buf3 , format_int ) iy
    write ( buf4 , format_int ) fy
    write ( buf5 , format_int ) iz
    write ( buf6 , format_int ) fz

    write ( buf , format_100 ) ' <RectilinearGrid WholeExtent="'                                  &
                               // trim ( adjustl(buf1) ) // ' ' // trim ( adjustl(buf2) ) // ' '  &
                               // trim ( adjustl(buf3) ) // ' ' // trim ( adjustl(buf4) ) // ' '  &
                               // trim ( adjustl(buf5) ) // ' ' // trim ( adjustl(buf6) ) // '">'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf1 , format_int ) ix
    write ( buf2 , format_int ) fx
    write ( buf3 , format_int ) iy
    write ( buf4 , format_int ) fy
    write ( buf5 , format_int ) iz
    write ( buf6 , format_int ) fz

    write ( buf , format_100 ) '  <Piece Extent="'                                                &
                               // trim ( adjustl(buf1) ) // ' ' // trim ( adjustl(buf2) ) // ' '  &
                               // trim ( adjustl(buf3) ) // ' ' // trim ( adjustl(buf4) ) // ' '  &
                               // trim ( adjustl(buf5) ) // ' ' // trim ( adjustl(buf6) ) // '">'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '   <Coordinates>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)


    offset = 0
    write ( buf1 , format_int ) offset


    write ( buf , format_100 ) '    <DataArray type="Float32" '               , &
                               'Name="X_COORDINATES" NumberOfComponents="1" ' , &
                               'format="appended" offset="'                     &
                               // trim ( adjustl(buf1) ) // '">'
    write ( unit_plot ) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '    </DataArray>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)


    offset = offset + (nx+1) * nbitsimple
    write ( buf1 , format_int   ) offset


    write ( buf  , format_100 ) '    <DataArray type="Float32" '               , &
                                'Name="Y_COORDINATES" NumberOfComponents="1" ' , &
                                'format="appended" offset="'                     &
                                // trim ( adjustl(buf1) ) // '">'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '    </DataArray>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)


    offset = offset + (ny+1) * nbitsimple
    write ( buf1 , format_int ) offset


    write ( buf , format_100 ) '    <DataArray type="Float32" '               , &
                               'Name="Z_COORDINATES" NumberOfComponents="1" ' , &
                               'format="appended" offset="'                     &
                               // trim ( adjustl(buf1) ) // '">'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '    </DataArray>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '   </Coordinates>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)


    offset = offset + (nz+1) * nbitsimple
    write( buf1 , format_int ) offset


    write ( buf , format_100 ) '   <PointData>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    do iv = 1 , inp % nvars

      write ( buf , format_100 ) '         <DataArray type="Float32" ' ,                 &
                                 'Name="' // trim ( adjustl ( inp % var_name (iv) ) ) // &
                                 '" NumberOfComponents="1" format="appended" ' ,         &
                                 'offset="' // trim ( adjustl(buf1) ) // '">'
      write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

      write ( buf , format_100 ) '         </DataArray>'
      write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)


      offset = offset + (nx*ny*nz+1) * nbitsimple
      write ( buf1 , format_int ) offset

    end do

    write ( buf , format_100 ) '   </PointData>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '   </Piece>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '  </RectilinearGrid>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '  <AppendedData encoding="raw">'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write (unit_plot) '_' ! necessary

    ! writing the mesh
    if ( inp % nondim_grid ) then
       write (unit_plot) nx*nbitsimple , real ( grid % x (ix:fx) * adi % L_ref / inp % length_ref )
       write (unit_plot) ny*nbitsimple , real ( grid % y (iy:fy) * adi % L_ref / inp % length_ref )
       write (unit_plot) nz*nbitsimple , real ( grid % z (iz:fz) * adi % L_ref / inp % length_ref )
    else
       write (unit_plot) nx*nbitsimple , real ( grid % x (ix:fx) * adi % L_ref , sp )
       write (unit_plot) ny*nbitsimple , real ( grid % y (iy:fy) * adi % L_ref , sp )
       write (unit_plot) nz*nbitsimple , real ( grid % z (iz:fz) * adi % L_ref , sp )
    end if

    ! plotting variables
    if ( .not. inp % temp_avg .and. .not. inp % spat_avg ) then ! instantaneous plot
       do iv = 1 , inp % nvars
          call plot_inst_var ( ifile , iv , inp , thd , adi , sim , grid , inst )
       end do
    else ! statistics plot
       do iv = 1 , inp % nvars
          call plot_stat_var ( iv , inp , thd , adi , grid , reyavg , favavg , stat )
       end do
    end if

    write ( buf , format_100 ) char(10) // '  </AppendedData>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    write ( buf , format_100 ) '</VTKFile>'
    write (unit_plot) buf ( 1:len_trim(buf) ) // char(10)

    close (unit_plot)


  end subroutine vtr_paraview


!> \brief Create a PVTR file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine pvtr_paraview ( ifile , inp )


    integer (ip)     , intent (in)                                   :: ifile !< file number
    type (inp_type)  , intent (in)                                   :: inp   !< input derived type


    integer (ip) , parameter                   :: nbitsimple = 4
    character (len_default) , parameter        :: format_int = ' ( I5 ) ' , &
                                                  format_100 = ' ( 100A ) '
    integer (ip)                               :: ix , fx , &
                                                  iy , fy , &
                                                  iz , fz
    integer (ip)                               :: ok , iv
    character ( len = 5 )                      :: buf1 , buf2 , buf3 , buf4 , buf5 , buf6
    character (len_default)                    :: rank_ascii , file_ascii
    integer (ip) , dimension (:) , allocatable :: seq_index , par_index


    call domain_paraview ( ix , fx , iy , fy , iz , fz )

    allocate ( seq_index (6)       , &
               par_index (6*nproc) , &
               stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate pvtr_paraview')

    seq_index (1) = ix
    seq_index (2) = fx
    seq_index (3) = iy
    seq_index (4) = fy
    seq_index (5) = iz
    seq_index (6) = fz

    call mpi_gather ( seq_index , 6 , MPI_INTEGER , &
                      par_index , 6 , MPI_INTEGER , &
                      0 , MPI_COMM_WORLD , mpicode )

    if ( rank == rank_default ) then


       write ( file_ascii , format_restart ) ifile


       ! writing the PVTR file
       open ( unit = unit_plot ,                            &
              file = trim (dir_plot) // trim (file_plot) // &
                     trim (file_ascii) // '.pvtr' ,         &
                     form   = 'formatted' ,                 &
                     status = 'unknown'   ,                 &
                     iostat = ok          ,                 &
                     action = 'write' )
       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_plot) // trim (file_plot) // file_ascii)

       write ( unit_plot , format_100 ) '<?xml version="1.0"?>'
       write ( unit_plot , format_100 ) '<VTKFile type="PRectilinearGrid" version="0.1" byte_order="BigEndian">'

       write ( buf1 , format_int ) 1
       write ( buf2 , format_int ) ntx
       write ( buf3 , format_int ) 1
       write ( buf4 , format_int ) nty
       write ( buf5 , format_int ) 1
       write ( buf6 , format_int ) ntz

       write ( unit_plot , format_100 ) ' <PRectilinearGrid WholeExtent="'                                &
                                        // trim ( adjustl(buf1) ) // ' ' // trim ( adjustl(buf2) ) // ' ' &
                                        // trim ( adjustl(buf3) ) // ' ' // trim ( adjustl(buf4) ) // ' ' &
                                        // trim ( adjustl(buf5) ) // ' ' // trim ( adjustl(buf6) ) // '"' &
                                        // ' GhostLevel="1">'

       write ( unit_plot , format_100 ) '   <PCoordinates>'

       write ( unit_plot , format_100 ) '    <PDataArray type="Float32" '               , &
                                        'Name="X_COORDINATES" NumberOfComponents="1"/>'

       write ( unit_plot , format_100 ) '    <PDataArray type="Float32" '               , &
                                        'Name="Y_COORDINATES" NumberOfComponents="1"/>'

       write ( unit_plot , format_100 ) '    <PDataArray type="Float32" '               , &
                                        'Name="Z_COORDINATES" NumberOfComponents="1"/>'

       write ( unit_plot , format_100 ) '   </PCoordinates>'
       write ( unit_plot , format_100 ) '   <PPointData>'
       do iv = 1 , inp % nvars
          write ( unit_plot , format_100 ) '         <PDataArray type="Float32" ' ,                      &
                                           'Name="' // trim ( adjustl ( inp % var_name (iv) ) ) // '"/>'
       end do
       write ( unit_plot , format_100 ) '   </PPointData>'

       do iv = 0 , nproc-1
          write ( rank_ascii , format_restart ) iv
          write ( file_ascii , format_restart ) ifile

          write ( buf1 , format_int ) par_index (iv*6+1)
          write ( buf2 , format_int ) par_index (iv*6+2)
          write ( buf3 , format_int ) par_index (iv*6+3)
          write ( buf4 , format_int ) par_index (iv*6+4)
          write ( buf5 , format_int ) par_index (iv*6+5)
          write ( buf6 , format_int ) par_index (iv*6+6)

          write ( unit_plot , * ) '   <Piece Extent="' , trim(buf1) , trim(buf2) ,         &
                                  trim(buf3) , trim(buf4) , trim(buf5) , trim(buf6) ,      &
                                  '" Source="' // trim (file_plot) // trim (file_ascii) // &
                                  '_' // trim (rank_ascii) // '.vtr"/>'
       end do

       write ( unit_plot , format_100 ) ' </PRectilinearGrid>'
       write ( unit_plot , format_100 ) '</VTKFile>'


       close (unit_plot)


       write (*,*) 'file ' , trim (dir_parent) // trim (dir_plot) // trim (file_plot) // trim (file_ascii) , ' written'


    end if


    deallocate ( seq_index , par_index )


  end subroutine pvtr_paraview


!> \brief Create a VTR file.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine pvd_paraview ( inp , adi , sim , t )


    type (inp_type) , intent (in)                               :: inp !< file number
    type (adi_type) , intent (in)                               :: adi !< non-dimensional derived type
    type (sim_type) , intent (in)                               :: sim !< simulation derived type
    real (dp) , allocatable , dimension (:) , intent (in)       :: t   !< time


    character (len_default) , parameter :: format_int = ' ( I8 ) ' , &
                                           format_100 = ' ( 100A ) '

    integer (ip)                        :: ok , ifile

    integer (ip)                        :: start_file , end_file , skip_file

    character (len_default)             :: file_ascii

    character ( len = 150 )             :: time , buf


    if ( rank == rank_default ) then


       skip_file  = inp % skip_file
       if  ( inp % read_stat ) then
          start_file = inp % end_file
          end_file   = inp % end_file
       else
          if ( inp % temp_avg ) then
             start_file = inp % end_file
             end_file   = inp % end_file
          else
             start_file = inp % start_file
             end_file   = inp % end_file
          end if
       end if


       write ( file_ascii , format_restart ) end_file

       open ( unit = unit_plot ,                            &
              file = trim (dir_plot) // trim (file_plot) // &
                     trim (file_ascii) // '.pvd' ,          &
              form   = 'formatted' ,                        &
              status = 'unknown'   ,                        &
              iostat = ok          ,                        &
              action = 'write' )
       if ( ok /= 0 ) call abort_mpi ('error opening ' // trim (dir_plot) // trim (file_plot) // file_ascii)


       write ( unit_plot , format_100 ) '<?xml version="1.0"?>'
       write ( unit_plot , format_100 ) '<VTKFile type="Collection" version="0.1">'
       write ( unit_plot , format_100 ) '  <Collection>'


       do ifile = start_file , end_file , skip_file

          write ( time , format_restart ) ifile

          write ( buf , format_100 ) trim (file_plot) // time ( 1:len_trim(time) )

          write ( time , * ) t (ifile) * adi % time_ref

          write ( unit_plot , format_100 ) '    <DataSet timestep="' // &
                                           trim ( adjustl(time) )    // &
                                           '" file="'                // &
                                           trim (adjustl(buf) )      // &
                                           '.pvtr"/>'

       end do


       write ( unit_plot , format_100 ) '  </Collection>'
       write ( unit_plot , format_100 ) '</VTKFile>'

       close (unit_plot)


    end if


  end subroutine pvd_paraview


end module paraview
