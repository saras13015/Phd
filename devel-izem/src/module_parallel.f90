!------------------------------------------------------------------------------
! MODULE: parallel
!------------------------------------------------------------------------------
!> \brief Parallel handling.
!!
!! This module contains all the operations related to parallel
!! communications through the Message Passing Interface (MPI) library.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module parallel

  use parameters

  use MPI

  implicit none


  real (dp)                                     :: tcomm !< time spent in communications

  integer (ip)                                  :: rank !< rank

  integer (ip)                                  :: nproc !< total number of processors

  integer (ip)                                  :: comm3d                         !< 3D communicator
  integer (ip)                                  :: comm2dxy , comm2dxz , comm2dyz !< 2D communicators
  integer (ip)                                  :: comm1dx , comm1dy , comm1dz    !< 1D communicators

  integer (ip) , parameter                      :: ndimcom = 3 !< dimension of comunicator: force to be always 3D

  integer (ip) , dimension (ndimcom)            :: dims !< parallel dimensions

  integer (ip) , dimension (ndimcom)            :: subdom_size !< Subdomain dimensions

  integer (ip) , dimension (ndimcom)            :: subdom_offset !< Subdomain offset (index in global grid)

  logical , dimension (ndimcom)                 :: period !< periodic BCs handling through MPI

  integer (ip) , dimension (ndimcom)            :: coords , dcoords !< coordinates MPI variable

  integer (ip) , parameter                      :: dimarray = 4 , dimone = 3         !< parameters for the conserved and the single array

  integer (ip) , parameter                      :: nneigh = 6 , nneighd = 4          !< number of neighbours: cartesian (6), diagonal (4)
  integer (ip) , parameter                      :: W=1 , E=2 , S=3 , N=4 , B=5 , F=6 !< direction parameters
  integer (ip) , dimension (nneigh)             :: neigh                             !< cartesian neighbours array
  integer (ip) , dimension (nneighd)            :: neighdxy , neighdxz , neighdyz    !< diagonal neighbours array

  character (len_default) , dimension (nneigh)  :: bc !< boundary conditions at neighbours

  !> MPI communication types for a conserved array
  integer (ip)                                  :: type_send_W , type_recv_W , &
  type_send_E , type_recv_E , &
  type_send_S , type_recv_S , &
  type_send_N , type_recv_N , &
  type_send_B , type_recv_B , &
  type_send_F , type_recv_F

  !> MPI communication types for a single array
  integer (ip)                                  :: type_send_one_W , type_recv_one_W , &
  type_send_one_E , type_recv_one_E , &
  type_send_one_S , type_recv_one_S , &
  type_send_one_N , type_recv_one_N , &
  type_send_one_B , type_recv_one_B , &
  type_send_one_F , type_recv_one_F

  !> MPI communication types for the viscous derivatives array
  integer (ip)                                  :: type_send_deriv_W , type_recv_deriv_W , &
  type_send_deriv_E , type_recv_deriv_E , &
  type_send_deriv_S , type_recv_deriv_S , &
  type_send_deriv_N , type_recv_deriv_N , &
  type_send_deriv_B , type_recv_deriv_B , &
  type_send_deriv_F , type_recv_deriv_F

  !> MPI communication types for the LES viscous additional derivatives array
  integer (ip)                                  :: type_send_derivSGS_W , type_recv_derivSGS_W , &
  type_send_derivSGS_E , type_recv_derivSGS_E , &
  type_send_derivSGS_S , type_recv_derivSGS_S , &
  type_send_derivSGS_N , type_recv_derivSGS_N , &
  type_send_derivSGS_B , type_recv_derivSGS_B , &
  type_send_derivSGS_F , type_recv_derivSGS_F

  !> MPI communication types for a single array through diagonals
  integer (ip)                                  :: type_send_edge_xyp1 , type_recv_edge_xyp1 , &
  type_send_edge_xyn1 , type_recv_edge_xyn1 , &
  type_send_edge_xzp1 , type_recv_edge_xzp1 , &
  type_send_edge_xzn1 , type_recv_edge_xzn1 , &
  type_send_edge_yzp1 , type_recv_edge_yzp1 , &
  type_send_edge_yzn1 , type_recv_edge_yzn1

  !> MPI communication types for a single array through diagonals
  integer (ip)                                  :: type_send_edge_xyp2 , type_recv_edge_xyp2 , &
  type_send_edge_xyn2 , type_recv_edge_xyn2 , &
  type_send_edge_xzp2 , type_recv_edge_xzp2 , &
  type_send_edge_xzn2 , type_recv_edge_xzn2 , &
  type_send_edge_yzp2 , type_recv_edge_yzp2 , &
  type_send_edge_yzn2 , type_recv_edge_yzn2


  !> MPI communication types for velocity derivatives array
  integer (ip) :: type_send_derivVel_W , type_recv_derivVel_W , &
  &               type_send_derivVel_E , type_recv_derivVel_E , &
  &               type_send_derivVel_S , type_recv_derivVel_S , &
  &               type_send_derivVel_N , type_recv_derivVel_N , &
  &               type_send_derivVel_B , type_recv_derivVel_B , &
  &               type_send_derivVel_F , type_recv_derivVel_F


  !> More communication types
  integer (ip) , dimension (dimarray)           :: coor_send_W , coor_recv_W , &
  coor_send_E , coor_recv_E , &
  coor_send_S , coor_recv_S , &
  coor_send_N , coor_recv_N , &
  coor_send_B , coor_recv_B , &
  coor_send_F , coor_recv_F

  integer (ip) , dimension (dimone)             :: coor_send_one_W , coor_recv_one_W , &
  coor_send_one_E , coor_recv_one_E , &
  coor_send_one_S , coor_recv_one_S , &
  coor_send_one_N , coor_recv_one_N , &
  coor_send_one_B , coor_recv_one_B , &
  coor_send_one_F , coor_recv_one_F

  integer (ip) , dimension (dimone)             :: coor_send_edge_xyp1 , coor_recv_edge_xyp1 , &
  coor_send_edge_xyn1 , coor_recv_edge_xyn1 , &
  coor_send_edge_xzp1 , coor_recv_edge_xzp1 , &
  coor_send_edge_xzn1 , coor_recv_edge_xzn1 , &
  coor_send_edge_yzp1 , coor_recv_edge_yzp1 , &
  coor_send_edge_yzn1 , coor_recv_edge_yzn1

  integer (ip) , dimension (dimone)             :: coor_send_edge_xyp2 , coor_recv_edge_xyp2 , &
  coor_send_edge_xyn2 , coor_recv_edge_xyn2 , &
  coor_send_edge_xzp2 , coor_recv_edge_xzp2 , &
  coor_send_edge_xzn2 , coor_recv_edge_xzn2 , &
  coor_send_edge_yzp2 , coor_recv_edge_yzp2 , &
  coor_send_edge_yzn2 , coor_recv_edge_yzn2


  integer (ip) , dimension (dimarray)           :: size , subsize1 , subsize2 , subsize3
  integer (ip) , dimension (dimone)             :: sizeone , subsizeone1 , subsizeone2 , subsizeone3

  integer (ip)                                  :: mpicode !< MPI default variable


  contains


  !> \brief Integer to string conversion
  !!
  !!
  !> \author
  !! Modeling, Simulation and Data Analysis (MSDA)\n
  !! Université Mohammed VI Polytechnique (UM6P)\n
  !! Benguerir, Morocco

  character(len=20) function str(k)
  integer, intent(in) :: k
  write (str, *) k
  str = adjustl(str)
end function str


!> \brief MPI initialization.
!!
!! Starts the MPI environment.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine init_mpi


 call MPI_INIT (mpicode)

 call MPI_COMM_RANK ( MPI_COMM_WORLD , rank , mpicode )

 call MPI_COMM_SIZE ( MPI_COMM_WORLD , nproc , mpicode )


end subroutine init_mpi


!> \brief MPI abortion.
!!
!! Aborts the MPI environment.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine abort_mpi (message)


 character (len=*) , intent (in) :: message
 integer :: dummy


 write (*,*) message
 call mpi_abort ( MPI_COMM_WORLD , dummy , mpicode )


end subroutine abort_mpi


!> \brief MPI termination.
!!
!! Finalizes the MPI environment.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine end_mpi


 call MPI_FINALIZE (mpicode)

 if ( mpicode /= 0 ) then
    write (*,675) 'Processor',rank,'/',nproc,'has NOT been closed correctly'
    675 format (1X,A9,I6,A2,I6,1X,A20)
 end if


end subroutine end_mpi


!> \brief MPI basic topology.
!!
!! Constructs the basic topology for parallel simulations.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine topology


 logical , parameter          :: reorganisation = .false.
 character (len_default)      :: word
 integer (ip)                 :: i , j , k , ok , ndir

 character (len_default) , parameter :: format_exit = '(3(1X,A,I7))'

 ! dims(:) = 0

  write(*,*) 'unit_grid', unit_grid
  write(*,*) 'file_grid', file_grid

 open ( unit_grid , file = file_grid , status = 'old' , iostat = ok )
 if ( ok /= 0 ) call abort_mpi ( 'error opening ' // trim (file_grid) )

 ! read (unit_grid,*) dims(3) ! number of processus in x-direction (0 if automatic)
 ! read (unit_grid,*) dims(2) ! number of processus in y-direction (0 if automatic)
 ! read (unit_grid,*) dims(1) ! number of processus in z-direction (0 if automatic)
 ! do i = 1 , nneigh
 !    read (unit_grid,*) bc(i) ! +/-(x,y,z) boundary
 ! end do
 ! read (unit_grid,*)

 ! count number of point in each direction
 do
    read (unit_grid,*) word
    if (word == 'X-DIRECTION') exit
 end do

 ! x-direction
 i = 0
 ndir = 1
 do
    read (unit_grid,*) word
    if (word == 'Y-DIRECTION') exit
    i = i + 1
 end do
 ntx = i ! count number of point in x-direction

 ! y-direction
 j = 0
 ndir = ndir +1
 do
    read (unit_grid,*) word
    if (word == 'Z-DIRECTION') exit
    j = j + 1
 end do
 nty = j ! count number of point in y-direction

 ! z-direction
 k = 0
 ndir = ndir +1
 do
    read (unit_grid,*) word
    if (word == 'END') exit
    k = k + 1
 end do
 ntz = k ! count number of point in z-direction

 close (unit_grid)

 if (ndir /= 3) call abort_mpi ( 'MISSING DIRECTION IN GRID.DAT : NDIR = ' // trim(str(ndir)) )

 ! if we decide to manually specify the processes in each direction
 ! verify that the total number of the agree
 if ( dims(1) == 0 .and. dims(2) == 0 .and. dims(3) == 0 ) then
    ! here the decomposition is automatic
 else
    if ( dims(1) * dims(2) * dims(3) /= nproc .and. nproc == 1 ) then
       ! Do sequential post-processing without always changing grid.dat file

       if ( rank == rank_default ) then
          write (*,*) 'nproc incompatible: correction nproc = 1'
          write (*,*) 'MPI division: ' , dims(3),' * ',dims(2),' * ',dims(1)
          dims = 1
          write (*,*) ' -----------> ' , dims(3),' * ',dims(2),' * ',dims(1)
          write (*,*) 'nproc = ', nproc
       end if

    end if

    if ( dims(1) * dims(2) * dims(3) /= nproc ) then

       if ( rank == rank_default ) then
          write (*,*) 'nproc incompatible'
          write (*,*) 'MPI division: ', &
          dims(3),' * ',dims(2),' * ',dims(1),' = ', &
          dims(1) * dims(2) * dims(3)
          write (*,*) 'nproc = ', nproc
       end if

       call abort_mpi (' ')

    end if
 end if


 ! print*, '00000000000000000000000',dims 
 ! call exit(0)
 ! checking for periodic boundary conditions
 do i = 1 , 5 , 2
    if ( bc (i) == periodic ) then

       if ( bc(i+1) == bc(i) ) then
       else
          call abort_mpi ('boundary conditions must appear in pairs')
       end if

    end if
 end do


 period(:) = .false.
 if ( ( bc(W) == bc(E) ) .and. ( bc(W) == periodic ) ) then
    period(3) = .true.
 end if
 if ( ( bc(S) == bc(N) ) .and. ( bc(S) == periodic ) ) then
    period(2) = .true.
 end if
 if ( ( bc(B) == bc(F) ) .and. ( bc(B) == periodic ) ) then
    period(1) = .true.
 end if


 call MPI_DIMS_CREATE ( nproc , ndimcom , dims , mpicode )

 ! the vector dims gives the number of blocks in each direction.
 ! to first start splitting blocks in x-direction use this : x -> dims(3) , y -> dims(2) , z -> dims(1)
 ! by default the program starts splitting in z-direction:   x -> dims(1) , y -> dims(2) , z -> dims(3)

 comm3d = MPI_COMM_NULL
 call MPI_CART_CREATE ( MPI_COMM_WORLD , ndimcom , dims , period , reorganisation , comm3d , mpicode )



end subroutine topology


!> \brief MPI topology: degeneration.
!! Degenerates the basic topology, e.g. 3D->2D (parallelepipeds to planes) and 2D->1D (planes to lines).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine degenerate_topology


 logical , dimension (3) :: subdiv3


 comm2dxy = MPI_COMM_NULL ; comm2dxz = MPI_COMM_NULL ; comm2dyz = MPI_COMM_NULL
 comm1dx  = MPI_COMM_NULL ; comm1dy  = MPI_COMM_NULL ; comm1dz  = MPI_COMM_NULL


 if ( dims(3) == 1 .and. dims(2) == 1 .and. dims(1) == 1 ) then ! sequential problem

    ! do anything

 else if ( dims(3) > 1 .and. dims(2) > 1 .and. dims(1) > 1 ) then ! 3D pure topology

    ! degenerate into 2Dxy
    subdiv3 (3) = .true.  ! X
    subdiv3 (2) = .true.  ! Y
    subdiv3 (1) = .false. ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm2dxy , mpicode )

    ! degenerate into 2Dxz
    subdiv3 (3) = .true.  ! X
    subdiv3 (2) = .false. ! Y
    subdiv3 (1) = .true.  ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm2dxz , mpicode )

    ! degenerate into 2Dyz
    subdiv3 (3) = .false. ! X
    subdiv3 (2) = .true.  ! Y
    subdiv3 (1) = .true.  ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm2dyz , mpicode )

    ! degenerate into 1Dx
    subdiv3 (3) = .true.  ! X
    subdiv3 (2) = .false. ! Y
    subdiv3 (1) = .false. ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm1dx , mpicode )

    ! degenerate into 1Dy
    subdiv3 (3) = .false. ! X
    subdiv3 (2) = .true.  ! Y
    subdiv3 (1) = .false. ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm1dy , mpicode )

    ! degenerate into 1Dz
    subdiv3 (3) = .false. ! X
    subdiv3 (2) = .false. ! Y
    subdiv3 (1) = .true.  ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm1dz , mpicode )

 else if ( dims(3) > 1 .and. dims(2) > 1 .and. dims(1) == 1 ) then ! XY planes

    ! degenerate into 2Dxy
    subdiv3 (3) = .true.  ! X
    subdiv3 (2) = .true.  ! Y
    subdiv3 (1) = .false. ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm2dxy , mpicode )

    ! degenerate into 1Dx
    subdiv3 (3) = .true.  ! X
    subdiv3 (2) = .false. ! Y
    subdiv3 (1) = .false. ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm1dx , mpicode )

    ! degenerate into 1Dy
    subdiv3 (3) = .false. ! X
    subdiv3 (2) = .true.  ! Y
    subdiv3 (1) = .false. ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm1dy , mpicode )

 else if ( dims(3) > 1 .and. dims(2) == 1 .and. dims(1) > 1 ) then ! XZ planes

    ! degenerate into 2Dxz
    subdiv3 (3) = .true.  ! X
    subdiv3 (2) = .false. ! Y
    subdiv3 (1) = .true.  ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm2dxz , mpicode )

    ! degenerate into 1Dx
    subdiv3 (3) = .true.  ! X
    subdiv3 (2) = .false. ! Y
    subdiv3 (1) = .false. ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm1dx , mpicode )

    ! degenerate into 1Dz
    subdiv3 (3) = .false. ! X
    subdiv3 (2) = .false. ! Y
    subdiv3 (1) = .true.  ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm1dz , mpicode )

 else if ( dims(3) == 1 .and. dims(2) > 1 .and. dims(1) > 1 ) then ! YZ planes

    ! degenerate into 2Dyz
    subdiv3 (3) = .false. ! X
    subdiv3 (2) = .true.  ! Y
    subdiv3 (1) = .true.  ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm2dyz , mpicode )

    ! degenerate into 1Dy
    subdiv3 (3) = .false. ! X
    subdiv3 (2) = .true.  ! Y
    subdiv3 (1) = .false. ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm1dy , mpicode )

    ! degenerate into 1Dz
    subdiv3 (3) = .false. ! X
    subdiv3 (2) = .false. ! Y
    subdiv3 (1) = .true.  ! Z

    call MPI_CART_SUB ( comm3d , subdiv3 , comm1dz , mpicode )

 else

    if ( rank == rank_default ) write (*,*) 'WARNING: no subcartesian topology to be created!'

 end if


end subroutine degenerate_topology


!> \brief Domain creation.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine domain


 call MPI_CART_COORDS ( comm3d , rank , ndimcom , coords , mpicode )

 sx = ( coords(3)*ntx ) / dims(3) + 1
 ex = ( ( coords(3)+1 )*ntx ) / dims(3)

 sy = ( coords(2)*nty ) / dims(2) + 1
 ey = ( ( coords(2)+1 )*nty ) / dims(2)

 sz = ( coords(1)*ntz ) / dims(1) + 1
 ez = ( ( coords(1)+1 )*ntz ) / dims(1)

 nxmax = ceiling ( 1.0_dp * ntx / dims (3) )
 nymax = ceiling ( 1.0_dp * nty / dims (2) )
 nzmax = ceiling ( 1.0_dp * ntz / dims (1) )

 subdom_size(1) = ex - sx + 1
 subdom_size(2) = ey - sy + 1
 subdom_size(3) = ez - sz + 1

 subdom_offset(1) = coords(3) * ntx / dims(3)
 subdom_offset(2) = coords(2) * nty / dims(2)
 subdom_offset(3) = coords(1) * ntz / dims(1)

end subroutine domain


!> \brief Neighbour calculation for each process.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine neighborhood


 neigh(:) = MPI_PROC_NULL
 neighdxy(:) = MPI_PROC_NULL ; neighdxz(:) = MPI_PROC_NULL ; neighdyz(:) = MPI_PROC_NULL

 call MPI_CART_SHIFT ( comm3d , 2 , 1 , neigh(W) , neigh(E) , mpicode )
 call MPI_CART_SHIFT ( comm3d , 1 , 1 , neigh(S) , neigh(N) , mpicode )
 call MPI_CART_SHIFT ( comm3d , 0 , 1 , neigh(B) , neigh(F) , mpicode )

 ! looking for diagonal neighboors

 ! XY plane
 ! diagonal 1
 dcoords = (/ coords(1) , coords(2)+1 , coords(3)+1 /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(2) <= dims(2)-1 ) .and. ( dcoords(3) <= dims(3)-1 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdxy (1) , mpicode )
 end if
 ! diagonal 2
 dcoords = (/ coords(1) , coords(2)-1 , coords(3)+1 /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(2) >= 0 ) .and. ( dcoords(3) <= dims(3)-1 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdxy (2) , mpicode )
 end if
 ! diagonal 3
 dcoords = (/ coords(1) , coords(2)-1 , coords(3)-1 /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(2) >= 0 ) .and. ( dcoords(3) >= 0 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdxy (3) , mpicode )
 end if
 ! diagonal 4
 dcoords = (/ coords(1) , coords(2)+1 , coords(3)-1 /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(2) <= dims(2)-1 ) .and. ( dcoords(3) >= 0 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdxy (4) , mpicode )
 end if

 ! XZ plane
 ! diagonal 1
 dcoords = (/ coords(1)+1 , coords(2) , coords(3)+1 /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(1) <= dims(1)-1 ) .and. ( dcoords(3) <= dims(3)-1 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdxz (1) , mpicode )
 end if
 ! diagonal 2
 dcoords = (/ coords(1)-1 , coords(2) , coords(3)+1 /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(1) >= 0 ) .and. ( dcoords(3) <= dims(3)-1 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdxz (2) , mpicode )
 end if
 ! diagonal 3
 dcoords = (/ coords(1)-1 , coords(2) , coords(3)-1 /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(1) >= 0 ) .and. ( dcoords(3) >= 0 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdxz (3) , mpicode )
 end if
 ! diagonal 4
 dcoords = (/ coords(1)+1 , coords(2) , coords(3)-1 /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(1) <= dims(1)-1 ) .and. ( dcoords(3) >= 0 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdxz (4) , mpicode )
 end if

 ! YZ plane
 ! diagonal 1
 dcoords = (/ coords(1)+1 , coords(2)+1 , coords(3) /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(1) <= dims(1)-1 ) .and. ( dcoords(2) <= dims(2)-1 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdyz (1) , mpicode )
 end if
 ! diagonal 2
 dcoords = (/ coords(1)-1 , coords(2)+1 , coords(3) /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(1) >= 0 ) .and. ( dcoords(2) <= dims(2)-1 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdyz (2) , mpicode )
 end if
 ! diagonal 3
 dcoords = (/ coords(1)-1 , coords(2)-1 , coords(3) /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(1) >= 0 ) .and. ( dcoords(2) >= 0 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdyz (3) , mpicode )
 end if
 ! diagonal 4
 dcoords = (/ coords(1)+1 , coords(2)-1 , coords(3) /) ! (/ zdirection , ydirection, xdirection /)
 if ( ( dcoords(1) <= dims(1)-1 ) .and. ( dcoords(2) >= 0 ) ) then
    call MPI_CART_RANK ( comm3d , dcoords , neighdyz (4) , mpicode )
 end if


 ! redefine the boundaries
 ! for peridicity the OpenMPI approach is not enough
 if ( period(3) ) then
    if ( sx /= 1   ) bc (W) = inner
    if ( ex /= ntx ) bc (E) = inner
 else
    if ( neigh (W) /= MPI_PROC_NULL ) bc (W) = inner
    if ( neigh (E) /= MPI_PROC_NULL ) bc (E) = inner
 end if

 if ( period(2) ) then
    if ( sy /= 1   ) bc (S) = inner
    if ( ey /= nty ) bc (N) = inner
 else
    if ( neigh (S) /= MPI_PROC_NULL ) bc (S) = inner
    if ( neigh (N) /= MPI_PROC_NULL ) bc (N) = inner
 end if

 if ( period(1) ) then
    if ( sz /= 1   ) bc (B) = inner
    if ( ez /= ntz ) bc (F) = inner
 else
    if ( neigh (B) /= MPI_PROC_NULL ) bc (B) = inner
    if ( neigh (F) /= MPI_PROC_NULL ) bc (F) = inner
 end if


end subroutine neighborhood


!> \brief Defintion of the communication types.
!!
!! Only conserved (independent) variables.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine type_conserved


 size = (/ ng+ex-sx+1+ng , ng+ey-sy+1+ng , ng+ez-sz+1+ng , nv /)

 ! create matrix data type to communicate on vertical Oxy plane

 subsize1 = (/ ex-sx+1 , ey-sy+1 , ng , nv /)

 coor_send_B = (/ ng , ng , ng         , 0 /) ! (/ sx   , sy , sz      , 1 /)
 coor_recv_F = (/ ng , ng , ez-sz+1+ng , 0 /) ! (/ sx   , sy , ez+1    , 1 /)

 coor_send_F = (/ ng , ng , ez-sz+1    , 0 /) ! (/ sx   , sy , ez-ng+1 , 1 /)
 coor_recv_B = (/ ng , ng , 0          , 0 /) ! (/ sx   , sy , sz-ng   , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_send_F , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_F , mpicode )
 call MPI_TYPE_COMMIT ( type_send_F , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_recv_B , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_B , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_send_B , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_B , mpicode )
 call MPI_TYPE_COMMIT ( type_send_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_recv_F , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_F , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_F , mpicode )


 ! create matrix data type to communicate on vertical Oxz plane

 subsize2 = (/ ex-sx+1 , ng , ez-sz+1 , nv /)

 coor_send_S = (/ ng , ng         , ng , 0 /) ! (/ sx , sy      , sz , 1 /)
 coor_recv_N = (/ ng , ey-sy+1+ng , ng , 0 /) ! (/ sx , ey+1    , sz , 1 /)

 coor_send_N = (/ ng , ey-sy+1    , ng , 0 /) ! (/ sx , ey-ng+1 , sz , 1 /)
 coor_recv_S = (/ ng , 0          , ng , 0 /) ! (/ sx , sy-ng   , sz , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_send_N , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_N , mpicode )
 call MPI_TYPE_COMMIT ( type_send_N , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_recv_S , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_S , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_send_S , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_S , mpicode )
 call MPI_TYPE_COMMIT ( type_send_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_recv_N , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_N , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_N , mpicode )


 ! create matrix data type to communicate on vertical Oyz plane

 subsize3 = (/ ng , ey-sy+1 , ez-sz+1 , nv /)

 coor_send_W = (/ ng         , ng , ng , 0 /) ! (/ sx      , sy , sz , 1 /)
 coor_recv_E = (/ ex-sx+1+ng , ng , ng , 0 /) ! (/ ex+1    , sy , sz , 1 /)

 coor_send_E = (/ ex-sx+1    , ng , ng , 0 /) ! (/ ex-ng+1 , sy , sz , 1 /)
 coor_recv_W = (/ 0          , ng , ng , 0 /) ! (/ sx-ng   , sy , sz , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_send_E , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_E , mpicode )
 call MPI_TYPE_COMMIT ( type_send_E , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_recv_W , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_W , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_send_W , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_W , mpicode )
 call MPI_TYPE_COMMIT ( type_send_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_recv_E , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_E , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_E , mpicode )


end subroutine type_conserved


!> \brief Defintion of the communication types.
!!
!! First-derivative variables for viscous problems.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine type_derivative


 size = (/ ng+ex-sx+1+ng , ng+ey-sy+1+ng , ng+ez-sz+1+ng , nderiv /)

 ! create matrix data type to communicate on vertical Oxy plan

 subsize1 = (/ ex-sx+1 , ey-sy+1 , ng , nderiv /)

 coor_send_B = (/ ng , ng , ng         , 0 /) ! (/ sx   , sy , sz      , 1 /)
 coor_recv_F = (/ ng , ng , ez-sz+1+ng , 0 /) ! (/ sx   , sy , ez+1    , 1 /)

 coor_send_F = (/ ng , ng , ez-sz+1    , 0 /) ! (/ sx   , sy , ez-ng+1 , 1 /)
 coor_recv_B = (/ ng , ng , 0          , 0 /) ! (/ sx   , sy , sz-ng   , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_send_F , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_deriv_F , mpicode )
 call MPI_TYPE_COMMIT ( type_send_deriv_F , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_recv_B , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_deriv_B , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_deriv_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_send_B , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_deriv_B , mpicode )
 call MPI_TYPE_COMMIT ( type_send_deriv_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_recv_F , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_deriv_F , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_deriv_F , mpicode )


 ! create matrix data type to communicate on vertical Oxz plan

 subsize2 = (/ ex-sx+1 , ng , ez-sz+1 , nderiv /)

 coor_send_S = (/ ng , ng         , ng , 0 /) ! (/ sx , sy      , sz , 1 /)
 coor_recv_N = (/ ng , ey-sy+1+ng , ng , 0 /) ! (/ sx , ey+1    , sz , 1 /)

 coor_send_N = (/ ng , ey-sy+1    , ng , 0 /) ! (/ sx , ey-ng+1 , sz , 1 /)
 coor_recv_S = (/ ng , 0          , ng , 0 /) ! (/ sx , sy-ng   , sz , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_send_N , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_deriv_N , mpicode )
 call MPI_TYPE_COMMIT ( type_send_deriv_N , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_recv_S , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_deriv_S , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_deriv_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_send_S , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_deriv_S , mpicode )
 call MPI_TYPE_COMMIT ( type_send_deriv_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_recv_N , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_deriv_N , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_deriv_N , mpicode )


 ! create matrix data type to communicate on vertical Oyz plan

 subsize3 = (/ ng , ey-sy+1 , ez-sz+1 , nderiv /)

 coor_send_W = (/ ng         , ng , ng , 0 /) ! (/ sx      , sy , sz , 1 /)
 coor_recv_E = (/ ex-sx+1+ng , ng , ng , 0 /) ! (/ ex+1    , sy , sz , 1 /)

 coor_send_E = (/ ex-sx+1    , ng , ng , 0 /) ! (/ ex-ng+1 , ey , sz , 1 /)
 coor_recv_W = (/ 0          , ng , ng , 0 /) ! (/ sx-ng   , sy , sz , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_send_E , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_deriv_E , mpicode )
 call MPI_TYPE_COMMIT ( type_send_deriv_E , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_recv_W , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_deriv_W , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_deriv_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_send_W , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_deriv_W , mpicode )
 call MPI_TYPE_COMMIT ( type_send_deriv_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_recv_E , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_deriv_E , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_deriv_E , mpicode )


end subroutine type_derivative


!> \brief Defintion of the communication types.
!!
!! First-derivative variables of velocity.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine type_derivVel


 size = (/ ng+ex-sx+1+ng , ng+ey-sy+1+ng , ng+ez-sz+1+ng , ndim * ndim  /)

 ! create matrix data type to communicate on vertical Oxy plan

 subsize1 = (/ ex-sx+1 , ey-sy+1 , ng , ndim * ndim  /)

 coor_send_B = (/ ng , ng , ng         , 0 /) ! (/ sx   , sy , sz      , 1 /)
 coor_recv_F = (/ ng , ng , ez-sz+1+ng , 0 /) ! (/ sx   , sy , ez+1    , 1 /)

 coor_send_F = (/ ng , ng , ez-sz+1    , 0 /) ! (/ sx   , sy , ez-ng+1 , 1 /)
 coor_recv_B = (/ ng , ng , 0          , 0 /) ! (/ sx   , sy , sz-ng   , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_send_F , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_send_derivVel_F , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivVel_F , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_recv_B , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_recv_derivVel_B , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivVel_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_send_B , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_send_derivVel_B , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivVel_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_recv_F , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_recv_derivVel_F , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivVel_F , mpicode )


 ! create matrix data type to communicate on vertical Oxz plan

 subsize2 = (/ ex-sx+1 , ng , ez-sz+1 , ndim * ndim  /)

 coor_send_S = (/ ng , ng         , ng , 0 /) ! (/ sx , sy      , sz , 1 /)
 coor_recv_N = (/ ng , ey-sy+1+ng , ng , 0 /) ! (/ sx , ey+1    , sz , 1 /)

 coor_send_N = (/ ng , ey-sy+1    , ng , 0 /) ! (/ sx , ey-ng+1 , sz , 1 /)
 coor_recv_S = (/ ng , 0          , ng , 0 /) ! (/ sx , sy-ng   , sz , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_send_N , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_send_derivVel_N , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivVel_N , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_recv_S , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_recv_derivVel_S , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivVel_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_send_S , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_send_derivVel_S , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivVel_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_recv_N , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_recv_derivVel_N , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivVel_N , mpicode )


 ! create matrix data type to communicate on vertical Oyz plan

 subsize3 = (/ ng , ey-sy+1 , ez-sz+1 , ndim * ndim  /)

 coor_send_W = (/ ng         , ng , ng , 0 /) ! (/ sx      , sy , sz , 1 /)
 coor_recv_E = (/ ex-sx+1+ng , ng , ng , 0 /) ! (/ ex+1    , sy , sz , 1 /)

 coor_send_E = (/ ex-sx+1    , ng , ng , 0 /) ! (/ ex-ng+1 , ey , sz , 1 /)
 coor_recv_W = (/ 0          , ng , ng , 0 /) ! (/ sx-ng   , sy , sz , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_send_E , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_send_derivVel_E , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivVel_E , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_recv_W , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_recv_derivVel_W , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivVel_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_send_W , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_send_derivVel_W , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivVel_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_recv_E , &
   &                              MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   &                              type_recv_derivVel_E , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivVel_E , mpicode )


end subroutine type_derivVel

!> \brief Defintion of the communication types.
!!
!! Additional third-derivative variables for LES viscous problems.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine type_derivative_SGS


 size = (/ ng+ex-sx+1+ng , ng+ey-sy+1+ng , ng+ez-sz+1+ng , nderivSGS /)

 ! create matrix data type to communicate on vertical Oxy plan

 subsize1 = (/ ex-sx+1 , ey-sy+1 , ng , nderivSGS /)

 coor_send_B = (/ ng , ng , ng         , 0 /) ! (/ sx   , sy , sz      , 1 /)
 coor_recv_F = (/ ng , ng , ez-sz+1+ng , 0 /) ! (/ sx   , sy , ez+1    , 1 /)

 coor_send_F = (/ ng , ng , ez-sz+1    , 0 /) ! (/ sx   , sy , ez-ng+1 , 1 /)
 coor_recv_B = (/ ng , ng , 0          , 0 /) ! (/ sx   , sy , sz-ng   , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_send_F , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_derivSGS_F , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivSGS_F , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_recv_B , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_derivSGS_B , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivSGS_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_send_B , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_derivSGS_B , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivSGS_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize1 , coor_recv_F , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_derivSGS_F , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivSGS_F , mpicode )


 ! create matrix data type to communicate on vertical Oxz plan

 subsize2 = (/ ex-sx+1 , ng , ez-sz+1 , nderivSGS /)

 coor_send_S = (/ ng , ng         , ng , 0 /) ! (/ sx , sy      , sz , 1 /)
 coor_recv_N = (/ ng , ey-sy+1+ng , ng , 0 /) ! (/ sx , ey+1    , sz , 1 /)

 coor_send_N = (/ ng , ey-sy+1    , ng , 0 /) ! (/ sx , ey-ng+1 , sz , 1 /)
 coor_recv_S = (/ ng , 0          , ng , 0 /) ! (/ sx , sy-ng   , sz , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_send_N , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_derivSGS_N , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivSGS_N , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_recv_S , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_derivSGS_S , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivSGS_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_send_S , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_derivSGS_S , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivSGS_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize2 , coor_recv_N , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_derivSGS_N , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivSGS_N , mpicode )


 ! create matrix data type to communicate on vertical Oyz plan

 subsize3 = (/ ng , ey-sy+1 , ez-sz+1 , nderivSGS /)

 coor_send_W = (/ ng         , ng , ng , 0 /) ! (/ sx      , sy , sz , 1 /)
 coor_recv_E = (/ ex-sx+1+ng , ng , ng , 0 /) ! (/ ex+1    , sy , sz , 1 /)

 coor_send_E = (/ ex-sx+1    , ng , ng , 0 /) ! (/ ex-ng+1 , ey , sz , 1 /)
 coor_recv_W = (/ 0          , ng , ng , 0 /) ! (/ sx-ng   , sy , sz , 1 /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_send_E , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_derivSGS_E , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivSGS_E , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_recv_W , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_derivSGS_W , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivSGS_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_send_W , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_send_derivSGS_W , mpicode )
 call MPI_TYPE_COMMIT ( type_send_derivSGS_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimarray , size , subsize3 , coor_recv_E , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION , &
   type_recv_derivSGS_E , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_derivSGS_E , mpicode )


end subroutine type_derivative_SGS


!> \brief Communication of an array of velocity tensor
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine comm_derivVel (veltensor)


 real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: veltensor !< array of derivative variables to be communicated

 integer (ip) , parameter  :: label=100 , nb_request = ndimmax * dimarray
 integer (ip) , dimension (MPI_STATUS_SIZE,nb_request) :: tab_stat
 integer (ip) , dimension (nb_request) :: request

 real (dp) :: tcomm_tmp

 tcomm_tmp = MPI_WTIME()


 !********* Est/West communication ************************************
 ! send my boundary to Est and receive from West
 call MPI_ISEND ( veltensor , 1 , type_send_derivVel_E , neigh(E) , &
   &               label , comm3d , request(1) , mpicode )
 call MPI_IRECV ( veltensor , 1 , type_recv_derivVel_W , neigh(W) , &
   &               label , comm3d , request(2) , mpicode )

 call MPI_ISEND ( veltensor , 1 , type_send_derivVel_W , neigh(W) , &
   &               label , comm3d , request(3) , mpicode )
 call MPI_IRECV ( veltensor , 1 , type_recv_derivVel_E , neigh(E) , &
   &               label , comm3d , request(4) , mpicode )


 !********* North/South communication ************************************
 ! send my boundary to North and receive from South
 call MPI_ISEND ( veltensor , 1 , type_send_derivVel_N , neigh(N) , &
   &               label , comm3d , request(5) , mpicode )
 call MPI_IRECV ( veltensor , 1 , type_recv_derivVel_S , neigh(S) , &
   &               label , comm3d , request(6) , mpicode )

 call MPI_ISEND ( veltensor , 1 , type_send_derivVel_S , neigh(S) , &
   &               label , comm3d , request(7) , mpicode )
 call MPI_IRECV ( veltensor , 1 , type_recv_derivVel_N , neigh(N) , &
   &               label , comm3d , request(8) , mpicode )


 !********* Behind/Front communication ************************************
 ! send my boundary to Front and receive from Behind
 call MPI_ISEND ( veltensor , 1 , type_send_derivVel_F , neigh(F) , &
   &               label , comm3d , request (9) , mpicode )
 call MPI_IRECV ( veltensor , 1 , type_recv_derivVel_B , neigh(B) , &
   &               label , comm3d , request(10) , mpicode )

 call MPI_ISEND ( veltensor , 1 , type_send_derivVel_B , neigh(B) , &
   &               label , comm3d , request(11) , mpicode )
 call MPI_IRECV ( veltensor , 1 , type_recv_derivVel_F , neigh(F) , &
   &               label , comm3d , request(12) , mpicode )

 call MPI_WAITALL ( nb_request , request , tab_stat , mpicode )

 tcomm = tcomm + MPI_WTIME() - tcomm_tmp


end subroutine comm_derivVel
!> \brief one-variable communication type.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine type_one


 sizeone = (/ ng+ex-sx+1+ng , ng+ey-sy+1+ng , ng+ez-sz+1+ng /)

 ! create matrix data type to communicate on vertical Oxy plan

 subsizeone1     = (/ ex-sx+1 , ey-sy+1 , ng /)

 coor_send_one_B = (/ ng , ng , ng         /) ! (/ sx   , sy , sz      /)
 coor_recv_one_F = (/ ng , ng , ez-sz+1+ng /) ! (/ sx   , sy , ez+1    /)

 coor_send_one_F = (/ ng , ng , ez-sz+1    /) ! (/ sx   , sy , ez-ng+1 /)
 coor_recv_one_B = (/ ng , ng , 0          /) ! (/ sx   , sy , sz-ng   /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_send_one_F , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_send_one_F , mpicode )
 call MPI_TYPE_COMMIT ( type_send_one_F , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_recv_one_B , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_recv_one_B , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_one_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_send_one_B , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_send_one_B , mpicode )
 call MPI_TYPE_COMMIT ( type_send_one_B , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_recv_one_F , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_recv_one_F , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_one_F , mpicode )


 ! create matrix data type to communicate on vertical Oxz plan

 subsizeone2     = (/ ex-sx+1 , ng , ez-sz+1 /)

 coor_send_one_S = (/ ng , ng         , ng /) ! (/ sx , sy      , sz /)
 coor_recv_one_N = (/ ng , ey-sy+1+ng , ng /) ! (/ sx , ey+1    , sz /)

 coor_send_one_N = (/ ng , ey-sy+1    , ng /) ! (/ sx , ey-ng+1 , sz /)
 coor_recv_one_S = (/ ng , 0          , ng /) ! (/ sx , sy-ng   , sz /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_send_one_N , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_send_one_N , mpicode )
 call MPI_TYPE_COMMIT ( type_send_one_N , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_recv_one_S , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_recv_one_S , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_one_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_send_one_S , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_send_one_S , mpicode )
 call MPI_TYPE_COMMIT ( type_send_one_S , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_recv_one_N , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_recv_one_N , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_one_N , mpicode )


 ! create matrix data type to communicate on vertical Oyz plan

 subsizeone3     = (/ ng , ey-sy+1 , ez-sz+1 /)

 coor_send_one_W = (/ ng         , ng , ng /) ! (/ sx      , sy , sz /)
 coor_recv_one_E = (/ ex-sx+1+ng , ng , ng /) ! (/ ex+1    , sy , sz /)

 coor_send_one_E = (/ ex-sx+1    , ng , ng /) ! (/ ex-ng+1 , sy , sz /)
 coor_recv_one_W = (/ 0          , ng , ng /) ! (/ sx-ng   , sy , sz /)

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_send_one_E , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_send_one_E , mpicode )
 call MPI_TYPE_COMMIT ( type_send_one_E , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_recv_one_W , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_recv_one_W , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_one_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_send_one_W , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_send_one_W , mpicode )
 call MPI_TYPE_COMMIT ( type_send_one_W , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_recv_one_E , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,         &
   type_recv_one_E , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_one_E , mpicode )


end subroutine type_one


!> \brief edge-variable communication type.
!!
!! Required for continuous plot when using parallel files in ParaView.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine type_edges


 sizeone = (/ ng+ex-sx+1+ng , ng+ey-sy+1+ng , ng+ez-sz+1+ng /)

 ! create matrix data type to communicate Z edge in the XY plane

 ! first diagonal : +x * +y and -x * -y

 subsizeone1 = (/ 1 , 1 , ez-sz+1+2 /)

 coor_send_edge_xyp1 = (/ ng         , ng         , ng-1  /) ! (/ sx   , sy   , sz-1  /)
 coor_recv_edge_xyn1 = (/ ex-sx+1+ng , ey-sy+1+ng , ng-1  /) ! (/ ex+1 , ey+1 , sz-1  /)

 coor_send_edge_xyn1 = (/ ex-sx+ng   , ey-sy+ng   , ng-1  /) ! (/ ex   , ey   , sz-1  /)
 coor_recv_edge_xyp1 = (/ ng-1       , ng-1       , ng-1  /) ! (/ sx-1 , sy-1 , sz-1  /)


 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_send_edge_xyp1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_xyp1 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_xyp1 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_send_edge_xyn1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_xyn1 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_xyn1 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_recv_edge_xyp1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_xyp1 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_xyp1 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_recv_edge_xyn1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_xyn1 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_xyn1 , mpicode )


 ! second diagonal : +x * -y and -x * +y

 subsizeone1 = (/ 1 , 1 , ez-sz+1+2 /)

 coor_send_edge_xyp2 = (/ ex-sx+ng   , ng         , ng-1  /) ! (/ ex   , sy   , sz-1  /)
 coor_recv_edge_xyn2 = (/ ng-1       , ey-sy+1+ng , ng-1  /) ! (/ sx-1 , ey+1 , sz-1  /)

 coor_send_edge_xyn2 = (/ ng         , ey-sy+ng   , ng-1  /) ! (/ sx   , ey   , sz-1  /)
 coor_recv_edge_xyp2 = (/ ex-sx+1+ng , ng-1       , ng-1  /) ! (/ ex+1 , sy-1 , sz-1  /)


 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_send_edge_xyp2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_xyp2 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_xyp2 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_send_edge_xyn2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_xyn2 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_xyn2 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_recv_edge_xyp2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_xyp2 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_xyp2 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone1 , coor_recv_edge_xyn2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_xyn2 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_xyn2 , mpicode )


 ! create matrix data type to communicate Y edge in the XZ plane

 ! first diagonal : +x * +z and -x * -z

 subsizeone2 = (/ 1 , ey-sy+1+2 , 1 /)

 coor_send_edge_xzp1 = (/ ng         , ng-1 , ng            /) ! (/ sx   , sy-1   , sz    /)
 coor_recv_edge_xzn1 = (/ ex-sx+1+ng , ng-1 , ez-sz+1+ng    /) ! (/ ex+1 , sy-1   , ez+1  /)

 coor_send_edge_xzn1 = (/ ex-sx+ng   , ng-1     , ez-sz+ng  /) ! (/ ex   , sy-1   , ez    /)
 coor_recv_edge_xzp1 = (/ ng-1       , ng-1     , ng-1      /) ! (/ sx-1 , sy-1   , sz-1  /)


 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_send_edge_xzp1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_xzp1 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_xzp1 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_send_edge_xzn1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_xzn1 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_xzn1 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_recv_edge_xzp1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_xzp1 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_xzp1 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_recv_edge_xzn1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_xzn1 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_xzn1 , mpicode )


 ! second diagonal : +x * -z and -x * +z

 subsizeone2 = (/ 1 , ey-sy+1+2 , 1 /)

 coor_send_edge_xzp2 = (/ ex-sx+ng   , ng-1 , ng            /) ! (/ ex   , sy-1   , sz    /)
 coor_recv_edge_xzn2 = (/ ng-1       , ng-1 , ez-sz+1+ng    /) ! (/ sx-1 , sy-1   , ez+1  /)

 coor_send_edge_xzn2 = (/ ng         , ng-1     , ez-sz+ng  /) ! (/ sx   , sy-1   , ez    /)
 coor_recv_edge_xzp2 = (/ ex-sx+1+ng , ng-1     , ng-1      /) ! (/ ex+1 , sy-1   , sz-1  /)


 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_send_edge_xzp2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_xzp2 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_xzp2 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_send_edge_xzn2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_xzn2 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_xzn2 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_recv_edge_xzp2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_xzp2 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_xzp2 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone2 , coor_recv_edge_xzn2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_xzn2 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_xzn2 , mpicode )


 ! create matrix data type to communicate X edge in the YZ plane

 ! first diagonal : +y * +z and -y * -z

 subsizeone3 = (/ ex-sx+1+2 , 1 , 1 /)

 coor_send_edge_yzp1 = (/ ng-1   , ng         , ng            /) ! (/ sx-1   , sy   , sz    /)
 coor_recv_edge_yzn1 = (/ ng-1   , ey-sy+1+ng , ez-sz+1+ng    /) ! (/ sx-1   , ey+1 , ez+1  /)

 coor_send_edge_yzn1 = (/ ng-1   , ey-sy+ng , ez-sz+ng        /) ! (/ sx-1   , ey   , ez    /)
 coor_recv_edge_yzp1 = (/ ng-1   , ng-1     , ng-1            /) ! (/ sx-1   , sy-1 , sz-1  /)


 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_send_edge_yzp1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_yzp1 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_yzp1 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_send_edge_yzn1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_yzn1 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_yzn1 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_recv_edge_yzp1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_yzp1 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_yzp1 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_recv_edge_yzn1 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_yzn1 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_yzn1 , mpicode )


 ! second diagonal : +y * -z and -y * +z

 subsizeone3 = (/ ex-sx+1+2 , 1 , 1 /)

 coor_send_edge_yzp2 = (/ ng-1   , ey-sy+ng   , ng            /) ! (/ sx-1   , ey   , sz    /)
 coor_recv_edge_yzn2 = (/ ng-1   , ng-1       , ez-sz+1+ng    /) ! (/ sx-1   , sy-1 , ez+1  /)

 coor_send_edge_yzn2 = (/ ng-1   , ng         , ez-sz+ng      /) ! (/ sx-1   , sy   , ez    /)
 coor_recv_edge_yzp2 = (/ ng-1   , ey-sy+1+ng , ng-1          /) ! (/ sx-1   , ey+1 , sz-1  /)


 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_send_edge_yzp2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_yzp2 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_yzp2 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_send_edge_yzn2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_send_edge_yzn2 , mpicode )
 call MPI_TYPE_COMMIT ( type_send_edge_yzn2 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_recv_edge_yzp2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_yzp2 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_yzp2 , mpicode )

 call MPI_TYPE_CREATE_SUBARRAY ( dimone , sizeone , subsizeone3 , coor_recv_edge_yzn2 , &
   MPI_ORDER_FORTRAN , MPI_DOUBLE_PRECISION ,             &
   type_recv_edge_yzn2 , mpicode )
 call MPI_TYPE_COMMIT ( type_recv_edge_yzn2 , mpicode )


end subroutine type_edges


!> \brief One-variable communication.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine comm_one (v)


 real (dp) , allocatable , dimension (:,:,:) , intent (inout) :: v !< variable to be communicated

 integer (ip) , parameter                               :: label=100 , nb_request = ndimmax * dimarray
 integer (ip) , dimension (MPI_STATUS_SIZE,nb_request)  :: tab_stat
 integer (ip) , dimension (nb_request)                  :: request


 real (dp) :: tcomm_tmp


 tcomm_tmp = MPI_WTIME()


 !********* Est/West communication ************************************
 ! send my boundary to Est and receive from West
 call MPI_ISEND ( v , 1 , type_send_one_E , neigh(E) , &
   label , comm3d , request(1) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_one_W , neigh(W) , &
   label , comm3d , request(2) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_one_W , neigh(W) , &
   label , comm3d , request(3) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_one_E , neigh(E) , &
   label , comm3d , request(4) , mpicode )


 !********* North/South communication ************************************
 ! send my boundary to North and receive from South
 call MPI_ISEND ( v , 1 , type_send_one_N , neigh(N) , &
   label , comm3d , request(5) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_one_S , neigh(S) , &
   label , comm3d , request(6) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_one_S , neigh(S) , &
   label , comm3d , request(7) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_one_N , neigh(N) , &
   label , comm3d , request(8) , mpicode )


 !********* Behind/Front communication ************************************
 ! send my boundary to Front and receive from Behind
 call MPI_ISEND ( v , 1 , type_send_one_F , neigh(F) , &
   label , comm3d , request (9) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_one_B , neigh(B) , &
   label , comm3d , request(10) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_one_B , neigh(B) , &
   label , comm3d , request(11) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_one_F , neigh(F) , &
   label , comm3d , request(12) , mpicode )


 call MPI_WAITALL ( nb_request , request , tab_stat , mpicode )


 tcomm = tcomm + MPI_WTIME() - tcomm_tmp


end subroutine comm_one


!> \brief Edge communication.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine comm_edges (v)


 real (dp) , allocatable , dimension (:,:,:) , intent (inout) :: v !< variable to be communicated


 integer (ip) , parameter                               :: label=100 , nb_request = ndimmax * dimarray
 integer (ip) , dimension (MPI_STATUS_SIZE,nb_request)  :: tab_stat
 integer (ip) , dimension (nb_request)                  :: request


 real (dp) :: tcomm_tmp


 tcomm_tmp = MPI_WTIME()


 ! first diagonal

 ! Z edge in XY plane
 call MPI_ISEND ( v , 1 , type_send_edge_xyn1 , neighdxy(1) , &
   label , comm3d , request(1) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_xyp1 , neighdxy(3) , &
   label , comm3d , request(2) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_edge_xyp1 , neighdxy(3) , &
   label , comm3d , request(3) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_xyn1 , neighdxy(1) , &
   label , comm3d , request(4) , mpicode )


 ! Y edge in XZ plane
 call MPI_ISEND ( v , 1 , type_send_edge_xzn1 , neighdxz(1) , &
   label , comm3d , request(5) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_xzp1 , neighdxz(3) , &
   label , comm3d , request(6) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_edge_xzp1 , neighdxz(3) , &
   label , comm3d , request(7) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_xzn1 , neighdxz(1) , &
   label , comm3d , request(8) , mpicode )


 ! X edge in YZ plane
 call MPI_ISEND ( v , 1 , type_send_edge_yzn1 , neighdyz(1) , &
   label , comm3d , request(9) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_yzp1 , neighdyz(3) , &
   label , comm3d , request(10) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_edge_yzp1 , neighdyz(3) , &
   label , comm3d , request(11) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_yzn1 , neighdyz(1) , &
   label , comm3d , request(12) , mpicode )


 call MPI_WAITALL ( nb_request , request , tab_stat , mpicode )


 ! second diagonal

 ! Z edge in XY plane
 call MPI_ISEND ( v , 1 , type_send_edge_xyn2 , neighdxy(4) , &
   label , comm3d , request(1) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_xyp2 , neighdxy(2) , &
   label , comm3d , request(2) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_edge_xyp2 , neighdxy(2) , &
   label , comm3d , request(3) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_xyn2 , neighdxy(4) , &
   label , comm3d , request(4) , mpicode )


 ! Y edge in XZ plane
 call MPI_ISEND ( v , 1 , type_send_edge_xzn2 , neighdxz(4) , &
   label , comm3d , request(5) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_xzp2 , neighdxz(2) , &
   label , comm3d , request(6) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_edge_xzp2 , neighdxz(2) , &
   label , comm3d , request(7) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_xzn2 , neighdxz(4) , &
   label , comm3d , request(8) , mpicode )


 ! X edge in YZ plane
 call MPI_ISEND ( v , 1 , type_send_edge_yzn2 , neighdyz(4) , &
   label , comm3d , request(9) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_yzp2 , neighdyz(2) , &
   label , comm3d , request(10) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_edge_yzp2 , neighdyz(2) , &
   label , comm3d , request(11) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_edge_yzn2 , neighdyz(4) , &
   label , comm3d , request(12) , mpicode )


 call MPI_WAITALL ( nb_request , request , tab_stat , mpicode )


 tcomm = tcomm + MPI_WTIME() - tcomm_tmp


end subroutine comm_edges


!> \brief Communication of an array of conserved variables.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine comm_cons (v)


 real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v !< array of conserved variables to be communicated


 integer (ip) , parameter                               :: label=100 , nb_request = ndimmax * dimarray
 integer (ip) , dimension (MPI_STATUS_SIZE,nb_request)  :: tab_stat
 integer (ip) , dimension (nb_request)                  :: request


 real (dp) :: tcomm_tmp


 tcomm_tmp = MPI_WTIME()


 !********* Est/West communication ************************************
 ! send my boundary to Est and receive from West
 call MPI_ISEND ( v , 1 , type_send_E , neigh(E) , &
   label , comm3d , request(1) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_W , neigh(W) , &
   label , comm3d , request(2) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_W , neigh(W) , &
   label , comm3d , request(3) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_E , neigh(E) , &
   label , comm3d , request(4) , mpicode )


 !********* North/South communication ************************************
 ! send my boundary to North and receive from South
 call MPI_ISEND ( v , 1 , type_send_N , neigh(N) , &
   label , comm3d , request(5) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_S , neigh(S) , &
   label , comm3d , request(6) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_S , neigh(S) , &
   label , comm3d , request(7) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_N , neigh(N) , &
   label , comm3d , request(8) , mpicode )


 !********* Behind/Front communication ************************************
 ! send my boundary to Front and receive from Behind
 call MPI_ISEND ( v , 1 , type_send_F , neigh(F) , &
   label , comm3d , request (9) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_B , neigh(B) , &
   label , comm3d , request(10) , mpicode )

 call MPI_ISEND ( v , 1 , type_send_B , neigh(B) , &
   label , comm3d , request(11) , mpicode )
 call MPI_IRECV ( v , 1 , type_recv_F , neigh(F) , &
   label , comm3d , request(12) , mpicode )


 call MPI_WAITALL ( nb_request , request , tab_stat , mpicode )


 tcomm = tcomm + MPI_WTIME() - tcomm_tmp


end subroutine comm_cons


!> \brief Communication of an array of derivative variables.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine comm_deriv (dv)


 real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: dv !< array of derivative variables to be communicated


 integer (ip) , parameter                               :: label=100 , nb_request = ndimmax * dimarray
 integer (ip) , dimension (MPI_STATUS_SIZE,nb_request)  :: tab_stat
 integer (ip) , dimension (nb_request)                  :: request


 real (dp) :: tcomm_tmp


 tcomm_tmp = MPI_WTIME()


 !********* Est/West communication ************************************
 ! send my boundary to Est and receive from West
 call MPI_ISEND ( dv , 1 , type_send_deriv_E , neigh(E) , &
   label , comm3d , request(1) , mpicode )
 call MPI_IRECV ( dv , 1 , type_recv_deriv_W , neigh(W) , &
   label , comm3d , request(2) , mpicode )

 call MPI_ISEND ( dv , 1 , type_send_deriv_W , neigh(W) , &
   label , comm3d , request(3) , mpicode )
 call MPI_IRECV ( dv , 1 , type_recv_deriv_E , neigh(E) , &
   label , comm3d , request(4) , mpicode )


 !********* North/South communication ************************************
 ! send my boundary to North and receive from South
 call MPI_ISEND ( dv , 1 , type_send_deriv_N , neigh(N) , &
   label , comm3d , request(5) , mpicode )
 call MPI_IRECV ( dv , 1 , type_recv_deriv_S , neigh(S) , &
   label , comm3d , request(6) , mpicode )

 call MPI_ISEND ( dv , 1 , type_send_deriv_S , neigh(S) , &
   label , comm3d , request(7) , mpicode )
 call MPI_IRECV ( dv , 1 , type_recv_deriv_N , neigh(N) , &
   label , comm3d , request(8) , mpicode )


 !********* Behind/Front communication ************************************
 ! send my boundary to Front and receive from Behind
 call MPI_ISEND ( dv , 1 , type_send_deriv_F , neigh(F) , &
   label , comm3d , request (9) , mpicode )
 call MPI_IRECV ( dv , 1 , type_recv_deriv_B , neigh(B) , &
   label , comm3d , request(10) , mpicode )

 call MPI_ISEND ( dv , 1 , type_send_deriv_B , neigh(B) , &
   label , comm3d , request(11) , mpicode )
 call MPI_IRECV ( dv , 1 , type_recv_deriv_F , neigh(F) , &
   label , comm3d , request(12) , mpicode )


 call MPI_WAITALL ( nb_request , request , tab_stat , mpicode )


 tcomm = tcomm + MPI_WTIME() - tcomm_tmp


end subroutine comm_deriv


!> \brief Communication of an array of additional SGS derivative variables.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine comm_derivSGS (dvSGS)


 real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: dvSGS !< array of derivative variables to be communicated


 integer (ip) , parameter                               :: label=100 , nb_request = ndimmax * dimarray
 integer (ip) , dimension (MPI_STATUS_SIZE,nb_request)  :: tab_stat
 integer (ip) , dimension (nb_request)                  :: request


 real (dp) :: tcomm_tmp


 tcomm_tmp = MPI_WTIME()


 !********* Est/West communication ************************************
 ! send my boundary to Est and receive from West
 call MPI_ISEND ( dvSGS , 1 , type_send_derivSGS_E , neigh(E) , &
   label , comm3d , request(1) , mpicode )
 call MPI_IRECV ( dvSGS , 1 , type_recv_derivSGS_W , neigh(W) , &
   label , comm3d , request(2) , mpicode )

 call MPI_ISEND ( dvSGS , 1 , type_send_derivSGS_W , neigh(W) , &
   label , comm3d , request(3) , mpicode )
 call MPI_IRECV ( dvSGS , 1 , type_recv_derivSGS_E , neigh(E) , &
   label , comm3d , request(4) , mpicode )


 !********* North/South communication ************************************
 ! send my boundary to North and receive from South
 call MPI_ISEND ( dvSGS , 1 , type_send_derivSGS_N , neigh(N) , &
   label , comm3d , request(5) , mpicode )
 call MPI_IRECV ( dvSGS , 1 , type_recv_derivSGS_S , neigh(S) , &
   label , comm3d , request(6) , mpicode )

 call MPI_ISEND ( dvSGS , 1 , type_send_derivSGS_S , neigh(S) , &
   label , comm3d , request(7) , mpicode )
 call MPI_IRECV ( dvSGS , 1 , type_recv_derivSGS_N , neigh(N) , &
   label , comm3d , request(8) , mpicode )


 !********* Behind/Front communication ************************************
 ! send my boundary to Front and receive from Behind
 call MPI_ISEND ( dvSGS , 1 , type_send_derivSGS_F , neigh(F) , &
   label , comm3d , request (9) , mpicode )
 call MPI_IRECV ( dvSGS , 1 , type_recv_derivSGS_B , neigh(B) , &
   label , comm3d , request(10) , mpicode )

 call MPI_ISEND ( dvSGS , 1 , type_send_derivSGS_B , neigh(B) , &
   label , comm3d , request(11) , mpicode )
 call MPI_IRECV ( dvSGS , 1 , type_recv_derivSGS_F , neigh(F) , &
   label , comm3d , request(12) , mpicode )


 call MPI_WAITALL ( nb_request , request , tab_stat , mpicode )


 tcomm = tcomm + MPI_WTIME() - tcomm_tmp


end subroutine comm_derivSGS


!> \brief Domain selection.
!!
!! E.g. of subdomains: inner domain or ghost points in W, E, S, N, B or F directions (it allows to define values on ghost points).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


 integer (ip) , intent (in)    :: domain_id !< subdomain selection
 integer (ip) , intent (inout) :: i0 , i1   !< start/end points in x-direction
 integer (ip) , intent (inout) :: j0 , j1   !< start/end points in y-direction
 integer (ip) , intent (inout) :: k0 , k1   !< start/end points in zo-direction


 select case (domain_id)
 case ( -30 ) ! Dirichlet -z
    i0 = sx             ; i1 = ex
    j0 = sy             ; j1 = ey
    k0 = sz             ; k1 = sz
 case ( -20 ) ! Dirichlet -y
    i0 = sx             ; i1 = ex
    j0 = sy             ; j1 = sy
    k0 = sz             ; k1 = ez
 case ( -10 ) ! Dirichlet -x
    i0 = sx             ; i1 = sx
    j0 = sy             ; j1 = ey
    k0 = sz             ; k1 = ez
 case ( -3 ) ! -z boundaries
    i0 = sx             ; i1 = ex
    j0 = sy             ; j1 = ey
    k0 = sz-ng          ; k1 = sz-1
 case ( -2 ) ! -y boundaries
    i0 = sx             ; i1 = ex
    j0 = sy-ng          ; j1 = sy-1
    k0 = sz             ; k1 = ez
 case ( -1 ) ! -x boundaries
    i0 = sx-ng          ; i1 = sx-1
    j0 = sy             ; j1 = ey
    k0 = sz             ; k1 = ez
 case ( 0 ) ! inner domain
    i0 = sx             ; i1 = ex
    j0 = sy             ; j1 = ey
    k0 = sz             ; k1 = ez
 case ( 1 ) ! +x boundaries
    i0 = ex+1           ; i1 = ex+ng
    j0 = sy             ; j1 = ey
    k0 = sz             ; k1 = ez
 case ( 2 ) ! +y boundaries
    i0 = sx             ; i1 = ex
    j0 = ey+1           ; j1 = ey+ng
    k0 = sz             ; k1 = ez
 case ( 3 ) ! +z boundaries
    i0 = sx             ; i1 = ex
    j0 = sy             ; j1 = ey
    k0 = ez+1           ; k1 = ez+ng
 case ( 10 ) ! Dirichlet +x
    i0 = ex             ; i1 = ex
    j0 = sy             ; j1 = ey
    k0 = sz             ; k1 = ez
 case ( 20 ) ! Dirichlet +y
    i0 = sx             ; i1 = ex
    j0 = ey             ; j1 = ey
    k0 = sz             ; k1 = ez
 case ( 30 ) ! Dirichlet +z
    i0 = sx             ; i1 = ex
    j0 = sy             ; j1 = ey
    k0 = ez             ; k1 = ez
 case ( 6 ) ! inner domain + all boundaries
    ! This is very dangeurous
    write (*,*) 'WARNING: this option is very slow'
    i0 = sx-ng          ; i1 = ex+ng
    j0 = sy-ng          ; j1 = ey+ng
    k0 = sz-ng          ; k1 = ez+ng
 case default
    call abort_mpi ('thd_domain not defined')
 end select


end subroutine domain_select


!> \brief Domain selection for ParaView.
!!
!! Takes into account overlaping (necessary for continous parallel plots).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

subroutine domain_paraview ( i0 , i1 , j0 , j1 , k0 , k1 )


 integer (ip) , intent (inout) :: i0 , i1 !< start/end points in x-direction
 integer (ip) , intent (inout) :: j0 , j1 !< start/end points in y-direction
 integer (ip) , intent (inout) :: k0 , k1 !< start/end points in z-direction


 if ( bc (W) == inner ) then
    i0 = sx-1
 else
    i0 = sx
 end if

 if ( bc (E) == inner ) then
    i1 = ex+1
 else
    i1 = ex
 end if


 if ( bc (S) == inner ) then
    j0 = sy-1
 else
    j0 = sy
 end if

 if ( bc (N) == inner ) then
    j1 = ey+1
 else
    j1 = ey
 end if


 if ( bc (B) == inner ) then
    k0 = sz-1
 else
    k0 = sz
 end if

 if ( bc (F) == inner ) then
    k1 = ez+1
 else
    k1 = ez
 end if


end subroutine domain_paraview


!=================================================================
subroutine print_warning(awarning)
  implicit none

  ! ----------------------
  character(len = *), intent(in) :: awarning
  ! ----------------------
  if ( rank == rank_default )  then
     write(6,'(A)') "Warning :"//trim(awarning)
     call flush(6)
  end if
end subroutine print_warning
!=================================================================


!=================================================================
subroutine error_and_exit(amessage)
  implicit none

  ! ----------------------
  character(len = *), intent(in), optional :: amessage
  ! ----------------------
  character(len=len_default) :: themessage
  ! ----------------------

  ! message
  themessage="Stop"
  if (present(amessage)) then
     themessage = amessage
  end if

  write(6,'(A)') trim(themessage)
  write(6,'(A)') "Exiting..."
  call flush(6)

  call MPI_Abort(MPI_COMM_WORLD,1,mpicode)

end subroutine error_and_exit
!=================================================================

end module parallel
