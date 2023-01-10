!-------------------------------------------------------|
!     FileName <- Name_Fsub.fix
!-------------------------------------------------------|
subroutine filefix_getname(FileName,SimulName,Fsub)
  use asphodele
  use commonfile
#ifdef MPI
  use mpi_mod
#endif

  implicit none
  character(len=*), intent(in)   :: SimulName,Fsub
  character(len=*), intent(out)  :: FileName
  character(len=4) :: x1,x2,x3

  x1='0000'
  x2='0000'
  x3='0000'

#ifdef MPI
 write(x1,'(I4.4)') Cmp%Cartesian_position_node(1)
 if(Ndim_grid > 1) write(x2,'(I4.4)') Cmp%Cartesian_position_node(2)
 if(Ndim_grid > 2) write(x3,'(I4.4)') Cmp%Cartesian_position_node(3)
 write(FileName,'(A)') trim(adjustl(Simulname))//'_'//trim(adjustl(Fsub)) &
  & //'.'//x1//'.'//x2//'.'//x3//'.fix'
#else
  write(FileName,'(A)') trim(adjustl(Simulname))//'_'//trim(adjustl(Fsub)) &
  & //'_'//x1//'_'//x2//'_'//x3//'.fix'

#endif
end subroutine filefix_getname

! !-------------------------------------------------------|
! !     FileName <- Name_Fsub.Number
! !-------------------------------------------------------|
! subroutine file_getname(FileName,SimulName,Fsub,SimulNumber)
!   use asphodele
!   use commonfile
!   use paramgrid
! #ifdef MPI
!   use mpi_mod
! #endif

!   implicit none
!   integer(IA), intent(in)        :: SimulNumber
!   character(len=*), intent(in)   :: SimulName,Fsub
!   character(len=*), intent(out)  :: FileName
! !-------------------------------------------------------|
!   integer(IA), parameter  :: LStr = 16
!   character(len=LSTr) :: Str
!   character(len=3)    :: Str3
! !  integer(IA)             :: numproc = 0 ! provisoire
!   character(len=4) :: x1,x2,x3

!   x1='0000'
!   x2='0000'
!   x3='0000'

!   write(Str,'(I4)')  SimulNumber
! #ifdef MPI
!   write(x1,'(I4.4)') Cmp%Cartesian_position_node(1)
!   if(Ndim_grid > 1) then
!      write(x2,'(I4.4)') Cmp%Cartesian_position_node(2)
!   end if
!   if(Ndim_grid > 2) then
!      write(x3,'(I4.4)') Cmp%Cartesian_position_node(3)
!   end if
! #endif
!   write(FileName,'(A)') trim(adjustl(Simulname))//'_'//trim(adjustl(Fsub)) &
!   & //'_'//x1//'_'//x2//'_'//x3//'.'// trim(adjustl(Str))
! end subroutine file_getname


! !-------------------------------------------------------|
! !     FileName <- Name_Fsub.Number
! !-------------------------------------------------------|
! subroutine file_getname_real(FileName,SimulName,Fsub,SimulTime)
!   use asphodele
!   use commonfile
!   use paramgrid
! #ifdef MPI
!   use mpi_mod
! #endif

!   implicit none
!   real(RA), intent(in)           :: SimulTime
!   character(len=*), intent(in)   :: SimulName,Fsub
!   character(len=*), intent(out)  :: FileName
! !-------------------------------------------------------|
!   integer(IA), parameter  :: LStr = 16
!   character(len=LSTr) :: Str
!   character(len=3)    :: Str3
! !  integer(IA)             :: numproc = 0 ! provisoire
!   character(len=4) :: x1,x2,x3

!   x1='0000'
!   x2='0000'
!   x3='0000'

! #ifdef MPI
!   write(x1,'(I4.4)') Cmp%Cartesian_position_node(1)
!   if(Ndim_grid > 1) then
!      write(x2,'(I4.4)') Cmp%Cartesian_position_node(2)
!   end if
!   if(Ndim_grid > 2) then
!      write(x3,'(I4.4)') Cmp%Cartesian_position_node(3)
!   end if
! #endif

! !-
!   write(Str,'(F10.5)')  SimulTime
! !  call split_unit(numproc,Str3)
! !  write(FileName,'(A)') trim(adjustl(Simulname))//'_'//trim(adjustl(Fsub))//   &
! !                        '_P' // trim(adjustl(Str3)) // '.' // trim(adjustl(Str))
!   write(FileName,'(A)') trim(adjustl(Simulname))//'_'//trim(adjustl(Fsub)) &
!   & //'_'//x1//'_'//x2//'_'//x3//'_T_'// trim(adjustl(Str)) // '.ana'
! end subroutine file_getname_real
