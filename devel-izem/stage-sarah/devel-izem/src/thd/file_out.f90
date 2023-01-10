!-------------------------------------------------------|
!-------------------------------------------------------|
!- To plot one field (mainly with matlab type code)
!-------------------------------------------------------|
!-------------------------------------------------------|
subroutine filegraph_gridonefield(Grid,Field,filename)
 ! modules
 use asphodele
 use commonfile
 use typegrid
 implicit none
 type (grid_type), intent(in)  :: Grid
 character (len=*), intent(in) :: filename
 real(RA), dimension(Grid%Nx(1),Grid%Nx(2),Grid%Nx(3)), intent(in) :: Field
 ! local
 integer(IA)   :: i,j,k,i1,i2,i3,formatfile=20061
 integer(IGO) :: MedData=-2,FinalData=-1,OneField=1
 integer(RGO) :: Int32
 real(RGO)    :: Real32
 i1=1;i2=1;i3=1
 !
 call exist_file(FileName)
 open(unit=unit_data_file,file=filename,form='unformatted')
  Int32=int(formatfile,IGO); write(unit_data_file) Int32
  Int32=int(Grid%Ndim,IGO);  write(unit_data_file) Int32
  if(Grid%Ndim > 0) then
   Int32=int(Grid%Nx(1),IGO);write(unit_data_file) Int32
   do i = 1,Grid%Nx(1)
      Real32 =real(Grid%XP(i,i1),RGO);write(unit_data_file) Real32
   enddo
  endif
  if(Grid%Ndim > 1) then
   i2 = 2
   Int32=int(Grid%Nx(2),IGO);write(unit_data_file) Int32
   do j = 1,Grid%Nx(2)
      Real32 = real(Grid%XP(j,i2),RGO);write(unit_data_file) Real32
   enddo
  endif
  if(Grid%Ndim > 2) then
   i3 = 3
   Int32=int(Grid%Nx(3),IGO);write(unit_data_file) Int32
   do k = 1,Grid%Nx(3)
      Real32 = real(Grid%XP(k,i3),RGO);write(unit_data_file) Real32
   enddo
  endif
  write(unit_data_file) MedData
  write(unit_data_file) OneField
  do i = 1,Grid%Nx(1)
   do j = 1,Grid%Nx(2)
    do k = 1,Grid%Nx(3)
       Real32 = real(Field(i,j,k),RGO);write(unit_data_file) Real32
    enddo
   enddo
  enddo
  write(unit_data_file) FinalData
 close(unit_data_file)
 !
end subroutine filegraph_gridonefield
!-------------------------------------------------------|

!-------------------------------------------------------|
!-------------------------------------------------------|
!- To plotxy1y2 file
!-------------------------------------------------------|
!-------------------------------------------------------|
subroutine filegraph_xyy(nx,x,y1,y2,filename,cflag)
 ! modules
 use asphodele
 use commonfile
 use filetools
 implicit none
 character (len=*), intent(in)           :: filename
 character, intent(in)                   :: cflag
 integer(IA),intent(in)                  :: nx
 real(RA), dimension(nx), intent(in)     :: x,y1,y2
 ! local
 integer(IA)  :: i,j,formatfile=20061
 integer(IGO) :: MedData=-2,FinalData=-1,OneField=1
 integer(RGO) :: Int32
 real(RGO)    :: Real32
 !
 if(cflag=='b') then
  call exist_file(FileName)
   open(unit=unit_data_file,file=filename,form='unformatted')
    Int32=formatfile; write(unit_data_file) Int32
    Int32=nx;         write(unit_data_file) Int32
    do i = 1,nx
       Real32 = real(x(i),RGO); write(unit_data_file) Real32
       Real32 = real(y1(i),RGO);write(unit_data_file) Real32
       Real32 = real(y2(i),RGO);write(unit_data_file) Real32
    enddo
    write(unit_data_file) FinalData
    close(unit_data_file)
    call crline('[x y y] file',filename)
 elseif(cflag=='t') then
   call exist_file(FileName)
   open(unit=unit_data_file,file=filename)
   do i = 1,nx
      write(unit_data_file,5003) x(i),y1(i),y2(i)
   enddo
   close(unit_data_file)
   5003 format(3(E15.6))
   call crline('[x y y] file',filename)
  else
    stop ' Error cflagformat filegraph_xyy in file_out.f90'
  endif

 !
end subroutine filegraph_xyy
!-------------------------------------------------------|
