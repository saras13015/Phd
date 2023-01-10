module filetools
use paramthd
use commonfile
use parallel , only : rank
implicit none
contains

!-------------------------------------------------------|
!-------------------------------------------------------|
!- Mark detection in a data file
!-------------------------------------------------------|
!-------------------------------------------------------|
! iflag = -1 : file not foundth
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
         call stopline('SECURITY : '// filename // ' must contain # END FILE #')
       endif
!        endread = .false.
!        do while (.not.endread)
!           read (unit_data_file,'(a50)') readmark
!           if (readmark == '# END FILE #') then
!              flag    = -2
!              close(unit_data_file)
!              endread = .true.
!              exit
!           endif
!           if (readmark == mark) then
!              read (unit_data_file,'(a50)')
!              flag    = unit_data_file
!              endread =.true.
!              exit
!           endif
!        enddo
    endif
!-
end subroutine file_mark
!-------------------------------------------------------|

!-------------------------------------------------------|
!-------------------------------------------------------|
!- Mark detection in a data file
!-------------------------------------------------------|
!-------------------------------------------------------|
! iflag = -1 : file not found
! iflag = -2 : marker not found
! iflag > 0  : unit number where to read the data
!-------------------------------------------------------|
subroutine file_marks(filename,mark,flag)
 implicit none
 character (len=*), intent(in) :: filename,mark
 integer(IA), intent(out)          :: flag
!-
 character (len=50) :: readmark
 logical            :: endread
 integer(IA)        :: err
!-
   call exist_file(filename)
   open(unit=unit_data_file,file=filename, form='formatted',status='old',iostat&
                                                                           =err)
    if(err /= 0) then
       flag = -1
    else
       endread = .false.
       do while (.not.endread)
          read (unit_data_file,'(a50)') readmark
          if (readmark == '# END FILE #') then
             flag    = -2
             close(unit_data_file)
             endread = .true.
             exit
          endif
          if (readmark == mark) then
             flag    = unit_data_file
             endread =.true.
             exit
          endif
       enddo
    endif
!-
end subroutine file_marks
!-------------------------------------------------------|
!-------------------------------------------------------|

!-------------------------------------------------------|
!     FileName <- Name_Fsub.Number
!-------------------------------------------------------|
subroutine file_tecplot_getname(FileName,SimulName,Fsub,SimulNumber)
   use paramthd
   use commonfile
   implicit none
   integer(IA), intent(in)        :: SimulNumber
   character(len=*), intent(in)   :: SimulName,Fsub
   character(len=*), intent(out)  :: FileName
!-------------------------------------------------------|
   character(len=4)    :: Str
   character(len=3)    :: Str3
   integer(IA)         :: numproc = 0 ! provisoire
!-
   write(Str,'(I4)') SimulNumber
   write(Str,'(A)') repeat('0',4-len(trim(adjustl(Str)))) // trim(adjustl(Str))
   call split_unit(numproc,Str3)
   write(FileName,'(A)') trim(adjustl(Simulname)) // '_' // trim(adjustl(Fsub)) //&
                  '_P' // trim(adjustl(Str3)) // '-' // trim(adjustl(Str)) //'.dat'
!-
end subroutine file_tecplot_getname

subroutine fileana_getname(FileName,SimulName,Fsub)
   use paramthd
   use commonfile
   implicit none
   character(len=*), intent(in)   :: SimulName,Fsub
   character(len=*), intent(out)  :: FileName
!-------------------------------------------------------|
  write(FileName,'(A)') trim(adjustl(Simulname))//'_'//trim(adjustl(Fsub))//'.ana'

end subroutine fileana_getname
!-------------------------------------------------------|
!- Clean output
!-------------------------------------------------------|
subroutine trline(ch)
  !use parallel
  implicit none
  character(len=*), intent(in) :: ch
  if(rank==0) then
   write(iu_out,'(a42)') 'x......................................asp'
   if (len(ch)<25) then
     130 format('x.......',a25,'     .asp') !a32,'     .asp')
     write(iu_out,130) ch !char(27)//'[1m'//ch//char(27)//'[m'
   else
     131 format('x.......',a25,'     .asp')
     write(iu_out,131) ch
   endif
   write(iu_out,'(a42)') 'x......................................asp'
  endif
end subroutine trline

subroutine crline(ch,c)
  !use parallel
   implicit none
   character(len=*), intent(in) :: ch,c
!-
if (rank==0) then
   130 format('x...',a25,' : ',a80)
  write(iu_out,130) ch,c
endif
end subroutine crline

subroutine frline(ch,f)
  !use parallel
   implicit none
   character(len=*), intent(in) :: ch
   real(RA), intent(in) :: f
!-
  if (rank==0) then
   130 format('x...',a25,' : ',e15.6)
   write(iu_out,130) ch,f
  endif
end subroutine frline

subroutine irline(ch,i)
  !use parallel
 implicit none
   character(len=*), intent(in) :: ch
   integer(IA), intent(in) :: i
!-
  if (rank==0) then
   130 format('x...',a25,' : ',i15)
   write(iu_out,130)  ch,i
  endif
end subroutine irline

subroutine wrline(c)
  !use parallel
 implicit none
   character, intent(in) :: c
!-
  if (rank==0) then
   130 format('x...',75A)
   write(iu_out,130) repeat(c,75)
  endif
end subroutine wrline

subroutine srline
  !use parallel
 implicit none
!-
  if (rank==0) then
   130 format('x...')
   write(iu_out,130)
  endif
end subroutine srline

subroutine mrline(ch)
  !use parallel
 implicit none
   character(len=*), intent(in) :: ch
!-
  if (rank==0) then
   130 format('x...     ',a70)
   call wrline('.')
   write(iu_out,130) ch
   call wrline('.')
  endif
end subroutine mrline

subroutine stopline(c)
  !use parallel
 implicit none
   character(len=*), intent(in) :: c
!-
  if (rank==0) then
   130 format('x...',a20,' : ',a80)
   write(iu_out,130) 'stop',c
  endif
   stop

end subroutine stopline

subroutine iterline(i,t,dt)
  !use parallel
 implicit none
   integer(IA), intent(in) :: i
   real(RA), intent(in)    :: t,dt
!-
  !call flush()
  if (rank==0) then
    130 format('x... Ite  :',I8)
    131 format('x... Time :',F15.8,',  Tstep:',F15.8)
    write(iu_out,130) i
    write(iu_out,131) t,dt
  endif
end subroutine iterline

subroutine itercomb(reaction_rate)
  !use parallel
  implicit none
   real(RA), intent(in)    :: reaction_rate
  !call flush()
  if (rank==0) then
    130 format('x... max reaction rate  :',F15.10)
    write(iu_out,130) reaction_rate
  endif
end subroutine itercomb

subroutine clear_screen()
  !use parallel
  implicit none
  if (rank==0) then
    write(iu_out,'(a)',advance='no')char(27)//'c'
    write(iu_out,'(a)',advance='no')char(27)//'[H'
  endif
end subroutine clear_screen

subroutine progress_bar(pc,larg,init,save_pos)
   implicit none
   integer(IA), intent(in) :: init,larg,save_pos
   real(RA), intent(in)    :: pc
   character(len=larg)     :: progress
   character(len=10)        :: str
   if ((pc > 1.0_RA).or.(pc < 0.0_RA)) return
   if (iu_out /= 6) return
!-
   if (init == 1) then
     write(iu_out,'(a)',advance='no')char(27)//'[1B' ! go to end of terminal / &
                                                     ! file
     write(iu_out,*) ' '
     if (save_pos == 1) write(iu_out,'(a)',advance='no')char(27)//'7'
   else
     if (save_pos == 1) write(iu_out,'(a)',advance='no')char(27)//'8'
     write(iu_out,'(a)',advance='no')char(27)//'[1A'
     progress(:) = ' '
     progress(1:int(pc*(larg-8_IA))) = repeat('=',int(pc*(larg-8_IA)))
!      do j=1, int(pc*(larg-8_IA))
!        progress(j:j) = '='
!      enddo
     select case(modulo(int(100*pc*(larg-8_IA)),4))
       case(0)
         progress(larg-7:larg-7) = char(45)
       case(1)
         progress(larg-7:larg-7) = char(92)
       case(2)
         progress(larg-7:larg-7) = char(124)
       case(3)
         progress(larg-7:larg-7) = char(47)
     end select
     write(Str,'(i10)') int(1000.0_RA*pc,IA)
     progress(larg-5:larg-3) = Str(7:9)
     if (progress(larg-3:larg-3) == ' ') progress(larg-3:larg-3) = '0'
     progress(larg-2:larg-2) = '.'
     progress(larg-1:larg-1) = Str(10:10)
     progress(larg:larg) = '%'
     write(iu_out,*) progress
      write(iu_out,'(a)',advance='no')char(27)//'[2B'  ! go to end of      &
                                                           ! terminal / file
!      write(iu_out,'(a)',advance='no')char(27)//'[999999B'  ! go to end of      &
                                                           ! terminal / file
   endif
end subroutine progress_bar

!-------------------------------------------------------|
!-------------------------------------------------------|
subroutine split_unit(I_unit,Cunit) ! from InPetto
  use paramthd
  implicit none
  integer(IA), intent(in)          :: I_unit
  character(len=3),intent(out)     :: Cunit
!-------------------------------------------------------|
  write(Cunit,'(I3)') I_unit
  write(Cunit,'(A)') repeat('0',3-len(trim(adjustl(Cunit)))) // trim(adjustl(  &
                                                                         Cunit))
end subroutine split_unit
!-------------------------------------------------------|
!-------------------------------------------------------|
 subroutine exist_file(filename)
  implicit none
  character (len=*), intent(in)     :: filename
  logical                           :: ex
  inquire(file = fileName, exist=ex)
   if (ex.eqv..FALSE.) then
    write(iu_out,*) char(27)//'[5m'//char(27)//'[1m'//                         &
                        'no such file or directory : ', filename//char(27)//'[m'
!     write(*,*) 'no such file or directory : ', filename
    stop
   endif
 end subroutine exist_file

end module filetools
