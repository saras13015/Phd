program vectorsum
#ifdef _OPENACC
  use openacc
#endif
  implicit none
  integer, parameter :: rk = selected_real_kind(12)
  integer, parameter :: ik = selected_int_kind(9)
  integer, parameter :: nx = 102400
  integer :: ok

  real(kind=rk), dimension(nx) :: vecA, vecB, vecC
  real(kind=rk)    :: sum
  integer(kind=ik) :: i
  real    :: t1, t2, t3, t4, t5, t6, t7, t8 
  real    :: percent_utilization, pu1, pu2, pu3

  call cpu_time(t1)
  ! Initialization of vectors
  call cpu_time(t5)
  ! !$acc kernels
  do i = 1, nx
     vecA(i) = 1.0_rk/(real(nx - i + 1, kind=rk))
     vecB(i) = vecA(i)**2
  end do
  !  !$acc end kernels
  call cpu_time(t6)

  ! TODO
  ! Implement vector addition on device with OpenACC 
  ! vecC = vecA + vecB
  call cpu_time(t3)
!   !$acc data copy(vecA, vecB, vecC)
!  !$acc kernels
  do i = 1, nx
        vecC(i) = vecA(i) + vecB(i)
  end do
!  !$acc end kernels
! !$acc end data
  call cpu_time(t4)

  ! Compute the check value
  call cpu_time(t7)
  write(*,*) 'Reduction sum: ', sum(vecC)
  call cpu_time(t8)

  call cpu_time(t2)
!  write(*,*) 'Time Elapsed ', t2-t1, ' sec'
!  write(*,*) 'Time Elapsed in ms ', (t2-t1)*1000, ' ms'
!  write(*,*) '----------------------------'
!  write(*,*) '         In details'
!  write(*,*) '----------------------------'
!  write(*,*) 'Vector addition time: ', (t4-t3)*1000, ' ms'
!  write(*,*) 'Time rate of the vector addition: ',((t4-t3)/(t2-t1))*100, ' %'
!  write(*,*)
!  write(*,*) 'Initalization time: ', (t6-t5)*1000, ' ms'
!  write(*,*) 'Time rate of the initialization: ',((t6-t5)/(t2-t1))*100, ' %'
!  write(*,*)
!  write(*,*) 'Sum operator time: ', (t8-t7)*1000, ' ms'
!  write(*,*) 'Time rate of the sum: ',((t8-t7)/(t2-t1))*100, ' %'
  !------------ TIME PROFILING USING ROUTINES ------------!
  call time_profiling('Initialization', t6-t5, t2-t1, pu1)
  call time_profiling('Vector addition', t4-t3, t2-t1, pu2)
  call time_profiling('Sum', t8-t7, t2-t1, pu3)

!---------------- Saving data into a file --------------!
open(unit=10, file="serial.dat", iostat=ok)
if (ok/=0) stop
write(10,*) '-----------------------------------------------------------------------------' 
write(10,*) '                 |   Initialization  |  Vector addition |        Sum         |'
write(10,*) '-----------------------------------------------------------------------------'
write(10,"(A19,f12.3,A8,f12.3,A8,f12.3,A8)") 'Time consumed   |',(t6-t5)*1000,'|', &
        (t4-t3)*1000, ' | ', (t8-t7)*1000, ' |' 
write(10,*) '      (ms)'
write(10,*) '-----------------------------------------------------------------------------'
write(10,"(A19,f12.3,A8,f12.3,A8,f12.3,A8)") 'Percentage   |',pu1,'|', &
        pu2, ' | ', pu3, ' |' 
write(10,*) '-----------------------------------------------------------------------------'
close(10)

  contains

          subroutine time_profiling(string, elapsed_time, execution_time, percent_utilization)
                  implicit none
                  character (len=*), intent(in) :: string
                  character (len=len(string)+5) :: temp
                  real, intent(in):: elapsed_time, execution_time
                  real, intent(out) :: percent_utilization

                  write(*,*) '=================================='
                  write(*,*) '      Time consumed by ', string, ' routine'
                  write(*,*)
                  write(*,*) elapsed_time, ' sec (', elapsed_time*1000, ' ms).'
                  write(*,*) 
                  write(*,*) 'Percent of time consumed by this routine'
                  percent_utilization = (elapsed_time/execution_time)*100
                  write(*,*) percent_utilization, ' %'
                  write(*,*) '=================================='
                  write(*,*)

          end subroutine time_profiling
  
end program vectorsum
