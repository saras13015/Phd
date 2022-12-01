program vectorsum
#ifdef _OPENACC
  use openacc
#endif
  implicit none
  integer, parameter :: rk = selected_real_kind(12)
  integer, parameter :: ik = selected_int_kind(9)
  integer, parameter :: nx = 102400

  real(kind=rk), dimension(nx) :: vecA, vecB, vecC
  real(kind=rk)    :: sum
  integer(kind=ik) :: i
  real(kind=rk)    :: t1, t2, t3, t4, t5, t6, t7, t8 

  call cpu_time(t1)
  ! Initialization of vectors
  call cpu_time(t5)
  do i = 1, nx
     vecA(i) = 1.0_rk/(real(nx - i + 1, kind=rk))
     vecB(i) = vecA(i)**2
  end do
  call cpu_time(t6)

  ! TODO
  ! Implement vector addition on device with OpenACC 
  ! vecC = vecA + vecB
  call cpu_time(t3)
  do i = 1, nx
        vecC(i) = vecA(i) + vecB(i)
  end do
  call cpu_time(t4)

  ! Compute the check value
  call cpu_time(t7)
  write(*,*) 'Reduction sum: ', sum(vecC)
  call cpu_time(t8)

  call cpu_time(t2)
  write(*,*) 'Time Elapsed ', t2-t1, ' sec'
  write(*,*) 'Time Elapsed in ms ', (t2-t1)*1000, ' ms'
  write(*,*) '----------------------------'
  write(*,*) '         In details'
  write(*,*) '----------------------------'
  write(*,*) 'Vector addition time: ', (t4-t3)*1000, ' ms'
  write(*,*) 'Time rate of the vector addition: ',((t4-t3)/(t2-t1))*100, ' %'
  write(*,*)
  write(*,*) 'Initalization time: ', (t6-t5)*1000, ' ms'
  write(*,*) 'Time rate of the initialization: ',((t6-t5)/(t2-t1))*100, ' %'
  write(*,*)
  write(*,*) 'Sum operator time: ', (t8-t7)*1000, ' ms'
  write(*,*) 'Time rate of the sum: ',((t8-t7)/(t2-t1))*100, ' %'
  
end program vectorsum
