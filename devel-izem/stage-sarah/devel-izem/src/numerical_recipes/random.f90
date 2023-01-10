module random_num

  use nrtype

  implicit none


contains


!>  Function involving the generation of random numbers.
!!  This Function involves a "long" period (> 2.e+18) random number generator of l'ecuyer with
!!  bays-durham shuffle and add safeguards. It returns a uniform random deviate between 0.0 and 1
!!  (exclusive of the endpoint values). Call with 'idum', a make negative integer to initialize.
!!  Thereafter, do not alter idum between succesive deviates in a sequence. 'rnmx' should approximate
!!  the largest floating value that is less than 1. 
!!
!!  Ref: Numerical recepies - pp 271
!!
!!  @param idum Negative integer used for initialization (IN)
!!
!!  @author J. M. Vedovoto
function rand2 (idum)

    integer(I4B) :: idum

    integer(I4B) :: im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
    real(DP) ::  rand2,am,eps,rnmx
    parameter (im1=2147483563 , im2=2147483399 , am=1.0_dp/im1)
    parameter (imm1=im1-1 , ia1=40014 , ia2=40692)
    parameter (iq1=53668 , iq2=52774 , ir1=12211 , ir2=3791)
    parameter (ntab=64,ndiv=1+imm1/ntab,eps=1.2e-7_dp,rnmx=1.0_dp-eps)
    integer(I4B) :: idum2,j,k,iv(ntab),iy
    save iv,iy,idum2
    data idum2 /123456789/, iv /ntab*0/ , iy /0/


    if (idum.le.0) then    ! initialize.
           idum  = max(-idum,1) ! Be sure to prevent idum = 0.
           idum2 = idum
           do j=ntab+8,1,-1 ! Load the shuffle table (after 8 warm-ups).
                   k    = idum/iq1
                   idum = ia1*(idum - k*iq1)  - k*ir1
                   if (idum.lt.0)   idum = idum + im1
                   if (j.le.ntab)   iv(j) = idum
           enddo
           iy=iv(1)
    endif

    k    = idum/iq1                        ! Start here when not initializing.
    idum = ia1*(idum - k*iq1) - k*ir1      ! Compute idum=mod(IA1*idum,IM1) without over-
    if (idum.lt.0)  idum = idum + im1      !     flows by Schrage’s method.
    k     = idum2/iq2
    idum2 = ia2*(idum2 - k*iq2) - k*ir2    ! Compute idum2=mod(IA2*idum2,IM2) likewise.
    if (idum2.lt.0) idum2 = idum2 + im2
    j     = 1 + iy/ndiv                    ! Will be in the range 1:NTAB.
    iy    = iv(j) - idum2                  ! Here idum is shuffled, idum and idum2 are com-
    iv(j) = idum                           !     bined to generate output.
    if (iy.lt.1) iy = iy + imm1


    rand2 = min(am*iy , rnmx)              ! Because users don’t expect endpoint values.

    return

end function rand2


end module random_num
