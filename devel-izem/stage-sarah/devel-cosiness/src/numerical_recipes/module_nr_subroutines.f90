module nr_subroutines

contains

!========================================================
FUNCTION gammln(xx)
!========================================================
USE nrtype; USE nrutil, ONLY : arth,assert
IMPLICIT NONE
REAL(DP), INTENT(IN) :: xx
REAL(DP) :: gammln
!Returns the value ln[Γ(xx)] for xx > 0.0_dp
REAL(DP) :: tmp,x
!Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
!accuracy is good enough.
REAL(DP) :: stp = 2.5066282746310005_dp
REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
-86.50532032941677_dp,24.01409824083091_dp,&
-1.231739572450155_dp,0.1208650973866179e-2_dp,&
-0.5395239384953e-5_dp/)
call assert(xx > 0.0, 'gammln arg')
x=xx
tmp=x+5.5_dp
tmp=(x+0.5_dp)*log(tmp)-tmp
gammln=tmp+log(stp*(1.000000000190015_dp+&
sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
END FUNCTION gammln


!========================================================
FUNCTION beta_s(z,w)
!========================================================
USE nrtype
!USE nr, ONLY : gammln
IMPLICIT NONE
REAL(DP), INTENT(IN) :: z,w
REAL(DP) :: beta_s
!Returns the value of the beta function B(z, w).
beta_s = 1.0_dp
beta_s=exp(gammln(z)+gammln(w)-gammln(z+w))
END FUNCTION beta_s


!========================================================
FUNCTION betai(a,b,x)
!========================================================
use nrtype
use nrutil ,  ONLY : assert 
implicit none
real(dp) , intent(in) :: a,b,x 
real(dp) :: betai , x1 
real(dp) :: bt 
!Returns the incomplete beta function I_x (a,b).

x1=x
if(x1<0.0_dp) x1=0.0_dp
if(x1>1.0_dp) x1=1.0_dp
call assert(x1 >= 0.0_dp ,  x1 <= 1.0_dp ,   'betai arg' ) 

if (x1 == 0.0_dp .or. x1 == 1.0_dp) then 
  bt=0.0_dp
else  !Factors in front of the continued fraction. 
  bt=exp(gammln(a+b)-gammln(a)-gammln(b)+& 
         a*dlog(x1)+b*dlog(1.0_dp-x1)) 
end if 

if (x1 < (a+1.0_dp)/(a+b+2.0_dp)) then   
  betai=bt*betacf(a,b,x1)/a   !Use continued fraction directly.
else   !Use continued fraction after making the symmetry transformation. 
  betai=1.0_dp-bt*betacf(b,a,1.0_dp-x1)/b 
end if 
END FUNCTION betai


!========================================================
FUNCTION betacf(a,b,x)
!========================================================
use nrtype
use nrutil ,  ONLY : nrerror 

implicit none

real(dp) , intent(in) :: a,b,x 
real(dp)  :: betacf
integer(I4B) , parameter :: MAXIT=500
real(dp)     , parameter :: EPS=epsilon(x), FPMIN=tiny(x)/EPS 
!Used by betai: Evaluates continued fraction 
!for incomplete beta function by modified Lentz s method (§5.2). 

real(dp) :: aa , c , d , del , h , qab , qam , qap 
integer(I4B) :: it , m2 
!write(*,*) 'a,b,x :' , a,b,x
!These q's will be used in factors that occur in the coeffcients (6.4.6). 
  qab=a+b
  qap=a+1.0_dp
  qam=a-1.0_dp 
!first step of Lentz s method.
  c=1.0
  d=1.0_dp-qab*x/qap 

if (abs(d) < FPMIN) d=FPMIN 
d=1.0_dp/d 
h=d 
betacf = 0.0_dp ! Initialization (A.Techer for MIL)

do it=1 , MAXIT 
  m2=2*it 
  aa=it*(b-it)*x/((qam+m2)*(a+m2)) 
  d=1.0_dp+aa*d         !One step (the even one) of the recurrence. 
  if (abs(d) < FPMIN) d=FPMIN 
  c=1.0_dp+aa/c 
  if (abs(c) < FPMIN) c=FPMIN 
  d=1.0_dp/d 
  h=h*d*c 
  aa=-(a+it)*(qab+it)*x/((a+m2)*(qap+m2)) 
  d=1.0_dp+aa*d         !Next step of the recurrence (the odd one). 
  if (abs(d) < FPMIN) d=FPMIN 
  c=1.0_dp+aa/c 
  if (abs(c) < FPMIN) c=FPMIN 
  d=1.0_dp/d 
  del=d*c 
  h=h*del 
  if (abs(del-1.0_dp) <= EPS) exit        !Are we done? 
end do 

if (it > MAXIT) then
!  write (*,'(1X,A,3(1X,1PE10.3),2(1X,I5))') 'a,b,x,it,MAXIT =' , a,b,x,it,MAXIT
!  call nrerror('a or b too big, or MAXIT too small in betacf')
  return ! (A.Techer for MIL)
end if

betacf=h 
!        write(* , *) 'betacf' ,  h
END FUNCTION betacf


!========================================================
SUBROUTINE jacobi(a,d,v,nrot)
!========================================================
  USE nrtype; USE nrutil, ONLY : assert_eq,get_diag,nrerror,unit_matrix,&
       upper_triangle
  IMPLICIT NONE
  INTEGER(I4B), INTENT(OUT) :: nrot
  REAL(SP), DIMENSION(:), INTENT(OUT) :: d
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
  REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
  ! Computes all eigenvalues and eigenvectors of a real symmetric N × N matrix a. On output,
  ! elements of a above the diagonal are destroyed. d is a vector of length N that returns the
  ! eigenvalues of a. v is an N × N matrix whose columns contain, on output, the normalized
  ! eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
  INTEGER(I4B) :: i,ip,iq,n
  REAL(SP) :: c,g,h,s,sm,t,tau,theta,tresh
  REAL(SP), DIMENSION(size(d)) :: b,z
  n=assert_eq((/size(a,1),size(a,2),size(d),size(v,1),size(v,2)/),'jacobi')
  call unit_matrix(v(:,:))
  ! Initialize v to the identity matrix.
  b(:)=get_diag(a(:,:))
  d(:)=b(:)
  z(:)=0.0
  ! This vector will accumulate terms of
  ! the form tapq as in eq. (11.1.14).
  nrot=0
  do i=1,50
     sm=sum(abs(a),mask=upper_triangle(n,n))
     ! Sum off-diagonal elements.
     if (sm == 0.0) RETURN
     ! The normal return, which relies on quadratic convergence to machine underflow.
     tresh=merge(0.2_sp*sm/n**2,0.0_sp, i < 4 )
     ! On the first three sweeps, we will rotate only if tresh exceeded.
     do ip=1,n-1
        do iq=ip+1,n
           g=100.0_sp*abs(a(ip,iq))
           ! After four sweeps, skip the rotation if the off-diagonal element is small.
           if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
                .and. (abs(d(iq))+g == abs(d(iq)))) then
              a(ip,iq)=0.0
           else if (abs(a(ip,iq)) > tresh) then
              h=d(iq)-d(ip)
              if (abs(h)+g == abs(h)) then
                 t=a(ip,iq)/h                 ! t = 1/(2θ)
              else
                 theta=0.5_sp*h/a(ip,iq)      ! Equation (11.1.10).
                 t=1.0_sp/(abs(theta)+sqrt(1.0_sp+theta**2))
                 if (theta < 0.0) t=-t
              end if
              c=1.0_sp/sqrt(1+t**2)
              s=t*c
              tau=s/(1.0_sp+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.0
              call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
              ! Case of rotations 1 ≤ j < p.
              call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
              ! Case of rotations p < j < q.
              call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
              ! Case of rotations q < j ≤ n.
              call jrotate(v(:,ip),v(:,iq))
              nrot=nrot+1
           end if
        end do
     end do
     b(:)=b(:)+z(:)
     d(:)=b(:)
     ! Update d with the sum of tapq ,
     z(:)=0.0
     ! and reinitialize z.
  end do
  call nrerror('too many iterations in jacobi')
CONTAINS
  SUBROUTINE jrotate(a1,a2)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a1,a2
    REAL(SP), DIMENSION(size(a1)) :: wk1
    wk1(:)=a1(:)
    a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
    a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
  END SUBROUTINE jrotate
END SUBROUTINE jacobi

end module
