!------------------------------------------------------------------------------
! MODULE: thermodynamics
!------------------------------------------------------------------------------
!> \brief Thermodynamic relations.
!!
!! This module contains all the subroutines related to the
!! thermophysical mixture.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module thermodynamics

  use parameters
  use parallel
  use type_thd


  implicit none


  integer (ip) , parameter , private        :: itmax_mach    = 100 !< maximum number of iterations (start simulation from restart file)
  integer (ip) , parameter , private        :: itmax         = 200 !< maximum number of iterations (during simulation)
  integer (ip) , parameter , private        :: nptnegT_max   = 40000 !< maximum number of points where T is negative

  real (dp) , parameter , private           :: eps_temp_mach = 1.0e-10_dp !< maximum temperature error (start simulation from restart file)
  real (dp) , parameter , private           :: eps_temp      = 1.0e-5_dp  !< maximum temperature error (during simulation)



contains


!> \brief Convert molar faction to mass fraction.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine X_to_Y ( thd , Ya )


    type (thd_type) , intent (in)                                  :: thd !< thermodynamic derived type
    real (dp) , dimension (:) , intent (inout)                     :: Ya  !< molar (in) to mass(out) fraction


    integer (ip) :: l
    real (dp)    :: wrk


    wrk = 0.0_dp
    do l = 1 , nrv
       wrk = wrk + Ya (l) * thd % Wc (l)
    end do
    wrk = 1.0_dp / wrk


    do l = 1 , nrv
       Ya (l) = Ya (l) * thd % Wc (l) * wrk
    end do


  end subroutine X_to_Y


!> \brief Calculate primitive inviscid variables.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine prim_inv_var ( domain_id , thd , v , W_i , T , cp , ha )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha        !< especies enthalpy


    integer (ip)                                      :: ok , it , i , j , k , l , ipt
    integer (ip)                                      :: ivar
    integer (ip) , dimension (nptnegT_max)            :: ir , jr , kr
    integer (ip)                                      :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                         :: error
    real (dp) , allocatable , dimension (:,:,:)       :: DF
    logical                                           :: correct_T


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
    call Wmix_i ( domain_id , thd , v , W_i )


    allocate ( DF ( i0:i1 , j0:j1 , k0:k1 ) , stat = ok )
    if ( ok > 0 ) call abort_mpi ('error allocate prim_inv_var')

    correct_T = .false.

    it = 0
    ipt = 0
    do

       it = it + 1
       if ( it > itmax ) exit

       call newton ( domain_id , thd , v , T , cp , ha , DF )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                T (i,j,k)  = T (i,j,k) + DF (i,j,k)
                DF (i,j,k) = DF (i,j,k) * DF (i,j,k)

                if ( T (i,j,k) < 0.0_dp ) then

                   write (*,'(1X,A,10(1X,I10))') 'WARNING: negative temperature at (rank,domain_id,i,j,k) =' , &
                                                  rank , domain_id , i , j , k

                   ipt = ipt + 1
                   if ( ipt > nptnegT_max ) call abort_mpi ('WARNING: too many points with negative temperature')

                   ir(ipt) = i ; jr(ipt) = j ; kr(ipt) = k
                   it = 0 ! Reset Newton routine

                   ! Negative temperature correction consists in averaging the conservatives variables
                   ! on the left and right.
                   correct_T = .true.
                   do ivar = 1,nv
                      v (i,j,k,ivar) = 0.5_dp * ( v (i-1,j,k,ivar) + v (i+1,j,k,ivar) )
                   end do

                end if

             end do
          end do
       end do

       error = sqrt ( maxval ( DF (i0:i1,j0:j1,k0:k1) ) )
       if ( error < eps_temp ) exit

       if ( correct_T ) then

          ! Averaging the temperature on the left and right
          do l = 1 , ipt
             T (ir(l),jr(l),kr(l)) = 0.5_dp * ( T (ir(l)-1,jr(l),kr(l)) + T (ir(l)+1,jr(l),kr(l)) )
          end do

          ! Write a report in a file
          open (unit=13,file='report_negative_temp.out',status='unknown')

             write (13,'(1X,A,10(1X,I10))') '(rank,domain_id,i,j,k)=' , &
                                             rank , domain_id , ir(ipt) , jr(ipt) , kr(ipt)

             write (13,'(20(A10))') 'i','rho','rhou','rhov','rhow','rhoet',&
                                    'rhoY1','rhoY2','rhoY3','rhoY4','rhoY5','rhoY6','rhoY7','rhoY8','rhoY9','T'
             do i = i0 , i1
                write (13,'(I5,20(1PE13.3))') i , v(i,jr(ipt),kr(ipt),1:nv) , T (i,jr(ipt),kr(ipt)) * thd % Tinf
             end do

          close (13)

       end if

    end do


    if ( it > itmax ) write (*,*) 'WARNING: failure in temperatures convergence' , error

    if ( error /= error ) call abort_mpi ('error is a NaN')


    deallocate (DF)


  end subroutine prim_inv_var


!> \brief Calculate primitive inviscid variables (machine precision).
!!
!! Use only for initialization and restart simulations. Machine
!! precision allows to increase the accuracy of the solution for the
!! first time step.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine prim_inv_var_err_mach ( domain_id , thd , v , W_i , T , cp , ha )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha        !< especies enthalpy


    integer (ip)                                      :: ok , it , i , j , k
    integer (ip)                                      :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                         :: error
    real (dp) , allocatable , dimension (:,:,:)       :: DF


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
    call Wmix_i ( domain_id , thd , v , W_i )


    allocate ( DF ( i0:i1 , j0:j1 , k0:k1 ) , stat = ok )
    if ( ok > 0 ) stop 'error allocate prim_inv_var_err_mach'


    it = 0
    do

       it = it + 1
       if ( it > itmax_mach ) exit

       call newton ( domain_id , thd , v , T , cp , ha , DF )

       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                T (i,j,k)  = T (i,j,k) + DF (i,j,k)
                DF (i,j,k) = DF (i,j,k) * DF (i,j,k)
             end do
          end do
       end do

       error = sqrt ( maxval ( DF (i0:i1,j0:j1,k0:k1) ) )
       if ( error < eps_temp_mach ) exit

    end do


    if ( minval ( T (i0:i1,j0:j1,k0:k1) ) < 0.0_dp )        &
         write (*,*) 'WARNING: negative temperature at =' , &
         minloc ( T (i0:i1,j0:j1,k0:k1) )

    if ( it > itmax_mach ) write (*,*) 'WARNING: failure in temperatures convergence' , error


    deallocate (DF)


  end subroutine prim_inv_var_err_mach


!> \brief Newton-Ramphson iterative procedure.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine newton ( domain_id , thd , v , T , cp , ha , F )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha        !< especies enthalpy
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: F         !< error function


    integer (ip) , parameter                      :: npoly = 5 ! order of the polynom
    integer (ip) , parameter                      :: io1=1 , io2=2 , io3=3 , io4=4 , io5=5 , io6=6

    integer (ip)                                  :: ok , is
    integer (ip)                                  :: i , j , k
    integer (ip)                                  :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                     :: Tinf_i , onethird , gam2_i , wrk
    real (dp)                                     :: wrk1
    real (dp) , dimension (npoly)                 :: tempkelvin

    real (dp) , allocatable , dimension (:,:,:)   :: rho_i
    real (dp) , allocatable , dimension (:,:,:,:) :: Ya

    real (dp)                                     :: ha0  ,& ! specific enthalpy at Tmin or Tmax
                                                     dha0    ! specific enthalpy derivative at Tmin or Tmax
    real (dp)                                     :: cp0  ,& ! heat capacity at Tmin or Tmax
                                                     dcp0    ! heat capacity derivative at Tmin or Tmax

    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    allocate ( rho_i      ( i0:i1 , j0:j1 , k0:k1 )       , &
               Ya         ( i0:i1 , j0:j1 , k0:k1 , nrv ) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate newton'


    Tinf_i   = 1.0_dp / thd % Tinf
    onethird = 1.0_dp / 3.0_dp
    gam2_i   = 1.0_dp / thd % gam2


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             rho_i (i,j,k) = 1.0_dp / v ( i,j,k, 1 )
          end do
       end do
    end do


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             F ( i,j,k ) = - v ( i,j,k, 5 ) * rho_i (i,j,k) + 0.5_dp * &
                           ( v ( i,j,k, 2 ) * v ( i,j,k, 2 )         + &
                             v ( i,j,k, 3 ) * v ( i,j,k, 3 )         + &
                             v ( i,j,k, 4 ) * v ( i,j,k, 4 ) )       * &
                           ( rho_i (i,j,k) * rho_i (i,j,k) )
          end do
       end do
    end do


    do is = 1 , nrv
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                Ya (i,j,k,is) = v (i,j,k,niv+is) * rho_i (i,j,k)
             end do
          end do
       end do
    end do


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1

             cp (i,j,k) = 0.0_dp

             do is = 1 , nrv

                tempkelvin (io1) = T (i,j,k) * thd % Tinf
                tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1)
                tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1)
                tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1)
                tempkelvin (io5) = tempkelvin (io4) * tempkelvin (io1)

                wrk1 = - T (i,j,k)

                ha (i,j,k,is) = 0.0_dp

                if ( tempkelvin (io1) < thd % Tpolymin (is) ) then ! T < Tmin linear fit from Tmin

                   tempkelvin (io1) = thd % Tpolymin (is)                 ! T
                   tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
                   tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
                   tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4
                   tempkelvin (io5) = tempkelvin (io4) * tempkelvin (io1) ! T^5

                   ha0=    ( thd % HeatCap (io1,1,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io2,1,is) * tempkelvin (io2) * 0.50_dp  + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io3) * onethird + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io4) * 0.25_dp  + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io5) * 0.20_dp  + &
                             thd % HeatCap (io6,1,is) )

                   dha0=   ( thd % HeatCap (io1,1,is)                               + &
                             thd % HeatCap (io2,1,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io2)            + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io3)            + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io4)            )  

                   ha (i,j,k , is ) = ha (i,j,k , is )  + thd % Wc_i (is) * Tinf_i  * &
                           ( dha0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymin (is) ) + ha0 )

                   cp0 =   dha0

                   dcp0 =  ( thd % HeatCap (io2,1,is)                                     + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io1) * 2.0_dp         + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io2) * 3.0_dp         + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io3) * 4.0_dp )

                   cp (i,j,k) = cp (i,j,k)  + Ya (i,j,k,is) * thd% gam2 * thd % Wc_i (is) * &
                           ( dcp0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymin (is) ) + cp0 )

                   F ( i,j,k ) = F ( i,j,k ) + Ya (i,j,k,is) * &
                               ( ha ( i,j,k, is ) + wrk1 * thd % Wc_i (is) )

                else if ( tempkelvin (io1) < thd % Tpolylim (is) ) then ! T < Tlim

                   ha (i,j,k , is ) = ha (i,j,k , is )  + thd % Wc_i (is) * Tinf_i  * &
                           ( thd % HeatCap (io1,1,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io2,1,is) * tempkelvin (io2) * 0.50_dp  + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io3) * onethird + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io4) * 0.25_dp  + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io5) * 0.20_dp  + &
                             thd % HeatCap (io6,1,is) )

                   cp (i,j,k) = cp (i,j,k)  + Ya (i,j,k,is) * thd% gam2 * thd % Wc_i (is) * &
                           ( thd % HeatCap (io1,1,is)                                     + &
                             thd % HeatCap (io2,1,is) * tempkelvin (io1)                  + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io2)                  + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io3)                  + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io4) )

                   F ( i,j,k ) = F ( i,j,k ) + Ya (i,j,k,is) * &
                               ( ha ( i,j,k, is ) + wrk1 * thd % Wc_i (is) )

                else if ( tempkelvin (io1) < thd % Tpolymax (is) ) then ! T < Tmax

                   ha (i,j,k , is ) = ha (i,j,k , is )  + thd % Wc_i (is) * Tinf_i  * &
                           ( thd % HeatCap (io1,2,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io2,2,is) * tempkelvin (io2) * 0.50_dp  + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io3) * onethird + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io4) * 0.25_dp  + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io5) * 0.20_dp  + &
                             thd % HeatCap (io6,2,is) )

                   cp (i,j,k) = cp (i,j,k)  + Ya (i,j,k,is) * thd% gam2 * thd % Wc_i (is) * &
                           ( thd % HeatCap (io1,2,is)                                     + &
                             thd % HeatCap (io2,2,is) * tempkelvin (io1)                  + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io2)                  + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io3)                  + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io4) )

                   F ( i,j,k ) = F ( i,j,k ) + Ya (i,j,k,is) * &
                               ( ha ( i,j,k, is ) + wrk1 * thd % Wc_i (is) )

                else ! T > Tmax linear fit from Tmax

                   tempkelvin (io1) = thd % Tpolymax (is)                 ! T
                   tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
                   tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
                   tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4
                   tempkelvin (io5) = tempkelvin (io4) * tempkelvin (io1) ! T^5

                   ha0=    ( thd % HeatCap (io1,2,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io2,2,is) * tempkelvin (io2) * 0.50_dp  + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io3) * onethird + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io4) * 0.25_dp  + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io5) * 0.20_dp  + &
                             thd % HeatCap (io6,2,is) )

                   dha0=   ( thd % HeatCap (io1,2,is)                               + &
                             thd % HeatCap (io2,2,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io2)            + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io3)            + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io4)            )  

                   ha (i,j,k , is ) = ha (i,j,k , is )  + thd % Wc_i (is) * Tinf_i  * &
                           ( dha0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymax (is) ) + ha0 )

                   cp0 =   dha0

                   dcp0 =  ( thd % HeatCap (io2,2,is)                                     + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io1) * 2.0_dp         + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io2) * 3.0_dp         + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io3) * 4.0_dp )

                   cp (i,j,k) = cp (i,j,k)  + Ya (i,j,k,is) * thd% gam2 * thd % Wc_i (is) * &
                           ( dcp0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymax (is) ) + cp0 )

                   F ( i,j,k ) = F ( i,j,k ) + Ya (i,j,k,is) * &
                               ( ha ( i,j,k, is ) + wrk1 * thd % Wc_i (is) )

                end if

             end do

          end do
       end do
    end do


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1

             wrk = 0.0_dp
             do is = 1 , nrv
                wrk = wrk + Ya (i,j,k,is) * thd % Wc_i (is)
             end do

             F (i,j,k) = F (i,j,k) / ( wrk - cp (i,j,k) * gam2_i ) ! DF <- F

          end do
       end do
    end do


    deallocate ( rho_i , Ya )


  end subroutine newton


!> \brief Calculate primitive inviscid variables (except temperature).
!!
!! It provides optimization to the code. 
!! For reactive cases : recalculate W_i(Ya) , cp(T) and ha(T)
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine upd_prim_var_woT ( thd , domain_id , v , T , W_i , cp , ha )


    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: ha        !< especies enthalpy


    integer (ip) , parameter                      :: npoly = 4 ! order of the polynom
    integer (ip) , parameter                      :: io1=1 , io2=2 , io3=3 , io4=4 , io5=5 , io6=6

    integer (ip)                                  :: is , ns
    integer (ip)                                  :: i , j , k
    integer (ip)                                  :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)                                     :: rho_i , Tinf_i , onethird , Ya , gam2_i
    real (dp) , dimension (npoly+1)               :: tempkelvin

    real (dp)                                     :: ha0  ,& ! specific enthalpy at Tmin or Tmax
                                                     dha0    ! specific enthalpy derivative at Tmin or Tmax
    real (dp)                                     :: cp0  ,& ! heat capacity at Tmin or Tmax
                                                     dcp0    ! heat capacity derivative at Tmin or Tmax


    ns       = nrv
    Tinf_i   = 1.0_dp / thd % Tinf
    onethird = 1.0_dp / 3.0_dp
    gam2_i   = 1.0_dp / thd % gam2


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )
    call Wmix_i ( domain_id , thd , v , W_i )


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1

             cp (i,j,k) = 0.0_dp
             rho_i = 1.0_dp / v ( i,j,k, 1 )

             do is = 1 , ns

                tempkelvin (io1) = T (i,j,k) * thd % Tinf
                tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1)
                tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1)
                tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1)
                tempkelvin (io5) = tempkelvin (io4) * tempkelvin (io1)

                Ya = v ( i,j,k, is+niv ) * rho_i
                ha (i,j,k,is) = 0.0_dp

                if ( tempkelvin (io1) < thd % Tpolymin (is) ) then ! T < Tmin linear fit of cp from Tmin

                   tempkelvin (io1) = thd % Tpolymin (is)                 ! T
                   tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
                   tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
                   tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4
                   tempkelvin (io5) = tempkelvin (io4) * tempkelvin (io1) ! T^5

                   ha0=    ( thd % HeatCap (io1,1,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io2,1,is) * tempkelvin (io2) * 0.50_dp  + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io3) * onethird + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io4) * 0.25_dp  + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io5) * 0.20_dp  + &
                             thd % HeatCap (io6,1,is) )

                   dha0=   ( thd % HeatCap (io1,1,is)                               + &
                             thd % HeatCap (io2,1,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io2)            + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io3)            + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io4)            )  

                   ha (i,j,k , is ) = ha (i,j,k , is )  + thd % Wc_i (is) * Tinf_i  * &
                           ( dha0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymin (is) ) + ha0 )

                   cp0 =   dha0

                   dcp0 =  ( thd % HeatCap (io2,1,is)                                     + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io1) * 2.0_dp         + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io2) * 3.0_dp         + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io3) * 4.0_dp )

                   cp (i,j,k) = cp (i,j,k)  +            Ya * thd% gam2 * thd % Wc_i (is) * &
                           ( dcp0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymin (is) ) + cp0 )

                else if ( tempkelvin (io1) < thd % Tpolylim (is) ) then ! T < Tlim

                   ha (i,j,k , is ) = ha (i,j,k , is )  + thd % Wc_i (is) * Tinf_i  * &
                           ( thd % HeatCap (io1,1,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io2,1,is) * tempkelvin (io2) * 0.50_dp  + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io3) * onethird + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io4) * 0.25_dp  + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io5) * 0.20_dp  + &
                             thd % HeatCap (io6,1,is) )

                   cp (i,j,k) = cp (i,j,k)  +            Ya * thd% gam2 * thd % Wc_i (is) * &
                           ( thd % HeatCap (io1,1,is)                                     + &
                             thd % HeatCap (io2,1,is) * tempkelvin (io1)                  + &
                             thd % HeatCap (io3,1,is) * tempkelvin (io2)                  + &
                             thd % HeatCap (io4,1,is) * tempkelvin (io3)                  + &
                             thd % HeatCap (io5,1,is) * tempkelvin (io4) )

                else if ( tempkelvin (io1) < thd % Tpolymax (is) ) then ! T < Tmax

                   ha (i,j,k , is ) = ha (i,j,k , is )  + thd % Wc_i (is) * Tinf_i  * &
                           ( thd % HeatCap (io1,2,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io2,2,is) * tempkelvin (io2) * 0.50_dp  + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io3) * onethird + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io4) * 0.25_dp  + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io5) * 0.20_dp  + &
                             thd % HeatCap (io6,2,is) )

                   cp (i,j,k) = cp (i,j,k)  +            Ya * thd% gam2 * thd % Wc_i (is) * &
                           ( thd % HeatCap (io1,2,is)                                     + &
                             thd % HeatCap (io2,2,is) * tempkelvin (io1)                  + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io2)                  + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io3)                  + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io4) )

                else ! T > Tmax linear fit of h from Tmax

                   tempkelvin (io1) = thd % Tpolymax (is)                 ! T
                   tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
                   tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
                   tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4
                   tempkelvin (io5) = tempkelvin (io4) * tempkelvin (io1) ! T^5

                   ha0=    ( thd % HeatCap (io1,2,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io2,2,is) * tempkelvin (io2) * 0.50_dp  + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io3) * onethird + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io4) * 0.25_dp  + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io5) * 0.20_dp  + &
                             thd % HeatCap (io6,2,is) )

                   dha0=   ( thd % HeatCap (io1,2,is)                               + &
                             thd % HeatCap (io2,2,is) * tempkelvin (io1)            + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io2)            + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io3)            + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io4)            )  

                   ha (i,j,k , is ) = ha (i,j,k , is )  + thd % Wc_i (is) * Tinf_i  * &
                           ( dha0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymax (is) ) + ha0 )

                   cp0 =   dha0

                   dcp0 =  ( thd % HeatCap (io2,2,is)                                     + &
                             thd % HeatCap (io3,2,is) * tempkelvin (io1) * 2.0_dp         + &
                             thd % HeatCap (io4,2,is) * tempkelvin (io2) * 3.0_dp         + &
                             thd % HeatCap (io5,2,is) * tempkelvin (io3) * 4.0_dp )

                   cp (i,j,k) = cp (i,j,k)  +            Ya * thd% gam2 * thd % Wc_i (is) * &
                           ( dcp0 * ( T (i,j,k) * thd % Tinf - thd % Tpolymax (is) ) + cp0 )

                end if

             end do


          end do
       end do
    end do


  end subroutine upd_prim_var_woT


!> \brief Calculate inverted molar mass.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine Wmix_i ( domain_id , thd , v , W_i )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: W_i       !< inverted molar mass


    integer (ip) :: i , j , k , l
    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    W_i (i0:i1,j0:j1,k0:k1) = 0.0_dp
    do l = 1 , nrv
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                W_i (i,j,k) = W_i (i,j,k) + ( v (i,j,k, niv+l ) / v (i,j,k, 1 ) ) * thd % Wc_i (l)
             end do
          end do
       end do
    end do


  end subroutine Wmix_i


!> \brief Calculate inverted molar mass (scalar version).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine Wmix_i_scalar ( thd , Ya , W_i )


    type (thd_type) , intent (in)                                  :: thd !< thermodynamic derived type
    real (dp) , dimension (:) , intent (in)                        :: Ya  !< species array
    real (dp) , intent (inout)                                     :: W_i !< inverted molar mass scalar


    integer (ip) :: l


    W_i = 0.0_dp
    do l = 1 , nrv
       W_i = W_i + Ya (l) * thd % Wc_i (l)
    end do


  end subroutine Wmix_i_scalar


!> \brief Calculate heat capacity (scalar version).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine cp_scalar ( thd , temp , Y , cp )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , intent (in)                                        :: temp !< temperature
    real (dp) , dimension (:) , intent (in)                        :: Y    !< species array
    real (dp) , intent (inout)                                     :: cp   !< heat capacity scalar


    integer (ip) , parameter                                       :: npoly = 4 ! order of the polynom
    integer (ip) , parameter                                       :: io1=1 , io2=2 , io3=3 , io4=4 , io5=5

    integer (ip)                                                   :: is , ns
    real (dp) , dimension (npoly)                                  :: tempkelvin

    real (dp)                      :: cp0  ,& ! heat capacity at Tmin or Tmax
                                      dcp0    ! heat capacity derivative at Tmin or Tmax


    ns = thd % Nspc

    cp = 0.0_dp

    do is = 1 , ns

    tempkelvin (io1) = temp             * thd % Tinf       ! T
    tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
    tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
    tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4

       if ( tempkelvin (io1) < thd % Tpolymin (is) ) then ! T < Tmin linear fit of cp from Tmin

          tempkelvin (io1) = thd % Tpolymin (is)                 ! T
          tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
          tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
          tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4

          cp0 =   ( thd % HeatCap (io1,1,is)                                     + &
                    thd % HeatCap (io2,1,is) * tempkelvin (io1)                  + &
                    thd % HeatCap (io3,1,is) * tempkelvin (io2)                  + &
                    thd % HeatCap (io4,1,is) * tempkelvin (io3)                  + &
                    thd % HeatCap (io5,1,is) * tempkelvin (io4) )

          dcp0 =  ( thd % HeatCap (io2,1,is)                                     + &
                    thd % HeatCap (io3,1,is) * tempkelvin (io1) * 2.0_dp         + &
                    thd % HeatCap (io4,1,is) * tempkelvin (io2) * 3.0_dp         + &
                    thd % HeatCap (io5,1,is) * tempkelvin (io3) * 4.0_dp )

          cp = cp                         + Y (is) * thd% gam2 * thd % Wc_i (is) * &
                  ( dcp0 * ( temp * thd % Tinf - thd % Tpolymin (is) ) + cp0 )

       else if ( tempkelvin (io1) < thd % Tpolylim (is) ) then ! T < Tlim

          cp = cp                         + Y (is) * thd% gam2 * thd % Wc_i (is) * &
                  ( thd % HeatCap (io1,1,is)                                     + &
                    thd % HeatCap (io2,1,is) * tempkelvin (io1)                  + &
                    thd % HeatCap (io3,1,is) * tempkelvin (io2)                  + &
                    thd % HeatCap (io4,1,is) * tempkelvin (io3)                  + &
                    thd % HeatCap (io5,1,is) * tempkelvin (io4) )

       else if ( tempkelvin (io1) < thd % Tpolymax (is) ) then ! T < Tmax

          cp = cp                         + Y (is) * thd% gam2 * thd % Wc_i (is) * &
                  ( thd % HeatCap (io1,2,is)                                     + &
                    thd % HeatCap (io2,2,is) * tempkelvin (io1)                  + &
                    thd % HeatCap (io3,2,is) * tempkelvin (io2)                  + &
                    thd % HeatCap (io4,2,is) * tempkelvin (io3)                  + &
                    thd % HeatCap (io5,2,is) * tempkelvin (io4) )

       else ! T > Tmax linear fit of cp from Tmax

          tempkelvin (io1) = thd % Tpolymax (is)                 ! T
          tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
          tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
          tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4

          cp0 =   ( thd % HeatCap (io1,2,is)                                     + &
                    thd % HeatCap (io2,2,is) * tempkelvin (io1)                  + &
                    thd % HeatCap (io3,2,is) * tempkelvin (io2)                  + &
                    thd % HeatCap (io4,2,is) * tempkelvin (io3)                  + &
                    thd % HeatCap (io5,2,is) * tempkelvin (io4) )

          dcp0 =  ( thd % HeatCap (io2,2,is)                                     + &
                    thd % HeatCap (io3,2,is) * tempkelvin (io1) * 2.0_dp         + &
                    thd % HeatCap (io4,2,is) * tempkelvin (io2) * 3.0_dp         + &
                    thd % HeatCap (io5,2,is) * tempkelvin (io3) * 4.0_dp )

          cp = cp                         + Y (is) * thd% gam2 * thd % Wc_i (is) * &
                  ( dcp0 * ( temp * thd % Tinf - thd % Tpolymax (is) ) + cp0 )

       end if

    end do


  end subroutine cp_scalar


!> \brief Calculate species enthalpy (scalar version).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine ha_scalar ( thd , temp , ha )


    type (thd_type) , intent (in)                                  :: thd  !< thermodynamic derived type
    real (dp) , intent (in)                                        :: temp !< temperature
    real (dp) , dimension (:) , intent (inout)                     :: ha   !< species enthalpy scalar


    integer (ip) , parameter                      :: npoly = 4 ! order of the polynom
    integer (ip) , parameter                      :: io1=1 , io2=2 , io3=3 , io4=4 , io5=5 , io6=6

    integer (ip)                                  :: is , ns
    real (dp)                                     :: Tinf_i , onethird
    real (dp) , dimension (npoly+1)               :: tempkelvin

    real (dp)                      :: ha0 ,& ! specific enthalpy at Tmin or Tmax
                                      dha0   ! specific enthalpy derivative at Tmin or Tmax

    ns       = thd % Nspc
    Tinf_i   = 1.0_dp / thd % Tinf
    onethird = 1.0_dp / 3.0_dp

    ha = 0.0_dp

    do is = 1 , ns

    tempkelvin (io1) = temp             * thd % Tinf       ! T
    tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
    tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
    tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4
    tempkelvin (io5) = tempkelvin (io4) * tempkelvin (io1) ! T^5

       if ( tempkelvin (io1) < thd % Tpolymin (is) ) then ! T < Tmin linear fit of h from Tmin

          tempkelvin (io1) = thd % Tpolymin (is)                 ! T
          tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
          tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
          tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4
          tempkelvin (io5) = tempkelvin (io4) * tempkelvin (io1) ! T^5

          ha0=    ( thd % HeatCap (io1,1,is) * tempkelvin (io1)            + &
                    thd % HeatCap (io2,1,is) * tempkelvin (io2) * 0.50_dp  + &
                    thd % HeatCap (io3,1,is) * tempkelvin (io3) * onethird + &
                    thd % HeatCap (io4,1,is) * tempkelvin (io4) * 0.25_dp  + &
                    thd % HeatCap (io5,1,is) * tempkelvin (io5) * 0.20_dp  + &
                    thd % HeatCap (io6,1,is) )

          dha0=   ( thd % HeatCap (io1,1,is)                               + &
                    thd % HeatCap (io2,1,is) * tempkelvin (io1)            + &
                    thd % HeatCap (io3,1,is) * tempkelvin (io2)            + &
                    thd % HeatCap (io4,1,is) * tempkelvin (io3)            + &
                    thd % HeatCap (io5,1,is) * tempkelvin (io4)            )  

          ha (is) = ha (is)                    + thd % Wc_i (is) * Tinf_i  * &
                  ( dha0 * ( temp * thd % Tinf - thd % Tpolymin (is) ) + ha0 )

       else if ( tempkelvin (io1) < thd % Tpolylim (is) ) then ! T < Tlim

          ha (is) = ha (is)                    + thd % Wc_i (is) * Tinf_i  * &
                  ( thd % HeatCap (io1,1,is) * tempkelvin (io1)            + &
                    thd % HeatCap (io2,1,is) * tempkelvin (io2) * 0.50_dp  + &
                    thd % HeatCap (io3,1,is) * tempkelvin (io3) * onethird + &
                    thd % HeatCap (io4,1,is) * tempkelvin (io4) * 0.25_dp  + &
                    thd % HeatCap (io5,1,is) * tempkelvin (io5) * 0.20_dp  + &
                    thd % HeatCap (io6,1,is) )

       else if ( tempkelvin (io1) < thd % Tpolymax (is) ) then ! T < Tmax

          ha (is) = ha (is)                    + thd % Wc_i (is) * Tinf_i  * &
                  ( thd % HeatCap (io1,2,is) * tempkelvin (io1)            + &
                    thd % HeatCap (io2,2,is) * tempkelvin (io2) * 0.50_dp  + &
                    thd % HeatCap (io3,2,is) * tempkelvin (io3) * onethird + &
                    thd % HeatCap (io4,2,is) * tempkelvin (io4) * 0.25_dp  + &
                    thd % HeatCap (io5,2,is) * tempkelvin (io5) * 0.20_dp  + &
                    thd % HeatCap (io6,2,is) )

       else ! T > Tmax linear fit of h from Tmax

          tempkelvin (io1) = thd % Tpolymax (is)                 ! T
          tempkelvin (io2) = tempkelvin (io1) * tempkelvin (io1) ! T^2
          tempkelvin (io3) = tempkelvin (io2) * tempkelvin (io1) ! T^3
          tempkelvin (io4) = tempkelvin (io3) * tempkelvin (io1) ! T^4
          tempkelvin (io5) = tempkelvin (io4) * tempkelvin (io1) ! T^5

          ha0=    ( thd % HeatCap (io1,2,is) * tempkelvin (io1)            + &
                    thd % HeatCap (io2,2,is) * tempkelvin (io2) * 0.50_dp  + &
                    thd % HeatCap (io3,2,is) * tempkelvin (io3) * onethird + &
                    thd % HeatCap (io4,2,is) * tempkelvin (io4) * 0.25_dp  + &
                    thd % HeatCap (io5,2,is) * tempkelvin (io5) * 0.20_dp  + &
                    thd % HeatCap (io6,2,is) )

          dha0=   ( thd % HeatCap (io1,2,is)                               + &
                    thd % HeatCap (io2,2,is) * tempkelvin (io1)            + &
                    thd % HeatCap (io3,2,is) * tempkelvin (io2)            + &
                    thd % HeatCap (io4,2,is) * tempkelvin (io3)            + &
                    thd % HeatCap (io5,2,is) * tempkelvin (io4)            )  

          ha (is) = ha (is)                    + thd % Wc_i (is) * Tinf_i  * &
                  ( dha0 * ( temp * thd % Tinf - thd % Tpolymax (is) ) + ha0 )

       end if

    end do


  end subroutine ha_scalar


!> \brief Calculate primitive viscous variables.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine prim_vis_vars ( domain_id , thd , v , W_i , T , Xa , dm , mu , ct )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: Xa        !< molar mass array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: dm        !< diffusion matrix
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mu        !< dynamic viscosity
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: ct        !< thermal conductivity


    integer (ip) , parameter :: npoly = 4 ! order of the polynom
    integer (ip) , parameter :: io1=1 , io2=2 , io3=3 , io4=4

    integer (ip) :: ok
    integer (ip) :: n , ns
    integer (ip) :: i , j , k , l , m
    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1


    real (dp) :: P , Kl , Ku , c1 , c2 , c3 , c4 , Xa_l , Xa_m
    real (dp) , dimension (nrv)                    :: Kd

    real (dp) , allocatable , dimension (:,:,:)     :: wrk1
    real (dp) , allocatable , dimension (:,:,:,:)   :: logtempkelvin
    real (dp) , allocatable , dimension (:,:,:,:)   :: wrk2
    real (dp) , allocatable , dimension (:,:,:,:,:) :: DiffKJ


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    ns        = nrv
    Kl        = 0.5_dp * 1.0e-5_dp / thd % laminf
    Ku        = 1.0e-1_dp / thd % muinf
    Kd (1:ns) = 1.0e-4_dp * thd % p0 / ( thd % pinf * thd % Dinf (1:ns) )


    allocate ( DiffKJ        ( i0:i1 , j0:j1 , k0:k1 , 1:ns , 1:ns ) , &
               logtempkelvin ( 1:npoly , i0:i1 , j0:j1 , k0:k1     ) , &
               wrk1          ( i0:i1 , j0:j1 , k0:k1               ) , &
               wrk2          ( i0:i1 , j0:j1 , k0:k1 , 1:ns        ) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate prim_vis_vars'


    logtempkelvin ( io1 , i0:i1,j0:j1,k0:k1) = 1.0_dp                                        ! (log T)^0
    logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1) = log ( T (i0:i1,j0:j1,k0:k1) * thd % Tinf )    ! (log T)^1
    logtempkelvin ( io3 , i0:i1,j0:j1,k0:k1) = logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1 ) * &
                                               logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1 )     ! (log T)^2
    logtempkelvin ( io4 , i0:i1,j0:j1,k0:k1) = logtempkelvin ( io3 , i0:i1,j0:j1,k0:k1 ) * &
                                               logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1 )     ! (log T)^3


    ! Diffusion matrix


    DiffKJ (i0:i1,j0:j1,k0:k1,1:ns,1:ns) = 0.0_dp
    do m = 1 , ns
       do l = 1 , ns
          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   do n = 1 , npoly
                      DiffKJ (i,j,k,l,m) = DiffKJ (i,j,k,l,m)      + &
                                           thd % CofD (n,l,m)      * &
                                           logtempkelvin (n,i,j,k)
                   end do
                end do
             end do
          end do
       end do
    end do


    dm (i0:i1,j0:j1,k0:k1,1:ns) = 0.0_dp
    do m = 1 , ns

       wrk1 (i0:i1,j0:j1,k0:k1) = 0.0_dp

       do l = 1 , ns

          if ( m /= l ) then

             wrk1 (i0:i1,j0:j1,k0:k1) = wrk1 (i0:i1,j0:j1,k0:k1)                + &
                                      ( Xa (i0:i1,j0:j1,k0:k1,l) + epsm12_asp ) * &
                                        thd % Wc (l) * W_i (i0:i1,j0:j1,k0:k1)

             dm (i0:i1,j0:j1,k0:k1,m) = dm (i0:i1,j0:j1,k0:k1,m)                         + &
                                      ( Xa (i0:i1,j0:j1,k0:k1,l) + epsm12_asp )          * &
                                      (  1.0_dp / exp ( DiffKJ (i0:i1,j0:j1,k0:k1,m,l) ) )

          end if

       end do


       ! Kd(m)
       ! CGS units: cm2/s
       ! conversion into SI units: m2/s
       ! pressure adjustement in bar (1)


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                P            = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
                dm (i,j,k,m) = wrk1 (i,j,k) * Kd (m) / ( dm (i,j,k,m) * P )
             end do
          end do
       end do


    end do


    deallocate (DiffKJ)


    ! Conductivity


    wrk2 (i0:i1,j0:j1,k0:k1,1:ns) = 0.0_dp
    do l = 1 , ns
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                do m = 1 , npoly
                   wrk2 (i,j,k,l) = wrk2 (i,j,k,l)          + &
                                    thd % CofLam (m,l)      * &
                                    logtempkelvin (m,i,j,k)
                end do
             end do
          end do
       end do
    end do


    wrk1 (i0:i1,j0:j1,k0:k1) = 0.0_dp
    ct   (i0:i1,j0:j1,k0:k1) = 0.0_dp
    do l = 1 , ns
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                c1           = exp ( wrk2 (i,j,k,l) )
                c2           = Xa (i,j,k,l)
                ct (i,j,k)   = ct (i,j,k)   + c2 / c1
                wrk1 (i,j,k) = wrk1 (i,j,k) + c2 * c1
             end do
          end do
       end do
    end do


    ! Kl
    ! CGS units: erg/(cm*K*s)
    ! conversion into SI units: J/(m*K*s)


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             ct (i,j,k) = Kl * ( 1.0_dp / ct (i,j,k) + wrk1 (i,j,k) )
          end do
       end do
    end do


    deallocate (wrk1)


    ! Viscosity


    do l = 1 , ns
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                wrk2 (i,j,k,l) = 0.0_dp

                do m = 1 , npoly
                   wrk2 (i,j,k,l) = wrk2 (i,j,k,l)          + &
                                    thd % CofEta (m,l)      * &
                                    logtempkelvin (m,i,j,k)
                end do

                wrk2 (i,j,k,l) = exp ( wrk2 (i,j,k,l) )

             end do
          end do
       end do
    end do


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1

             mu (i,j,k) = 0.0_dp

             do m = 1 , ns

                Xa_m = Xa (i,j,k,m)
                c4 = 0.0_dp

                do l = 1 , ns

                   Xa_l = Xa (i,j,k,l)

                      c1 = ( thd % W (l) / thd % W (m) ) ** 0.25_dp
                      c2 = 1.0_dp / sqrt ( ( 1.0_dp + thd % W (m) / thd % W (l) ) * 8.0_dp )
                      c3 = 1.0_dp + c1 * sqrt ( wrk2 (i,j,k,m) / wrk2 (i,j,k,l) )

                      c2 = c3 * c3 * c2

                      c4 = c4 + c2 * Xa_l

                end do


                ! Ku
                ! CGS units: g/(cm/s)
                ! conversion into SI units: kg/(m*s)
                mu (i,j,k) = mu (i,j,k) + Ku * ( wrk2 (i,j,k,m) * Xa_m / c4 )


             end do


          end do
       end do
    end do


    deallocate ( logtempkelvin , wrk2 )


  end subroutine prim_vis_vars


!> \brief Calculate primitive inviscid variables (except viscosity).
!!
!! It provides optimization to the code.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine prim_vis_vars_wo_mu ( domain_id , thd , v , W_i , T , Xa , dm , ct )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: Xa        !< molar mass array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: dm        !< diffusion matrix
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: ct        !< thermal conductivity


    integer (ip) , parameter :: npoly = 4 ! order of the polynom
    integer (ip) , parameter :: io1=1 , io2=2 , io3=3 , io4=4

    integer (ip) :: ok
    integer (ip) :: n , ns
    integer (ip) :: i , j , k , l , m
    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1

    real (dp)                   :: P , Kl , c1 , c2
    real (dp) , dimension (nrv) :: Kd

    real (dp) , allocatable , dimension (:,:,:)     :: wrk1
    real (dp) , allocatable , dimension (:,:,:,:)   :: logtempkelvin
    real (dp) , allocatable , dimension (:,:,:,:)   :: wrk2
    real (dp) , allocatable , dimension (:,:,:,:,:) :: DiffKJ


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    ns        = nrv
    Kl        = 0.5_dp * 1.0e-5_dp / thd % laminf
    Kd (1:ns) = 1.0e-4_dp * thd % p0 / ( thd % pinf * thd % Dinf (1:ns) )


    allocate ( DiffKJ        ( i0:i1 , j0:j1 , k0:k1 , 1:ns , 1:ns ) , &
               logtempkelvin ( 1:npoly , i0:i1 , j0:j1 , k0:k1     ) , &
               wrk1          ( i0:i1 , j0:j1 , k0:k1               ) , &
               wrk2          ( i0:i1 , j0:j1 , k0:k1 , 1:ns        ) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate prim_vis_vars_wo_mu'


    logtempkelvin ( io1 , i0:i1,j0:j1,k0:k1) = 1.0_dp                                        ! (log T)^0
    logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1) = log ( T (i0:i1,j0:j1,k0:k1) * thd % Tinf )    ! (log T)^1
    logtempkelvin ( io3 , i0:i1,j0:j1,k0:k1) = logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1 ) * &
                                               logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1 )     ! (log T)^2
    logtempkelvin ( io4 , i0:i1,j0:j1,k0:k1) = logtempkelvin ( io3 , i0:i1,j0:j1,k0:k1 ) * &
                                               logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1 )     ! (log T)^3


    ! Diffusion matrix


    DiffKJ (i0:i1,j0:j1,k0:k1,1:ns,1:ns) = 0.0_dp
    do m = 1 , ns
       do l = 1 , ns
          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   do n = 1 , npoly
                      DiffKJ (i,j,k,l,m) = DiffKJ (i,j,k,l,m)      + &
                                           thd % CofD (n,l,m)      * &
                                           logtempkelvin (n,i,j,k)
                   end do
                end do
             end do
          end do
       end do
    end do


    dm (i0:i1,j0:j1,k0:k1,1:ns) = 0.0_dp
    do m = 1 , ns

       wrk1 (i0:i1,j0:j1,k0:k1) = 0.0_dp

       do l = 1 , ns

          if ( m /= l ) then

             wrk1 (i0:i1,j0:j1,k0:k1) = wrk1 (i0:i1,j0:j1,k0:k1)                + &
                                      ( Xa (i0:i1,j0:j1,k0:k1,l) + epsm12_asp ) * &
                                        thd % Wc (l) * W_i (i0:i1,j0:j1,k0:k1)

             dm (i0:i1,j0:j1,k0:k1,m) = dm (i0:i1,j0:j1,k0:k1,m)                         + &
                                      ( Xa (i0:i1,j0:j1,k0:k1,l) + epsm12_asp )          * &
                                      (  1.0_dp / exp ( DiffKJ (i0:i1,j0:j1,k0:k1,m,l) ) )

          end if

       end do


       ! Kd(m)
       ! cgs unit : cm2/s
       ! passage to m2/s
       ! pressure adjustement pres en bar (1)


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                P            = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
                dm (i,j,k,m) = wrk1 (i,j,k) * Kd (m) / ( dm (i,j,k,m) * P )
             end do
          end do
       end do


    end do


    deallocate (DiffKJ)


    ! Conductivity


    wrk2 (i0:i1,j0:j1,k0:k1,1:ns) = 0.0_dp
    do l = 1 , ns
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                do m = 1 , npoly
                   wrk2 (i,j,k,l) = wrk2 (i,j,k,l)          + &
                                    thd % CofLam (m,l)      * &
                                    logtempkelvin (m,i,j,k)
                end do
             end do
          end do
       end do
    end do


    deallocate (logtempkelvin)


    wrk1 (i0:i1,j0:j1,k0:k1) = 0.0_dp
    ct   (i0:i1,j0:j1,k0:k1) = 0.0_dp
    do l = 1 , ns
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                c1           = exp ( wrk2 (i,j,k,l) )
                c2           = Xa (i,j,k,l)
                ct (i,j,k)   = ct (i,j,k)   + c2 / c1
                wrk1 (i,j,k) = wrk1 (i,j,k) + c2 * c1
             end do
          end do
       end do
    end do


    deallocate (wrk2)


    ! Kl
    ! unit = erg/(cm*K*s)
    ! passage to SI : J/(m.K.s)


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             ct (i,j,k) = Kl * ( 1.0_dp / ct (i,j,k) + wrk1 (i,j,k) )
          end do
       end do
    end do


    deallocate (wrk1)


  end subroutine prim_vis_vars_wo_mu


!> \brief Calculate thermal conductivity
!!
!! It provides optimization to the code.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine prim_vis_ct ( domain_id , thd , T , Xa , ct )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: Xa        !< molar mass array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: ct        !< thermal conductivity


    integer (ip) , parameter :: npoly = 4 ! order of the polynom
    integer (ip) , parameter :: io1=1 , io2=2 , io3=3 , io4=4

    integer (ip) :: ok
    integer (ip) :: n , ns
    integer (ip) :: i , j , k , l , m
    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1

    real (dp)                   :: Kl , c1 , c2

    real (dp) , allocatable , dimension (:,:,:)     :: wrk1
    real (dp) , allocatable , dimension (:,:,:,:)   :: logtempkelvin
    real (dp) , allocatable , dimension (:,:,:,:)   :: wrk2


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    ns        = nrv
    Kl        = 0.5_dp * 1.0e-5_dp / thd % laminf


    allocate ( logtempkelvin ( 1:npoly , i0:i1 , j0:j1 , k0:k1     ) , &
               wrk1          ( i0:i1 , j0:j1 , k0:k1               ) , &
               wrk2          ( i0:i1 , j0:j1 , k0:k1 , 1:ns        ) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate prim_vis_ct'


    logtempkelvin ( io1 , i0:i1,j0:j1,k0:k1) = 1.0_dp                                        ! (log T)^0
    logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1) = log ( T (i0:i1,j0:j1,k0:k1) * thd % Tinf )    ! (log T)^1
    logtempkelvin ( io3 , i0:i1,j0:j1,k0:k1) = logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1 ) * &
                                               logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1 )     ! (log T)^2
    logtempkelvin ( io4 , i0:i1,j0:j1,k0:k1) = logtempkelvin ( io3 , i0:i1,j0:j1,k0:k1 ) * &
                                               logtempkelvin ( io2 , i0:i1,j0:j1,k0:k1 )     ! (log T)^3


    ! Conductivity


    wrk2 (i0:i1,j0:j1,k0:k1,1:ns) = 0.0_dp
    do l = 1 , ns
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                do m = 1 , npoly
                   wrk2 (i,j,k,l) = wrk2 (i,j,k,l)          + &
                                    thd % CofLam (m,l)      * &
                                    logtempkelvin (m,i,j,k)
                end do
             end do
          end do
       end do
    end do


    deallocate (logtempkelvin)


    wrk1 (i0:i1,j0:j1,k0:k1) = 0.0_dp
    ct   (i0:i1,j0:j1,k0:k1) = 0.0_dp
    do l = 1 , ns
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                c1           = exp ( wrk2 (i,j,k,l) )
                c2           = Xa (i,j,k,l)
                ct (i,j,k)   = ct (i,j,k)   + c2 / c1
                wrk1 (i,j,k) = wrk1 (i,j,k) + c2 * c1
             end do
          end do
       end do
    end do


    deallocate (wrk2)


    ! Kl
    ! unit = erg/(cm*K*s)
    ! passage to SI : J/(m.K.s)


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             ct (i,j,k) = Kl * ( 1.0_dp / ct (i,j,k) + wrk1 (i,j,k) )
          end do
       end do
    end do


    deallocate (wrk1)


  end subroutine prim_vis_ct


!> \brief Calculate diffusion matrix.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine prim_vis_Dm ( domain_id , thd , v , W_i , T , Xa , dm )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: Xa        !< molar mass array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: dm        !< diffusion matrix


    integer (ip) , parameter :: npoly = 4 ! order of the polynom
    integer (ip) , parameter :: io1=1 , io2=2 , io3=3 , io4=4

    integer (ip) :: ok
    integer (ip) :: n , ns
    integer (ip) :: i , j , k , l , m
    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1

    real (dp)                       :: P
    real (dp) , dimension (nrv)     :: Kd
    real (dp) , dimension (1:npoly) :: logtempkelvin

    real (dp) , allocatable , dimension (:,:,:)     :: wrk1
    real (dp) , allocatable , dimension (:,:,:,:,:) :: DiffKJ


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    ns        = nrv
    Kd (1:ns) = 1.0e-4_dp * thd % p0 / ( thd % pinf * thd % Dinf (1:ns) )


    allocate ( DiffKJ        ( i0:i1 , j0:j1 , k0:k1 , 1:ns , 1:ns ) , &
               wrk1          ( i0:i1 , j0:j1 , k0:k1               ) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate prim_vis_Dm'


    ! Diffusion matrix

    logtempkelvin (io1) = 1.0_dp
    DiffKJ (i0:i1,j0:j1,k0:k1,1:ns,1:ns) = 0.0_dp
    do m = 1 , ns
       do l = 1 , ns
          do k = k0 , k1
             do j = j0 , j1
                do i = i0 , i1
                   logtempkelvin (io2) = log ( T (i,j,k) * thd % Tinf )
                   logtempkelvin (io3) = logtempkelvin (io2) * logtempkelvin (io2)
                   logtempkelvin (io4) = logtempkelvin (io3) * logtempkelvin (io2)
                   do n = 1 , npoly
                      DiffKJ (i,j,k,l,m) = DiffKJ (i,j,k,l,m) + &
                                           thd % CofD (n,l,m) * &
                                           logtempkelvin (n)
                   end do
                end do
             end do
          end do
       end do
    end do


    dm (i0:i1,j0:j1,k0:k1,1:ns) = 0.0_dp
    do m = 1 , ns

       wrk1 (i0:i1,j0:j1,k0:k1) = 0.0_dp

       do l = 1 , ns

          if ( m /= l ) then

             wrk1 (i0:i1,j0:j1,k0:k1) = wrk1 (i0:i1,j0:j1,k0:k1)                + &
                                      ( Xa (i0:i1,j0:j1,k0:k1,l) + epsm12_asp ) * &
                                        thd % Wc (l) * W_i (i0:i1,j0:j1,k0:k1)

             dm (i0:i1,j0:j1,k0:k1,m) = dm (i0:i1,j0:j1,k0:k1,m)                         + &
                                      ( Xa (i0:i1,j0:j1,k0:k1,l) + epsm12_asp )          * &
                                      (  1.0_dp / exp ( DiffKJ (i0:i1,j0:j1,k0:k1,m,l) ) )

          end if

       end do


       ! Kd(m)
       ! cgs unit : cm2/s
       ! passage to m2/s
       ! pressure adjustement pres en bar (1)


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                P            = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
                dm (i,j,k,m) = wrk1 (i,j,k) * Kd (m) / ( dm (i,j,k,m) * P )
             end do
          end do
       end do


    end do


    deallocate ( DiffKJ , wrk1 )


  end subroutine prim_vis_Dm


!> \brief Calculate molar fraction.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine molefrac ( domain_id , thd , W_i , v , Xa )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: Xa        !< molar mass array


    integer (ip) :: i , j , k , l
    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    do l = 1 , nrv
       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1
                Xa (i,j,k,l) = v (i,j,k, l+niv ) / &
                             ( v (i,j,k,1) * W_i (i,j,k) * thd % Wc (l) )
             end do
          end do
       end do
    end do


  end subroutine molefrac


!> \brief Force sum of mass fractions to be equal to unity.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine force_sumY ( domain_id , v )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    real (dp) , allocatable , dimension (:,:,:,:) , intent (inout) :: v         !< conserved variables array


    integer (ip) :: i , j , k , l
    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1
    real (dp)    :: rho , rho_i , sum , Ya , ps , Ya_max , Ya_min


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    if ( npv == 0 ) then


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho   = v (i,j,k,1)
                rho_i = 1.0_dp / rho
                sum   = 0.0_dp

                do l = 1 , nrv-1
                   Ya              = v (i,j,k,niv+l) * rho_i
                   Ya              = min ( max ( Ya , 0.0_dp ) , 1.0_dp )
                   v (i,j,k,niv+l) = rho * Ya
                   sum             = sum + Ya
                end do

                Ya = 1.0_dp - sum
                Ya = min ( max ( Ya , 0.0_dp ) , 1.0_dp )

                v (i,j,k,niv+nrv) = rho * Ya

             end do
          end do
       end do


    else


       do k = k0 , k1
          do j = j0 , j1
             do i = i0 , i1

                rho   = v (i,j,k,1)
                rho_i = 1.0_dp / rho
                sum   = 0.0_dp

                do l = 1 , nrv-1
                   Ya              = v (i,j,k,niv+l) * rho_i
                   Ya              = min ( max ( Ya , 0.0_dp ) , 1.0_dp )
                   v (i,j,k,niv+l) = rho * Ya
                   sum             = sum + Ya
                end do

                Ya = 1.0_dp - sum
                Ya = min ( max ( Ya , 0.0_dp ) , 1.0_dp )

                v (i,j,k,niv+nrv) = rho * Ya

                ! passive scalars
                do l = nrv+1 , nrv+npv
                   Ya              = v (i,j,k,niv+l) * rho_i
                   Ya              = min ( max ( Ya , 0.0_dp ) , 1.0_dp )
                   v (i,j,k,niv+l) = rho * Ya
                end do

             end do
          end do
       end do


    end if


    if ( nvv > 0 ) then

        do k = k0 , k1
           do j = j0 , j1
              do i = i0 , i1

                 rho   = v (i,j,k,1)
                 rho_i = 1.0_dp / rho
                 ps    = v (i,j,k,niv+5) * rho_i


                 ! Squared filtered scalar
                 !=======================================================
                    l = 8
                    Ya              = v (i,j,k,niv+l) * rho_i
                    Ya_min          = 0.0_dp
                    Ya_max          = 1.0_dp

                    Ya              = min ( max ( Ya , Ya_min ) , Ya_max )
                    v (i,j,k,niv+l) = rho * Ya


                 ! SGS scalar variance
                 !=======================================================
                    l = 6
                    Ya              = v (i,j,k,niv+l) * rho_i
                    Ya_min          = 0.0_dp
                    Ya_max          = ps * ( 1.0_dp - ps )

                    Ya              = min ( max ( Ya , Ya_min ) , Ya_max )
                    v (i,j,k,niv+l) = rho * Ya


                 ! Filtered squared scalar
                 !=======================================================
                    l = 7
                    Ya              = v (i,j,k,niv+l) * rho_i
                    Ya_min          = ps * ps
                    Ya_max          = ps * 1.0_dp

                    Ya              = min ( max ( Ya , Ya_min ) , Ya_max )
                    v (i,j,k,niv+l) = rho * Ya


                 ! Departure from maximal variance
                 !=======================================================
                    l = 9
                    Ya              = v (i,j,k,niv+l) * rho_i
                    Ya_min          = 0.0_dp
                    Ya_max          = ps * ( 1.0_dp - ps )

                    Ya              = min ( max ( Ya , Ya_min ) , Ya_max )
                    v (i,j,k,niv+l) = rho * Ya

              end do
           end do
        end do

    end if


  end subroutine force_sumY


!> \brief Calculate ratio of specific heats.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine gammamix ( domain_id , thd , W_i , cp , gamma )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: gamma     !< ratio of specific heats


    integer (ip) :: i , j , k

    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             gamma (i,j,k) = thd % gam2 * W_i (i,j,k) / cp (i,j,k)
             gamma (i,j,k) = 1.0_dp / ( 1.0_dp - gamma (i,j,k) )
          end do
       end do
    end do


  end subroutine gammamix


!> \brief Calculate speed of sound.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine soundmix ( domain_id , thd , W_i , cp , T , v , cs )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: cs        !< speed of sound


    integer (ip) :: i , j , k

    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1

    real (dp)    :: gamma , P


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )



    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             P          = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             gamma      = thd % gam2 * W_i (i,j,k) / cp (i,j,k)
             gamma      = 1.0_dp / ( 1.0_dp - gamma )
             cs (i,j,k) = sqrt ( gamma * P / v (i,j,k,1) )
          end do
       end do
    end do


  end subroutine soundmix


!> \brief Calculate Mach number.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine machmix ( domain_id , thd , W_i , cp , T , v , mach )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mach      !< Mach number


    integer (ip) :: i , j , k

    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1

    real (dp)    :: P , gamma , cs


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             P            = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             gamma        = thd % gam2 * W_i (i,j,k) / cp (i,j,k)
             gamma        = 1.0_dp / ( 1.0_dp - gamma )
             cs           = sqrt ( gamma * P / v (i,j,k,1) )
             mach (i,j,k) = sqrt ( v (i,j,k,2) * v (i,j,k,2)   + &
                                   v (i,j,k,3) * v (i,j,k,3)   + &
                                   v (i,j,k,4) * v (i,j,k,4) ) / &
                                   v (i,j,k,1)
             mach (i,j,k) = mach (i,j,k) / cs
          end do
       end do
    end do


  end subroutine machmix


!> \brief Calculate a _relative_ Mach number.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine machmix_rel ( domain_id , thd , W_i , cp , T , v , rux , rvy , rwz , mach )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: rux       !< relative x-component of the velocity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: rvy       !< relative y-component of the velocity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: rwz       !< relative z-component of the velocity
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: mach      !< Mach number


    integer (ip) :: i , j , k

    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1

    real (dp)    :: P , gamma , cs


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             P            = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             gamma        = thd % gam2 * W_i (i,j,k) / cp (i,j,k)
             gamma        = 1.0_dp / ( 1.0_dp - gamma )
             cs           = sqrt ( gamma * P / v (i,j,k,1) )
             mach (i,j,k) = sqrt ( rux (i,j,k) * rux (i,j,k)    + &
                                   rvy (i,j,k) * rvy (i,j,k)    + &
                                   rwz (i,j,k) * rwz (i,j,k) )  / &
                                   v (i,j,k,1)
             mach (i,j,k) = mach (i,j,k) / cs
          end do
       end do
    end do


  end subroutine machmix_rel


!> \brief Calculate total temperature.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine Ttmix ( domain_id , thd , W_i , cp , T , v , Tt )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    type (thd_type) , intent (in)                                  :: thd       !< thermodynamic derived type
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: W_i       !< inverted molar mass
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: T         !< temperature
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)    :: v         !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: Tt        !< total temperature


    integer (ip) :: i , j , k

    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1

    real (dp)    :: P , gamma , cs , mach , wrk


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1

             P          = v (i,j,k,1) * T (i,j,k) * W_i (i,j,k)
             gamma      = thd % gam2 * W_i (i,j,k) / cp (i,j,k)
             gamma      = 1.0_dp / ( 1.0_dp - gamma )
             cs         = sqrt ( gamma * P / v (i,j,k,1) )
             mach       = sqrt ( v (i,j,k,2) * v (i,j,k,2)   + &
                                 v (i,j,k,3) * v (i,j,k,3)   + &
                                 v (i,j,k,4) * v (i,j,k,4) ) / &
                                 v (i,j,k,1)
             mach       = mach / cs

             wrk        = 0.5_dp * ( gamma - 1.0_dp )
             Tt (i,j,k) = T (i,j,k) * ( mach * mach * wrk + 1.0_dp )

          end do
       end do
    end do


  end subroutine Ttmix


!> \brief Calculate Prandtl number.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine prandtlmix ( domain_id , cp , mu , ct , prandtl )


    integer (ip) , intent (in)                                     :: domain_id !< subdomain selection
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: cp        !< heat capacity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: mu        !< dynamic viscosity
    real (dp) , allocatable , dimension (:,:,:) , intent (in)      :: ct        !< thermal conductivity
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)   :: prandtl   !< Prandtl number


    integer (ip) :: i , j , k

    integer (ip) :: i0 , i1 , j0 , j1 , k0 , k1


    call domain_select ( domain_id , i0 , i1 , j0 , j1 , k0 , k1 )


    do k = k0 , k1
       do j = j0 , j1
          do i = i0 , i1
             prandtl (i,j,k) = cp (i,j,k) * mu (i,j,k) / ct (i,j,k)
          end do
       end do
    end do


  end subroutine prandtlmix


end module thermodynamics
