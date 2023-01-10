module thdtools
!
!
  use parallel
  implicit none
!
!
contains
!
!
!
  subroutine thd_init ( Thd , adi , thdname )
!
!
!
    use paramthd
!
    use type_thd
!
    use adim
!
    use commonfile
!
    use filetools
!
!
!
    implicit none
!
!
!
    type (adi_type)   , intent ( in )          :: adi
    type (thd_type)   , intent ( inout )       :: Thd
    character (len=*) , intent ( in )          :: thdname ! = 'tranfit.out'
    !- local declaration
    !-
    integer(IA) :: i , l , npoly=4
    character (len=Length_thd) :: thdfilename , cpfilename
    ! integer(IA) :: i , io , n , l , nt , npoly=4
    ! real(RA)    :: Tmin , T , CpT
    ! character (len=Length_thd) :: thdfilename , cpfilename
!
!
!-
!- read polynomial coefficient from tranfit code
    call thd_cofd_read ( thdname    , Thd % Nreac , Thd % Nspc     ,  &
                         thd % Nele , Thd % NameSpecy              ,  &
                         thd % Cofd , Thd % CofLam , Thd % CofEta  ,  &
                         Thd % W    , thd % AWT    , thd % KNCF       )
!write (*,*) 'thd_cofd_read... OK'
!
!
!-
! 5000 format(A15)
! 5005 format(5(E15.0))
! 5004 format(4(E15.0))
      !
      ! - Reference parameters
      !
    Thd % R  = adi % R_ref
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DO NOT MODIFY
    Thd % T0 = 298.15_RA    ! The reference to calculate the formation enthalpy
    Thd % P0 = 101325.0_RA  ! 1 atmosphere
    ! DO NOT MODIFY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    Thd % T0I (1) = Thd % T0
    do i = 2 , Npolymax_thd
       Thd % T0I (i) = Thd % T0 ** i
    enddo
    !
    ! - Unresolved species
    !
    !
    ! - Infinite mixture
    !
    Thd % Tinf = adi % t_ref
    Thd % Pinf = adi % p_ref
    !
    ! - End data reading
    !
    ! call trline ( 'thd init' )
    ! call frline ( 'Reference Temp' , Thd % T0 )
    ! call frline ( 'Reference Pres' , Thd % P0 )
    !   call crline('Unresolved specy',Thd % NameSpecy(Thd % UnrSpc))
    !
    !  Physical and Thermodynamical data
    !
!    call trline ( 'Thermodynamical data' )
    call thd_file_read_therm ( 'therm.dat'   , Thd % Nspc     , Thd % NameSpecy , &
                               Thd % HeatCap , Thd % Phase    , Thd % Tpolymin  , &
                               Thd % Tpolymax, Thd % Tpolylim )
    !
    do i = 1 , Thd % Nspc
       call GetSpeciesFileName(thdfilename,cpfilename,Thd % NameSpecy(i))
       ! call crline('Specy',Thd % NameSpecy(i))
       ! call crline('Phase',Thd % Phase(i))
       ! call frline('Molar mass (kg/mol)',Thd % W(i))
       !
      ! Break temperature and infinite mixture heat capacity
      !
       Thd % IntCpInf (i) = Thd % HeatCap (1,1,i) *         &
                          ( Thd % Tpolylim (i) - Thd % T0 )
       Thd % TLI(1,i)     = Thd % Tpolylim (i)
       do l = 2 , npoly+1
          Thd % TLI (l,i)    = Thd % Tpolylim (i) ** l
          Thd % IntCpInf (i) = Thd % IntCpInf (i) + Thd % HeatCap (l,1,i) *        &
                             ( Thd % TLI (l,i) - Thd % T0I (l) ) / real ( l , RA )
       enddo
       !
       ! Formation enthalpy
       !
       ! VERIFER POUR LES CAS CP=CTE SI OK
       Thd % h0f (i)     = Thd % HeatCap (npoly+2,1,i) + &
                           Thd % HeatCap (1,1,i) * Thd % T0
       do l = 2 , npoly+1
          Thd % h0f (i) = Thd % h0f (i) + Thd % HeatCap (l,1,i) * &
                          Thd % T0I (l) / real ( l , RA )
       enddo
       Thd % h0f(i)     = Thd % h0f (i) * Thd % R / Thd % W (i)
       !call frline('Form. enth (J/kg)',Thd % h0f(i))
       !
       ! Plot Cp versus Temp
       !
       ! nt = 200
       ! Tmin = Thd % Tpolymin (i) - 50.0_RA
       ! open( unit = unit_data_file , file = cpfilename , form = 'formatted' )
       ! 5002 format( 2(E15.6) )
       ! do n = 1 , nt
       !    T = Tmin + real( n-1 , RA ) * ( Thd % Tpolymax (i) - Tmin ) / &
       !        real( nt-1 , RA )
       !    if( T < Thd % Tpolylim (i) ) then
       !       CpT = Thd % HeatCap (1,1,i) + Thd % HeatCap (2,1,i) * T
       !       do io = 2 , npoly
       !          CpT = CpT + Thd % HeatCap (io+1,1,i) * T ** io
       !       enddo
       !    else
       !       CpT = Thd % HeatCap (1,2,i) + Thd % HeatCap (2,2,i) * T
       !       do io = 2 , npoly
       !          CpT = CpT + Thd % HeatCap (io+1,2,i) * T ** io
       !       enddo
       !    endif
       !    CpT = CpT * Thd % R / Thd % W (i)
       !    write ( unit_data_file , 5002 ) T , CpT
       ! enddo
       ! close ( unit_data_file )
    enddo
    !
    ! call file_mark('reacteur.data','# INITIALS CONDITIONS #',iunit_tic)
    ! if(iunit_tic > 0) then
    !    read(iunit_tic,*) Thd % Tempic
    ! else
    !    Thd % Tempic=0.0_RA
    ! endif
    Thd % Tempic = 0.0_RA   !!! A supprimer de ce module !
    !
    !
    ! Infinite mixture, normalization data
    !
    ! CpPhys  = 0.0_RA
    ! SumW    = 0.0_RA
    ! do i = 1, Thd % Nspc
    !   l       = 1 ! if we are below the limit temperature
    !   if(Thd % Tinf > Thd % Tpolylim(i)) l = 2 ! otherwise
    !   CpPhysk = Thd % HeatCap(1,l,i)+Thd % HeatCap(2,l,i)*Thd % Tinf
    !   do io = 2,npoly
    !         CpPhysk   = CpPhysk + Thd % HeatCap(io+1,l,i)*Thd % Tinf**io
    !   enddo
    !   !CpPhysk = (J.kg-1.mol-1) // *R because database given in Cp/R
    !   CpPhysk         = CpPhysk * Thd % R/Thd % W(i)
    !   Thd % CpConsc(i)  = CpPhysk
    !   CpPhys          = CpPhys  + CpPhysk * Thd % InfMixY(i)
    !   SumW            = SumW    + Thd % InfMixY(i)/Thd % W(i)
    ! enddo
    Thd % Winf        = adi % W_ref ! 1.0_RA/SumW ! Winf (kg.mol-1) masse molaire melange infini
    Thd % Cpinf       = adi % cp_ref !CpPhys      ! Cp (J.kg-1.K-1) : capacite cal melange infini
    ! (R/W) <=> Cpinf and (W) <=> Winf then (R) <=> Cpinf*Winf
    Thd % Rinf        = adi % R_ref !Thd%Cpinf*Thd%Winf
    Thd % Gamma       = adi % gamma_inf !1.0_RA/(1.0_RA-Thd%R/Thd%Cpinf/Thd%Winf)
    Thd % GaM1        = Thd % Gamma - 1.0_RA
    Thd % GaM2        = Thd % GaM1 / thd % Gamma
    ! if(Runflag%Acoustic == 0) then
    Thd % Tc2p          = 1.0_RA  ! LMN formulation
    Thd % Tp2c          = 1.0_RA  ! LM Nformulation
    Thd % Pc2p          = 1.0_RA  ! LMN formulation
    Thd % Pp2c          = 1.0_RA  ! LMN formulation
    ! elseif(Runflag % Acoustic == 1) then
    !  Thd % Tc2p          = Thd % GaM1        ! Ac formulation
    !  Thd % Tp2c          = 1.0_RA/Thd % GaM1 ! Ac formulation
    !  Thd % Pc2p          = Thd % Gamma       ! Ac formulation
    !  Thd % Pp2c          = 1.0_RA/Thd % Gamma! Ac formulation
    ! else
    !  call abort_mpi ('Thd_init needs to know Runflag.Acoustic')
    ! endif
    ! if((Runflag % VarConsEnt /= 0).and.(Runflag % VarConsEnt /= 1))                   &
    !  call abort_mpi ('Thd_init needs to know Runflag.VarConsEnt')
    ! if((Runflag % VarConsSca /= 0).and.(Runflag % VarConsSca /= 1))                   &
    !  call abort_mpi ('Thd_init needs to know Runflag.VarConsSca')
    ! !
    Thd % Rhoinf                  = adi % rho_ref !Thd%Pinf/((Thd%R/Thd%Winf)*Thd%Tinf)
    Thd % Soundinf                = adi % c_inf !sqrt(Thd%Gamma*(Thd%R/Thd%Winf)*Thd%Tinf)
    Thd % T0c                     = Thd % Tp2c * Thd % T0 / Thd % Tinf
    Thd % Rc                      = Thd % R / Thd % Rinf
    Thd % Wc ( 1 : Thd % Nspc )   = Thd % W ( 1 : Thd % Nspc ) / Thd % Winf
    Thd % Wc_i ( 1 : Thd % Nspc ) = 1.0_dp / Thd % Wc ( 1 : Thd % Nspc )
    Thd % AWTc ( 1 : Thd % Nele ) = Thd % AWT ( 1 : Thd % Nele ) / Thd % Winf
    ! Cp constante adimension code (Cpk/Cpinf)
    Thd % CpConsc ( 1 : Thd % Nspc ) = Thd % CpConsc ( 1 : Thd % Nspc ) / Thd % Cpinf
    ! Enthalpie de formation adimension code
    Thd % h0fc ( 1 : Thd % Nspc ) = Thd % h0f ( 1 : Thd % Nspc ) / &
                                  ( adi % u_ref * adi % u_ref )
    !Thd % h0fc ( 1 : Thd % Nspc )= Thd % h0f ( 1 : Thd % Nspc ) * Thd%Tp2c / (Thd%Cpinf*Thd%Tinf)
    !
    !call thddt_Transport_inf(Thd)
    !
    Thd % Muinf                   = adi % mu_ref
    Thd % Dinf ( 1:Thd % Nspc )   = adi % D_ref
    Thd % Laminf                  = adi % lbda_ref
    !
    Thd % Prinf                   = adi % pr
    Thd % Leinf ( 1:Thd % Nspc )  = adi % le
    Thd % Scinf ( 1:Thd % Nspc )  = thd % Leinf ( 1 : Thd % Nspc ) * thd % Prinf
    !
    !
    !
  contains
!
!
    subroutine GetSpeciesFileName ( thdfn , cpfn , Sname )
!
!
      implicit none
!
!
      character(len=*) , intent (in)  :: Sname
      character(len=*) , intent (out) :: thdfn , cpfn
!
!
      cpfn  = 'thd/thd_'//trim(adjustl(Sname))//'.cpte'
      thdfn = 'thd_'//trim(adjustl(Sname))//'.data'
!
!
    end subroutine GetSpeciesFileName
!
!
  end subroutine thd_init



  subroutine thd_display (Thd)
!
!
    use paramthd
!
    use type_thd
!
    use commonfile
!
    use filetools
!
!
    implicit none
!
!
    type (thd_type) , intent (in)  :: Thd
    !- local declaration
    !-
    integer (IA) :: i
    !-
    !-
    call trline ( 'Thermodynamical data' )
    do i = 1 , Thd % Nspc
       ! call crline ( 'Specy' , Thd % NameSpecy (i) )
       ! call crline ( 'Phase' , Thd % Phase (i) )
       ! call frline ( 'Molar mass (kg/mol)' , Thd % W (i) )
       ! call frline ( 'Form. enth (J/kg)' , Thd % h0f (i) )
       ! call frline ( 'Scinf (-)' , Thd % Scinf (i) )
       ! call frline ( 'Leinf (-)' , Thd % Leinf (i) )
    end do
    do i = 1 , thd % nele
!       call frline ( 'Atom Weig (kg/mol)' , Thd % AWT (i) )
    enddo
!    call trline ( 'Infinite Mixture' )
    do i = 1 , Thd % Nspc
       if ( Thd % InfMixY (i) > 0 ) then
!          call crline ( 'Specy:' , Thd % NameSpecy (i) )
!          call frline ( 'infinit Y' ,Thd % InfMixY (i) )
       endif
    enddo
    ! call trline ( 'Unresolved specy' )
    ! call irline ( 'Number' , Thd % UnrSpc )
    ! ! call crline('name',Thd%NameSpecy(Thd%UnrSpc))
    ! call trline ( 'dim thermo. data')
    ! call frline ( 'Ref Temp 0 (K)' , Thd % T0 )
    ! call frline ( 'Ref Pres 0 (Pa)' , Thd % P0 )
    ! call frline ( 'Temp inf (K)' , Thd % Tinf )
    ! call frline ( 'Pres inf (Pa)' , Thd % Pinf )
    ! call frline ( 'R (J/mol.K)' , Thd % R )
    ! call frline ( 'Cp inf (J/kg.K)' , Thd % Cpinf )
    ! call frline ( 'W inf (kg/mol)' , Thd % Winf )
    ! call frline ( 'R inf (J/kg.K)' , Thd % Rinf )
    ! call frline ( 'Gamma inf (-)' , Thd % Gamma )
    ! call frline ( 'Rho inf (kg/m^3)' , Thd % RhoInf )
    ! call frline ( 'Sound inf (m/s)' , Thd % Soundinf )
    ! call frline ( 'Visc Dyn inf (kg/m.s)' , Thd % Muinf )
    if ( Thd % RhoInf /= 0.0_RA ) then
       ! call frline ( 'Visc Cin inf (m^2/s)' , Thd % Muinf / Thd % RhoInf )
       ! call frline ( 'Cd. therm inf (W/K.s)' , Thd % Laminf )
       ! call frline ( 'Prandtl inf (-)' , Thd % Prinf )
       ! call trline ( 'code thermo. data' )
       ! call frline ( 'T0c' , Thd % T0c )
       ! call frline ( 'Rc' , Thd % Rc )
       ! call frline ( 'Tc2p' , Thd % Tc2p )
       ! call frline ( 'Pc2p' , Thd % Pc2p )
    end if
 !
 !
  end subroutine thd_display



  subroutine thd_species_index ( Thd , name , number )
!
!
!
    use paramthd

    use type_thd
!
!
!
    implicit none
!
!
!
    type (thd_type) , intent (in)   :: Thd
    character (len=*) , intent (in) :: name
    integer(IA) , intent (out)      :: number
!
    integer(IA) :: i , flag
!
!
!
    i     = 0
    flag  = 0
    do while ( flag == 0 )
       i = i + 1
       if ( name == Thd % NameSpecy (i) ) then
          flag   = 1
          number = i
       end if
       if ( ( i == Thd % Nspc ) .and. ( flag==0 ) ) flag = -1
    end do
    if ( flag == -1 ) then
       write (*,*) 'WARNING undefined species: ' , name
       number = -1
!       call abort_mpi ('I will stop now the program')
    end if
!
!
  end subroutine thd_species_index



  subroutine thd_file_read_therm ( filename , Nspc     , NameSpecy , HeatCap  , &
                                   Phase    , Tpolymin , Tpolymax  , Tpolylim )
!
!
    use paramthd
!
    use commonfile
!
    use filetools
!
!
    implicit none
!
!
    integer(IA) , intent (in)       :: Nspc
    character (len=*) , intent (in) :: filename
    character (len=Length_thd) , dimension (Nspcmax_thd) , intent (in) :: NameSpecy
    character , dimension (Nspcmax_thd) , intent (out) :: Phase
    real(RA)  , dimension (Nspcmax_thd) , intent (out) :: Tpolymin , Tpolymax , Tpolylim
    real(RA)  , dimension (Npolymax_thd,2,Nspcmax_thd) , intent(out) :: HeatCap
!
!
    character (len=Length_thd) :: readmark , localNameSpecy
    character (len=18)         :: readspec
    character (len=6)          :: readdate
    character (len=20)         :: readjunk
    character                  :: readphase
    real(RA)    :: readTmin , readTmax , readTlim
    logical     :: endread , flagspec
    integer(IA) :: err , l , flag , i , check
    real(RA)    :: Tmin , Tmax , Tbreak
!-
    Tpolymin = 0.0_RA ; Tpolymax = 0.0_RA ; Tpolylim = 0.0_RA ; HeatCap = 0.0_RA ;
!-
    do i = 1 , Nspc
       localNameSpecy = NameSpecy (i)
       open ( unit = unit_data_file , file = filename , form = 'formatted' , &
              status = 'old' , iostat = err )
       if ( err /= 0 ) then
          flag = -1
       else
          endread = .false.
          read ( unit_data_file , '(a50)' ) readmark
          if ( readmark /= 'THERMO' .and. readmark /= 'THERMO ALL' ) then
             flag    = -3
             close ( unit_data_file )
             endread = .true.
          else
             read ( unit_data_file , '(3(F10.0))' ) Tmin , Tbreak , Tmax
          endif
          do while ( .not. endread )
             read ( unit_data_file , '(a18,a6,a20,a1,e10.0,e10.0,e10.0,I5)' ) &
                  readspec , readdate , readjunk , readphase ,                &
                  readTmin , readTmax , readTlim , check
             !! prevoir le cas où readTmin,readTmax,readTlim
             !! n'existent pas et utiliser le defaut
             if ( check /= 1 )                                                         &
                  call stopline ( 'format error while reading therm.dat : thd_file_read_therm' )
             if ( readspec == 'END' ) then
                flag    = -2 ! did not found the specy
                endread = .true.
                exit
             else
                flagspec = .true.
                l        = 1
                do while ( (readspec (l:l) /= ' ' ) .and. flagspec )
                   if ( readspec (l:l) /= localNameSpecy (l:l) ) then
                      flagspec = .false.
                   else
                      l = l + 1
                      if ( readspec (l:l) == ' ' .and. localNameSpecy (l:l) /= '' ) &
                           flagspec = .false.
                   endif
                enddo
                if ( flagspec ) then
                   flag         = unit_data_file
                   endread      = .true.
                   Phase (i)    = readphase
                   Tpolymin (i) = readTmin
                   Tpolymax (i) = readTmax
                   Tpolylim (i) = readTlim
                   read ( unit_data_file , 5005 ) HeatCap (1,2,i) , HeatCap (2,2,i) , &
                                HeatCap (3,2,i) , HeatCap (4,2,i) , HeatCap (5,2,i)
                   read ( unit_data_file , 5005 ) HeatCap (6,2,i) , HeatCap (7,2,i) , &
                                HeatCap (1,1,i) , HeatCap (2,1,i) , HeatCap (3,1,i)
                   read ( unit_data_file , 5004 ) HeatCap (4,1,i) , HeatCap (5,1,i) , &
                                                  HeatCap (6,1,i) , HeatCap (7,1,i)
                   exit
                else
                   read ( unit_data_file , * )
                   read ( unit_data_file , * )
                   read ( unit_data_file , * )
                endif
             endif
          enddo
       endif
!
       close ( unit_data_file )
!
    enddo
!
!
5005 format(5(E15.0))
5004 format(4(E15.0))
!
!
  end subroutine thd_file_read_therm




 subroutine thd_cofd_read ( filename , Nreac , Nspc , Nele , NameSpecy , Cofd , Coflam , Cofeta , W , &
                             AWT , KNCF )
!
!
!
    use paramthd
!
    use commonfile
!
    use filetools
!
!
!
    implicit none
!
!
!
    character (len=*) , intent (in)  :: filename
    character (len=Length_thd) , dimension (Nspcmax_thd), intent (out)  :: NameSpecy
    real (RA) , dimension (NTranspolymax_thd,Nspcmax_thd,Nspcmax_thd) ,         &
                                                           intent (out) :: Cofd
    real(RA) , dimension (NTranspolymax_thd,Nspcmax_thd) , intent (out) :: Coflam
    real(RA) , dimension (NTranspolymax_thd,Nspcmax_thd) , intent (out) :: Cofeta
    real(RA) , dimension (Nspcmax_thd) , intent (out)                   :: W
    integer(IA) , intent (out)                                          :: Nspc , Nele , Nreac
    integer (IA) , dimension (Nelemax_thd,Nspcmax_thd) , intent (out)   :: KNCF
    real(RA) , dimension (Nelemax_thd) , intent (out)                   :: AWT
!
!
    character (len=16) , allocatable , dimension (:)   :: SpcA
    double precision , allocatable , dimension (:,:,:) :: cofd_lu
    double precision , allocatable , dimension (:,:)   :: coflam_lu
    double precision , allocatable , dimension (:)     :: W_lu
    integer (IA) :: err , j , n , k , j_lu
    integer (IA) :: Nspec_lu , Npoly_lu , Nele_lu , Nreac_lu
    real (RA)    :: ter
!
!
    call exist_file (filename)
    open ( unit = unit_data_file , file = filename , & !, form = 'unformatted' , &
           status = 'old' , iostat = err )
    if ( err /= 0 ) then
       call abort_mpi ('File not found thd_cofd_read')
    else
       read (unit_data_file,*) Nreac_lu
       read (unit_data_file,*) Nspec_lu
       read (unit_data_file,*) Npoly_lu
       read (unit_data_file,*) Nele_lu
       if ( Nelemax_thd < Nele_lu ) call stopline ( 'error number max of atoms' )
       if ( NTranspolymax_thd < Npoly_lu ) call stopline ( 'error poly order thd_cofd_read' )
       if ( Nspcmax_thd < Nspec_lu ) then
          call irline ( 'Nspcmax_thd' , Nspcmax_thd )
          call irline ( 'Nspec_lu' , Nspec_lu )
          call stopline ( 'error Nspcmax_thd<Nspec_lu thd_cofd_read' )
       endif
       Nreac = Nreac_lu
       Nspc = Nspec_lu
       Nele = Nele_lu
       allocate ( cofd_lu (Npoly_lu,Nspec_lu,Nspec_lu) , &
                  coflam_lu (Npoly_lu,Nspec_lu) )
       allocate ( SpcA (Nspec_lu) , W_lu (Nspec_lu) )
       ! call trline ( 'Open Thd Data File' )
       ! call irline ( 'Number Of Species' , Nspc )
       ! call irline ( 'Polynomial order' , Npoly_lu )
       !
       ! Diffusion coefficient
       !
!       call trline( 'Diff coefficients' )
       do j = 1 , Nspec_lu
          read (unit_data_file,*) j_lu
          if ( j_lu /= j ) call abort_mpi ('error file format thd_cofd_read')
          read( unit_data_file,* ) SpcA (j)
          do k = 1 , j
             read(unit_data_file,*) SpcA (k)
             do n = 1 , Npoly_lu
                read (unit_data_file,*) ter ; cofd_lu (n,j,k) = real ( ter , RA )
             enddo
          enddo
          !write(*,8110) (SpcA(j), SpcA(k), (cofd_lu(n,j,k),n=1,Npoly_lu), k=1,j)
       enddo
       read (unit_data_file,*) j_lu
       if ( j_lu /= -2 ) call abort_mpi ('error intermediate file thd_cofd_read')
       do j = 1 , Nspec_lu
          NameSpecy (j)                    = SpcA (j)
          Cofd (1:Npoly_lu,j,1:j)          = cofd_lu (1:Npoly_lu,j,1:j)
          Cofd (1:Npoly_lu,j,j+1:Nspec_lu) = cofd_lu (1:Npoly_lu,j+1:Nspec_lu,j)
       enddo
       !
       ! coflam
       !
!       call trline('Lambda coefficients')
       do k = 1 , Nspec_lu
          read (unit_data_file,*) SpcA (k)
          do n = 1 , Npoly_lu
             read (unit_data_file,*) ter ; coflam_lu (n,k) = real ( ter , RA )
          enddo
       enddo
       Coflam (1:Npoly_lu,1:Nspec_lu) = coflam_lu (1:Npoly_lu,1:Nspec_lu)
       !write(*,8100) (SpcA(K), (Coflam(N,K),N=1,Npoly_lu), K=1,Nspec_lu)
       read (unit_data_file,*) j_lu
       if ( j_lu /= -2 ) call abort_mpi ('error intermediate file thd_cofd_read')
       !
       ! cofeta
       !
!       call trline ( 'Eta coefficients' )
       Coflam_lu (:,:) = 0.0_RA
       do k = 1 , Nspec_lu
          read (unit_data_file,*) SpcA (k)
          do n = 1 , Npoly_lu
             read (unit_data_file,*) ter ; coflam_lu (n,k) = real ( ter , RA )
          enddo
       enddo
       Cofeta (1:Npoly_lu,1:Nspec_lu) = coflam_lu (1:Npoly_lu,1:Nspec_lu)
       !write(*,8100)  (SpcA(k), (Cofeta(n,k),n=1,Npoly_lu), k=1,Nspec_lu)
       read (unit_data_file,*) j_lu
       if ( j_lu /= -2 ) call abort_mpi ('error intermediate file thd_cofd_read')
       !
       ! cofeta
       !
!       call trline ( 'Molecular Weigths' )
       W (:) = 0.0_RA
       do k = 1 , Nspec_lu
          read (unit_data_file,*) SpcA (k)
          read (unit_data_file,*) ter ; W_lu (k) = real ( ter , RA )
       enddo
       W (1:Nspec_lu) = W_lu (1:Nspec_lu) / 1000.0_RA    !  en kg/mol
       !write(*,8120)  (SpcA(k), W_lu(k), k=1,Nspec_lu)
       read (unit_data_file,*) j_lu
       deallocate ( cofd_lu , coflam_lu , SpcA , W_lu )
       !
       ! AWT atomic weights
       !
!       call trline ( 'Atomic Weigths' )
       do k = 1 , Nele_lu
          read (unit_data_file,*) ter ;  AWT (k) = real ( ter , RA )
       enddo
       AWT (:) = AWT (:) / 1000.0_RA   ! en kg/mol
       read (unit_data_file,*) j_lu
       if ( j_lu /= -2 ) call abort_mpi ('error intermediate file thd_cofd_read')
       !
       ! KNCF matrix
       !
!       call trline ( 'KNCF matrix' )
       do k = 1 , Nspec_lu
          do j = 1 , Nele_lu
             read (unit_data_file,*) KNCF (j,k)
          end do
       enddo
       read (unit_data_file,*) j_lu
       if ( j_lu /= -1 ) call abort_mpi ('error end file thd_cofd_read')
!       call trline ( 'Thd File OK' )
    endif
    close (unit_data_file)
    !
    !
! 8100 FORMAT (2X, A10, 4E12.3)
! 8110 FORMAT (2X, A10, A10, 4E12.3)
! 8120 FORMAT (2X, A10, E12.3, ' g/mol')
!
!
  end subroutine thd_cofd_read


  subroutine thd_write ( thd , WriteFileName )
!
!
    use commonfile
!
    use filetools
!
    use paramthd
!
    use type_thd
!
!
    implicit none
!
!
    character (len=50) , intent (in) :: WriteFileName
    type (thd_type) , intent (in)        :: Thd
    ! local declaration
    integer (IA)       :: formatfile = 20061 , Lastdata = -1 ,  &
                          ns , nele , nt , np , is , iss , it , ipp

 !
 ! call filefix_getname(FileName,simulname,'thd')
!
!
!     ! PM, avoiding problems if istore = 0
!     if ( istore == 0 ) then
!        istore = 1
! !       WriteFileName = 'thd'//na(istore)//'_'//na(myb)//'.fix'
!        WriteFileName = 'thd.dat'
!        istore = 0
!     else
! !       WriteFileName = 'thd'//na(istore)//'_'//na(myb)//'.fix'
!        WriteFileName = 'thd.dat'
!     end if
    !
 !
    open ( unit = unit_data_file , file = WriteFileName )! , form = 'unformatted' )
 ! general data
    write (unit_data_file,*) formatfile
 !
    write (unit_data_file,*) Ndim_thd
    write (unit_data_file,*) Nspcmax_thd
    write (unit_data_file,*) Nelemax_thd
    write (unit_data_file,*) Npolymax_thd
    write (unit_data_file,*) NTranspolymax_thd
    write (unit_data_file,*) Length_thd
 !
    write (unit_data_file,*) Thd % Nspc
    write (unit_data_file,*) Thd % Nele
 !
    ns = Thd % Nspc
    nele = Thd % Nele
    nt = NTranspolymax_thd
    np = Npolymax_thd
 !
    write (unit_data_file,*) Thd % R
    write (unit_data_file,*) Thd % T0
    write (unit_data_file,*) Thd % P0
    write (unit_data_file,*) Thd % Tinf
    write (unit_data_file,*) Thd % Pinf
    write (unit_data_file,*) Thd % Cpinf
    write (unit_data_file,*) Thd % Winf
    write (unit_data_file,*) Thd % Rinf
    write (unit_data_file,*) Thd % Rhoinf
    write (unit_data_file,*) Thd % Soundinf
    write (unit_data_file,*) Thd % Muinf
    write (unit_data_file,*) ( Thd % Dinf (is) , is = 1,ns )
    write (unit_data_file,*) Thd % Laminf
    write (unit_data_file,*) Thd % Prinf
    write (unit_data_file,*) ( Thd % Scinf (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % Leinf (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % IntCpInf (is) , is = 1,ns )
 !
    write (unit_data_file,*) Thd % Gamma
    write (unit_data_file,*) Thd % GaM1
    write (unit_data_file,*) Thd % GaM2
    write (unit_data_file,*) Thd % Tc2p
    write (unit_data_file,*) Thd % Tp2c
    write (unit_data_file,*) Thd % Pc2p
    write (unit_data_file,*) Thd % Pp2c
    write (unit_data_file,*) Thd % Xc2p
    write (unit_data_file,*) Thd % Xp2c
    write (unit_data_file,*) Thd % Uc2p
    write (unit_data_file,*) Thd % Up2c
 !
    write (unit_data_file,*) ( Thd % T0I (ipp) , ipp = 1,np )
    write (unit_data_file,*) ( ( Thd % TLI (ipp,is) , ipp = 1,np ) , is = 1,ns )
    write (unit_data_file,*) ( Thd % NameSpecy (is) , is = 1,ns )
    do is=1,ns
       write (unit_data_file,*) Thd % Phase (is)
    end do
    write (unit_data_file,*) Thd % UnrSpc
    write (unit_data_file,*) ( Thd % Tpolymin (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % Tpolymax (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % Tpolylim (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % h0f (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % W (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % Wc (is), is = 1,ns )
    write (unit_data_file,*) ( Thd % Wc_i (is), is = 1,ns )
    write (unit_data_file,*) ( ( Thd % HeatCap (ipp,1,is) , ipp = 1,np ) , is = 1,ns )
    write (unit_data_file,*) ( ( Thd % HeatCap (ipp,2,is) , ipp = 1,np ) , is = 1,ns )
 !
    write (unit_data_file,*) ( ( ( Thd % CofD (it,is,iss) , it = 1,nt ) , is = 1,ns ) , iss = 1,ns )
    write (unit_data_file,*) ( ( Thd % CofLam (it,is) , it = 1,nt ) , is = 1,ns )
    write (unit_data_file,*) ( ( Thd % CofEta (it,is) , it = 1,nt ) , is = 1,ns )
 !
    write (unit_data_file,*) Thd % Rc
    write (unit_data_file,*) Thd % T0c
    write (unit_data_file,*) ( Thd % Sc (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % Le (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % InfMixY (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % CpConsc (is) , is = 1,ns )
    write (unit_data_file,*) ( Thd % h0fc (is) , is= 1,ns )
    write (unit_data_file,*) Thd % Tempic
!
    write (unit_data_file,*)   ( Thd % AWT  (is) , is = 1,nele )
    write (unit_data_file,*)   ( Thd % AWTc (is) , is = 1,nele )
    write (unit_data_file,*) ( ( Thd % KNCF (is,it) , is = 1,nele ) , it = 1,ns )

!
    write (unit_data_file,*) Lastdata
!
    close (unit_data_file)
 !
    ! call srline
    ! call wrline ('.')
    ! call crline ( 'Ecriture fichier' , WriteFileName )
    ! call irline ( 'Species actives' , Thd % Nspc )
    ! call wrline ('.')
    ! call srline
!
!
  end subroutine thd_write
!
!
!
!
!
  subroutine thd_read ( Thd , ReadFileName )
!
!
    use commonfile
!
    use filetools
!
    use paramthd
!
    use type_thd
!
!   use commonrun
    implicit none
    character (len=50) , intent (in)     :: ReadFileName
    type (thd_type) , intent (out)       :: Thd
    ! local declaration
    integer (IA)       :: formatfile , Lastdata
    integer (IA)       :: Ndim_lu , Nspcmax_lu , Nelemax_lu , Npolymax_lu
    integer (IA)       :: NTranspolymax_lu , Length_lu
    integer (IA)       :: ns , nele , nt , np , is , iss , ipp , it

!
!
!     ! PM, avoiding problems if istore = 0
!     ! call filefix_getname(FileName,Rundata%SimulName,'thd')
!     if ( istore == 0 ) then
!        istore = 1
! !       ReadFileName = 'thd'//na(istore)//'_'//na(myb)//'.fix'
!        ReadFileName = 'thd.dat'
!        istore = 0
!     else
! !       ReadFileName = 'thd'//na(istore)//'_'//na(myb)//'.fix'
!        ReadFileName = 'thd.dat'
!     end if

 !
    call thd_reset (Thd)
 !
    ! call srline
    ! call wrline ('.')
    ! call crline ( 'Lecture fichier' , ReadFileName )
    ! call wrline ('.')
    ! call srline
 !
    call exist_file (Readfilename)
    open ( unit = unit_data_file , file = ReadFileName )! , form = 'unformatted' )
    ! general data
    read (unit_data_file,*) formatfile
    !
    if(formatfile == 20061) then
       read (unit_data_file,*) Ndim_lu
       if ( Ndim_lu /= Ndim_thd )                   call abort_mpi ('error Ndim thd_read')
       read (unit_data_file,*) Nspcmax_lu
       if ( Nspcmax_lu /= Nspcmax_thd )             call abort_mpi ('error Nspcmax thd_read')
       read (unit_data_file,*) Nelemax_lu
       if ( Nelemax_lu /= Nelemax_thd )             call abort_mpi ('error Nelemax thd_read')
       read (unit_data_file,*) Npolymax_lu
       if ( Npolymax_lu /= Npolymax_thd )           call abort_mpi ('error Npolymax thd_read')
       read (unit_data_file,*) NTranspolymax_lu
       if ( NTranspolymax_lu /= NTranspolymax_thd ) call abort_mpi ('error NTranspolymax thd_read')
       read (unit_data_file,*) Length_lu
       if ( Length_lu /= Length_thd )               call abort_mpi ('error Length thd_read')
  !
       read(unit_data_file,*) Thd % Nspc
       read(unit_data_file,*) Thd % Nele
  !
       ns = Thd % Nspc
       nele = Thd % Nele
       nt = NTranspolymax_thd
       np = Npolymax_thd
  !
       read (unit_data_file,*) Thd % R
       read (unit_data_file,*) Thd % T0
       read (unit_data_file,*) Thd % P0
       read (unit_data_file,*) Thd % Tinf
       read (unit_data_file,*) Thd % Pinf
       read (unit_data_file,*) Thd % Cpinf
       read (unit_data_file,*) Thd % Winf
       read (unit_data_file,*) Thd % Rinf
       read (unit_data_file,*) Thd % Rhoinf
       read (unit_data_file,*) Thd % Soundinf
       read (unit_data_file,*) Thd % Muinf
       read (unit_data_file,*) ( Thd % Dinf (is) , is = 1,ns )
       read (unit_data_file,*) Thd % Laminf
       read (unit_data_file,*) Thd % Prinf
       read (unit_data_file,*) ( Thd % Scinf (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % Leinf (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % IntCpInf (is) , is=1,ns )
  !
       read (unit_data_file,*) Thd % Gamma
       read (unit_data_file,*) Thd % GaM1
       read (unit_data_file,*) Thd % GaM2
       read (unit_data_file,*) Thd % Tc2p
       read (unit_data_file,*) Thd % Tp2c
       read (unit_data_file,*) Thd % Pc2p
       read (unit_data_file,*) Thd % Pp2c
       read (unit_data_file,*) Thd % Xc2p
       read (unit_data_file,*) Thd % Xp2c
       read (unit_data_file,*) Thd % Uc2p
       read (unit_data_file,*) Thd % Up2c
  !
       read (unit_data_file,*) ( Thd % T0I (ipp) , ipp = 1,np )
       read (unit_data_file,*) ( ( Thd % TLI (ipp,is) , ipp = 1,np ) , is = 1,ns )
  !
       read (unit_data_file,*) ( Thd % NameSpecy (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % Phase (is) , is = 1,ns )
       read (unit_data_file,*) Thd % UnrSpc
       read (unit_data_file,*) ( Thd % Tpolymin (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % Tpolymax (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % Tpolylim (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % h0f (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % W (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % Wc (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % Wc_i (is) , is = 1,ns )
       read (unit_data_file,*) ( ( Thd % HeatCap(ipp,1,is) , ipp = 1,np ) , is = 1,ns )
       read (unit_data_file,*) ( ( Thd % HeatCap(ipp,2,is) , ipp = 1,np ) , is = 1,ns )
  !
       read (unit_data_file,*) ( ( ( Thd % CofD(it,is,iss),it=1,nt) , is = 1,ns ),iss=1,ns )
       read (unit_data_file,*) ( ( Thd % CofLam(it,is),it=1,nt) , is = 1,ns )
       read (unit_data_file,*) ( ( Thd % CofEta(it,is),it=1,nt) , is = 1,ns )
  !
       read (unit_data_file,*) Thd % Rc
       read (unit_data_file,*) Thd % T0c
       read (unit_data_file,*) ( Thd % Sc (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % Le (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % InfMixY (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % CpConsc (is) , is = 1,ns )
       read (unit_data_file,*) ( Thd % h0fc (is) , is = 1,ns )
       read (unit_data_file,*) Thd % Tempic
!
       read (unit_data_file,*)   ( Thd % AWT (is)  , is = 1,nele)
       read (unit_data_file,*)   ( Thd % AWTc (is) , is = 1,nele)
       read (unit_data_file,*) ( ( Thd % KNCF (is,it) , is = 1,nele ) , it = 1,ns )

!
       read (unit_data_file,*) Lastdata
       if(Lastdata /= -1) call abort_mpi ('error LastData thd_read')
    else
       call abort_mpi ('error formatfile thd_read')
    endif
!
    close (unit_data_file)

 !
  end subroutine thd_read




 subroutine thd_reset (Thd)
!
!
   use paramthd
!
   use type_thd
!
!
   implicit none
!
!
   type (thd_type)    :: Thd
!- local declaration
   integer (IA) :: i
!-
   Thd % Nspc          = 0
   thd % Nele          = 0
 !
   Thd % R             = 0.0_RA
   Thd % T0            = 0.0_RA
   Thd % P0            = 0.0_RA
   Thd % Tinf          = 0.0_RA
   Thd % Pinf          = 0.0_RA
   Thd % Cpinf         = 0.0_RA
   Thd % Winf          = 0.0_RA
   Thd % Rinf          = 0.0_RA
   Thd % Rhoinf        = 0.0_RA
   Thd % Soundinf      = 0.0_RA
   Thd % Muinf         = 0.0_RA
   Thd % Dinf (:)      = 0.0_RA
   Thd % Scinf (:)     = 0.0_RA
   Thd % Leinf (:)     = 0.0_RA
   Thd % Prinf         = 0.0_RA
   Thd % Laminf        = 0.0_RA
   Thd % Gamma         = 0.0_RA
   Thd % GaM1          = 0.0_RA
   Thd % GaM2          = 0.0_RA
   Thd % Tc2p          = 0.0_RA
   Thd % Tp2c          = 0.0_RA
   Thd % Pc2p          = 0.0_RA
   Thd % Pp2c          = 0.0_RA
   Thd % Xc2p          = 0.0_RA
   Thd % Xp2c          = 0.0_RA
   Thd % Uc2p          = 0.0_RA
   Thd % Up2c          = 0.0_RA
 !
   Thd % T0I (:)       = 0.0_RA
   Thd % TLI (:,:)     = 0.0_RA
 !
 do i = 1 , Nspcmax_thd
    Thd % NameSpecy (i)  = 'None'
    Thd % Phase (i)      = '0'
 enddo
 Thd % UnrSpc        = 0
 !
 Thd % Tpolymin (:)    = 0.0_RA
 Thd % Tpolymax (:)    = 0.0_RA
 Thd % Tpolylim (:)    = 0.0_RA
 Thd % h0f (:)         = 0.0_RA
 Thd % W (:)           = 0.0_RA
 Thd % Wc (:)          = 0.0_RA
 Thd % Wc_i (:)        = 0.0_RA
 Thd % HeatCap (:,:,:) = 0.0_RA
 !
 Thd % CofD (:,:,:)   = 0.0_RA
 Thd % CofLam (:,:)   = 0.0_RA
 Thd % CofEta (:,:)   = 0.0_RA
 !
 Thd % Rc            = 0.0_RA
 Thd % T0c           = 0.0_RA
 Thd % Sc (:)        = 0.0_RA
 Thd % Le (:)        = 0.0_RA
 Thd % InfMixY (:)   = 0.0_RA
 Thd % Cpconsc (:)   = 0.0_RA
 Thd % h0fc (:)      = 0.0_RA
 Thd % IntCpInf (:)  = 0.0_RA
 !
 !
 thd % AWT  (:)      = 0.0_RA
 thd % AWTc (:)      = 0.0_RA
 thd % KNCF (:,:)    = 0
 !
!
 end subroutine thd_reset

!
!
end module thdtools
