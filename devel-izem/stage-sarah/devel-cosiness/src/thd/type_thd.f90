module type_thd


  use paramthd


  implicit none


  type thd_type


     integer(IA) :: Nreac    ! number of elementary reactions
     integer(IA) :: Nspc     ! number of effective species
     integer(IA) :: Nele     ! number of effective elements

     real(RA)   :: R        ! perfect gas constant (8.314D0 J/mol/K)
     real(RA)   :: T0       ! reference temperature (es=INT_T0^T Cp dT (K))
     real(RA)   :: P0       ! reference pressure (Pa = pascal)
     real(RA)   :: Tinf     ! infinite mixture : temperature (=T0)
     real(RA)   :: Pinf     ! infinite mixture : pressure (= P0)
     real(RA)   :: Cpinf    ! infinite mixture : heat capacity (= T0)
     real(RA)   :: Winf     ! infinite mixture : molar mass
     real(RA)   :: Rinf     ! infinite mixture : Rinf = Winf*Cpinf
     real(RA)   :: Rhoinf   ! infinite mixture : density
     real(RA)   :: Soundinf ! infinite mixture : sound velocity
     real(RA)   :: Muinf    ! infinite mixture : viscosity (kg/(m.s))
     real(RA)   :: Laminf   ! infinite mixture : Thermal conductivity W/(K.s)
     real(RA)   :: Prinf    ! infinite mixture : Prandtl number      (-)
     real(RA),dimension(Nspcmax_thd)::IntCpInf    ! infinite mixture :
     real(RA),dimension(Nspcmax_thd)::Dinf  ! infinite mixture : Diffusion
     real(RA),dimension(Nspcmax_thd)::Scinf ! infinite mixture : Schmidt number
     real(RA),dimension(Nspcmax_thd)::Leinf ! infinite mixture : Lewis number

     real(RA)    :: Tc2p         ! code2physical adim
     real(RA)    :: Tp2c         !
     real(RA)    :: Pc2p         !
     real(RA)    :: Pp2c         !
     real(RA)    :: Xc2p
     real(RA)    :: Xp2c
     real(RA)    :: Uc2p
     real(RA)    :: Up2c

     real(RA)    :: Gamma        ! Gamma=Cp/Cv, reference mixture
     real(RA)    :: GaM1         ! = Gamma-1
     real(RA)    :: GaM2         ! = ( Gamma -1 ) / Gamma

     real(RA), dimension(Npolymax_thd)             :: T0I
     ! T0I : powers of the reference temperature
     real(RA), dimension(Npolymax_thd,Nspcmax_thd) :: TLI
     ! TLI : powers of the limit  temperature

     character(len=Length_thd), dimension(Nspcmax_thd) :: NameSpecy
     character, dimension(Nspcmax_thd)                 :: Phase
     integer(IA)                                       :: UnrSpc
     ! UnrSpc : number of the unresolved specy if necessary

     real(RA), dimension(Nspcmax_thd) :: Tpolymin   ! T minimum polynom
     real(RA), dimension(Nspcmax_thd) :: Tpolymax   ! T maximum polynom
     real(RA), dimension(Nspcmax_thd) :: Tpolylim   ! T limite pol change
     real(RA), dimension(Nspcmax_thd) :: h0f        ! formation enthalpy
     real(RA), dimension(Nspcmax_thd) :: W          ! molar mass
     real(RA), dimension(Nspcmax_thd) :: Wc         ! molar mass dim code
     real(RA), dimension(Nspcmax_thd) :: Wc_i       ! molar mass dim code (inv)
     real(RA), dimension(Npolymax_thd,2,Nspcmax_thd) :: HeatCap
     ! HeatCap : Polynom coeff for Cp_molar/R (CHEMKIN FORMAT)
     real(RA), dimension(NTranspolymax_thd,Nspcmax_thd,Nspcmax_thd) :: CofD
     ! CofD   : Polynom coeff for Diffusion
     real(RA), dimension(NTranspolymax_thd,Nspcmax_thd)             :: CofLam
     ! CofLam : Polynom coeff for Lambda (energy diffusion⁾
     real(RA), dimension(NTranspolymax_thd,Nspcmax_thd)             :: CofEta
     ! CofEta : Polynom coeff for Eta, viscosity

     ! normalized data (InPetto Code)
     real(RA)                       ::Rc  ! perfect gas constant: (8.314D0 J/mol/K)/CpInf
     real(RA)                       ::T0c ! reference temperature
     real(RA),dimension(Nspcmax_thd)::Sc
     real(RA),dimension(Nspcmax_thd)::Le
     real(RA),dimension(Nspcmax_thd)::InfMixY!composition of infinite mixture
     real(RA),dimension(Nspcmax_thd)::CpConsc! constant Cp, code format
     real(RA),dimension(Nspcmax_thd)::h0fc   ! normalized formation enthalpy

     ! Initials conditions (Reactor)
     real(RA)                       :: Tempic  ! Initial Temperature (K)

     ! added by PM (30-04-2012) to account for mixture tracking scalar (Inga Mahle)
     real (RA)    , dimension (Nelemax_thd)             :: AWT , AWTc ! atomic weights
     integer (IA) , dimension (Nelemax_thd,Nspcmax_thd) :: KNCF       ! element count matrix


   end type thd_type


end module type_thd
