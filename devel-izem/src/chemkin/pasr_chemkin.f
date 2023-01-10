      SUBROUTINE PASR_CHEMKIN ( TINIT , RHOINIT , PINIT ,
     1                          YINIT , DTINIT , SZINIT ,
     2                          TSTARINIT , YSTARINIT ,
     3                          GAMMASTARINIT , TAUMINIT ,
     4                          HFINIT , EPS0INIT , ISTATE )
!
!
      IMPLICIT NONE
!
!
      INCLUDE 'param.common'
      INCLUDE 'chemkin.common'
!
!
      DOUBLE PRECISION TINIT                       !> INOUT - Initial/Final temperature (K)
      DOUBLE PRECISION RHOINIT                     !> IN    - Density for constant volume reaction (g/cm^3)
      DOUBLE PRECISION PINIT                       !> IN    - Pressure for constant pressure reaction (g/cm/s2 = dynes/cm2)
      DOUBLE PRECISION YINIT ( 2 * KMAX + 10 )     !> INOUT - Species mass fraction (g/g)
      DOUBLE PRECISION DTINIT                      !> IN    - Integration time step (s)
      DOUBLE PRECISION SZINIT                      !> IN    - Scalar segregation rate (-)
      DOUBLE PRECISION TSTARINIT                   !> INOUT - Initial/Final temperature at the fine-scale structure (K)
      DOUBLE PRECISION YSTARINIT ( 2 * KMAX + 10 ) !> INOUT - Species mass fraction at the fine-scale structure (g/g)
      DOUBLE PRECISION GAMMASTARINIT               !> IN    - Fine-scale structure volume fraction (-)
      DOUBLE PRECISION TAUMINIT                    !> IN    - Subgrid mixing time scale (s)
      DOUBLE PRECISION HFINIT ( 2 * KMAX + 10 )    !> IN    - Species formation enthalpy (erg/g)
      DOUBLE PRECISION EPS0INIT                    !> IN    - Tolerance for zero value
      INTEGER          ISTATE                      !> INOUT - Index to specify the state of the calculation
!
!
      INTEGER K
!
!
! Common parameters
      RHO_     = RHOINIT  ! Constant Volume reaction
      P        = PINIT    ! Constant Pressure reaction
      OMS      = 1.0D0 - SZINIT
      Z (1)    = TINIT      ! T: [K]
      Z (2+KK) = TSTARINIT  ! T*: [K]
      DO K = 1 , KK
         Z (1+K)    = YINIT (K)      ! Yk: [-]
         Z (2+KK+K) = YSTARINIT (K)  ! Yk*: [-]
         HF (K)     = HFINIT (K)     ! Formation enthalpy: [erg / g]
      ENDDO
      GAMMASTAR = GAMMASTARINIT  ! [-]
      TAUM      = TAUMINIT       ! [s]
      EPS0      = EPS0INIT
!
!
      IF ( ABS(GAMMASTAR) >= EPS0 ) THEN  ! PSR / U-PaSR
         CALL PASR_REACTION_V ( DTINIT , ISTATE ) ! Constant volume reaction
!         CALL PASR_REACTION_P ( DTINIT , ISTATE ) ! Constant pressure reaction
!      ELSE  ! No chemical reaction
!         Z(2+KK) = 0.0D0  ! No T*
!         DO K = 1 , KK
!            Z(2+KK+K) = 0.0D0  ! No Yk*
!         END DO
      END IF
!
!
      TINIT     = Z (1)     ! New temperature
      TSTARINIT = Z (2+KK)  ! New temperature at the fine-scale structure
      DO K = 1 , KK
         YINIT (K)     = Z (1+K)     ! New species
         YSTARINIT (K) = Z (2+KK+K)  ! New species at the fine-scale structure
      ENDDO
!
!
      END SUBROUTINE PASR_CHEMKIN
