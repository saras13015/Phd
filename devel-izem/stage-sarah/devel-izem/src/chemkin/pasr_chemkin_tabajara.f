      SUBROUTINE PASR_CHEMKIN_TABAJARA ( TINIT , RHOINIT , PINIT ,
     *                         YINIT , DTINIT , SZINIT , ISTATE )
!
!
      IMPLICIT NONE
!
!
      INCLUDE 'param.common'
      INCLUDE 'chemkin.common'
!
!
      DOUBLE PRECISION TINIT                    !> INOUT - Initial/Final temperature (K)
      DOUBLE PRECISION RHOINIT                  !> IN    - Density for constant volume reaction (g/cm^3)
      DOUBLE PRECISION PINIT                    !> IN    - Pressure for constant pressure reaction (g/cm/s2 = dynes/cm2)
      DOUBLE PRECISION YINIT ( 2 * KMAX + 10 )  !> INOUT - Species mass fraction (g/g)
      DOUBLE PRECISION DTINIT                   !> IN    - Integration time step (s)
      DOUBLE PRECISION SZINIT                   !> IN    - Scalar segregation rate (-)
      INTEGER          ISTATE                   !> INOUT - Index to specify the state of the calculation
!
!
      INTEGER K
!
!
! Common parameters
      RHO_  = RHOINIT ! Constant volume reaction
      P     = PINIT   ! Constant pressure reaction
      OMS   = 1.0D0 - SZINIT
      Z (1) = TINIT  ! T: [K]
      DO K = 2 , KK + 1
         Z (K) = YINIT (K-1)  ! Yk: [-]
      ENDDO
!
!
      CALL PASR_REACTION_V_TABAJARA ( DTINIT , ISTATE ) ! Constant volume reaction
!      CALL PASR_REACTION_P_TABAJARA ( DTINIT , ISTATE ) ! Constant pressure reaction
!
!
      TINIT = Z (1) ! New temperature
      DO K = 1 , KK
         YINIT (K) = Z ( K + 1 ) ! New species
      ENDDO
!
!
      END SUBROUTINE PASR_CHEMKIN_TABAJARA
