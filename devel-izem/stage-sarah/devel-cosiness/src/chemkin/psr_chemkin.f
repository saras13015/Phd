      SUBROUTINE PSR_CHEMKIN ( TINIT , RHOINIT , PINIT ,
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
!                                                          if ( ISTATE .LE. -2 ) the 
!
!
      INTEGER K
!
!
      ! these are a common parameters
      RHO_  = RHOINIT ! Constant Volume reaction
      P     = PINIT   ! Constant Pressure reaction
      OMS   = 1.0D0 - SZINIT
      Z (1) = TINIT
      DO K = 2 , KK + 1
         Z (K) = YINIT (K-1)
      ENDDO
      ! these are a common parameters
!
!
      CALL REACTION_V (DTINIT,ISTATE) ! Constant Volume reaction
!      CALL REACTION_P (DTINIT,ISTATE) ! Constant Pressure reaction
!
!
      TINIT = Z (1) ! new temperature
      DO K = 1 , KK
         YINIT (K) = Z ( K + 1 ) ! new species
      ENDDO
!
!
      END SUBROUTINE PSR_CHEMKIN
