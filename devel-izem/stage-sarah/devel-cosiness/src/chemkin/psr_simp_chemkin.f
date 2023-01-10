      SUBROUTINE PSR_SIMP_CHEMKIN ( TINIT , RHOINIT , PINIT ,
     *                              YINIT , DTINIT )
!
!
      IMPLICIT NONE
!
!
      INCLUDE 'param.common'
      INCLUDE 'chemkin.common'
!
!
      DOUBLE PRECISION TINIT                    !> IN    - Initial/Final temperature (K)
      DOUBLE PRECISION RHOINIT                  !> IN    - Density for constant volume reaction (g/cm^3)
      DOUBLE PRECISION PINIT                    !> IN    - Pressure for constant pressure reaction (g/cm/s2 = dynes/cm2)
      DOUBLE PRECISION YINIT ( 2 * KMAX + 10 )  !> INOUT - IN  -> Species mass fraction (g/g) , 
!                                                          OUT -> Mass fractions departure: Yfinal-Yinitial (g/g)
      DOUBLE PRECISION DTINIT                   !> IN    - Integration time step (s)
!
!
      INTEGER K
      DOUBLE PRECISION MYCONSTANT
!
!
      MYCONSTANT = DTINIT / RHOINIT
      Z (1) = TINIT
      DO K = 2 , KK + 1
         Z (K) = YINIT (K-1)
      ENDDO
!
!
      ! CALCULATING PRODUCTION RATES SUPER FAST WITHOUT USING DVODE
      CALL CKWYR (RHOINIT, Z(1), Z(2), IWORK, RWORK, RWORK(NWDOT)) ! Constant Volume reaction
!      CALL CKWYP (PINIT   , Z(1), Z(2), IWORK, RWORK, RWORK(NWDOT)) ! Constant Pressure reaction


      ! MASS FRACTIONS DEPARTURE
      DO K = 1 , KK
         YINIT (K) = RWORK ( NWDOT + K - 1 ) * ! CGS units
     .               RWORK ( NWT   + K - 1 ) * ! Molar mass of species K
     .               MYCONSTANT                ! YINIT is the nondimensional increment of the species K
!         write (*,*) K , YINIT (K)
      ENDDO
!
!
      END SUBROUTINE PSR_SIMP_CHEMKIN
