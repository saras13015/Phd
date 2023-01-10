      SUBROUTINE PSR_SIMP_CHEMKIN_QI
     *           ( TINIT , RHOINIT , PINIT , WINIT , YINIT ,
     *             DTINIT , MYINDEX , RKR , RKF ,
     *             ROP , WSCAL , YASS )
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
      DOUBLE PRECISION WINIT                    !> IN    - Mean molar mass (g/mol)
      DOUBLE PRECISION YINIT ( 2 * KMAX + 10 )  !> INOUT - IN  -> Species mass fraction (g/g) ; 
!                                                          OUT -> Chemical production rates of each species (1/s)
      DOUBLE PRECISION DTINIT                   !> IN    - Reference time (s)
      INTEGER          MYINDEX                  !> IN    - Index scalar (-)
      DOUBLE PRECISION RKF(*),RKR(*)            !> OUT   - Forward/Reversed rate of progress for the reaction of the species (1/s)
      DOUBLE PRECISION ROP(*)                   !> OUT   - Rate of progress variable for the reactions (1/s)
      DOUBLE PRECISION WSCAL (*)                !> OUT   - Rate of progress of the MYINDEX-th species on each elementary step (1/s)
      DOUBLE PRECISION YASS  (*)                !> OUT   - Mass fraction of steady state species (g/g)
!
!
      INTEGER K
      DOUBLE PRECISION MYCONSTANT , OTHERCONSTANT
!
!
      MYCONSTANT    = DTINIT / RHOINIT
      OTHERCONSTANT = MYCONSTANT * WINIT
!
      Z (1) = TINIT
      DO K = 2 , KK + 1
         Z (K) = YINIT (K-1)
      ENDDO
!
!
      ! CALCULATING PRODUCTION RATES SUPER FAST WITHOUT USING DVODE
      CALL CKWYR_qi ( RHOINIT , Z(1) , Z(2) , IWORK , RWORK ,
     *                RWORK(NWDOT) , MYINDEX , RKR , RKF , ROP , WSCAL ,
     *                YASS )

!      CALL CKWYR ( RHOINIT , Z(1),Z(2) , IWORK , RWORK , RWORK (NWDOT) )
!
!
      DO K = 1 , KK ! KK is the number of species
         YINIT (K) = RWORK ( NWDOT + K - 1 ) * ! CGS units
     *               RWORK ( NWT   + K - 1 ) * ! Molar mass of species K
     *               MYCONSTANT                ! YINIT is the nondimensional increment of the species K
      ENDDO
!
!
      DO K = 1 , II ! II is the number of elementary reaction steps
         WSCAL (K) = WSCAL (K)                   * ! CGS units
     *               RWORK ( NWT + MYINDEX - 1 ) * ! Molar mass of species MYINDEX
     *               MYCONSTANT                    ! MYCONSTANT is nondimensional
      END DO
!
!
      DO K = 1 , II ! II is the number of elementary reaction steps
         ROP (K) = ROP (K)       * ! CGS units
     *             OTHERCONSTANT   ! multiply by the molar mass of the mixture
         RKR (K) = RKR (K)       * ! CGS units
     *             OTHERCONSTANT   ! multiply by the molar mass of the mixture
         RKF (K) = RKF (K)       * ! CGS units
     *             OTHERCONSTANT   ! multiply by the molar mass of the mixture
      END DO
!
!
      END SUBROUTINE PSR_SIMP_CHEMKIN_QI
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE COSINESS_to_CKCDRKFR ( RKF , RKR , CDOT , DDOT )
!
!
      IMPLICIT NONE
!
!
      INCLUDE 'param.common'
      INCLUDE 'chemkin.common'

!
      DOUBLE PRECISION RKF(*) , RKR(*)   ! IN    - forward/reversed local rate of progress
      DOUBLE PRECISION CDOT(*) , DDOT(*) ! INOUT - creation/destruction rate
!
!
      ! CALCULATING CREATION AND DESTRUCTION RATES SUPER FAST WITHOUT USING DVODE
      CALL CKCDRKFR ( RKF , RKR , IWORK , RWORK , CDOT , DDOT )
!
!
      END SUBROUTINE COSINESS_to_CKCDRKFR
