      SUBROUTINE PASR_REACTION_V_TABAJARA (TF,ISTATE)
!
!
      IMPLICIT NONE
!
!
      INCLUDE 'chemkin.common'
      INCLUDE 'param.common'
!
!
      DOUBLE PRECISION TI , TF
!
      INTEGER ISTATE , IBOUCLE , NBOUCLEMAX , MF
      PARAMETER ( NBOUCLEMAX = 20 )
      EXTERNAL PASR_FUN_TABAJARA , PASR_JAC_TABAJARA
!
!
      TI      = 0.0D0
      MF      = 22
      IBOUCLE = 0
      ISTATE  = 1
!
!
C$$$      write (*,*) 'before dvode'
C$$$      write (*,*) RTOL , ATOL
C$$$      WRITE (*,*) Z(1) , Z(2) , Z(3) , Z(4) , ISTATE , TI , TF
!
!
      DO WHILE ( (ISTATE .NE. 2) .AND. (IBOUCLE .LT. NBOUCLEMAX) )
!
         ISTATE  = 1 ! OVERWRITING ISTATE=-1 TO AVOID STOP THE PROGRAM
         IBOUCLE = 1 + IBOUCLE
!
         CALL DVODE ( PASR_FUN_TABAJARA , NEQ , Z , TI    , TF     ,
     *                ITOL , RTOL          , ATOL , ITASK , ISTATE ,
     *                IOPT , RWORK (NVODE) , LRW  , IWORK (IVODE)  ,
     *                LIW  , PASR_JAC_TABAJARA , MF , RWORK , IWORK )
!
C$$$        WRITE (*,*) 'boucle' , IBOUCLE
C$$$        WRITE (*,*) Z(1) , Z(2) , Z(3) , Z(4) , ISTATE , TI , TF
!
      END DO
!
!
!      IF ( IBOUCLE .GE. NBOUCLEMAX ) THEN
!         WRITE (*,'(1X,A,I5,3X,A,I5)') 'ISTATE = ' , ISTATE ,
!     *   'NO CHEMICAL SOLUTION IN PASR_REACTION_V_TABAJARA FUNCTION, loop number= ' ,
!     *   IBOUCLE
!      END IF
!
!
!      IF ( IBOUCLE .GT. 1 ) THEN
!!         WRITE ( LOUT , * )
!         WRITE (*,*)
!     *   'CHEMICAL SOLUTION IN PASR_REACTION_V_TABAJARA FUNCTION, loop number = ' , IBOUCLE
!      ENDIF
!
!
      IF ( ISTATE .LE. -2 ) THEN
         WRITE (*,'(A,I5)') 'PSR: ISTATE = ' , ISTATE
         RETURN ! to stop the simulation with MPI_ABORT in CFD solver
      END IF
!
!
      END SUBROUTINE PASR_REACTION_V_TABAJARA
!
!===============================================================================
!
      SUBROUTINE PASR_FUN_TABAJARA ( N , TIME , Z , ZP , RPAR , IPAR )
!
!
      IMPLICIT NONE
!
!
      INCLUDE 'param.common'
!
!
      DOUBLE PRECISION RHO_ , RU   , Z  , ZP  , RPAR ,
     1                 RHO  , SUM  , H  , CPB , WDOT ,
     2                 WT   , TIME , P_ , CVB , P    ,
     3                 OMS  , EPS0
!
      INTEGER IPAR , KK , NWT , NH , NWDOT , K , N
!
      COMMON /RCONS/ P , RU , RHO_ , OMS , EPS0
      COMMON /ICONS/ KK , NWT , NH , NWDOT
      DIMENSION Z(*) , ZP(*) , RPAR(*) , IPAR(*)
!
      DOUBLE PRECISION RHO_I , MRU
      DOUBLE PRECISION GAMMA_STAR
!
!
      CALL CKCVBS ( Z(1) , Z(2) , IPAR , RPAR , CVB )                  ! cv: [erg / (g * K)]
      CALL CKWYR  ( RHO_ , Z(1) , Z(2) , IPAR , RPAR , RPAR (NWDOT) )  ! wdot: [mol / (cm^3 * s)]
      CALL CKHMS  ( Z(1) , IPAR , RPAR , RPAR (NH) )                   ! hk: [erg / g]
!
      RHO_I      = 1.0D0 / RHO_ * OMS  ! rho^{-1}: [cm^3 / g]
      MRU        = - RU                ! -R: [erg / (mol * K)]
      GAMMA_STAR = 0.9D0
!
      SUM = 0.0D0
      DO K = 1 , KK
!
         H            = RPAR ( NH    + K - 1 )               ! hk: [erg / g]
         WDOT         = GAMMA_STAR * RPAR ( NWDOT + K - 1 )  ! gammastar*wdot: [mol / (cm^3 * s)]
         WT           = RPAR ( NWT   + K - 1 )               ! MWk: [g / mol]
         ZP   ( K+1 ) = WDOT * WT * RHO_I                    ! dYk/dt: [1 / s]
         SUM          = SUM + H * WDOT * WT + MRU * Z(1) * WDOT
!
      END DO
!
!
      ZP   (1) = - SUM / ( RHO_ * CVB ) * OMS  ! dT/dt: [K / s]
!
!
      END SUBROUTINE PASR_FUN_TABAJARA
!
!===============================================================================
!
      SUBROUTINE PASR_JAC_TABAJARA
!
!
      IMPLICIT NONE
!
!
      END SUBROUTINE PASR_JAC_TABAJARA
