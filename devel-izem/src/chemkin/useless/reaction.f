      SUBROUTINE REACTION ( TF , DT )
!
!
      IMPLICIT NONE
!
!
      INCLUDE 'chemkin.common'
      INCLUDE 'param.common'
!
!
      DOUBLE PRECISION TI , TF , DT
!
      INTEGER ISTATE , IBOUCLE , NBOUCLEMAX , MF
!
      PARAMETER ( NBOUCLEMAX = 20 )
!
      EXTERNAL FUN , JAC
!
!
      TI      = 0.0D0
      MF      = 22
      IBOUCLE = 0
      ISTATE  = 1
!
!
c$$$      write (*,*) 'before dvode'
c$$$      write (*,*) RTOL , ATOL
c$$$      WRITE (*,*) Z(1) , Z(2) , Z(3) , Z(4) , ISTATE , TI , TF
!
      DO WHILE ( ( ISTATE .NE. 2 ) .AND. ( IBOUCLE .LT. NBOUCLEMAX ) )
!
         ISTATE  = 1 ! OVERWRITING ISTATE=-1 TO AVOID STOP THE PROGRAM
         IBOUCLE = 1 + IBOUCLE
!
         CALL DVODE ( FUN  , NEQ             , Z    , TI    , TF     ,
     *                ITOL , RTOL            , ATOL , ITASK , ISTATE ,
     *                IOPT , RWORK (NVODE)   , LRW  , IWORK (IVODE)  ,
     *                LIW  , JAC             , MF   , RWORK  , IWORK )
!
c$$$        WRITE (*,*) 'boucle' , IBOUCLE
c$$$        WRITE (*,*) Z(1) , Z(2) , Z(3) , Z(4) , ISTATE , TI , TF
!
      ENDDO
!
!
      IF ( IBOUCLE .GE. NBOUCLEMAX ) THEN
         WRITE (*,*) 'ISTATE = ' , ISTATE ,
     *   'PAS DE SOLUTION CHIMIQUE dans REACTION.F, Nbre boucle= ' ,
     *   IBOUCLE
      ENDIF
!
!
      IF ( IBOUCLE .GT. 1 ) THEN
         WRITE ( LOUT , * )
     *   'SOLUTION CHIMIQUE dans REACTION.F, Nbre boucle = ' , IBOUCLE
      ENDIF
!
!
      IF ( ISTATE .LE. -2 ) THEN
         WRITE ( LOUT , * ) 'ISTATE = ' , ISTATE
         STOP
      ENDIF
!
!
!
      END SUBROUTINE REACTION
!
!--------------------------------------------------------
!--------------------------------------------------------
!--------------------------------------------------------
!--------------------------------------------------------
!
!--------------------------------------------------------
!--------------------------------------------------------
!--------------------------------------------------------
!--------------------------------------------------------
      SUBROUTINE FUN ( N , TIME , Z , ZP , RPAR , IPAR )
!
!
      IMPLICIT NONE
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
!
!
      INCLUDE 'param.common'
!
!
      DOUBLE PRECISION P   , RU   , Z , ZP  , RPAR ,
     1                 RHO , SUM  , H , CPB , WDOT ,
     2                 WT  , TIME
!
      INTEGER IPAR , KK , NWT , NH , NWDOT , K , N
!
      COMMON /RCONS/ P  , RU , RHO_
      COMMON /ICONS/ KK , NWT , NH , NWDOT
      DIMENSION Z(*) , ZP(*) , RPAR(*) , IPAR(*)
!
!
      CALL CKRHOY ( P , Z(1) , Z(2) , IPAR , RPAR , RHO )
!
      CALL CKCPBS ( Z(1) , Z(2) , IPAR , RPAR , CPB )
!
      CALL CKWYP  ( P , Z(1) , Z(2) , IPAR , RPAR , RPAR (NWDOT) )
!
      CALL CKHMS  ( Z(1) , IPAR , RPAR , RPAR (NH) )
!
!
      SUM = 0.0D0
!
      DO 100 K = 1 , KK
!
         H            = RPAR ( NH    + K - 1 )
         WDOT         = RPAR ( NWDOT + K - 1 )
         WOME ( K+1 ) = WDOT
         WT           = RPAR ( NWT   + K - 1 )
         ZP   ( K+1 ) = WDOT * WT / RHO
         SUM          = SUM + H * WDOT * WT
!
 100  CONTINUE
!
!
      ZP   (1) = - SUM / ( RHO * CPB )
      WOME (1) = - SUM / ( RHO * CPB )
!
!
!
      END SUBROUTINE FUN
!
!
!
      SUBROUTINE JAC
!
!
      IMPLICIT NONE
!
!
      END SUBROUTINE JAC
