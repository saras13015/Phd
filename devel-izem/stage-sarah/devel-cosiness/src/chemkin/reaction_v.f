      SUBROUTINE REACTION_V (TF,ISTATE)
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
!
      DO WHILE ( ( ISTATE .NE. 2 ) .AND. ( IBOUCLE .LT. NBOUCLEMAX ) )
!
         ISTATE  = 1 ! OVERWRITING ISTATE=-1 TO AVOID STOP THE PROGRAM
         IBOUCLE = 1 + IBOUCLE
!
         CALL DVODE ( FUN  , NEQ           , Z    , TI    , TF     ,
     *                ITOL , RTOL          , ATOL , ITASK , ISTATE ,
     *                IOPT , RWORK (NVODE) , LRW  , IWORK (IVODE)  ,
     *                LIW  , JAC           , MF   , RWORK  , IWORK )
!
c$$$        WRITE (*,*) 'boucle' , IBOUCLE
c$$$        WRITE (*,*) Z(1) , Z(2) , Z(3) , Z(4) , ISTATE , TI , TF
!
      ENDDO
!
!
!      IF ( IBOUCLE .GE. NBOUCLEMAX ) THEN
!         WRITE (*,100) 'ISTATE = ' , ISTATE ,
!     *   'PAS DE SOLUTION CHIMIQUE dans REACTION_V.F, Nbre boucle= ' ,
!     *   IBOUCLE
!  100 FORMAT(1X,A,I5,3X,A,I5)
!      ENDIF
!
!
!      IF ( IBOUCLE .GT. 1 ) THEN
!!         WRITE ( LOUT , * )
!!     *   'SOLUTION CHIMIQUE dans REACTION_V.F, Nbre boucle = ' , IBOUCLE
!         WRITE (*,*)
!     *   'SOLUTION CHIMIQUE dans REACTION_V.F, Nbre boucle = ' , IBOUCLE
!      ENDIF
!
!
      IF ( ISTATE .LE. -2 ) THEN
         WRITE (*,'(A,I5)') 'ISTATE = ' , ISTATE
         RETURN ! to stop the simulation with MPI_ABORT in CFD solver
      ENDIF
!
!
      END SUBROUTINE REACTION_V
!
!--------------------------------------------------------
!--------------------------------------------------------
!--------------------------------------------------------
!--------------------------------------------------------
!
      SUBROUTINE FUN ( N , TIME , Z , ZP , RPAR , IPAR )
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
     3                 OMS
!
      INTEGER IPAR , KK , NWT , NH , NWDOT , K , N
!
      COMMON /RCONS/ P , RU , RHO_ , OMS
      COMMON /ICONS/ KK , NWT , NH , NWDOT
      DIMENSION Z(*) , ZP(*) , RPAR(*) , IPAR(*)
!
      DOUBLE PRECISION RHO_i , MRU
!
!
      CALL CKCVBS ( Z(1) , Z(2) , IPAR , RPAR , CVB )
      CALL CKWYR  ( RHO_ , Z(1) , Z(2) , IPAR , RPAR , RPAR (NWDOT) )
      CALL CKHMS  ( Z(1) , IPAR , RPAR , RPAR (NH) )
!
      RHO_I = 1.0D0 / RHO_ * OMS
      MRU   = - RU
!
      SUM = 0.0D0
      DO 100 K = 1 , KK
!
         H            = RPAR ( NH    + K - 1 )
         WDOT         = RPAR ( NWDOT + K - 1 )
         WT           = RPAR ( NWT   + K - 1 )
         ZP   ( K+1 ) = WDOT * WT * RHO_I
         SUM          = SUM + H * WDOT * WT + MRU * Z(1) * WDOT
!
 100  CONTINUE
!
!
      ZP   (1) = - SUM / ( RHO_ * CVB ) * OMS
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
