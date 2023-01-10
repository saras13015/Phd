      SUBROUTINE PASR_REACTION_V (TF,ISTATE)
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
      EXTERNAL PASR_FUN , PASR_JAC
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
         CALL DVODE ( PASR_FUN , NEQ           , Z    , TI    , TF     ,
     *                ITOL     , RTOL          , ATOL , ITASK , ISTATE ,
     *                IOPT     , RWORK (NVODE) , LRW  , IWORK (IVODE)  ,
     *                LIW      , PASR_JAC      , MF   , RWORK  , IWORK )
!
C$$$        WRITE (*,*) 'boucle' , IBOUCLE
C$$$        WRITE (*,*) Z(1) , Z(2) , Z(3) , Z(4) , ISTATE , TI , TF
!
      END DO
!
!
!      IF ( IBOUCLE .GE. NBOUCLEMAX ) THEN
!         WRITE (*,'(1X,A,I5,3X,A,I5)') 'ISTATE = ' , ISTATE ,
!     *   'NO CHEMICAL SOLUTION IN PASR_REACTION_V FUNCTION, loop number= ' ,
!     *   IBOUCLE
!      END IF
!
!
!      IF ( IBOUCLE .GT. 1 ) THEN
!!         WRITE ( LOUT , * )
!         WRITE (*,*)
!     *   'CHEMICAL SOLUTION IN PASR_REACTION_V FUNCTION, loop number = ' , IBOUCLE
!      ENDIF
!
!
      IF ( ISTATE .LE. -2 ) THEN
         WRITE (*,'(A,I5)') 'U-PaSR: ISTATE = ' , ISTATE
         RETURN ! to stop the simulation with MPI_ABORT in CFD solver
      END IF
!
!
      END SUBROUTINE PASR_REACTION_V
!
!===============================================================================
!
      SUBROUTINE PASR_FUN ( N , TIME , Z , ZP , RPAR , IPAR )
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
!
      DOUBLE PRECISION HF , GAMMASTAR , TAUM
      COMMON /PASRCONS/ HF , GAMMASTAR , TAUM
      DIMENSION HF(110) ! chemkin.common (2*KMAX+10)
!
      DOUBLE PRECISION HKF(KK) , OMEGADOT(KK)
      DOUBLE PRECISION TAV , YKAV(KK) , HAV , HKAV(KK) , CPAV , CVAV
      DOUBLE PRECISION TS  , YKS(KK)  , HS  , HKS(KK)  , CPS  , CVS
      DOUBLE PRECISION DYKAVDT(KK) , DYKSDT(KK)
      DOUBLE PRECISION DHAVDT , DHSDT
      DOUBLE PRECISION DTAVDT , DTSDT
!
      DOUBLE PRECISION T_YKAV(2) , T_HAV(2) , T_TAV
      DOUBLE PRECISION T_YKS(2)  , T_HS(2)  , T_TS
!
      TAV = Z(1)     ! T: [K]
      TS  = Z(2+KK)  ! T*: [K]
      DO K = 1 , KK
         YKAV(K) = Z(1+K)     ! Yk: [-]
         YKS(K)  = Z(2+KK+K)  ! Yk*: [-]
         HKF(K)  = HF(K)      ! hkf: [erg / g]
      END DO
!
      CALL CKHMS  ( TAV , IPAR , RPAR , HKAV )         ! hk: [erg / g]
      CALL CKHBMS ( TAV , YKAV , IPAR , RPAR , HAV )   ! h: [erg / g]
      CALL CKCPBS ( TAV , YKAV , IPAR , RPAR , CPAV )  ! cp: [erg / (g * K)]
!      CALL CKCVBS ( TAV , YKAV , IPAR , RPAR , CVAV )  ! cv: [erg / (g * K)]
!
      CALL CKHMS  ( TS , IPAR , RPAR , HKS )        ! hk*: [erg / g]
      CALL CKHBMS ( TS , YKS , IPAR , RPAR , HS )   ! h*: [erg / g]
      CALL CKCPBS ( TS , YKS , IPAR , RPAR , CPS )  ! cp*: [erg / (g * K)]
!      CALL CKCVBS ( TS , YKS , IPAR , RPAR , CVS )  ! cv*: [erg / (g * K)]
!
      CALL CKWYR  ( RHO_ , TS , YKS , IPAR , RPAR , OMEGADOT )  ! wdot: [mol / (cm^3 * s)]
!
      MRU   = - RU                ! -R: [erg / (mol * K)]
      RHO_I = 1.0D0 / RHO_ * OMS  ! rho^{-1}: [cm^3 / g]
!
      T_HAV(1) = 0.0D0
      T_HS(1)  = 0.0D0
      T_TAV    = 0.0D0
      T_TS     = 0.0D0
      DO K = 1 , KK
!
         WT   = RPAR    ( NWT   + K - 1 )   ! MW: [g / mol]
         WDOT = OMEGADOT( K ) * WT * RHO_I  ! dYk/dt: [1 / s]
!
         T_YKAV(1) = GAMMASTAR * WDOT  ! dYk/dt: [1 / s]
         T_YKAV(2) = 0.0D0             ! [-]
!
         T_YKS(1) = WDOT   ! dYk*/dt: [1 / s]
         T_YKS(2) = 0.0D0  ! [-]
         IF ( ABS(1.0D0 - GAMMASTAR) >= EPS0 ) THEN  ! if U-PaSR
            T_YKS(2) = (YKS(K) - YKAV(K)) / (TAUM * (1.0D0 - GAMMASTAR))  ! U-PaSR term: [1 / s]
         END IF
!
         DYKAVDT(K) = T_YKAV(1) - T_YKAV(2)  ! dYk/dt: [1 / s]
         DYKSDT (K) = T_YKS(1)  - T_YKS(2)   ! dYk*/dt: [1 / s]
!
         ZP(1+K)    = DYKAVDT(K)  ! dYk/dt: [1 / s]
         ZP(2+KK+K) = DYKSDT(K)   ! dYk*/dt: [1 / s]
!
         T_HAV(1) = T_HAV(1) + ( GAMMASTAR * HKF(K) * WDOT )  ! [erg / (g * s)]
         T_HS(1)  = T_HS(1)  + ( HKF(K) * WDOT )              ! [erg / (g * s)]
!
         T_TAV = T_TAV + ( HKAV(K) * DYKAVDT(K) )          ! [erg / (g * s)]
!          T_TAV = ( T_TAV + HKAV(K) * OMEGADOT(K) * WT + 
!     *              MRU * TAV * OMEGADOT(K) ) * GAMMASTAR  ! [erg / (g * s)]
         T_TS  = T_TS  + ( HKS(K)  * DYKSDT(K)  )          ! [erg / (g * s)]
!
      END DO
!
      T_HAV(2) = 0.0D0  ! [-]
      T_HS(2)  = 0.0D0  ! [-]
      IF ( ABS(1.0D0 - GAMMASTAR) >= EPS0 ) THEN  ! if U-PaSR
         T_HS(2)  = (HS - HAV) / (TAUM * (1.0D0 - GAMMASTAR))  ! U-PaSR term: [erg / (g * s)]
      END IF
!
      DHAVDT = T_HAV(1) - T_HAV(2)  ! dh/dt: [erg / (g * s)]
      DHSDT  = T_HS(1)  - T_HS(2)   ! dh*/dt: [erg / (g * s)]
!
      DTAVDT = (DHAVDT - T_TAV) / CPAV  ! dT/dt: [K / s]
!      DTAVDT = - T_TAV / (RHO_ * CVAV)  ! dT/dt: [K / s]
      DTSDT  = (DHSDT  - T_TS)  / CPS   ! dT*/dt: [K / s]
!
      ZP(1)    = DTAVDT * OMS  ! dT/dt: [K / s]
      ZP(2+KK) = DTSDT  * OMS  ! dT*/dt: [K / s]
!
      END SUBROUTINE PASR_FUN
!
!===============================================================================
!
      SUBROUTINE PASR_JAC
!
!
      IMPLICIT NONE
!
!
      END SUBROUTINE PASR_JAC
