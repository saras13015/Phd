      SUBROUTINE CKWYR_qi  (RHO, T, Y, ICKWRK, RCKWRK, WDOT,
     *                      MYINDEX , RKR , RKF , RROP , WSCAL , YASS )
C
C  START PROLOGUE
C
C  SUBROUTINE CKWYR  (RHO, T, Y, ICKWRK, RCKWRK, WDOT)
C     Returns the molar production rates of the species given the
C     mass density, temperature and mass fractions;  see Eq. (49).
C
C  INPUT
C     P      - Pressure. (for subroutine CKWYP)
C                   cgs units - gm/cm**2
C                   Data type - real scalar
C     RHO    - Density. (for subroutine CKWYR)
C                   cgs units - gm/cm**3
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension Y(*) at least  KK=6 (number of species)
C                   Species must be arranged in the order of:
C                   Y(1:6) = H2  O2  H2O  H  HO2  N2
C     ICKWRK - Dummy Array of integer workspace.
C                   Data type - integer array
C                   (Not used; simply to be Compatible with Chemkin)
C     RCKWRK - Dummy Array of real work space.
C                   Data type - real array
C                   (Not used; simply to be Compatible with Chemkin)
C     MYINDEX- Scalar index you want to calculate the Rate of progress
C                   Data type - integer
C
C  OUTPUT
C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C     RKF    - Dummy Forward molar rate of reaction of the species.
C                   In reduced chemistry, it's not possible to access 
C                   only to the forward molar rate.
C     RKR    - Dummy Reversed molar rate of reaction of the species.
C                   In reduced chemistry, it's not possible to access 
C                   only to the reversed molar rate.
C     RROP   - Rate of progress variable for the reactions.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension RROP(*) at least II, the total number of
C                   reactions.
C     WSCAL  - Rate of progress of the MYINDEX-th specie for the reactions.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WSCAL(*) at least II, the total number of
C                   reactions.
C     YASS   - Mass fraction of the species in steady-state.
C                   cgs units - g/g
C                   Data type - real array
C                   Dimension YASS(*) at least KSS=2, the total number 
C                   of species in steady-state.
C                   YASS(1:KSS) = OH O
C
C  END PROLOGUE
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C
      PARAMETER (KK=7,KSS=2,IREAC=12,II=4)

      DIMENSION ICKWRK(*), RCKWRK(*), Y(*), WDOT(*)
      DIMENSION XCON(KK+KSS), WT(KK+KSS)
      DIMENSION RF(IREAC), RB(IREAC), WF(IREAC), WB(IREAC)
      DIMENSION AF(IREAC), EXPTF(IREAC), EAF(IREAC)
      DIMENSION AB(IREAC), EXPTB(IREAC), EAB(IREAC)
      DIMENSION RKR(II),RKF(II),RROP(II),WSCAL(II)
      DIMENSION YASS(KSS)

!      INCLUDE 'myvar.h'                                                 ! include for SENKIN: if you want to extract SS_k and SLAMBDA

      DATA SMALL/1.0D-50/
      DATA TOLER/1.0D-20/
      DATA ITEMAX/10/ ! Good value, checked by P.Boivin and A.Techer
      DATA RU/8.314510D7/  ! ideal gas constant in (g*cm^2/s^-2)/mol/K
C-- Molar mass of the kth species (g/mol)
C-- 6+1+2 species: H2     O2      H2O     H    HO2     H2O2    (N2)     ! 4-Step
      DATA WT/  2.016, 31.999, 18.015, 1.008, 33.007, 34.015, 28.013, 
C                OH        O
     &         17.008, 15.999/

C-- Pre-exponential factor AForward and ABackward (cm^3, s, K, mol)
      DATA AF/  3.52D16, 5.06D04, 1.17D09, 5.75D19, 7.08D13, 1.66D13, 
     &          2.89D13, 4.00D22, 1.30D18, 3.02D12, 1.62D11, 8.15D23 /
      DATA AB/  7.04D13, 3.03D04, 1.28D10, 4.65D12, 0.00   , 2.69D12, 
     &          0.00   , 1.03D23, 3.04D17, 0.00   , 0.00   , 2.62D19 /

C-- Temperature exponent Forward and Backward (without unity)
      DATA EXPTF/  -0.7, 2.67, 1.3, -1.4, 0.0, 0.0, 0.0, -2.0, -1.0, 
     &              0.0, 0.61, -1.9/
      DATA EXPTB/  -0.26, 2.63, 1.19, 0.44, 0.0, 0.36, 0.0, -1.75, 
     &             -0.65, 0.0, 0.0, -1.39/

C-- Activation energy EaForward and EaBackward (kcal/mol)
      DATA EAF/  17.07, 6.291, 3.635,  0.00, 0.294, 0.822, -0.497, 0.00, 
     &            0.00,  1.39, 23.93, 49.62/
      DATA EAB/  0.143,  4.84, 18.70,  0.00, 0.00, 55.42, 0.00, 118.58, 
     &           103.5,  0.00,  0.00, 51.32/

C

C*** INPUT: PRESSURE P --> Section for routine compatible with CKWYP
!      SUMYOW = 0.0
!      DO K = 1, KK
!        SUMYOW = SUMYOW + Y(K)/WT(K)
!      END DO ! SUMYOW = 1/W (masse molaire du mélange en g/mol)
!      SUMYOW = SUMYOW*T*RU  ! SUMYOW = RT/W (en (g*cm^2/s^-2) /g)
!      DO K = 1, KK
!        XCON(K) = P*Y(K) / ( SUMYOW*WT(K) ) ! = P*W*Y_k / (R*T*W_k) = molar concentration (mol/cm^-3)
!      END DO
C*** END INPUT: PRESSURE P

C*** INPUT: DENSITY RHO --> Section for routine compatible with CKWYR
      DO K = 1, KK
        XCON(K) = RHO*Y(K) / WT(K) ! = rho*Y_k/W_k = molar concentration (mol/cm^-3)
      END DO
C*** END INPUT: DENSITY RHO

      SUMXCON = 0.0d0
      DO K=1, KK
         WDOT(K) = 0.0D0
         SUMXCON = SUMXCON + XCON (K)
      END DO
      DO I = 1, II
         WSCAL(I) = 0.0D0
      END DO
C
C     SET LOCAL VALUES FOR THE CONCENTRATIONS
C         H2  O2  H2O  H  X=HO2  N2
C
      CH2   =  MAX( XCON(1), SMALL )
      CO2   =  MAX( XCON(2), SMALL )
      CH2O  =  MAX( XCON(3), SMALL )
      CH    =  MAX( XCON(4), SMALL )
      CHO2  =  MAX( XCON(5), SMALL )
      CH2O2 =  MAX( XCON(6), SMALL )
C
C     EFFECTIVE THIRD BODY FOR ALL REACTIONS
C       3rd body concentration for the footnotes pages b, c and d
C       using the Chaperon efficiencies
C
      CMb = SUMXCON + 1.5d0*CH2 + 15.0d0*CH2O
      CMc = SUMXCON + 1.5d0*CH2 + 11.0d0*CH2O
      CMd = SUMXCON + 1.5d0*CH2 +  5.0d0*CH2O
C
      RUC= 1.987D0        ! ideal gas constant in cal/mol/K
      ALOGT=LOG(T)        ! ln(T)
      RTR=1.0D3/(RUC*T)   ! 10^3 / (cal/mol) ==> 10^3 because Ea is given in kcal
C
C
C     SET THE ELEMENTARY RATES
C       RF = k_forward = A_k * exp ( n_k * ln(T) - Ea / RT )
C       with A_k  = Frequency factor in cm, s, K et mol
C            n_k  = Temperature exponent (-)
C            Ea_k = Activation energy in kcal/mol --> ATTENTION in KILOcal !
C
      DO I = 1, IREAC
         RF(I) = AF(I) * EXP( EXPTF(I)*ALOGT - EAF(I)*RTR )
         RB(I) = AB(I) * EXP( EXPTB(I)*ALOGT - EAB(I)*RTR )
      END DO 
C
C
C     SET THE ELEMENTARY RATES INCLUDING THIRD BODIES: TROE FALLOFF CORRECTION
C
      RKinf = RB(4)
      RK0   = RF(4)
      PR = RK0 * CMb / RKinf
      PRLOG = DLOG10(MAX(PR,SMALL))
      FC = 0.5D0
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(4) = RKinf*PCOR

      RKinf = RB(12)
      RK0   = RF(12)
      PR = RK0 * CMd / RKinf
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     0.735D0
      F5 =     94.0D0
      F6 =     1756.0D0
      F7 =     5182.0D0
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(12) = RKinf*PCOR
C
C
C     ELEMENTARY FORWARD AND REVERSE RATE OF PROGRESS
      WF1_CH = RF(1)*CO2
      WF(1) = WF1_CH*CH
      WF2_CO = RF(2)*CH2
      WB2_COH = RB(2)*CH
      WF3_COH = RF(3)*CH2
      WB(3) = RB(3)*CH*CH2O
      WF4_CH = RF(4)*CO2
      WF(4) = WF4_CH*CH ! *CMb ! DON'T USE the 3rd body concentration !
      WF5_CHO2 = RF(5)*CH
      WF(5) = WF5_CHO2*CHO2
      WF(6) = RF(6)*CHO2*CH
      WB(6) = RB(6)*CH2*CO2
      WF7_COH = RF(7)*CHO2
      WF8_COH = RF(8)*CH*CMc
      WB(8) = RB(8)*CH2O*CMc
      WF(9) = RF(9)*CH*CH*CMc
      WB(9) = RB(9)*CH2*CMc
      WF(10) = RF(10)*CHO2*CHO2
      WF11_CHO2 = RF(11)*CH2
      WF(11) = WF11_CHO2*CHO2
      WF(12) = RF(12)*CH2O2


C         STEADY STATE OF OH (CF_2013,eq.4,6) 
C           => COH=f(CHO2)

!        IF ( CH > 1.0d-3 * SUMXCON ) then
          WRK1 = WB(3) + WF(5)+WF(5) + WF(12)+WF(12) + WB(8)
          WRK2 = WF3_COH + WF7_COH + WF8_COH

          A0 = WF2_CO * ( WF(1)+WF(1) + WRK1 )
          A1 = WF2_CO * WRK2 - RB(1) * WRK1
          A2 = RB(1) * ( WB2_COH+WB2_COH + WRK2 )                         ! already added EPSC to avoid divisions by zero
          AA = A1*A1 + 4.0D0*A0*A2                                        ! >= 0 by definition

          COH = ( DSQRT(AA) - A1 ) / ( A2+A2 )
          COH = MAX ( 0.0D0 , COH )
!        ELSE
!          COH = 0.0d0
!        END IF


C       STEADY STATE OF O (CF_2013,eq.5 = PCI_2011,eq.4) 
C         => CO=f(COH) but COH doesn't depend on CO !

!        IF ( COH > 1.0d-3 * SUMXCON .OR. CH2 > 1.0d-2 * SUMXCON ) THEN
          DENOM = RB(1)*COH + WF2_CO                                        ! already added EPSC to avoid divisions by zero

          CO = ( WF(1) + WB2_COH*COH ) / DENOM
          CO = MAX ( 0.0D0 , CO )
!        ELSE
!          CO = 0.0d0
!        END IF


C       COMPLETION OF THE WHOLE XCON VECTOR
        XCON(KK+1) = COH
        XCON(KK+2) = CO


C       UPD ELEMENTARY FORWARD AND REVERSE RATES OF PROGRESS
        WB(1) = RB(1)*COH*CO
        WF(7) = WF7_COH*COH
        WF(8) = WF8_COH*COH


C       GLOBAL FORWARD AND REVERSE REACTION RATES OF PROGRESS (CF_2013,eq.2) 

        RKF(1) = WF(1)+WF(5)+WF(12)
        RKR(1) = WB(1)

        RKF(2) = WF(4)+WF(8)+WF(9)
        RKR(2) = WB(8)+WB(9)+WF(10)+WF(11)

        RKF(3) = WF(4)+WB(6)
        RKR(3) = WF(5)+WF(6)+WF(7)+WF(10)+WF(10)+WF(11)

        RKF(4) = WF(10)+WF(11)
        RKR(4) = WF(12)


C       STEADY STATE PARAMETER FOR HO2
        RPROD_HO2 = RKF(3)
        RDEST_HO2 = RKR(3)
        SS_HO2 = CKSSparam ( RPROD_HO2 , RDEST_HO2 )


C       STEADY STATE PARAMETER FOR H
        RPROD_H = RKF(1) + RKF(1) + RKR(2) + RKR(2) + RKF(3)
        RDEST_H = RKR(1) + RKR(1) + RKF(2) + RKF(2) + RKR(3)
        SS_H = CKSSparam ( RPROD_H , RDEST_H )


        delta_corr = 0.03D0
        ss_corr    = 0.05D0


C       The correction FLAMBDA is to be applied only when SS_HO2 > 5%, and H2 is 
C       sufficiently close to steady state (that is necessary not to 
C       trigger the correction in very lean deflagrations)
C       SLAMBDA = small_lambda (CF_2012,eq.5 = PCI_2011,eq.6 + 4th_reaction)

        METHODE = 2
C       Meth 1: SLAMBDA from PCI_2011,eq.6
C       Meth 2: SLAMBDA from CF_2012,eq.5

!        IF ( METHODE == 1 ) WF4_CH = 0.0d0
        CALL SLAMBDA_COEFF (WF1_CH, WF2_CO, WF3_COH, WF4_CH, 
     *                      A0, A1, A2, AA)
        ! A2 /= 0 normally, it's zero if you have CO2 .AND CH2 = 0
        ! AA CAN BE < 0, if T<Tc. Cf. CF_2013 for correction

        CALL CK_SLAMBDA (AA < 0.0D0, WF1_CH, WF2_CO, WF3_COH, WF4_CH, 
     *                   A0, A1, A2, AA, SLAMBDA)
        CALL CK_SLAMBDA (SLAMBDA < 0.0D0, 
     *                   WF1_CH, WF2_CO, WF3_COH, WF4_CH, 
     *                   A0, A1, A2, AA, SLAMBDA)

!        IF ( SLAMBDA < 0.0d0 ) 
!     &   write(*,*) "negative lambda 2 = ", SLAMBDA
        SLAMBDA = max ( 0.0d0 , SLAMBDA )

!        SLAMBDA = ( WF1_CH + WF1_CH - WF4_CH )                            ! Decomment if you want to post-process lambda=f(phi) WITHOUT correction


!        IF ( SS_HO2 > ss_corr-delta_corr .AND. SS_H > 0.0D0 ) THEN       ! by P.Boivin
        IF ( SS_HO2 > ss_corr .AND. SS_H > ss_corr ) THEN                 ! modified by A.Techer
C       HO2 is out of SS, apply corrections for 

C         The correction for autoignition is needed to correct errors of 
C         the order of 40%. Then, there is no need to let FLAMBDA be below 0.25% 
C         (or 75% of error in autoignition time).
C         FLAMBDA = FACTOR LAMBDA (PCI_2011,eq.6 and using CF_2012,eq.5)

          FLAMBDA = SLAMBDA / ( WF1_CH + WF1_CH - WF4_CH )
          IF ( FLAMBDA < 0.0d0 ) write(*,*)'FLAMBDA < 0:' , FLAMBDA
          FLAMBDA = min ( 1.0d0 , FLAMBDA )
!          FLAMBDA = max ( 0.1d0 , FLAMBDA )              ! add by P.Boivin. A.Techer: I've never seen FLAMBDA exceed this bounds

C         Correction is only to be used for SS_HO2 > 5%,
C           =>  linear smoothing around ss_corr=5%
!          WRK1    = min ( 1.0d0 , 1.0d0 + (ss_corr-SS_HO2)/delta_corr )
!          FLAMBDA = max ( FLAMBDA , WRK1 )                               ! add by P.Boivin. A.Techer: I've never seen FLAMBDA exceed this bounds

C         CORRECTION FOR GLOBAL REACTION RATES OF PROGRESS (PCI_2011,eq.9)
          DO I = 1, II
             RKF(I) = FLAMBDA * RKF(I)
             RKR(I) = FLAMBDA * RKR(I)
          END DO

        END IF


C                                                                             Decomment if you want to post-process / extract X_kss
C     STEADY STATE MASS FRACTIONS
         SUMCW_i = 0.0D0
         DO K = 1, KK+KSS
           SUMCW_i = SUMCW_i + XCON(K) * WT(K)
         END DO
         SUMCW_i = 1.0D0 / SUMCW_i
         YASS(1) = COH * WT(KK+1) * SUMCW_i
         YASS(2) = CO  * WT(KK+2) * SUMCW_i
C     OR... STEADY STATE MOLAR FRACTIONS
!         SUMC_i = 0.0D0
!         DO K = 1, KK+KSS
!           SUMC_i = SUMC_i + XCON(K)
!         END DO
!         SUMC_i = 1.0D0 / SUMC_i
!         YASS(1) = COH * SUMC_i
!         YASS(2) = CO  * SUMC_i

!      XOHss = YASS(1)
!      XOss  = YASS(2)


C     GLOBAL RATE OF PROGRESS (CF_2013,eq.2)

      DO I = 1, II
         RROP(I) = RKF(I) - RKR(I)
      END DO


C     SPECIES PRODUCTION RATE for    H2  O2  H2O  H  HO2  H2O2  N2
C         3H2 + O2   = 2H2O + 2H    (R1)
C         H + H + M  = H2 + M       (R2)
C         H2 + O2    = HO2 + H      (R3)
C         H2 + O2    = H2O2         (R4)

      WDOT(1) = -RROP(1) -RROP(1) -RROP(1) +RROP(2) -RROP(3) -RROP(4)   ! H2
      WDOT(2) = -RROP(1) -RROP(3) -RROP(4)                              ! O2
      WDOT(3) = +RROP(1) +RROP(1)                                       ! H2O
      WDOT(4) = +RROP(1) +RROP(1) -RROP(2) -RROP(2) +RROP(3)            ! H
      WDOT(5) = +RROP(3)                                                ! HO2
      WDOT(6) = +RROP(4)                                                ! H2O2
!     WDOT(6) = 0.0D0                                                   ! N2


C     Rate of progress of the MYINDEX-th specie for the reactions.
      IF     ( MYINDEX .EQ. 1 ) THEN ! H2
         WSCAL(1) = -3.0D0*RROP(1)
         WSCAL(2) = +1.0D0*RROP(2)
         WSCAL(3) = -1.0D0*RROP(3)
         WSCAL(4) = -1.0D0*RROP(4)
      ELSEIF ( MYINDEX .EQ. 2 ) THEN ! O2
         WSCAL(1) = -1.0D0*RROP(1)
         WSCAL(3) = -1.0D0*RROP(3)
         WSCAL(4) = -1.0D0*RROP(4)
      ELSEIF ( MYINDEX .EQ. 3 ) THEN ! H2O
         WSCAL(1) = +2.0D0*RROP(1)
      ELSEIF ( MYINDEX .EQ. 4 ) THEN ! H
         WSCAL(1) = +2.0D0*RROP(1)
         WSCAL(2) = -2.0D0*RROP(2)
         WSCAL(3) = +1.0D0*RROP(3)
      ELSEIF ( MYINDEX .EQ. 5 ) THEN ! HO2
         WSCAL(3) = +1.0D0*RROP(3)
      ELSEIF ( MYINDEX .EQ. 6 ) THEN ! H2O2
         WSCAL(4) = +1.0D0*RROP(4)
      ELSEIF ( MYINDEX .EQ. 7 ) THEN ! N2
C        WSCAL(:) =  0.0D0
      ENDIF

      RETURN
      END
