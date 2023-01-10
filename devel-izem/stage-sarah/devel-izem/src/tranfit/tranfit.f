C///////////////////////////////////////////////////////////////////
C
      SUBROUTINE TRANFIT
C
C     WRITTEN BY:
C         ROBERT J. KEE
C         COMPUTATIONAL MECHANICS DIVISION
C         SANDIA NATIONAL LABORATORIES
C         LIVERMORE, CA  94550
C         (415) 294-3272
C
C/////////////////////////////////////////////////////////////////////
C
C     VERSION 3.0
C
C     CHANGES FROM PREVIOUS VERSION:
C     1.  Restructured data in OMEG1236
C     2.  Changed REAL*8 to DOUBLE PRECISION
C     3.  Changed POLFIT and PCOEF to call single and double precision
C         SLATEC subroutines POLFIT,DPOLFT and PCOEF,DPCOEF
C     4.  Change vms open statements
C     CHANGES FROM VERSION 1.3
C     1.  Change THIGH to 3500
C     2.  Change name to TRANFT to conform to ANSI standard
C     3.  Add "unix" and "snla" change blocks
C     CHANGES FROM VERSION 1.4
C     1.  modify OPEN statements for unix
C     CHANGES FOR VERSION 1.6
C     1.  Find THIGH and TLOW from species thermodynamic data
C     CHANGES FOR VERSION 1.9
C     1.  Change binary file to include version, precision, error
C         status, and required length of work arrays
C     2.  Allow user input from a file.
C     3.  Implement non-formatted transport data input
C     CHANGES FOR VERSION 3.0 (3/15/94 F. Rupley)
C     1.  DOS/PC compatibility effort includes adding file names to
C         OPEN statements, removing unused variables in CALL lists,
C         unusued but possibly initialized variables.
C
C/////////////////////////////////////////////////////////////////////
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (LINKCK=25, LINKTP=35, LTRAN=31, LDATA=10, 
     1           MDIM=50, KDIM=99, NFDIM=165, NO=4, NT=50, NCHECK=1000,
     2           LENICK=10000, LENRCK=10000, LENCCK=100, MAXFIT=7,
     3           NRANGE=2, MAXTMP=3)
C
      DIMENSION EPS(KDIM), SIG(KDIM), DIP(KDIM), POL(KDIM), ZROT(KDIM),
     1          NLIN(KDIM), WT(KDIM), CV(KDIM), FITWT(NT), FITRES(NT),
     2          ALOGT(NT), XLA(NT), XETA(NT), XD(NT), FITWRK(NFDIM),
     3          COFLAM(NO,KDIM), COFETA(NO,KDIM), COFD(NO,KDIM,KDIM),
     4          COFTD(NO,KDIM,10), ICKWRK(LENICK), RCKWRK(LENRCK),
     5          VALUE(6), KTDIF(KDIM), NTEMP(KDIM), TEMP(MAXTMP,KDIM),
     6          THERM(MAXFIT, NRANGE, KDIM) ,
     7          AWT(MDIM), KNCF(MDIM,KDIM)
C
      CHARACTER CCKWRK(LENCCK)*16, KSYM(KDIM)*16, LINE*80, VERS*16,
     1          PREC*16
      LOGICAL IERR, KERR, LEXIST
      COMMON /TPPARS/ PI, RU, PATM, BOLTZ, EPSIL, FDTCGS, FATCM,
     1                DIPMIN
C
      DATA NLIN/KDIM*NCHECK/, EMAXL,EMAXE,EMAXD,EMAXTD/4*0.0/,
     1     KERR/.FALSE./
C
C
      PI     = 3.1415926535
      EPSIL  = 0.0
      FDTCGS = 1.0E-18
      FATCM  = 1.0E+08
      DIPMIN = 1.0E-20
C
      CALL CKINTP36(KDIM2,AWT,KNCF)
c$$$      write (*,*) 'hola'
c$$$      write (*,*) ( awt (k) , k = 1,3)
c$$$
c$$$      do k =1,4
c$$$      write (*,*) kncf (1:4,k)
c$$$      end do
c$$$      read (*,*)
c$$$
c$$$      read (*,*)

      OPEN (LINKCK, FORM='UNFORMATTED', STATUS='UNKNOWN',
     1              FILE='chem.bin')
      OPEN (LTRAN,  FORM='FORMATTED', STATUS='UNKNOWN',
     1              FILE='tran.dat')

!      OPEN (LOUT, FORM='FORMATTED', STATUS='UNKNOWN', FILE='chem2.out')


C      OPEN (LINKTP, FORM='UNFORMATTED', STATUS='UNKNOWN',
C     1              FILE='tran.bin')
C      OPEN (LDATA,FORM='UNFORMATTED',STATUS='UNKNOWN',
C     1              FILE='tran.inp')
C
C          WRITE VERSION NUMBER
C
C      VERS = '3.0'
C      WRITE (LOUT, 15) VERS(:3)
C   15 FORMAT(
C     1/' TRANFT:  Transport property fitting,',
C     2/'           CHEMKIN-II Version ',A,', November 1990',
C*****precision > double
C     3/'           DOUBLE PRECISION')
C      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C     3/'           SINGLE PRECISION')
C      PREC = 'SINGLE'
C*****END precision > single
C
C      WRITE (LOUT, 8310)
C
C         INITIALIZE CHEMKIN
C
      CALL CKINIT36 (LENICK, LENRCK, LENCCK, LINKCK, LOUT, ICKWRK,
     1             RCKWRK, CCKWRK)
      CLOSE (LINKCK)
      CALL CKINDX36 (ICKWRK, RCKWRK, MM, KK, II, NFIT)
!      write(*,*) 'toto'
!      write(*,*) 'KK',KK,'KDIM',KDIM
      IF (KK .GT. KDIM) THEN
         WRITE (LOUT, 8300) KDIM
         STOP
      ENDIF
C
      CALL CKSYMS36 (CCKWRK, LOUT, KSYM, IERR)
      KERR = KERR.OR.IERR
      CALL CKWT36   (ICKWRK, RCKWRK, WT)
      CALL CKRP36   (ICKWRK, RCKWRK, RU, RUC, PATM)
      P = PATM
      BOLTZ = RU/6.02217E23
C
      CALL CKATHM36 (MAXFIT, NRANGE, ICKWRK, RCKWRK, MAXTMP,
     1               NTEMP, TEMP, THERM)
      TLOW  = TEMP(1,1)
      THIGH = TEMP(NTEMP(1),1)
      DO 10 K = 2, KK
         THIGH = MIN (THIGH, TEMP (NTEMP(K), K))
         TLOW  = MAX (TLOW,  TEMP (1,K))
   10 CONTINUE
      DT = (THIGH-TLOW) / (NT-1)
C
C      WRITE (LOUT, 20) TLOW, THIGH
C   20 FORMAT (/,' DATA HAS BEEN FIT OVER THE RANGE:  TLOW=',F9.2,
C     1        ', THIGH=',F9.2)
C
C       READ THE TRANSPORT DATA BASE
C
C      WRITE (LOUT, 7020)
      LIN = LTRAN
C
   50 CONTINUE
      READ (LIN, '(A)', END=500) LINE
      ILEN = IPPLEN_T36(LINE)
      IF (ILEN .GT. 0) THEN
         K1 = IFIRCH_T36(LINE(:ILEN))
         K2 = K1 + INDEX(LINE(K1:ILEN),' ') - 1
         CALL CKCOMP_T36 (LINE(K1:K2-1), KSYM, KK, KFIND)
C
         IF (KFIND .GT. 0) THEN
            CALL CKXNUM36 (LINE(K2:ILEN), 6, LOUT, NVAL, VALUE, IERR)
            KERR = KERR.OR.IERR
            IF (IERR) WRITE (LOUT, 8000)
C            WRITE (LOUT, 7010) LINE
C
            NLIN(KFIND) = INT(VALUE(1))
            EPS(KFIND)  = VALUE(2)
            SIG(KFIND)  = VALUE(3)
            DIP(KFIND)  = VALUE(4)
            POL(KFIND)  = VALUE(5)
            ZROT(KFIND) = VALUE(6)
         ENDIF
      ENDIF
      GO TO 50
C
  500 CONTINUE
C      IF (LEXIST .AND. LIN.EQ.LTRAN) THEN
C         LIN = LDATA
C         GO TO 50
C      ENDIF
C
      DO 600 K = 1, KK
         IF (NLIN(K) .EQ. NCHECK) THEN
            DO 750 J = K, KK
               IF (NLIN(J) .EQ. NCHECK) WRITE (LOUT, 8010) KSYM(J)
  750       CONTINUE
            KERR = .TRUE.
         ENDIF
  600 CONTINUE
C
      IF (.NOT. KERR) THEN
C
C        FIT THE CONDUCTIVITIES AND VISCOSITIES
C
         CALL LAMFIT36 (KK, NT, NO, LOUT, WT, SIG, EPS, DIP, ZROT, NLIN,
     1                  P, TLOW, DT, ALOGT, FITRES, FITWT, FITWRK, XLA,
     2                  XETA, CV, ICKWRK, RCKWRK, COFLAM, COFETA, EMAXL,
     3                  EMAXE)
C
C        FIT THE DIFFUSION COEFFICIENTS
C
         CALL DIFFIT36 (KK, NT, NO, KDIM, LOUT, WT, SIG, EPS, DIP, POL,
     1                  P, TLOW, DT, ALOGT, FITRES, FITWT, FITWRK, XD,
     2                  COFD, EMAXD)
C
C        FIT THE THERMAL DIFFUSION RATIOS
C
         CALL THMFIT36 (KK, NT, NO, KDIM, LOUT, WT, SIG, EPS, DIP, POL,
     1                  TLOW, DT, ALOGT, FITRES, FITWT, FITWRK, XD,
     2                  NLITE, KTDIF, COFTD, EMAXTD)
C
C        PRINT THE FITS
C
C         WRITE (LOUT, 7030)
C         WRITE (LOUT, 7035) EMAXL
C         WRITE (LOUT, 8200)
C         WRITE (LOUT, 8100) (KSYM(K), (COFLAM(N,K),N=1,NO), K=1,KK)
C
C         WRITE (LOUT, 7050)
C         WRITE (LOUT, 7035) EMAXE
C         WRITE (LOUT, 8200)
C         WRITE (LOUT, 8100) (KSYM(K), (COFETA(N,K),N=1,NO), K=1,KK)
C
C         WRITE (LOUT, 7060)
C         WRITE (LOUT, 7035) EMAXD
          LastData = -1
          MedData  = -2
          IFORMFILE = 1
          IF(IFORMFILE.EQ.0) THEN
C FORMER COMPACT VERSION NON COMPATIBLE WITH DEBUG ABSOFT
           OPEN(unit=20,file='tranfit.out')!,form='unformatted')
           WRITE (20,*) II
           WRITE (20,*) KK
           WRITE (20,*) NO
           WRITE (20,*) MM
           DO J = 1, KK
              WRITE (LOUT, 8200)
              WRITE (LOUT, 8110)
     1             (KSYM(J), KSYM(K), (COFD(N,J,K),N=1,NO), K=1,J)
               WRITE (20,*) J
               WRITE (20,*)
     1             (KSYM(J), KSYM(K), (COFD(N,J,K),N=1,NO), K=1,J)
           ENDDO
           WRITE (20,*) MedData
           WRITE (20,*) (KSYM(K), (COFLAM(N,K),N=1,NO), K=1,KK)
           WRITE (20,*) MedData
           WRITE (20,*) (KSYM(K), (COFETA(N,K),N=1,NO), K=1,KK)
           WRITE (20,*) MedData
           WRITE (20,*) (KSYM(K), WT(K), K=1,KK)
           WRITE (20,*) MedData
           WRITE (20,*) (AWT(K), K=1,MM)
           WRITE (20,*) MedData
           WRITE (20,*) ((KNCF(N,K), M=1,MM), K=1,KK)
           WRITE (20,*) LastData
           CLOSE(20)
          ELSE
C NEW LESS COMPACT VERSION
           OPEN(unit=20,file='tranfit.out')!,form='unformatted')
           WRITE (20,*) II
           WRITE (20,*) KK
           WRITE (20,*) NO
           WRITE (20,*) MM
           DO J = 1, KK
              WRITE (LOUT, 8200)
              WRITE (LOUT, 8110)
     1             (KSYM(J), KSYM(K), (COFD(N,J,K),N=1,NO), K=1,J)
              WRITE (20,*) J
              WRITE (20,*) KSYM(J)
              DO K = 1,J
                 WRITE (20,*) KSYM(K)
                 DO N = 1,NO
                    WRITE (20,*) COFD(N,J,K)
                 ENDDO
              ENDDO
           ENDDO
           WRITE (20,*) MedData
           DO K = 1,KK
              WRITE (20,*) KSYM(K)
              DO N = 1,NO
                 WRITE (20,*) COFLAM(N,K)
              ENDDO
           ENDDO
           WRITE (20,*) MedData
           DO K = 1,KK
              WRITE (20,*) KSYM(K)
              DO N = 1,NO
                 WRITE (20,*) COFETA(N,K)
              ENDDO
           ENDDO
           WRITE (20,*) MedData
           DO K = 1,KK
              WRITE (20,*) KSYM(K)
              WRITE (20,*) WT(K)
           ENDDO
           WRITE (20,*) MedData
           DO K = 1,MM
              WRITE (20,*) AWT(K) ! AWT(M) = RCKWRK(NcAW+M-1)
           ENDDO
           WRITE (20,*) MedData
           DO K = 1,KK
              DO N = 1,MM
                 WRITE (20,*) KNCF(N,K)
              END DO
           ENDDO
           WRITE (20,*) LastData
           CLOSE(20)
          ENDIF
C
C         WRITE (LOUT, 7070)
C         WRITE (LOUT, 7035) EMAXTD
C         DO 2400 M = 1, NLITE
C            K = KTDIF(M)
C            WRITE (LOUT, 8200)
C            WRITE (LOUT, 8110)
C     1            (KSYM(K), KSYM(J), (COFTD(N,J,M),N=1,NO), J=1,KK)
C 2400    CONTINUE
      ELSE
         WRITE (LOUT, '(/A)')
     1   ' WARNING...THERE IS AN ERROR IN THE TRANSPORT LINKING FILE'
      ENDIF
C
C        WRITE THE LINKING FILE
C
C      WRITE (LINKTP) VERS, PREC, KERR
C
C      LENIMC = 4*KK + NLITE
C      LENRMC = (19 + 2*NO + NO*NLITE)*KK + (15+NO)*KK**2
C
C      WRITE (LINKTP) LENIMC, LENRMC, NO, KK, NLITE
C      WRITE (LINKTP) PATM, (WT(K), EPS(K), SIG(K), DIP(K),
C     1               POL(K), ZROT(K), NLIN(K), K=1,KK),
C     2               ((COFLAM(N,K),N=1,NO),K=1,KK),
C     4               ((COFETA(N,K),N=1,NO),K=1,KK),
C     5               (((COFD(N,J,K),N=1,NO),J=1,KK),K=1,KK),
C     6               (KTDIF(N),N=1,NLITE),
C     7               (((COFTD(N,J,L),N=1,NO),J=1,KK),L=1,NLITE)
C      CLOSE (LINKTP)
C
!      STOP
C
C       FORMATS
C
 7000 FORMAT (A)
 7010 FORMAT (1X,A)
 7020 FORMAT (///,' TRANSPORT PARAMETERS FROM DATA BASE:',/)
 7030 FORMAT (///,' COEFFICIENTS FOR SPECIES CONDUCTIVITIES',/)
 7035 FORMAT ( '  MAXIMUM FITTING ERROR = ', 1PE12.3)
 7040 FORMAT (///,' COEFFICIENTS FOR MONOTOMIC PARTS OF CONDUCTIVITIES'
     1, /)
 7050 FORMAT (///,' COEFFICIENTS FOR SPECIES VISCOSITIES',/)
 7060 FORMAT (///,' COEFFICIENTS FOR SPECIES DIFFUSION COEFFICIENTS',
     1  /)
 7070 FORMAT (///,' COEFFICIENTS FOR THERMAL DIFFUSION RATIOS',/)
 8000 FORMAT (10X, 'ERROR IN CKXNUM READING FROM TRANSPORT DATA BASE')
 8010 FORMAT (10X, 'ERROR...TRANSPORT DATA NOT GIVEN FOR ',A10)
C
 8100 FORMAT (2X, A10, 4E12.3)
 8110 FORMAT (2X, A10, A10, 4E12.3)
 8200 FORMAT (1X, ' ')
 8300 FORMAT (10X, 'ERROR...THE TRANSPORT FITTING CODE IS DIMENSIONED
     1FOR ONLY', I3, ' SPECIES')
 8310 FORMAT(///,' OUTPUT FROM THE TRANSPORT PROPERTY FITTING PACKAGE'
     1,//)

      CLOSE (LTRAN)
!      CLOSE (LOUT)

      END
C
      SUBROUTINE LAMFIT36 (KK,NT,NO, NOUT, WT, SIG, EPS, DIP, ZROT,NLIN,
     1                   P, TLOW, DT, ALOGT, FITRES, FITWT, FITWRK, XLA,
     2                   XETA, CV, ICKWRK, RCKWRK, COFLAM, COFETA,EMAXL,
     3                   EMAXE)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WT(*), SIG(*), EPS(*), DIP(*), ZROT(*), NLIN(*), CV(*),
     1          FITRES(*), FITWRK(*), ALOGT(*), FITWT(*), XLA(*),
     2          XETA(*), ICKWRK(*), RCKWRK(*), COFLAM(NO,*),
     3          COFETA(NO,*)
C
      COMMON /TPPARS/ PI, RU, PATM, BOLTZ, EPSIL, FDTCGS, FATCM, DIPMIN
C
      FCON = 0.5 * FDTCGS**2 * FATCM**3 / BOLTZ
C
      PFAC = PATM / P
      DO 1000 K = 1, KK
         DST = FCON * DIP(K)**2 / (EPS(K) * SIG(K)**3)
         HELPE = 2.6693E-5 * SQRT(WT(K)) / SIG(K)**2
         HELPD = 2.6280E-3 * PFAC / (SQRT(WT(K)) * SIG(K)**2)
C
         DO 300 N = 1, NT
            T = TLOW + (N-1)*DT
            PRUT = P / (RU * T)
            TR = T / EPS(K)
            CALL CKCVML36 (T, ICKWRK, RCKWRK, CV)
            ALOGT(N) = LOG(T)
            DII = SQRT(T**3)*HELPD / OMEG1236(1,TR,DST)
            XETA(N) = SQRT(T)*HELPE / OMEG1236(2,TR,DST)
            RODET = DII * WT(K) * PRUT / XETA(N)
            AA = 2.5 - RODET
            CALL PARKER36 (T, EPS(K), PARK)
            IF (NLIN(K) .EQ. 2) THEN
               BB = ZROT(K)*PARK + 2.0*(5.0/2.0 + RODET)/PI
               CROCT = 1.0
            ELSE
               BB = ZROT(K)*PARK + 2.0*(5.0/3.0 + RODET)/PI
               CROCT = 2.0/3.0
            ENDIF
            FTRA = 2.5 * (1.0 - 2.0 * CROCT * AA / (BB*PI))
            FROT = RODET * (1.0 + 2.0 * AA / (BB*PI))
            FVIB = RODET
            IF (NLIN(K) .EQ. 0) THEN
               FAKTOR = 5.0/2.0 * 1.5*RU
            ELSEIF (NLIN(K) .EQ. 1) THEN
               FAKTOR = (FTRA*1.5 + FROT)*RU + FVIB*(CV(K)-2.5*RU)
            ELSEIF (NLIN(K) .EQ. 2) THEN
               FAKTOR = (FTRA+FROT)*1.5*RU + FVIB*(CV(K)-3.0*RU)
            ENDIF
            XLA(N)  = LOG( XETA(N)/WT(K) * FAKTOR)
            XETA(N) = LOG(XETA(N))
300      CONTINUE
C
C      FIT CONDUCTIVITY
C
         FITWT(1) = -1.0
         EPSL = EPSIL
C
C*****precision > double
         CALL DPOLFT36 (NT, ALOGT, XLA, FITWT, NO-1, NORD, EPSL, FITRES,
     1                IERR, FITWRK)
C*****END precision > double
C*****precision > single
C         CALL POLFIT (NT, ALOGT, XLA, FITWT, NO-1, NORD, EPSL, FITRES,
C     1                IERR, FITWRK)
C*****END precision > single
C
         EMAXL = MAX (EMAXL, EPSL)
         IF (IERR .NE. 1) THEN
            WRITE (NOUT,7000)
            STOP
         ENDIF
C
         CCC = 0.0
C*****precision > double
         CALL DPCOEF36 (NORD, CCC, COFLAM(1,K), FITWRK)
C*****END precision > double
C*****precision > single
C         CALL PCOEF (NORD, CCC, COFLAM(1,K), FITWRK)
C*****END precision > single
C
C      FIT VISCOSITY
C
         FITWT(1) = -1.0
         EPSL = EPSIL
C*****precision > double
         CALL DPOLFT36
C*****END precision > double
C*****precision > single
C         CALL POLFIT
C*****END precision > single
C
     1 (NT, ALOGT, XETA, FITWT, NO-1, NORD, EPSL, FITRES, IERR, FITWRK)
C
         EMAXE = MAX (EMAXE, EPSL)
         IF (IERR .NE. 1) THEN
            WRITE (NOUT,7000)
            STOP
         ENDIF
C
         CCC = 0.0
C*****precision > double
         CALL DPCOEF36 (NORD, CCC, COFETA(1,K), FITWRK)
C*****END precision > double
C*****precision > single
C         CALL PCOEF (NORD, CCC, COFETA(1,K), FITWRK)
C*****END precision > single
C
1000  CONTINUE
C
7000  FORMAT(9X,'ERROR IN POLFIT WHILE FITTING VISCOSITY OR ',
     1'CONDUCTIVITY')
      RETURN
      END
C
      SUBROUTINE DIFFIT36 (KK, NT, NO, KDIM, NOUT, WT, SIG, EPS,
     1                     DIP, POL, P, TLOW, DT, ALOGT, FITRES, FITWT,
     2                     FITWRK, XD, COFD, EMAXD)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WT(*), SIG(*), EPS(*), DIP(*), POL(*), FITRES(*),
     1          FITWT(*), FITWRK(*), ALOGT(*), XD(*), COFD(NO,KDIM,*)
C
      COMMON /TPPARS/ PI, RU, PATM, BOLTZ, EPSIL, FDTCGS, FATCM, DIPMIN
C
      DO 50 N = 1, NT
         T = TLOW + (N-1)*DT
         ALOGT(N) = LOG(T)
   50 CONTINUE
C
      FCON = 0.5 * FDTCGS**2 * FATCM**3 / BOLTZ
C
      DO 1000 J = 1, KK
         DO 1000 K = 1, J
            SIGM = 0.5 * (SIG(J)+SIG(K))
            EPSM = SQRT(EPS(J)*EPS(K))
            IF (DIP(J).LT.DIPMIN .AND. DIP(K).GT.DIPMIN) THEN
C
C        K IS POLAR, J IS NONPOLAR
C
               DST = 0.
               XI = 1.0 +
     1              POL(J) * FCON * DIP(K)**2 * SQRT(EPS(K)/EPS(J)) /
     2              ( 2.0 * SIG(J)**3 * EPS(K) * SIG(K)**3)
C
            ELSEIF (DIP(J).GT.DIPMIN .AND. DIP(K).LT.DIPMIN) THEN
C
C        J IS POLAR, K IS NONPOLAR
C
               DST = 0.
               XI = 1.0 +
     1              POL(K) * FCON * DIP(J)**2 * SQRT(EPS(J)/EPS(K)) /
     2              (2.0 * SIG(K)**3 * EPS(J) * SIG(J)**3)
C
C
            ELSE
C
C         NORMAL CASE, EITHER BOTH POLAR OR BOTH NONPOLAR
C
               DST = FCON * DIP(J) * DIP(K) / (EPSM * SIGM**3)
               XI = 1.0
            ENDIF
C
            SIGM = SIGM * XI**(-1.0/6.0)
            EPSM = EPSM * XI**2
            HELP1 = (WT(J)+WT(K)) / (2.0*WT(J)*WT(K))
            DO 500 N = 1, NT
               T = TLOW + (N-1)*DT
               TR = T / EPSM
               HELP2 = 0.002628 * SQRT(HELP1*T**3)
               XD(N) = HELP2*PATM / (OMEG1236(1,TR,DST) * SIGM**2 * P)
               XD(N) = XD(N) *
     1            D1236(WT(J),WT(K),T,EPS(J),EPS(K),SIG(J),SIG(K),DST)
               XD(N) = LOG(XD(N))
500         CONTINUE
C
            FITWT(1) = -1.0
            EPSL = EPSIL
C
C*****precision > double
            CALL DPOLFT36
C*****END precision > double
C*****precision > single
C            CALL POLFIT
C*****END precision > single
     1 (NT, ALOGT, XD, FITWT, NO-1, NORD, EPSL, FITRES, IERR, FITWRK)
C
            EMAXD = MAX (EMAXD, EPSL)
            IF (IERR .NE. 1) THEN
               WRITE (NOUT,7000) J,K
               STOP
            ENDIF
C
            CCC = 0.0
C
C*****precision > double
            CALL DPCOEF36 (NORD, CCC, COFD(1,K,J), FITWRK)
C*****END precision > double
C*****precision > single
C            CALL PCOEF (NORD, CCC, COFD(1,K,J), FITWRK)
C*****END precision > single
1000  CONTINUE
C
C          FILL OUT THE LOWER TRIANGLE
C
      DO 2000 K = 1, KK-1
         DO 2000 J = K+1, KK
            DO 2000 N = 1, NO
C MATTHIEU
C Coefficients de diffusion de l'argon = coefficients de diffusion du CH4
C Coefficient force dans oppdif.f mais pas premix !
C Dans oppdif les coefficients n etaient pas egaux malgre cette ligne:
               COFD(N,32,K) = COFD(N,1,K)
               COFD(N,J,K) = COFD(N,K,J)
 2000 CONTINUE
      RETURN
7000  FORMAT (10X, 'ERROR IN FITTING', I3, ',', I3,
     1                              'DIFFUSION COEFFICIENT')
      END
C
      SUBROUTINE THMFIT36 (KK, NT, NO, KDIM, NOUT, WT, SIG, EPS, DIP,
     1                     POL, TLOW, DT, AT, FITRES, FITWT, FITWRK, TD,
     2                     NLITE, KTDIF, COFTD, EMAXTD)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WT(*), SIG(*), EPS(*), DIP(*), POL(*), KTDIF(*),
     1          COFTD(NO,KDIM,*), FITRES(*), FITWT(*), FITWRK(*),
     2          AT(*), TD(*), FITAST(7), FITBST(7), FITCST(7)
C
      COMMON /TPPARS/ PI, RU, PATM, BOLTZ, EPSIL, FDTCGS, FATCM, DIPMIN
C
      DATA FITAST/ .1106910525E+01, -.7065517161E-02,-.1671975393E-01,
     1             .1188708609E-01,  .7569367323E-03,-.1313998345E-02,
     2             .1720853282E-03/
C
      DATA FITBST/ .1199673577E+01, -.1140928763E+00,-.2147636665E-02,
     1             .2512965407E-01, -.3030372973E-02,-.1445009039E-02,
     2             .2492954809E-03/
C
      DATA FITCST/ .8386993788E+00,  .4748325276E-01, .3250097527E-01,
     1            -.1625859588E-01, -.2260153363E-02, .1844922811E-02,
     2            -.2115417788E-03/
C
      NLITE = 0
      DO 50 N = 1, NT
         T = TLOW + (N-1)*DT
         AT(N) = T
   50 CONTINUE
C
      DO 1100 J = 1, KK
C
         IF (WT(J) .LE. 5.0) THEN
            NLITE = NLITE + 1
            KTDIF(NLITE) = J
            EPSJ = EPS(J) * BOLTZ
            SIGJ = SIG(J) * 1.0E-8
            POLJ = POL(J) * 1.0E-24
            POLJST = POLJ / SIGJ**3
C
            DO 1000 K = 1, KK
               EPSK = EPS(K) * BOLTZ
               SIGK = SIG(K) * 1.0E-8
               DIPK = DIP(K) * 1.0E-18
               DIPKST = DIPK / SQRT(EPSK*SIGK**3)
               EKOEJ = EPSK / EPSJ
               TSE = 1.0 + 0.25*POLJST*DIPKST**2*SQRT(EKOEJ)
               EOK = TSE**2 * SQRT(EPS(J)*EPS(K))
               WTKJ= (WT(K)-WT(J)) / (WT(K)+WT(J))
C
               DO 500 N = 1, NT
                  TSLOG = LOG(AT(N) / EOK)
                  CALL HORNER36 (7, FITAST, TSLOG, ASTAR)
                  CALL HORNER36 (7, FITBST, TSLOG, BSTAR)
                  CALL HORNER36 (7, FITCST, TSLOG, CSTAR)
                  TD(N) = 7.5 * WTKJ * (2.0*ASTAR + 5.0) *
     1                    (6.0*CSTAR - 5.0)/
     2               (ASTAR * (16.0*ASTAR - 12.0*BSTAR + 55.0))
  500          CONTINUE
C
               FITWT(1) = -1.0
               EPSL = EPSIL
C
C*****precision > double
               CALL DPOLFT36
C*****END precision > double
C*****precision > single
C               CALL POLFIT
C*****END precision > single
     1   (NT, AT, TD, FITWT, NO-1, NORD, EPSL, FITRES, IERR, FITWRK)
C
               EMAXTD = MAX (EMAXTD, EPSL)
               IF (IERR .NE. 1) THEN
                  WRITE (NOUT, 7000) J,K
                  STOP
               ENDIF
C
               CCC = 0.0
C
C*****precision > double
               CALL DPCOEF36 (NORD, CCC, COFTD(1,K,NLITE), FITWRK)
C*****END precision > double
C*****precision > single
C               CALL PCOEF (NORD, CCC, COFTD(1,K,NLITE), FITWRK)
C*****END precision > single
C
 1000       CONTINUE
         ENDIF
 1100 CONTINUE
 7000 FORMAT (10X, 'ERROR IN FITTING THE ', I3, ',', I3,
     1                ' THERMAL DIFFUSION RATIO')
      RETURN
      END
C
      SUBROUTINE HORNER36 (N, A, X, P)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION A(*)
      NM1 = N-1
      B = A(N)
      DO 10 I = 1, NM1
         B = A(N-I) + B*X
   10 CONTINUE
      P = B
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION D1236
     1                (W1, W2, T, EPS1, EPS2, SIG1, SIG2, DST)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      REAL FUNCTION D1236 (W1, W2, T, EPS1, EPS2, SIG1, SIG2, DST)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C        CORRECTION OF D(1,2).
C         REFERENCE: MARRERO AND MASON,J.PHYS CHEM REF DAT.1,3(1972)
C
      SUMW = W1+W2
      SIG12 = 0.5 * (SIG1 + SIG2)
      TR11 = T / EPS1
      TR22 = T / EPS2
      TR12 = T / SQRT(EPS1*EPS2)
C
      CALL OMEGXX36 (2, TR11, OM2F11, DST)
      CALL OMEGXX36 (2, TR22, OM2F22, DST)
      CALL OMEGXX36 (1, TR12, OM1F12, DST)
      CALL OMEGXX36 (2, TR12, OM2F12, DST)
      CALL OMEGXX36 (3, TR12, BST12,  DST)
      CALL OMEGXX36 (4, TR12, CST12,  DST)
      AST12 = OM2F12 / OM1F12
C
      H = (SIG1/SIG2)**2 * SQRT(2.0*W2/SUMW) * 2.0*W1**2/(SUMW*W2)
      P1 = H * OM2F11/OM1F12
C
      H = (SIG2/SIG1)**2 * SQRT(2.0*W1/SUMW) * 2.0*W2**2/(SUMW*W1)
      P2 = H * OM2F22/OM1F12
C
      P12 = 15.0 * ((W1-W2)/SUMW)**2 + 8.0*W1*W2*AST12/SUMW**2
C
      H = 2.0/(W2*SUMW) * SQRT(2.0*W2/SUMW) * OM2F11/OM1F12 *
     1    (SIG1/SIG12)**2
      Q1 = ((2.5-1.2*BST12)*W1**2 +
     1                           3.0*W2**2 + 1.6*W1*W2*AST12)*H
C
      H = 2.0/(W1*SUMW) * SQRT(2.0*W1/SUMW) * OM2F22/OM1F12 *
     1    (SIG2/SIG12)**2
      Q2 = ((2.5-1.2*BST12)*W2**2 +
     1                           3.0*W1**2 + 1.6*W1*W2*AST12)*H
C
      H = ((W1-W2)/SUMW)**2 * (2.5-1.2*BST12)*15.0 +
     1      4.0*W1*W2*AST12/SUMW**2 * (11.0 - 2.4*BST12)
      Q12 = H + 1.6*SUMW*OM2F11*OM2F22 / (SQRT(W1*W2)*OM1F12**2) *
     1      SIG1**2*SIG2**2/SIG12**4
      D1236 = ((6.0*CST12-5.0)**2/10.0) * (P1+P2+P12) /
     1      (Q1+Q2+Q12) + 1.0
C
      RETURN
      END
C
      SUBROUTINE PARKER36 (T, EPS, P)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C        TEMPERATURE DEPENDENCE OF THE ROTATIONAL COLLISION NUMBERS
C
C         REF: BRAU, C.A., JONKMAN, R.M., "CLASSICAL THEORY
C            OF ROTATIONAL RELAXATION IN DIATOMIC GASES",
C            JCP,VOL52,P477(1970)
C
      DATA PI32O2/2.7842/, P2O4P2/4.4674/, PI32/5.5683/
      D  = EPS / T
      DR = EPS / 298.0
      P = (1.0 + PI32O2*SQRT(DR) + P2O4P2*DR + PI32*DR**1.5) /
     1    (1.0 + PI32O2*SQRT(D)  + P2O4P2*D  + PI32*D**1.5)
      RETURN
      END
C
      SUBROUTINE INTERP36 (ARG, X, Y, VAL)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C//// QUADRATIC INTERPOLATION //////////////////////////////////////////
C
C
      DIMENSION X(3),Y(3)
      VAL1 = Y(1) + (ARG-X(1))*(Y(2)-Y(1)) / (X(2)-X(1))
      VAL2 = Y(2) + (ARG-X(2))*(Y(3)-Y(2)) / (X(3)-X(2))
      FAC1 = (ARG-X(1)) / (X(2)-X(1)) / 2.0
      FAC2 = (X(3)-ARG) / (X(3)-X(2)) / 2.0
      IF (ARG .GE. X(2)) THEN
         VAL = (VAL1*FAC2+VAL2) / (1.0+FAC2)
      ELSE
         VAL = (VAL1+VAL2*FAC1) / (1.0+FAC1)
      ENDIF
      RETURN
      END
C
C*****precision > double
      DOUBLE PRECISION FUNCTION OMEG1236 (N, TR, DR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      REAL FUNCTION OMEG1236 (N, TR, DR)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C//// CALC. OF COLLISION INTEGRALS FOR A KRIEGER-12-6-3-POTENTIAL //////
C
      DIMENSION VERT(3), ARG(3), VAL(3), TSTERN(37), DELTA(8), O(37,8),
     1          P(37,8)
      DATA TSTERN/.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,1.4,1.6,1.8,2.,2.5,
     1            3.,3.5,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,25.,
     2            30.,35.,40.,50.,75.,100./
      DATA DELTA/0.,.25,.5,.75,1.,1.5,2.,2.5/
      DATA ((O(I,L),L=1,8),I=1,10) /
     1   4.008 , 4.002 , 4.655 , 5.52  , 6.454 , 8.214 , 9.824 ,11.31  ,
     2   3.130 , 3.164 , 3.355 , 3.721 , 4.198 , 5.23  , 6.225 , 7.160 ,
     3   2.649 , 2.657 , 2.77  , 3.002 , 3.319 , 4.054 , 4.785 , 5.483 ,
     4   2.314 , 2.32  , 2.402 , 2.572 , 2.812 , 3.386 , 3.972 , 4.539 ,
     5   2.066 , 2.073 , 2.14  , 2.278 , 2.472 , 2.946 , 3.437 , 3.918 ,
     6   1.877 , 1.885 , 1.944 , 2.06  , 2.225 , 2.628 , 3.054 , 3.747 ,
     7   1.729 , 1.738 , 1.79  , 1.893 , 2.036 , 2.388 , 2.763 , 3.137 ,
     8   1.6122, 1.622 , 1.67  , 1.76  , 1.886 , 2.198 , 2.535 , 2.872 ,
     9   1.517 , 1.527 , 1.572 , 1.653 , 1.765 , 2.044 , 2.35  , 2.657 ,
     *   1.44  , 1.45  , 1.49  , 1.564 , 1.665 , 1.917 , 2.196 , 2.4780/
      DATA ((O(I,L),L=1,8),I=11,20) /
     1   1.3204, 1.33  , 1.364 , 1.425 , 1.51  , 1.72  , 1.956 , 2.199 ,
     2   1.234 , 1.24  , 1.272 , 1.324 , 1.394 , 1.573 , 1.777 , 1.99  ,
     3   1.168 , 1.176 , 1.202 , 1.246 , 1.306 , 1.46  , 1.64  , 1.827 ,
     4   1.1166, 1.124 , 1.146 , 1.185 , 1.237 , 1.372 , 1.53  , 1.7   ,
     5   1.075 , 1.082 , 1.102 , 1.135 , 1.181 , 1.3   , 1.441 , 1.592 ,
     6   1.0006, 1.005 , 1.02  , 1.046 , 1.08  , 1.17  , 1.278 , 1.397 ,
     7    .95  ,  .9538,  .9656,  .9852, 1.012 , 1.082 , 1.168 , 1.265 ,
     8    .9131,  .9162,  .9256,  .9413,  .9626, 1.019 , 1.09  , 1.17  ,
     9    .8845,  .8871,  .8948,  .9076,  .9252,  .972 , 1.03  , 1.098 ,
     *    .8428,  .8446,  .850 ,  .859 ,  .8716,  .9053,  .9483,  .9984/
      DATA ((O(I,L),L=1,8),I=21,30) /
     1    .813 ,  .8142,  .8183,  .825 ,  .8344,  .8598,  .8927,  .9316,
     2    .7898,  .791 ,  .794 ,  .7993,  .8066,  .8265,  .8526,  .8836,
     3    .7711,  .772 ,  .7745,  .7788,  .7846,  .8007,  .822 ,  .8474,
     4    .7555,  .7562,  .7584,  .7619,  .7667,  .78  ,  .7976,  .8189,
     5    .7422,  .743 ,  .7446,  .7475,  .7515,  .7627,  .7776,  .796 ,
     6    .72022, .7206,  .722 ,  .7241,  .7271,  .7354,  .7464,  .76  ,
     7    .7025,  .703 ,  .704 ,  .7055,  .7078,  .7142,  .7228,  .7334,
     8    .68776, .688,   .6888,  .6901,  .6919,  .697 ,  .704 ,  .7125,
     9    .6751,  .6753,  .676 ,  .677 ,  .6785,  .6827,  .6884,  .6955,
     *    .664 ,  .6642,  .6648,  .6657,  .6669,  .6704,  .6752,  .681 /
      DATA ((O(I,L),L=1,8),I=31,37) /
     1    .6414,  .6415,  .6418,  .6425,  .6433,  .6457,  .649 ,  .653 ,
     2    .6235,  .6236,  .6239,  .6243,  .6249,  .6267,  .629 ,  .632 ,
     3    .60882, .6089,  .6091,  .6094,  .61  ,  .6112,  .613 ,  .6154,
     4    .5964,  .5964,  .5966,  .597 ,  .5972,  .5983,  .600 ,  .6017,
     5    .5763,  .5763,  .5764,  .5766,  .5768,  .5775,  .5785,  .58  ,
     6    .5415,  .5415,  .5416,  .5416,  .5418,  .542 ,  .5424,  .543 ,
     7    .518 ,  .518 ,  .5182,  .5184,  .5184,  .5185,  .5186,  .5187/
C
      DATA ((P(I,L),L=1,8),I=1,10) /
     1   4.1   , 4.266 , 4.833 , 5.742 , 6.729 , 8.624 ,10.34  ,11.890 ,
     2   3.263 , 3.305 , 3.516 , 3.914 , 4.433 , 5.57  , 6.637 , 7.618 ,
     3   2.84  , 2.836 , 2.936 , 3.168 , 3.511 , 4.329 , 5.126 , 5.874 ,
     4   2.531 , 2.522 , 2.586 , 2.749 , 3.004 , 3.64  , 4.282 , 4.895 ,
     5   2.284 , 2.277 , 2.329 , 2.46  , 2.665 , 3.187 , 3.727 , 4.249 ,
     6   2.084 , 2.081 , 2.13  , 2.243 , 2.417 , 2.862 , 3.329 , 3.786 ,
     7   1.922 , 1.924 , 1.97  , 2.072 , 2.225 , 2.641 , 3.028 , 3.435 ,
     8   1.7902, 1.795 , 1.84  , 1.934 , 2.07  , 2.417 , 2.788 , 3.156 ,
     9   1.682 , 1.689 , 1.733 , 1.82  , 1.944 , 2.258 , 2.596 , 2.933 ,
     *   1.593 , 1.60  , 1.644 , 1.725 , 1.84  , 2.124 , 2.435 , 2.746 /
      DATA ((P(I,L),L=1,8),I=11,20) /
     1   1.455 , 1.465 , 1.504 , 1.574 , 1.67  , 1.913 , 2.181 , 2.45  ,
     2   1.355 , 1.365 , 1.4   , 1.461 , 1.544 , 1.754 , 1.989 , 2.228 ,
     3   1.28  , 1.289 , 1.321 , 1.374 , 1.447 , 1.63  , 1.838 , 2.053 ,
     4   1.222 , 1.231 , 1.26  , 1.306 , 1.37  , 1.532 , 1.718 , 1.912 ,
     5   1.176 , 1.184 , 1.209 , 1.25  , 1.307 , 1.45  , 1.618 , 1.795 ,
     6   1.0933, 1.1   , 1.119 , 1.15  , 1.193 , 1.304 , 1.435 , 1.578 ,
     7   1.039 , 1.044 , 1.06  , 1.083 , 1.117 , 1.204 , 1.31  , 1.428 ,
     8    .9996, 1.004 , 1.016 , 1.035 , 1.062 , 1.133 , 1.22  , 1.32  ,
     9    .9699,  .9732,  .983 ,  .9991, 1.021 , 1.08  , 1.153 , 1.236 ,
     *    .9268,  .9291,  .936 ,  .9473,  .9628, 1.005 , 1.058 , 1.12  /
      DATA ((P(I,L),L=1,8),I=21,30) /
     1    .8962,  .8979,  .903 ,  .9114,  .923 ,  .9545,  .9955, 1.044 ,
     2    .8727,  .8741,  .878 ,  .8845,  .8935,  .918 ,  .9505,  .9893,
     3    .8538,  .8549,  .858 ,  .8632,  .8703,  .890 ,  .9164,  .9482,
     4    .8379,  .8388,  .8414,  .8456,  .8515,  .868 ,  .8895,  .916 ,
     5    .8243,  .8251,  .8273,  .8308,  .8356,  .8493,  .8676,  .89  ,
     6    .8018,  .8024,  .8039,  .8065,  .810 ,  .820 ,  .8337,  .8504,
     7    .7836,  .784 ,  .7852,  .7872,  .7899,  .7976,  .808 ,  .8212,
     8    .7683,  .7687,  .7696,  .771 ,  .7733,  .7794,  .788 ,  .7983,
     9    .7552,  .7554,  .7562,  .7575,  .7592,  .764 ,  .771 ,  .7797,
     *    .7436,  .7438,  .7445,  .7455,  .747 ,  .7512,  .757 ,  .7642/
      DATA ((P(I,L),L=1,8),I=31,37) /
     1    .71982, .72  ,  .7204,  .7211,  .7221,  .725 ,  .7289,  .7339,
     2    .701 ,  .7011,  .7014,  .702 ,  .7026,  .7047,  .7076,  .7112,
     3    .68545, .6855,  .686 ,  .686 ,  .6867,  .6883,  .6905,  .693 ,
     4    .6723,  .6724,  .6726,  .673 ,  .6733,  .6745,  .676 ,  .6784,
     5    .651 ,  .651 ,  .6512,  .6513,  .6516,  .6524,  .6534,  .6546,
     6    .614 ,  .614 ,  .6143,  .6145,  .6147,  .6148,  .6148,  .6147,
     7    .5887,  .5889,  .5894,  .59  ,  .5903,  .5901,  .5895,  .5885/
C
      IF (DR.LT.-.00001 .OR. DR.GT.2.5 .OR. TR.LT..09 .OR. TR.GT.500.
     1 .OR. (ABS(DR).GT.1.0E-5 .AND. TR.GT.75.0)) THEN
         WRITE (6,'(A)') 'COLLISION INTEGRAL UNDEFINED'
         RETURN
      ENDIF
C
      IF (TR .GT. 75.0)  THEN
         IF (N .EQ. 1) THEN
            OMEG1236 = 0.623 - 0.136E-2*TR + 0.346E-5*TR*TR
     1                         -0.343E-8*TR*TR*TR
         ELSEIF (N .EQ. 2) THEN
            OMEG1236 = 0.703 - 0.146E-2*TR + 0.357E-5*TR*TR
     1                         -0.343E-8*TR*TR*TR
         ENDIF
         RETURN
      ENDIF
C
      IF (TR .LE. 0.2) THEN
         II = 2
      ELSE
         II = 37
         DO 30 I = 2, 36
            IF (TSTERN(I-1).LT.TR .AND. TSTERN(I).GE.TR)  II=I
   30    CONTINUE
      ENDIF
C
      IF (ABS(DR) .GE. 1.0E-5)  THEN
         IF (DR .LE. 0.25) THEN
            KK = 2
         ELSE
            KK = 7
            DO 10 K = 2, 7
               IF (DELTA(K-1).LT.DR .AND. DELTA(K) .GE. DR)  KK=K
   10       CONTINUE
         ENDIF
         DO 50 I = 1, 3
            DO 40 K = 1, 3
               ARG(K) = DELTA(KK-2+K)
               IF (N .EQ. 1)  THEN
                  VAL(K) = O(II-2+I, KK-2+K)
               ELSEIF (N .EQ. 2)  THEN
                  VAL(K) = P(II-2+I, KK-2+K)
               ENDIF
   40       CONTINUE
            CALL INTERP36 (DR, ARG, VAL, VERT(I))
   50    CONTINUE
         DO 60 I = 1, 3
            ARG(I) = TSTERN(II-2+I)
   60    CONTINUE
         CALL INTERP36 (TR, ARG, VERT, OMEG1236)
C
      ELSE
         DO 80 I = 1, 3
            ARG(I) = TSTERN(II-2+I)
            IF (N .EQ. 1) THEN
               VAL(I) = O(II-2+I, 1)
            ELSEIF (N .EQ. 2) THEN
               VAL(I) = P(II-2+I, 1)
            ENDIF
   80    CONTINUE
         CALL INTERP36 (TR, ARG, VAL, OMEG1236)
      ENDIF
      RETURN
      END
C
      SUBROUTINE OMEGXX36 (N, TR, OM, DST)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
C//// CALCULATION OF THE REDUCED COLLISION INTEGRALS ///////////////////
C
      IF (N.EQ.1 .OR. N.EQ.2) THEN
         OM = OMEG1236 (N, TR, DST)
      ELSEIF (N .EQ. 3) THEN
         IF (TR .LE. 5.0) THEN
            OM = 1.36 - 0.223*TR + 0.0613*TR**2 -0.00554*TR**3
         ELSE
            OM = 1.095
         ENDIF
      ELSEIF (N .EQ. 4) THEN
         IF (TR .LE. 5.0) THEN
            OM =0.823 + 0.0128*TR + 0.0112*TR**2 -0.00193*TR**3
         ELSE
            OM = .9483
         ENDIF
      ELSE
         WRITE (6,'(10X,A)') 'OMEGXX IS CALLED WITH WRONG PARAMETERS'
      ENDIF
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKCVML36 (T, ICKWRK, RCKWRK, CVML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCVML (T, ICKWRK, RCKWRK, CVML)
C     Returns the specific heats in constant volume in molar units;
C     see Eq. (22).
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     CVML   - Specific heats at constant volume in molar units
C              for the species.
C                   cgs units - ergs/(mole*K)
C                   Data type - real array
C                   Dimension CVML(*) at least KK, the total number of
C                   species.
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
      DIMENSION ICKWRK(*), RCKWRK(*), CVML(*)
C
      INCLUDE 'ckstrt36.common'
C
      CALL CKCPML36 (T, ICKWRK, RCKWRK, CVML)
C
      DO 150 K = 1, NKK
         CVML(K) = CVML(K) - RCKWRK(NcRU)
150   CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKINIT36 (LENIWK, LENRWK, LENCWK, LINC, LOUT, ICKWRK,
     1                     RCKWRK, CCKWRK)
C
C  A.Techer: This is the subroutine CKINIT called by the code. With 
C     also the subroutine CKINIT in chemkin/vode.f
C
C  START PROLOGUE
C
C  SUBROUTINE CKINIT (LENIWK, LENRWK, LENCWK, LINC, LOUT, ICKWRK,
C                     RCKWRK, CCKWRK)*
C     Reads the binary file and creates the internal work arrays
C     ICKWRK, CCKWRK, and RCKWORK.  CKINIT must be called before any
C     other CHEMKIN subroutine is called.  The work arrays must then
C     be made available as input to the other CHEMKIN subroutines.
C
C  INPUT
C     LENIWK - Length of the integer work array, ICKWRK.
C                   Data type - integer scalar
C     LENCWK - Length of the character work array, CCKWRK.
C              The minimum length of CCKWRK(*) is MM + KK.
C                   Data type - integer scalar
C     LENRWK - Length of the real work array, WORK.
C                   Data type - integer scalar
C     LINC  -  Logical file number for the binary file.
C                   Data type - integer scalar
C     LOUT  -  Output file for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C     CCKWRK - Array of character work space.
C                   Data type - CHARACTER*16 array
C                   Dimension CCKWRK(*) at least LENCWK.
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
      INTEGER LOUT
      DIMENSION ICKWRK(*), RCKWRK(*)
      CHARACTER CCKWRK(*)*(*), VERS*16, PREC*16
      LOGICAL IOK, ROK, COK, KERR
      COMMON /CKCONS/ PREC, VERS, KERR, LENI, LENR, LENC
C
      INCLUDE 'ckstrt36.common'
C
      COMMON /MACH/ SMALL,BIG,EXPARG
      COMMON /CMIN/ CKMIN
C
C     Data about the machine dependent constants is carried in
C
C     COMMON/MACH/SMALL,BIG,EXPARG
C
      DATA RU,RUC,PA /8.314E7, 1.987, 1.01325E6/
C
C      THIS STATEMENT WILL NOT COMPILE, MACHINE-DEPENDENT CONSTANTS
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
      EXPARG = LOG(BIG)
C
      LOUT=6 ! problem with gfortran/ifort
      WRITE (LOUT,15)
   15 FORMAT (/1X,' CKLIB:  Chemical Kinetics Library')
C*****END precision > double
C*****precision > single
C     2       /1X,'         SINGLE PRECISION')
C*****END precision > single
C
      CALL CKLEN36 (LINC, LOUT, LI, LR, LC)
C
      IOK = (LENIWK .GE. LI)
      ROK = (LENRWK .GE. LR)
      COK = (LENCWK .GE. LC)
      IF (.NOT.IOK .OR. .NOT.ROK .OR. .NOT.COK) THEN
         IF (.NOT. IOK) WRITE (LOUT, 300) LI
         IF (.NOT. ROK) WRITE (LOUT, 350) LR
         IF (.NOT. COK) WRITE (LOUT, 375) LC
         STOP
      ENDIF
C
      REWIND LINC
      READ (LINC, ERR=110) VERS, PREC, KERR
      READ (LINC, ERR=110) LENI, LENR, LENC, MM, KK, II,
     1                     MAXSP, MAXTB, MAXTP, NTHCF, NIPAR, NITAR,
     2                     NIFAR, NRV, NFL, NTB, NLT, NRL, NW, NCHRG,
     3                     NSTO, NOR, MAXORD, CKMN
C
      IF (LEN(CCKWRK(1)) .LT. 16) THEN
         WRITE (LOUT,475)
         STOP
      ENDIF
C
      NMM = MM
      NKK = KK
      NII = II
      MXSP = MAXSP
      MXTB = MAXTB
      MXTP = MAXTP
      MXOR = MAXORD
      NCP  = NTHCF
      NCP1 = NTHCF+1
      NCP2 = NTHCF+2
      NCP2T = NCP2*(MAXTP-1)
      NPAR = NIPAR
      NLAR = NITAR
      NFAR = NIFAR
      NTHB = NTB
      NLAN = NLT
      NFAL = NFL
      NREV = NRV
      NRLT = NRL
      NWL  = NW
      NRNU= NSTO
      NORD = NOR
      MXORD= MAXORD
      CKMIN= CKMN
C
C             APPORTION work arrays
C
C             SET  ICKWRK(*)=1  TO FLAG THAT CKINIT HAS BEEN CALLED
C
      ICKWRK(1) = 1
C
C             STARTING LOCATIONS OF INTEGER SPACE
C
C! elemental composition of species
      IcNC = 2
C! species phase array
      IcPH = IcNC + KK*MM
C! species charge array
      IcCH = IcPH + KK
C! # of temperatures for fit
      IcNT = IcCH + KK
C! stoichiometric coefficients
      IcNU = IcNT + KK
C! species numbers for the coefficients
      IcNK = IcNU + MAXSP*II
C! # of non-zero coefficients  (<0=reversible, >0=irreversible)
      IcNS = IcNK + MAXSP*II
C! # of reactants
      IcNR = IcNS + II
C! Landau-Teller reaction numbers
      IcLT = IcNR + II
C! Reverse Landau-Teller reactions
      IcRL = IcLT + NLAN
C! Fall-off reaction numbers
      IcFL = IcRL + NRLT
C! Fall-off option numbers
      IcFO = IcFL + NFAL
C! Fall-off enhanced species
      IcKF = IcFO + NFAL
C! Third-body reaction numbers
      IcTB = IcKF + NFAL
C! number of 3rd bodies for above
      IcKN = IcTB + NTHB
C! array of species #'s for above
      IcKT = IcKN + NTHB
C! Reverse parameter reaction numbers
      IcRV = IcKT + MAXTB*NTHB
C! Radiation wavelength reactions
      IcWL = IcRV + NREV
C! Real stoichometry reactions
      IcRNU= IcWL + NWL
C! Change of order reactions
      IcORD= IcRNU + NRNU
C! Species for which there is a change of order
      IcKOR= IcORD + NORD
C
      ITOT = IcKOR + NORD*MXORD - 1
C
C             STARTING LOCATIONS OF CHARACTER SPACE
C
C! start of element names
      IcMM = 1
C! start of species names
      IcKK = IcMM + MM
      ITOC = IcKK + KK - 1
C
C             STARTING LOCATIONS OF REAL SPACE
C
C! atomic weights
      NcAW = 1
C! molecular weights
      NcWT = NcAW + MM
C! temperature fit array for species
      NcTT = NcWT + KK
C! thermodynamic coefficients
      NcAA = NcTT + MAXTP*KK
C! Arrhenius coefficients (3)
      NcCO = NcAA + (MAXTP-1)*NCP2*KK
C! Reverse coefficients
      NcRV = NcCO + (NPAR+1)*II
C! Landau-Teller #'s for NLT reactions
      NcLT = NcRV + (NPAR+1)*NREV
C! Reverse Landau-Teller #'s
      NcRL = NcLT + NLAR*NLAN
C! Fall-off parameters for NFL reactions
      NcFL = NcRL + NLAR*NRLT
C! 3rd body coef'nts for NTHB reactions
      NcKT = NcFL + NFAR*NFAL
C! wavelength
      NcWL = NcKT + MAXTB*NTHB
C! real stoichometric coefficients
      NcRNU= NcWL + NWL
C! change of order for species/reactions
      NcKOR= NcRNU + NRNU*MXSP
C! universal gas constant
      NcRU = NcKOR + NORD*MXORD
C! universal gas constant in units
      NcRC = NcRU + 1
C! pressure of one atmosphere
      NcPA = NcRC + 1
C! intermediate temperature-dependent forward rates
      NcKF = NcPA + 1
C! intermediate temperature-dependent reverse rates
      NcKR = NcKF + II
C! internal work space of length kk
      NcK1 = NcKR + II
C!          'ditto'
      NcK2 = NcK1 + KK
C!          'ditto'
      NcK3 = NcK2 + KK
C!          'ditto'
      NcK4 = NcK3 + KK
      NcI1 = NcK4 + KK
      NcI2 = NcI1 + II
      NcI3 = NcI2 + II
      NcI4 = NcI3 + II
      NTOT = NcI4 + II - 1
C
C        SET UNIVERSAL CONSTANTS IN CGS UNITS
C
      RCKWRK(NcRU) = RU
      RCKWRK(NcRC) = RUC
      RCKWRK(NcPA) = PA
C
C!element names, !atomic weights
      READ (LINC,err=111) (CCKWRK(IcMM+M-1), RCKWRK(NcAW+M-1), M=1,MM)
C
C!species names, !composition, !phase, !charge, !molec weight,
C!# of fit temps, !array of temps, !fit coeff'nts
      READ (LINC,err=222) (CCKWRK(IcKK+K-1),
     1     (ICKWRK(IcNC+(K-1)*MM + M-1),M=1,MM),
     2     ICKWRK(IcPH+K-1),
     3     ICKWRK(IcCH+K-1),
     4     RCKWRK(NcWT+K-1),
     5     ICKWRK(IcNT+K-1),
     6     (RCKWRK(NcTT+(K-1)*MAXTP + L-1),L=1,MAXTP),
     7     ((RCKWRK(NcAA+(L-1)*NCP2+(K-1)*NCP2T+N-1),
     8     N=1,NCP2), L=1,(MAXTP-1)),    K = 1,KK)
C
      IF (II .EQ. 0) RETURN
C
C!# spec,reactants, !Arr. coefficients, !stoic coef, !species numbers
      READ (LINC,end=100,err=333)
     1     (ICKWRK(IcNS+I-1), ICKWRK(IcNR+I-1),
     2      (RCKWRK(NcCO+(I-1)*(NPAR+1)+N-1), N=1,NPAR),
     3      (ICKWRK(IcNU+(I-1)*MAXSP+N-1),
     4       ICKWRK(IcNK+(I-1)*MAXSP+N-1), N=1,MAXSP),
     5      I = 1,II)
C
C     PERTURBATION FACTOR
C
      DO 10 I = 1, II
         RCKWRK(NcCO + (I-1)*(NPAR+1) + NPAR) = 1.0
   10 CONTINUE
C
      IF (NREV .GT. 0) READ (LINC,err=444)
     1   (ICKWRK(IcRV+N-1), (RCKWRK(NcRV+(N-1)*(NPAR+1)+L-1),
     1   L=1,NPAR), N = 1,NREV)
C
      IF (NFAL .GT. 0) READ (LINC,err=555)
     1   (ICKWRK(IcFL+N-1), ICKWRK(IcFO+N-1), ICKWRK(IcKF+N-1),
     2   (RCKWRK(NcFL+(N-1)*NFAR+L-1),L=1,NFAR),N=1,NFAL)
C
      IF (NTHB .GT. 0) READ (LINC,err=666)
     1   (ICKWRK(IcTB+N-1), ICKWRK(IcKN+N-1),
     2   (ICKWRK(IcKT+(N-1)*MAXTB+L-1),
     3     RCKWRK(NcKT+(N-1)*MAXTB+L-1),L=1,MAXTB),N=1,NTHB)
C
      IF (NLAN .GT. 0) READ (LINC,err=777)
     1   (ICKWRK(IcLT+N-1), (RCKWRK(NcLT+(N-1)*NLAR+L-1),L=1,NLAR),
     2    N=1,NLAN)
C
      IF (NRLT .GT. 0) READ (LINC,err=888)
     1   (ICKWRK(IcRL+N-1), (RCKWRK(NcRL+(N-1)*NLAR+L-1),L=1,NLAR),
     2    N=1,NRLT)
C
      IF (NWL .GT. 0) READ (LINC,err=999)
     1   (ICKWRK(IcWL+N-1), RCKWRK(NcWL+N-1), N=1,NWL)
C
      IF (NRNU .GT. 0) READ (LINC,err=1111)
     1   (ICKWRK(IcRNU+N-1), (RCKWRK(NcRNU+(N-1)*MAXSP+L-1),L=1,MAXSP),
     2    N=1,NRNU)
C
      IF (NORD .GT. 0) READ (LINC,err=2222)
     1   (ICKWRK(IcORD+N-1), (ICKWRK(IcKOR+(N-1)*MXORD+L-1),
     2                        RCKWRK(NcKOR+(N-1)*MXORD+L-1),
     3   L = 1, MXORD), N=1, NORD)
C
  100 CONTINUE
      RETURN
C
  110 WRITE (LOUT,*) ' Error reading binary file...'
      STOP
  111 WRITE (LOUT,*) ' Error reading element data...'
      STOP
  222 WRITE (LOUT,*) ' Error reading species data...'
      STOP
  333 WRITE (LOUT,*) ' Error reading reaction data...'
      STOP
  444 WRITE (LOUT,*) ' Error reading reverse Arrhenius parameters...'
      STOP
  555 WRITE (LOUT,*) ' Error reading Fall-off data...'
      STOP
  666 WRITE (LOUT,*) ' Error reading third-body data...'
      STOP
  777 WRITE (LOUT,*) ' Error reading Landau-Teller data...'
      STOP
  888 WRITE (LOUT,*) ' Error reading reverse Landau-Teller data...'
      STOP
  999 WRITE (LOUT,*) ' Error reading Wavelength data...'
      STOP
 1111 WRITE (LOUT,*) ' Error reading real stoichometric data...'
      STOP
 2222 WRITE (LOUT,*) ' Error reading order data...'
      STOP
C
  300 FORMAT (10X,'ICKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  350 FORMAT (10X,'RCKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  375 FORMAT (10X,'CCKWRK MUST BE DIMENSIONED AT LEAST ',I5)
  475 FORMAT (10X,'CHARACTER LENGTH OF CCKWRK MUST BE AT LEAST 16 ')
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKLEN36 (LINC, LOUT, LI, LR, LC)
C
C  START PROLOGUE
C
C  SUBROUTINE CKLEN (LINC, LOUT, LENI, LENR, LENC)
C     Returns the lengths required for the work arrays.
C
C  INPUT
C
C     LINC  -  Logical file number for the binary file.
C                   Data type - integer scalar
C     LOUT  -  Output file for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     LENI  -  Minimum length required for the integer work array.
C                   Data type - integer scalar
C     LENR  -  Minimum length required for the real work array.
C                   Data type - integer scalar
C     LENC  -  Minimum length required for the character work array.
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
      PARAMETER (NLIST = 3)
      LOGICAL KERR, VOK, POK
      CHARACTER LIST(NLIST)*16, PREC*16, VERS*16
      COMMON /CKCONS/ PREC, VERS, KERR, LENI, LENR, LENC
C      DATA LIST/'1.9','2.0','2.1','2.2','2.3','2.4','2.5','2.6',
C     1          '2.7','2.8','2.9','3.0','3.1','3.2','3.3'/
      DATA LIST /'3.4','3.5','3.6'/
C
      VERS = ' '
      PREC = ' '
      LENI = 0
      LENR = 0
      LENC = 0
C
      KERR = .FALSE.
      REWIND LINC
      READ (LINC, ERR=999) VERS, PREC, KERR
C
      VOK = .FALSE.
      DO 5 N = 1, NLIST
         IF (VERS .EQ. LIST(N)) VOK = .TRUE.
    5 CONTINUE
C
      POK = .FALSE.
C*****precision > double
      IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.
C*****END precision > double
C*****precision > single
C      IF (INDEX(PREC, 'SING') .GT. 0) POK = .TRUE.
C*****END precision > single
C
      IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
         IF (KERR) THEN
            WRITE (LOUT,'(/A,/A)')
     1      ' There is an error in the Chemkin binary file...',
     2      ' Check CHEMKIN INTERPRETER output for error conditions.'
         ENDIF
         IF (.NOT. VOK) THEN
            WRITE (LOUT,'(/A,A)')
     1      ' Chemkin binary file is incompatible with Chemkin',
     2      ' Library Version 4.9'
         ENDIF
         IF (.NOT. POK) THEN
            WRITE (LOUT,'(/A,A)')
     1      ' Precision of Chemkin binary file does not agree with',
     2      ' precision of Chemkin library'
         ENDIF
         STOP
      ENDIF
C
      READ (LINC, ERR=999) LENICK, LENRCK, LENCCK, MM, KK, II,
     1                     MAXSP, MAXTB, MAXTP, NTHCF, NIPAR, NITAR,
     2                     NIFAR, NRV, NFL, NTB, NLT, NRL, NW, NCHRG,
     3                     NSTO, NOR, MAXORD, CKMN
      REWIND LINC
C
      LENI = LENICK
      LENR = LENRCK
      LENC = LENCCK
      LI   = LENI
      LR   = LENR
      LC   = LENC
      RETURN
C
  999 CONTINUE
      WRITE (LOUT, 50)
   50 FORMAT (' Error reading Chemkin binary file.')
      STOP
      END
C----------------------------------------------------------------------C
C
      SUBROUTINE CKINDX36 (ICKWRK, RCKWRK, MM, KK, II, NFIT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKINDX (ICKWRK, RCKWRK, MM, KK, II, NFIT)*
C     Returns a group of indices defining the size of the particular
C     reaction mechanism
C
C  INPUT
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     MM     - Total number of elements in mechanism.
C                   Data type - integer scalar
C     KK     - Total number of species in mechanism.
C                   Data type - integer scalar
C     II     - Total number of reactions in mechanism.
C                   Data type - integer scalar
C     NFIT   - number of coefficients in fits to thermodynamic data
C              for one temperature range; NFIT = number of
C              coefficients in polynomial fits to CP/R  +  2.
C                   Data type - integer scalar
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
      DIMENSION ICKWRK(*), RCKWRK(*)
C
      INCLUDE 'ckstrt36.common'
C
      MM = NMM
      KK = NKK
      II = NII
      NFIT = NCP2
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKSYMS36 (CCKWRK, LOUT, KNAME, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSYMS (CCKWRK, LOUT, KNAME, KERR)*
C     Returns the character strings of species names
C
C  INPUT
C     CCKWRK - Array of character work space.
C                   Data type - CHARACTER*16 array
C                   Dimension CCKWRK(*) at least LENCWK.
C     LOUT   - Output unit for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     KNAME  - Species names.
C                   Data type - CHARACTER*(*) array
C                   Dimension KNAME(*) at least KK,
C                   the total number of species.
C     KERR   - Error flag; character length errors will result in
C              KERR = .TRUE.
C                   Data type - logical
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
      CHARACTER*(*) CCKWRK(*), KNAME(*)
      LOGICAL KERR
C
      INCLUDE 'ckstrt36.common'
C
      KERR = .FALSE.
      ILEN = LEN(KNAME(1))
      DO 150 K = 1, NKK
         LT = ILASCH_T36(CCKWRK(IcKK + K - 1))
         KNAME(K) = ' '
         IF (LT .LE. ILEN) THEN
            KNAME(K) = CCKWRK(IcKK+K-1)
         ELSE
            WRITE (LOUT,*)
     1      ' Error in CKSYM...character string length too small '
            KERR = .TRUE.
         ENDIF
150   CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKWT36   (ICKWRK, RCKWRK, WT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWT   (ICKWRK, RCKWRK, WT)
C     Returns the molecular weights of the species
C
C  INPUT
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     WT     - Molecular weights of the species.
C                   cgs units - gm/mole
C                   Data type - real array
C                   Dimension WT(*) at least KK, the total number of
C                   species.
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
      DIMENSION ICKWRK(*), RCKWRK(*), WT(*)
C
      INCLUDE 'ckstrt36.common'
C
      DO 100 K = 1, NKK
         WT(K) = RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKRP36   (ICKWRK, RCKWRK, RU, RUC, PA)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRP   (ICKWRK, RCKWRK, RU, RUC, PA)
C     Returns universal gas constants and the pressure of one standard
C     atmosphere
C
C  INPUT
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     RU     - Universal gas constant.
C                   cgs units - 8.314E7 ergs/(mole*K)
C                   Data type - real scalar
C     RUC    - Universal gas constant used only in conjuction with
C              activation energy.
C                   preferred units - 1.987 cal/(mole*K)
C                   Data type - real scalar
C     PA     - Pressure of one standard atmosphere.
C                   cgs units - 1.01325E6 dynes/cm**2
C                   Data type - real scalar
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
      DIMENSION ICKWRK(*), RCKWRK(*)
C
      INCLUDE 'ckstrt36.common'
C
      RU  = RCKWRK(NcRU)
      RUC = RCKWRK(NcRC)
      PA  = RCKWRK(NcPA)
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKATHM36 (NDIM1, NDIM2, ICKWRK, RCKWRK, MAXTP, NT, TMP,
     1                     A)
C
C  START PROLOGUE
C
C  SUBROUTINE CKATHM (NDIM1, NDIM2, ICKWRK, RCKWRK, MAXTP, NT, TMP,
C                     A)
C     Returns the coefficients of the fits for thermodynamic properties
C     of the species; see Eqns. (19)-(21).
C
C  INPUT
C     NDIM1  - First dimension of the three-dimensional array of
C              thermodynamic fit coefficients, A; NDIM1 must be at
C              least NCP2, the total number of coefficients for one
C              temperature range.
C     NDIM2  - Second dimension of the three-dimensionalarray of
C              thermodynamic fit coefficients, A; NDIM2 must be at
C              least MXPT-1, the total number of temperature ranges.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     NT     - Number of temperatures used for fitting coefficients of
C              thermodynamic properties for the species.
C                   Data type - integer array
C                   Dimension NT(*) at least KK, the total number of
C                   species.
C     TMP    - Common temperatures dividing the thermodynamic fits for
C              the species.
C                   cgs units - K
C                   Data type - real array
C                   Dimension TMP(MAXT,*) exactly MAXT for the first
C                   dimension (the maximum number of temperatures
C                   allowed for a species) , and at least KK for the
C                   second dimension (the total number of species)
C     A      - Three dimensional array of fit coefficients to the
C              thermodynamic data for the species.
C              The indicies in  A(N,L,K) mean-
C              N = 1,NN are polynomial coefficients in CP/R
C                CP/R(K)=A(1,L,K) + A(2,L,K)*T + A(3,L,K)*T**2 + ...
C              N = NN+1 is a6 in Eq. (20)
C              N = NN+2 is a7 in Eq. (21)
C              L = 1..MXTP-1 is for each temperature range.
C              K  is  the  species index
C                   Data type - real array
C                   Dimension A(NPCP2,NDIM2,*) exactly NPCP2 and MXTP-1
C                   for the first and second dimensions and at least
C                   KK for the third.
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
      DIMENSION NT(*), TMP(MAXTP,*), A(NDIM1,NDIM2,*),
     1          ICKWRK(*), RCKWRK(*)
C
      INCLUDE 'ckstrt36.common'
C
      DO 100 K = 1, NKK
         NT(K) = ICKWRK(IcNT + K - 1)
  100 CONTINUE
C
      DO 140 L = 1, MXTP
         DO 140 K = 1, NKK
            TMP(L,K) = RCKWRK(NcTT + (K-1)*MXTP + L - 1)
  140 CONTINUE
C
      DO 150 K = 1, NKK
         DO 150 L = 1, MXTP-1
            NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
            DO 150 M = 1, NCP2
               A(M, L, K) = RCKWRK(NA1 + M - 1)
150   CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKCOMP_T36 (IST, IRAY, II, I)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCOMP_T (IST, IRAY, II, I)*
C     Returns the index of an element of a reference character
C     string array which corresponds to a character string;
C     leading and trailing blanks are ignored.
C
C
C  INPUT
C     IST   - A character string.
C                  Data type - CHARACTER*(*)
C     IRAY  - An array of character strings;
C                  Data type - CHARACTER*(*)
C                  Dimension IRAY(*) at least II
C     II    - The length of IRAY.
C                  Data type - integer scalar.
C
C  OUTPUT
C     I     - The first integer location in IRAY in which IST
C             corresponds to IRAY(I); if IST is not also an
C             entry in IRAY, I=0.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) IST, IRAY(*)
C
      I = 0
      DO 10 N = II, 1, -1
         IS1 = IFIRCH_T36(IST)
         IS2 = ILASCH_T36(IST)
         IR1 = IFIRCH_T36(IRAY(N))
         IR2 = ILASCH_T36(IRAY(N))
         IF ( IS2.GE.IS1 .AND. IS2.GT.0 .AND.
     1        IR2.GE.IR1 .AND. IR2.GT.0 .AND.
     2        IST(IS1:IS2).EQ.IRAY(N)(IR1:IR2) ) I = N
   10 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKXNUM36 (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKXNUM (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
C     This subroutine is called to parse a character string, LINE,
C     that is composed of several blank-delimited substrings.
C     Each substring is expected to represent a number, which
C     is converted to entries in the array of real numbers, RVAL(*).
C     NEXP is the number of values expected, and NVAL is the
C     number of values found.  This allows format-free input of
C     numerical data.  For example:
C
C     input:  LINE    = " 0.170E+14 0 47780.0"
C             NEXP    = 3, the number of values requested
C             LOUT    = 6, a logical unit number on which to write
C                       diagnostic messages.
C     output: NVAL    = 3, the number of values found
C             RVAL(*) = 1.700E+13, 0.000E+00, 4.778E+04
C             KERR    = .FALSE.
C
C  INPUT
C     LINE   - A character string.
C                   Data type - CHARACTER*80
C     NEXP   - Number of real values to be found in character string.
C              If NEXP is negative, then IABS(NEXP) values are
C              expected.  However, it is not an error condition,
C              if less values are found.
C                   Data type - integer scalar
C     LOUT   - Output unit for printed diagnostics.
C                   Data type - integer scalar
C
C  OUTPUT
C     NVAL   - Number of real values found in character string.
C                   Data type - integer scalar
C     RVAL   - Array of real values found.
C                   Data type - real array
C                   Dimension RVAL(*) at least NEXP
C     KERR   - Error flag;  syntax or dimensioning error results
C              in KERR = .TRUE.
C                   Data type - logical
C
C  END PROLOGUE
C
C     A '!' will comment out a line, or remainder of the line.
C
C*****precision > double
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*), ITEMP*80
      DIMENSION RVAL(*), RTEMP(80)
      LOGICAL KERR
C
C----------Find Comment String (! signifies comment)
C
      ILEN = IPPLEN_T36(LINE)
      NVAL = 0
      KERR = .FALSE.
C
      IF (ILEN .LE. 0) RETURN
      IF (ILEN .GT. 80) THEN
         WRITE (LOUT,*)     ' Error in CKXNUM...line length > 80 '
         WRITE (LOUT,'(A)') LINE
         KERR = .TRUE.
         RETURN
      ENDIF
C
      ITEMP = LINE(:ILEN)
      IF (NEXP .LT. 0) THEN
         CALL IPPAR_TR36 (ITEMP, -1, NEXP, RTEMP, NVAL, IERR, LOUT)
      ELSE
         CALL IPPAR_TR36 (ITEMP, -1, -NEXP, RTEMP, NVAL, IERR, LOUT)
         IF (IERR .EQ. 1) THEN
            WRITE (LOUT, *)    ' Syntax errors in CKXNUM...'
            WRITE (LOUT,'(A)') LINE
            KERR = .TRUE.
         ELSEIF (NVAL .NE. NEXP) THEN
            WRITE (LOUT,*) ' Error in CKXNUM...'
            WRITE (LOUT,'(A)') LINE
            KERR = .TRUE.
            WRITE (LOUT,*) NEXP,' values expected, ',
     1                     NVAL,' values found.'
         ENDIF
      ENDIF
      IF (NVAL .LE. IABS(NEXP)) THEN
         DO 20 N = 1, NVAL
            RVAL(N) = RTEMP(N)
   20    CONTINUE
      ENDIF
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKCPML36 (T, ICKWRK, RCKWRK, CPML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCPML (T, ICKWRK, RCKWRK, CPML)
C     Returns the specific heats at constant pressure in molar units
C
C  INPUT
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     CPML   - Specific heats at constant pressure in molar units
C              for the species.
C                   cgs units - ergs/(mole*K)
C                   Data type - real array
C                   Dimension CPML(*) at least KK, the total number of
C                   species.
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
      DIMENSION ICKWRK(*), RCKWRK(*), CPML(*), TN(10)
C
      INCLUDE 'ckstrt36.common'
C
      TN(1) = 1.0
      DO 150 N = 2, NCP
         TN(N) = T**(N-1)
150   CONTINUE
C
      DO 250 K = 1, NKK
         L = 1
         DO 220 N = 2, ICKWRK(IcNT + K - 1)-1
            TEMP = RCKWRK(NcTT + (K-1)*MXTP + N - 1)
            IF (T .GT. TEMP) L = L+1
 220     CONTINUE
C
         NA1 = NcAA + (L-1)*NCP2 + (K-1)*NCP2T
         CPML(K) = 0.0
         DO 250 N = 1, NCP
            CPML(K) = CPML(K) + RCKWRK(NcRU)*TN(N)*RCKWRK(NA1 + N - 1)
250   CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE IPPAR_TR36 (STRING,ICARD,NEXPEC,RVAL,NFOUND,IERR,LOUT)
C   BEGIN PROLOGUE  IPPAR_TR
C   REFER TO  IPGETR
C   DATE WRITTEN  850625   (YYMMDD)
C   REVISION DATE 851625   (YYMMDD)
C   CATEGORY NO.  J3.,J4.,M2.
C   KEYWORDS  PARSE
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Parses real variables from a character variable.  Called
C            by IPGETR, the IOPAK routine used for interactive input.
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  IPPAR_TR may be used for parsing an input record that contains real
C  values, but was read into a character variable instead of directly
C  into real variables.
C  The following benefits are gained by this approach:
C    - specification of only certain elements of the array is allowed,
C      thus letting the others retain default values
C    - variable numbers of values may be input in a record, up to a
C      specified maximum
C    - control remains with the calling program in case of an input
C      error
C    - diagnostics may be printed by IPPAR_TR to indicate the nature
C      of input errors
C
C   The contents of STRING on input indicate which elements of RVAL
C   are to be changed from their entry values, and values to which
C   they should be changed on exit.  Commas and blanks serve as
C   delimiters, but multiple blanks are treated as a single delimeter.
C   Thus, an input record such as:
C     '   1.,   2,,4.e-5   , ,6.e-6'
C   is interpreted as the following set of instructions by IPGETR:
C
C     (1) set RVAL(1) = 1.0
C     (2) set RVAL(2) = 2.0
C     (3) leave RVAL(3) unchanged
C     (4) set RVAL(4) = 4.0E-05
C     (5) leave RVAL(5) unchanged
C     (6) set RVAL(6) = 6.0E-06
C
C   IPPAR_TR will print diagnostics on the default output device, if
C   desired.
C
C   IPPAR_TR is part of IOPAK, and is written in ANSI FORTRAN 77
C
C   Examples:
C
C      Assume RVAL = (0., 0., 0.) and NEXPEC = 3 on entry:
C
C   input string           RVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '  2.34e-3,  3 45.1'    (2.34E-03, 3.0, 45.1)     0       3
C  '2,,3.-5'               (2.0, 0.0, 3.0E-05)       0       3
C  ',1.4,0.028E4'          (0.0, 1.4, 280.0)         0       3
C  '1.0, 2.a4, 3.0'        (1.0, 0.0, 0.0)           1       1
C  '1.0'                   (1.0, 0.0, 0.0)           2       1
C
C      Assume RVAL = (0.,0.,0.,0.) and NEXPEC = -4 on entry:
C
C   input string           RVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '1.,2.'                 (1.0, 2.0)                0       2
C  ',,3  4.0'              (0.0, 0.0, 3.0, 4.0)      0       4
C  '1,,3,,5.0'             (0.0, 0.0, 3.0, 0.0)      3       4
C
C  arguments: (I=input,O=output)
C  -----------------------------
C  STRING (I) - the character string to be parsed.
C
C  ICARD  (I) - data statement number, and error processing flag
C         < 0 : no error messages printed
C         = 0 : print error messages, but not ICARD
C         > 0 : print error messages, and ICARD
C
C  NEXPEC (I) - number of real variables expected to be input.  If
C         < 0, the number is unknown, and any number of values
C         between 0 and abs(nexpec) may be input.  (see NFOUND)
C
C  PROMPT (I) - prompting string, character type.  A question
C         mark will be added to form the prompt at the screen.
C
C  RVAL (I,O) - the real value or values to be modified.  On entry,
C       the values are printed as defaults.  The formal parameter
C       corresponding to RVAL must be dimensioned at least NEXPEC
C       in the calling program if NEXPEC > 1.
C
C  NFOUND (O) - the number of real values represented in STRING,
C         only in the case that there were as many or less than
C         NEXPEC.
C
C  IERR (O) - error flag:
C       = 0 if no errors found
C       = 1 syntax errors or illegal values found
C       = 2 for too few values found (NFOUND < NEXPEC)
C       = 3 for too many values found (NFOUND > NEXPEC)
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  IFIRCH_T36,ILASCH_T36
C   END PROLOGUE  IPPAR_TR
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STRING*(*), ITEMP*80
      DIMENSION RVAL(*)
      CHARACTER *8 FMT(16)
      LOGICAL OKINCR
C
C   FIRST EXECUTABLE STATEMENT  IPPAR_TR
      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = ILASCH_T36(STRING)
      IF (IE .EQ. 0) GO TO 500
      NC = 1
C
C--- OKINCR is a flag that indicates it's OK to increment
C--- NFOUND, the index of the array into which the value
C--- should be read.  It is set negative when a space follows
C--- a real value substring, to keep incrementing from
C--- occurring if a comma should be encountered before the
C--- next value.
C
      OKINCR = .TRUE.
C
C--- begin overall loop on characters in string
C
100   CONTINUE
C
      IF (STRING(NC:NC) .EQ. ',') THEN
         IF (OKINCR) THEN
            NFOUND = NFOUND + 1
         ELSE
            OKINCR = .TRUE.
         ENDIF
C
         GO TO 450
      ENDIF
      IF (STRING(NC:NC) .EQ. ' ') GO TO 450
C
C--- first good character (non-delimeter) found - now find
C--- last good character
C
      IBS = NC
160   CONTINUE
      NC = NC + 1
      IF (NC .GT. IE) GO TO 180
      IF (STRING(NC:NC) .EQ. ' ')THEN
         OKINCR = .FALSE.
      ELSEIF (STRING(NC:NC) .EQ. ',')THEN
         OKINCR = .TRUE.
      ELSE
         GO TO 160
      ENDIF
C
C--- end of substring found - read value into real array
C
180   CONTINUE
      NFOUND = NFOUND + 1
      IF (NFOUND .GT. NEXP) THEN
         IERR = 3
         GO TO 500
      ENDIF
C
      DATA FMT/     ' (E1.0)', ' (E2.0)', ' (E3.0)', ' (E4.0)',
     1   ' (E5.0)', ' (E6.0)', ' (E7.0)', ' (E8.0)', ' (E9.0)',
     2   '(E10.0)', '(E11.0)', '(E12.0)', '(E13.0)', '(E14.0)',
     3   '(E15.0)', '(E16.0)'/
      IES = NC - 1
      NCH = IES - IBS + 1
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
      READ (ITEMP(:NCH), FMT(NCH), ERR = 400) RVAL(NFOUND)
      GO TO 450
400   CONTINUE
      IERR = 1
      GO TO 510
450   CONTINUE
      NC = NC + 1
      IF (NC .LE. IE) GO TO 100
C
500   CONTINUE
      IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510   CONTINUE
C
      IF (IERR .EQ. 0 .OR. ICARD .LT. 0) RETURN
      IF (ICARD .NE. 0) WRITE(LOUT,'(A,I3)')
     1   '!! ERROR IN DATA STATEMENT NUMBER', ICARD
      IF (IERR .EQ. 1) WRITE(LOUT,'(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
      IF (IERR .EQ. 2) WRITE(LOUT,'(A,I2, A, I2)')
     1   ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND,
     2   '  NUMBER EXPECTED = ', NEXPEC
      IF (IERR .EQ. 3) WRITE(LOUT,'(A,I2)')
     1   ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
      END
C
C----------------------------------------------------------------------C
C
      FUNCTION IFIRCH_T36   (STRING)
C   BEGIN PROLOGUE  IFIRCH_T36
C   DATE WRITTEN   850626
C   REVISION DATE  850626
C   CATEGORY NO.  M4.
C   KEYWORDS  CHARACTER STRINGS,SIGNIFICANT CHARACTERS
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Determines first significant (non-blank) character
C            in character variable
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  IFIRCH_T36 locates the first non-blank character in a string of
C  arbitrary length.  If no characters are found, IFIRCH_T36 is set = 0.
C  When used with the companion routine ILASCH_T36, the length of a string
C  can be determined, and/or a concatenated substring containing the
C  significant characters produced.
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  (NONE)
C   END PROLOGUE IFIRCH_T36
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER* (*)STRING
C
C   FIRST EXECUTABLE STATEMENT IFIRCH_T36
      NLOOP = LEN(STRING)
C
      IF (NLOOP.EQ.0 .OR. STRING.EQ.' ') THEN
         IFIRCH_T36 = 0
         RETURN
      ENDIF
C
      DO 100 I = 1, NLOOP
         IF (STRING(I:I) .NE. ' ') GO TO 120
100   CONTINUE
C
      IFIRCH_T36 = 0
      RETURN
120   CONTINUE
      IFIRCH_T36 = I
      END
C
C----------------------------------------------------------------------C
C
      FUNCTION ILASCH_T36   (STRING)
C   BEGIN PROLOGUE  ILASCH_T36
C   DATE WRITTEN   850626
C   REVISION DATE  850626
C   CATEGORY NO.  M4.
C   KEYWORDS  CHARACTER STRINGS,SIGNIFICANT CHARACTERS
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Determines last significant (non-blank) character
C            in character variable
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  IFIRCH_T36 locates the last non-blank character in a string of
C  arbitrary length.  If no characters are found, ILASCH_T36 is set = 0.
C  When used with the companion routine IFIRCH_T36, the length of a string
C  can be determined, and/or a concatenated substring containing the
C  significant characters produced.
C  Note that the FORTRAN intrinsic function LEN returns the length
C  of a character string as declared, rather than as filled.  The
C  declared length includes leading and trailing blanks, and thus is
C  not useful in generating 'significant' substrings.
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  (NONE)
C   END PROLOGUE IFIRCH_T36
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) STRING
C
C   FIRST EXECUTABLE STATEMENT ILASCH_T36
      NLOOP = LEN(STRING)
      IF (NLOOP.EQ.0 .OR. STRING.EQ.' ') THEN
         ILASCH_T36 = 0
         RETURN
      ENDIF
C
      DO 100 I = NLOOP, 1, -1
         ILASCH_T36 = I
         IF (STRING(I:I) .NE. ' ') RETURN
100   CONTINUE
C
      END
C
C----------------------------------------------------------------------C
C
      FUNCTION IPPLEN_T36 (LINE)
C
C  BEGIN PROLOGUE
C
C  FUNCTION IPPLEN_T36 (LINE)
C     Returns the effective length of a character string, i.e.,
C     the index of the last character before an exclamation mark (!)
C     indicating a comment.
C
C  INPUT
C     LINE  - A character string.
C                  Data type - CHARACTER*(*)
C
C  OUTPUT
C     IPPLEN_T36 - The effective length of the character string.
C                   Data type - integer scalar
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*)
C
      IN = IFIRCH_T36(LINE)
      IF (IN.EQ.0 .OR. LINE(IN:IN) .EQ. '!') THEN
         IPPLEN_T36 = 0
      ELSE
         IN = INDEX(LINE,'!')
         IF (IN .EQ. 0) THEN
            IPPLEN_T36 = ILASCH_T36(LINE)
         ELSE
            IPPLEN_T36 = ILASCH_T36(LINE(:IN-1))
         ENDIF
      ENDIF
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE DPOLFT36 (N, X, Y, W, MAXDEG, NDEG, EPS, R, IERR, A)
C***BEGIN PROLOGUE  DPOLFT
C***PURPOSE  Fit discrete data in a least squares sense by polynomials
C            in one variable.
C***LIBRARY   SLATEC
C***CATEGORY  K1A1A2
C***TYPE      DOUBLE PRECISION (POLFIT-S, DPOLFT-D)
C***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
C***AUTHOR  Shampine, L. F., (SNLA)
C           Davenport, S. M., (SNLA)
C           Huddleston, R. E., (SNLL)
C***DESCRIPTION
C
C     Abstract
C
C     Given a collection of points X(I) and a set of values Y(I) which
C     correspond to some function or measurement at each of the X(I),
C     subroutine  DPOLFT  computes the weighted least-squares polynomial
C     fits of all degrees up to some degree either specified by the user
C     or determined by the routine.  The fits thus obtained are in
C     orthogonal polynomial form.  Subroutine  DP1VLU  may then be
C     called to evaluate the fitted polynomials and any of their
C     derivatives at any point.  The subroutine  DPCOEF  may be used to
C     express the polynomial fits as powers of (X-C) for any specified
C     point C.
C
C     The parameters for  DPOLFT  are
C
C     Input -- All TYPE REAL variables are DOUBLE PRECISION
C         N -      the number of data points.  The arrays X, Y and W
C                  must be dimensioned at least  N  (N .GE. 1).
C         X -      array of values of the independent variable.  These
C                  values may appear in any order and need not all be
C                  distinct.
C         Y -      array of corresponding function values.
C         W -      array of positive values to be used as weights.  If
C                  W(1) is negative,  DPOLFT  will set all the weights
C                  to 1.0, which means unweighted least squares error
C                  will be minimized.  To minimize relative error, the
C                  user should set the weights to:  W(I) = 1.0/Y(I)**2,
C                  I = 1,...,N .
C         MAXDEG - maximum degree to be allowed for polynomial fit.
C                  MAXDEG  may be any non-negative integer less than  N.
C                  Note -- MAXDEG  cannot be equal to  N-1  when a
C                  statistical test is to be used for degree selection,
C                  i.e., when input value of  EPS  is negative.
C         EPS -    specifies the criterion to be used in determining
C                  the degree of fit to be computed.
C                  (1)  If  EPS  is input negative,  DPOLFT  chooses the
C                       degree based on a statistical F test of
C                       significance.  One of three possible
C                       significance levels will be used:  .01, .05 or
C                       .10.  If  EPS=-1.0 , the routine will
C                       automatically select one of these levels based
C                       on the number of data points and the maximum
C                       degree to be considered.  If  EPS  is input as
C                       -.01, -.05, or -.10, a significance level of
C                       .01, .05, or .10, respectively, will be used.
C                  (2)  If  EPS  is set to 0.,  DPOLFT  computes the
C                       polynomials of degrees 0 through  MAXDEG .
C                  (3)  If  EPS  is input positive,  EPS  is the RMS
C                       error tolerance which must be satisfied by the
C                       fitted polynomial.  DPOLFT  will increase the
C                       degree of fit until this criterion is met or
C                       until the maximum degree is reached.
C
C     Output -- All TYPE REAL variables are DOUBLE PRECISION
C         NDEG -   degree of the highest degree fit computed.
C         EPS -    RMS error of the polynomial of degree  NDEG .
C         R -      vector of dimension at least NDEG containing values
C                  of the fit of degree  NDEG  at each of the  X(I) .
C                  Except when the statistical test is used, these
C                  values are more accurate than results from subroutine
C                  DP1VLU  normally are.
C         IERR -   error flag with the following possible values.
C             1 -- indicates normal execution, i.e., either
C                  (1)  the input value of  EPS  was negative, and the
C                       computed polynomial fit of degree  NDEG
C                       satisfies the specified F test, or
C                  (2)  the input value of  EPS  was 0., and the fits of
C                       all degrees up to  MAXDEG  are complete, or
C                  (3)  the input value of  EPS  was positive, and the
C                       polynomial of degree  NDEG  satisfies the RMS
C                       error requirement.
C             2 -- invalid input parameter.  At least one of the input
C                  parameters has an illegal value and must be corrected
C                  before  DPOLFT  can proceed.  Valid input results
C                  when the following restrictions are observed
C                       N .GE. 1
C                       0 .LE. MAXDEG .LE. N-1  for  EPS .GE. 0.
C                       0 .LE. MAXDEG .LE. N-2  for  EPS .LT. 0.
C                       W(1)=-1.0  or  W(I) .GT. 0., I=1,...,N .
C             3 -- cannot satisfy the RMS error requirement with a
C                  polynomial of degree no greater than  MAXDEG .  Best
C                  fit found is of degree  MAXDEG .
C             4 -- cannot satisfy the test for significance using
C                  current value of  MAXDEG .  Statistically, the
C                  best fit found is of order  NORD .  (In this case,
C                  NDEG will have one of the values:  MAXDEG-2,
C                  MAXDEG-1, or MAXDEG).  Using a higher value of
C                  MAXDEG  may result in passing the test.
C         A -      work and output array having at least 3N+3MAXDEG+3
C                  locations
C
C     Note - DPOLFT  calculates all fits of degrees up to and including
C            NDEG .  Any or all of these fits can be evaluated or
C            expressed as powers of (X-C) using  DP1VLU  and  DPCOEF
C            after just one call to  DPOLFT .
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  DP1VLU, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900911  Added variable YP to DOUBLE PRECISION declaration.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920527  Corrected erroneous statements in DESCRIPTION.  (WRB)
C***END PROLOGUE  DPOLFT
      INTEGER I,IDEGF,IERR,J,JP1,JPAS,K1,K1PJ,K2,K2PJ,K3,K3PI,K4,
     * K4PI,K5,K5PI,KSIG,M,MAXDEG,MOP1,NDEG,NDER,NFAIL
      DOUBLE PRECISION TEMD1,TEMD2
      DOUBLE PRECISION A(*),DEGF,DEN,EPS,ETST,F,FCRIT,R(*),SIG,SIGJ,
     * SIGJM1,SIGPAS,TEMP,X(*),XM,Y(*),YP,W(*),W1,W11
      DOUBLE PRECISION CO(4,3)
      SAVE CO
      DATA  CO(1,1), CO(2,1), CO(3,1), CO(4,1), CO(1,2), CO(2,2),
     1      CO(3,2), CO(4,2), CO(1,3), CO(2,3), CO(3,3),
     2  CO(4,3)/-13.086850D0,-2.4648165D0,-3.3846535D0,-1.2973162D0,
     3          -3.3381146D0,-1.7812271D0,-3.2578406D0,-1.6589279D0,
     4          -1.6282703D0,-1.3152745D0,-3.2640179D0,-1.9829776D0/
C***FIRST EXECUTABLE STATEMENT  DPOLFT
      M = ABS(N)
      IF (M .EQ. 0) GO TO 30
      IF (MAXDEG .LT. 0) GO TO 30
      A(1) = MAXDEG
      MOP1 = MAXDEG + 1
      IF (M .LT. MOP1) GO TO 30
      IF (EPS .LT. 0.0D0 .AND.  M .EQ. MOP1) GO TO 30
      XM = M
      ETST = EPS*EPS*XM
      IF (W(1) .LT. 0.0D0) GO TO 2
      DO 1 I = 1,M
        IF (W(I) .LE. 0.0D0) GO TO 30
 1      CONTINUE
      GO TO 4
 2    DO 3 I = 1,M
 3      W(I) = 1.0D0
 4    IF (EPS .GE. 0.0D0) GO TO 8
C
C DETERMINE SIGNIFICANCE LEVEL INDEX TO BE USED IN STATISTICAL TEST FOR
C CHOOSING DEGREE OF POLYNOMIAL FIT
C
      IF (EPS .GT. (-.55D0)) GO TO 5
      IDEGF = M - MAXDEG - 1
      KSIG = 1
      IF (IDEGF .LT. 10) KSIG = 2
      IF (IDEGF .LT. 5) KSIG = 3
      GO TO 8
 5    KSIG = 1
      IF (EPS .LT. (-.03D0)) KSIG = 2
      IF (EPS .LT. (-.07D0)) KSIG = 3
C
C INITIALIZE INDEXES AND COEFFICIENTS FOR FITTING
C
 8    K1 = MAXDEG + 1
      K2 = K1 + MAXDEG
      K3 = K2 + MAXDEG + 2
      K4 = K3 + M
      K5 = K4 + M
      DO 9 I = 2,K4
 9      A(I) = 0.0D0
      W11 = 0.0D0
      IF (N .LT. 0) GO TO 11
C
C UNCONSTRAINED CASE
C
      DO 10 I = 1,M
        K4PI = K4 + I
        A(K4PI) = 1.0D0
 10     W11 = W11 + W(I)
      GO TO 13
C
C CONSTRAINED CASE
C
 11   DO 12 I = 1,M
        K4PI = K4 + I
 12     W11 = W11 + W(I)*A(K4PI)**2
C
C COMPUTE FIT OF DEGREE ZERO
C
 13   TEMD1 = 0.0D0
      DO 14 I = 1,M
        K4PI = K4 + I
        TEMD1 = TEMD1 + W(I)*Y(I)*A(K4PI)
 14     CONTINUE
      TEMD1 = TEMD1/W11
      A(K2+1) = TEMD1
      SIGJ = 0.0D0
      DO 15 I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = TEMD1*A(K4PI)
        R(I) = TEMD2
        A(K5PI) = TEMD2 - R(I)
 15     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
      J = 0
C
C SEE IF POLYNOMIAL OF DEGREE 0 SATISFIES THE DEGREE SELECTION CRITERION
C
      IF (EPS) 24,26,27
C
C INCREMENT DEGREE
C
 16   J = J + 1
      JP1 = J + 1
      K1PJ = K1 + J
      K2PJ = K2 + J
      SIGJM1 = SIGJ
C
C COMPUTE NEW B COEFFICIENT EXCEPT WHEN J = 1
C
      IF (J .GT. 1) A(K1PJ) = W11/W1
C
C COMPUTE NEW A COEFFICIENT
C
      TEMD1 = 0.0D0
      DO 18 I = 1,M
        K4PI = K4 + I
        TEMD2 = A(K4PI)
        TEMD1 = TEMD1 + X(I)*W(I)*TEMD2*TEMD2
 18     CONTINUE
      A(JP1) = TEMD1/W11
C
C EVALUATE ORTHOGONAL POLYNOMIAL AT DATA POINTS
C
      W1 = W11
      W11 = 0.0D0
      DO 19 I = 1,M
        K3PI = K3 + I
        K4PI = K4 + I
        TEMP = A(K3PI)
        A(K3PI) = A(K4PI)
        A(K4PI) = (X(I)-A(JP1))*A(K3PI) - A(K1PJ)*TEMP
 19     W11 = W11 + W(I)*A(K4PI)**2
C
C GET NEW ORTHOGONAL POLYNOMIAL COEFFICIENT USING PARTIAL DOUBLE
C PRECISION
C
      TEMD1 = 0.0D0
      DO 20 I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = W(I)*((Y(I)-R(I))-A(K5PI))*A(K4PI)
 20     TEMD1 = TEMD1 + TEMD2
      TEMD1 = TEMD1/W11
      A(K2PJ+1) = TEMD1
C
C UPDATE POLYNOMIAL EVALUATIONS AT EACH OF THE DATA POINTS, AND
C ACCUMULATE SUM OF SQUARES OF ERRORS.  THE POLYNOMIAL EVALUATIONS ARE
C COMPUTED AND STORED IN EXTENDED PRECISION.  FOR THE I-TH DATA POINT,
C THE MOST SIGNIFICANT BITS ARE STORED IN  R(I) , AND THE LEAST
C SIGNIFICANT BITS ARE IN  A(K5PI) .
C
      SIGJ = 0.0D0
      DO 21 I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = R(I) + A(K5PI) + TEMD1*A(K4PI)
        R(I) = TEMD2
        A(K5PI) = TEMD2 - R(I)
 21     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
C
C SEE IF DEGREE SELECTION CRITERION HAS BEEN SATISFIED OR IF DEGREE
C MAXDEG  HAS BEEN REACHED
C
      IF (EPS) 23,26,27
C
C COMPUTE F STATISTICS  (INPUT EPS .LT. 0.)
C
 23   IF (SIGJ .EQ. 0.0D0) GO TO 29
      DEGF = M - J - 1
      DEN = (CO(4,KSIG)*DEGF + 1.0D0)*DEGF
      FCRIT = (((CO(3,KSIG)*DEGF) + CO(2,KSIG))*DEGF + CO(1,KSIG))/DEN
      FCRIT = FCRIT*FCRIT
      F = (SIGJM1 - SIGJ)*DEGF/SIGJ
      IF (F .LT. FCRIT) GO TO 25
C
C POLYNOMIAL OF DEGREE J SATISFIES F TEST
C
 24   SIGPAS = SIGJ
      JPAS = J
      NFAIL = 0
      IF (MAXDEG .EQ. J) GO TO 32
      GO TO 16
C
C POLYNOMIAL OF DEGREE J FAILS F TEST.  IF THERE HAVE BEEN THREE
C SUCCESSIVE FAILURES, A STATISTICALLY BEST DEGREE HAS BEEN FOUND.
C
 25   NFAIL = NFAIL + 1
      IF (NFAIL .GE. 3) GO TO 29
      IF (MAXDEG .EQ. J) GO TO 32
      GO TO 16
C
C RAISE THE DEGREE IF DEGREE  MAXDEG  HAS NOT YET BEEN REACHED  (INPUT
C EPS = 0.)
C
 26   IF (MAXDEG .EQ. J) GO TO 28
      GO TO 16
C
C SEE IF RMS ERROR CRITERION IS SATISFIED  (INPUT EPS .GT. 0.)
C
 27   IF (SIGJ .LE. ETST) GO TO 28
      IF (MAXDEG .EQ. J) GO TO 31
      GO TO 16
C
C RETURNS
C
 28   IERR = 1
      NDEG = J
      SIG = SIGJ
      GO TO 33
 29   IERR = 1
      NDEG = JPAS
      SIG = SIGPAS
      GO TO 33
 30   IERR = 2
C      CALL XERMSG ('SLATEC', 'DPOLFT', 'INVALID INPUT PARAMETER.', 2,
C     +   1)
      GO TO 37
 31   IERR = 3
      NDEG = MAXDEG
      SIG = SIGJ
      GO TO 33
 32   IERR = 4
      NDEG = JPAS
      SIG = SIGPAS
C
 33   A(K3) = NDEG
C
C WHEN STATISTICAL TEST HAS BEEN USED, EVALUATE THE BEST POLYNOMIAL AT
C ALL THE DATA POINTS IF  R  DOES NOT ALREADY CONTAIN THESE VALUES
C
      IF(EPS .GE. 0.0  .OR.  NDEG .EQ. MAXDEG) GO TO 36
      NDER = 0
      DO 35 I = 1,M
        CALL DP1VLU36 (NDEG,NDER,X(I),R(I),YP,A)
 35     CONTINUE
 36   EPS = SQRT(SIG/XM)
 37   RETURN
      END

      SUBROUTINE DPCOEF36 (L, C, TC, A)
C***BEGIN PROLOGUE  DPCOEF
C***PURPOSE  Convert the DPOLFT coefficients to Taylor series form.
C***LIBRARY   SLATEC
C***CATEGORY  K1A1A2
C***TYPE      DOUBLE PRECISION (PCOEF-S, DPCOEF-D)
C***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
C***AUTHOR  Shampine, L. F., (SNLA)
C           Davenport, S. M., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C     DPOLFT  computes the least squares polynomial fit of degree  L  as
C     a sum of orthogonal polynomials.  DPCOEF  changes this fit to its
C     Taylor expansion about any point  C , i.e. writes the polynomial
C     as a sum of powers of (X-C).  Taking  C=0.  gives the polynomial
C     in powers of X, but a suitable non-zero  C  often leads to
C     polynomials which are better scaled and more accurately evaluated.
C
C     The parameters for  DPCOEF  are
C
C     INPUT -- All TYPE REAL variables are DOUBLE PRECISION
C         L -      Indicates the degree of polynomial to be changed to
C                  its Taylor expansion.  To obtain the Taylor
C                  coefficients in reverse order, input  L  as the
C                  negative of the degree desired.  The absolute value
C                  of L  must be less than or equal to NDEG, the highest
C                  degree polynomial fitted by  DPOLFT .
C         C -      The point about which the Taylor expansion is to be
C                  made.
C         A -      Work and output array containing values from last
C                  call to  DPOLFT .
C
C     OUTPUT -- All TYPE REAL variables are DOUBLE PRECISION
C         TC -     Vector containing the first LL+1 Taylor coefficients
C                  where LL=ABS(L).  If  L.GT.0 , the coefficients are
C                  in the usual Taylor series order, i.e.
C                    P(X) = TC(1) + TC(2)*(X-C) + ... + TC(N+1)*(X-C)**N
C                  If L .LT. 0, the coefficients are in reverse order,
C                  i.e.
C                    P(X) = TC(1)*(X-C)**N + ... + TC(N)*(X-C) + TC(N+1)
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  DP1VLU
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPCOEF
C
      INTEGER I,L,LL,LLP1,LLP2,NEW,NR
      DOUBLE PRECISION A(*),C,FAC,SAVE,TC(*)
C***FIRST EXECUTABLE STATEMENT  DPCOEF
      LL = ABS(L)
      LLP1 = LL + 1
      CALL DP1VLU36 (LL,LL,C,TC(1),TC(2),A)
      IF (LL .LT. 2) GO TO 2
      FAC = 1.0D0
      DO 1 I = 3,LLP1
        FAC = FAC*(I-1)
 1      TC(I) = TC(I)/FAC
 2    IF (L .GE. 0) GO TO 4
      NR = LLP1/2
      LLP2 = LL + 2
      DO 3 I = 1,NR
        SAVE = TC(I)
        NEW = LLP2 - I
        TC(I) = TC(NEW)
 3      TC(NEW) = SAVE
 4    RETURN
      END

      SUBROUTINE DP1VLU36 (L, NDER, X, YFIT, YP, A)
C***BEGIN PROLOGUE  DP1VLU
C***PURPOSE  Use the coefficients generated by DPOLFT to evaluate the
C            polynomial fit of degree L, along with the first NDER of
C            its derivatives, at a specified point.
C***LIBRARY   SLATEC
C***CATEGORY  K6
C***TYPE      DOUBLE PRECISION (PVALUE-S, DP1VLU-D)
C***KEYWORDS  CURVE FITTING, LEAST SQUARES, POLYNOMIAL APPROXIMATION
C***AUTHOR  Shampine, L. F., (SNLA)
C           Davenport, S. M., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C     The subroutine  DP1VLU  uses the coefficients generated by  DPOLFT
C     to evaluate the polynomial fit of degree  L , along with the first
C     NDER  of its derivatives, at a specified point.  Computationally
C     stable recurrence relations are used to perform this task.
C
C     The parameters for  DP1VLU  are
C
C     Input -- ALL TYPE REAL variables are DOUBLE PRECISION
C         L -      the degree of polynomial to be evaluated.  L  may be
C                  any non-negative integer which is less than or equal
C                  to  NDEG , the highest degree polynomial provided
C                  by  DPOLFT .
C         NDER -   the number of derivatives to be evaluated.  NDER
C                  may be 0 or any positive value.  If NDER is less
C                  than 0, it will be treated as 0.
C         X -      the argument at which the polynomial and its
C                  derivatives are to be evaluated.
C         A -      work and output array containing values from last
C                  call to  DPOLFT .
C
C     Output -- ALL TYPE REAL variables are DOUBLE PRECISION
C         YFIT -   value of the fitting polynomial of degree  L  at  X
C         YP -     array containing the first through  NDER  derivatives
C                  of the polynomial of degree  L .  YP  must be
C                  dimensioned at least  NDER  in the calling program.
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DP1VLU
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER I,IC,ILO,IN,INP1,IUP,K1,K1I,K2,K3,K3P1,K3PN,K4,K4P1,K4PN,
     * KC,L,LM1,LP1,MAXORD,N,NDER,NDO,NDP1,NORD
      DOUBLE PRECISION A(*),CC,DIF,VAL,X,YFIT,YP(*)
      CHARACTER*8 XERN1, XERN2
C***FIRST EXECUTABLE STATEMENT  DP1VLU
C      IF (L .LT. 0) GO TO 12
      NDO = MAX(NDER,0)
      NDO = MIN(NDO,L)
      MAXORD = A(1) + 0.5D0
      K1 = MAXORD + 1
      K2 = K1 + MAXORD
      K3 = K2 + MAXORD + 2
      NORD = A(K3) + 0.5D0
      IF (L .GT. NORD) GO TO 11
      K4 = K3 + L + 1
      IF (NDER .LT. 1) GO TO 2
      DO 1 I = 1,NDER
 1      YP(I) = 0.0D0
 2    IF (L .GE. 2) GO TO 4
      IF (L .EQ. 1) GO TO 3
C
C L IS 0
C
      VAL = A(K2+1)
      GO TO 10
C
C L IS 1
C
 3    CC = A(K2+2)
      VAL = A(K2+1) + (X-A(2))*CC
      IF (NDER .GE. 1) YP(1) = CC
      GO TO 10
C
C L IS GREATER THAN 1
C
 4    NDP1 = NDO + 1
      K3P1 = K3 + 1
      K4P1 = K4 + 1
      LP1 = L + 1
      LM1 = L - 1
      ILO = K3 + 3
      IUP = K4 + NDP1
      DO 5 I = ILO,IUP
 5      A(I) = 0.0D0
      DIF = X - A(LP1)
      KC = K2 + LP1
      A(K4P1) = A(KC)
      A(K3P1) = A(KC-1) + DIF*A(K4P1)
      A(K3+2) = A(K4P1)
C
C EVALUATE RECURRENCE RELATIONS FOR FUNCTION VALUE AND DERIVATIVES
C
      DO 9 I = 1,LM1
        IN = L - I
        INP1 = IN + 1
        K1I = K1 + INP1
        IC = K2 + IN
        DIF = X - A(INP1)
        VAL = A(IC) + DIF*A(K3P1) - A(K1I)*A(K4P1)
        IF (NDO .LE. 0) GO TO 8
        DO 6 N = 1,NDO
          K3PN = K3P1 + N
          K4PN = K4P1 + N
 6        YP(N) = DIF*A(K3PN) + N*A(K3PN-1) - A(K1I)*A(K4PN)
C
C SAVE VALUES NEEDED FOR NEXT EVALUATION OF RECURRENCE RELATIONS
C
        DO 7 N = 1,NDO
          K3PN = K3P1 + N
          K4PN = K4P1 + N
          A(K4PN) = A(K3PN)
 7        A(K3PN) = YP(N)
 8      A(K4P1) = A(K3P1)
 9      A(K3P1) = VAL
C
C NORMAL RETURN OR ABORT DUE TO ERROR
C
 10   YFIT = VAL
      RETURN
C
   11 WRITE (XERN1, '(I8)') L
      WRITE (XERN2, '(I8)') NORD
C      CALL XERMSG ('SLATEC', 'DP1VLU',
C     *   'THE ORDER OF POLYNOMIAL EVALUATION, L = ' // XERN1 //
C     *   ' REQUESTED EXCEEDS THE HIGHEST ORDER FIT, NORD = ' // XERN2 //
C     *   ', COMPUTED BY DPOLFT -- EXECUTION TERMINATED.', 8, 2)
      RETURN
C
C   12 CALL XERMSG ('SLATEC', 'DP1VLU',
C     +   'INVALID INPUT PARAMETER.  ORDER OF POLYNOMIAL EVALUATION ' //
C     +   'REQUESTED IS NEGATIVE.', 2, 2)
      RETURN
      END
