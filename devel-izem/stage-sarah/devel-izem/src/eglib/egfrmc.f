      SUBROUTINE EGFRMC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*16 CDUMMY
      LOGICAL LDUMMY
C-----------------------------------------------------------------------
      COMMON /ARAYS/ W(20000), IW(10000)
C-----------------------------------------------------------------------
C
C     Test des sous programmes de eglib.f
C
C-----------------------------------------------------------------------
      LLTP   = 10
      LLMC   = 10
C-----------------------------------------------------------------------
      LENWK  = 20000
      LENIWK = 10000
C-----------------------------------------------------------------------
C     LENGTH OF THE TOTAL WORK SPACE
C-----------------------------------------------------------------------
C------------ READ INFORMATION FROM A CHEMKIN-I LINKTP FILE ------------
C------------ OR AN OLD CHEMKIN-II LINKMC FILE -------------------------
C      OPEN(UNIT=LLMC,STATUS='OLD',FORM='UNFORMATTED',FILE='Linktp')
C         READ (LLMC) NO, NS, IDUMMY
C      CLOSE(UNIT=LLMC)
C------------ READ INFORMATION FROM A CHEMKIN LINKMC FILE --------------
      OPEN(UNIT=LLMC,STATUS='OLD',FORM='UNFORMATTED',FILE='tran.bin')
         READ (LLMC) CDUMMY, CDUMMY, LDUMMY
         READ (LLMC) IDUMMY, IDUMMY, NO, NS, IDUMMY
      CLOSE(UNIT=LLMC)
C------------- LENGTH OF WORK ARRAYS -----------------------------------
      LW   = NS*NS*NO + 2*NS*NO + 6*NS + 2
      LIW  = NS
C-----------------------------------------------------------------------
C     DETERMINE IF ENOUGH WORK SPACE HAS BEEN ALLOCATED
C-----------------------------------------------------------------------
      IF (LENWK .LT. LW  .OR.  LENIWK .LT. LIW) THEN
         WRITE (6, '(1X,''ERROR...THE LENGTH OF W MUST BE '',
     &               ''AT LEAST'', 2X, I10, /, 1X, 4X,
     &               ''AND THE LENGTH OF IW MUST BE AT LEAST '',
     &               2X, I10)' )  LW, LIW
         STOP
      ELSE
c$$$         WRITE (6, '(//,1X, ''THE PROBLEM ONLY REQUIRES'', 2X, I10, 3X,
c$$$     &               ''STORAGE LOCATIONS FOR W AND'', 2X, I10, 3X,
c$$$     &               ''STORAGE LOCATIONS FOR IW '',//)')
c$$$     &               LW, LIW
      ENDIF
C-----------------------------------------------------------------------
C    SETS THE POINTERS
C-----------------------------------------------------------------------
      CALL POINTR (NS, NO)
C-----------------------------------------------------------------------
C     CALL TRANSFER
C-----------------------------------------------------------------------
      CALL TRANSF (NS, NO)
C-----------------------------------------------------------------------
C      STOP
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE POINTR (NS, NO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      COMMON /IRTRF/ IWT,    IPATM,  IRU,    IPOLA,  IDIPO,  ISIGM,
     &               IEPSI,  IZROT,  ICOFE,  ICOFL,  ICOFD
      COMMON /IITRF/ ILINA
C-----------------------------------------------------------------------
C    MAXIMUM DIMENSIONS
C-----------------------------------------------------------------------
      LO  = NO
      LS  = NS
C-----------------------------------------------------------------------
C    POINTERS FOR RTST
C-----------------------------------------------------------------------
      IWT    = 1
      IPATM  = IWT   +  LS
      IRU    = IPATM +  1
      IPOLA  = IRU   +  1
      IDIPO  = IPOLA +  LS
      ISIGM  = IDIPO +  LS
      IEPSI  = ISIGM +  LS
      IZROT  = IEPSI +  LS
      ICOFE  = IZROT +  LS
      ICOFL  = ICOFE +  LO*LS
      ICOFD  = ICOFL +  LO*LS
      IRNEXT = ICOFD +  LO*LS*LS
C-----------------------------------------------------------------------
C    POINTERS FOR ITST
C-----------------------------------------------------------------------
      ILINA  = 1
      IINEXT = ILINA +  LS
C-----------------------------------------------------------------------
C      WRITE(6, '(1X,''LW = '',I10, 2X,''LIW = '',I10)' ) IRNEXT,IINEXT
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE TRANSF (NS, NO)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      COMMON /ARAYS/ W(20000), IW(10000)
C-----------------------------------------------------------------------
      COMMON /IRTRF/ IWT,    IPATM,  IRU,    IPOLA,  IDIPO,  ISIGM,
     &               IEPSI,  IZROT,  ICOFE,  ICOFL,  ICOFD
      COMMON /IITRF/ ILINA
C-----------------------------------------------------------------------
      NONS   = NO*NS
      NONSNS = NO*NS*NS
C-----------------------------------------------------------------------
      CALL LTRA (NS, NO, NONS, NONSNS,
     &           W(IWT), W(IPATM), W(IRU),
     &           W(IPOLA), W(IDIPO), W(ISIGM), W(IEPSI), W(IZROT),
     &           W(ICOFE), W(ICOFL), W(ICOFD), IW(ILINA))
C-----------------------------------------------------------------------
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE LTRA  (NS, NO, NONS, NONSNS,
     &                  WT, PATM, RGAZ, POLA, DIPO, SIGM, EPSI,
     &                  ZROT, COFE, COFL, COFD, LINA)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*16 CDUMMY
      LOGICAL LDUMMY
C-----------------------------------------------------------------------
      DIMENSION WT(NS), POLA(NS), DIPO(NS), SIGM(NS), EPSI(NS),
     &          ZROT(NS), LINA(NS),
     &          COFE(NONS), COFL(NONS), COFD(NONSNS)
C-----------------------------------------------------------------------
      LLTP  = 10
      LLMC  = 10
      LLEG  = 11
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C     INITIALIZATION
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C-----------------------------------------------------------------------
C     CHEMKIN-I TRANSPORT LINKING FILE
C     OR OLD VERSIONS OF THE CHEMKIN-II TRANSPORT LINKING FILE
C-----------------------------------------------------------------------
C      OPEN (UNIT=LLTP,STATUS='OLD',FORM='UNFORMATTED',FILE='Linktp')
C         READ (LLTP) IDUMMY, IDUMMY, IDUMMY
C         READ (LLTP) PATM, (WT(K), K=1, NS), (EPSI(K), K=1, NS),
C     &               (SIGM(K), K=1, NS), (DIPO(K), K=1, NS),
C     &               (POLA(K), K=1, NS), (ZROT(K), K=1, NS),
C     &               (LINA(K), K=1, NS), (COFL(N), N=1, NONS),
C     &               (COFE(N), N=1, NONS), (COFD(N), N=1, NONSNS)
C      CLOSE(UNIT=LLTP)
C-----------------------------------------------------------------------
C     CHEMKIN-II TRANSPORT LINKING FILE
C-----------------------------------------------------------------------
      OPEN (UNIT=LLMC,STATUS='OLD',FORM='UNFORMATTED',FILE='tran.bin')
         READ (LLMC) CDUMMY, CDUMMY, LDUMMY
         READ (LLMC) IDUMMY, IDUMMY, IDUMMY, IDUMMY, IDUMMY
         READ (LLMC) PATM, (WT(K), EPSI(K),  SIGM(K),
     &               DIPO(K), POLA(K), ZROT(K), LINA(K), K=1, NS),
     &               (COFL(N), N=1, NONS),
     &               (COFE(N), N=1, NONS), (COFD(N), N=1, NONSNS)
      CLOSE(UNIT=LLMC)
C-----------------------------------------------------------------------
C     MODIFIED TRANSPORT LINKING FILE
C-----------------------------------------------------------------------
      OPEN (UNIT=LLEG,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE='Linkeg')
         WRITE (LLEG) NS, NO, (WT(K), K=1, NS), (EPSI(K), K=1, NS),
     &                (SIGM(K), K=1, NS), (DIPO(K), K=1, NS),
     &                (POLA(K), K=1, NS), (ZROT(K), K=1, NS),
     &                (LINA(K), K=1, NS), (COFE(N), N=1, NONS),
     &                (COFL(N), N=1, NONS), (COFD(N), N=1, NONSNS)
      CLOSE(UNIT=LLEG)
C-----------------------------------------------------------------------
C      STOP
      END
