      SUBROUTINE INITCHEMKIN
!
!
      IMPLICIT NONE
!
!
      INCLUDE 'chemkin.common'
!
!
      DATA KERR/.FALSE./, X/KMAX*0.0/, KSYM/KMAX*' '/
!
!
!
!     Initialize CHEMKIN
      OPEN ( LINKCK, FORM = 'UNFORMATTED' , STATUS = 'UNKNOWN' ,
     *               FILE = 'chem.bin' )
!
      CALL CKLEN ( LINKCK , LOUT , LENI , LENR , LENC , IFLAG )
!
      IF ( IFLAG.GT.0 ) STOP '1'
!
      CALL CKINIT ( LENIWK , LENRWK , LENCWK , LINKCK , LOUT , IWORK ,
     *              RWORK , CWORK , IFLAG )
!
      IF ( IFLAG.GT.0 ) STOP '2'
!
      CLOSE ( LINKCK )
!
      CALL CKINDX ( IWORK , RWORK , MM , KK , II , NFIT )
!
!
      NEQ   = KK + 1
      LRW   = 22 + 9*NEQ + 2*NEQ**2
      NVODE = LENR + 1
      NWT   = NVODE + LRW
      NH    = NWT  + KK
      NWDOT = NH   + KK
      NTOT  = NWDOT+ KK - 1
!
      LIW   = 30 + NEQ
      IVODE = LENI + 1
      ITOT  = IVODE + LIW - 1
!
!
      IF ( KK.GT.KMAX .OR. LENRWK.LT.NTOT .OR. LENIWK.LT.ITOT ) THEN
!
         IF ( KK.GT.KMAX ) WRITE (LOUT, *)
     *   'Error...KMAX too small...must be at least' , KK
!
         IF ( LENRWK.LT.NTOT ) WRITE (LOUT, *)
     *   'Error...LENRWK too small...must be at least' , NTOT
!
         IF ( LENIWK.LT.ITOT ) WRITE (LOUT, *)
     *   'Error...LENIWK too small...must be at least' , ITOT
!
         STOP
!
      ENDIF
!
!
      CALL CKSYMS ( CWORK , LOUT , KSYM , IERR )
!
      IF (IERR) KERR = .TRUE.
!
      CALL CKWT   ( IWORK , RWORK , RWORK (NWT) )
!
      CALL CKRP   ( IWORK , RWORK , RU , RUC , PATM )
!
!
!
      END SUBROUTINE INITCHEMKIN
