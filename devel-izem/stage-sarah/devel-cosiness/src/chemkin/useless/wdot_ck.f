      SUBROUTINE WDOT_CK ( PRE , TEM , Y_ , WDOT_ )
!
!
      IMPLICIT NONE
!
!
      INCLUDE 'chemkin.common'
!
!
      DOUBLE PRECISION PRE    , TEM     , Y_     , WDOT_
!
!
      DIMENSION Y_(*)   , WDOT_(*) 
!
!
!
      ! Calculates the MOLAR production rates of each species
      CALL CKWYP ( PRE , TEM , Y_ , IWORK , RWORK , WDOT_ )
!
!
!
c$$$      WRITE (*,*) 'WDOT_CK'
c$$$      WRITE (*,*) WDOT_(1)
!
!
      END SUBROUTINE WDOT_CK
