!------------------------------------------------------------------------------
! MODULE: adim
!------------------------------------------------------------------------------
!> \brief Non-dimensional parameters.
!!
!! This module contains defines all the non-dimensional variables
!! necessary to solve the non-dimensional system of equations or
!! post-treat non-dimensional data.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module adim

  use parameters

  implicit none


!> \brief Non-dimensional derived type declaration.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  type adi_type

     !> Adimensional numbers
     real (dp) :: ma
     real (dp) :: re
     real (dp) :: pr
     real (dp) :: le
     real (dp) :: Sc

     !> Reference state
     real (dp) :: t_ref
     real (dp) :: p_ref
     real (dp) :: rho_ref
     real (dp) :: R_ref
     real (dp) :: r_m_ref
     real (dp) :: W_ref
     real (dp) :: cp_ref
     real (dp) :: mu_ref
     real (dp) :: lbda_ref
     real (dp) :: D_ref
     real (dp) :: time_ref
     real (dp) :: L_ref
     real (dp) :: u_ref

     !> Infinity state
     real (dp) :: t_inf
     real (dp) :: p_inf
     real (dp) :: rho_inf
     real (dp) :: R_inf
     real (dp) :: r_m_inf
     real (dp) :: W_inf
     real (dp) :: cp_inf
     real (dp) :: gamma_inf
     real (dp) :: mu_inf
     real (dp) :: lbda_inf
     real (dp) :: D_inf
     real (dp) :: c_inf
     real (dp) :: u_inf

     !> Constants
     real (dp) :: gamma
     real (dp) :: gm1
     real (dp) :: gm1_i
     real (dp) :: ttrd
     real (dp) :: ggmopr
     real (dp) :: sqgmr

  end type adi_type


contains


!> \brief Print at the screen the adi derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine adi_display (adi)


    type (adi_type) , intent (in)      :: adi !< non-dimensional derived type
    character ( len = 50 )             :: format

    format  = ' ( A19 , 3X , F19.12 ) '

    write (*,*) '########################################'
    write (*,*) 'Non-dimensional numbers'
    write (*,*) '########################################'

    write ( * , format ) 'Mach'     , adi % ma
    write ( * , format ) 'Reynolds' , adi % re
    write ( * , format ) 'Prandtl'  , adi % pr
    write ( * , format ) 'Lewis'    , adi % Le
    write ( * , format ) 'Schmidt'  , adi % Sc


    write (*,*) '########################################'
    write (*,*) 'Infinity state, iptype'
    write (*,*) '########################################'

    write ( * , format ) 'Temperature'         , adi % t_inf
    write ( * , format ) 'Pressure'            , adi % p_inf
    write ( * , format ) 'Density'             , adi % rho_inf
    write ( * , format ) 'Univ. Gas Constant'  , adi % R_inf
    write ( * , format ) 'Mass Gas Constant'   , adi % r_m_inf
    write ( * , format ) 'Molar mass'          , adi % W_inf
    write ( * , format ) 'Spec. heat capacity' , adi % cp_inf
    write ( * , format ) 'Ratio spec. heats'   , adi % gamma_inf
    write ( * , format ) 'Viscosity'           , adi % mu_inf
    write ( * , format ) 'Conductivity'        , adi % lbda_inf
    write ( * , format ) 'Diffusivity'         , adi % D_inf
    write ( * , format ) 'Speed of sound'      , adi % c_inf
    write ( * , format ) 'Velocity'            , adi % u_inf


    write (*,*) '########################################'
    write (*,*) 'Reference state, iptype'
    write (*,*) '########################################'

    write ( * , format ) 'Temperature'         , adi % t_ref
    write ( * , format ) 'Pressure'            , adi % p_ref
    write ( * , format ) 'Density'             , adi % rho_ref
    write ( * , format ) 'Univ. Gas Constant'  , adi % R_ref
    write ( * , format ) 'Mass Gas Constant'   , adi % r_m_ref
    write ( * , format ) 'Molar mass'          , adi % W_ref
    write ( * , format ) 'Spec. heat capacity' , adi % cp_ref
    write ( * , format ) 'Viscosity'           , adi % mu_ref
    write ( * , format ) 'Conductivity'        , adi % lbda_ref
    write ( * , format ) 'Diffusivity'         , adi % D_ref
    write ( * , format ) 'Time'                , adi % time_ref
    write ( * , format ) 'Lenght'              , adi % L_ref
    write ( * , format ) 'Velocity'            , adi % u_ref

  end subroutine adi_display


!> \brief Read from file the adi derived type. (NOT USED, call nowhere... ?!)
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
!  subroutine adi_read (adi)

!    type (adi_type) , intent (inout) :: adi !< non-dimensional derived type

!    open ( unit = unit_adi , file = file_adi )

!        read (unit_adi,*)  adi % ma
!        read (unit_adi,*)  adi % re
!        read (unit_adi,*)  adi % pr
!        read (unit_adi,*)  adi % Le
!        read (unit_adi,*)  adi % Sc

!        read (unit_adi,*)  adi % t_inf
!        read (unit_adi,*)  adi % p_inf
!        read (unit_adi,*)  adi % rho_inf
!        read (unit_adi,*)  adi % R_inf
!        read (unit_adi,*)  adi % r_m_inf
!        read (unit_adi,*)  adi % W_inf
!        read (unit_adi,*)  adi % cp_inf
!        read (unit_adi,*)  adi % gamma_inf
!        read (unit_adi,*)  adi % mu_inf
!        read (unit_adi,*)  adi % lbda_inf
!        read (unit_adi,*)  adi % D_inf
!        read (unit_adi,*)  adi % c_inf
!        read (unit_adi,*)  adi % u_inf

!        read (unit_adi,*)  adi % t_ref
!        read (unit_adi,*)  adi % p_ref
!        read (unit_adi,*)  adi % rho_ref
!        read (unit_adi,*)  adi % R_ref
!        read (unit_adi,*)  adi % r_m_ref
!        read (unit_adi,*)  adi % W_ref
!        read (unit_adi,*)  adi % cp_ref
!        read (unit_adi,*)  adi % mu_ref
!        read (unit_adi,*)  adi % lbda_ref
!        read (unit_adi,*)  adi % D_ref
!        read (unit_adi,*)  adi % time_ref
!        read (unit_adi,*)  adi % L_ref
!        read (unit_adi,*)  adi % u_ref

!    close (unit_adi)

!  end subroutine adi_read


!> \brief Write to file the adi derived type.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine adi_write (adi)

    type (adi_type) , intent (in) :: adi !< non-dimensional derived type

    open ( unit = unit_adi , file = file_adi )

        write (unit_adi,*)  adi % ma
        write (unit_adi,*)  adi % re
        write (unit_adi,*)  adi % pr
        write (unit_adi,*)  adi % Le
        write (unit_adi,*)  adi % Sc

        write (unit_adi,*)  adi % t_inf
        write (unit_adi,*)  adi % p_inf
        write (unit_adi,*)  adi % rho_inf
        write (unit_adi,*)  adi % R_inf
        write (unit_adi,*)  adi % r_m_inf
        write (unit_adi,*)  adi % W_inf
        write (unit_adi,*)  adi % cp_inf
        write (unit_adi,*)  adi % gamma_inf
        write (unit_adi,*)  adi % mu_inf
        write (unit_adi,*)  adi % lbda_inf
        write (unit_adi,*)  adi % D_inf
        write (unit_adi,*)  adi % c_inf
        write (unit_adi,*)  adi % u_inf

        write (unit_adi,*)  adi % t_ref
        write (unit_adi,*)  adi % p_ref
        write (unit_adi,*)  adi % rho_ref
        write (unit_adi,*)  adi % R_ref
        write (unit_adi,*)  adi % r_m_ref
        write (unit_adi,*)  adi % W_ref
        write (unit_adi,*)  adi % cp_ref
        write (unit_adi,*)  adi % mu_ref
        write (unit_adi,*)  adi % lbda_ref
        write (unit_adi,*)  adi % D_ref
        write (unit_adi,*)  adi % time_ref
        write (unit_adi,*)  adi % L_ref
        write (unit_adi,*)  adi % u_ref

    close (unit_adi)

  end subroutine adi_write


end module adim
