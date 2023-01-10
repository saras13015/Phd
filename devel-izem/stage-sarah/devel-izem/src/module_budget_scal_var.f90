!------------------------------------------------------------------------------
! MODULE: budget_scal_var
!------------------------------------------------------------------------------
!> \brief Budget equations for the scalar variance transport equations.
!!
!! ### Each term has its own sign!
!!
!! In module_stat_plot.f90 specific signs and other customizations are
!! applied to each term to make the _addition_ of all the terms in the
!! equation equal to zero.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

module budget_scal_var


  use parameters

  use parallel

  use input

  use type_thd

  use deriv

  use thermodynamics

  use tools_post


  implicit none


contains


!> \brief Calculate the mean convection term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine scal_var_mean_conv ( dx_i , dy_i , dz_i , ux , vy , wz , rYa , mean_conv )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux        !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vy        !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: wz        !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rYa       !< density times species array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: mean_conv !< mean convection term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3 , &
                                                           wrk4 , wrk5 , wrk6


    allocate ( wrk1 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk2 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk3 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk4 (sx:ex,sy:ey,sz:ez)                   , &
               wrk5 (sx:ex,sy:ey,sz:ez)                   , &
               wrk6 (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate scal_var_mean_conv'


    if ( ndim == 2 ) then ! 2D problem


       wrk1 (sx:ex,sy:ey,sz:ez) = rYa (sx:ex,sy:ey,sz:ez) * ux (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = rYa (sx:ex,sy:ey,sz:ez) * vy (sx:ex,sy:ey,sz:ez)

       call comm_one (wrk1)           ; call comm_one (wrk2)
       call dx ( dx_i , wrk1 , wrk4 ) ; call dy ( dy_i , wrk2 , wrk5 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                mean_conv (i,j,k) = wrk4 (i,j,k) + wrk5 (i,j,k)
                mean_conv (i,j,k) = 0.5_dp * mean_conv (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       wrk1 (sx:ex,sy:ey,sz:ez) = rYa (sx:ex,sy:ey,sz:ez) * ux (sx:ex,sy:ey,sz:ez)
       wrk2 (sx:ex,sy:ey,sz:ez) = rYa (sx:ex,sy:ey,sz:ez) * vy (sx:ex,sy:ey,sz:ez)
       wrk3 (sx:ex,sy:ey,sz:ez) = rYa (sx:ex,sy:ey,sz:ez) * wz (sx:ex,sy:ey,sz:ez)

       call comm_one (wrk1)           ; call comm_one (wrk2)           ; call comm_one (wrk3)
       call dx ( dx_i , wrk1 , wrk4 ) ; call dy ( dy_i , wrk2 , wrk5 ) ; call dz ( dz_i , wrk3 , wrk6 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                mean_conv (i,j,k) = wrk4 (i,j,k) + wrk5 (i,j,k) + wrk6 (i,j,k)
                mean_conv (i,j,k) = 0.5_dp * mean_conv (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( wrk1 , wrk2 , wrk3 , &
                 wrk4 , wrk5 , wrk6 )


  end subroutine scal_var_mean_conv


!> \brief Calculate the turbulent transport term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine scal_var_turb_transp ( rho , ux , vy , wz , Ya , turb_transp )


    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho         !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux          !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vy          !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: wz          !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: Ya          !< species array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: turb_transp !< turbulent transport


    integer (ip)                                        :: i , j , k
    real (dp)                                           :: wrk


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk                   = rho (i,j,k) * Ya (i,j,k) * Ya (i,j,k)
                turb_transp (i,j,k,1) = wrk * ux (i,j,k)
                turb_transp (i,j,k,2) = wrk * vy (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk                   = rho (i,j,k) * Ya (i,j,k) * Ya (i,j,k)
                turb_transp (i,j,k,1) = wrk * ux (i,j,k)
                turb_transp (i,j,k,2) = wrk * vy (i,j,k)
                turb_transp (i,j,k,3) = wrk * wz (i,j,k)
             end do
          end do
       end do


    end if


  end subroutine scal_var_turb_transp


!> \brief Calculate the turbulent transport term (last call).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine scal_var_turb_transp_last ( dx_i , dy_i , dz_i , turb_transp , res )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i        !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i        !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i        !< inverted dx array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: turb_transp !< previous turbulent transport term
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: res         !< final turbulent transport term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3 , &
                                                           wrk4 , wrk5 , wrk6


    allocate ( wrk1 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk2 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk3 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk4 (sx:ex,sy:ey,sz:ez)                   , &
               wrk5 (sx:ex,sy:ey,sz:ez)                   , &
               wrk6 (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate scal_var_turb_transp_last'


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk1 (i,j,k) = turb_transp (i,j,k,1)
                wrk2 (i,j,k) = turb_transp (i,j,k,2)
             end do
          end do
       end do

       call comm_one (wrk1)           ; call comm_one (wrk2)
       call dx ( dx_i , wrk1 , wrk4 ) ; call dy ( dy_i , wrk2 , wrk5 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                res (i,j,k) = wrk4 (i,j,k) + wrk5 (i,j,k)
                res (i,j,k) = 0.5_dp * res (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk1 (i,j,k) = turb_transp (i,j,k,1)
                wrk2 (i,j,k) = turb_transp (i,j,k,2)
                wrk3 (i,j,k) = turb_transp (i,j,k,3)
             end do
          end do
       end do

       call comm_one (wrk1)           ; call comm_one (wrk2)           ; call comm_one (wrk3)
       call dx ( dx_i , wrk1 , wrk4 ) ; call dy ( dy_i , wrk2 , wrk5 ) ; call dz ( dz_i , wrk3 , wrk6 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                res (i,j,k) = wrk4 (i,j,k) + wrk5 (i,j,k) + wrk6 (i,j,k)
                res (i,j,k) = 0.5_dp * res (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( wrk1 , wrk2 , wrk3 , &
                 wrk4 , wrk5 , wrk6 )


  end subroutine scal_var_turb_transp_last


!> \brief Calculate the turbulent production term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine scal_var_turb_prod ( rho , ux , vy , wz , Ya , turb_prod )


    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: rho       !< density
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: ux        !< x-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vy        !< y-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: wz        !< z-component of the velocity
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: Ya        !< species array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: turb_prod !< turbulent production term


    integer (ip)                                        :: i , j , k
    real (dp)                                           :: wrk


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk                 = rho (i,j,k) * Ya (i,j,k)
                turb_prod (i,j,k,1) = wrk * ux (i,j,k)
                turb_prod (i,j,k,2) = wrk * vy (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk                 = rho (i,j,k) * Ya (i,j,k)
                turb_prod (i,j,k,1) = wrk * ux (i,j,k)
                turb_prod (i,j,k,2) = wrk * vy (i,j,k)
                turb_prod (i,j,k,3) = wrk * wz (i,j,k)
             end do
          end do
       end do


    end if


  end subroutine scal_var_turb_prod


!> \brief Calculate the turbulent production term (last call).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine scal_var_turb_prod_last ( dx_i , dy_i , dz_i , turb_prod , Ya , res )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i      !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i      !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i      !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: turb_prod !< previous turbulent production term
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: Ya        !< species array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: res       !< final turbulent production term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3


    allocate ( wrk1 (sx:ex,sy:ey,sz:ez) , &
               wrk2 (sx:ex,sy:ey,sz:ez) , &
               wrk3 (sx:ex,sy:ey,sz:ez) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate scal_var_turb_prod_last'


    if ( ndim == 2 ) then ! 2D problem


       call comm_one (Ya)
       call dx ( dx_i , Ya , wrk1 ) ; call dy ( dy_i , Ya , wrk2 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                res (i,j,k) = turb_prod (i,j,k,1) * wrk1 (i,j,k) + &
                              turb_prod (i,j,k,2) * wrk2 (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       call comm_one (Ya)
       call dx ( dx_i , Ya , wrk1 ) ; call dy ( dy_i , Ya , wrk2 ) ; call dz ( dz_i , Ya , wrk3 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                res (i,j,k) = turb_prod (i,j,k,1) * wrk1 (i,j,k) + &
                              turb_prod (i,j,k,2) * wrk2 (i,j,k) + &
                              turb_prod (i,j,k,3) * wrk3 (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( wrk1 , wrk2 , wrk3 )


  end subroutine scal_var_turb_prod_last


!> \brief Calculate the molecular diffusion term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine scal_var_mol_diff ( index_scalar , thd , dx_i , dy_i , dz_i , T , W_i , dm , tdr , rd , v , &
                                 D_equiv , fluc_Ya , Ya , mol_diff )


    integer (ip) , intent (in)                                           :: index_scalar !< species index
    type (thd_type) , intent (in)                                        :: thd          !< thermodynamic derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i         !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i         !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i         !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: T            !< temperature
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: W_i          !< inverted molar mass
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: dm           !< diffusion matrix
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: tdr          !< thermal diffusivity
    real (dp) , dimension (:,:,:,:,:) , allocatable , intent (in)        :: rd           !< density times diffusion
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: v            !< conserved variables array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: D_equiv      !< "equivalent" diffusion matrix
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: fluc_Ya      !< species fluctuation
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: Ya           !< species array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: mol_diff     !< molecular diffusion term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: wrk
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3
    real (dp) , dimension (:,:,:,:) , allocatable       :: rYVx , rYVy , rYVz


    allocate ( wrk1 (sx:ex,sy:ey,sz:ez)     , &
               wrk2 (sx:ex,sy:ey,sz:ez)     , &
               wrk3 (sx:ex,sy:ey,sz:ez)     , &
               rYVx (sx:ex,sy:ey,sz:ez,nrv) , &
               rYVy (sx:ex,sy:ey,sz:ez,nrv) , &
               rYVz (sx:ex,sy:ey,sz:ez,nrv) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate scal_var_mol_diff'


    if ( ndim == 2 ) then ! 2D problem


       if ( index_scalar == 0 .or. index_scalar == nrv+npv ) then

          call comm_one (Ya)
          call dx ( dx_i , Ya , wrk1 ) ; call dy ( dy_i , Ya , wrk2 )

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   wrk                = v (i,j,k,1) * D_equiv (i,j,k) * fluc_Ya (i,j,k)
                   mol_diff (i,j,k,1) = wrk * wrk1 (i,j,k) ! these terms are positive
                   mol_diff (i,j,k,2) = wrk * wrk2 (i,j,k)
                end do
             end do
          end do

       else

          if ( vis .and. .not. eglib ) then
             call rhoYaVa_mix ( thd , dx_i , dy_i , dx_i , v , W_i , dm , rYVx , rYVy , rYVx )
          else if ( vis .and. eglib ) then
             call rhoYaVa_eglib ( thd , dx_i , dy_i , dx_i , T , W_i , tdr , rd , v , rYVx , rYVy , rYVx )
          end if

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   mol_diff (i,j,k,1) = - rYVx (i,j,k,index_scalar) * fluc_Ya (i,j,k) ! these terms are negative
                   mol_diff (i,j,k,2) = - rYVy (i,j,k,index_scalar) * fluc_Ya (i,j,k)
                end do
             end do
          end do

       end if


    else if ( ndim == 3 ) then ! 3D problem


       if ( index_scalar == 0 .or. index_scalar == nrv+npv ) then

          call comm_one (Ya)
          call dx ( dx_i , Ya , wrk1 ) ; call dy ( dy_i , Ya , wrk2 ) ; call dz ( dz_i , Ya , wrk3 )

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   wrk                = v (i,j,k,1) * D_equiv (i,j,k) * fluc_Ya (i,j,k)
                   mol_diff (i,j,k,1) = wrk * wrk1 (i,j,k) ! these terms are positive
                   mol_diff (i,j,k,2) = wrk * wrk2 (i,j,k)
                   mol_diff (i,j,k,3) = wrk * wrk3 (i,j,k)
                end do
             end do
          end do

       else

          if ( vis .and. .not. eglib ) then
             call rhoYaVa_mix ( thd , dx_i , dy_i , dz_i , v , W_i , dm , rYVx , rYVy , rYVz )
          else if ( vis .and. eglib ) then
             call rhoYaVa_eglib ( thd , dx_i , dy_i , dz_i , T , W_i , tdr , rd , v , rYVx , rYVy , rYVz )
          end if

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   mol_diff (i,j,k,1) = - rYVx (i,j,k,index_scalar) * fluc_Ya (i,j,k) ! these terms are negative
                   mol_diff (i,j,k,2) = - rYVy (i,j,k,index_scalar) * fluc_Ya (i,j,k)
                   mol_diff (i,j,k,3) = - rYVz (i,j,k,index_scalar) * fluc_Ya (i,j,k)
                end do
             end do
          end do

       end if


    end if


    deallocate ( wrk1 , wrk2 , wrk3 , &
                 rYVx , rYVy , rYVz )


  end subroutine scal_var_mol_diff


!> \brief Calculate the molecular diffusion term (last call).
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine scal_var_mol_diff_last ( dx_i , dy_i , dz_i , mol_diff , res )


    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i     !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i     !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i     !< inverted dz array
    real (dp) , dimension (:,:,:,:) , allocatable , intent (inout)       :: mol_diff !< previous molecular diffusion term
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: res      !< final molecular diffusion term


    integer (ip)                                        :: ok , i , j , k
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3 , &
                                                           wrk4 , wrk5 , wrk6


    allocate ( wrk1 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk2 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk3 (sx-ng:ex+ng,sy-ng:ey+ng,sz-ng:ez+ng) , &
               wrk4 (sx:ex,sy:ey,sz:ez)                   , &
               wrk5 (sx:ex,sy:ey,sz:ez)                   , &
               wrk6 (sx:ex,sy:ey,sz:ez)                   , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate scal_var_mol_diff_last'


    if ( ndim == 2 ) then ! 2D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk1 (i,j,k) = mol_diff (i,j,k,1)
                wrk2 (i,j,k) = mol_diff (i,j,k,2)
             end do
          end do
       end do

       call comm_one (wrk1)           ; call comm_one (wrk2)
       call dx ( dx_i , wrk1 , wrk4 ) ; call dy ( dy_i , wrk2 , wrk5 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                res (i,j,k) = wrk4 (i,j,k) + wrk5 (i,j,k)
             end do
          end do
       end do


    else if ( ndim == 3 ) then ! 3D problem


       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                wrk1 (i,j,k) = mol_diff (i,j,k,1)
                wrk2 (i,j,k) = mol_diff (i,j,k,2)
                wrk3 (i,j,k) = mol_diff (i,j,k,3)
             end do
          end do
       end do

       call comm_one (wrk1)           ; call comm_one (wrk2)           ; call comm_one (wrk3)
       call dx ( dx_i , wrk1 , wrk4 ) ; call dy ( dy_i , wrk2 , wrk5 ) ; call dz ( dz_i , wrk3 , wrk6 )

       do k = sz,ez
          do j = sy,ey
             do i = sx,ex
                res (i,j,k) = wrk4 (i,j,k) + wrk5 (i,j,k) + wrk6 (i,j,k)
             end do
          end do
       end do


    end if


    deallocate ( wrk1 , wrk2 , wrk3 , &
                 wrk4 , wrk5 , wrk6 )


  end subroutine scal_var_mol_diff_last


!> \brief Calculate the dissipation term.
!!
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)\n
!! Benguerir, Morocco

  subroutine scal_var_dissip ( index_scalar , thd , dx_i , dy_i , dz_i , T , W_i , dm , tdr , rd , v , &
                               D_equiv , fluc_Ya , Ya , dissip )


    integer (ip) , intent (in)                                           :: index_scalar !< species index
    type (thd_type) , intent (in)                                        :: thd          !< thermodynamic derived type
    real (dp) , dimension (:) , allocatable , intent (in)                :: dx_i         !< inverted dx array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dy_i         !< inverted dy array
    real (dp) , dimension (:) , allocatable , intent (in)                :: dz_i         !< inverted dz array
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: T            !< temperature
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: W_i          !< inverted molar mass
    real (dp) , dimension (:,:,:,:) , allocatable , intent (in)          :: dm           !< diffusion matrix
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)          :: tdr          !< thermal diffusivity
    real (dp) , allocatable , dimension (:,:,:,:,:) , intent (in)        :: rd           !< density times diffusion
    real (dp) , allocatable , dimension (:,:,:,:) , intent (in)          :: v            !< conserved variables array
    real (dp) , allocatable , dimension (:,:,:) , intent (in)            :: D_equiv      !< "equivalent" diffusion matrix
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)         :: fluc_Ya      !< species fluctuation
    real (dp) , allocatable , dimension (:,:,:) , intent (inout)         :: Ya           !< species array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: dissip       !< dissipation term


    integer (ip)                                        :: ok , i , j , k
    real (dp)                                           :: wrk
    real (dp) , dimension (:,:,:) , allocatable         :: wrk1 , wrk2 , wrk3 , &
                                                           wrk4 , wrk5 , wrk6
    real (dp) , dimension (:,:,:,:) , allocatable       :: rYVx , rYVy , rYVz

    allocate ( wrk1 (sx:ex,sy:ey,sz:ez)     , &
               wrk2 (sx:ex,sy:ey,sz:ez)     , &
               wrk3 (sx:ex,sy:ey,sz:ez)     , &
               wrk4 (sx:ex,sy:ey,sz:ez)     , &
               wrk5 (sx:ex,sy:ey,sz:ez)     , &
               wrk6 (sx:ex,sy:ey,sz:ez)     , &
               rYVx (sx:ex,sy:ey,sz:ez,nrv) , &
               rYVy (sx:ex,sy:ey,sz:ez,nrv) , &
               rYVz (sx:ex,sy:ey,sz:ez,nrv) , &
               stat = ok )
    if ( ok > 0 ) stop 'error allocate scal_var_dissip'


    if ( ndim == 2 ) then ! 2D problem


       if ( index_scalar == 0 .or. index_scalar == nrv+npv ) then

          call comm_one (Ya) ; call comm_one (fluc_Ya)
          call dx ( dx_i , Ya , wrk1 )      ; call dy ( dy_i , Ya , wrk2 )
          call dx ( dx_i , fluc_Ya , wrk4 ) ; call dy ( dy_i , fluc_Ya , wrk5 )

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   wrk            = v (i,j,k,1) * D_equiv (i,j,k)
                   dissip (i,j,k) = wrk * ( wrk1 (i,j,k) * wrk4 (i,j,k) + &
                                            wrk2 (i,j,k) * wrk5 (i,j,k) )
                   dissip (i,j,k) = - dissip (i,j,k) ! this term is negative
                end do
             end do
          end do

       else

          call comm_one (fluc_Ya)
          call dx ( dx_i , fluc_Ya , wrk4 ) ; call dy ( dy_i , fluc_Ya , wrk5 )

          if ( vis .and. .not. eglib ) then
             call rhoYaVa_mix ( thd , dx_i , dy_i , dx_i , v , W_i , dm , rYVx , rYVy , rYVx )
          else if ( vis .and. eglib ) then
             call rhoYaVa_eglib ( thd , dx_i , dy_i , dx_i , T , W_i , tdr , rd , v , rYVx , rYVy , rYVx )
          end if

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   dissip (i,j,k) = rYVx (i,j,k,index_scalar) * wrk4 (i,j,k) + & ! this term is positive
                                    rYVy (i,j,k,index_scalar) * wrk5 (i,j,k)
                end do
             end do
          end do

       end if


    else if ( ndim == 3 ) then ! 3D problem


       if ( index_scalar == 0 .or. index_scalar == nrv+npv ) then

          call comm_one (Ya) ; call comm_one (fluc_Ya)
          call dx ( dx_i , Ya , wrk1 )      ; call dy ( dy_i , Ya , wrk2 )      ; call dz ( dz_i , Ya , wrk3 )
          call dx ( dx_i , fluc_Ya , wrk4 ) ; call dy ( dy_i , fluc_Ya , wrk5 ) ; call dz ( dz_i , fluc_Ya , wrk6 )

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   wrk            = v (i,j,k,1) * D_equiv (i,j,k)
                   dissip (i,j,k) = wrk * ( wrk1 (i,j,k) * wrk4 (i,j,k) + &
                                            wrk2 (i,j,k) * wrk5 (i,j,k) + &
                                            wrk3 (i,j,k) * wrk6 (i,j,k) )
                   dissip (i,j,k) = - dissip (i,j,k) ! this term is negative
                end do
             end do
          end do

       else

          call comm_one (fluc_Ya)
          call dx ( dx_i , fluc_Ya , wrk4 ) ; call dy ( dy_i , fluc_Ya , wrk5 ) ; call dz ( dz_i , fluc_Ya , wrk6 )

          if ( vis .and. .not. eglib ) then
             call rhoYaVa_mix ( thd , dx_i , dy_i , dz_i , v , W_i , dm , rYVx , rYVy , rYVz )
          else if ( vis .and. eglib ) then
             call rhoYaVa_eglib ( thd , dx_i , dy_i , dz_i , T , W_i , tdr , rd , v , rYVx , rYVy , rYVz )
          end if

          do k = sz,ez
             do j = sy,ey
                do i = sx,ex
                   dissip (i,j,k) = rYVx (i,j,k,index_scalar) * wrk4 (i,j,k) + & ! this term is positive
                                    rYVy (i,j,k,index_scalar) * wrk5 (i,j,k) + &
                                    rYVz (i,j,k,index_scalar) * wrk6 (i,j,k)
                end do
             end do
          end do

       end if


    end if


    deallocate ( wrk1 , wrk2 , wrk3 , &
                 wrk4 , wrk5 , wrk6 , &
                 rYVx , rYVy , rYVz )


  end subroutine scal_var_dissip


end module budget_scal_var
