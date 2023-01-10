!------------------------------------------------------------------------------
! MODULE: autocorrelations
!------------------------------------------------------------------------------
!> \brief Autocorrelations functions.
!!
!! This module provides autocorrelations functions used when
!! post-treating the data.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
module autocorrelations


  use parameters
  use parallel
  use input


  implicit none


contains


!> \brief Calculate the correlation in z-direction (intermediate operation).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine correlation_Z ( inp , grid , vv , corr )


    type (inp_type) , intent (in)                                        :: inp  !< input derived type
    type (inp_grid) , intent (in)                                        :: grid !< grid derived type
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: vv   !< variable to correlate
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: corr !< intermediate correlation result


    integer (ip)                                :: length
    integer (ip)                                :: np , npl , npg , ix , fx , iy , fy , iz , fz
    integer (ip)                                :: ok , i , j , k , px , py , pz , n_z2
    integer (ip) , dimension (nproc)            :: length_array , displ_array
    integer (ip) , dimension (ndimmax)          :: mycoords
    real (dp) , dimension (ndimmax)             :: pt
    real (dp) , dimension (:) , allocatable     :: wrk
    real (dp) , dimension (:) , allocatable     :: global_array , local_array
    real (dp) , dimension (:,:,:) , allocatable :: global_vv


    if ( rank == rank_default ) then
       allocate  ( wrk (ntz)               , &
                   global_vv (ntx,nty,ntz) , &
                   stat = ok )
       if ( ok > 0 ) stop 'error allocate correlationZ'
    end if


    ! first: sequentialize global array


    if ( nproc > 1 ) then ! parallel treatment


       npl = (ex-sx+1)*(ey-sy+1)*(ez-sz+1)
       npg = ntx*nty*ntz


       allocate ( local_array (npl)  , &
                  global_array (npg) , &
                  stat = ok )
       if ( ok > 0 ) stop 'error allocate correlationZ_2'


       np = 0
       call domain_select ( 0 , ix , fx , iy , fy , iz , fz )
       do k = iz , fz
          do j = iy , fy
             do i = ix , fx
                np               = np + 1
                local_array (np) = vv (i,j,k)
             end do
          end do
       end do

       call MPI_ALLGATHER ( npl          , 1 , MPI_INTEGER , &
                            length_array , 1 , MPI_INTEGER , &
                            MPI_COMM_WORLD , mpicode )

       length = 0
       do np = 1 , nproc
          displ_array (np) = length
          length           = length + length_array (np)
       end do

       call MPI_GATHERV ( local_array , npl , MPI_DOUBLE_PRECISION  ,           &
                          global_array , length_array , displ_array ,           &
                          MPI_DOUBLE_PRECISION , 0 , MPI_COMM_WORLD , mpicode )


       if ( rank == rank_default ) then ! this part is only made by the first process

          length     = 0
          np         = 0
          mycoords (:) = 0
          do pz = 1 , dims (1)

             iz = ( mycoords(1)*ntz ) / dims(1) + 1
             fz = ( ( mycoords(1)+1 )*ntz ) / dims(1)

             mycoords (1) = mycoords (1) + 1
             mycoords (2) = 0
             mycoords (3) = 0

             do py = 1 , dims (2)

                iy = ( mycoords(2)*nty ) / dims(2) + 1
                fy = ( ( mycoords(2)+1 )*nty ) / dims(2)

                mycoords (2) = mycoords (2) + 1
                mycoords (3) = 0

                do px = 1 , dims (3)

                   ix = ( mycoords(3)*ntx ) / dims(3) + 1
                   fx = ( ( mycoords(3)+1 )*ntx ) / dims(3)

                   mycoords (3) = mycoords (3) + 1

                   do k = iz,fz
                      do j = iy,fy
                         do i = ix,fx
                            length = length+1
                            global_vv (i,j,k) = global_array (length)
                         end do
                      end do
                   end do

                   np = np+1

                end do
             end do
          end do

       end if


    else ! sequential treatement


       global_vv (1:ntx,1:nty,1:ntz) = vv (1:ntx,1:nty,1:ntz)


    end if


    ! second: calculate the auto-correlation functions


    if ( rank == rank_default ) then


       px = ntx ; py = nty ; pz = ntz ! initialize this to avoid bugs (verified)


       pt (1) = inp % corrspec_coord (1)
       pt (2) = inp % corrspec_coord (2)
       pt (3) = inp % corrspec_coord (3)


       ! calculate points
       do i = ntx , 1 , -1
          if ( grid % xt (i) > pt (1) ) px = i
       end do
       do j = nty , 1 , -1
          if ( grid % yt (j) > pt (2) ) py = j
       end do
       if ( ndim == 3 ) then
          do k = ntz , 1 , -1
             if ( grid % zt (k) > pt (3) ) pz = k
          end do
       else
          pz = sz
       end if


       ! calculate their correlation
       wrk (:) = 0.0_dp
       n_z2 = ntz/2-1
       do k = 0 , n_z2
          wrk (pz+k) = global_vv (px,py,pz) * global_vv (px,py,pz+k)
       end do


    end if


    ! third: parallelize global array


    if ( nproc > 1 ) then ! parallel treatment


       if ( rank == rank_default ) then ! this part is only made by the first process

          do k = 1,ntz
             do j = 1,nty
                do i = 1,ntx
                   global_vv (i,j,k) = wrk (k)
                end do
             end do
          end do

          length       = 0
          np           = 0
          mycoords (:) = 0
          do pz = 1 , dims (1)

             iz = ( mycoords(1)*ntz ) / dims(1) + 1
             fz = ( ( mycoords(1)+1 )*ntz ) / dims(1)

             mycoords (1) = mycoords (1) + 1
             mycoords (2) = 0
             mycoords (3) = 0

             do py = 1 , dims (2)

                iy = ( mycoords(2)*nty ) / dims(2) + 1
                fy = ( ( mycoords(2)+1 )*nty ) / dims(2)

                mycoords (2) = mycoords (2) + 1
                mycoords (3) = 0

                do px = 1 , dims (3)

                   ix = ( mycoords(3)*ntx ) / dims(3) + 1
                   fx = ( ( mycoords(3)+1 )*ntx ) / dims(3)

                   mycoords (3) = mycoords (3) + 1

                   do k = iz,fz
                      do j = iy,fy
                         do i = ix,fx
                            length = length+1
                            global_array (length) = global_vv (i,j,k)
                         end do
                      end do
                   end do

                   np = np+1

                end do
             end do
          end do

       end if

       call MPI_SCATTERV ( global_array , length_array , displ_array ,           &
                           MPI_DOUBLE_PRECISION , local_array , npl ,            &
                           MPI_DOUBLE_PRECISION , 0 , MPI_COMM_WORLD , mpicode )

       np = 0
       call domain_select ( 0 , ix , fx , iy , fy , iz , fz )
       do k = iz , fz
          do j = iy , fy
             do i = ix , fx
                np           = np + 1
                corr (i,j,k) = local_array (np)
             end do
          end do
       end do

       deallocate ( global_array , local_array )


    else ! sequential treatment


       do k = 1,ntz
          do j = 1,nty
             do i = 1,ntx
                corr (i,j,k) = wrk (k)
             end do
          end do
       end do


    end if


    if ( rank == rank_default ) deallocate ( wrk , global_vv )


  end subroutine correlation_Z


!> \brief Calculate the correlation in z-direction (final operation).
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine last_correlation_Z ( inp , grid , corr )


    type (inp_type) , intent (in)                                        :: inp  !< input derived type
    type (inp_grid) , intent (in)                                        :: grid !< grid derived type
!    real (dp) , dimension (:) , allocatable , intent (in)                :: xt   !< absolute x-coordinate array
!    real (dp) , dimension (:) , allocatable , intent (in)                :: yt   !< absolute y-coordinate array
!    real (dp) , dimension (:) , allocatable , intent (in)                :: zt   !< absolute z-coordinate array
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: corr !< final correlation result


    integer (ip)                        :: myrank , np , n_z2 , i , j , k , px , py , pz
    integer (ip)                        :: ix , fx , iy , fy , iz , fz
    real (dp)                           :: tmp
    integer (ip) , dimension (ndimmax)  :: mycoords
    real (dp) , dimension (ndimmax)     :: pt


    myrank = 0


    pt (1) = inp % corrspec_coord (1)
    pt (2) = inp % corrspec_coord (2)
    pt (3) = inp % corrspec_coord (3)


    np           = 0
    mycoords (:) = 0
    do pz = 1 , dims (1)

       iz = ( mycoords(1)*ntz ) / dims(1) + 1
       fz = ( ( mycoords(1)+1 )*ntz ) / dims(1)

       mycoords (1) = mycoords (1) + 1
       mycoords (2) = 0
       mycoords (3) = 0

       do py = 1 , dims (2)

          iy = ( mycoords(2)*nty ) / dims(2) + 1
          fy = ( ( mycoords(2)+1 )*nty ) / dims(2)

          mycoords (2) = mycoords (2) + 1
          mycoords (3) = 0

          do px = 1 , dims (3)

             ix = ( mycoords(3)*ntx ) / dims(3) + 1
             fx = ( ( mycoords(3)+1 )*ntx ) / dims(3)

             mycoords (3) = mycoords (3) + 1

             if ( pt (1) >= grid % xt (ix) .and. pt (1) <= grid % xt (fx) .and. &
                  pt (2) >= grid % yt (iy) .and. pt (2) <= grid % yt (fy) .and. &
                  pt (3) >= grid % zt (iz) .and. pt (3) <= grid % zt (fz) ) then
                myrank = np
             end if

             np = np+1

          end do
       end do
    end do


    ! calculate points
    px = ntx ; py = nty ; pz = ntz
    do i = ntx , 1 , -1
       if ( grid % xt (i) > pt (1) ) px = i
    end do
    do j = nty , 1 , -1
       if ( grid % yt (j) > pt (2) ) py = j
    end do
    if ( ndim == 3 ) then
       do k = ntz , 1 , -1
          if ( grid % zt (k) > pt (3) ) pz = k
       end do
    else
       pz = sz
    end if

    tmp = 0.0_dp
    if ( rank == myrank ) then
       tmp = 1.0_dp / corr (px,py,pz)
    end if

    call MPI_BCAST ( tmp , 1 , MPI_DOUBLE_PRECISION ,    &
                     myrank , MPI_COMM_WORLD , mpicode )

    n_z2 = ntz/2-1
    do k = sz,min(ez,pz+n_z2)
       do j = sy,ey
          do i = sx,ex
             if ( k < pz ) then
                corr (i,j,k) = 1.0_dp
             else
                corr (i,j,k) = corr (i,j,k) * tmp
             end if
          end do
       end do
    end do


  end subroutine last_correlation_Z


!> \brief Calculate the integral length in z-direction.
!!
!> \author
!! Modeling, Simulation and Data Analysis (MSDA)\n
!! Université Mohammed VI Polytechnique (UM6P)
  subroutine integral_length_Z ( inp , grid , corr , length )


    type (inp_type) , intent (in)                                        :: inp    !< input derived type
    type (inp_grid) , intent (in)                                        :: grid   !< grid derived type
    real (dp) , dimension (:,:,:) , allocatable , intent (in)            :: corr   !< correlation (z-direction)
    real (dp) , dimension (:,:,:) , allocatable , intent (inout)         :: length !< integral length (z-direction)


    integer (ip)                        :: i , j , k , px , py , pz , n_z2
    real (dp)                           :: tmp1 , tmp2
    real (dp) , dimension (ndimmax)     :: pt


    length = 0.0_dp


    px = ntx ; py = nty ; pz = ntz ! initialize this to avoid bugs (verified)


    pt (1) = inp % corrspec_coord (1)
    pt (2) = inp % corrspec_coord (2)
    pt (3) = inp % corrspec_coord (3)


    ! calculate points
    do i = ntx , 1 , -1
       if ( grid % xt (i) > pt (1) ) px = i
    end do
    do j = nty , 1 , -1
       if ( grid % yt (j) > pt (2) ) py = j
    end do
    if ( ndim == 3 ) then
       do k = ntz , 1 , -1
          if ( grid % zt (k) > pt (3) ) pz = k
       end do
    else
       pz = sz
    end if


    tmp1 = 0.0_dp
    n_z2 = ntz/2-1
    do k = sz,min(ez,pz+n_z2)
       if ( k >= pz ) tmp1 = tmp1 + corr (sx,sy,k) / grid % dz_i (k) ! sx and sy index are valid
    end do


    if ( comm1dz /= MPI_COMM_NULL ) then
       call mpi_allreduce ( tmp1 , tmp2 , 1 , MPI_DOUBLE_PRECISION , &
                            MPI_SUM , comm1dz , mpicode )
    else
       tmp2 = tmp1
    end if


    length (sx:ex,sy:ey,sz:ez) = tmp2


  end subroutine integral_length_Z


end module autocorrelations
