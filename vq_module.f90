module vq_module
  use parameters_module
  use hdf5_module
  use aux
  implicit none


contains

!===============================================================================================================================
subroutine read_vq(iq, v)
  use aux
  implicit none
  integer, intent(in) :: iq
  complex(kind=8), intent(out) :: v(ndim2,ndim2)
  complex(kind=8) :: vq(ndim,ndim,ndim,ndim)
  integer(hid_t) :: vq_file_id, grp_id, iq_id, iq_space_id
  integer :: err, ind, i, j, k, l, i1, i2
  integer :: nmembers, imembers, itype
  character(len=20) :: name_buffer
  integer(hsize_t), dimension(2) :: iq_dims, iq_maxdims
  integer(hsize_t), dimension(1) :: vq_dims
  double precision :: vq_tmp_r(nqp), vq_tmp_i(nqp)

  call h5fopen_f(filename_vq, h5f_acc_rdonly_f, vq_file_id, err)

  call h5dopen_f(vq_file_id, ".axes/Q-points", iq_id, err)
  call h5dget_space_f(iq_id, iq_space_id, err)
  call h5sget_simple_extent_dims_f(iq_space_id, iq_dims, iq_maxdims, err)
  if(iq_dims(2) .ne. nqp) then
     write(*,*) 'Inconsistent number of q-points in V^q!', iq_dims(2),'/',nqp
     stop
  endif

  vq = 0.d0

  call h5gn_members_f(vq_file_id, "/", nmembers, err)
  do imembers = 1,nmembers - 1
     call h5gget_obj_info_idx_f(vq_file_id, "/", imembers, name_buffer, itype, err)

     read(name_buffer,'(I5.5)') ind
     call index2component_band(ndim,ind,i,j,k,l)
     call h5dopen_f(vq_file_id, name_buffer, grp_id, err)
     call h5dread_f(grp_id, type_r_id, vq_tmp_r, vq_dims, err)
     call h5dread_f(grp_id, type_i_id, vq_tmp_i, vq_dims, err)

     vq(i,j,k,l) = vq_tmp_r(iq)+ci*vq_tmp_i(iq)

     call h5dclose_f(grp_id, err)
  enddo
  call h5fclose_f(vq_file_id, err)

  v = 0.d0
  i2 = 0
  do l=1,ndim
     do j=1,ndim
        i2 = i2+1
        i1 = 0
        do i=1,ndim
           do k=1,ndim
              i1 = i1+1
              v(i1,i2) = vq(i,j,k,l)
           enddo
        enddo
     enddo
  enddo

end subroutine read_vq
!========================================================================================================


!========================================================================================================
subroutine read_u(u, u_tilde)
  implicit none
  real(kind=8) :: u_tmp(ndim,ndim,ndim,ndim), u_tilde_tmp(ndim,ndim,ndim,ndim)
  real(kind=8) :: u_value
  complex(kind=8), intent(out) :: u(ndim**2, ndim**2), u_tilde(ndim**2, ndim**2)
  integer :: n,i,j,k,l,i1,i2


  open(21,file=filename_umatrix,status='old')
  read(21,*)
  do n=1,ndim**4
     read(21,*) i, j, k, l, u_value
     u_tmp(i,j,k,l) = u_value
     u_tilde_tmp(i,j,l,k) = u_value
  enddo
  close(21)

  !go into compound index:
  u = 0.d0
  u_tilde = 0.d0
  i2 = 0
  do l=1,ndim
     do j=1,ndim
        i2 = i2+1
        i1 = 0
        do i=1,ndim
           do k=1,ndim
              i1 = i1+1
              u(i1,i2) = u_tmp(i,j,k,l)
              u_tilde(i1,i2) = u_tilde_tmp(i,j,k,l)
           enddo
        enddo
     enddo
  enddo

end subroutine read_u


subroutine create_u(u, u_tilde)
  implicit none
  real(kind=8) :: u_tmp(ndim,ndim,ndim,ndim), u_tilde_tmp(ndim,ndim,ndim,ndim)
  real(kind=8) :: u_value
  complex(kind=8), intent(out) :: u(ndim**2, ndim**2), u_tilde(ndim**2, ndim**2)
  integer :: n,i,j,k,l,i1,i2,ineq

  allocate(Umat(ndim,ndim,ndim,ndim))
  Umat=0.d0

  do i=1,ndim
  do j=1,ndim
  do k=1,ndim
  do l=1,ndim

    ineq=index2ineq(nineq,ndims,i,j,k,l) 
    if (ineq .eq. 0) then ! not on the same impurity
      goto 100
    endif
  
    ! DD - VALUES
    if (index2cor(nineq,ndims,i,j,k,l)) then
      if (i .eq. j .and. k .eq. l)  then
        if (i .eq. k) then
          Umat(i,j,k,l) = Udd(ineq) ! onsite density-density -- 1 1 1 1
          goto 100
        else if (interaction_mode(ineq) .eq. 1) then
          Umat(i,j,k,l) = Jdd(ineq) ! kanamori double hopping -- 1 1 2 2
          goto 100
        endif
      endif

      if (i .eq. k .and. j .eq. l) then
        Umat(i,j,k,l) = Vdd(ineq) ! screened density density -- 1 2 1 2
        goto 100
      endif
      if (i .eq. l .and. j .eq. k .and. interaction_mode(ineq) .eq. 1) then
        Umat(i,j,k,l) = Jdd(ineq) ! kanamori spin flip -- 1 2 2 1
        goto 100
      endif

    ! PP - VALUES
    else if(index2uncor(nineq,ndims,i,j,k,l)) then
      if (i .eq. j .and. k .eq. l)  then
        if (i .eq. k) then
          Umat(i,j,k,l) = Upp(ineq) ! onsite density-density -- 1 1 1 1
          goto 100
        else if (interaction_mode(ineq) .eq. 1) then
          Umat(i,j,k,l) = Jpp(ineq) ! kanamori double hopping -- 1 1 2 2
          goto 100
        endif
      endif

      if (i .eq. k .and. j .eq. l) then
        Umat(i,j,k,l) = Vpp(ineq) ! screened density density -- 1 2 1 2
        goto 100
      endif
      if (i .eq. l .and. j .eq. k .and. interaction_mode(ineq) .eq. 1) then
        Umat(i,j,k,l) = Jpp(ineq) ! kanamori spin flip -- 1 2 2 1
        goto 100
      endif

    ! DP - VALUES
    else
      if (i .eq. j .and. k .eq. l)  then
        if (i .eq. k) then
          Umat(i,j,k,l) = Udp(ineq) ! onsite density-density -- 1 1 1 1
          goto 100
        else if (interaction_mode(ineq) .eq. 1) then
          Umat(i,j,k,l) = Jdp(ineq) ! kanamori double hopping -- 1 1 2 2
          goto 100
        endif
      endif

      if (i .eq. k .and. j .eq. l) then
        Umat(i,j,k,l) = Vdp(ineq) ! screened density density -- 1 2 1 2
        goto 100
      endif
      if (i .eq. l .and. j .eq. k .and. interaction_mode(ineq) .eq. 1) then
        Umat(i,j,k,l) = Jdp(ineq) ! kanamori spin flip -- 1 2 2 1
        goto 100
      endif
    endif

  100 continue
  ! 100 WRITE(10,'(4I10,F15.8)') I,J,K,L,UMAT(I,J,K,L)

  enddo
  enddo
  enddo
  enddo

  deallocate(interaction_mode)
  deallocate(Udd,Vdd,Jdd,Upp,Vpp,Jpp,Udp,Vdp,Jdp)

  !go into compound index:
  u = 0.d0
  u_tilde = 0.d0
  i2 = 0
  do l=1,ndim
     do j=1,ndim
        i2 = i2+1
        i1 = 0
        do i=1,ndim
           do k=1,ndim
              i1 = i1+1
              u(i1,i2) = Umat(i,j,k,l)
              u_tilde(i1,i2) = Umat(i,j,k,l)
           enddo
        enddo
     enddo
  enddo

end subroutine create_u
!======================================================================================================


!subroutine read_v_r(v_r,r_data)
!  implicit none
!  real(kind=8) v_r(:,:,:)!(ndim**2,ndim**2,nr)
!  real(kind=8) r_data(3,nr),v_r_real(ndim2,ndim2)
!  integer :: nr_file,ir,i,j,nd

!  open(unit=2,file=filename_vr)
!  read(2,*) nr_file,nd,a,b,c
!  if (nr_file .ne. nr) then
!    write(*,*) 'V(r) file says there are',nr_file,'r points. '
!    write(*,*) 'Please adapt config file.'
!    stop
!  end if
!  if (nd .ne. ndim) then
!    write(*,*) 'V(r) file says there are',nd,'orbitals. '
!    write(*,*) 'Please adapt config file.'
!    stop
!  end if

!  do ir=1,nr
!    read(2,*) (r_data(i,ir),i=1,3)
!! TODO: correctly read multi-band components and go to compound index.
!    do i=1,nd**2
!       read(2,*) (v_r(i,j,ir),j=1,nd**2)
!    enddo
!  enddo

!  close(2)

!end subroutine read_v_r

!subroutine get_vq(v,q,v_r,r_data)
!  implicit none
!  complex(kind=8),intent(out) :: v(ndim**2,ndim**2)
!  real(kind=8),intent(in) :: q(3)
!  real(kind=8) :: v_r(:,:,:),r_data(:,:)
!  integer :: i

!  v=cmplx(0.d0,0.d0,kind=8)
!  if (nr.eq.0) then
!    v=u
!  end if

!  do i=1,nr
!    v = v + v_r(:,:,i)*exp(2.d0*pi*ci*(r_data(1,i)*q(1)+r_data(2,i)*q(2)+r_data(3,i)*q(3)))
!  end do

!end subroutine get_vq

end module vq_module
