module kq_tools
  use parameters_module
  implicit none

  interface k_vector
    module procedure k_vector_1, k_vector_3
  end interface k_vector
  interface k_index
    module procedure k_index_1, k_index_3
  end interface k_index

contains

subroutine generate_q_vol(nqpx,nqpy,nqpz,qdata)
  implicit none
  integer :: nqpx,nqpy,nqpz
  integer :: qdata(nqpx*nqpy*nqpz),i,j,k,i1

  i1=0
  do i=0,nqpx-1
    do j=0,nqpy-1
      do k=0,nqpz-1
        i1 = i1+1
        qdata(i1)=k_index(i*nkpx/nqpx,j*nkpy/nqpy,k*nkpz/nqpz)
      enddo
    enddo
  enddo

end subroutine generate_q_vol

subroutine index_kq(ind)
  implicit none
  integer, intent(out) :: ind(nkp,nqp)

  integer ikp,jkp
  ind = 0

  do ikp=1,nkp
    do jkp=1,nqp
      ind(ikp,jkp)=k_minus_q(ikp,q_data(jkp))
    end do
  end do
end subroutine index_kq

subroutine index_kq_eom(ind)
  implicit none
  integer, intent(out) :: ind(nkp_eom,nqp)

  integer ikp,jkp
  ind = 0

  do ikp=1,nkp_eom
    do jkp=1,nqp
      ind(ikp,jkp)=k_minus_q(k_data_eom(ikp),q_data(jkp))
    end do
  end do
end subroutine index_kq_eom

! The following function calculates the index of \vec{k} - \vec{q}.
! It uses only integers
! \vec{k} is associated to (ix,iy,iz)
! \vec{q} is associated to (lx,ly,lz)
! k-space is assumed to have nkpx*nkpy*nkpz points
! q-space is assumed to have nqpx*nqpy*nqpz points,
! where each element of the q-space has to be an element of the k-space.
! subtractions are done in integers,
! fold-back to BZ is achieved by modulo division.
function k_minus_q(ik,iq)
  implicit none
  integer :: ik,iq,k_minus_q
  integer :: ix,iy,iz,lx,ly,lz

  call k_vector(ik,ix,iy,iz)
  call k_vector(iq,lx,ly,lz)

  k_minus_q=1+mod(nkpz+iz-lz,nkpz) + &
              mod(nkpy+iy-ly,nkpy)*nkpz + &
              mod(nkpx+ix-lx,nkpx)*nkpy*nkpz

end function k_minus_q

function k_plus_q(ik,iq)
  implicit none
  integer :: ik,iq,k_plus_q
  integer :: ix,iy,iz,lx,ly,lz

  call k_vector(ik,ix,iy,iz)
  call k_vector(iq,lx,ly,lz)

  k_plus_q=1+mod(nkpz+iz+lz,nkpz) + &
              mod(nkpy+iy+ly,nkpy)*nkpz + &
              mod(nkpx+ix+lx,nkpx)*nkpy*nkpz

end function k_plus_q

function k_index_1(k)
  implicit none
  integer,intent(in) :: k(3)
  integer :: k_index_1

  k_index_1 = 1 + k(3) + k(2)*nkpz + k(1)*nkpy*nkpz

end function k_index_1

function k_index_3(kx,ky,kz)
  implicit none
  integer :: k_index_3,kx,ky,kz

  k_index_3 = 1 + kz + ky*nkpz + kx*nkpy*nkpz

end function k_index_3

subroutine k_vector_1(ik,k)
  implicit none
  integer,intent(in) :: ik
  integer,intent(out) :: k(3)

  k(3)=mod(ik-1,nkpz)
  k(2)=mod((ik-1)/nkpz,nkpy)
  k(1)=(ik-1)/(nkpy*nkpz)
end subroutine k_vector_1

subroutine k_vector_3(ik,kx,ky,kz)
  implicit none
  integer :: ik,kx,ky,kz

  kz=mod(ik-1,nkpz)
  ky=mod((ik-1)/nkpz,nkpy)
  kx=(ik-1)/(nkpy*nkpz)
end subroutine k_vector_3

subroutine qdata_from_file()
  use parameters_module

  implicit none
  integer :: iostatus,iq
  character(100) :: str_tmp
  real(kind=8) :: qx,qy,qz

  iostatus=0
  open(unit=101,file=filename_qdata)

  nqp=-1
  do while (iostatus.eq.0)
    read(101,*,iostat=iostatus) str_tmp
    nqp=nqp+1
  end do
  close(101)
  !write(*,*) nqp,' q points'

  allocate(q_data(nqp))
  open(unit=101,file=filename_qdata)
  do iq=1,nqp
    ! We read three real numbers.
    ! If all of them are zero, it is the gamma point.
    ! If one of them is larger or equal to 1, the coordinates are cast to integers
    ! and assumed to be given in integer basis [0,nkpi-1]
    ! If neither of above is true, the coordinates are assumed to lie in the interval [0,1).
    read(101,*) qx,qy,qz
    if (qx .eq. 0 .and. qy .eq. 0 .and. qz .eq. 0) then
      q_data(iq) = k_index(0,0,0) ! gamma point
    else if (qx .ge. 1 .or. qy .ge. 1 .or. qz .ge. 1) then
      q_data(iq) = k_index(int(qx),int(qy),int(qz)) ! cast to integers
    else
      q_data(iq) = k_index(nint(qx*nkpx),nint(qy*nkpy),nint(qz*nkpz)) ! round to nearest integers
    end if
  end do
  close(101)
  !write(*,*) 'q data',q_data

end subroutine qdata_from_file

subroutine qdata_susc_from_file()
  use parameters_module

  implicit none
  integer :: iostatus,iq
  character(100) :: str_tmp
  real(kind=8) :: qx,qy,qz

  logical :: warning
  warning = .false.

  iostatus=0
  open(unit=101,file=filename_qdata_susc)

  nqpphbar=-1
  do while (iostatus.eq.0)
    read(101,*,iostat=iostatus) str_tmp
    nqpphbar=nqpphbar+1
  end do
  close(101)
  !write(*,*) nqpphbar,' q points'

  allocate(q_data_phbar(nqpphbar))
  allocate(q_half_data_phbar(nqpphbar))
  open(unit=101,file=filename_qdata_susc)
  do iq=1,nqpphbar
    ! We read three real numbers.
    ! If all of them are zero, it is the gamma point.
    ! If one of them is larger or equal to 1, the coordinates are cast to integers
    ! and assumed to be given in integer basis [0,nkpi-1]
    ! If neither of above is true, the coordinates are assumed to lie in the interval [0,1).
    read(101,*) qx,qy,qz
    if (qx .eq. 0 .and. qy .eq. 0 .and. qz .eq. 0) then
      q_data_phbar(iq) = k_index(0,0,0) ! gamma point
      q_half_data_phbar(iq) = q_data_phbar(iq) ! the same point
    else if (qx .ge. 1 .or. qy .ge. 1 .or. qz .ge. 1) then
      q_data_phbar(iq) = k_index(int(qx),int(qy),int(qz)) ! cast to integers
      if (mod(int(qx),2) /= 0 .or. mod(int(qy),2) /= 0 .or. mod(int(qz),2) /= 0) then
        warning = .true.
      endif
      q_half_data_phbar(iq) = k_index(floor(qx/2),floor(qy/2),floor(qz/2)) ! cast to integers
    else
      q_data_phbar(iq) = k_index(nint(qx*nkpx),nint(qy*nkpy),nint(qz*nkpz)) ! round to nearest integers
    end if
  end do
  close(101)

  if (warning) then
    write(*,*) 'WARNING: only use q-vectors that are divisible by 2 for jj-correlators'
  endif
  !write(*,*) 'q data',q_data

end subroutine qdata_susc_from_file

subroutine kdata_from_file()
  ! defining k_data_eom by mapping to the contigous k_data array analogous to q_data
  use parameters_module

  implicit none
  integer :: iostatus,ik
  character(100) :: str_tmp
  real(8) :: kx,ky,kz

  iostatus=0
  open(unit=101,file=filename_kdata)

  nkp_eom = -1
  do while (iostatus.eq.0)
    read(101,*,iostat=iostatus) str_tmp
    nkp_eom=nkp_eom+1
  end do
  close(101)

  allocate(k_data_eom(nkp_eom))
  open(unit=101,file=filename_kdata)
  do ik=1,nkp_eom
    read(101,*) kx,ky,kz
    if (kx .eq. 0 .and. ky .eq. 0 .and. kz .eq. 0) then
      k_data_eom(ik) = k_index(0,0,0) ! gamma point
    else if (kx .ge. 1 .or. ky .ge. 1 .or. kz .ge. 1) then
      k_data_eom(ik) = k_index(int(kx),int(ky),int(kz)) ! cast to integers
    else
      k_data_eom(ik) = k_index(nint(kx*nkpx),nint(ky*nkpy),nint(kz*nkpz)) ! round to nearest integers
    end if
  end do
  close(101)

end subroutine kdata_from_file

end module kq_tools
