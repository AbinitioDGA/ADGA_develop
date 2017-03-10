module vq_module
  use parameters_module
  implicit none


contains


subroutine get_vq(v,q,v_r,r_data)
  implicit none
  complex(kind=8),intent(out) :: v(ndim**2,ndim**2)
  real(kind=8),intent(in) :: q(3)
  real(kind=8) :: v_r(:,:,:),r_data(:,:)
  integer :: i

  v=cmplx(0.d0,0.d0,kind=8)
  if (nr.eq.0) then
    v=u
  end if

  do i=1,nr
    v = v + v_r(:,:,i)*exp(2.d0*pi*ci*(r_data(1,i)*q(1)+r_data(2,i)*q(2)+r_data(3,i)*q(3)))
  end do

end subroutine get_vq

end module vq_module
