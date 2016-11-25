module vq_module
  use parameters_module
  implicit none


contains


subroutine get_vq(v,q)
  implicit none
  complex(kind=8),intent(out) :: v(ndim2,ndim2)
  real(kind=8),intent(in) :: q(3)
  integer :: i

  v=0.d0

  do i=1,nr
    v = v + v_r(:,:,i)*exp(ci*(r_data(1,i)*q(1)+r_data(2,i)*q(2)+r_data(3,i)*q(3)))
  end do
end subroutine get_vq

end module vq_module
