module parameters_module
  implicit none
  public

  complex(kind=8) ci
  parameter (ci=(0.d0,1.d0))
  integer :: nkp, ndim, ndims, ndim2, maxdim
  integer :: iwmax, iwbmax, iwfmax, iwbmax_small, iwfmax_small
  double precision :: mu
  integer :: nqp
  

end module parameters_module
