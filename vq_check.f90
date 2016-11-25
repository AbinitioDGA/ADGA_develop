program vq_check
  use parameters_module
  use vq_module
  use kq_tools

  implicit none

  integer,parameter           :: nq1=10
  complex(kind=8),allocatable :: vq(:,:,:)
  character(len=150)          :: cmd_arg
  real(kind=8),allocatable    :: rdata(:,:),qdata(:,:),vr(:,:,:)
  integer :: i,j,iq,ir

  
  call getarg(1,cmd_arg)
  filename_vr=trim(cmd_arg) ! declared in parameters_module

  open(unit=20,file=filename_vr)

  read(20,*) nr,ndim,a,b,c
  allocate(vq(ndim**2,ndim**2,nq1**3))
  allocate(vr(ndim**2,ndim**2,nr))
  allocate(rdata(3,nr),qdata(3,nq1**3))

  call read_v_r(vr,rdata)

  close(20)

  call init_k_grid_cubic(qdata,nq1,nq1,nq1,a,b,c)
 

  open(unit=30,file='V_q.dat')
  do iq=1,nq1**3
    do i=1,ndim**2
      do j=1,ndim**2
        vq(i,j,iq) = cmplx(0.,0.,kind=8)
        do ir=1,nr
          vq(i,j,iq)=vq(i,j,iq)+vr(i,j,ir)*exp(-ci*dot_product(qdata(:,iq),rdata(:,ir)))
        end do
      end do
    end do
    write(30,'(F10.4  F10.4  F10.4  )',advance='no')  qdata(1,iq),qdata(2,iq),qdata(3,iq)
    write(30,*) ((real(vq(i,j,iq)),i=1,ndim),j=1,ndim)
  end do 
  close(30)

end program vq_check
