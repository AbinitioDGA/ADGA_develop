program vq_check
  use parameters_module
  use vq_module
  use kq_tools

  implicit none

  integer,parameter           :: nqx=1000,nqy=1,nqz=1
  complex(kind=8),allocatable :: vq(:,:),vq_sum(:,:)
  character(len=150)          :: cmd_arg
  real(kind=8),allocatable    :: rdata(:,:),qdata(:,:),vr(:,:,:)
  integer                     :: i,j,iq,ir

  
  call getarg(1,cmd_arg)
  filename_vr=trim(cmd_arg) ! declared in parameters_module

  open(unit=20,file=filename_vr)

  read(20,*) nr,ndim,a,b,c
  allocate(vq(ndim**2,ndim**2),vq_sum(ndim**2,ndim**2))
  allocate(vr(ndim**2,ndim**2,nr))
  allocate(rdata(3,nr),qdata(3,nqx*nqy*nqz))

  call read_v_r(vr,rdata)

  close(20)

  call init_k_grid_cubic(qdata,nqx,nqy,nqz,a,b,c)
 
  vq_sum=0.d0
  open(unit=30,file='V_q_exp.dat')
  open(unit=31,file='V_q_exp_ind.dat')
  do iq=1,nqx*nqy*nqz !nq1**3
    call get_vq(vq,qdata(:,iq),vr,rdata)
    vq_sum=vq_sum+vq
    write(30,'(F10.4  F10.4  F10.4  )',advance='no')  qdata(1,iq),qdata(2,iq),qdata(3,iq)
    write(30,*) ((real(vq(i,j)),i=1,ndim**2),j=1,ndim**2)
    write(31,'(I6)',advance='no') iq
    write(31,*) ((real(vq(i,j)),i=1,ndim**2),j=1,ndim**2)
  end do 
  close(30)
  close(31)

  write(*,*) vq_sum/(nqx*nqy*nqz)
end program vq_check
