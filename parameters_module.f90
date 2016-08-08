module parameters_module
  implicit none
  public

  complex(kind=8) ci
  parameter (ci=(0.d0,1.d0))
  double precision,parameter :: pi=3.1415926535897932385d0
  integer :: nkp, ndim, ndims, ndim2, maxdim,i1,i2,i3,i4
  integer :: nkpx,nkpy,nkpz,nkp1,nqpx,nqpy,nqpz,nqp1
  integer :: iwmax, iwbmax, iwfmax, iwbmax_small, iwfmax_small,nk_frac
  integer :: iwstart,iwstop
  double precision :: mu, beta
  double precision, allocatable :: k_data(:,:)
  integer :: nqp,nkp_eom
  integer,allocatable :: q_data(:),k_data_eom(:)
  character(len=150) :: filename,filename_vertex,filename_umatrix,output_dir,filename_q_path
  logical :: orb_sym,small_freq_box,full_chi0
  logical :: do_eom,do_chi
  logical :: q_path_susc,k_path_eom,q_vol
  integer :: vertex_type
  integer,parameter :: full_g4=0,connected_g4=1,chi_g4=2

contains
subroutine read_config()
  implicit none
  character(len=100) :: cmd_arg
  character(len=100) :: config_file
  character(len=150) :: str_tmp
  integer :: int_tmp_1,int_tmp_2,int_tmp_3

  call getarg(1,cmd_arg)
  config_file=trim(cmd_arg)
  write(*,*) 'Reading config: ',config_file

  open(unit=1,file=config_file)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_vertex=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_umatrix=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2,int_tmp_3
  orb_sym=int_tmp_1
  small_freq_box=int_tmp_2
  full_chi0=int_tmp_3

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2
  iwfmax_small=int_tmp_1
  iwbmax_small=int_tmp_2

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2,int_tmp_3
  q_vol=int_tmp_1
  k_path_eom=int_tmp_2
  q_path_susc=int_tmp_3

  read(1,*)
  read(1,*)
  read(1,*) str_tmp
  filename_q_path=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) nk_frac

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  output_dir=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2
  do_chi=int_tmp_1
  do_eom=int_tmp_2
  close(1)
end subroutine read_config


subroutine init() 
  implicit none
  maxdim = ndim*ndim*2*iwfmax_small
  ndim2 = ndim*ndim
  if (full_chi0) then
    iwstart=-iwmax+iwbmax
    iwstop=iwmax-iwbmax-1
  else
    iwstart=-iwfmax_small
    iwstop=iwfmax_small-1
  end if
!  nkp=nkpx*nkpy*nkpz

  !define k-grid:
  nkp1 = nkp**(1./3.)  
  nkpx=nkp1
  nkpy=nkp1
  nkpz=nkp1

  if(mod(nkp1,nk_frac).ne.0)then
    stop 'mismatch between k- and q-grid!'
  endif

  nqp1 = nkp1/nk_frac
  nqpx=nqp1
  nqpy=nqp1
  nqpz=nqp1 

  if (q_vol) then
    nqp = nqpx*nqpy*nqpz
  if (mod(nkpx,nqpx).ne.0 .or. mod(nkpy,nqpy).ne.0 .or. mod(nkpz,nqpz).ne.0) then
     stop 'mismatch between k- and q-grid!'
   endif
  end if
  if (q_path_susc .or. k_path_eom) then
    nkp1=nkpx
    nqp1=nqpx
  end if
end subroutine init

end module parameters_module
