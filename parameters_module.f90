module parameters_module
  implicit none
  public

  complex(kind=8) ci
  parameter (ci=(0.d0,1.d0))
  double precision,parameter :: pi=3.1415926535897932385d0
  integer :: nkp, ndim, ndim2, maxdim, nineq
  integer, allocatable :: ndims(:,:)
  integer :: nkpx,nkpy,nkpz,nkp1,nqpx,nqpy,nqpz,nqp1
  integer :: iwmax, iwbmax, iwfmax, iwbmax_small, iwfmax_small,nk_frac
  integer :: iwstart,iwstop
  integer :: ineq
  double precision :: mu, beta
  double precision, allocatable :: k_data(:,:), r_data(:,:)
  integer :: nqp,nkp_eom, idp
  integer,allocatable :: q_data(:),k_data_eom(:)
  character(len=150) :: filename, filename_umatrix, filename_vq, filename_hk, output_dir, filename_q_path
  logical :: orb_sym,full_chi0
  logical :: do_eom,do_chi,do_vq
  logical :: q_path_susc,k_path_eom,q_vol,read_ext_hk
  integer :: vertex_type
  integer,parameter :: full_g4=0,connected_g4=1,chi_g4=2
  integer :: nr ! number of r-points in extrapolated V(r)
  character(len=150) filename_vr ! filename of extrapolated V(r)
  real(kind=8) :: a,b,c ! lattice spacing
  character(len=150) ::  filename_vertex, filename_vertex_sym
  character(len=100) :: config_file


contains

subroutine read_config()
  implicit none
  character(len=100) :: cmd_arg
  character(len=150) :: str_tmp
  integer :: int_tmp_1,int_tmp_2,int_tmp_3,int_tmp_4

  call getarg(1,cmd_arg)
  config_file=trim(cmd_arg)

  open(unit=1,file=config_file)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2
  do_chi=int_tmp_1
  do_eom=int_tmp_2

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1
  nineq=int_tmp_1
  allocate(ndims(nineq,2))

  read(1,*)
  read(1,*)
  read(1,*) ((ndims(ineq, idp), idp=1,2), ineq=1,nineq)
  ndim=sum(ndims)

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_hk=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) nkpx, nkpy, nkpz, nqpx, nqpy, nqpz

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2,int_tmp_3,int_tmp_4
  q_vol=int_tmp_1
  k_path_eom=int_tmp_2
  q_path_susc=int_tmp_3
  read_ext_hk=int_tmp_4

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_q_path=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1
  do_vq = int_tmp_1

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_vq = trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_umatrix=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1
  orb_sym=int_tmp_1

  read(1,*)
  read(1,*)
  read(1,*) vertex_type

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_vertex_sym=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2
  iwfmax_small=int_tmp_1
  iwbmax_small=int_tmp_2

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  output_dir=trim(str_tmp)

  ! read(1,*)
  ! read(1,*)
  ! read(1,*) nr ! VR

  ! if (nr .gt. 0) then
  !   read(1,*)
  !   read(1,*)
  !   read(1,'(A)') str_tmp
  !   filename_vr = trim(str_tmp)
  ! end if

  close(1)

end subroutine read_config


subroutine init()
  implicit none
  integer :: i,j,k,l,n
  maxdim = ndim*ndim*2*iwfmax_small
  ndim2 = ndim*ndim
  if (full_chi0) then
    iwstart=-iwmax+iwbmax
    iwstop=iwmax-iwbmax-1
  else
    iwstart=-iwfmax_small
    iwstop=iwfmax_small-1
  end if

  ! Since currently the BZ sum has to go over all points of the Hamiltonian, we
  ! do a consistency check here.
  if (nkpx*nkpy*nkpz .ne. nkp) then
    stop 'Wrong number of k points in config file.'
  end if
  nkp=nkpx*nkpy*nkpz


  if (q_vol) then
    nqp = nqpx*nqpy*nqpz
    if (mod(nkpx,nqpx).ne.0 .or. mod(nkpy,nqpy).ne.0 .or. mod(nkpz,nqpz).ne.0) then
      stop 'mismatch between k- and q-grid!'
    endif
  end if

  if (q_path_susc .or. k_path_eom) then
    stop 'q paths currently not stable'
    nkp1=nkpx
    nqp1=nqpx
  end if


end subroutine init

end module parameters_module
