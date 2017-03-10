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
  double precision, allocatable :: k_data(:,:), r_data(:,:)
  integer :: nqp,nkp_eom
  integer,allocatable :: q_data(:),k_data_eom(:)
  character(len=150) :: filename,filename_vertex,filename_umatrix,filename_hk,output_dir,filename_q_path
  logical :: orb_sym,small_freq_box,full_chi0
  logical :: do_eom,do_chi
  logical :: q_path_susc,k_path_eom,q_vol,read_ext_hk
  integer :: vertex_type
  integer,parameter :: full_g4=0,connected_g4=1,chi_g4=2
  integer :: nr ! number of r-points in extrapolated V(r)
  character(len=150) filename_vr ! filename of extrapolated V(r)
  real(kind=8),allocatable :: v_r(:,:,:), u_tmp(:,:,:,:), u_tilde_tmp(:,:,:,:)
  real(kind=8) :: a,b,c ! lattice spacing
  real(kind=8) :: u_value
  complex(kind=8),allocatable :: u(:,:), u_tilde(:,:)

contains

subroutine read_config()
  implicit none
  character(len=100) :: cmd_arg
  character(len=100) :: config_file
  character(len=150) :: str_tmp
  integer :: int_tmp_1,int_tmp_2,int_tmp_3,int_tmp_4

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
  read(1,*) nkpx, nkpy, nkpz, nqpx, nqpy, nqpz, ndim

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  output_dir=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2
  do_chi=int_tmp_1
  do_eom=int_tmp_2

  read(1,*)
  read(1,*)
  read(1,*) vertex_type

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_hk=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) nr

  if (nr .gt. 0) then
    read(1,*)
    read(1,*)
    read(1,'(A)') str_tmp
    filename_vr = trim(str_tmp)
  end if

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
  !  stop 'q paths currently not stable'
    nkp1=nkpx
    nqp1=nqpx
  end if

  if (nr .eq. 0) then
    write(*,*) 'Run without V(q)'
  else
    allocate(r_data(3,nr))
    allocate(v_r(ndim2,ndim2,nr))

    write(*,*) 'Read V(r)'

    call read_v_r(v_r,r_data)
  end if

  ! read Umatrix from file:
  allocate(u_tmp(ndim, ndim, ndim, ndim))
  allocate(u_tilde_tmp(ndim, ndim, ndim, ndim))

  open(21,file=filename_umatrix,status='old')
  read(21,*)

  do n=1,ndim**4

     read(21,*) i, j, k, l, u_value   
     u_tmp(i,j,k,l) = u_value
     u_tilde_tmp(i,j,l,k) = u_value

  enddo

  close(21)

  
  !go into compound index:
  allocate(u(ndim**2, ndim**2))
  allocate(u_tilde(ndim**2, ndim**2))
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

  deallocate(u_tmp, u_tilde_tmp)
  

end subroutine init

subroutine read_v_r(v_r,r_data)
  implicit none
  real(kind=8) v_r(:,:,:)!(ndim**2,ndim**2,nr)
  real(kind=8) r_data(3,nr),v_r_real(ndim2,ndim2)
  integer :: nr_file,ir,i,j,nd

  open(unit=2,file=filename_vr)
  read(2,*) nr_file,nd,a,b,c
  if (nr_file .ne. nr) then
    write(*,*) 'V(r) file says there are',nr_file,'r points. '
    write(*,*) 'Please adapt config file.'
    stop
  end if
  if (nd .ne. ndim) then
    write(*,*) 'V(r) file says there are',nd,'orbitals. '
    write(*,*) 'Please adapt config file.'
    stop
  end if

  do ir=1,nr
    read(2,*) (r_data(i,ir),i=1,3)
! TODO: correctly read multi-band components and go to compound index.
    do i=1,nd**2
       read(2,*) (v_r(i,j,ir),j=1,nd**2)
    enddo
  enddo

  close(2)

end subroutine read_v_r

end module parameters_module
