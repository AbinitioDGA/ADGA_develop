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
  double precision, allocatable :: iw_data(:), iwb_data(:), iwf_data(:)
  integer :: iwstart,iwstop
  double precision :: mu, beta
  complex(kind=8), allocatable :: u(:,:), u_tilde(:,:)
  double precision, allocatable :: k_data(:,:), r_data(:,:)
  complex(kind=8), allocatable :: hk(:,:,:),dc(:,:)
  complex(kind=8), allocatable :: siw(:,:),giw(:,:) 
  complex(kind=8), allocatable :: n_dga(:), n_dmft(:), n_fock(:,:,:)
  integer :: nqp,nkp_eom, idp
  integer,allocatable :: q_data(:),k_data_eom(:)
  logical :: orb_sym,full_chi0
  logical :: do_eom,do_chi,do_vq
  logical :: q_path_susc,k_path_eom,q_vol,read_ext_hk
  integer :: vertex_type
  logical :: exist_p
  integer,parameter :: full_g4=0,connected_g4=1,chi_g4=2
  integer :: nr ! number of r-points in extrapolated V(r)
  real(kind=8) :: a,b,c ! lattice spacing
  character(len=150) filename_vr ! filename of extrapolated V(r)
  character(len=150) :: filename_vertex, filename_vertex_sym
  character(len=150) :: filename, filename_umatrix, filename_vq
  character(len=150) :: filename_hk, output_dir, filename_q_path
  character(len=100) :: config_file


contains

subroutine read_config()
  implicit none
  character(len=100) :: cmd_arg
  character(len=150) :: str_tmp
  integer :: int_tmp_1,int_tmp_2,int_tmp_3,int_tmp_4
  integer :: ineq

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
  integer :: i
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

  ! create arrays with Matsubara frequencies
  allocate(iw_data(-iwmax:iwmax-1),iwb_data(-iwbmax:iwbmax),iwf_data(-iwfmax:iwfmax-1))
  do i=-iwmax,iwmax-1
    iw_data(i)=pi*(2*i+1)/beta
  end do
  do i=-iwfmax,iwfmax-1
    iwf_data(i)=pi*(2*i+1)/beta
  end do
  do i=-iwbmax,iwbmax
    iwb_data(i)=pi*2*i/beta
  end do
  
  allocate(u(ndim**2,ndim**2), u_tilde(ndim**2,ndim**2))

end subroutine init

subroutine finalize()
  implicit none
  deallocate(iw_data,iwf_data,iwb_data)
end subroutine finalize


subroutine check_freq_range(mpi_wrank,master)
  implicit none
  integer :: mpi_wrank, master

  if (iwfmax_small .le. 0) then
    iwfmax_small = iwfmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Calculating with maximum number of fermionic frequencies &
                  =', iwfmax
    endif
  endif

  if (iwbmax_small .lt. 0) then
    iwbmax_small = iwbmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Calculating with maximum number of bosonic frequencies &
                  =', iwbmax
    endif
  endif

  if (iwfmax_small .gt. iwfmax) then
    iwfmax_small = iwfmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Error: Wrong input for fermionic frequencies'
      write(*,*) 'Calculating with maximum number of fermionic frequencies &
                  =', iwfmax
    endif
  endif

  if (iwbmax_small .gt. iwbmax) then
    iwbmax_small = iwbmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Error: Wrong input for bosonic frequencies'
      write(*,*) 'Calculating with maximum number of bosonic frequencies &
                  =', iwbmax
    endif
  endif

end subroutine check_freq_range


subroutine check_config()
  implicit none
  logical :: there
  integer :: ineq

  if (do_eom .lt. 0 .or. do_eom .gt. 1 .or. do_chi .lt. 0 .or. do_chi .gt. 1) then
    stop "Error: Choose appropriate calculation mode (config line 5)"
  endif

  exist_p = .false.
  do ineq=1,nineq
    if(ndims(ineq,2) .ne. 0) exist_p=.true.
  enddo

  inquire (file=trim(filename_hk), exist=there)
  if (.not. there) then
    stop "Error: Hamiltonian file does not exist"
  endif

  if (q_vol .lt. 0 .or. q_vol .gt. 1 .or. &
    k_path_eom .lt. 0 .or. k_path_eom .gt. 1 .or. &
    q_path_susc .lt. 0 .or. q_path_susc .gt. 1 .or. &
    read_ext_hk .lt. 0 .or. read_ext_hk .gt. 1 ) then
    stop "Error: Choose appropriate settings (config line 20)"
  endif

  if (q_vol .eq. 0) then
    inquire (file=trim(filename_q_path), exist=there)
    if (.not. there) then
      stop "Error: Q-Path file does not exist"
    endif
  endif

  if (do_vq .lt. 0 .or. do_vq .gt. 1) then
    stop "Error: Choose appropriate V(q) mode (config line 29)"
  endif

  if (do_vq .eq. 1) then
    inquire (file=trim(filename_vq), exist=there)
    if (.not. there) then
      stop "Error: V(Q) file does not exist"
    endif
  endif

  inquire (file=trim(filename_vertex), exist=there)
  if (.not. there) then
    stop "Error: One-Particle data file does not exist"
  endif

  inquire (file=trim(filename_umatrix), exist=there)
  if (.not. there) then
    stop "Error: Umatrix file does not exist"
  endif

  if (orb_sym .lt. 0 .or. orb_sym .gt. 1) then
    stop "Error: Choose appropriate symmetrization mode (config line 38)"
  endif

  if (vertex_type .lt. 0 .or. vertex_type .gt. 2) then
    stop "Error: Choose appropriate vertex type (config line 41)"
  endif

  inquire (file=trim(filename_vertex_sym), exist=there)
  if (.not. there) then
    stop "Error: Two-Particle data file does not exist"
  endif

  end subroutine check_config

end module parameters_module
