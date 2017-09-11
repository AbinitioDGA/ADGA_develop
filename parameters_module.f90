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

  integer :: lines
  character(len=1), parameter :: cmnt = '#'
  character(len=1), parameter :: seperator = '='
  character(len=1), parameter :: multseperator = ' ' ! space
  character(len=150), allocatable :: file_temp(:), file_save(:)

  real(8), allocatable :: Umat(:,:,:,:)
  real(8), allocatable :: Udd(:),Vdd(:),Jdd(:)
  real(8), allocatable :: Udp(:),Vdp(:),Jdp(:)
  real(8), allocatable :: Upp(:),Vpp(:),Jpp(:)
  character(len=150), allocatable :: interaction(:)
  integer,allocatable :: interaction_mode(:)

end module parameters_module
