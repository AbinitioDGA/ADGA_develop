module parameters_module
  implicit none
  public

  integer   :: ounit ! output unit
  logical   :: verbose ! User defined verbose
  character(len=200) :: verbstr
  logical   :: debug ! User defined debug
  character(len=200) :: dbgstr

  !complex(kind=8),parameter :: ci=(0.d0,1.d0)
  !double precision,parameter :: pi=3.1415926535897932385d0

  ! different loop variables and information about box sizes
  integer :: iwstart,iwstop
  integer :: nqp, nkp_eom
  integer :: nqpphbar
  integer :: nkp, ndim, ndim2, maxdim, nineq
  integer :: nkpx,nkpy,nkpz,nqpx,nqpy,nqpz
  integer :: iwmax, iwbmax, iwfmax, iwbmax_small, iwfmax_small,n2iwb,n3iwf,n3iwb

  integer :: iwbcond ! number of bosonic frequencies for the conductivity
  integer :: iwfcond ! number of fermionic frequencies (internal loop) for the conductivity
  integer,allocatable :: q_data(:), k_data_eom(:)
  integer,allocatable :: q_data_phbar(:)
  integer,allocatable :: q_half_data_phbar(:)


  ! input data created from dmft
  integer, allocatable :: ndims(:,:)
  real(kind=8) :: mu, beta
  real(kind=8), allocatable :: iw_data(:), iwb_data(:), iwf_data(:)
  real(kind=8), allocatable :: k_data(:,:)
  complex(kind=8), allocatable :: u(:,:), u_tilde(:,:) ! u matrices in compound indices
  complex(kind=8), allocatable :: hk(:,:,:),dc(:,:)
  complex(kind=8), allocatable :: siw(:,:),giw(:,:)
  complex(kind=8), allocatable :: n_dmft(:), n_fock(:,:,:)


  ! for the optical conductivity
  double precision, allocatable :: hkder(:,:,:)
  integer :: iwcstart, iwcstop
  logical :: cond_dmftlegs
  logical :: do_cond_phbar
  logical :: do_cond_ph

  ! conductivity legs
  complex(kind=8), allocatable :: gkiwfull(:,:,:,:)
  complex(kind=8), allocatable :: gkiwfullbubble(:,:,:,:)

  ! input data created from previous ADGA in SC cycle
  complex(kind=8), allocatable :: skiw(:,:,:,:)

  ! output data created from the dga selfenergy
  complex(kind=8), allocatable :: n_dga(:), n_dga_k(:,:,:)

  ! run parameters and flags
  logical :: orb_sym
  logical :: do_eom,do_chi,do_cond,do_vq, do_ph
  logical :: do_chi_phbar, do_chi_ph
  logical :: extend_cond_bubble
  logical :: q_path_susc,q_vol,read_ext_hk,read_ext_u
  logical :: q_path_suscphbar
  logical :: k_path_eom
  logical :: external_chi_loc,external_threelegs
  logical :: exist_p
  logical :: susc_full_output
  logical :: text_output
  integer :: gzip_compression
  integer :: vertex_type
  integer,parameter :: full_g4=0,connected_g4=1,chi_g4=2
  integer :: summation_order ! order of geometric summation
  logical :: bse_inversion
  logical :: sc_mode

  ! filenames
  character(len=150) :: filename_vertex, filename_vertex_sym
  character(len=150) :: filename_chi_loc,filename_threelegs
  character(len=150) :: filename_umatrix, filename_vq
  character(len=150) :: filename_1p
  character(len=150) :: filename_hk, output_dir, filename_qdata
  character(len=150) :: filename_kdata, filename_qdata_susc
  character(len=150) :: output_filename, config_file

  character(len=150) :: filename_hkder
  character(len=150) :: filename_condlegs


  ! hdf5 iters
  character(len=150)  :: dmft_iter

  ! config file auxiliary variables
  integer :: lines
  character(len=1), parameter :: cmnt = '#'
  character(len=1), parameter :: seperator = '='
  character(len=1), parameter :: multseperator = ' ' ! space
  character(len=150), allocatable :: file_temp(:), file_save(:)

  ! umatrix parameters
  real(8), allocatable :: Umat(:,:,:,:) ! umatrix in band indices
  real(8), allocatable :: Udd(:),Vdd(:),Jdd(:)
  real(8), allocatable :: Udp(:),Vdp(:),Jdp(:)
  real(8), allocatable :: Upp(:),Vpp(:),Jpp(:)
  character(len=150), allocatable :: interaction(:)
  integer,allocatable :: interaction_mode(:)

  ! eigenvalues parameters and arrays
  logical :: calc_eigen
  integer :: number_eigenvalues
  integer :: number_eigenvectors

  contains

   subroutine dumpdata(ain,str)
   complex(KIND=8),intent(in)    :: ain(:,:)
   character(len=*),intent(in)   :: str
   integer                       :: i,j

   write(ounit,*) TRIM(ADJUSTL(str))
   write(ounit,'(1x,"Real:")')
   do i=1,size(ain,1)
      write(ounit,'(9999f12.8)') (dble(ain(i,j)),j=1,size(ain,2))
   enddo
   write(ounit,'(1x,"Imag:")')
   do i=1,size(ain,1)
      write(ounit,'(9999f12.8)') (dimag(ain(i,j)),j=1,size(ain,2))
   enddo

   return
   end subroutine dumpdata

end module parameters_module
