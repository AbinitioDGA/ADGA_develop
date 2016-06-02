program main

#ifdef MPI
  use mpi
#endif

  use hdf5
  use lapack_module
  use parameters_module

  implicit none


!  character(len=41), parameter :: filename = "SVO-Vertex-2016-02-19-Fri-17-50-41.hdf5" 
!  character(len=15), parameter :: filename_vertex = "vertex_sym.hdf5"

  integer(hid_t) :: file_id, file_vert_id, iw_id, iwb_id, iwf_id, siw_id, giw_id, plist_id, compound_id, r_id, k_id, hk_id, mu_id, dc_id, config_id, beta_id
  integer(hid_t) :: iw_space_id, iwb_space_id, iwf_space_id, siw_space_id, giw_space_id, k_space_id, hk_space_id, dc_space_id

  character(len=20) :: grpname_magn, grpname_dens, name_buffer
  character(len=30) :: name_buffer_dset
  integer(hid_t) :: grp_magn_id, grp_dens_id, nmembers, itype, dset_magn_id, dset_dens_id
  integer(hid_t) :: type_i_id, type_r_id
  integer(hsize_t), dimension(2) :: tmp_dims

  integer :: error, ierr
  integer(size_t) :: compound_size, type_sized 
  integer(hsize_t), dimension(1) :: iw_dims, iw_maxdims, iwb_dims, iwb_maxdims, iwf_dims, iwf_maxdims 
  integer(hsize_t), dimension(3) :: siw_dims, siw_maxdims, giw_dims, giw_maxdims
  integer(hsize_t), dimension(3) :: hk_dims, hk_maxdims
  integer(hsize_t), dimension(2) :: k_dims, k_maxdims, dc_dims, dc_maxdims
  integer(hsize_t), dimension(0) :: mu_dims, beta_dims
  
  double precision, allocatable :: iw_data(:), iwb_data(:)
  double precision, allocatable :: siw_data(:,:,:,:), giw_data(:,:,:,:)  
  complex(kind=8), allocatable :: siw(:,:)
  double precision, allocatable :: k_data(:,:), q_data(:,:) 
  double precision, allocatable :: hk_data(:,:,:,:)
  double precision, allocatable :: dc(:,:)
  complex(kind=8), allocatable :: hk(:,:,:)
  
  integer :: iw, ik, iq, ikq, iwf, iwb, iv, i, j, k, l, n, dum, dum1, nk, nq, ind_iwb, ind_grp, iwf1, iwf2
  integer :: imembers
  complex(kind=8), allocatable :: giw(:,:), gkiw(:,:)
  complex(kind=8), allocatable :: g4iw_magn(:,:,:,:,:,:), g4iw_dens(:,:,:,:,:,:) 
  double precision, allocatable :: tmp_r(:,:), tmp_i(:,:)
  complex(kind=8), allocatable :: chi0_loc(:,:), chi0_loc_inv(:,:,:), chi0(:,:,:), chi0_sum(:,:,:), chi_loc_slice_dens(:,:), sum_chi0_loc(:,:)
  complex(kind=8), allocatable ::  chi_loc_slice_magn(:,:), chi_loc_dens(:,:), chi_loc_magn(:,:), chi_loc(:,:)
  complex(kind=8), allocatable ::  chi_loc_magn_sum_left(:,:), chi_loc_dens_sum_left(:,:), gamma_dmft_dens(:,:), gamma_dmft_magn(:,:)
  complex(kind=8), allocatable :: chi_qw_dens(:,:,:),chi_qw_magn(:,:,:)
  integer, allocatable :: kq_ind(:,:), qw(:,:)
  complex(kind=8), allocatable :: interm2_magn(:,:), interm3_dens(:,:), interm3_magn(:,:), c(:,:)
  complex(kind=8), allocatable :: interm1(:,:), interm1_v(:,:), interm2_dens(:,:)

  real(kind=8 ):: start, finish, start1, finish1
  complex(kind=8) :: alpha, delta
  integer :: iqw, qwstart, qwstop
  logical :: update_chi_loc_flag!, small_freq_box, orb_sym
  integer :: b1, b2, b3, b4

  double precision :: u_value, kx, ky, kz
  double precision, allocatable :: u_tmp(:,:,:,:), u_tilde_tmp(:,:,:,:), hr(:,:), hi(:,:)
  complex(kind=8), allocatable :: u(:,:), u_tilde(:,:), u_work(:,:), m_work(:,:)
  complex(kind=8), allocatable :: m_tot(:,:), m_tot_array(:,:,:,:,:), gamma_loc(:,:), gamma_loc_sum_left(:,:), v(:,:)

  complex(kind=8), allocatable :: sigma(:,:,:,:), sigma_sum(:,:,:,:), sigma_loc(:,:,:)
  double precision :: beta
  integer(hsize_t) ::  inull
  integer :: iwb_zero, iband, ispin

  double precision :: iw_val, giw_r, giw_i, siw_r, siw_i


#ifdef MPI
  integer :: mpi_wrank
  integer :: mpi_wsize
  integer :: master
#endif


  ! read command line argument -> file name of config file
  if (iargc().eq.0 .or. iargc().gt.1) then
    write(*,*) 'The program has to be executed with exactly one argument. (Name of config file)'
    stop
  end if

  call read_config()
  call init()
  
#ifdef MPI
  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD,mpi_wrank,ierr)
  call MPI_comm_size(MPI_COMM_WORLD,mpi_wsize,ierr)
  master = 0
#endif

 
!  orb_sym = .true.

!  small_freq_box = .false.
!  iwfmax_small = 60
!  iwbmax_small = 15 

!  nk_frac = 1   !number of q-points in each direction nq=nk/nk_frac (cubic case assumed)

  !THE FOLLOWING PARAMETERS ARE READ FROM THE W2DYNAMICS OUTPUT-FILE:
  !iwmax or iw_dims(1)/2    number of fermionic Matsubara frequencies for single particle quantities 
  !nkp or hk_dims(3)    number of k-points in H(k)
  !ndim or hk_dims(1(2)) total number of bands in H(k) (d+p)
  !ndims or siw_dims(3)   number of d-bands

!#################################################################
 
  
!##################  READ W2DYNAMICS HDF5 OUTPUT FILE  #####################################

  call h5open_f(error)

  inull = 0

! Set dataset transfer property to preserve partially initialized fields during write/read to/from dataset with compound datatype (necessary?)
  call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  call h5pset_preserve_f(plist_id, .true., error)

  call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, error)
  call h5fopen_f(filename_vertex, h5f_acc_rdonly_f, file_vert_id, error) 

! create compound datatype for complex arrays:
  call h5tget_size_f(h5t_native_double, type_sized, error)
  compound_size = 2*type_sized
  call h5tcreate_f(h5t_compound_f, compound_size, compound_id, error)
  call h5tinsert_f(compound_id, "r", inull, h5t_native_double, error)
  call h5tinsert_f(compound_id, "i", type_sized, h5t_native_double, error)

  !complex type to write real and imaginary individually:
  call h5tcreate_f(h5t_compound_f, type_sized, type_r_id, error)
  call h5tinsert_f(type_r_id, "r", inull, h5t_native_double, error)
  call h5tcreate_f(h5t_compound_f, type_sized, type_i_id, error)
  call h5tinsert_f(type_i_id, "i", inull, h5t_native_double, error)

! read Matsubara frequencies iw (big range): 
  call h5dopen_f(file_id, ".axes/iw", iw_id, error)
  call h5dget_space_f(iw_id, iw_space_id, error)
  call h5sget_simple_extent_dims_f(iw_space_id, iw_dims, iw_maxdims, error)
  iwmax = iw_dims(1)/2
  allocate(iw_data(-iwmax:iwmax-1))
  call h5dread_f(iw_id, h5t_native_double, iw_data, iw_dims, error)
  call h5dclose_f(iw_id, error)

! read dimension of fermionic Matsubara frequencies iwf-g4 (small range):
  call h5dopen_f(file_vert_id, ".axes/iwf-g4", iwf_id, error)
  call h5dget_space_f(iwf_id, iwf_space_id, error)
  call h5sget_simple_extent_dims_f(iwf_space_id, iwf_dims, iwf_maxdims, error)
  iwfmax = iwf_dims(1)/2

  if (small_freq_box == .false.) iwfmax_small = iwfmax
  if (iwfmax_small .gt. iwfmax) then
     write(*,*) 'Error: Maximum number of fermionic frequencies =', iwfmax
  endif
  write(*,*) 'iwfmax=', iwfmax, 'iwfmax_small=', iwfmax_small
! read bosonic Matsubara frequencies iwb-g4:
  call h5dopen_f(file_vert_id, ".axes/iwb-g4", iwb_id, error)
  call h5dget_space_f(iwb_id, iwb_space_id, error)
  call h5sget_simple_extent_dims_f(iwb_space_id, iwb_dims, iwb_maxdims, error)
  iwbmax = iwb_dims(1)/2
  allocate(iwb_data(-iwbmax:iwbmax))
  call h5dread_f(iwb_id, h5t_native_double, iwb_data, iwb_dims, error)
  call h5dclose_f(iwb_id, error)
  ! "Find" the zero 
  iwb_zero = 0

  if (small_freq_box == .false.) iwbmax_small = iwbmax
  if (iwbmax_small .gt. iwbmax) then
     write(*,*) 'Error: Maximum number of bosonic frequencies =', iwbmax
  endif
  write(*,*)'iwbmax=',iwbmax, 'iwbmax_small=', iwbmax_small
 

! read siw:
  call h5dopen_f(file_id, "stat-001/ineq-001/siw/value", siw_id, error)
  call h5dget_space_f(siw_id, siw_space_id, error)
  call h5sget_simple_extent_dims_f(siw_space_id, siw_dims, siw_maxdims, error)
  ndims = siw_dims(3)
  allocate(siw_data(2,-iwmax:iwmax-1,siw_dims(2),siw_dims(3))) !indices: real/imag iw spin band
  call h5dread_f(siw_id, compound_id, siw_data, siw_dims, error)
  allocate(siw(-iwmax:iwmax-1,siw_dims(3)))

  !paramagnetic:
  siw = 0.d0
  siw(:,:) = siw_data(1,:,1,:)+siw_data(1,:,2,:)+ci*siw_data(2,:,1,:)+ci*siw_data(2,:,2,:)
  siw = siw/2.d0

  call h5dclose_f(siw_id, error)

  if (orb_sym) then

     ! enforce orbital symmetry:
     do iband=2,ndims
        siw(:,1) = siw(:,1)+siw(:,iband)
     enddo

     do iband=1,ndims
        siw(:,iband) = siw(:,1)
     enddo
     siw = siw/dble(ndims)

  endif

  ! test siw:
  open(34, file="siw.dat", status='unknown')
  do iw=-iwmax,iwmax-1   
     write(34,'(100F12.6)')iw_data(iw), (real(siw(iw,i)),aimag(siw(iw,i)), i=1,ndims)
  enddo
  close(34)


  ! read giw:
  call h5dopen_f(file_id, "stat-001/ineq-001/giw/value", giw_id, error)
  call h5dget_space_f(giw_id, giw_space_id, error)
  call h5sget_simple_extent_dims_f(giw_space_id, giw_dims, giw_maxdims, error)
  allocate(giw_data(2,-iwmax:iwmax-1,giw_dims(2),giw_dims(3))) !indices: real/imag iw spin band
  call h5dread_f(giw_id, compound_id, giw_data, giw_dims, error)
  allocate(giw(-iwmax:iwmax-1,giw_dims(3)))

  !paramagnetic:
  giw = 0.d0
  giw(:,:) = giw_data(1,:,1,:)+giw_data(1,:,2,:)+ci*giw_data(2,:,1,:)+ci*giw_data(2,:,2,:)
  giw = giw/2.d0

  call h5dclose_f(giw_id, error)
  
  if (orb_sym) then

     ! enforce orbital symmetry:
     do iband=2,ndims
        giw(:,1) = giw(:,1)+giw(:,iband)
     enddo

     do iband=1,ndims
        giw(:,iband) = giw(:,1)
     enddo
     giw = giw/dble(ndims)

  endif 

! test giw:
  open(54, file="giw.dat", status='unknown')
  do iw=-iwmax,iwmax-1
     write(54,'(100F12.6)')iw_data(iw), (real(giw(iw,1)),aimag(giw(iw,1)),i=1,ndims)
  enddo
  close(54)

 
! read k-points:
  call h5dopen_f(file_id, ".axes/k-points", k_id, error)
  call h5dget_space_f(k_id, k_space_id, error)
  call h5sget_simple_extent_dims_f(k_space_id, k_dims, k_maxdims, error)
  nkp = k_dims(2)
  allocate(k_data(k_dims(1),k_dims(2))) !indices: 3 ik
  call h5dread_f(k_id, h5t_native_double, k_data, k_dims, error)
  call h5dclose_f(k_id, error)

! write k-points:
!  open(37, file='k_points.dat', status='unknown')
!  do ik=1,100
!    write(37,'(100F12.6)') k_data(2,ik), k_data(3,ik)
!  enddo
!  close(37)

! read Hamiltonian H(k):
  call h5dopen_f(file_id, "start/hk/value", hk_id, error)
  call h5dget_space_f(hk_id, hk_space_id, error)
  call h5sget_simple_extent_dims_f(hk_space_id, hk_dims, hk_maxdims, error)
  ndim = hk_dims(1)
  allocate(hk_data(2,hk_dims(1),hk_dims(2),hk_dims(3)))
  call h5dread_f(hk_id, compound_id, hk_data, hk_dims, error)
  allocate(hk(hk_dims(1),hk_dims(2),hk_dims(3))) !indices: band band ik
  hk = 0.d0
  hk(:,:,:) = hk_data(1,:,:,:)+ci*hk_data(2,:,:,:)
  call h5dclose_f(hk_id, error)

! test hk:
!  open(34, file="hk.dat", status='unknown')
!  do ik=1,hk_dims(3)
!     write(34,*)k_data(:,ik)
!     do i=1,hk_dims(2)
!        write(34,'(100F12.6)')hk(:,i,ik)
!     enddo
!  enddo
!  close(34)


! read chemical potential:
  call h5dopen_f(file_id, "stat-001/mu/value", mu_id, error)
  call h5dread_f(mu_id, h5t_native_double, mu, mu_dims, error)
  call h5dclose_f(mu_id, error)

! read double counting:
  call h5dopen_f(file_id, "stat-001/ineq-001/dc/value", dc_id, error)
  call h5dget_space_f(dc_id, dc_space_id, error)
  call h5sget_simple_extent_dims_f(dc_space_id, dc_dims, dc_maxdims, error)
  allocate(dc(dc_dims(1),dc_dims(2))) !indices: spin band
  call h5dread_f(dc_id, h5t_native_double, dc, dc_dims, error)
  call h5dclose_f(dc_id, error)

! read inverse temperature beta:
  call h5gopen_f(file_id, ".config", config_id, error)
  call h5aopen_f(config_id, "general.beta", beta_id, error)
  call h5aread_f(beta_id, h5t_native_double, beta, beta_dims, error)
  call h5aclose_f(beta_id, error)
  call h5gclose_f(config_id,error)
 
  call h5fclose_f(file_vert_id, error)
  call h5fclose_f(file_id, error)
  
  call h5close_f(error)

  write(*,*) 'beta=', beta

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
  
  
!################################################################################################

! compute local single-particle Greens function:
  call cpu_time(start)

!  allocate(giw(-iwmax:iwmax-1),ndim) 
!  call get_giw(iw_data, hk, siw, dc, giw)
!  call cpu_time(finish)

  !write(*,*)'computing giw:', finish-start
 
! test giw:
  open(35, file="giw_calc.dat", status='unknown')
  do iw=-iwmax,iwmax-1
     write(35,'(100F12.6)') iw_data(iw), (real(giw(iw,i)),aimag(giw(iw,i)),i=1,1)
  enddo
  close(35)


  !allocate(chi0_loc(ndim2,ndim2)) !only for test reasons
  allocate(sum_chi0_loc(ndim2,ndim2)) !test

  allocate(chi0_loc_inv(ndim2,ndim2,-iwfmax:iwfmax-1))
  allocate(chi0(ndim2,ndim2,-iwfmax_small:iwfmax_small-1))
  allocate(chi0_sum(ndim2,ndim2,-iwfmax_small:iwfmax_small-1)) 
  
  allocate(interm1(ndim2,ndim2))
  allocate(interm1_v(ndim2,ndim2))
  
  allocate(gamma_dmft_dens(ndim2,maxdim), gamma_dmft_magn(ndim2,maxdim))
  allocate(chi_loc_slice_dens(ndim2,maxdim))
  allocate(chi_loc_slice_magn(ndim2,maxdim))
  allocate(chi_loc_dens_sum_left(ndim2,maxdim))
  allocate(chi_loc_magn_sum_left(ndim2,maxdim))
  allocate(chi_loc_magn(maxdim,maxdim))
  allocate(chi_loc_dens(maxdim,maxdim))
  allocate(c(ndim2,maxdim))

  allocate(interm3_dens(ndim2,maxdim))!allocate somewhere else?
  allocate(interm3_magn(ndim2,maxdim))
  allocate(gamma_loc_sum_left(ndim2,maxdim))

  allocate(gkiw(ndim,ndim))
  allocate(sigma(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp))

  allocate(v(ndim2,ndim2))


  !define q-grid:
  nk = nkp**(1./3.)  

  if(mod(nk,nk_frac).ne.0)then
     stop 'mismatch between k- and q-grid!'
   endif

  nq = nk/nk_frac
  nqp = nq**3  
 
  allocate(q_data(3,nqp))

  i1=0
  do i=0,nq-1
     do j=0,nq-1
        do k=0,nq-1
           i1 = i1+1
           q_data(1,i1) = i*1.d0/(dble(nq))
           q_data(2,i1) = j*1.d0/(dble(nq)) 
           q_data(3,i1) = k*1.d0/(dble(nq))
        enddo
     enddo
  enddo

  !search for k+q - index:
  call cpu_time(start)
  allocate(kq_ind(nkp,nqp))
  call index_kq(k_data, q_data, kq_ind) !assumes cubic case
  call cpu_time(finish)
  !write(*,*)'finding k-q index:', finish-start



!##################### parallel code ##################################################

  write(*,*)'nqp=', nqp !test

!define qw compound index for mpi:
  allocate(qw(2,nqp*(2*iwbmax+1)))
  i1=0
  do iwb=-iwbmax_small,iwbmax_small
     do iq=1,nqp
        i1 = i1+1
        qw(1,i1) = iwb
        qw(2,i1) = iq
     enddo
  enddo

  allocate(chi_qw_dens(ndim2,ndim2,nqp*(2*iwbmax_small+1)),chi_qw_magn(ndim2,ndim2,nqp*(2*iwbmax_small+1)))
  chi_qw_dens=0.d0
  chi_qw_magn=0.d0


!distribute the qw compound index:
#ifdef MPI
  qwstop = 0
  do i=0,mpi_wrank
     j = (nqp*(2*iwbmax_small+1) - qwstop)/(mpi_wsize-i)
     qwstart = qwstop + 1
     qwstop = qwstop + j
  enddo
#else 
  qwstart = 1
  qwstop = nqp*(2*iwbmax_small+1)
#endif

#ifdef MPI
write(*,*)'rank=',mpi_wrank, 'qwstart=', qwstart, 'qwstop=', qwstop
start = mpi_wtime()
#endif



  sigma = 0.d0

  open(55, file="chi0_loc_sum.dat", status='unknown')

  iwb = iwbmax_small+3
  do iqw=qwstart,qwstop
     update_chi_loc_flag = qw(1,iqw) .ne. iwb

     iq = qw(2,iqw)
     iwb = qw(1,iqw)

     !to be done here: read nonlocal interaction v and go into compound index
     v = 0.d0

     !update chi_loc only if iwb is different than the previous one:
     if(update_chi_loc_flag) then    
  
        !call cpu_time(start)
        do iwf=-iwfmax,iwfmax-1

           ! compute local bubble chi0_loc^{-1}(i1,i2)(orbital compound index i1,i2):
           call get_chi0_loc_inv(beta, iwf, iwb, giw, chi0_loc_inv(:,:,iwf))

           ! test chi0_loc_inv:
           !  open(36, file="chi0_loc_inv.dat", status='unknown')
           !  do i1=1,ndim2
           !     do i2=1,ndim2
           !        write(36,*)i1,i2,chi0_loc_inv(i1,i2,iwf)
           !     enddo
           !  enddo

        enddo
        !call cpu_time(finish)
        !write(*,*) 'calculating chi0_loc_inv', finish-start
        
        !call cpu_time(start)

        !get iwb-slice of w2dynamics vertex:
        allocate(g4iw_magn(ndims, ndims, -iwfmax:iwfmax-1, ndims, ndims, -iwfmax:iwfmax-1))
        allocate(g4iw_dens(ndims, ndims, -iwfmax:iwfmax-1, ndims, ndims, -iwfmax:iwfmax-1))
        allocate(tmp_r(-iwfmax:iwfmax-1, -iwfmax:iwfmax-1))
        allocate(tmp_i(-iwfmax:iwfmax-1, -iwfmax:iwfmax-1))
        tmp_dims = (/2*iwfmax, 2*iwfmax/)

        g4iw_magn = 0.d0
        g4iw_dens = 0.d0
 
        ind_iwb = iwb+iwbmax
        write(grpname_magn, '(A5,(I5.5),A1,(I5.5))'), "magn/", ind_iwb
        write(grpname_dens, '(A5,(I5.5),A1,(I5.5))'), "dens/", ind_iwb

        call h5fopen_f(filename_vertex, h5f_acc_rdonly_f, file_vert_id, error)
        call h5gopen_f(file_vert_id, grpname_magn, grp_magn_id, error) 
        call h5gopen_f(file_vert_id, grpname_dens, grp_dens_id, error)

        call h5gn_members_f(file_vert_id, grpname_magn, nmembers, error)
        
        do imembers=0,nmembers-1
           
           call h5gget_obj_info_idx_f(file_vert_id, grpname_magn, imembers, name_buffer, itype, error)
           read(name_buffer,'(I5.5)')ind_grp
           
           call index2component_band(ndims, ind_grp, b1, b2, b3, b4)
           
           write(name_buffer_dset, '(A5,(I5.5),A1,(I5.5),A6)')"magn/", ind_iwb, "/", ind_grp, "/value"
           call h5dopen_f(file_vert_id, name_buffer_dset, dset_magn_id, error)
           call h5dread_f(dset_magn_id, type_r_id, tmp_r, tmp_dims, error)
           call h5dread_f(dset_magn_id, type_i_id, tmp_i, tmp_dims, error)

           g4iw_magn(b1,b2,:,b3,b4,:) = (tmp_r(:,:)+ci*tmp_i(:,:))*beta

           call h5dclose_f(dset_magn_id, error)

           write(name_buffer_dset, '(A5,(I5.5),A1,(I5.5),A6)')"dens/", ind_iwb, "/", ind_grp, "/value"
           call h5dopen_f(file_vert_id, name_buffer_dset, dset_dens_id, error)
           call h5dread_f(dset_dens_id, type_r_id, tmp_r, tmp_dims, error)
           call h5dread_f(dset_dens_id, type_i_id, tmp_i, tmp_dims, error)

           g4iw_dens(b1,b2,:,b3,b4,:) = (tmp_r(:,:)+ci*tmp_i(:,:))*beta
           
           call h5dclose_f(dset_dens_id, error)
           
           
        enddo

        call h5gclose_f(grp_dens_id, error)
        call h5gclose_f(grp_magn_id, error)
        call h5fclose_f(file_vert_id, error)

     
        !compute chi_loc (go into compound index and subtract straight term):
        chi_loc_magn = 0.d0
        chi_loc_dens = 0.d0

        !open(73, file="chi_loc_magn_before.dat", status='unknown') !test
        !open(74, file="chi_loc_magn_after.dat", status='unknown') !test
        
        i2 = 0
        do iwf2=-iwfmax_small,iwfmax_small-1
           do l=1,ndim
              do k=1,ndim
                 i2 = i2+1
                 i1 = 0
                 do iwf1=-iwfmax_small,iwfmax_small-1
                    do i=1,ndim
                       do j=1,ndim
                          i1 = i1+1
                          chi_loc_magn(i1,i2) = g4iw_magn(i,j,iwf1,k,l,iwf2)
                          chi_loc_dens(i1,i2) = g4iw_dens(i,j,iwf1,k,l,iwf2)

                           !if (iwf2==0 .and. iwb==1)then !test
                           !     write(73,'(100F12.6)') iw_data(iwf1), chi_loc_magn(i1,i2)
                           !  endif

                          !straight term is subtracted (twice) only in the dens channel and only for iw=0:
                          if((iwb .eq. iwb_zero) .and. i==j .and. k==l)then
                             chi_loc_dens(i1,i2) = chi_loc_dens(i1,i2)-2.d0*beta*giw(iwf1,i)*giw(iwf2,l) 
                          endif

                         !if (iwf2==0 .and. iwb==1)then !test
                         !       write(74,'(100F12.6)') iw_data(iwf1), chi_loc_magn(i1,i2)
                         !endif

                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo

        !close(73) !test
        !close(74)

        deallocate(g4iw_magn, g4iw_dens, tmp_r, tmp_i)
       
        !time reversal symmetry (which is simply a transpose in our compound index)
        do i1=1,maxdim
           do i2=i1+1,maxdim

              chi_loc_magn(i1,i2) = 0.5d0*(chi_loc_magn(i1,i2)+chi_loc_magn(i2,i1))
              chi_loc_magn(i2,i1) = chi_loc_magn(i1,i2)

              chi_loc_dens(i1,i2) = 0.5d0*(chi_loc_dens(i1,i2)+chi_loc_dens(i2,i1))
              chi_loc_dens(i2,i1) = chi_loc_dens(i1,i2)

           enddo
        enddo

        

        !compute chi_loc*chi0_loc_inv (chi0_loc_inv is diagonal in compound index):
        !use blas-routine instead?
        dum = 0
        do iwf=-iwfmax_small,iwfmax_small-1
           do i2=1,ndim2
              do i1=1,maxdim
           
                 chi_loc_magn(i1,i2+dum*ndim2) = chi_loc_magn(i1,i2+dum*ndim2)*chi0_loc_inv(i2,i2,iwf)
                 chi_loc_dens(i1,i2+dum*ndim2) = chi_loc_dens(i1,i2+dum*ndim2)*chi0_loc_inv(i2,i2,iwf)

              enddo
           enddo
           dum = dum+1
        enddo

           
        !sum up the left fermionic frequency of chi_loc (quantity needed afterwards)
        chi_loc_magn_sum_left = 0.d0
        chi_loc_dens_sum_left = 0.d0

        do i1=1,ndim2
           do dum=0,2*iwfmax_small-1
              
              chi_loc_magn_sum_left(i1,:) = chi_loc_magn_sum_left(i1,:)+chi_loc_magn(i1+dum*ndim2,:)
              chi_loc_dens_sum_left(i1,:) = chi_loc_dens_sum_left(i1,:)+chi_loc_dens(i1+dum*ndim2,:)
              
           enddo
        enddo
        
        !call cpu_time(finish)
        gamma_dmft_dens = 0.d0
        gamma_dmft_magn = 0.d0

        do i1=1,ndim2
           do dum=0,2*iwfmax_small-1
              gamma_dmft_dens(i1,i1+dum*ndim2) = chi_loc_dens_sum_left(i1,i1+dum*ndim2)-1.d0
           enddo
        enddo

        do i1=1,ndim2
           do dum=0,2*iwfmax_small-1
              gamma_dmft_magn(i1,i1+dum*ndim2) = chi_loc_magn_sum_left(i1,i1+dum*ndim2)-1.d0
           enddo
        enddo
        !write(*,*)'update local quantities:', finish-start
      
     endif


     dum1 = 0 !index for second matrix multiplication to get interm2
     
     allocate(interm2_dens(maxdim,maxdim))
     allocate(interm2_magn(maxdim,maxdim))
     allocate(gamma_loc(maxdim,maxdim))
     interm2_dens = 0.d0
     interm2_magn = 0.d0
     gamma_loc = 0.d0

     !chi_loc_magn = 0.d0 !only for test reason (when chi_loc_magn is replaced with chi0_loc) 
     sum_chi0_loc = 0.d0  !only to test summed up local bubble

     !call cpu_time(start)
     chi0_sum=0.d0
     do iwf=-iwfmax_small,iwfmax_small-1
       
        ! compute k-summed (but still q-dependent) bubble chi0(i1,i2):
        !replace chi0_sum with chi0_loc for dmft test!!!!!!
!        chi0_sum = 0.d0
        do ik=1,nkp
           ikq = kq_ind(ik,iq) !Index of G(k+q)
           call get_chi0(beta, ik, ikq, iwf, iwb, iw_data, siw, hk, dc, chi0(:,:,iwf)) 

           !call get_chi0_loc(beta, iwf,iwb,giw,chi0_sum)
           
           chi0_sum(:,:,iwf) = chi0_sum(:,:,iwf)+chi0(:,:,iwf)
        enddo
        chi0_sum(:,:,iwf) = chi0_sum(:,:,iwf)/dble(nkp)

        !test local bubble:
        sum_chi0_loc = sum_chi0_loc + chi0_sum(:,:,iwf)

        
        ! compute intermediate quantity chi0*chi0_loc_inv (chi0_loc_inv is diagonal):
        do j=1,ndim2
           do i=1,ndim2

              interm1(i,j) = chi0_sum(i,j,iwf)*chi0_loc_inv(j,j,iwf) !remove one index for chi0_loc_inv
          
           enddo
        enddo
        

        ! compute part containing nonlocal interaction v: chi0*v
        ! do it with matrix multiplication routine?
        interm1_v = 0.d0

        do i=1,ndim2
           do j=1,ndim2
              do k=1,ndim2
                    
                 interm1_v(i,j) = interm1_v(i,j)+chi0_sum(i,k,iwf)*v(k,j)
                 
              enddo
           enddo
        enddo
                    
        interm1_v = interm1_v/(beta**2)

        !get horizontal slice of chi_loc for matmul (interm2):
        chi_loc_slice_magn(:,:) = chi_loc_magn((dum1*ndim2)+1 : (dum1+1)*ndim2, :)
        chi_loc_slice_dens(:,:) = chi_loc_dens((dum1*ndim2)+1 : (dum1+1)*ndim2, :)
        
             
!########################################################################

        !test, replace chi_loc by bubble:
        !chi_loc_slice_magn = 0.d0
      
 
        !call get_chi0_loc(beta, iwf, iwb, giw, chi0_loc)

        !do i=1,ndim2
        !   do j=1,ndim2
        !      chi_loc_slice_magn(i,j+(dum1*ndim2)) = chi0_loc(i,j)
        !   enddo
        !enddo

        !do i=1,ndim2
        !   do j=1,ndim2
        !      chi_loc_magn_sum_left(i,j+(dum1*ndim2)) = chi0_loc(i,j)
        !   enddo
        !enddo

        !do i=1,ndim2
        !   do j=1,ndim2
        !      chi_loc_magn(i+(dum1*ndim2),j+(dum1*ndim2)) = chi0_loc(i,j)
        !   enddo
        !enddo

!########################################################################
     
              
        ! compute intermediate quantity (chi0*chi0_loc_inv)*chi_loc:
        c = 0.d0
        alpha = 1.d0
        delta = 0.d0
        call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, interm1, ndim2, chi_loc_slice_dens, ndim2, delta, c, ndim2)
             
        do i=1,ndim2
           do j=1,maxdim
              interm2_dens(i+(dum1*ndim2),j) = c(i,j)
           enddo
        enddo
        
        !same for the magn channel:
        c = 0.d0
        alpha = 1.d0
        delta = 0.d0
        call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, interm1, ndim2, chi_loc_slice_magn, ndim2, delta, c, ndim2)
             
        do i=1,ndim2
           do j=1,maxdim
              interm2_magn(i+(dum1*ndim2),j) = c(i,j)
           enddo
        enddo
            
             
        ! subtract interm1(diagonal in iwf) from interm2:
        ! (for 'chi0 test' replace interm1 with chi0_sum)
        do i=1,ndim2
           do j=1,ndim2
              interm2_dens(i+(dum1*ndim2),j+(dum1*ndim2)) = interm2_dens(i+(dum1*ndim2),j+(dum1*ndim2))-interm1(i,j)
              interm2_magn(i+(dum1*ndim2),j+(dum1*ndim2)) = interm2_magn(i+(dum1*ndim2),j+(dum1*ndim2))-interm1(i,j)
              !interm2_magn(i+(dum1*ndim2),j+(dum1*ndim2)) = interm2_magn(i+(dum1*ndim2),j+(dum1*ndim2))-chi0_sum(i,j)
           enddo
        enddo

        !store interm2 here since this quantity will be used in the EOM 
        do i=1,ndim2
           do j=1,maxdim
              gamma_loc(i+(dum1*ndim2),j) = interm2_dens(i+(dum1*ndim2),j)
           enddo
        enddo

        ! compute the part of interm2 containing v: (chi0*v)*chi_loc_sum_left
        c = 0.d0
        alpha = 1.d0
        delta = 0.d0
        call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, interm1_v, ndim2, chi_loc_dens_sum_left, ndim2, delta, c, ndim2)
             
        ! add it to interm2:
        do i=1,ndim2
           do j=1,maxdim
              interm2_dens(i+(dum1*ndim2),j) = interm2_dens(i+(dum1*ndim2),j)+c(i,j)
           enddo
        enddo

        !same for the magnetic channel:
        c = 0.d0
        alpha = 1.d0
        delta = 0.d0
        call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, interm1_v, ndim2, chi_loc_magn_sum_left, ndim2, delta, c, ndim2)
             
        do i=1,ndim2
           do j=1,maxdim
              interm2_magn(i+(dum1*ndim2),j) = interm2_magn(i+(dum1*ndim2),j)+c(i,j)
           enddo
        enddo
        
                    
        dum1 = dum1+1
     enddo !iwf

     !call cpu_time(finish)
     !write(*,*)'doing the iwf-loop:', finish-start

     !call cpu_time(start)
     
     !sum up left iwf-index of interm2_dens and store the quantity (is afterwards directly used in the EOM):
     gamma_loc_sum_left = 0.d0

     do i1=1,ndim2
        do dum=0,2*iwfmax_small-1
              
           gamma_loc_sum_left(i1,:) = gamma_loc_sum_left(i1,:)+gamma_loc(i1+dum*ndim2,:)
              
        enddo
     enddo

     deallocate(gamma_loc)


     !construct quantity that is inverted right afterwards:
     interm2_dens = chi_loc_dens-interm2_dens
     interm2_magn = chi_loc_magn-interm2_magn

     !call cpu_time(finish)
     !write(*,*)'do steps in between:', finish-start

#ifdef MPI
  start1 = mpi_wtime()
#endif

     call inverse_matrix(interm2_dens)
     call inverse_matrix(interm2_magn)

#ifdef MPI
  finish1 = mpi_wtime()
  !if(mpi_wrank .eq. master) then
  !  write(*,*)'doing inversion:',finish1-start1
  !endif
#endif
       

     !check:
     !do i=1,maxdim
     !   do j=1,maxdim
     !      write(*,*) i, j, interm2_dens(i,j)
     !   enddo
     !enddo

     !call cpu_time(finish)
     !write(*,*)'doing inversion:', finish-start
     
     !call cpu_time(start)

     !do matrix multiplication chi_loc_sum_left*interm2
     interm3_dens = 0.d0
     interm3_magn = 0.d0

     alpha = 1.d0
     delta = 0.d0
     call zgemm('n', 'n', ndim2, maxdim, maxdim, alpha, chi_loc_dens_sum_left &
             , ndim2, interm2_dens, maxdim, delta, interm3_dens, ndim2)

     alpha = 1.d0
     delta = 0.d0
     call zgemm('n', 'n', ndim2, maxdim, maxdim, alpha, chi_loc_magn_sum_left &
               , ndim2, interm2_magn, maxdim, delta, interm3_magn, ndim2)



     deallocate(interm2_dens, interm2_magn)

     !output for 'chi0 test' (unitary matrix in orbital compound index):
     !do i=1,ndim2
     !   do j=i,maxdim
     !      write(*,*) i, j, interm3_magn(i,j)
     !   enddo
     !enddo
     !stop




! Calculation of q dependent susceptibility by multiplication with chi0
     call calc_chi_qw(chi_qw_dens(:,:,iqw),interm3_dens,chi0_sum)
     call calc_chi_qw(chi_qw_magn(:,:,iqw),interm3_magn,chi0_sum)

! from here on: equation of motion     
     call calc_eom(interm3_dens,interm3_magn,gamma_dmft_dens,gamma_dmft_magn,gamma_loc_sum_left,sigma,kq_ind,iwb,iq,iw_data)

     !call cpu_time(finish)
     !write(*,*)'equation of motion:', finish-start

    
  enddo !iqw

#ifdef MPI
  finish = mpi_wtime()
  if (mpi_wrank .eq. master) then
    write(*,*)'doing one qw-loop:', finish-start
  endif
#endif
  
#ifdef MPI
  allocate(sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp))
  allocate(sigma_loc(ndim, ndim, -iwfmax_small:iwfmax_small-1))
 
  call MPI_reduce(sigma, sigma_sum, ndim*ndim*2*iwfmax_small*nkp, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)

  sigma_sum = -sigma_sum/(beta*nqp)
 
   ! local contribution is replaced by the DMFT self energy for better asymptotics
  do ik=1,nkp
     do iwf=-iwfmax_small,iwfmax_small-1
        do iband=1,ndim
           sigma_sum(iband, iband, iwf, ik) = sigma_sum(iband, iband, iwf, ik) + siw(iwf, iband)
        enddo
     enddo
  enddo

  sigma_loc = 0.d0
  do ik=1,nkp
    sigma_loc(:,:,:) = sigma_loc(:,:,:)+sigma_sum(:,:,:,ik)
  enddo
  sigma_loc = sigma_loc/dble(nkp)

  call MPI_finalize(ierr)

  !TEST:
  if (mpi_wrank .eq. master) then
    open(34, file="siw_0_0_0_nostraight.dat", status='unknown')
    open(35, file="siw_0_0_0.5_nostraight.dat", status='unknown')
    open(36, file="siw_0_0.5_0.5_nostraight.dat", status='unknown')
    open(37, file="siw_loc_nostraight.dat", status='unknown')
    open(38, file="siw_iwf_0_nostraight.dat", status='unknown')
    open(39, file="siw_0.5_0.5_0.5_nostraight.dat", status='unknown')

    do ik=1,100
       write(38,'(100F12.6)') k_data(2,ik), k_data(3,ik), (real(sigma_sum(i,i,0,ik)), aimag(sigma_sum(i,i,0,ik)), i=1,3)
    enddo 

    do iwf=-iwfmax_small,iwfmax_small-1
       write(34,'(100F12.6)')iw_data(iwf), (real(sigma_sum(i,i,iwf,1)), aimag(sigma_sum(i,i,iwf,1)), i=1,3)
       write(35,'(100F12.6)')iw_data(iwf), (real(sigma_sum(i,i,iwf,6)), aimag(sigma_sum(i,i,iwf,6)),i=1,3)
       write(36,'(100F12.6)')iw_data(iwf), (real(sigma_sum(i,i,iwf,56)), aimag(sigma_sum(i,i,iwf,56)),i=1,3)
       write(37,'(100F12.6)')iw_data(iwf), (real(sigma_loc(i,i,iwf)), aimag(sigma_loc(i,i,iwf)),i=1,3)
       write(39,'(100F12.6)')iw_data(iwf), (real(sigma_sum(i,i,iwf,556)),aimag(sigma_sum(i,i,iwf,556)),i=1,3)

    enddo

    close(34)
    close(35)
    close(36)
    close(37)
    close(38)
    close(39)

    call output_chi_qw(chi_qw_dens,q_data,qw,'Output/chi_qw_dens.dat')
    call output_chi_qw(chi_qw_magn,q_data,qw,'Output/chi_qw_magn.dat')

      
  endif
  
#endif

     
end program main


! subroutine to output susceptibility
! for now not adapted to more than 1 band
! but i wrote it already somewhere for more bands...
subroutine output_chi_qw(chi_qw,q_data,qw,filename_output)
  use parameters_module
  implicit none
  character(len=*) :: filename_output
  real*8 :: q_data(3,nqp)
  complex(kind=8) :: chi_qw(ndim2,ndim2,nqp*(2*iwbmax_small+1))
  integer :: iwb,iq,qw(2,nqp*(2*iwbmax+1))
    open(unit=10,file=filename_output)
    write(10,*) '#iwb  ','iq  ','      (q)      ','chi_qw'
    do i1=1,nqp*(2*iwbmax_small+1)
      iq = qw(2,i1)
      iwb = qw(1,i1)
      write(10,'(I5,2X,I5,2X,5(E14.7E2,2X))') iwb,iq,q_data(:,iq),real(chi_qw(1,1,i1),8),dimag(chi_qw(1,1,i1))
      if (mod(i1,nqp).eq.0) then
        write(10,*) ' '
      end if
    end do
    close(10)
end subroutine output_chi_qw


subroutine calc_eom(interm3_dens,interm3_magn,gamma_dmft_dens,gamma_dmft_magn,gamma_loc_sum_left,sigma,kq_ind,iwb,iq,iw_data)
  use parameters_module
  use lapack_module

  implicit none
  complex(kind=8), allocatable :: u(:,:), u_tilde(:,:), u_work(:,:), m_work(:,:)
  complex(kind=8) :: interm3_dens(ndim2,maxdim),interm3_magn(ndim2,maxdim)
  complex(kind=8) :: gamma_dmft_dens(ndim2,maxdim), gamma_dmft_magn(ndim2,maxdim)
  complex(kind=8) :: gamma_loc_sum_left(ndim2,maxdim)
  complex(kind=8) :: alpha, delta
  complex(kind=8) :: v(ndim2,ndim2),sigma(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp)
  complex(kind=8), allocatable :: m_tot_array(:,:,:,:,:),m_tot(:,:)
  integer :: dum,i,j,iwf,iwb,iwf2,l,k,ik,iq,ikq,kq_ind(nkp,nqp)
  real*8 :: iw_data(-iwmax:iwmax-1)
  complex(kind=8) :: siw(-iwmax:iwmax-1,ndims)
  complex(kind=8) :: hk(ndim,ndim,nkp)
  double precision :: dc(2,ndim)
  complex(kind=8) :: gkiw(ndim,ndim)


  !subtract -1 in the diagonal of the orbital blocks:
  do i1=1,ndim2
     do dum=0,2*iwfmax_small-1
        
        interm3_dens(i1,i1+dum*ndim2) = interm3_dens(i1,i1+dum*ndim2)-1.d0
        interm3_magn(i1,i1+dum*ndim2) = interm3_magn(i1,i1+dum*ndim2)-1.d0

     enddo
  enddo
  
  

  !call cpu_time(finish)
  !write(*,*)'doing matrix multiplication:', finish-start


  
  !call cpu_time(start)

  !EQUATION OF MOTION:
  allocate(m_work(ndim2, maxdim), m_tot(ndim2, maxdim), u_work(ndim2,ndim2))
  m_work = 0.d0
  m_tot = 0.d0
  u_work = 0.d0

  !density part:
  u_work = v + u - 0.5d0*u_tilde

  !subtract dmft part (dmft self energy will be added in the end):
  interm3_dens = interm3_dens - gamma_dmft_dens
  
  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, interm3_dens, ndim2, delta, m_work, ndim2)

  m_tot = m_work

  !magnetic part:
  u_work = -1.5d0*u_tilde
  m_work = 0.d0

  !subtract dmft part:
  interm3_magn = interm3_magn - gamma_dmft_magn

  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, interm3_magn, ndim2, delta, m_work, ndim2)

  m_tot = m_tot + m_work

  !local part:
  u_work = v + u
  m_work = 0.d0

  gamma_loc_sum_left = gamma_loc_sum_left - gamma_dmft_dens
  
  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, gamma_loc_sum_left, ndim2, delta, m_work, ndim2)

  m_tot = m_tot - m_work

  !break up the compound index:
  allocate(m_tot_array(ndim,ndim,ndim,ndim,-iwfmax_small:iwfmax_small))
  m_tot_array = 0.d0

  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1 = i1+1
        i2 = 0
        do iwf2=-iwfmax_small,iwfmax_small-1
           do l=1,ndim
              do k=1,ndim
                 i2 = i2+1
                 
                 m_tot_array(i,j,k,l,iwf2) = m_tot(i1,i2)

              enddo
           enddo
        enddo
     enddo
  enddo

  
  !compute k-dependent self energy:
  do ik=1,nkp
     ikq = kq_ind(ik,iq) 
     do iwf=-iwfmax_small,iwfmax_small-1

        call get_gkiw(ikq, iwf, iwb, iw_data, siw, hk, dc, gkiw)
        
        do i=1,ndim
           do l=1,ndim
              do j=1,ndim
                 do k=1,ndim

                    sigma(i,l,iwf,ik) = sigma(i,l,iwf,ik)+m_tot_array(i,j,k,l,iwf)*gkiw(k,j)
                    !sigma(i,l,iwf,ik) = sigma(i,l,iwf,ik)+m_tot_array(i,j,j,l,iwf)*gkiw(j,j) !test

                 enddo
              enddo
           enddo
        enddo

     enddo
  enddo
             

  deallocate(m_tot, m_tot_array, m_work, u_work)

end subroutine calc_eom

subroutine calc_chi_qw(chi_qw,interm3,chi0_sum)
  use parameters_module
  implicit none
  complex(kind=8) :: chi_qw(ndim2,ndim2)
  complex(kind=8) :: interm3(ndim2,maxdim)
  complex(kind=8) :: chi0_sum(ndim2,ndim2,-iwfmax_small:iwfmax_small-1)

  do i1=1,2*iwfmax_small
    do i2=1,ndim2
      do i3=1,ndim2
        do i4=1,ndim2
          chi_qw(i2,i3)=chi_qw(i2,i3)+interm3(i2,(i1-1)*ndim2+i4)*chi0_sum(i4,i3,i1-1-iwfmax_small)
        end do
      end do
    end do
  end do

end subroutine calc_chi_qw

subroutine get_giw(iw_data, hk, siw, dc, giw)
  use lapack_module
  use parameters_module
  implicit none
 
  double precision :: iw_data(-iwmax:iwmax-1)
  complex(kind=8) :: hk(ndim,ndim,nkp)
  complex(kind=8) :: siw(-iwmax:iwmax-1,ndims)
  double precision :: dc(2,ndim)
  integer :: ik, iw, i
  complex(kind=8) :: g(ndim,ndim), g2(ndim,ndim)
  complex(kind=8), intent (out) :: giw(-iwmax:iwmax-1,ndim)

  giw = 0.d0
  do ik=1,nkp
     g(:,:) = -hk(:,:,ik)
     do iw=0,iwmax-1 !use symmetry of giw(-w)=giw^*(w) 
        do i=1,ndim
           g(i,i) = ci*iw_data(iw)+mu-hk(i,i,ik)-dc(1,i)
        enddo
        do i=1,ndims
           g(i,i) = g(i,i)-siw(iw,i) !no spin dependence in single particle Greens function
        enddo
        g2 = g(:,:)
        call inverse_matrix(g2)
        do i=1,ndim
           giw(iw,i) = giw(iw,i)+g2(i,i)
        enddo
     enddo
  enddo

  do iw=0,iwmax-1
     do i=1,ndim
        giw(-iw-1,i) = real(giw(iw,i),kind=8)-ci*aimag(giw(iw,i))
     enddo
  enddo

  giw = giw/dble(nkp)

end subroutine get_giw



subroutine get_gkiw(ikq, iwf, iwb, iw_data, siw, hk, dc, gkiw)
  use lapack_module
  use parameters_module
  implicit none
  integer :: i
  integer :: iwf, iwb, ikq
  double precision :: iw_data(-iwmax:iwmax-1)
  complex(kind=8) :: hk(ndim,ndim,nkp)
  complex(kind=8) :: siw(-iwmax:iwmax-1,ndims)
  double precision :: dc(2,ndim)
  complex(kind=8), intent(out) :: gkiw(ndim,ndim)

  
  gkiw(:,:) = -hk(:,:,ikq)
  do i=1,ndim
     gkiw(i,i) = ci*iw_data(iwf-iwb)+mu-hk(i,i,ikq)-dc(1,i)
  enddo
  do i=1,ndims
  gkiw(i,i) = gkiw(i,i)-siw(iwf-iwb,i) 
  enddo
  call inverse_matrix(gkiw)
  
  
end subroutine get_gkiw



subroutine get_chi0_loc(beta, iwf, iwb, giw, chi0_loc)
  use parameters_module
  implicit none
  integer :: i, j, k, l
  integer :: iwf, iwb
  double precision :: beta
  complex(kind=8) :: giw(-iwmax:iwmax-1,ndim)
  complex(kind=8), intent(out) :: chi0_loc(ndim*ndim,ndim*ndim)
  
  chi0_loc = 0.d0
  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1=i1+1
        i2=0
        do l=1,ndim
           do k=1,ndim
              i2=i2+1
              if(i==l .and. j==k)then
                 chi0_loc(i1,i2) = - beta*giw(iwf,i)*giw(iwf-iwb,j)
              endif
           enddo
        enddo
     enddo
  enddo

end subroutine get_chi0_loc
  
  
subroutine get_chi0_loc_inv(beta, iwf, iwb, giw, chi0_loc)
  use parameters_module
  implicit none
  integer :: i, j, k, l
  integer :: iwf, iwb
  double precision :: beta
  complex(kind=8) :: giw(-iwmax:iwmax-1,ndim)
  complex(kind=8), intent(out) :: chi0_loc(ndim*ndim,ndim*ndim)
  
  chi0_loc = 0.d0
  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1=i1+1
        i2=0
        do l=1,ndim
           do k=1,ndim
              i2=i2+1
              if(i==l .and. j==k)then
                 chi0_loc(i1,i2) = -1.d0/(beta*giw(iwf,i)*giw(iwf-iwb,j))
              endif
           enddo
        enddo
     enddo
  enddo

end subroutine get_chi0_loc_inv


subroutine get_chi0(beta, ik, ikq, iwf, iwb, iw_data, siw, hk, dc, chi0)
  use lapack_module
  use parameters_module
  implicit none
  integer :: i, j, k, l
  integer :: iwf, iwb, ik, ikq
  complex(kind=8) :: g1(ndim,ndim), g2(ndim,ndim)
  double precision :: iw_data(-iwmax:iwmax-1)
  complex(kind=8) :: hk(ndim,ndim,nkp)
  complex(kind=8) :: siw(-iwmax:iwmax-1,ndims)
  double precision :: dc(2,ndim), beta
  complex(kind=8), intent(out) :: chi0(ndim*ndim,ndim*ndim)

  g1(:,:) = -hk(:,:,ik)
  do i=1,ndim
     g1(i,i) = ci*iw_data(iwf)+mu-hk(i,i,ik)-dc(1,i)
  enddo
  do i=1,ndims
  g1(i,i) = g1(i,i)-siw(iwf,i) 
  enddo
  call inverse_matrix(g1)

  g2(:,:) = -hk(:,:,ikq)
  do i=1,ndim
     g2(i,i) = ci*iw_data(iwf-iwb)+mu-hk(i,i,ikq)-dc(1,i)
  enddo
  do i=1,ndims
  g2(i,i) = g2(i,i)-siw(iwf-iwb,i) 
  enddo
  call inverse_matrix(g2)
  
  chi0 = 0.d0
  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1=i1+1
        i2=0
        do l=1,ndim
           do k=1,ndim
              i2=i2+1
              chi0(i1,i2) = - beta*g1(i,l)*g2(k,j)
           enddo
        enddo
     enddo
  enddo

end subroutine get_chi0



subroutine component2index_band(Nbands, ind, b1, b2, b3, b4)

  implicit none

  integer,intent(in) :: Nbands
  integer,intent(in) :: b1, b2, b3, b4
  integer,intent(out) :: ind

  ind =  Nbands**3*(b1-1) + Nbands**2*(b2-1) + Nbands*(b3-1) + b4

end subroutine component2index_band



! converting an index into a band pattern
subroutine index2component_band(Nbands, ind, b1, b2, b3, b4)

  implicit none

  integer,intent(in) :: Nbands,ind
  integer,intent(out) :: b1, b2, b3, b4
  integer :: tmp1,tmp2,tmp3,ind_tmp
  integer :: g1,g2,g3,g4

  ! the proposed back conversion assumes the indices are
  ! given form 0 to max-1  
  ind_tmp = ind - 1
  tmp1 = Nbands**3
  tmp2 = Nbands**2
  tmp3 = Nbands

  b1 = ind_tmp/tmp1 + 1
  b2 = (ind_tmp-tmp1*(b1-1))/tmp2 + 1
  b3 = (ind_tmp-tmp1*(b1-1)-tmp2*(b2-1))/tmp3 + 1
  b4 = (ind_tmp-tmp1*(b1-1)-tmp2*(b2-1)-tmp3*(b3-1)) + 1


end subroutine index2component_band

