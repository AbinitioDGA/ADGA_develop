program main

#ifdef MPI
  use mpi
#endif

  use hdf5
  use lapack_module
  use parameters_module
  use one_particle_quant_module
  use eom_module
  use susc_module
  use kq_tools

  implicit none

  integer(hid_t) :: file_id, file_vert_id, iw_id, iwb_id, iwf_id, siw_id, giw_id, plist_id, compound_id, r_id, k_id, hk_id, mu_id
  integer :: dc_id, config_id, beta_id
  integer(hid_t) :: iw_space_id, iwb_space_id, iwf_space_id, siw_space_id, giw_space_id, k_space_id, hk_space_id, dc_space_id

  character(len=20) :: grpname_magn, grpname_dens, name_buffer
  character(len=30) :: name_buffer_dset
  character(len=100) :: out_tmp,rank_str
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
  double precision, allocatable :: hk_data(:,:,:,:)
  double precision, allocatable :: dc(:,:)
  complex(kind=8), allocatable :: hk(:,:,:)
  
  integer :: iw, ik, iq, ikq, iwf, iwb, iv, i, j, k, l, n, dum, dum1, ind_iwb, ind_grp, iwf1, iwf2
  integer :: imembers
  complex(kind=8), allocatable :: giw(:,:)
  complex(kind=8), allocatable :: g4iw_magn(:,:,:,:,:,:), g4iw_dens(:,:,:,:,:,:) 
  double precision, allocatable :: tmp_r(:,:), tmp_i(:,:)
  complex(kind=8), allocatable :: chi0_loc(:,:,:), chi0_loc_inv(:,:,:), chi0(:,:), chi0_sum(:,:,:)
  complex(kind=8), allocatable :: chi_loc_slice_dens(:,:), sum_chi0_loc(:,:)
  complex(kind=8), allocatable ::  chi_loc_slice_magn(:,:), chi_loc_dens_full(:,:), chi_loc_magn_full(:,:), chi_loc(:,:)
  complex(kind=8), allocatable ::  chi_loc_magn_sum_left(:,:), chi_loc_dens_sum_left(:,:), gamma_dmft_dens(:,:), gamma_dmft_magn(:,:)
  complex(kind=8), allocatable :: chi_qw_dens(:,:,:),chi_qw_magn(:,:,:),bubble(:,:,:),chi_qw_full(:,:,:)
  complex(kind=8), allocatable :: chi_loc_dens(:,:,:),chi_loc_magn(:,:,:),bubble_loc(:,:,:)
  integer, allocatable :: kq_ind(:,:), qw(:,:)
  complex(kind=8), allocatable :: interm2_magn(:,:), interm3_dens(:,:), interm3_magn(:,:), c(:,:)
  complex(kind=8), allocatable :: interm1(:,:), interm1_v(:,:), interm2_dens(:,:)

  real(kind=8 ):: start, finish, start1, finish1
  complex(kind=8) :: alpha, delta
  integer :: iqw, qwstart, qwstop
  logical :: update_chi_loc_flag
  integer :: b1, b2, b3, b4

  double precision :: u_value, kx, ky, kz
  double precision, allocatable :: u_tmp(:,:,:,:), u_tilde_tmp(:,:,:,:), hr(:,:), hi(:,:)
  complex(kind=8), allocatable :: u(:,:), u_tilde(:,:) 
  complex(kind=8), allocatable :: gamma_loc(:,:), gamma_loc_sum_left(:,:), v(:,:)

  complex(kind=8), allocatable :: sigma(:,:,:,:), sigma_sum(:,:,:,:), sigma_loc(:,:,:)
  integer(hsize_t) ::  inull
  integer :: iwb_zero, iband, ispin

  double precision :: iw_val, giw_r, giw_i, siw_r, siw_i


#ifdef MPI
  integer :: mpi_wrank
  integer :: mpi_wsize
  integer :: master
  integer,allocatable :: rct(:),disp(:)
#endif


  ! read command line argument -> file name of config file
  if (iargc().eq.0 .or. iargc().gt.1) then
    write(*,*) 'The program has to be executed with exactly one argument. (Name of config file)'
    stop
  end if

  call read_config()
  

#ifdef MPI
  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD,mpi_wrank,ierr)
  call MPI_comm_size(MPI_COMM_WORLD,mpi_wsize,ierr)
  master = 0
  allocate(rct(mpi_wsize),disp(mpi_wsize))
#endif


  !THE FOLLOWING PARAMETERS ARE READ FROM THE W2DYNAMICS OUTPUT-FILE:
  !iwmax or iw_dims(1)/2    number of fermionic Matsubara frequencies for single particle quantities 
  !nkp or hk_dims(3)    number of k-points in H(k)
  !ndim or hk_dims(1(2)) total number of bands in H(k) (d+p)
  !ndims or siw_dims(3)   number of d-bands

!#################################################################


! read  external w2wannier Hamitonian:
  open(21, file="/home/lv70808/jkaufmann/Data/Hubbard/1band_beta8/Hk.dat", status='unknown')

  read(21,*) nkp,ndim
  allocate(hr(ndim,ndim),hi(ndim,ndim))
  allocate(hk(ndim,ndim,nkp))
  allocate(k_data(3,nkp))

  do ik=1,nkp

     read(21,*)kx,ky,kz
     k_data(1,ik) = kx
     k_data(2,ik) = ky
     k_data(3,ik) = kz

     do i=1,ndim
        read(21,*) (hr(i,j),hi(i,j),j=1,ndim)
     enddo

     hk(:,:,ik)=hr(:,:)+ci*hi(:,:)
  
  enddo

  close(21)


 
  
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

  if (small_freq_box .eqv. .false.) iwfmax_small = iwfmax
  if (iwfmax_small .gt. iwfmax) then
     write(*,*) 'Error: Maximum number of fermionic frequencies =', iwfmax
  endif
  write(*,*) 'iwfmax=', iwfmax, 'iwfmax_small=', iwfmax_small
  write(*,*) 'iwmax=', iwmax

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

  if (small_freq_box .eqv. .false.) iwbmax_small = iwbmax
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
  deallocate(siw_data)

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
  open(34, file=trim(output_dir)//"siw.dat", status='unknown')
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
  deallocate(giw_data)
  
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
  open(54, file=trim(output_dir)//"giw.dat", status='unknown')
  do iw=-iwmax,iwmax-1
     write(54,'(100F12.6)')iw_data(iw), (real(giw(iw,1)),aimag(giw(iw,1)),i=1,ndims)
  enddo
  close(54)

  write(*,*) 'ok' 
 
! read k-points:
!  call h5dopen_f(file_id, ".axes/k-points", k_id, error)
!  call h5dget_space_f(k_id, k_space_id, error)
!  call h5sget_simple_extent_dims_f(k_space_id, k_dims, k_maxdims, error)
!  nkp = k_dims(2)
!  allocate(k_data(k_dims(1),k_dims(2))) !indices: 3 ik
!  call h5dread_f(k_id, h5t_native_double, k_data, k_dims, error)
!  call h5dclose_f(k_id, error)

write(*,*) nkp,'k points'
! write k-points:
  open(37, file=trim(output_dir)//'k_points.dat', status='unknown')
  do ik=1,100
    write(37,'(100F12.6)') k_data(2,ik), k_data(3,ik)
  enddo
  close(37)

! read Hamiltonian H(k):
!  call h5dopen_f(file_id, "start/hk/value", hk_id, error)
!  call h5dget_space_f(hk_id, hk_space_id, error)
!  call h5sget_simple_extent_dims_f(hk_space_id, hk_dims, hk_maxdims, error)
!  ndim = hk_dims(1)
!  allocate(hk_data(2,hk_dims(1),hk_dims(2),hk_dims(3)))
!  call h5dread_f(hk_id, compound_id, hk_data, hk_dims, error)
!  allocate(hk(hk_dims(1),hk_dims(2),hk_dims(3))) !indices: band band ik
!  hk = 0.d0
!  hk(:,:,:) = hk_data(1,:,:,:)+ci*hk_data(2,:,:,:)
!  call h5dclose_f(hk_id, error)
!  deallocate(hk_data)

! test hk:
!  open(34, file=trim(output_dir)//"hk.dat", status='unknown')
!  do ik=1,hk_dims(3)
!     write(34,*)k_data(:,ik)
!     do i=1,hk_dims(2)
!        write(34,'(100F12.6)')hk(:,i,ik)
!     enddo
!  enddo
!  close(34)

  call init()
!  maxdim = ndim*ndim*2*iwfmax_small
!  ndim2 = ndim*ndim

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

!  allocate(dc(2,ndim))
!  dc=0.d0

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
  write(*,*) 'mu=', mu
  write(*,*) 'dc=', dc

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
!  allocate(giw(-iwmax:iwmax-1),ndim) 
!  call get_giw(iw_data, hk, siw, dc, giw)
 
! test giw:
  open(35, file=trim(output_dir)//"giw_calc.dat", status='unknown')
  do iw=-iwmax,iwmax-1
     write(35,'(100F12.6)') iw_data(iw), (real(giw(iw,i)),aimag(giw(iw,i)),i=1,1)
  enddo
  close(35)


  allocate(chi0_loc(ndim2,ndim2,iwstart:iwstop)) 
  allocate(chi0_sum(ndim2,ndim2,iwstart:iwstop)) 
!  allocate(chi0_sum(ndim2,ndim2,-iwfmax_small:iwfmax_small-1)) 

  allocate(chi0_loc_inv(ndim2,ndim2,-iwfmax:iwfmax-1))
  
  allocate(interm1(ndim2,ndim2))
  allocate(interm1_v(ndim2,ndim2))
  
  allocate(gamma_dmft_dens(ndim2,maxdim), gamma_dmft_magn(ndim2,maxdim))
  allocate(chi_loc_dens_sum_left(ndim2,maxdim))
  allocate(chi_loc_magn_sum_left(ndim2,maxdim))
  allocate(chi_loc_magn_full(maxdim,maxdim))
  allocate(chi_loc_dens_full(maxdim,maxdim))
  allocate(chi_loc_slice_dens(ndim2,maxdim))
  allocate(chi_loc_slice_magn(ndim2,maxdim))
  allocate(interm3_dens(ndim2,maxdim))!allocate somewhere else?
  allocate(interm3_magn(ndim2,maxdim))
  allocate(c(ndim2,maxdim))

  allocate(gamma_loc_sum_left(ndim2,maxdim))
  allocate(v(ndim2,ndim2))

  if (q_path_susc .and. do_chi .and. (.not. q_vol)) then
    if (do_eom) then
      write(*,*) 'Error: it is impossible to use do_eom and q_path_susc'
      stop
    end if
    nqp=n_segments()*nqp1/2+1
    allocate(q_data(nqp))
    call generate_q_path(nqp1,q_data)
    write(*,*) 'q path with',n_segments(),'segments and',nqp,'points.'
  else 
    if (q_path_susc .and. q_vol) then
      write(*,*) 'Warning: q_path_susc .and. q_vol currently has the same effect as only q_vol.'
    end if
    nqp=nqp1**3
    allocate(q_data(nqp))
    call generate_q_vol(nqp1,q_data)
  end if

  if (mpi_wrank .eq. 1) then
    write(*,*) nkp1,'k points in one direction'
    write(*,*) nkp,'k points total'
    write(*,*) nqp1,'k points in one direction'
    write(*,*) nqp,'q points total'
  end if

!  else if (.not. do_eom .and. do_chi .and. q_path) then
!    nqp=n_segments()*nqp1/2+1
!    allocate(q_data(nqp))
!    call generate_q_path(q_data)
!    write(*,*) 'q path'
!    write(*,*)'nqp=', nqp !test
!    write(*,*) q_data
!  else if (do_eom .and. q_path) then
!    nqp = nqp1**3
!    allocate(q_data(nqp))
!    call generate_q_vol(nqp1,q_data)
!    nkp_eom=n_segments()*nqp1/2+1
!    allocate(k_data_eom(nkp_eom))
!    call generate_q_path(k_data_eom)
!  else
!    write(*,*) 'only (q_vol) and (.not. do_eom .and. q_path) are implemented yet.'
!    stop
!  end if
  
  !search for k+q - index:
  call cpu_time(start)
  allocate(kq_ind(nkp,nqp))
!  call index_kq_search(k_data, q_data, kq_ind) ! old method, assumes cubic case
  call index_kq(kq_ind) ! new method
  call cpu_time(finish)
  write(*,*)'finding k-q index:', finish-start
!##################### parallel code ##################################################

  write(*,*)'nqp=', nqp !test

  !stop
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

!distribute the qw compound index:
#ifdef MPI
  rct=0
  disp=0
  qwstop = 0
  do i=0,mpi_wrank
     j = (nqp*(2*iwbmax_small+1) - qwstop)/(mpi_wsize-i)
     qwstart = qwstop + 1
     qwstop = qwstop + j
  enddo
  rct(mpi_wrank+1)=(qwstop-qwstart+1)*ndim2**2
  if (mpi_wrank .eq. master) then
    call MPI_reduce(MPI_IN_PLACE,rct,mpi_wsize,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  else
    call MPI_reduce(rct,rct,mpi_wsize,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  end if
  if (mpi_wrank .eq. master) then
    do i=2,mpi_wsize
      disp(i)=sum(rct(1:i-1))! the first displacing has to be 0
    end do
    write(*,*) 'receive ct',rct
    write(*,*) 'displacing',disp
  end if
#else 
  qwstart = 1
  qwstop = nqp*(2*iwbmax_small+1)
#endif

if (do_chi) then
  allocate(chi_qw_dens(ndim2,ndim2,qwstart:qwstop),chi_qw_magn(ndim2,ndim2,qwstart:qwstop))
  allocate(bubble(ndim2,ndim2,qwstart:qwstop))
  allocate(chi_loc_dens(ndim2,ndim2,-iwbmax_small:iwbmax_small),chi_loc_magn(ndim2,ndim2,-iwbmax_small:iwbmax_small))
  allocate(bubble_loc(ndim2,ndim2,-iwbmax_small:iwbmax_small))
  chi_qw_dens=0.d0
  chi_qw_magn=0.d0
  bubble=0.d0
  chi_loc_dens=0.d0
  chi_loc_magn=0.d0
  bubble_loc=0.d0
end if

#ifdef MPI
write(*,*)'rank=',mpi_wrank, 'qwstart=', qwstart, 'qwstop=', qwstop
start = mpi_wtime()
#endif


if (do_eom) then
  allocate(sigma(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp))
  sigma = 0.d0
end if

!  open(55, file=trim(output_dir)//"chi0_loc_sum.dat", status='unknown')

  iwb = iwbmax_small+3
  do iqw=qwstart,qwstop
     call cpu_time(start)
     update_chi_loc_flag = qw(1,iqw) .ne. iwb
     iq = qw(2,iqw)
     iwb = qw(1,iqw)
  
     !to be done here: read nonlocal interaction v and go into compound index
     v = 0.d0

     !update chi_loc only if iwb is different than the previous one:
     if(update_chi_loc_flag) then    
        
        chi0_loc=0.d0
        ! compute local bubble chi0_loc^{-1}(i1,i2)(orbital compound index i1,i2):
        do iwf=-iwfmax,iwfmax-1
           call get_chi0_loc_inv(iwf, iwb, giw, chi0_loc_inv(:,:,iwf))
        enddo
        if (do_chi) then
           do iwf=iwstart,iwstop
             call get_chi0_loc(  iwf, iwb, giw, chi0_loc(:,:,iwf))
           enddo
        end if
       
!        if (do_chi .and. update_chi_loc_flag) then

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
        chi_loc_magn_full = 0.d0
        chi_loc_dens_full = 0.d0

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
                          chi_loc_magn_full(i1,i2) = g4iw_magn(i,j,iwf1,k,l,iwf2)
                          chi_loc_dens_full(i1,i2) = g4iw_dens(i,j,iwf1,k,l,iwf2)

                          !straight term is subtracted (twice) only in the dens channel and only for iw=0:
                          if((iwb .eq. iwb_zero) .and. i==j .and. k==l)then
                             chi_loc_dens_full(i1,i2) = chi_loc_dens_full(i1,i2)-2.d0*beta*giw(iwf1,i)*giw(iwf2,l) 
                          endif

                          !HACK: for reading vertex_sym.hdf5 from g4iw_conn.hdf5 instead of
                          !vertex_full.hdf5
                          !cross term is subtracted once for each channel
                          !if((iwf2 .eq. iwf1) .and. i==l .and. j==k)then
                          !   chi_loc_dens_full(i1,i2) = chi_loc_dens_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j) 
                          !   chi_loc_magn_full(i1,i2) = chi_loc_magn_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j) 
                          !endif

                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo

        deallocate(g4iw_magn, g4iw_dens, tmp_r, tmp_i)
       
        !time reversal symmetry (which is simply a transpose in our compound index)
        do i1=1,maxdim
           do i2=i1+1,maxdim
              chi_loc_magn_full(i1,i2) = 0.5d0*(chi_loc_magn_full(i1,i2)+chi_loc_magn_full(i2,i1))
              chi_loc_magn_full(i2,i1) = chi_loc_magn_full(i1,i2)

              chi_loc_dens_full(i1,i2) = 0.5d0*(chi_loc_dens_full(i1,i2)+chi_loc_dens_full(i2,i1))
              chi_loc_dens_full(i2,i1) = chi_loc_dens_full(i1,i2)
           enddo
        enddo

        !compute chi_loc*chi0_loc_inv (chi0_loc_inv is diagonal in compound index):
        !use blas-routine instead?
        dum = 0
        do iwf=-iwfmax_small,iwfmax_small-1
           do i2=1,ndim2
              do i1=1,maxdim
                 chi_loc_magn_full(i1,i2+dum*ndim2) = chi_loc_magn_full(i1,i2+dum*ndim2)*chi0_loc_inv(i2,i2,iwf)
                 chi_loc_dens_full(i1,i2+dum*ndim2) = chi_loc_dens_full(i1,i2+dum*ndim2)*chi0_loc_inv(i2,i2,iwf)
              enddo
           enddo
           dum = dum+1
        enddo

        !sum up the left fermionic frequency of chi_loc (quantity needed afterwards)
        chi_loc_magn_sum_left = 0.d0!von rechts mit chi0_loc mutliplizieren
        chi_loc_dens_sum_left = 0.d0

        do i1=1,ndim2
           do dum=0,2*iwfmax_small-1
              chi_loc_magn_sum_left(i1,:) = chi_loc_magn_sum_left(i1,:)+chi_loc_magn_full(i1+dum*ndim2,:)
              chi_loc_dens_sum_left(i1,:) = chi_loc_dens_sum_left(i1,:)+chi_loc_dens_full(i1+dum*ndim2,:)
           enddo
        enddo
        
       
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
      
     endif !update local quantities

     dum1 = 0 !index for second matrix multiplication to get interm2
     
     allocate(interm2_dens(maxdim,maxdim))
     allocate(interm2_magn(maxdim,maxdim))
     allocate(gamma_loc(maxdim,maxdim))
     allocate(chi0(ndim2,ndim2))
     interm2_dens = 0.d0
     interm2_magn = 0.d0
     gamma_loc = 0.d0
     chi0_sum=0.d0
     
     do iwf=iwstart,iwstop
        ! compute k-summed (but still q-dependent) bubble chi0(i1,i2):
        do ik=1,nkp
           ikq = kq_ind(ik,iq) !Index of G(k+q)
           call get_chi0(ik, ikq, iwf, iwb, iw_data, siw, hk, dc, chi0(:,:)) 
           chi0_sum(:,:,iwf) = chi0_sum(:,:,iwf)+chi0(:,:)
        enddo
        chi0_sum(:,:,iwf) = chi0_sum(:,:,iwf)/dble(nkp)
     end do
!     call cpu_time(finish)
!     write(*,*) finish-start

     do iwf=-iwfmax_small,iwfmax_small-1
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
        chi_loc_slice_magn(:,:) = chi_loc_magn_full((dum1*ndim2)+1 : (dum1+1)*ndim2, :)
        chi_loc_slice_dens(:,:) = chi_loc_dens_full((dum1*ndim2)+1 : (dum1+1)*ndim2, :)
        

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
        do i=1,ndim2
           do j=1,ndim2
              interm2_dens(i+(dum1*ndim2),j+(dum1*ndim2)) = interm2_dens(i+(dum1*ndim2),j+(dum1*ndim2))-interm1(i,j)
              interm2_magn(i+(dum1*ndim2),j+(dum1*ndim2)) = interm2_magn(i+(dum1*ndim2),j+(dum1*ndim2))-interm1(i,j)
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
     
     deallocate(chi0)

     !sum up left iwf-index of interm2_dens and store the quantity (is afterwards directly used in the EOM):
     gamma_loc_sum_left = 0.d0

     do i1=1,ndim2
        do dum=0,2*iwfmax_small-1
           gamma_loc_sum_left(i1,:) = gamma_loc_sum_left(i1,:)+gamma_loc(i1+dum*ndim2,:)
        enddo
     enddo

     deallocate(gamma_loc)


     !construct quantity that is inverted right afterwards:
     interm2_dens = chi_loc_dens_full-interm2_dens
     interm2_magn = chi_loc_magn_full-interm2_magn
    
 
     call inverse_matrix(interm2_dens)
     call inverse_matrix(interm2_magn)


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

     if (do_chi) then
        ! Calculation of q dependent susceptibility by multiplication with chi0
        call calc_chi_qw(chi_qw_dens(:,:,iqw),interm3_dens,chi0_sum(:,:,-iwfmax_small:iwfmax_small-1))
        call calc_chi_qw(chi_qw_magn(:,:,iqw),interm3_magn,chi0_sum(:,:,-iwfmax_small:iwfmax_small-1))
        call calc_bubble(bubble(:,:,iqw),chi0_sum)

        ! Calculation of local susceptibility and bubble
        if (iq.eq.1) then !this should be calculated only once, otherwise wrong result due to mpi_reduce sum.
          call calc_chi_qw(chi_loc_dens(:,:,iwb),chi_loc_dens_sum_left,chi0_loc(:,:,-iwfmax_small:iwfmax_small-1))
          call calc_chi_qw(chi_loc_magn(:,:,iwb),chi_loc_magn_sum_left,chi0_loc(:,:,-iwfmax_small:iwfmax_small-1))
          call calc_bubble(bubble_loc(:,:,iwb),chi0_loc)
        end if
     end if
     if (do_eom) then
        !equation of motion     
        call calc_eom(interm3_dens,interm3_magn,gamma_dmft_dens,gamma_dmft_magn,gamma_loc_sum_left,sigma,kq_ind,iwb,iq,iw_data,u,v,u_tilde,hk,dc,siw)
     end if
     call cpu_time(finish)

     !Output the calculation progress
!     if (mpi_wrank .eq. master) then
       write(*,'((A),I5,2X,I6,2X,I6,2X,I6,2X,(A),F8.4)') 'iqw/qwstart/qwstop on rank',mpi_wrank,iqw,qwstart,qwstop,'time ',finish-start
!     end if
  enddo !iqw

  deallocate(interm3_dens,interm3_magn) 
  deallocate(hk,giw,dc,u,u_tilde,chi0_loc,chi0_loc_inv,chi0_sum,interm1,interm1_v,gamma_dmft_dens,gamma_dmft_magn)
  deallocate(chi_loc_dens_sum_left,chi_loc_magn_sum_left)
  deallocate(chi_loc_magn_full,chi_loc_dens_full,gamma_loc_sum_left,v,c)
  deallocate(chi_loc_slice_dens,chi_loc_slice_magn)
  
#ifdef MPI
  ! MPI reduction and output
  if (do_eom) then
     allocate(sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp))
     allocate(sigma_loc(ndim, ndim, -iwfmax_small:iwfmax_small-1))
     call MPI_reduce(sigma, sigma_sum, ndim*ndim*2*iwfmax_small*nkp, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
     sigma_sum = -sigma_sum/(beta*nqp)
     call add_siw_dmft(siw, sigma_sum, sigma_loc)

     if (mpi_wrank .eq. master) then
       call output_eom(iw_data, k_data, sigma_sum, sigma_loc)
     end if
     deallocate(sigma,sigma_sum,sigma_loc)
  end if

  if (do_chi) then
    if (mpi_wrank.eq.master) then
       call MPI_reduce(MPI_IN_PLACE,chi_loc_dens,ndim2*ndim2*(2*iwbmax_small+1),MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
       call MPI_reduce(MPI_IN_PLACE,chi_loc_magn,ndim2*ndim2*(2*iwbmax_small+1),MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
       call MPI_reduce(MPI_IN_PLACE,bubble_loc,ndim2*ndim2*(2*iwbmax_small+1),MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    else 
       call MPI_reduce(chi_loc_dens,chi_loc_dens,ndim2*ndim2*(2*iwbmax_small+1),MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
       call MPI_reduce(chi_loc_magn,chi_loc_magn,ndim2*ndim2*(2*iwbmax_small+1),MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
       call MPI_reduce(bubble_loc,bubble_loc,ndim2*ndim2*(2*iwbmax_small+1),MPI_DOUBLE_COMPLEX,MPI_SUM,master,MPI_COMM_WORLD,ierr)
       deallocate(chi_loc_dens,chi_loc_magn,bubble_loc)
    end if
    allocate(chi_qw_full(1,1,1))
    if (mpi_wrank .eq. master) then 
      deallocate(chi_qw_full)
      allocate(chi_qw_full(ndim2,ndim2,nqp*(2*iwbmax+1))) 
    end if
    chi_qw_full=0.d0

    call MPI_gatherv(chi_qw_dens,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,chi_qw_full,rct,disp,MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
    deallocate(chi_qw_dens)
    if (mpi_wrank .eq. master) then 
      call output_chi_qw(chi_qw_full,iwb_data,qw,'chi_qw_dens.dat')
    end if

    call MPI_gatherv(chi_qw_magn,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,chi_qw_full,rct,disp,MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
    deallocate(chi_qw_magn)
    if (mpi_wrank .eq. master) then 
      call output_chi_qw(chi_qw_full,iwb_data,qw,'chi_qw_magn.dat')
    end if

    call MPI_gatherv(bubble,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,chi_qw_full,rct,disp,     MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
    deallocate(bubble)
    if (mpi_wrank .eq. master) then 
      call output_chi_qw(chi_qw_full,iwb_data,qw,'bubble.dat')
    end if

    deallocate(chi_qw_full)
  end if
  
  call MPI_finalize(ierr)

! Output
  if (mpi_wrank .eq. master) then
     if (do_chi) then
        call output_chi_loc(chi_loc_dens,iwb_data,'chi_loc_dens.dat')
        call output_chi_loc(chi_loc_magn,iwb_data,'chi_loc_magn.dat')
        call output_chi_loc(bubble_loc,iwb_data,'bubble_loc.dat')
        deallocate(chi_loc_dens,chi_loc_magn,bubble_loc)
     end if
  endif
  
#endif
deallocate(iw_data,iwb_data,siw,k_data,q_data,kq_ind,qw)
end program main




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

