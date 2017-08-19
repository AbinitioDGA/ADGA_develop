program main

#ifdef MPI
  use mpi
#endif

  use hdf5_module
  use lapack_module
  use parameters_module
  use one_particle_quant_module
  use eom_module
  use susc_module
  use kq_tools
  use vq_module

  implicit none

  integer(hid_t) :: file_id, file_vert_id, iw_id, iwb_id, iwf_id, r_id, k_id, hk_id
  integer :: config_id, beta_id
  integer(hid_t) :: iw_space_id, iwb_space_id, iwf_space_id, k_space_id, hk_space_id

  character(len=20) :: grpname_magn, grpname_dens, name_buffer
  character(len=100) :: name_buffer_dset
  character(len=100) :: out_tmp,rank_str
  integer(hid_t) :: grp_magn_id, grp_dens_id, nmembers, itype, dset_magn_id, dset_dens_id
  integer(hsize_t), dimension(2) :: tmp_dims

  integer :: error, ierr
  integer(hsize_t), dimension(0) :: beta_dims


  integer :: iw, ik, iq, ikq, iwf, iwb, iv, dum, dum1, ind_iwb, ind_grp, iwf1, iwf2
  integer :: i, j, k, l, n, i1, i2, i3, i4
  integer :: ineq,dimstart, dimend
  integer :: imembers
  complex(kind=8), allocatable :: gloc(:,:,:), gkiw(:,:)
  complex(kind=8), allocatable :: g4iw_magn(:,:,:,:,:,:), g4iw_dens(:,:,:,:,:,:)
  double precision, allocatable :: tmp_r(:,:), tmp_i(:,:)
  complex(kind=8), allocatable :: chi0_loc(:,:,:), chi0_loc_inv(:,:,:), chi0(:,:), chi0_sum(:,:,:)
  complex(kind=8), allocatable :: chi_loc_slice_dens(:,:), sum_chi0_loc(:,:)
  complex(kind=8), allocatable ::  chi_loc_slice_magn(:,:), chi_loc_dens_full(:,:), chi_loc_magn_full(:,:), chi_loc(:,:)
  complex(kind=8), allocatable ::  chi_loc_magn_sum_left(:,:), chi_loc_dens_sum_left(:,:), gamma_dmft_dens(:,:), gamma_dmft_magn(:,:)
  complex(kind=8), allocatable :: chi_qw_dens(:,:,:),chi_qw_magn(:,:,:),bubble(:,:,:),chi_qw_full(:,:,:)
  complex(kind=8), allocatable :: chi_loc_dens(:,:,:),chi_loc_magn(:,:,:),bubble_loc(:,:,:)
  integer, allocatable :: kq_ind(:,:), qw(:,:)
  complex(kind=8), allocatable :: interm2_magn(:,:), interm3_dens(:,:), interm3_magn(:,:), c_interm(:,:)
  complex(kind=8), allocatable :: interm1(:,:), interm1_v(:,:), interm2_dens(:,:)

  real(kind=8 ):: start, finish, start1, finish1
  complex(kind=8) :: alpha, delta
  integer :: iqw, qwstart, qwstop
  logical :: update_chi_loc_flag
  integer :: b1, b2, b3, b4

  complex(kind=8), allocatable :: gamma_loc(:,:), gamma_loc_sum_left(:,:), v(:,:)

  complex(kind=8), allocatable :: sigma(:,:,:,:), sigma_hf(:,:,:,:), sigma_dmft(:,:,:)
  complex(kind=8), allocatable :: sigma_sum(:,:,:,:), sigma_sum_hf(:,:,:,:), sigma_sum_dmft(:,:,:), sigma_loc(:,:,:)
  complex(kind=8), allocatable :: giw_sum(:), n_dga(:), n_dmft(:), n_fock(:,:,:)
  integer :: iband, ispin

  double precision :: iw_val, giw_r, giw_i


#ifdef MPI
  integer :: mpi_wrank
  integer :: mpi_wsize
  integer :: master
  integer,allocatable :: rct(:),disp(:)
#endif

! variables for date-time string
  character(20) :: date,time,zone
  character(200) :: output_filename
  integer,dimension(8) :: values


  ! read command line argument -> file name of config file
  if (iargc() .ne. 1) then
    write(*,*) 'The program has to be executed with exactly one argument. (Name of config file)'
    stop
  end if


#ifdef MPI
  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD,mpi_wrank,ierr)
  call MPI_comm_size(MPI_COMM_WORLD,mpi_wsize,ierr)
  master = 0
  allocate(rct(mpi_wsize),disp(mpi_wsize))
#endif

  call read_config()

! master id creates output folder if necessary
if (mpi_wrank .eq. master) call system("mkdir -p "//trim(output_dir))

if (mpi_wrank .eq. master) then
  write(*,*)
  write(*,*) '/-----------------------------------------------------------\'
  write(*,*) '|  Ab initio dynamical vertex approximation program (ADGA)  |'
  write(*,*) '|  Running on ',mpi_wsize,' core(s).                        |'
  write(*,*) '\-----------------------------------------------------------/'
  write(*,*)
end if

  !THE FOLLOWING PARAMETERS ARE READ FROM THE W2DYNAMICS OUTPUT-FILE:
  !iwmax or iw_dims(1)/2    number of fermionic Matsubara frequencies for single particle quantities
  !nkp or hk_dims(3)    number of k-points in H(k)
  !ndim or hk_dims(1(2)) total number of bands in H(k) (d+p)
  !ndims or siw_dims(3)   number of d-bands

!#################################################################


! read  external w2wannier Hamitonian:
  if(read_ext_hk) then
    call read_hk_w2w()
  end if


!##################  READ W2DYNAMICS HDF5 OUTPUT FILE  #####################################

  call init_h5() ! open the hdf5-fortran interface

! read bosonic and fermionic Matsubara axes iwf-g4,iwb-g4:
  call get_freq_range()


  if (iwfmax_small .le. 0 .or. iwfmax_small .gt. iwfmax) then
    iwfmax_small = iwfmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Error: Wrong input for fermionic frequencies'
      write(*,*) 'Calculating with maximum number of fermionic frequencies &
                  =', iwfmax
    endif
  endif

  if (iwbmax_small .le. 0 .or. iwbmax_small .gt. iwbmax) then
    iwbmax_small = iwbmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Error: Wrong input for bosonic frequencies'
      write(*,*) 'Calculating with maximum number of bosonic frequencies &
                  =', iwbmax
    endif
  endif

  if (mpi_wrank .eq. master) then
    write(*,*) 'iwmax=', iwmax, ' (number of fermionic matsubara frequencies of one-particle quantities)'
    write(*,*) 'iwfmax=', iwfmax, 'iwfmax_small=', iwfmax_small, ' (number of fermionic matsubara frequencies of two-particle quantities)'
    write(*,*)'iwbmax=',iwbmax, 'iwbmax_small=', iwbmax_small, ' (number of bosonic matsubara frequencies of two-particle quantities)'
  end if

! read in all inequivalent atoms
  ! siw from DMFT contains all possible bands
  ! siw at p (non-interacting bands) is set to 0
  call read_siw()

  call read_giw()

  if (mpi_wrank .eq. master) then
    write(*,*) 'orb_symmetry = ', orb_sym
  end if


  if(.not. read_ext_hk) then
    call read_hk_w2dyn()
  end if


  call init()

  if (.not. do_vq .and. mpi_wrank .eq. master) then
    write(*,*) 'Main: Run without V(q)'
  end if

  if (mpi_wrank .eq. master) then
    ! generate a date-time string for output file name
    call date_and_time(date,time,zone,values)
    output_filename=trim(output_dir)//'adga-'//trim(date)//'-'//trim(time)//'-output.hdf5'
    write(*,*) 'writing output to ',output_filename
    call init_h5_output(output_filename)
  end if

  call read_mu()
  call read_beta()
  call read_dc()


! read double counting:


  call finalize_h5() ! close the hdf5-fortran interface

  if (mpi_wrank .eq. master) then
    write(*,*) 'beta=', beta
    write(*,*) 'mu=', mu
    write(*,*) 'dc=', dc
  end if


  !read umatrix from separate file:
  call read_u(u,u_tilde)


!################################################################################################

! compute local single-particle Greens function:
! allocate(giw(-iwmax:iwmax-1,ndim))
  call get_giw()

! test giw:
!   open(35, file=trim(output_dir)//"giw_calc.dat", status='unknown')
!   do iw=-iwmax,iwmax-1
!      write(35,'(100F12.6)') iw_data(iw), (real(giw(iw,i)),aimag(giw(iw,i)),i=1,ndim)
!   enddo
!   close(35)

  !compute DMFT filling n_dmft
  allocate(gloc(-iwmax:iwmax-1,ndim,ndim), gkiw(ndim,ndim))
  allocate(giw_sum(ndim), n_dmft(ndim), n_fock(nkp,ndim,ndim), n_dga(ndim))
  giw_sum = 0.d0
  n_dmft = 0.d0
  do iw=0,iwmax-1
      giw_sum(:) = giw_sum(:)+giw(iw,:)
  enddo

  n_dmft = 0.d0
  n_dmft(:) = 2.d0*real(giw_sum(:))/beta+0.5d0
  open(56, file=trim(output_dir)//"n_dmft.dat", status='unknown')
  write(56,'(100F12.6)') (real(n_dmft(i)),i=1,ndim)

  !compute k-dependent filling for Fock-term (computed in the EOM):
  n_fock = 0.d0
  gkiw = 0.d0
  open(110, file=trim(output_dir)//"n_fock.dat", status='unknown')
  do ik=1,nkp
     do iw=0,iwmax-1
        call get_gkiw(ik, iw, 0,gkiw)
        n_fock(ik,:,:) = n_fock(ik,:,:)+real(gkiw(:,:))
     enddo
     n_fock = 2.d0*n_fock/beta
     do i=1,ndim
        n_fock(ik,i,i) = n_fock(ik,i,i)+0.5d0
     enddo
     write(110,'(100F12.6)') k_data(1,ik),k_data(2,ik),k_data(3,ik), (real(n_fock(ik,i,i)),i=1,ndim)
  enddo


  allocate(chi0_loc(ndim2,ndim2,iwstart:iwstop))
  allocate(chi0_sum(ndim2,ndim2,iwstart:iwstop))

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
  allocate(c_interm(ndim2,maxdim))

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
    nqp=nqpx*nqpy*nqpz
    allocate(q_data(nqp))
    call generate_q_vol(nqpx,nqpy,nqpz,q_data)
  end if


  if (mpi_wrank .eq. master) then
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
  if (mpi_wrank .eq. master) then
    write(*,*)'finding k-q index:', finish-start
  end if
!##################### parallel code ##################################################

  if (mpi_wrank .eq. master) then
    write(*,*)'nqp=', nqp !test
    write(*,*) maxval(kq_ind)
  end if


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
  allocate(sigma(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp), sigma_hf(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp))
  allocate(sigma_dmft(ndim,ndim,-iwfmax_small:iwfmax_small-1))
  sigma = 0.d0
  sigma_hf = 0.d0
  sigma_dmft = 0.d0
end if


if (mpi_wrank.eq.0) then
  open(unit=256,file=trim(output_dir)//'kdata',status='replace')
  do ik=1,nkp
    write(256,*) k_data(:,ik)
  end do
  close(256)
  open(unit=266,file=trim(output_dir)//'qdata',status='replace')
  do iq=1,nqp
    write(266,*) k_data(:,q_data(iq))
  end do
  close(266)
end if


  iwb = iwbmax_small+3
  do iqw=qwstart,qwstop
     call cpu_time(start)
     update_chi_loc_flag = qw(1,iqw) .ne. iwb
     iq = qw(2,iqw)
     iwb = qw(1,iqw)

     !read nonlocal interaction v and go into compound index:
     if(do_vq) then
        call read_vq(iq,v)
       ! v = v-u  !otherwise, local U would be included twice
     else
        v = 0.d0
     endif

     !update chi_loc only if iwb is different than the previous one:
     if(update_chi_loc_flag) then

        chi0_loc=0.d0
        ! compute local bubble chi0_loc^{-1}(i1,i2)(orbital compound index i1,i2):
        do iwf=-iwfmax,iwfmax-1
           call get_chi0_loc_inv(iwf, iwb, chi0_loc_inv(:,:,iwf))
        enddo
        if (do_chi) then
           do iwf=iwstart,iwstop
             call get_chi0_loc(  iwf, iwb, chi0_loc(:,:,iwf))
           enddo
        end if

!        if (do_chi .and. update_chi_loc_flag) then


        call read_vertex(chi_loc_dens_full,chi_loc_magn_full,iwb)

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
           call get_chi0(ik, ikq, iwf, iwb, chi0(:,:))
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
        c_interm = 0.d0
        alpha = 1.d0
        delta = 0.d0
        call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, interm1, ndim2, chi_loc_slice_dens, ndim2, delta, c_interm, ndim2)

        do i=1,ndim2
           do j=1,maxdim
              interm2_dens(i+(dum1*ndim2),j) = c_interm(i,j)
           enddo
        enddo

        !same for the magn channel:
        c_interm = 0.d0
        alpha = 1.d0
        delta = 0.d0
        call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, interm1, ndim2, chi_loc_slice_magn, ndim2, delta, c_interm, ndim2)

        do i=1,ndim2
           do j=1,maxdim
              interm2_magn(i+(dum1*ndim2),j) = c_interm(i,j)
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
        ! only in density channel
        c_interm = 0.d0
        alpha = 1.d0
        delta = 0.d0
        call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, interm1_v, ndim2, chi_loc_dens_sum_left, ndim2, delta, c_interm, ndim2)

        ! add it to interm2:
        do i=1,ndim2
           do j=1,maxdim
              interm2_dens(i+(dum1*ndim2),j) = interm2_dens(i+(dum1*ndim2),j)+c_interm(i,j)
           enddo
        enddo

!        !same for the magnetic channel:
!        c_interm = 0.d0
!        alpha = 1.d0
!        delta = 0.d0
!        call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, interm1_v, ndim2, chi_loc_magn_sum_left, ndim2, delta, c_interm, ndim2)
!
!        do i=1,ndim2
!           do j=1,maxdim
!              interm2_magn(i+(dum1*ndim2),j) = interm2_magn(i+(dum1*ndim2),j)+c_interm(i,j)
!           enddo
!        enddo


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
        call calc_eom(interm3_dens,interm3_magn,gamma_dmft_dens,gamma_dmft_magn,gamma_loc_sum_left,sigma,sigma_dmft,sigma_hf,kq_ind,iwb,iq,iw_data,u,v,u_tilde,n_dmft,n_fock)
     end if
     call cpu_time(finish)

     !Output the calculation progress
!     if (mpi_wrank .eq. master) then
       write(*,'((A),I5,2X,I6,2X,I6,2X,I6,2X,(A),F8.4)') 'iqw/qwstart/qwstop on rank',mpi_wrank,iqw,qwstart,qwstop,'time ',finish-start
!     end if
  enddo !iqw


!  close(267)

  deallocate(interm3_dens,interm3_magn)
  deallocate(giw,u,u_tilde,chi0_loc,chi0_loc_inv,chi0_sum,interm1,interm1_v,gamma_dmft_dens,gamma_dmft_magn)
  deallocate(chi_loc_dens_sum_left,chi_loc_magn_sum_left)
  deallocate(chi_loc_magn_full,chi_loc_dens_full,gamma_loc_sum_left,v,c_interm)
  deallocate(chi_loc_slice_dens,chi_loc_slice_magn)

#ifdef MPI
  ! MPI reduction and output
  write(*,*)'nkp=', nkp
  if (do_eom) then
     allocate(sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp), sigma_sum_dmft(ndim, ndim, -iwfmax_small:iwfmax_small-1))
     allocate(sigma_sum_hf(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp))
     allocate(sigma_loc(ndim, ndim, -iwfmax_small:iwfmax_small-1))
     call MPI_reduce(sigma, sigma_sum, ndim*ndim*2*iwfmax_small*nkp, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
     call MPI_reduce(sigma_hf, sigma_sum_hf, ndim*ndim*2*iwfmax_small*nkp, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
     call MPI_reduce(sigma_dmft, sigma_sum_dmft, ndim*ndim*2*iwfmax_small, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
     sigma_sum = -sigma_sum/(beta*nqp)
     sigma_sum_hf = -sigma_sum_hf/(beta*nqp)
     sigma_sum_dmft = -sigma_sum_dmft/(beta*nqp)
     call add_siw_dmft(sigma_sum)
     call get_sigma_g_loc(iw_data, sigma_sum, sigma_loc, gloc, n_dga)
     if (mpi_wrank .eq. master) then
       call output_eom(iw_data, k_data, sigma_sum, sigma_sum_dmft, sigma_sum_hf, sigma_loc, gloc, n_dga)
     end if
     deallocate(sigma, sigma_sum, sigma_sum_dmft, sigma_sum_hf, sigma_loc)
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
      call output_chi_qw(chi_qw_full,qw,'chi_qw_dens.dat')
      call output_chi_qw_h5(output_filename,'dens',chi_qw_full)
    end if

    call MPI_gatherv(chi_qw_magn,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,chi_qw_full,rct,disp,MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
    deallocate(chi_qw_magn)
    if (mpi_wrank .eq. master) then
      call output_chi_qw(chi_qw_full,qw,'chi_qw_magn.dat')
      call output_chi_qw_h5(output_filename,'magn',chi_qw_full)
    end if

    call MPI_gatherv(bubble,(qwstop-qwstart+1)*ndim2**2,MPI_DOUBLE_COMPLEX,chi_qw_full,rct,disp,     MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
    deallocate(bubble)
    if (mpi_wrank .eq. master) then
      call output_chi_qw(chi_qw_full,qw,'bubble.dat')
      call output_chi_qw_h5(output_filename,'bubble',chi_qw_full)
    end if

    deallocate(chi_qw_full)
  end if


! Output
  if (mpi_wrank .eq. master) then
     if (do_chi) then
        call output_chi_loc(chi_loc_dens,'chi_loc_dens.dat')
        call output_chi_loc(chi_loc_magn,'chi_loc_magn.dat')
        call output_chi_loc(bubble_loc,'bubble_loc.dat')
        call output_chi_loc_h5(output_filename,'dens',chi_loc_dens)
        call output_chi_loc_h5(output_filename,'magn',chi_loc_magn)
        call output_chi_loc_h5(output_filename,'bubble',bubble_loc)
        deallocate(chi_loc_dens,chi_loc_magn,bubble_loc)
     end if
  endif

  call MPI_finalize(ierr)

#endif
deallocate(iw_data,iwb_data,siw,k_data,q_data,kq_ind,qw)
write(*,*) 'end of program'
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

