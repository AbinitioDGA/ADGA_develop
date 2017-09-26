program main

#ifdef MPI
  use mpi
#endif

  use hdf5_module
  use lapack_module
  use parameters_module
  use config_module
  use one_particle_quant_module
  use eom_module
  use susc_module
  use kq_tools
  use vq_module
  use aux
  use mpi_org

  implicit none
  integer :: iw, ik, iq, ikq, iwf, iwb, iv, dum, dum1, iwf1, iwf2
  integer :: i, j, k, l, n, i1, i2, i3, i4
  complex(kind=8), allocatable :: gloc(:,:,:)
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
  integer :: iqw
  logical :: update_chi_loc_flag
  complex(kind=8), allocatable :: gamma_loc(:,:), gamma_loc_sum_left(:,:), v(:,:)
  complex(kind=8), allocatable :: sigma(:,:,:,:), sigma_hf(:,:,:,:), sigma_dmft(:,:,:)
  complex(kind=8), allocatable :: sigma_sum(:,:,:,:), sigma_sum_hf(:,:,:,:), sigma_sum_dmft(:,:,:), sigma_loc(:,:,:)
! variables for date-time string
  character(20) :: date,time,zone
  character(200) :: output_filename
  integer,dimension(8) :: time_date_values


  ! read command line argument -> file name of config file
  if (iargc() .ne. 1) then
    stop 'The program has to be executed with exactly one argument. (Name of config file)'
  end if

  ! mpi initialization
  call mpi_initialize()

  ! read config settings
  if (mpi_wrank .eq. master) write(*,*) 'Reading config file'
  call read_config()

  write(*,*) nkpx, nkpy, nkpz

  ! check config settings
  if (mpi_wrank .eq. master) write(*,*) 'Checking config file'
  call check_config() 

  ! create output folder if not yet existing
  if (mpi_wrank .eq. master) call system('mkdir -p ' // trim(adjustl(output_dir)))
  call mpi_barrier(mpi_comm_world,ierr)

  ! program introduction
  if (mpi_wrank .eq. master) then
    write(*,*)
    write(*,*) '/-----------------------------------------------------------\'
    write(*,*) '|  Ab initio dynamical vertex approximation program (ADGA)  |'
    write(*,*) '|  Running on ',mpi_wsize,' core(s).                        |'
    write(*,*) '\-----------------------------------------------------------/'
    write(*,*)
  end if

  ! creation of hdf5 output file
  if (mpi_wrank .eq. master) then
    ! generate a date-time string for output file name
    call date_and_time(date,time,zone,time_date_values)
    output_filename=trim(output_dir)//'adga-'//trim(date)//'-'//trim(time)//'-output.hdf5'
    ! while the name is generated here to get the correct starting time, 
    ! the file is created later, just before the beginning of the parallel loop.
    ! Thus, the parameters can be written immediately.
  end if

!##################  READ W2DYNAMICS HDF5 OUTPUT FILE ##################
! open the hdf5-fortran interface
  call init_h5()

! read bosonic and fermionic Matsubara axes iwf-g4,iwb-g4:
  call get_freq_range(mpi_wrank,master)
  call check_freq_range(mpi_wrank,master)

  if (mpi_wrank .eq. master) then
    write(*,*) 'iwmax=', iwmax, ' (number of fermionic matsubara frequencies of one-particle quantities)'
    write(*,*) 'iwfmax=', iwfmax, 'iwfmax_small=', iwfmax_small, ' (number of fermionic matsubara frequencies of two-particle quantities)'
    write(*,*)'iwbmax=',iwbmax, 'iwbmax_small=', iwbmax_small, ' (number of bosonic matsubara frequencies of two-particle quantities)'
  end if

! after frequencies and dimensions are obtained, arrays can be allocated 
  allocate(siw(-iwmax:iwmax-1,ndim))
  allocate(giw(-iwmax:iwmax-1,ndim))
  allocate(dc(2,ndim)) ! indices: spin band

! read Hamiltonian
  if(read_ext_hk) then
    call read_hk_w2w()
  else
    call read_hk_w2dyn()
  end if

  call read_siw()  ! w2d self energy
! this only works for no p-bands since read_giw reads the giw array which only
! contains correlated bands
  if(.not. exist_p) then
    call read_giw()  ! w2d greens function G_dmft
  endif

  call read_mu()   ! w2d chemical potential
  call read_beta() ! w2d inverse temperature
  call read_dc()   ! w2d double counting

  call finalize_h5() ! close the hdf5-fortran interface

  if (mpi_wrank .eq. master) then
    write(*,*) 'orb_symmetry = ', orb_sym
  end if

  if (.not. do_vq .and. mpi_wrank .eq. master) then
    write(*,*) 'Main: Run without V(q)'
  end if

  if (mpi_wrank .eq. master) then
    write(*,*) 'beta=', beta
    write(*,*) 'mu=', mu
    write(*,*) 'dc=', dc
  end if

  call init() ! this requires beta, so it has to be called after read_beta()

  !read umatrix from separate file:
  if (read_ext_u) then
    call read_u(u,u_tilde)
  else
    if (mpi_wrank .eq. master) write(*,*) 'Creating u matrix'
    call create_u(u,u_tilde)

    ! test umatrix
    if (mpi_wrank .eq. master) then
      open(unit=10,file=trim(output_dir)//"umatrix.dat")
      write(10,*) 'Umatrix File for the ADGA code : band,band,band,band,Uvalue'
      do i=1,ndim
      do j=1,ndim
      do k=1,ndim
      do l=1,ndim
        write(10,'(4I10,F15.8)') i,j,k,l,Umat(i,j,k,l)
      enddo
      enddo
      enddo
      enddo
    endif
    deallocate(Umat)
  endif

! COMPUTE local single-particle Greens function -- this overwrites the readin from w2d above
! This calculation is necessary if one wants to do calculations with p-bands
  call get_giw() ! writes giw_calc.dat

  allocate(gloc(-iwmax:iwmax-1,ndim,ndim))
  allocate(n_dmft(ndim), n_fock(nkp,ndim,ndim), n_dga(ndim))

!compute DMFT filling n_dmft
  call get_ndmft() ! writes n_dmft.dat
!compute k-dependent filling for Fock-term (computed in the EOM):
  call get_nfock() ! writes n_fock.dat


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


! define qw compound index for mpi:
  call mpi_distribute()



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



if (mpi_wrank .eq. master) then
  write(*,*) 'writing output to ',output_filename
  call init_h5_output(output_filename)
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
        chi_loc_magn_sum_left = 0.d0
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
        call calc_eom(interm3_dens,interm3_magn,gamma_dmft_dens,gamma_dmft_magn,gamma_loc_sum_left,sigma,sigma_dmft,sigma_hf,kq_ind,iwb,iq,v)
     end if
     call cpu_time(finish)

     !Output the calculation progress
     write(*,'("Core:",I5,"  Progress:",I7,"  of",I7,"  Time: ",F8.4,"  (working area:",I9,"  -",I9,")")') &
      mpi_wrank, iqw-qwstart+1, qwstop-qwstart+1, finish-start, qwstart, qwstop
       ! write(*,'((A),I5,2X,I6,2X,I6,2X,I6,2X,F8.4,(A),F8.4)') 'mpi_rank || qwstart -- iqw -- qwstop || % : ',&
                            ! mpi_wrank,qwstart,iqw,qwstop, (iqw-qwstart)/dble(qwstart-qwstop+1), 'time: ',finish-start
  enddo !iqw


!  close(267)

  deallocate(interm3_dens,interm3_magn)
  deallocate(giw,u,u_tilde,chi0_loc,chi0_loc_inv,chi0_sum,interm1,interm1_v,gamma_dmft_dens,gamma_dmft_magn)
  deallocate(chi_loc_dens_sum_left,chi_loc_magn_sum_left)
  deallocate(chi_loc_magn_full,chi_loc_dens_full,gamma_loc_sum_left,v,c_interm)
  deallocate(chi_loc_slice_dens,chi_loc_slice_magn)

#ifdef MPI
  ! MPI reduction and output
  if (do_eom) then
     allocate(sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp),sigma_sum_hf(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp))
     allocate(sigma_loc(ndim, ndim, -iwfmax_small:iwfmax_small-1),sigma_sum_dmft(ndim, ndim, -iwfmax_small:iwfmax_small-1))
     call MPI_reduce(sigma, sigma_sum, ndim*ndim*2*iwfmax_small*nkp, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
     call MPI_reduce(sigma_hf, sigma_sum_hf, ndim*ndim*2*iwfmax_small*nkp, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
     call MPI_reduce(sigma_dmft, sigma_sum_dmft, ndim*ndim*2*iwfmax_small, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
     sigma_sum = -sigma_sum/(beta*nqp)
     sigma_sum_hf = -sigma_sum_hf/(beta*nqp)
     sigma_sum_dmft = -sigma_sum_dmft/(beta*nqp)
     call add_siw_dmft(sigma_sum)
     call get_sigma_g_loc(sigma_sum, sigma_loc, gloc)
     if (mpi_wrank .eq. master) then
       call output_eom(sigma_sum, sigma_sum_dmft, sigma_sum_hf, sigma_loc, gloc)
       call output_eom_hdf5(output_filename,sigma_sum,sigma_sum_hf,sigma_loc,sigma_sum_dmft)
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

#endif
  deallocate(iw_data,iwb_data,siw,k_data,q_data,kq_ind,qw)
  if (mpi_wrank .eq. master) write(*,*) 'end of program'
  call mpi_close()
end program main
