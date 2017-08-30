module hdf5_module

  use hdf5
  use parameters_module
  implicit none

  integer(hid_t) :: plist_id
  integer(hid_t) :: dspace_id
  integer(hid_t) :: compound_id, type_r_id, type_i_id
  integer(size_t) :: compound_size, type_sized
  integer(hsize_t), dimension(2) :: dims
  double precision, allocatable :: tmp_r_1(:,:), tmp_i_1(:,:), tmp_err_1(:,:)
  integer :: err

 contains

!=====================================================================================
! wrapper function that can be called from outside without need to "use hdf5"
  subroutine init_h5()
    implicit none
    
    call h5open_f(err)
    call create_complex_datatype

  end subroutine init_h5

!=====================================================================================
! wrapper function that can be called from outside without need to "use hdf5"
  subroutine finalize_h5()
    implicit none
    
    call h5close_f(err)

  end subroutine finalize_h5

!=====================================================================================
   subroutine create_complex_datatype

     integer(size_t), parameter :: zero = 0

     ! Set dataset transfer property to preserve partially initialized fields during write/read to/from dataset with compound datatype (necessary?)
     call h5pcreate_f(h5p_dataset_xfer_f, plist_id, err)
     call h5pset_preserve_f(plist_id, .true., err)

     ! create compound datatype for complex arrays:
     call h5tget_size_f(h5t_native_double, type_sized, err)
     compound_size = 2*type_sized
     call h5tcreate_f(h5t_compound_f, compound_size, compound_id, err)
     call h5tinsert_f(compound_id, "r", zero, h5t_native_double, err)
     call h5tinsert_f(compound_id, "i", type_sized, h5t_native_double, err)

     !complex type to write real and imaginary individually:
     call h5tcreate_f(h5t_compound_f, type_sized, type_r_id, err)
     call h5tinsert_f(type_r_id, "r", zero, h5t_native_double, err)
     call h5tcreate_f(h5t_compound_f, type_sized, type_i_id, err)
     call h5tinsert_f(type_i_id, "i", zero, h5t_native_double, err)

   end subroutine create_complex_datatype
!=====================================================================================

!=====================================================================================
   subroutine read_axes(file_id, iwb_array, iwf_array, dspace_iwb_id, dspace_iwf_id, dim_iwb, dim_iwf)
     implicit none
     double precision, allocatable :: iwb_array(:),iwf_array(:)
     integer(hsize_t), dimension(1), intent(out) :: dim_iwb, dim_iwf
     integer(hsize_t), dimension(1) :: dim_iwf_max, dim_iwb_max
     integer(hid_t) :: axes_id,iwb_id,iwf_id
     integer(hid_t), intent(out) :: dspace_iwb_id, dspace_iwf_id
     integer(hid_t) :: file_id
 
     ! read fermionic Matsubara frequencies iwf:
     call h5dopen_f(file_id, ".axes/iwf-g4", iwf_id, err)
     call h5dget_space_f(iwf_id, dspace_iwf_id, err)
     call h5sget_simple_extent_dims_f(dspace_iwf_id, dim_iwf, dim_iwf_max, err)
     iwfmax = dim_iwf(1)/2
     allocate(iwf_array(-iwfmax:iwfmax-1))
     call h5dread_f(iwf_id, h5t_native_double, iwf_array, dim_iwf, err)
     call h5dclose_f(iwf_id, err)

     ! read bosonic Matsubara frequencies iwf:
     call h5dopen_f(file_id, ".axes/iwb-g4", iwb_id, err)
     call h5dget_space_f(iwb_id, dspace_iwb_id, err)
     call h5sget_simple_extent_dims_f(dspace_iwb_id, dim_iwb, dim_iwb_max, err)
     iwbmax = dim_iwb(1)/2
     allocate(iwb_array(-iwbmax:iwbmax))
     call h5dread_f(iwb_id, h5t_native_double, iwb_array, dim_iwb, err)
     call h5dclose_f(iwb_id, err)
   
   end subroutine read_axes


   subroutine write_axes(file_id, iwb_array, iwf_array, dspace_iwb_id, dspace_iwf_id, dim_iwb, dim_iwf)
     integer(hid_t) :: file_id
     double precision, intent(in):: iwb_array(-iwbmax:iwbmax), iwf_array(-iwfmax:iwfmax-1)
     integer(hsize_t),dimension(1) :: dim_iwb, dim_iwf
     integer(hid_t) :: axes_id, iwb_id, iwf_id
     integer(hid_t), intent(in) :: dspace_iwb_id, dspace_iwf_id

     !write Matsubara frequency axes:
     call h5gcreate_f(file_id, ".axes", axes_id, err)
     call h5dcreate_f(axes_id, "iwb-g4", h5t_native_double, dspace_iwb_id, iwb_id, err)
     call h5dcreate_f(axes_id, "iwf-g4", h5t_native_double, dspace_iwf_id, iwf_id, err)

     call h5dwrite_f(iwb_id, h5t_native_double, iwb_array, dim_iwb, err) 
     call h5dwrite_f(iwf_id, h5t_native_double, iwf_array, dim_iwf, err)
     
     call h5dclose_f(iwb_id, err)
     call h5dclose_f(iwf_id, err)
     call h5gclose_f(axes_id, err)

   end subroutine write_axes
!============================================================================================

!===========================================================================================
   subroutine create_channels(file_id, ineq)
     implicit none

     integer :: iwb,ineq
     integer(hid_t) :: grp_dens_id,grp_magn_id,iw_magn_id,iw_dens_id,ineq_id
     character(len=20) :: name_buffer
     character(len=20) :: grp_name
     integer(hid_t) :: file_id
     
     !create dens and magn groups:
     write(grp_name,'("ineq-",I3.3)') ineq
     call h5gcreate_f(file_id, grp_name, ineq_id, err)
     call h5gcreate_f(ineq_id, 'dens', grp_dens_id, err)
     call h5gcreate_f(ineq_id, 'magn', grp_magn_id, err)

     do iwb=0,2*iwbmax

        write(name_buffer, '(I5.5)'), iwb
        call h5gcreate_f(grp_dens_id, name_buffer, iw_magn_id, err)
        call h5gcreate_f(grp_magn_id, name_buffer, iw_dens_id, err)
        call h5gclose_f(iw_magn_id, err)
        call h5gclose_f(iw_dens_id, err)

     enddo

     call h5gclose_f(grp_dens_id, err)
     call h5gclose_f(grp_magn_id, err)

     return
   end subroutine create_channels
!=====================================================================================

!=====================================================================================
    subroutine create_component(file_id, ichannel, iwb, ind_orb, ineq)
      implicit none

      integer, intent(in) :: ineq ! inequivalent atom
      integer, intent(in) :: ichannel !1=magn, 2=dens
      integer, intent(in) :: iwb, ind_orb
      character(len=30) :: grpname
      integer(hid_t) :: grp_id, dset_id, dset_err_id
      integer(hid_t) :: file_id

      if (ichannel==1) then
         write(grpname, '("ineq-",I3.3,"/magn/",I5.5,"/",i5.5)') ineq, iwb, ind_orb
      else
         write(grpname, '("ineq-",I3.3,"/dens/",I5.5,"/",i5.5)') ineq, iwb, ind_orb
      endif

      call h5gcreate_f(file_id, trim(grpname), grp_id, err)
      call h5dcreate_f(grp_id, "value", compound_id, dspace_id, dset_id, err)
      call h5dcreate_f(grp_id, "error", h5t_native_double, dspace_id, dset_err_id, err)

      tmp_r_1 = 0.d0
      tmp_i_1 = 0.d0
      tmp_err_1 = 0.d0

      call h5dwrite_f(dset_id, type_r_id, tmp_r_1, dims, err)
      call h5dwrite_f(dset_id, type_i_id, tmp_i_1, dims, err)
      call h5dwrite_f(dset_err_id, h5t_native_double, tmp_err_1, dims, err)

      call h5dclose_f(dset_err_id, err)
      call h5dclose_f(dset_id, err)
      call h5gclose_f(grp_id, err)

    end subroutine create_component
!====================================================================================

!====================================================================================
   subroutine add_to_component(file_id, ichannel, iwb, ind_orb, g4iw_r, g4iw_i, g4err, ineq)
     implicit none
     integer, intent(in) :: ichannel, iwb, ind_orb, ineq
     double precision,intent(in) :: g4iw_r(2*iwfmax, 2*iwfmax, 2*iwbmax+1)
     double precision,intent(in) :: g4iw_i(2*iwfmax, 2*iwfmax, 2*iwbmax+1)
     double precision,intent(in) :: g4err(2*iwfmax, 2*iwfmax, 2*iwbmax+1)
     character(len=30) :: grpname
     integer(hid_t) :: grp_id, dset_id, dset_err_id
     integer(hid_t) :: file_id

     if (ichannel==1) then
        write(grpname, '("ineq-",I3.3,"/magn/",(I5.5),"/",(I5.5))') ineq, iwb ,ind_orb
     else
        write(grpname, '("ineq-",I3.3,"/dens/",(I5.5),"/",(I5.5))') ineq, iwb, ind_orb
     endif
        
     call h5gopen_f(file_id, trim(grpname), grp_id, err)
     call h5dopen_f(grp_id, "value", dset_id, err)
     call h5dopen_f(grp_id, "error", dset_err_id, err)
     
     call h5dread_f(dset_id, type_r_id, tmp_r_1, dims, err)
     call h5dread_f(dset_id, type_i_id, tmp_i_1, dims, err)
     call h5dread_f(dset_err_id, h5t_native_double, tmp_err_1, dims, err)

     tmp_r_1(:,:) = g4iw_r(:,:,iwb+1)+tmp_r_1(:,:)
     tmp_i_1(:,:) = g4iw_i(:,:,iwb+1)+tmp_i_1(:,:)
     tmp_err_1(:,:) = sqrt(g4err(:,:,iwb+1)**2+tmp_err_1(:,:)**2)

     call h5dwrite_f(dset_id, type_r_id, tmp_r_1, dims, err)
     call h5dwrite_f(dset_id, type_i_id, tmp_i_1, dims, err)
     call h5dwrite_f(dset_err_id, h5t_native_double, tmp_err_1, dims, err)

     call h5dclose_f(dset_err_id, err)
     call h5dclose_f(dset_id, err)
     call h5gclose_f(grp_id, err)

 end subroutine add_to_component
!===========================================================================================


 subroutine get_freq_range(mpi_wrank,master)
     implicit none
     integer :: mpi_wrank, master
     integer(hid_t) :: file_id
     integer(hsize_t), dimension(1) :: iw_dims,iw_maxdims,dim_iwb,dim_iwf,dim_iwf_max,dim_iwb_max
     integer(hid_t) :: iw_id,iwb_id,iwf_id,iw_space_id,dspace_iwb_id,dspace_iwf_id
 
     ! read frequency range of one-particle quantities
     if (mpi_wrank .eq. master) then
       write(*,*) 'one particle quantities in ',filename
     endif
     call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, err)
     call h5dopen_f(file_id, ".axes/iw", iw_id, err)
     call h5dget_space_f(iw_id, iw_space_id, err)
     call h5sget_simple_extent_dims_f(iw_space_id, iw_dims, iw_maxdims, err)
     iwmax = iw_dims(1)/2
     call h5dclose_f(iw_id, err)
     call h5fclose_f(file_id,err)

     if (mpi_wrank .eq. master) then
       write(*,*) 'two particle quantities in ',filename_vertex_sym
     endif
     call h5fopen_f(filename_vertex_sym, h5f_acc_rdonly_f, file_id, err)
     ! read fermionic Matsubara frequencies iwf:
     call h5dopen_f(file_id, ".axes/iwf-g4", iwf_id, err)
     call h5dget_space_f(iwf_id, dspace_iwf_id, err)
     call h5sget_simple_extent_dims_f(dspace_iwf_id, dim_iwf, dim_iwf_max, err)
     iwfmax = dim_iwf(1)/2
     call h5dclose_f(iwf_id, err)

     ! read bosonic Matsubara frequencies iwf:
     call h5dopen_f(file_id, ".axes/iwb-g4", iwb_id, err)
     call h5dget_space_f(iwb_id, dspace_iwb_id, err)
     call h5sget_simple_extent_dims_f(dspace_iwb_id, dim_iwb, dim_iwb_max, err)
     iwbmax = dim_iwb(1)/2
     call h5dclose_f(iwb_id, err)
     call h5fclose_f(file_id,err)

 end subroutine get_freq_range



 subroutine read_hk_w2w()
     implicit none
     double precision, allocatable :: hr(:,:), hi(:,:)
     integer :: ik,i,j
     integer :: ndim_test
     double precision :: kx, ky, kz
     
     open(21, file=filename_hk, status='unknown') ! the filename_hk is taken from parameters_module
     read(21,*) nkp,ndim_test
     if (ndim .ne. ndim_test) then
       stop "Error: number of dimensions in input file do not coincide with Hamiltonian file"
     endif
     allocate(hr(ndim,ndim),hi(ndim,ndim))
     allocate(hk(ndim,ndim,nkp))
     allocate(k_data(3,nkp))
 
     do ik=1,nkp
        read(21,*) kx,ky,kz
        k_data(1,ik) = kx
        k_data(2,ik) = ky
        k_data(3,ik) = kz
        do i=1,ndim
           read(21,*) (hr(i,j),hi(i,j),j=1,ndim)
        enddo
        hk(:,:,ik)=hr(:,:)+ci*hi(:,:)
     enddo
 
     close(21)

 end subroutine read_hk_w2w

 subroutine read_siw()
     implicit none
     integer :: i,iw,ineq,iband,dimstart, dimend
     integer(hid_t) :: file_id,siw_id, siw_space_id
     integer(hsize_t), dimension(3) :: siw_dims, siw_maxdims
     double precision, allocatable :: siw_data(:,:,:,:)
     double precision :: siw_r, siw_i
     character(len=200) :: name_buffer

     siw=0.d0
   
     do ineq=1,nineq
       dimstart=1
       do i=2,ineq
         dimstart=dimstart+ndims(i-1,1)+ndims(i-1,2)
       enddo
       dimend=dimstart+ndims(ineq,1)-1 ! here we are only interested in the interacting orbitals
   
       write(name_buffer,'("ineq-",I3.3)') ineq
       ! read siw:
       ! local self energy - only for interacting orbitals == d
       call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, err)
       call h5dopen_f(file_id, "stat-001/"//trim(name_buffer)//"/siw/value", siw_id, err)
       call h5dget_space_f(siw_id, siw_space_id, err)
       call h5sget_simple_extent_dims_f(siw_space_id, siw_dims, siw_maxdims, err)
       ! ndims = siw_dims(3)
       allocate(siw_data(2,-iwmax:iwmax-1,siw_dims(2),siw_dims(3))) !indices: real/imag iw spin band
       call h5dread_f(siw_id, compound_id, siw_data, siw_dims, err)
   
       !paramagnetic (spin average):
       do i=dimstart,dimend
         siw(:,i) = siw_data(1,:,1,i-dimstart+1)+siw_data(1,:,2,i-dimstart+1)+ci*siw_data(2,:,1,i-dimstart+1)+ci*siw_data(2,:,2,i-dimstart+1)
         siw(:,i) = siw(:,i)/2.d0
       enddo
   
       call h5dclose_f(siw_id, err)
       call h5fclose_f(file_id,err)
       deallocate(siw_data)
   
       if (orb_sym) then
          ! enforce orbital symmetry:
          do iband=dimstart+1,dimend
             siw(:,dimstart) = siw(:,dimstart)+siw(:,iband)
          enddo
          siw(:,dimstart)=siw(:,dimstart)/dble(dimend-dimstart+1)
          do iband=dimstart+1,dimend
             siw(:,iband) = siw(:,dimstart)
          enddo
       endif
     enddo ! loop over inequivalent atoms
   
     ! test siw:
     open(34, file=trim(output_dir)//"siw.dat", status='unknown')
     do iw=-iwmax,iwmax-1
        write(34,'(100F12.6)') (real(siw(iw,i)),aimag(siw(iw,i)), i=1,ndim)
     enddo
     close(34)
   

 end subroutine read_siw

 subroutine read_giw()
     implicit none
     integer :: i,iw,ineq,iband,dimstart, dimend
     integer(hid_t) :: file_id,giw_id, giw_space_id
     integer(hsize_t), dimension(3) :: giw_dims, giw_maxdims
     double precision, allocatable :: giw_data(:,:,:,:)
     double precision :: giw_r, giw_i
     character(len=200) :: name_buffer

     giw=0.d0
   
     do ineq=1,nineq
       dimstart=1
       do i=2,ineq
         dimstart=dimstart+ndims(i-1,1)+ndims(i-1,2)
       enddo
       dimend=dimstart+ndims(ineq,1)+ndims(ineq,2)-1 ! here we are only interested in the interacting orbitals
   
       write(name_buffer,'("ineq-",I3.3)') ineq
       !read giw
       call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, err)
       call h5dopen_f(file_id, "stat-001/"//trim(name_buffer)//"/giw/value", giw_id, err)
       call h5dget_space_f(giw_id, giw_space_id, err)
       call h5sget_simple_extent_dims_f(giw_space_id, giw_dims, giw_maxdims, err)
       allocate(giw_data(2,-iwmax:iwmax-1,giw_dims(2),giw_dims(3))) !indices: real/imag iw spin band
       call h5dread_f(giw_id, compound_id, giw_data, giw_dims, err)
   
       !paramagnetic:
       do i=dimstart,dimend
         giw(:,i) = giw_data(1,:,1,i-dimstart+1)+giw_data(1,:,2,i-dimstart+1)+ci*giw_data(2,:,1,i-dimstart+1)+ci*giw_data(2,:,2,i-dimstart+1)
         giw(:,i) = giw(:,i)/2.d0
       enddo
   
       call h5dclose_f(giw_id, err)
       call h5fclose_f(file_id,err)
       deallocate(giw_data)
   
   
       if (orb_sym) then
       ! enforce orbital symmetry:
           ! here we need to enforce symmetry over one type of band specifically
           dimstart=1
           do i=2,ineq
             dimstart=dimstart+ndims(i-1,1)+ndims(i-1,2)
           enddo
   
           do i=1,2 ! d and p bands
             if (ndims(ineq,i) .eq. 0) cycle ! do nothing
             if (i .eq. 1) then
               dimend = dimstart+ndims(ineq,1)-1
             endif
             if (i .eq. 2) then
               dimend = dimstart+ndims(ineq,1)+ndims(ineq,2)-1
               dimstart = dimend-ndims(ineq,2)+1
             endif
   
             do iband=dimstart+1,dimend
               giw(:,dimstart) = giw(:,dimstart)+giw(:,iband)
             enddo
             giw(:,dimstart)=giw(:,dimstart)/dble(dimend-dimstart+1)
             do iband=dimstart+1,dimend
               giw(:,iband) = giw(:,dimstart)
             enddo
           enddo
       endif
   
     enddo ! inequivalent atom loop

  ! test giw:
  open(54, file=trim(output_dir)//"giw.dat", status='unknown')
  do iw=-iwmax,iwmax-1
     write(54,'(100F12.6)') (real(giw(iw,i)),aimag(giw(iw,i)),i=1,ndim)
  enddo
  close(54)
 end subroutine read_giw

 subroutine read_hk_w2dyn()
     implicit none
     integer :: ndim_test
     integer :: i,ik
     integer(hid_t) :: file_id,k_id,k_space_id,hk_id,hk_space_id
     integer(hsize_t), dimension(2) :: k_dims, k_maxdims
     integer(hsize_t), dimension(3) :: hk_dims, hk_maxdims
     double precision, allocatable :: hk_data(:,:,:,:)

    ! read k-points:
    call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, err)
    call h5dopen_f(file_id, ".axes/k-points", k_id, err)
    call h5dget_space_f(k_id, k_space_id, err)
    call h5sget_simple_extent_dims_f(k_space_id, k_dims, k_maxdims, err)
    nkp = k_dims(2)
    allocate(k_data(k_dims(1),k_dims(2))) !indices: 3 ik
    call h5dread_f(k_id, h5t_native_double, k_data, k_dims, err)
    call h5dclose_f(k_id, err)

! write k-points:
!    open(37, file=trim(output_dir)//'k_points.dat', status='unknown')
!    do ik=1,100
!      write(37,'(100F12.6)') k_data(2,ik), k_data(3,ik)
!    enddo
!    close(37)

    call h5dopen_f(file_id, "start/hk/value", hk_id, err)
    call h5dget_space_f(hk_id, hk_space_id, err)
    call h5sget_simple_extent_dims_f(hk_space_id, hk_dims, hk_maxdims, err)
    ndim_test = hk_dims(1)
    if (ndim .ne. ndim_test) then
      stop "Error: number of dimensions in input file do not coincide with Hamiltonian file"
    endif
    allocate(hk_data(2,hk_dims(1),hk_dims(2),hk_dims(3)))
    call h5dread_f(hk_id, compound_id, hk_data, hk_dims, err)
    allocate(hk(hk_dims(1),hk_dims(2),hk_dims(3))) !indices: band band ik
    hk = 0.d0
    hk(:,:,:) = hk_data(1,:,:,:)+ci*hk_data(2,:,:,:)
    call h5dclose_f(hk_id, err)
    deallocate(hk_data)

  ! test hk:
    !open(34, file=trim(output_dir)//"hk.dat", status='unknown')
    !do ik=1,hk_dims(3)
    !   write(34,*)k_data(:,ik)
    !   do i=1,hk_dims(2)
    !      write(34,'(100F12.6)')hk(:,i,ik)
    !   enddo
    !enddo
    !close(34)
    call h5fclose_f(file_id,err)   
 
 end subroutine read_hk_w2dyn

 subroutine read_beta()
     implicit none
     integer(hid_t) :: file_id,config_id,beta_id
     integer(hsize_t),dimension(0) :: beta_dims

     call h5fopen_f(filename,h5f_acc_rdonly_f,file_id,err)
     call h5gopen_f(file_id, ".config", config_id, err)
     call h5aopen_f(config_id, "general.beta", beta_id, err)
     call h5aread_f(beta_id, h5t_native_double, beta, beta_dims, err)
     call h5aclose_f(beta_id, err)
     call h5gclose_f(config_id,err)
     call h5fclose_f(file_id, err)
 end subroutine read_beta

 subroutine read_mu()
     implicit none
     integer(hid_t) :: mu_id,file_id,mu_space_id
     integer(hsize_t),dimension(0) :: mu_dims

     call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, err)
     call h5dopen_f(file_id, "stat-001/mu/value", mu_id, err)
     call h5dread_f(mu_id, h5t_native_double, mu, mu_dims, err)
     call h5dclose_f(mu_id, err)
     call h5fclose_f(file_id,err)
 end subroutine read_mu


 subroutine read_dc()
     implicit none
     integer :: ineq,i,iband,dimstart,dimend
     integer(hid_t) :: file_id,dc_id, dc_space_id
     integer(hsize_t), dimension(2) :: dc_dims, dc_maxdims
     double precision, allocatable :: dc_data(:,:)
     character(len=200) :: name_buffer

     ! dc for noninteracting bands set to 0
     dc = 0.d0

!     do ineq=1,nineq
!       dimstart=1
!       do i=2,ineq
!         dimstart=dimstart+ndims(i-1,1)+ndims(i-1,2)
!       enddo
!       dimend=dimstart+ndims(ineq,1)-1 ! here we are only interested in the interacting orbitals
!       write(name_buffer,'("ineq-",I3.3)') ineq
!       call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, err)
!       call h5dopen_f(file_id, "stat-001/"//trim(name_buffer)//"/dc/value", dc_id, err)
!       call h5dget_space_f(dc_id, dc_space_id, err)
!       call h5sget_simple_extent_dims_f(dc_space_id, dc_dims, dc_maxdims, err)
!       allocate(dc_data(dc_dims(1),dc_dims(2))) !indices: spin band
!       call h5dread_f(dc_id, h5t_native_double, dc_data, dc_dims, err)
!       call h5dclose_f(dc_id, err)
!       call h5fclose_f(file_id,err)
!
!       do iband=dimstart,dimend
!         dc(:,iband) = dc_data(:,iband-dimstart+1)
!       enddo
!       deallocate(dc_data)
!     enddo
 end subroutine read_dc

 subroutine read_vertex(chi_loc_dens_full,chi_loc_magn_full,iwb)
     use aux
     implicit none
     integer :: ineq,dimstart,dimend,imembers,ind_grp,b1,b2,b3,b4,ind_iwb
     integer :: i1,i2,iwf1,iwf2,i,j,k,l
     integer,intent(in) :: iwb
     integer(hid_t) :: file_vert_id,grp_magn_id,grp_dens_id,nmembers,dset_magn_id,dset_dens_id,itype
     integer(hsize_t), dimension(2) :: tmp_dims
     character(len=100) :: grpname_magn,grpname_dens,name_buffer,name_buffer_dset
     complex(kind=8), allocatable :: g4iw_magn(:,:,:,:,:,:), g4iw_dens(:,:,:,:,:,:)
     complex(kind=8),intent(out) :: chi_loc_magn_full(maxdim,maxdim),chi_loc_dens_full(maxdim,maxdim)
     double precision, allocatable :: tmp_r(:,:), tmp_i(:,:)

        allocate(g4iw_magn(ndim, ndim, -iwfmax:iwfmax-1, ndim, ndim, -iwfmax:iwfmax-1))
        allocate(g4iw_dens(ndim, ndim, -iwfmax:iwfmax-1, ndim, ndim, -iwfmax:iwfmax-1))

        allocate(tmp_r(-iwfmax:iwfmax-1, -iwfmax:iwfmax-1))
        allocate(tmp_i(-iwfmax:iwfmax-1, -iwfmax:iwfmax-1))
        tmp_dims = (/2*iwfmax, 2*iwfmax/)

        g4iw_magn = 0.d0
        g4iw_dens = 0.d0
        ind_iwb = iwb+iwbmax

        do ineq=1,nineq
          write(grpname_magn, '("ineq-",I3.3,"/magn/",(I5.5))'), ineq, ind_iwb
          write(grpname_dens, '("ineq-",I3.3,"/dens/",(I5.5))'), ineq, ind_iwb

          dimstart=1
          do i=2,ineq
            dimstart=dimstart+ndims(i-1,1)+ndims(i-1,2)
          enddo
          dimend=dimstart+ndims(ineq,1)+ndims(ineq,2)-1

          call h5open_f(err)
          call h5fopen_f(filename_vertex_sym, h5f_acc_rdonly_f, file_vert_id, err)
          call h5gopen_f(file_vert_id, grpname_magn, grp_magn_id, err)
          call h5gopen_f(file_vert_id, grpname_dens, grp_dens_id, err)

          call h5gn_members_f(file_vert_id, grpname_magn, nmembers, err)

          do imembers=0,nmembers-1

             call h5gget_obj_info_idx_f(file_vert_id, grpname_magn, imembers, name_buffer, itype, err)
             read(name_buffer,'(I5.5)')ind_grp

             call index2component_band(dimend-dimstart+1, ind_grp, b1, b2, b3, b4)

             write(name_buffer_dset, '("ineq-",I3.3,"/magn/",(I5.5),"/",(I5.5),"/value")') ineq, ind_iwb, ind_grp
             call h5dopen_f(file_vert_id, name_buffer_dset, dset_magn_id, err)
             call h5dread_f(dset_magn_id, type_r_id, tmp_r, tmp_dims, err)
             call h5dread_f(dset_magn_id, type_i_id, tmp_i, tmp_dims, err)

             g4iw_magn(dimstart+b1-1,dimstart+b2-1,:,dimstart+b3-1,dimstart+b4-1,:) = (tmp_r(:,:)+ci*tmp_i(:,:))*beta

             call h5dclose_f(dset_magn_id, err)

             write(name_buffer_dset, '("ineq-",I3.3,"/dens/",(I5.5),"/",(I5.5),"/value")') ineq, ind_iwb, ind_grp
             call h5dopen_f(file_vert_id, name_buffer_dset, dset_dens_id, err)
             call h5dread_f(dset_dens_id, type_r_id, tmp_r, tmp_dims, err)
             call h5dread_f(dset_dens_id, type_i_id, tmp_i, tmp_dims, err)

             g4iw_dens(dimstart+b1-1,dimstart+b2-1,:,dimstart+b3-1,dimstart+b4-1,:) = (tmp_r(:,:)+ci*tmp_i(:,:))*beta

             call h5dclose_f(dset_dens_id, err)


          enddo ! members

          call h5gclose_f(grp_dens_id, err)
          call h5gclose_f(grp_magn_id, err)
          call h5fclose_f(file_vert_id, err)
          call h5close_f(err)

        enddo ! loop for inequivalent atoms

        write(*,*) "Reading Vertex complete"

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

                          ! Depending on the type of vertex read, disconnected contributions need to be added/removed in order to
                          ! obtain \chi
                          if (vertex_type .eq. full_g4) then
                            !full 2-particle GF:
                            !straight term G(\nu)G(\nu') is subtracted (twice) only in the dens channel and only for iw=0:
                            if((iwb .eq. 0) .and. i==j .and. k==l)then
                              if(index2ineq(nineq,ndims,i,j,k,l)) then ! substracted in the correlated subspace
                                chi_loc_dens_full(i1,i2) = chi_loc_dens_full(i1,i2)-2.d0*beta*giw(iwf1,i)*giw(iwf2,l)
                              endif
                            endif

                            if((iwf2 .eq. iwf1) .and. i==l .and. j==k)then
                              if(.not. index2ineq(nineq,ndims,i,j,k,l)) then ! add bubble term only if not in the same correlated subspace
                                chi_loc_dens_full(i1,i2) = chi_loc_dens_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j)
                                chi_loc_magn_full(i1,i2) = chi_loc_magn_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j)
                              endif
                            endif

                          else if (vertex_type .eq. connected_g4) then
                            !G_conn:
                            !bubble term -G(\nu)G(\nu-\omega) is added in both channels
                            if((iwf2 .eq. iwf1) .and. i==l .and. j==k)then ! add all possible bubble terms
                              chi_loc_dens_full(i1,i2) = chi_loc_dens_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j)
                              chi_loc_magn_full(i1,i2) = chi_loc_magn_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j)
                            endif
                          else if (vertex_type .eq. chi_g4) then

                            if((iwf2 .eq. iwf1) .and. i==l .and. j==k)then
                              if(.not. index2ineq(nineq,ndims,i,j,k,l)) then ! add bubble term only if not in the same correlated subspace
                                chi_loc_dens_full(i1,i2) = chi_loc_dens_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j)
                                chi_loc_magn_full(i1,i2) = chi_loc_magn_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j)
                              endif
                            endif

                          endif

                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo

        deallocate(g4iw_magn, g4iw_dens, tmp_r, tmp_i)

 end subroutine read_vertex



! This subroutine creates the HDF5 output file and initializes its structure
subroutine init_h5_output(filename_output)
  implicit none

  integer :: err,i1,i2,ikx,iky,ikz
  integer(hid_t) :: file_id,grp_id_chi_qw,grp_id_chi_loc,grp_id_siwk,grp_id_input,grp_id_susc,grp_id_se
  character(len=*) :: filename_output
  integer(kind=8) :: chi_qw_dims(3),chi_loc_dims(3)
  integer(hid_t) :: dspace_scalar_id,dset_id_mu,dset_id_beta,dset_id_1,dspace_vec_id
  integer(hid_t) :: dspace_dc_id,dset_dc_id,dspace_id_g,dspace_id_s,dset_id_g,dset_id_s,dset_id_hk,dspace_id_hk
  integer(hsize_t),dimension(0) :: dims_scalar
  integer(hsize_t),dimension(1) :: dims_vec
  integer(hsize_t),dimension(2) :: dims_dc,dims_g
  integer(hsize_t),dimension(5) :: dims_hk_arr
  complex(kind=8), dimension(:,:,:,:,:),allocatable :: hk_arr
  write(*,*) 'initialize output file'

  call h5open_f(err)
  call h5fcreate_f(filename_output, h5f_acc_trunc_f, file_id, err)

! create the groups for the ADGA results
  if (do_chi) then
    call h5gcreate_f(file_id,'susceptibility',grp_id_susc,err)
    call h5gcreate_f(grp_id_susc,'nonloc',grp_id_chi_qw,err)
    call h5gclose_f(grp_id_chi_qw,err)
    call h5gcreate_f(grp_id_susc,'loc',grp_id_chi_loc,err)
    call h5gclose_f(grp_id_chi_loc,err)
    call h5gclose_f(grp_id_susc,err)
  end if
  
  if (do_eom) then
    call h5gcreate_f(file_id,'selfenergy',grp_id_se,err)
    call h5gclose_f(grp_id_se,err)
  end if

! create the group for parameters and input
  call h5gcreate_f(file_id,'input',grp_id_input,err)

! create dataspace for scalar quantities
  call h5screate_simple_f(0,dims_scalar,dspace_scalar_id,err)
  
! write chemical potential
  call h5dcreate_f(grp_id_input,'mu',H5T_NATIVE_DOUBLE,dspace_scalar_id,dset_id_mu,err)
  call h5dwrite_f(dset_id_mu,H5T_NATIVE_DOUBLE,mu,dims_scalar,err)
  call h5dclose_f(dset_id_mu,err)

! write inverse temperature
  call h5dcreate_f(grp_id_input,'beta',H5T_NATIVE_DOUBLE,dspace_scalar_id,dset_id_beta,err)
  call h5dwrite_f(dset_id_beta,H5T_NATIVE_DOUBLE,beta,dims_scalar,err)
  call h5dclose_f(dset_id_beta,err)
  
! create dataspace for double counting
  dims_dc=(/2,ndim/)
  call h5screate_simple_f(2,dims_dc,dspace_dc_id,err)

! write double counting
  call h5dcreate_f(grp_id_input,'dc',compound_id,dspace_dc_id,dset_dc_id,err)
  call h5dwrite_f(dset_dc_id,type_r_id,real(dc),dims_dc,err)
  call h5dwrite_f(dset_dc_id,type_i_id,aimag(dc),dims_dc,err)
  call h5dclose_f(dset_dc_id,err)

! create dataspace for siw and giw
  dims_g=(/2*iwmax,ndim/)
  call h5screate_simple_f(2,dims_g,dspace_id_g,err)

! write DMFT one-particle Green's function
  call h5dcreate_f(grp_id_input,'giw',compound_id,dspace_id_g,dset_id_g,err)
  call h5dwrite_f(dset_id_g,type_r_id,real(giw),dims_g,err)
  call h5dwrite_f(dset_id_g,type_i_id,aimag(giw),dims_g,err)
  call h5dclose_f(dset_id_g,err)

! write DMFT self-energy
  call h5dcreate_f(grp_id_input,'siw',compound_id,dspace_id_g,dset_id_s,err)
  call h5dwrite_f(dset_id_s,type_r_id,real(siw),dims_g,err)
  call h5dwrite_f(dset_id_s,type_i_id,aimag(siw),dims_g,err)
  call h5dclose_f(dset_id_s,err)

! create dataspace for hamiltonian
  dims_hk_arr=(/ ndim,ndim,nkpz,nkpy,nkpx /)
  allocate(hk_arr(ndim,ndim,nkpz,nkpy,nkpx))
  do ikx=1,nkpx
    do iky=1,nkpy
      do ikz=1,nkpz
        do i1=1,ndim
          do i2=1,ndim
            hk_arr(i2,i1,ikz,iky,ikx)=hk(i1,i2,1+ikz+iky*nkpz+ikx*nkpy*nkpz)
          end do
        end do
      end do
    end do
  end do
  call h5screate_simple_f(5,dims_hk_arr,dspace_id_hk,err)

! write hamiltonian
  call h5dcreate_f(grp_id_input,'hk',compound_id,dspace_id_hk,dset_id_hk,err)
  call h5dwrite_f(dset_id_hk,type_r_id,real(hk_arr),dims_hk_arr,err)
  call h5dwrite_f(dset_id_hk,type_i_id,aimag(hk_arr),dims_hk_arr,err)
  call h5dclose_f(dset_id_hk,err)
  



! write iw,iwfmax,iwbmax,nkp,nqp etc.
  call h5dcreate_f(grp_id_input,'iwmax',H5T_NATIVE_INTEGER,dspace_scalar_id,dset_id_1,err)
  call h5dwrite_f(dset_id_1,H5T_NATIVE_INTEGER,iwmax,dims_scalar,err)
  call h5dclose_f(dset_id_1,err)
  call h5dcreate_f(grp_id_input,'iwfmax',H5T_NATIVE_INTEGER,dspace_scalar_id,dset_id_1,err)
  call h5dwrite_f(dset_id_1,H5T_NATIVE_INTEGER,iwfmax,dims_scalar,err)
  call h5dclose_f(dset_id_1,err)
  call h5dcreate_f(grp_id_input,'iwbmax',H5T_NATIVE_INTEGER,dspace_scalar_id,dset_id_1,err)
  call h5dwrite_f(dset_id_1,H5T_NATIVE_INTEGER,iwbmax,dims_scalar,err)
  call h5dclose_f(dset_id_1,err)
  call h5dcreate_f(grp_id_input,'iwfmax_small',H5T_NATIVE_INTEGER,dspace_scalar_id,dset_id_1,err)
  call h5dwrite_f(dset_id_1,H5T_NATIVE_INTEGER,iwfmax_small,dims_scalar,err)
  call h5dclose_f(dset_id_1,err)
  call h5dcreate_f(grp_id_input,'iwbmax_small',H5T_NATIVE_INTEGER,dspace_scalar_id,dset_id_1,err)
  call h5dwrite_f(dset_id_1,H5T_NATIVE_INTEGER,iwbmax_small,dims_scalar,err)
  call h5dclose_f(dset_id_1,err)

  call h5dcreate_f(grp_id_input,'nkp',H5T_NATIVE_INTEGER,dspace_scalar_id,dset_id_1,err)
  call h5dwrite_f(dset_id_1,H5T_NATIVE_INTEGER,nkp,dims_scalar,err)
  call h5dclose_f(dset_id_1,err)
  call h5dcreate_f(grp_id_input,'nqp',H5T_NATIVE_INTEGER,dspace_scalar_id,dset_id_1,err)
  call h5dwrite_f(dset_id_1,H5T_NATIVE_INTEGER,nqp,dims_scalar,err)
  call h5dclose_f(dset_id_1,err)
  
  dims_vec=(/ 3 /)
  call h5screate_simple_f(1,dims_vec,dspace_vec_id,err)
  call h5dcreate_f(grp_id_input,'nkpxyz',H5T_NATIVE_INTEGER,dspace_vec_id,dset_id_1,err)
  call h5dwrite_f(dset_id_1,H5T_NATIVE_INTEGER,(/nkpx,nkpy,nkpz/),dims_vec,err)
  call h5dclose_f(dset_id_1,err)
  call h5dcreate_f(grp_id_input,'nqpxyz',H5T_NATIVE_INTEGER,dspace_vec_id,dset_id_1,err)
  call h5dwrite_f(dset_id_1,H5T_NATIVE_INTEGER,(/nqpx,nqpy,nqpz/),dims_vec,err)
  call h5dclose_f(dset_id_1,err)


! close the group and the file
  call h5gclose_f(grp_id_input,err)
  call h5fclose_f(file_id,err)

  write(*,*) 'hdf5 output initialized'
end subroutine init_h5_output


subroutine output_chi_loc_h5(filename_output,channel,chi_loc)
  implicit none
  character(len=*) :: filename_output,channel
  integer :: err,rank_chi_loc
  integer :: i1,i2,i3,i4,iwb
  integer(kind=8),dimension(:),allocatable :: dims_chi_loc
  integer(hid_t) :: file_id,grp_id_chi_loc
  integer(hid_t) :: dspace_id_chi_loc
  integer(hid_t) :: dset_id_chi_loc
  complex(kind=8),dimension(ndim**2,ndim**2,2*iwbmax_small+1) :: chi_loc
  complex(kind=8),dimension(ndim,ndim,ndim,ndim,2*iwbmax_small+1) :: chi_loc_outputarray

  rank_chi_loc=5
  allocate(dims_chi_loc(rank_chi_loc))
  dims_chi_loc=(/ ndim,ndim,ndim,ndim,2*iwbmax_small+1 /)


  ! reshape and transpose the array, i.e. break up the compound index
  do i1=1,ndim
    do i2=1,ndim
      do i3=1,ndim
        do i4=1,ndim
          do iwb=1,2*iwbmax_small+1
            chi_loc_outputarray(i4,i3,i2,i1,iwb) = chi_loc(ndim*(i1-1)+i2,ndim*(i4-1)+i3,iwb)
          end do
        end do
      end do
    end do
  end do



  call h5open_f(err)
  call h5fopen_f(filename_output,H5F_ACC_RDWR_F,file_id,err)
  call h5gopen_f(file_id,'susceptibility/loc',grp_id_chi_loc,err)

  call h5screate_f(H5S_SIMPLE_F,dspace_id_chi_loc,err)
  call h5sset_extent_simple_f(dspace_id_chi_loc,rank_chi_loc,dims_chi_loc,dims_chi_loc,err)

  call h5dcreate_f(grp_id_chi_loc,channel,compound_id,dspace_id_chi_loc,dset_id_chi_loc,err)
  call h5dwrite_f(dset_id_chi_loc,type_r_id,real(chi_loc_outputarray),dims_chi_loc,err)
  call h5dwrite_f(dset_id_chi_loc,type_i_id,aimag(chi_loc_outputarray),dims_chi_loc,err)

  call h5dclose_f(dset_id_chi_loc,err)
  call h5gclose_f(grp_id_chi_loc,err)
  call h5close_f(err)
 
  write(*,*) 'local susceptibility ',channel,' written to h5'
end subroutine output_chi_loc_h5


subroutine output_chi_qw_h5(filename_output,channel,chi_qw)

  implicit none
  character(len=*) :: filename_output,channel
  integer :: err,rank_chi_qw
  integer :: iwb,iqx,iqy,iqz,i1,i2,i3,i4
  integer(kind=8),dimension(:),allocatable :: dims_chi_qw,dims_chi_slice,offset_chi_slice,stride,block
  integer(hid_t) :: file_id,grp_id_chi_qw
  integer(hid_t) :: dspace_id_chi_qw,dspace_id_chi_slice
  integer(hid_t) :: dset_id_chi_qw
  complex(kind=8),dimension(ndim**2,ndim**2,nqp*(2*iwbmax_small+1)) :: chi_qw
  complex(kind=8),dimension(ndim,ndim,ndim,ndim,nqpz,nqpy,nqpx) :: chi_slice

  rank_chi_qw=8 
  allocate(dims_chi_qw(rank_chi_qw))
  dims_chi_qw=(/ ndim,ndim,ndim,ndim,nqpz,nqpy,nqpx,2*iwbmax_small+1 /)
  allocate(dims_chi_slice(rank_chi_qw))
  dims_chi_slice=(/ ndim,ndim,ndim,ndim,nqpz,nqpy,nqpx,1/)
  allocate(stride(rank_chi_qw),block(rank_chi_qw),offset_chi_slice(rank_chi_qw))
  stride=(/1,1,1,1,1,1,1,1/)
  block=(/1,1,1,1,1,1,1,1/)
  chi_slice=1.d0

  call h5open_f(err)
  call h5fopen_f(filename_output,H5F_ACC_RDWR_F,file_id,err)
  call h5gopen_f(file_id,'susceptibility/nonloc',grp_id_chi_qw,err)

  call h5screate_f(H5S_SIMPLE_F,dspace_id_chi_qw,err)
  call h5sset_extent_simple_f(dspace_id_chi_qw,rank_chi_qw,dims_chi_qw,dims_chi_qw,err)

  call h5dcreate_f(grp_id_chi_qw,channel,compound_id,dspace_id_chi_qw,dset_id_chi_qw,err)
  call h5dwrite_f(dset_id_chi_qw,type_r_id,real(chi_qw),dims_chi_qw,err) ! this is overwritten anyway
  call h5dwrite_f(dset_id_chi_qw,type_i_id,aimag(chi_qw),dims_chi_qw,err)
  call h5dclose_f(dset_id_chi_qw,err)

  call h5gclose_f(grp_id_chi_qw,err)

  ! re-open group
  call h5gopen_f(file_id,'susceptibility/nonloc',grp_id_chi_qw,err)

  ! since fortran knows only 7-dimensional arrays, we cannot write the 8-dimensional chi directly
  ! instead it is written slice-by-slice in a loop over omega
  ! additionally, the order of the indices is reversed here, because fortran uses column-major memory layout
  ! Furthermore, the last two band indices are swapped back here to break up the compound index correctly.
  ! TODO: check, if really everything is ok with the indices!
  do iwb=1,2*iwbmax_small+1
    do iqz=1,nqpz
      do iqy=1,nqpy
        do iqx=1,nqpx
          do i1=1,ndim
            do i2=1,ndim
              do i3=1,ndim
                do i4=1,ndim
                  chi_slice(i4,i3,i2,i1,iqz,iqy,iqx) = chi_qw(ndim*(i1-1)+i2,ndim*(i4-1)+i3,(iwb-1)*nqp+(iqz-1)+(iqy-1)*nqpz+(iqx-1)*nqpy*nqpz+1)
                end do ! i4
              end do ! i3
            end do ! i2
          end do ! i1
        end do ! nqpx
      end do ! nqpy
    end do ! nqpz
    ! re-open dataset
    call h5dopen_f(grp_id_chi_qw,channel,dset_id_chi_qw,err)
    offset_chi_slice=(/0,0,0,0,0,0,0,iwb-1/)
    
    ! select hyperslab
    call h5sselect_hyperslab_f(dspace_id_chi_qw,H5S_SELECT_SET_F,offset_chi_slice,dims_chi_slice,err,stride,block)
    ! create data space for subset
    call h5screate_f(H5S_SIMPLE_F,dspace_id_chi_slice,err)
    call h5sset_extent_simple_f(dspace_id_chi_slice,rank_chi_qw,dims_chi_slice,dims_chi_slice,err)

    call h5dwrite_f(dset_id_chi_qw,type_r_id,real(chi_slice),dims_chi_slice,err,dspace_id_chi_slice,dspace_id_chi_qw)
    call h5dwrite_f(dset_id_chi_qw,type_i_id,aimag(chi_slice),dims_chi_slice,err,dspace_id_chi_slice,dspace_id_chi_qw)
    call h5dclose_f(dset_id_chi_qw,err)
  end do ! iwb

  call h5gclose_f(grp_id_chi_qw,err)
  call h5close_f(err)
 
  write(*,*) 'nonlocal susceptibility ',channel,' written to h5'
end subroutine output_chi_qw_h5


subroutine output_eom_hdf5(filename_output,sigma_sum,sigma_sum_hf,sigma_loc,sigma_sum_dmft)
  implicit none
  character(len=*) :: filename_output
  integer(hid_t) :: file_id,grp_id_siwk,dset_id_sigmasum,dspace_id_sigmasum,dset_id_sigmasumhf
  integer(hid_t) :: dspace_id_siw,dset_id_siwloc,dset_id_siwdmft
  integer :: rank_siwk,rank_siw
  integer(hsize_t),dimension(:),allocatable :: dims_siwk,dims_siw
  complex(kind=8),dimension(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp) :: sigma_sum,sigma_sum_hf
  complex(kind=8),dimension(ndim,ndim,-iwfmax_small:iwfmax_small-1) :: sigma_loc,sigma_sum_dmft
  complex(kind=8),dimension(ndim,ndim,nkpz,nkpy,nkpx,2*iwfmax_small) :: siwk_outputarray
  complex(kind=8),dimension(ndim,ndim,2*iwfmax_small) :: siw_outputarray

  integer :: i1,i2,ikx,iky,ikz,iw


  call h5open_f(err)


  ! create outputarray (essentially just reshape the original array)
  ! doing this in-place would be much more memory-efficient!
  do i1=1,ndim
    do i2=1,ndim
      do ikz=1,nkpz
        do iky=1,nkpy
          do ikx=1,nkpx
            do iw=1,2*iwfmax_small
              siwk_outputarray(i2,i1,ikz,iky,ikx,iw) = sigma_sum(i1,i2,iw-iwfmax_small-1,(ikz-1)+(iky-1)*nkpz+(ikx-1)*nkpy*nkpz+1)
            end do
          end do
        end do
      end do
    end do
  end do

  rank_siwk=6
  allocate(dims_siwk(rank_siwk))
  dims_siwk=(/ndim,ndim,nkpz,nkpy,nkpx,2*iwfmax_small/)
  call h5screate_f(H5S_SIMPLE_F,dspace_id_sigmasum,err)
  call h5sset_extent_simple_f(dspace_id_sigmasum,rank_siwk,dims_siwk,dims_siwk,err)

  call h5fopen_f(filename_output,H5F_ACC_RDWR_F,file_id,err)
  call h5gopen_f(file_id,'selfenergy',grp_id_siwk,err)
  call h5dcreate_f(grp_id_siwk,'dga',compound_id,dspace_id_sigmasum,dset_id_sigmasum,err)
  call h5dwrite_f(dset_id_sigmasum,type_r_id,real(siwk_outputarray),dims_siwk,err)
  call h5dwrite_f(dset_id_sigmasum,type_i_id,aimag(siwk_outputarray),dims_siwk,err)
  call h5dclose_f(dset_id_sigmasum,err)

  ! in order to save memory, the same array is used again, now for the Hartree-Fock contribution.
  do i1=1,ndim
    do i2=1,ndim
      do ikz=1,nkpz
        do iky=1,nkpy
          do ikx=1,nkpx
            do iw=1,2*iwfmax_small
              siwk_outputarray(i2,i1,ikz,iky,ikx,iw) = sigma_sum_hf(i1,i2,iw-iwfmax_small-1,(ikz-1)+(iky-1)*nkpz+(ikx-1)*nkpy*nkpz+1)
            end do
          end do
        end do
      end do
    end do
  end do

  ! the dimensions are the same, so we can use the same data space
  call h5dcreate_f(grp_id_siwk,'hartree_fock_nonloc',compound_id,dspace_id_sigmasum,dset_id_sigmasumhf,err)
  call h5dwrite_f(dset_id_sigmasumhf,type_r_id,real(siwk_outputarray),dims_siwk,err)
  call h5dwrite_f(dset_id_sigmasumhf,type_i_id,aimag(siwk_outputarray),dims_siwk,err)
  call h5dclose_f(dset_id_sigmasumhf,err)


  deallocate(dims_siwk)

  rank_siw=3
  allocate(dims_siw(rank_siw))
  dims_siw=(/ndim,ndim,2*iwfmax_small/)
  call h5screate_f(H5S_SIMPLE_F,dspace_id_siw,err)
  call h5sset_extent_simple_f(dspace_id_siw,rank_siw,dims_siw,dims_siw,err)

  do i1=1,ndim
    do i2=1,ndim
      do iw=1,2*iwfmax_small
        siw_outputarray(i2,i1,iw) = sigma_loc(i1,i2,iw-iwfmax_small-1)
      end do
    end do
  end do

  call h5dcreate_f(grp_id_siwk,'dga_loc',compound_id,dspace_id_siw,dset_id_siwloc,err)
  call h5dwrite_f(dset_id_siwloc,type_r_id,real(siw_outputarray),dims_siw,err)
  call h5dwrite_f(dset_id_siwloc,type_i_id,aimag(siw_outputarray),dims_siw,err)
  call h5dclose_f(dset_id_siwloc,err)


  do i1=1,ndim
    do i2=1,ndim
      do iw=1,2*iwfmax_small
        siw_outputarray(i2,i1,iw) = sigma_sum_dmft(i1,i2,iw-iwfmax_small-1)
      end do
    end do
  end do

  call h5dcreate_f(grp_id_siwk,'dmft',compound_id,dspace_id_siw,dset_id_siwdmft,err)
  call h5dwrite_f(dset_id_siwdmft,type_r_id,real(siw_outputarray),dims_siw,err)
  call h5dwrite_f(dset_id_siwdmft,type_i_id,aimag(siw_outputarray),dims_siw,err)
  call h5dclose_f(dset_id_siwdmft,err)

  deallocate(dims_siw)

  call h5gclose_f(grp_id_siwk,err)
  call h5fclose_f(file_id,err)
  call h5close_f(err)

end subroutine output_eom_hdf5

end module hdf5_module
