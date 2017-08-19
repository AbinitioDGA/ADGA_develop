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


 subroutine get_freq_range()
     use parameters_module
     implicit none
     integer(hid_t) :: file_id
     integer(hsize_t), dimension(1) :: iw_dims,iw_maxdims,dim_iwb,dim_iwf,dim_iwf_max,dim_iwb_max
     integer(hid_t) :: iw_id,iwb_id,iwf_id,iw_space_id,dspace_iwb_id,dspace_iwf_id
 
     ! read frequency range of one-particle quantities
     write(*,*) 'one particle quantities in ',filename
     call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, err)
     call h5dopen_f(file_id, ".axes/iw", iw_id, err)
     call h5dget_space_f(iw_id, iw_space_id, err)
     call h5sget_simple_extent_dims_f(iw_space_id, iw_dims, iw_maxdims, err)
     iwmax = iw_dims(1)/2
     call h5dclose_f(iw_id, err)
     call h5fclose_f(file_id,err)

     write(*,*) 'two particle quantities in ',filename_vertex_sym
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
     use parameters_module
     implicit none
     double precision, allocatable :: hr(:,:), hi(:,:)
     integer :: ik,i,j
     double precision :: kx, ky, kz
     
     open(21, file=filename_hk, status='unknown') ! the filename_hk is taken from parameters_module
     read(21,*) nkp,ndim
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
     use parameters_module
     implicit none
     integer :: i,ineq,iband,dimstart, dimend
     integer(hid_t) :: file_id,siw_id, siw_space_id
     integer(hsize_t), dimension(3) :: siw_dims, siw_maxdims
     double precision, allocatable :: siw_data(:,:,:,:)
     double precision :: siw_r, siw_i
     character(len=200) :: name_buffer

     allocate(siw(-iwmax:iwmax-1,ndim))
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
   !  open(34, file=trim(output_dir)//"siw.dat", status='unknown')
   !  do iw=-iwmax,iwmax-1
   !     write(34,'(100F12.6)')iw_data(iw), (real(siw(iw,i)),aimag(siw(iw,i)), i=1,ndim)
   !  enddo
   !  close(34)
   

 end subroutine read_siw

 subroutine read_giw()
     use parameters_module
     implicit none
     integer :: i,ineq,iband,dimstart, dimend
     integer(hid_t) :: file_id,giw_id, giw_space_id
     integer(hsize_t), dimension(3) :: giw_dims, giw_maxdims
     double precision, allocatable :: giw_data(:,:,:,:)
     double precision :: giw_r, giw_i
     character(len=200) :: name_buffer


! read in all inequivalent atoms
     allocate(giw(-iwmax:iwmax-1,ndim))
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
!  open(54, file=trim(output_dir)//"giw.dat", status='unknown')
!  do iw=-iwmax,iwmax-1
!     write(54,'(100F12.6)')iw_data(iw), (real(giw(iw,i)),aimag(giw(iw,i)),i=1,ndim)
!  enddo
!  close(54)
 end subroutine read_giw

 subroutine read_hk_w2dyn()
     use parameters_module
     implicit none
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
    ndim = hk_dims(1)
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
     use parameters_module
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
     use parameters_module
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
     use parameters_module
     implicit none
     integer :: ineq,i,iband,dimstart,dimend
     integer(hid_t) :: file_id,dc_id, dc_space_id
     integer(hsize_t), dimension(2) :: dc_dims, dc_maxdims
     double precision, allocatable :: dc_data(:,:)
     character(len=200) :: name_buffer

     allocate(dc(2,ndim)) ! indices: spin band
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
     use parameters_module
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
                               chi_loc_dens_full(i1,i2) = chi_loc_dens_full(i1,i2)-2.d0*beta*giw(iwf1,i)*giw(iwf2,l)
                            endif

                          else if (vertex_type .eq. connected_g4) then
                            !G_conn:
                            !bubble term -G(\nu)G(\nu-\omega) is added in both channels
                            if((iwf2 .eq. iwf1) .and. i==l .and. j==k)then
                               chi_loc_dens_full(i1,i2) = chi_loc_dens_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j)
                               chi_loc_magn_full(i1,i2) = chi_loc_magn_full(i1,i2)-beta*giw(iwf1,i)*giw(iwf2-iwb,j)
                            endif
                          end if

                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo

        deallocate(g4iw_magn, g4iw_dens, tmp_r, tmp_i)

 end subroutine read_vertex
end module hdf5_module
