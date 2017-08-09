!=========================================================================
module indexutils
!==============================================
contains

!==============================================================================
subroutine get_channel(s1, s2, s3, s4, ichannel)
  implicit none

  integer, intent(in) :: s1, s2, s3, s4
  integer, intent(out) :: ichannel

  if ((s1==s4) .and. (s2==s3) .and. (s1 .ne. s2)) then
     ichannel = 1
  else
     ichannel = 2
  endif

end subroutine get_channel
!================================================================================


!================================================================================
subroutine get_orb_sym(b1, b2, b3, b4, nbands, ntot, ind_band_list)
  implicit none
  
  integer, intent(in) :: b1, b2, b3, b4, nbands
  integer,intent(out) :: ntot
  integer,intent(out) :: ind_band_list(:)
  integer :: i, j, icount, ind_band

  if ((b1==b2) .and. (b2==b3) .and. (b3==b4)) then
     ntot = nbands
     
     do i=1,nbands
        call component2index_band(nbands, ind_band, i, i, i, i)
        ind_band_list(i) = ind_band
     enddo

  else
     ntot = nbands**2-nbands
     icount = 1
     do i=1,nbands
        do j=1,nbands
           if (i .ne. j) then
              ind_band = -1
        
              if((b1==b2) .and. (b3==b4) .and. (b1 .ne. b3)) then
                 call component2index_band(nbands, ind_band, i, i, j, j)
                 
              else if((b1==b3) .and. (b2==b4) .and. (b1 .ne. b2)) then
                 call component2index_band(nbands, ind_band, i, j, i, j)

              else if((b1==b4) .and. (b2==b3) .and. (b1 .ne. b2)) then
                 call component2index_band(nbands, ind_band, i, j, j, i)

              endif

              ind_band_list(icount) = ind_band
              icount = icount+1
           endif
        enddo
     enddo

  endif

end subroutine get_orb_sym
!================================================================================


!================================================================================
! converting a band-spin pattern into an index
subroutine component2index(Nbands, ind, b1, s1, b2, s2, b3, s3, b4, s4)
  implicit none

  integer,intent(in) :: Nbands
  integer,intent(in) :: b1, s1, b2, s2, b3, s3, b4, s4
  integer,intent(inout) :: ind
  integer :: g1, g2, g3, g4

  g1=2*(b1-1) + s1
  g2=2*(b2-1) + s2
  g3=2*(b3-1) + s3
  g4=2*(b4-1) + s4

  ind =  8*Nbands**3*(g1-1) + 4*Nbands**2*(g2-1) + 2*Nbands*(g3-1) + g4

end subroutine component2index
!=================================================================================

!================================================================================
! converting an index into a band-spin pattern
subroutine index2component(Nbands, ind, b1, s1, b2, s2, b3, s3, b4, s4)

  implicit none
  integer,intent(in) :: Nbands,ind
  integer,intent(inout) :: b1, s1, b2, s2, b3, s3, b4, s4
  integer :: tmp1,tmp2,tmp3,ind_tmp
  integer :: g1,g2,g3,g4

  ! the proposed back conversion assumes the indices are
  ! given form 0 to max-1  
  ind_tmp = ind - 1
  tmp1 = 8*Nbands**3
  tmp2 = 4*Nbands**2
  tmp3 = 2*Nbands

  g1 = ind_tmp/tmp1 + 1
  g2 = (ind_tmp-tmp1*(g1-1))/tmp2 + 1
  g3 = (ind_tmp-tmp1*(g1-1)-tmp2*(g2-1))/tmp3 + 1
  g4 = (ind_tmp-tmp1*(g1-1)-tmp2*(g2-1)-tmp3*(g3-1)) + 1

  s1=mod(g1-1,2)+1
  b1=(g1-s1)/2+1

  s2=mod(g2-1,2)+1
  b2=(g2-s2)/2+1

  s3=mod(g3-1,2)+1
  b3=(g3-s3)/2+1

  s4=mod(g4-1,2)+1
  b4=(g4-s4)/2+1

end subroutine index2component
!============================================================================

!=============================================================================
!converting a band-pattern into an index
subroutine component2index_band(Nbands, ind, b1, b2, b3, b4)
  implicit none

  integer,intent(in) :: Nbands
  integer,intent(in) :: b1, b2, b3, b4
  integer,intent(out) :: ind

  ind =  Nbands**3*(b1-1) + Nbands**2*(b2-1) + Nbands*(b3-1) + b4

end subroutine component2index_band

end module
!============================================================================







!====================================================================================
program symmetrize_vertex
!====================================
  use indexutils
  use hdf5
  use hdf5_module
  use parameters_module
  implicit none

  character(len=150) :: cmd_arg
  integer(hid_t) :: file_id, new_file_id
  character(len=40) :: grpname, name_buffer, name_buffer_value, name_buffer_error
  integer(hid_t) :: nmembers, imembers, itype, grp_id, g4iw_id, g4err_id
  integer(hid_t) :: dset_dens_id, dset_magn_id, dset_err_id
  integer(hid_t) :: dspace_iwb_id, dspace_iwf_id
  integer(hsize_t), dimension(1) :: dim_iwb, dim_iwf
  integer(hsize_t), dimension(3) :: g4iw_dims, g4iw_maxdims, g4err_dims
  integer, parameter :: rank = 2, rank_iw = 1
  
  double precision, allocatable :: g4iw_r(:,:,:), g4iw_i(:,:,:), g4err(:,:,:), diff_r(:,:), diff_i(:,:)
  integer :: ind, b1, s1, b2, s2, b3, s3, b4, s4
  integer, allocatable :: Nbands(:)
  integer :: iwb, iwf, iwf1, iwf2 
  logical :: su2_only
  logical, allocatable :: create_comp(:,:)
  double precision, allocatable :: iwb_array(:), iwf_array(:)
  integer :: ichannel, ntot, ind_band, icount
  integer, allocatable :: ind_band_list(:)
  character(len=150), allocatable :: filename_vertex_ineq(:)
  character(len=1) :: arg_sym

  real(kind=8) :: start, finish
!================================================================

! alternatively - read input
! this is definitely more clearer
  write(*,'(A)',advance='no') 'Number of inequivalent atoms: '
  read(*,*) nineq
  allocate(filename_vertex_ineq(nineq))
  allocate(Nbands(nineq))
  do ineq=1,nineq
    write(*,'(A,I1,A)',advance='no') 'Vertex of inequivalent atom ', ineq, ': '
    read(*,*) filename_vertex_ineq(ineq)
    write(*,'(A,I1,A)',advance='no') 'Number of correlated bands for inequivalent atom ', ineq, ': '
    read(*,*) Nbands(ineq)
  enddo
  write(*,'(A)',advance='no') 'Outputfile for symmetrized Vertex: '
  read(*,*) filename_vertex_sym
  write(*,*)
  write(*,'(A)',advance='no') 'SU2 symmetry only (s) or SU2 AND orbital symmetry (o)?: '
  read(*,*) arg_sym
  if (arg_sym .eq. 'o' .or. arg_sym .eq. 'O') then
    su2_only = .false.
  elseif (arg_sym .eq. 's' .or. arg_sym .eq. 'S') then
    su2_only = .true.
  else
    su2_only = .false.
    write(*,*) 'Wrong input - Using only SU2 symmetry.'
  endif


! read command line arguments -> number ineq atoms, input filenames and output filename, number of bands.
  ! call getarg(1,cmd_arg)
  ! read(cmd_arg,'(I1)') nineq

  ! if (.not. iargc() .eq. 2*nineq+2 ) then
  !   write(*,*) 'The program has to be executed with the following arguments'
  !   write(*,*) 'number of inequivalent atoms, names of input files (as many as inequivalent atoms)'
  !   write(*,*) 'name of output file, number of correlated (d) bands for each input file'
  !   stop
  ! endif

  ! allocate(filename_vertex_ineq(nineq))
  ! allocate(Nbands(nineq))
  
  ! do ineq=1,nineq
  !   call getarg(1+ineq,cmd_arg)
  !   filename_vertex_ineq(ineq)=trim(cmd_arg)
  ! enddo
  ! call getarg(2+nineq,cmd_arg)
  ! filename_vertex_sym = trim(cmd_arg)
  ! do ineq=1,nineq
  !   call getarg(2+nineq+ineq,cmd_arg)
  !   read(cmd_arg,'(I1)') Nbands(ineq)
  ! enddo

!================================================================
!Define orbital symmetry here:
  su2_only = .false. 
  write(*,*) 'Symmetrizing ',(filename_vertex_ineq(ineq),ineq=1,nineq),'>>>>>',filename_vertex_sym
  write(*,*) 'Total number of bands: ',sum(Nbands)
  if(su2_only) then
    write(*,*) 'Using only SU2 symmetry'
  else
    write(*,*) 'Using orbital and SU2 symmetry'
  endif
    
!=================================================================

  call cpu_time(start)
  
  call h5open_f(err)
 
  call create_complex_datatype
 


! loop over number of inequivalent atoms
  do ineq=1,nineq


! open vertex file 'vertex_full.hdf5':
    call h5fopen_f(trim(filename_vertex_ineq(ineq)), h5f_acc_rdonly_f, file_id, err)

    if (ineq .eq. 1) then ! only write the axes once
! create new file for the symmetrised vertex:
      call h5fcreate_f(filename_vertex_sym, h5f_acc_trunc_f, new_file_id, err)
! get fermionic and bosonic Matsubara axes: 
      call read_axes(file_id, iwb_array, iwf_array, dspace_iwb_id, dspace_iwf_id, dim_iwb, dim_iwf)
      call write_axes(new_file_id, iwb_array, iwf_array, dspace_iwb_id, dspace_iwf_id, dim_iwb, dim_iwf)
    endif

! write magn/iwb and dens/iwb in the output file vertex_sym.hdf5: 
    call create_channels(new_file_id, ineq)

!=========================================================================

! allocate quantities that are needed in the loop afterwards: 
    g4iw_dims = (/2*iwfmax,2*iwfmax,2*iwbmax+1/)

    allocate(g4iw_r(2*iwfmax,2*iwfmax,2*iwbmax+1)) 
    allocate(g4iw_i(2*iwfmax,2*iwfmax,2*iwbmax+1))
    allocate(g4err(2*iwfmax, 2*iwfmax, 2*iwbmax+1))

    allocate(tmp_r_1(g4iw_dims(1), g4iw_dims(2)), tmp_i_1(g4iw_dims(1), g4iw_dims(2)), tmp_err_1(g4iw_dims(1),g4iw_dims(2)))
    allocate( diff_r(g4iw_dims(1), g4iw_dims(2)), diff_i(g4iw_dims(1), g4iw_dims(2)))

    allocate(create_comp(2,Nbands(ineq)**4))
    create_comp = .true.

    allocate(ind_band_list(Nbands(ineq)**2))

! create dataspace:
    dims = (/g4iw_dims(1), g4iw_dims(2)/)
    call h5screate_simple_f(rank, dims, dspace_id, err)
 
!============================================================================
! iterate over all groups(band-spin combinations) in the old vertex file:
    call h5gn_members_f(file_id, "/", nmembers, err)

    do imembers = 1,nmembers - 1 
       write(*,*) imembers
       call h5gget_obj_info_idx_f(file_id, "/", imembers, name_buffer, itype, err)
       read(name_buffer,'(I5.5)') ind
     
       ! read the current group in the old vertex file:
       call h5gopen_f(file_id,name_buffer,grp_id,err)

       write(name_buffer_value, '((I5.5),A6)'), ind, "/value"
       call h5dopen_f(file_id, name_buffer_value, g4iw_id, err)

       write(name_buffer_error, '((I5.5),A6)'), ind, "/error"
       call h5dopen_f(file_id, name_buffer_error, g4err_id, err)

       call h5dread_f(g4iw_id, type_r_id, g4iw_r, g4iw_dims, err)
       call h5dread_f(g4iw_id, type_i_id, g4iw_i, g4iw_dims, err)
       call h5dread_f(g4err_id, h5t_native_double, g4err, g4err_dims, err)

       call h5dclose_f(g4iw_id, err)
       call h5dclose_f(g4err_id, err)
       call h5gclose_f(grp_id, err)
      
  ! get band and spin indices:
       call index2component(Nbands(ineq), ind, b1, s1, b2, s2, b3, s3, b4, s4)
      
  ! get the channel (magnetic: ichannel=1, density: ichannel=2)
       call get_channel(s1, s2, s3, s4, ichannel)
      
       ! get list of ind_band (components to be written in vertex_sym.hdf5)
       if (su2_only) then

          ! without orbital symmetry (only su2):
          call component2index_band(Nbands(ineq), ind_band, b1, b2, b3, b4)
          ntot = 1 
          ind_band_list(1) = ind_band

       else

          !with orbital symmetry: 
          call get_orb_sym(b1, b2, b3, b4, Nbands(ineq), ntot, ind_band_list)

       endif     

       !divisions needed for the channels:
       g4iw_r = g4iw_r/(2.d0*ntot)
       g4iw_i = g4iw_i/(2.d0*ntot)

       if (ichannel == 1) then  !is this correct??? 
          g4err = g4err/(2.d0*ntot)
       else
          g4err = g4err/(4.d0*ntot)
       endif
      
       ! iterates over all components that need to be written: 
       do icount = 1, ntot

         
          if (create_comp(ichannel, ind_band_list(icount))) then
      
             ! creates new groups and datasets if they are not there yet:
             do iwb = 0, 2*iwbmax
                call create_component(new_file_id, ichannel, iwb, ind_band_list(icount), ineq)
             enddo
             
             create_comp(ichannel, ind_band_list(icount)) = .false.

          endif
          
          !adds the data:
          do iwb = 0, 2*iwbmax

             call add_to_component(new_file_id, ichannel, iwb, ind_band_list(icount), g4iw_r, g4iw_i, g4err, ineq)

          enddo
         
       enddo !icount
   
    enddo ! nmembers
  
 
    call h5sclose_f(dspace_id, err)

    deallocate(g4iw_r, g4iw_i, g4err)
    deallocate(tmp_r_1, tmp_i_1, tmp_err_1, diff_r, diff_i) 
    deallocate(ind_band_list)
    deallocate(create_comp)

    call h5fclose_f(file_id, err)

  enddo

  call h5fclose_f(new_file_id, err)

  call cpu_time(finish)
  write(*,*)'cpu-time =', finish-start
  

end program symmetrize_vertex


!su(2) symmetry check to be added somewhere:

!open(21, file="su2_error.dat", status='unknown')
!write(21,*)'# b1, b2, b3, b4, iwb, iwf1, iwf2'

!diff_r(:,:) = abs(2.d0*g4iw_r(:,:,iwb+1)-2.d0*tmp_r(:,:))
!diff_i(:,:) = abs(2.d0*g4iw_i(:,:,iwb+1)-2.d0*tmp_i(:,:))

!do iwf1=1,2*iwfmax
!   do iwf2=1,2*iwfmax
!      if((diff_r(iwf1,iwf2) .gt. 1.d0*g4err(iwf1,iwf2,iwb+1)) .or. (diff_i(iwf1,iwf2) .gt. 1.d0*g4err(iwf1,iwf2,iwb+1))) then
!         write(21,*) ind, b1, b2, b3, b4, iwb+1, iwf1, iwf2
!      endif
!   enddo
!enddo
!close(21)
