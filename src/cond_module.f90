module cond_module
  use parameters_module
  use one_particle_quant_module, only: get_gkiw
  use hdf5
  use hdf5_wrapper
  use mpi_org, only: master, mpi_wrank
  implicit none


  contains

  subroutine read_hkder(er, erstr)
    implicit none

    integer, intent(out) :: er
    character(len=*), intent(out) :: erstr
    integer(hid_t) :: ifile
    integer, allocatable :: hdf_shape(:)
    integer :: itest

    ! reset error code
    er = 0
    erstr = ''

    call hdf5_open_file(filename_hkder, ifile, rdonly=.true.)
    call hdf5_get_shape(ifile, 'hkder', hdf_shape)

    if (size(hdf_shape) /= 3) then
      er = 1
      erstr = 'Hamiltonian derivative shape from file has not appropriate shape'
      return
    endif

    if (hdf_shape(1) /= 3 .or. hdf_shape(2) /= ndim .or. hdf_shape(3) /= nkp) then
      er = 2
      erstr = 'Hamiltonian derivative array has not appropriate size'
      return
    endif

    ! if no errors we load the array
    call hdf5_read_data(ifile, 'hkder', hkder)
    call hdf5_close_file(ifile)

  end subroutine read_hkder

  subroutine read_gkiw_cond(er, erstr)
    implicit none

    integer, intent(out) :: er
    character(len=*), intent(out) :: erstr
    integer(hid_t) :: ifile
    integer, allocatable :: hdf_shape(:)
    complex(kind=8), allocatable :: tempc4(:,:,:,:) ! temporary array to load the array
    integer :: gkiwshape(4)
    integer :: itest
    integer :: iwgkiwmax, iwgkiwextmax
    integer :: iwf

    ! reset error code
    er = 0
    erstr = ''

    call hdf5_open_file(filename_condlegs, ifile, rdonly=.true.)
    call hdf5_get_shape(ifile, 'giwk', hdf_shape)

    if (size(hdf_shape) /= 4) then
      er = 1
      erstr = 'Greens function from file has not appropriate shape'
      return
    endif

    if (hdf_shape(1) /= ndim .or. hdf_shape(2) /= ndim .or. hdf_shape(3) /= nkp &
        .or. mod(hdf_shape(4), 2) /= 0) then
      er = 2
      erstr = 'Greens function array has not appropriate size'
      return
    endif

    ! if no errors we load the array
    ! gkiwfull
    call hdf5_read_data(ifile, 'giwk', tempc4)
    gkiwshape = shape(tempc4)
    iwgkiwmax = gkiwshape(4) / 2 !
    allocate(gkiwfull(ndim,ndim,nkp,-iwgkiwmax:iwgkiwmax-1))
    gkiwfull = tempc4 ! all this workaround to have the proper offset inbuilt

    if (mpi_wrank .eq. master) then
      if (extend_cond_bubble) then
        deallocate(tempc4)
  ! gkiwfullbubble
        call hdf5_read_data(ifile, 'giwkext', tempc4)
        gkiwshape = shape(tempc4)
        iwgkiwextmax = gkiwshape(4) / 2
        allocate(gkiwfullbubble(ndim,ndim,nkp,-iwgkiwextmax:iwgkiwextmax-1))
        gkiwfullbubble = tempc4
        ! OVERWRITE THE start and stop index for the conductivity bubble
        iwcstart = -iwgkiwextmax+iwbcond
        iwcstop  = iwgkiwextmax-iwbcond-1
      else
        allocate(gkiwfullbubble(ndim,ndim,nkp,iwcstart-iwbcond:iwcstop+iwbcond))
        do iwf = -iwcstart-iwbcond, iwcstop+iwbcond
          gkiwfullbubble(:,:,:,iwf) = gkiwfull(:,:,:,iwf)
        enddo
      endif
    endif

    deallocate(tempc4)
    call hdf5_close_file(ifile)

    if (iwgkiwmax .lt. (iwfmax_small + iwbmax_small)) then
      er = 3
      if (ounit .gt. 0) then
        write(ounit,*) 'Error: N1iwf from external Greens function must be greater than N4iwf + N4iwb'
        write(ounit,*) 'N1iwf=',iwgkiwmax,'  N4iwf=',iwfmax_small,'  N4iwb=',iwbmax_small
      endif
      return
    endif

  end subroutine read_gkiw_cond

  subroutine create_gkiw_cond()
    implicit none

    integer :: iwf,ik
    integer :: iwfmax_loc
    complex(kind=8), allocatable :: gkiw(:,:)

    iwfmax_loc = iwfcond+iwbmax_small+iwbcond
    allocate(gkiw(ndim,ndim))

    ! fill the 'small' gkiw array
    do iwf = -iwfmax_loc,iwfmax_loc-1
      do ik=1,nkp
        call get_gkiw(ik,iwf,0,gkiw) ! g(ik,iwf)
        gkiwfull(:,:,ik,iwf) = gkiw
      enddo
    enddo

    ! fill the 'large gkiw array'
    if (mpi_wrank .eq. master) then
      do iwf = iwcstart-iwbcond, iwcstop+iwbcond
        do ik=1,nkp
          call get_gkiw(ik,iwf,0,gkiw) ! g(ik,iwf)
          gkiwfullbubble(:,:,ik,iwf) = gkiw
        enddo
      enddo
    endif

    deallocate(gkiw)

  end subroutine create_gkiw_cond

end module cond_module
