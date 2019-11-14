module cond_module
  use parameters_module
  use hdf5
  use hdf5_wrapper
  implicit none


  contains

  subroutine read_hkder(er, erstr)
    implicit none

    integer, intent(out) :: er
    character(len=*), intent(out) :: erstr
    integer(hid_t) :: ifile
    integer, allocatable :: hdf_shape(:)
    integer :: hkdershape(3)
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

end module cond_module
