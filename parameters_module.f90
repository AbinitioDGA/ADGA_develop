module parameters_module
  implicit none
  public

  complex(kind=8) ci
  parameter (ci=(0.d0,1.d0))
  integer :: nkp, ndim, ndims, ndim2, maxdim
  integer :: iwmax, iwbmax, iwfmax, iwbmax_small, iwfmax_small,nk_frac
  double precision :: mu
  integer :: nqp
  character(len=100) :: filename,filename_vertex,filename_umatrix
  logical :: orb_sym,small_freq_box

contains
subroutine read_config()
  implicit none
  character(len=100) :: cmd_arg
  character(len=100) :: config_file
  character(len=150) :: str_tmp
  integer :: int_tmp_1,int_tmp_2

  call getarg(1,cmd_arg)
  config_file=trim(cmd_arg)
  write(*,*) 'Reading config.',config_file

  open(unit=1,file=config_file)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename=trim(str_tmp)
!  write(*,*) 'filename one-particle quantities: ',filename
  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_vertex=trim(str_tmp)
!  write(*,*) 'filename two-particle quantities: ',filename_vertex
  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_umatrix=trim(str_tmp)
  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2
  orb_sym=int_tmp_1
  small_freq_box=int_tmp_2
  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2
  iwfmax_small=int_tmp_1
  iwbmax_small=int_tmp_2
!  write(*,*) iwfmax_small
  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1
  nk_frac = int_tmp_1
  close(1)
end subroutine read_config

  

end module parameters_module
