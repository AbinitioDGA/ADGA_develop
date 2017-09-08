module global
  implicit none
  public
  integer :: lines
  character(len=1), parameter :: cmnt = '#'
  character(len=1), parameter :: seperator = '='
  character(len=150), allocatable :: file_temp(:), file_save(:)
end module

module subspaces
  implicit none

  contains

  ! return true if all 4 legs are in the same correlated subspace
  logical function index2cor(nineq,ndims,m,n,o,p)
    implicit none
    integer, intent(in) :: nineq
    integer, intent(in) :: ndims(nineq,2)
    ! band indices from specific beginning or end point of a 1PG
    integer, intent(in) :: m,n,o,p
    ! inequivalent atom number for specific index
    integer :: a,b,c,d
    integer :: dimstart,dimend,ineq, i

    a=0;b=0;c=0;d=0

    do ineq=1,nineq
      dimstart=1
      do i=2,ineq
        dimstart=dimstart+ndims(i-1,1)+ndims(i-1,2)
      enddo
      dimend=dimstart+ndims(ineq,1)-1 ! only in correlated sub space
      if ( m .ge. dimstart .and. m .le. dimend ) a=ineq
      if ( n .ge. dimstart .and. n .le. dimend ) b=ineq
      if ( o .ge. dimstart .and. o .le. dimend ) c=ineq
      if ( p .ge. dimstart .and. p .le. dimend ) d=ineq
    enddo

    ! checking if everything is on the same atom
    ! AND on correlated bands (non correlated lines would have ineq=0)
    if ( (a .eq. b) .and. (c .eq. d) .and. (a .eq. d) .and. (a .ne. 0)) then
      index2cor = .true.
    else
      index2cor = .false.
    endif
  end function index2cor

  ! return true if all 4 legs are in the same uncorrelated subspace
  logical function index2uncor(nineq,ndims,m,n,o,p)
    implicit none
    integer, intent(in) :: nineq
    integer, intent(in) :: ndims(nineq,2)
    ! band indices from specific beginning or end point of a 1PG
    integer, intent(in) :: m,n,o,p
    ! inequivalent atom number for specific index
    integer :: a,b,c,d
    integer :: dimstart,dimend,ineq, i

    a=0;b=0;c=0;d=0

    do ineq=1,nineq
      dimstart=ndims(1,1)+1
      do i=2,ineq
        dimstart=dimstart+ndims(i,1)+ndims(i-1,2)
      enddo
      dimend=dimstart+ndims(ineq,2)-1
      if ( m .ge. dimstart .and. m .le. dimend ) a=ineq
      if ( n .ge. dimstart .and. n .le. dimend ) b=ineq
      if ( o .ge. dimstart .and. o .le. dimend ) c=ineq
      if ( p .ge. dimstart .and. p .le. dimend ) d=ineq
    enddo

    if ( (a .eq. b) .and. (c .eq. d) .and. (a .eq. d) .and. (a .ne. 0)) then
      index2uncor = .true.
    else
      index2uncor = .false.
    endif
  end function index2uncor

  ! returns number of impurity if all 4 legs on the same
  ! otherwise returns 0
  integer function index2ineq(nineq,ndims,m,n,o,p)
    implicit none
    integer, intent(in) :: nineq
    integer, intent(in) :: ndims(nineq,2)
    ! band indices from specific beginning or end point of a 1PG
    integer, intent(in) :: m,n,o,p
    ! inequivalent atom number for specific index
    integer :: a,b,c,d
    integer :: dimstart,dimend,ineq, i

    a=0;b=0;c=0;d=0

    do ineq=1,nineq
      dimstart=1
      do i=2,ineq
        dimstart=dimstart+ndims(i-1,1)+ndims(i-1,2)
      enddo
      dimend=dimstart+ndims(ineq,1)+ndims(ineq,2)-1
      if ( m .ge. dimstart .and. m .le. dimend ) a=ineq
      if ( n .ge. dimstart .and. n .le. dimend ) b=ineq
      if ( o .ge. dimstart .and. o .le. dimend ) c=ineq
      if ( p .ge. dimstart .and. p .le. dimend ) d=ineq
    enddo

    if ( (a .eq. b) .and. (c .eq. d) .and. (a .eq. d) ) then
      index2ineq = a
    else
      index2ineq = 0
    endif
  end function index2ineq

end module subspaces

module lookup
  use global
  implicit none
  private
  integer :: i, pst
  character(len=150) :: str_temp
  public :: string_find, int_find, float_find, group_find, subgroup_find

  contains

  subroutine string_find(search_string, save_string, search_start, search_end)
    character(*), intent(in)  :: search_string
    character(len=150), intent(inout) :: save_string ! keep default string
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        save_string=adjustl(trim(str_temp(pst+1:)))
      endif
    enddo
  end subroutine string_find

  subroutine int_find(search_string, save_int, search_start, search_end)
    character(*), intent(in)  :: search_string
    integer, intent(inout) :: save_int ! keep default values
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=adjustl(trim(str_temp(pst+1:)))
        read(str_temp,'(I5)') save_int
      endif
    enddo
  end subroutine int_find

  subroutine float_find(search_string, save_float, search_start, search_end)
    character(*), intent(in)  :: search_string
    real(8), intent(inout) :: save_float ! keep default values
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=adjustl(trim(str_temp(pst+1:)))
        read(str_temp,'(F13.8)') save_float
      endif
    enddo
  end subroutine float_find

  subroutine group_find(search_string, save_start, save_end)
    character(*), intent(in) :: search_string
    integer, intent(out) :: save_start, save_end
    save_start=0
    save_end=0

    do i=1,lines
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        save_start=i+1
        exit
      endif
    enddo

    do i=save_start, lines
      if (index(trim(file_save(i)),'[') .eq. 1) then
        if (index(trim(file_save(i)),'[[') .eq. 1) then ! skip subgroups
          cycle
        endif
        save_end=i-1 ! one above the next session
        exit
      endif
    enddo

    if(save_end .eq. 0) then
      save_end = lines ! if nothing else is found, until the end of the file
    endif
  end subroutine group_find

  subroutine subgroup_find(search_string, search_start, search_end, save_start, save_end)
    character(*), intent(in) :: search_string
    integer, intent(in) :: search_start, search_end
    integer, intent(out) :: save_start, save_end
    save_start=0
    save_end=0

    do i=search_start, search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        save_start=i+1
        exit
      endif
    enddo

    do i=save_start, search_end
      if (index(trim(file_save(i)),'[') .eq. 1) then
        save_end=i-1 ! one above the next session
        exit
      endif
    enddo

    if(save_end .eq. 0) then
      save_end = lines ! if nothing else is found, until the end of the file
    endif
  end subroutine subgroup_find

end module lookup


program umatrix
  use global
  use subspaces
  use lookup

  implicit none
  integer :: i,j,k,l,ineq,stat
  integer :: int_temp
  character(len=150) :: str_temp, str_ineq

  integer :: nineq, ndim
  integer, allocatable :: ndims(:,:)
  real(8), allocatable :: Umat(:,:,:,:)
  real(8), allocatable :: Udd(:),Vdd(:),Jdd(:)
  real(8), allocatable :: Udp(:),Vdp(:),Jdp(:)
  real(8), allocatable :: Upp(:),Vpp(:),Jpp(:)
  character(len=150), allocatable :: hamiltonian(:)
  integer,allocatable :: mode_ineq(:)

  integer :: search_start, search_end
  integer :: subsearch_start, subsearch_end
  character(len=150) :: config_file, output
  integer :: pst, empty

  ! Config File checks
  if (iargc() .ne. 1) then
    stop 'The program has to be executed with exactly one argument. (Name of config file)'
  end if
  call getarg(1,config_file)

  open(unit=10,file=trim(config_file),action='read',iostat=stat)
    if (stat .ne. 0) then
      write(*,'(A)') 'Input file cannot be opened'
      close(10)
      goto 99
    endif

    ! line counting
    lines=0
    read_count: do
      read(10,'(A)',END=200) str_temp ! read whole line as string, doesnt skip empty lines
      lines=lines+1
    enddo read_count

    200 continue
    rewind 10

    allocate(file_temp(lines))

    ! remove empty lines and comment strings
    empty=0
    read_temp: do i=1,lines
      read(10,'(A)') str_temp
        str_temp = trim(adjustl(str_temp))
        pst=scan(str_temp,cmnt) ! find out possible comment symbol
        if (pst .eq. 1) then ! whole line is commented
          file_temp(i) = ''
        elseif (pst .eq. 0) then ! no comment symbol found
          file_temp(i) = str_temp
        else  ! getting left side of comment
          file_temp(i) = trim(str_temp(1:pst-1))
        endif

        if (len_trim(file_temp(i)) .eq. 0) then ! filter out empty lines
          empty=empty+1
        endif
    enddo read_temp

    ! rewrite everything to a new clean string array
    allocate(file_save(lines-empty))
    j=1
    read_save: do i=1,lines
      if(len_trim(file_temp(i)) .ne. 0) then
        file_save(j)=file_temp(i)
        j=j+1
      endif
    enddo read_save
    deallocate(file_temp)

    lines=lines-empty

  close(unit=10)


  ! FREE FORMAT LOOKUP
  ! defining default values
  output='umatrix.dat'
  nineq=1

  ! search for General stuff + Allocation of values

  call group_find('[General]', search_start, search_end)
  call string_find('Output', output, search_start, search_end)
  call int_find('NAt', nineq, search_start, search_end)

  allocate(hamiltonian(nineq))
  allocate(mode_ineq(nineq))
  allocate(ndims(nineq,2))
  allocate(Udd(nineq),Vdd(nineq),Jdd(nineq))
  allocate(Udp(nineq),Vdp(nineq),Jdp(nineq))
  allocate(Upp(nineq),Vpp(nineq),Jpp(nineq))
  hamiltonian=''
  mode_ineq=0; ndims=0
  Udd=0.d0; Vdd=0.d0; Jdd=0.d0
  Udp=0.d0; Vdp=0.d0; Jdp=0.d0
  Upp=0.d0; Vpp=0.d0; Jpp=0.d0

  call group_find('[Atoms]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    stop 'Group not found'
  endif
  do ineq=1,nineq
    write(str_ineq,'(A2,I1,A2)') '[[',ineq,']]'
    call subgroup_find(str_ineq, search_start, search_end, subsearch_start, subsearch_end)
    if (subsearch_start .eq. 0) then ! group was not found
      stop 'Subgroup not found'
    endif

    call string_find('Hamiltonian',hamiltonian(ineq),subsearch_start,subsearch_end)
    call int_find('Nd',ndims(ineq,1),subsearch_start,subsearch_end)
    call int_find('Np',ndims(ineq,2),subsearch_start,subsearch_end)
    call float_find('Udd',Udd(ineq),subsearch_start,subsearch_end)
    call float_find('Vdd',Vdd(ineq),subsearch_start,subsearch_end)
    call float_find('Jdd',Jdd(ineq),subsearch_start,subsearch_end)
    call float_find('Upp',Upp(ineq),subsearch_start,subsearch_end)
    call float_find('Vpp',Vpp(ineq),subsearch_start,subsearch_end)
    call float_find('Jpp',Jpp(ineq),subsearch_start,subsearch_end)
    call float_find('Udp',Udp(ineq),subsearch_start,subsearch_end)
    call float_find('Vdp',Vdp(ineq),subsearch_start,subsearch_end)
    call float_find('Jdp',Jdp(ineq),subsearch_start,subsearch_end)

    select case (hamiltonian(ineq))
      case ('Kanamori')
        mode_ineq(ineq) = 1
      case default
        mode_ineq(ineq) = 0
    end select
  enddo


  ! creating umatrix file
  ndim=sum(ndims)
  write(*,'(A,I3)') 'Total number of dimensions ndim: ', ndim
  do ineq=1,nineq
    write(*,'(A,I3,A,A,I3,A,I3)') 'Inequivalent atom', ineq, ': ', 'd: ', ndims(ineq,1), '  p: ', ndims(ineq,2)
  enddo
  write(*,'(A,A)') 'Calculating and writing umatrix to: ', trim(output)

  allocate(Umat(ndim,ndim,ndim,ndim))
  open(unit=10,file=trim(output))

  write(10,*) 'Umatrix File for the ADGA code : band,band,band,band,Uvalue'
  Umat=0.d0

  do i=1,ndim
  do j=1,ndim
  do k=1,ndim
  do l=1,ndim

    ineq=index2ineq(nineq,ndims,i,j,k,l) 
    if (ineq .eq. 0) then ! not on the same impurity
      goto 100
    endif
  
    ! DD - VALUES
    if (index2cor(nineq,ndims,i,j,k,l)) then
      if (i .eq. j .and. k .eq. l)  then
        if (i .eq. k) then
          Umat(i,j,k,l) = Udd(ineq) ! onsite density-density -- 1 1 1 1
          goto 100
        else if (mode_ineq(ineq) .eq. 1) then
          Umat(i,j,k,l) = Jdd(ineq) ! kanamori double hopping -- 1 1 2 2
          goto 100
        endif
      endif

      if (i .eq. k .and. j .eq. l) then
        Umat(i,j,k,l) = Vdd(ineq) ! screened density density -- 1 2 1 2
        goto 100
      endif
      if (i .eq. l .and. j .eq. k .and. mode_ineq(ineq) .eq. 1) then
        Umat(i,j,k,l) = Jdd(ineq) ! kanamori spin flip -- 1 2 2 1
        goto 100
      endif

    ! PP - VALUES
    else if(index2uncor(nineq,ndims,i,j,k,l)) then
      if (i .eq. j .and. k .eq. l)  then
        if (i .eq. k) then
          Umat(i,j,k,l) = Upp(ineq) ! onsite density-density -- 1 1 1 1
          goto 100
        else if (mode_ineq(ineq) .eq. 1) then
          Umat(i,j,k,l) = Jpp(ineq) ! kanamori double hopping -- 1 1 2 2
          goto 100
        endif
      endif

      if (i .eq. k .and. j .eq. l) then
        Umat(i,j,k,l) = Vpp(ineq) ! screened density density -- 1 2 1 2
        goto 100
      endif
      if (i .eq. l .and. j .eq. k .and. mode_ineq(ineq) .eq. 1) then
        Umat(i,j,k,l) = Jpp(ineq) ! kanamori spin flip -- 1 2 2 1
        goto 100
      endif

    ! DP - VALUES
    else
      if (i .eq. j .and. k .eq. l)  then
        if (i .eq. k) then
          Umat(i,j,k,l) = Udp(ineq) ! onsite density-density -- 1 1 1 1
          goto 100
        else if (mode_ineq(ineq) .eq. 1) then
          Umat(i,j,k,l) = Jdp(ineq) ! kanamori double hopping -- 1 1 2 2
          goto 100
        endif
      endif

      if (i .eq. k .and. j .eq. l) then
        Umat(i,j,k,l) = Vdp(ineq) ! screened density density -- 1 2 1 2
        goto 100
      endif
      if (i .eq. l .and. j .eq. k .and. mode_ineq(ineq) .eq. 1) then
        Umat(i,j,k,l) = Jdp(ineq) ! kanamori spin flip -- 1 2 2 1
        goto 100
      endif
    endif

  100 write(10,'(4I10,F15.8)') i,j,k,l,Umat(i,j,k,l)

  enddo
  enddo
  enddo
  enddo

  deallocate(file_save,hamiltonian,mode_ineq)
  deallocate(Udd,Vdd,Jdd,Upp,Vpp,Jpp,Udp,Vdp,Jdp)

  99 continue
  write(*,'(A)') 'Done.'

end program
