program umatrix
  implicit none
  integer :: i,j,k,l
  integer :: ineq, nineq, ndim
  integer :: lines, stat
  integer, allocatable :: ndims(:,:)
  real(8), allocatable :: Umat(:,:,:,:)
  real(8), allocatable :: Udd(:),Vdd(:),Jdd(:)
  real(8), allocatable :: Udp(:),Vdp(:),Jdp(:)
  real(8), allocatable :: Upp(:),Vpp(:),Jpp(:)
  character(len=150), allocatable :: hamiltonian(:)
  integer,allocatable :: mode_ineq(:)

  character(len=150), allocatable :: file_temp(:), file_save(:)
  integer, allocatable :: pos_ineq_start(:), pos_ineq_end(:)
  character(len=150) :: str_temp, str_ineq
  character(len=150) :: config_file
  character(len=150) :: output
  integer :: int_temp
  integer :: index2ineq
  logical :: index2cor, index2uncor
  character(len=1) :: cmnt, seperator
  integer :: pst, empty

  ! Config File checks
  if (iargc() .ne. 1) then
    stop 'The program has to be executed with exactly one argument. (Name of config file)'
  end if
  call getarg(1,config_file)

  open(unit=10,file=trim(config_file),action='read',iostat=stat)
    if (stat .ne. 0) then
      write(*,*) 'input file cannot be opened'
      goto 99
    endif

    lines=0
    read_count: do
      read(10,'(A)',END=200) str_temp ! read whole line as string, doesnt skip empty lines
      lines=lines+1
    enddo read_count

    200 continue
    rewind 10

    allocate(file_temp(lines))

    cmnt = '#'

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

    ! file save now contains all the information without comments 
    ! or break lines

    ! free format detection

    seperator='='
    
    ! search for General stuff + Allocation of values
    do i=1,lines
      if (index(trim(file_save(i)),'Output') .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        output=adjustl(trim(str_temp(pst+1:)))
      endif

      if (index(trim(file_save(i)),'NAt') .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=adjustl(trim(str_temp(pst+1:)))
        read(str_temp,'(I1)') nineq
      endif
    enddo

    allocate(hamiltonian(nineq))
    allocate(mode_ineq(nineq))
    allocate(pos_ineq_start(nineq),pos_ineq_end(nineq))
    allocate(ndims(nineq,2))
    allocate(Udd(nineq),Vdd(nineq),Jdd(nineq))
    allocate(Udp(nineq),Vdp(nineq),Jdp(nineq))
    allocate(Upp(nineq),Vpp(nineq),Jpp(nineq))
    hamiltonian=''
    mode_ineq=0
    ndims=0
    Udd=0.d0; Vdd=0.d0; Jdd=0.d0
    Udp=0.d0; Vdp=0.d0; Jdp=0.d0
    Upp=0.d0; Vpp=0.d0; Jpp=0.d0
    pos_ineq_start=0; pos_ineq_end=0

    ! search for [[i]] positions
    do ineq=1,nineq
      write(str_ineq,'(A2,I1,A2)') '[[',ineq,']]'
      do i=1,lines
        if (index(trim(file_save(i)),trim(str_ineq)) .ne. 0) then
          pos_ineq_start(ineq)=i+1 ! one below that sign
          exit
        endif
      enddo
    enddo

    ! search for areas of ineq definitions
    do ineq=1,nineq
      do i=pos_ineq_start(ineq)+1, lines
        if (index(trim(file_save(i)),'[') .eq. 1) then
          pos_ineq_end(ineq)=i-1 ! one above the next session
          exit
        endif
      enddo
      if(pos_ineq_end(ineq) .eq. 0) then
        pos_ineq_end(ineq) = lines ! if nothing is found, until the end of the file
      endif
    enddo

    
    do ineq=1,nineq
      do i=pos_ineq_start(ineq), pos_ineq_end(ineq)
        if (index(trim(file_save(i)),'Hamiltonian') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          hamiltonian(ineq)=adjustl(trim(str_temp(pst+1:)))
          cycle
        endif
        if (index(trim(file_save(i)),'Nd') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(I1)') ndims(ineq,1)
          cycle
        endif
        if (index(trim(file_save(i)),'Np') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(I1)') ndims(ineq,2)
          cycle
        endif
        if (index(trim(file_save(i)),'Udd') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(F13.8)') Udd(ineq)
          cycle
        endif
        if (index(trim(file_save(i)),'Vdd') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(F13.8)') Vdd(ineq)
          cycle
        endif
        if (index(trim(file_save(i)),'Jdd') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(F13.8)') Jdd(ineq)
          cycle
        endif
        if (index(trim(file_save(i)),'Upp') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(F13.8)') Upp(ineq)
          cycle
        endif
        if (index(trim(file_save(i)),'Vpp') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(F13.8)') Vpp(ineq)
          cycle
        endif
        if (index(trim(file_save(i)),'Jpp') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(F13.8)') Jpp(ineq)
          cycle
        endif
        if (index(trim(file_save(i)),'Udp') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(F13.8)') Udp(ineq)
          cycle
        endif
        if (index(trim(file_save(i)),'Vdp') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(F13.8)') Vdp(ineq)
          cycle
        endif
        if (index(trim(file_save(i)),'Jdp') .ne. 0) then
          str_temp=file_save(i)
          pst=scan(str_temp,seperator)
          str_temp=adjustl(trim(str_temp(pst+1:)))
          read(str_temp,'(F13.8)') Jdp(ineq)
          cycle
        endif
      enddo
    enddo


    do ineq=1,nineq
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

  deallocate(file_save,hamiltonian,mode_ineq,pos_ineq_end,pos_ineq_start)
  deallocate(Udd,Vdd,Jdd,Upp,Vpp,Jpp,Udp,Vdp,Jdp)
  99 continue
  close(10)
  write(*,'(A)') 'Done.'


end program



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
