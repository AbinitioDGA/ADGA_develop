program umatrix
  ! program to create ADGA specific Umatrix File from Values of U and J
  ! considers interaction to be kanamori
  ! if inequivalent atoms are existent: block diagonal
  ! currently made for SrVO3 - 2 layer -- change as you see fit
  implicit none
  integer :: i,j,k,l
  integer :: ineq, idp
  integer :: nineq
  integer :: ndim
  integer, allocatable :: ndims(:,:)
  real(8), allocatable :: Umat(:,:,:,:)
  real(8) :: U,JJ
  logical :: index2ineq
  character(len=150) :: output
  character(len=1) :: mode_input
  integer :: mode

  write(*,*) 'Program to create umatrix file for the ADGA program'
  write(*,'(A)',advance='no') 'Kanamori (k) or Density-Density only (d): '
  read(*,*) mode_input
  if (mode_input .eq. 'k' .or. mode_input .eq. 'K') then
    mode=1
    write(*,*) 'Using Kanamori interaction'
  else
    mode=0
    write(*,*) 'Using Density-Density interaction'
  endif
  write(*,*)
  write(*,'(A)',advance='no') 'Outputfile: '
  read(*,*) output
  write(*,*)
  write(*,'(A)',advance='no') 'Number of inequivalent atoms: '
  read(*,*) nineq
  allocate(ndims(nineq,2))
  write(*,'(A)',advance='no') 'Band configuration (format= d p d p ...): '
  read(*,*) ((ndims(ineq, idp), idp=1,2), ineq=1,nineq)
  ndim=sum(ndims)
  allocate(Umat(ndim,ndim,ndim,ndim))
  write(*,*)
  write(*,'(A)',advance='no') 'Value of U: '
  read(*,*) U
  write(*,'(A)',advance='no') 'Value of J: '
  read(*,*) JJ

  open(unit=10,file=trim(output))
  write(10,*) 'Umatrix File for the ADGA code : band,band,band,band,Uvalue'

  Umat=0.d0

  do i=1,ndim
  do j=1,ndim
  do k=1,ndim
  do l=1,ndim

    if (i .eq. j .and. k .eq. l)  then
      if (i .eq. k) then
        Umat(i,j,k,l) = U ! onsite density-density -- 1 1 1 1
        goto 100
      else if (mode .eq. 1) then
        Umat(i,j,k,l) = JJ ! kanamori double hopping -- 1 1 2 2
        goto 100
      endif
    endif

    if (i .eq. k .and. j .eq. l) then
      Umat(i,j,k,l) = U-2*JJ  ! screened density density -- 1 2 1 2
      ! factor 2 because of cubic unit cell
      goto 100
    endif
    if (i .eq. l .and. j .eq. k .and. mode .eq. 1) then
      Umat(i,j,k,l) = JJ ! kanamori spin flip -- 1 2 2 1
      goto 100
    endif

    !enforce block diagonality and write
100 if (.not. index2ineq(nineq,ndims,i,j,k,l)) then
        Umat(i,j,k,l) = 0.d0
    endif

    write(10,'(4I10,F15.8)') i,j,k,l,Umat(i,j,k,l)

  enddo
  enddo
  enddo
  enddo




  close(unit=10)
  write(*,*) 'Done.'

end program

logical function index2ineq(nineq,ndims,m,n,o,p)
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
    index2ineq = .true.
  else
    index2ineq = .false.
  endif
end function index2ineq
