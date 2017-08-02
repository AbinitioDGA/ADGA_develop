program umatrix
  ! program to create ADGA specific Umatrix File from Values of U and J
  ! if inequivalent atoms are existent: block diagonal
  implicit none
  integer :: i,j,k,l
  integer, parameter :: ndim=6
  real(8) :: U,JJ
  real(8), dimension(ndim,ndim,ndim,ndim) :: Umat

  write(*,'(A)',advance='no') 'Value of U: '
  read(*,*) U
  write(*,'(A)',advance='no') 'Value of J: '
  read(*,*) JJ

  open(unit=10,file='umatrix.dat')
  write(10,*) 'Umatrix File for the ADGA code : band,band,band,band,Uvalue'

  Umat=0.d0

  do i=1,ndim
  do j=1,ndim
  do k=1,ndim
  do l=1,ndim

    if (i .eq. j .and. k .eq. l)  then
      if (i .eq. k) then
        Umat(i,j,k,l) = U ! onsite density-density
        goto 100
      else
        Umat(i,j,k,l) = JJ ! double hopping
        goto 100
      endif
    endif

    if (i .eq. k .and. j .eq. l) then
      Umat(i,j,k,l) = U
      goto 100
    endif
    if (i .eq. l .and. j .eq. k) then
      Umat(i,j,k,l) = JJ
      goto 100
    endif

    !enforce block diagonality and write
100 if (abs(i-j) .ge. 3 .or. abs(i-k) .ge. 3 .or. abs(i-l) .ge. 3 .or. &
        abs(j-k) .ge. 3 .or. abs(j-l) .ge. 3 .or. abs(k-l) .ge. 3) then
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
