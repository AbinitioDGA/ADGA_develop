module aux
  implicit none

  contains

  subroutine component2index_band(Nbands, ind, b1, b2, b3, b4)
    implicit none
    integer,intent(in) :: Nbands
    integer,intent(in) :: b1, b2, b3, b4
    integer,intent(out) :: ind

    ind =  Nbands**3*(b1-1) + Nbands**2*(b2-1) + Nbands*(b3-1) + b4
  end subroutine component2index_band



  ! converting an index into a band pattern
  subroutine index2component_band(Nbands, ind, b1, b2, b3, b4)
    implicit none
    integer,intent(in) :: Nbands,ind
    integer,intent(out) :: b1, b2, b3, b4
    integer :: tmp1,tmp2,tmp3,ind_tmp
    integer :: g1,g2,g3,g4

    ! the proposed back conversion assumes the indices are
    ! given form 0 to max-1
    ind_tmp = ind - 1
    tmp1 = Nbands**3
    tmp2 = Nbands**2
    tmp3 = Nbands

    b1 = ind_tmp/tmp1 + 1
    b2 = (ind_tmp-tmp1*(b1-1))/tmp2 + 1
    b3 = (ind_tmp-tmp1*(b1-1)-tmp2*(b2-1))/tmp3 + 1
    b4 = (ind_tmp-tmp1*(b1-1)-tmp2*(b2-1)-tmp3*(b3-1)) + 1
  end subroutine index2component_band



  ! function which checks blockdiagonality of 1PG functions
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
    if ( (a .eqv. b) .and. (c .eqv. d) .and. (a .eqv. d) .and. (a .neqv. 0)) then
      index2ineq = .true.
    else
      index2ineq = .false.
    endif
  end function index2ineq

  subroutine inverse_matrix_3(A)
    implicit none
    complex(kind=8) :: A(3,3),B(3,3),detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                 - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                 + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

    A=B
  end subroutine inverse_matrix_3

  subroutine inverse_matrix_2(A)
    implicit none
    complex(kind=8) :: A(2,2),B(2,2),detinv

   ! Calculate the inverse determinant of the matrix
    detinv = 1.d0/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)

    A=B
  end subroutine inverse_matrix_2

end module
