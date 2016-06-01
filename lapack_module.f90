module lapack_module
  implicit none
  private
  public inverse_matrix!,inverse_matrix_z
  
!
!   Calculation of the inverse matrix
!
  interface inverse_matrix
    module procedure inverse_matrix_d,inverse_matrix_z
  end interface inverse_matrix
  
  contains
  
  subroutine inverse_matrix_z( a )
    implicit none
    double complex, intent (inout) :: a(:,:)
    integer                     :: dim,info,lwork
    integer,        allocatable :: ipiv(:)
    double complex, allocatable :: work(:)
  
    dim = size( a,1 )
    lwork = 6*dim
!      write(*,*)'dim',dim
    allocate( ipiv(dim),work(lwork) )

    call zgetrf( dim,dim,a,dim,ipiv,info )
    if( info /= 0 )then
      write(6,"(//,'INVERSE_MATRIX: error in ZGETRF, info =',1x,i3)") info
      stop
    end if
    call zgetri( dim,a,dim,ipiv,work,lwork,info )
    if( info /= 0 )then
      write(6,"(//,'INVERSE_MATRIX: error in ZGETRI, info =',1x,i3)") info
      stop
    end if
	
    deallocate( ipiv,work )
  
  end subroutine inverse_matrix_z
  
  subroutine inverse_matrix_d( a )
    implicit none
    double precision, intent (inout) :: a(:,:)
    integer                          :: dim,info,lwork
    integer,          allocatable    :: ipiv(:)
    double precision, allocatable    :: work(:)
  
    dim = size( a,1 )
    lwork = 6*dim
    allocate( ipiv(dim),work(lwork) )

    call dgetrf( dim,dim,a,dim,ipiv,info )
    if( info /= 0 )then
      write(6,"(//,'INVERSE_MATRIX: error in DGETRF, info =',1x,i3)") info
      stop
    end if
    call dgetri( dim,a,dim,ipiv,work,lwork,info )
    if( info /= 0 )then
      write(6,"(//,'INVERSE_MATRIX: error in DGETRI, info =',1x,i3)") info
      stop
    end if
	
    deallocate( ipiv,work )
  
  end subroutine inverse_matrix_d
  
end module lapack_module
