module lapack_module
  implicit none
  private
  public inverse_matrix,geometric_summation

  interface inverse_matrix
    module procedure inverse_matrix_d, inverse_matrix_z
  end interface inverse_matrix

  contains

  subroutine inverse_matrix_z(M, erstr, ierr)
    implicit none
    double complex, intent (inout)  :: M(:,:)
    integer, intent(out)            :: ierr
    character(len=200), intent(out) :: erstr
    integer                         :: ndim,lwork
    integer,        allocatable     :: ipiv(:)
    double complex, allocatable     :: work(:)
    double complex                  :: work_query(1)

    ndim = size(M,1)
    allocate(ipiv(ndim))
    call zgetrf(ndim,ndim,M,ndim,ipiv,ierr)
    if (ierr .ne. 0) then
      erstr = "ERROR in ZGETRF"
      return
    endif
    call zgetri(ndim,M,ndim,ipiv,work_query,-1,ierr) ! query for optimal work space
    if (ierr .ne. 0) then
      erstr = "ERROR in ZGETRI at workspace query"
      return
    endif
    lwork = work_query(1)
    allocate(work(lwork))
    call zgetri(ndim,M,ndim,ipiv,work,lwork,ierr)
    if (ierr .ne. 0) then
      erstr = "ERROR in ZGETRI"
      return
    endif
    deallocate(ipiv,work)
  end subroutine inverse_matrix_z

  subroutine inverse_matrix_d(M, erstr, ierr)
    implicit none
    double precision, intent (inout) :: M(:,:)
    integer, intent(out)             :: ierr
    character(len=200), intent(out)  :: erstr
    integer                          :: ndim,lwork
    integer,          allocatable    :: ipiv(:)
    double precision, allocatable    :: work(:)
    double precision                 :: work_query(1)

    ndim = size(M,1)
    allocate(ipiv(ndim))
    call dgetrf(ndim,ndim,M,ndim,ipiv,ierr)
    if(ierr .ne. 0) then
      erstr = "ERROR in DGETRF"
      return
    end if
    call dgetri(ndim,M,ndim,ipiv,work_query,-1,ierr) ! query for optimal work space
    if (ierr .ne. 0) then
      erstr = "ERROR in DGETRI at workspace query"
      return
    endif
    lwork = work_query(1)
    allocate(work(lwork))
    call dgetri(ndim,M,ndim,ipiv,work,lwork,ierr)
    if(ierr .ne. 0) then
      erstr = "ERROR in DGETRI"
      return
    end if
    deallocate(ipiv,work)
  end subroutine inverse_matrix_d


  ! Calculate the sum M + M**2 + ... + M**ord
  subroutine geometric_summation(M,ord)
    implicit none
    double complex, intent(inout) :: M(:,:)
    double complex,allocatable :: M_power(:,:),M_tmp(:,:),M_tmp_2(:,:)
    double complex :: alpha, beta
    integer, intent(in) :: ord
    integer :: mat_size,i,j
    mat_size = size(M,1)
    allocate(M_power(mat_size,mat_size),M_tmp(mat_size,mat_size),M_tmp_2(mat_size,mat_size))
    M_power=M
    M_tmp=M
    M_tmp_2=0.d0
    alpha=1.d0
    beta=0.d0
    do j=2,ord
      call zgemm('N','N',mat_size,mat_size,mat_size,alpha,M,mat_size,M_power,mat_size,beta,M_tmp_2,mat_size)
      M_power = M_tmp_2
      M_tmp = M_tmp + M_power
    enddo
    M = M_tmp
    deallocate(M_power,M_tmp,M_tmp_2)
  end subroutine geometric_summation

end module lapack_module
