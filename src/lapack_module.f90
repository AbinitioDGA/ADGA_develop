module lapack_module
  implicit none
  private
  public :: inverse_matrix, eigenvalues_matrix

  interface inverse_matrix
    module procedure inverse_matrix_d, inverse_matrix_z
  end interface inverse_matrix

  ! right eigenvalues
  interface eigenvalues_matrix
    module procedure eigenvalues_matrix_z
  end interface eigenvalues_matrix

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

  subroutine eigenvalues_matrix_z(M, Ve, Va, erstr, ierr)
    implicit none
    double complex, intent(inout)   :: M(:,:), Ve(:,:), Va(:) ! matrix, vectors, values
    integer, intent(out)            :: ierr
    character(len=200), intent(out) :: erstr
    integer                         :: ndim, lwork
    double complex, allocatable     :: work(:)
    double complex                  :: work_query(1)
    double complex                  :: dummy
    double precision, allocatable   :: rwork(:)

    ndim = size(M,1) ! identical to maxdim
    allocate(rwork(2*ndim))
    call zgeev('n','v',ndim,M,ndim,Va,Ve,ndim,Ve,ndim,work_query,-1,rwork,ierr)
    if (ierr .ne. 0) then
      erstr = "ERROR in ZGEEV at workspace query"
      return
    endif
    lwork = work_query(1)
    allocate(work(lwork))
    call zgeev('n','v',ndim,M,ndim,Va,Ve,ndim,Ve,ndim,work,lwork,rwork,ierr)
    if (ierr .ne. 0) then
      erstr = "ERROR in ZGEEV"
      return
    endif
    deallocate(work,rwork)
  end subroutine eigenvalues_matrix_z

end module lapack_module
