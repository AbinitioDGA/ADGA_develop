module mpi_org

  use parameters_module, only: iwbmax_small,nqp,ndim2
#ifdef MPI
    use mpi
#endif
  implicit none
  public

  integer :: mpi_wrank
  integer :: mpi_wsize
  integer :: master
  integer, allocatable :: rct(:),disp(:)
  integer :: ierr

  contains

  subroutine mpi_initialize()
    implicit none
#ifdef MPI
      call MPI_init(ierr)
      call MPI_comm_rank(MPI_COMM_WORLD,mpi_wrank,ierr)
      call MPI_comm_size(MPI_COMM_WORLD,mpi_wsize,ierr)
      master = 0
      allocate(rct(mpi_wsize),disp(mpi_wsize))
#else
      mpi_wrank=0
      mpi_wsize=1
      master=0
#endif
  end subroutine mpi_initialize


  subroutine mpi_distribute()
    implicit none
    integer :: qwstart,qwstop
    integer :: i,j
#ifdef MPI
      rct=0
      disp=0
      qwstop = 0
      do i=0,mpi_wrank
        j = (nqp*(2*iwbmax_small+1) - qwstop)/(mpi_wsize-i)
        qwstart = qwstop + 1
        qwstop = qwstop + j
      enddo
      rct(mpi_wrank+1)=(qwstop-qwstart+1)*ndim2**2
      if (mpi_wrank .eq. master) then
        call MPI_reduce(MPI_IN_PLACE,rct,mpi_wsize,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierr)
      else
        call MPI_reduce(rct,rct,mpi_wsize,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierr)
      end if
      if (mpi_wrank .eq. master) then
        do i=2,mpi_wsize
          disp(i)=sum(rct(1:i-1))! the first displacing has to be 0
        end do
        write(*,*) 'receive ct',rct
        write(*,*) 'displacing',disp
      end if
#else
      qwstart = 1
      qwstop = nqp*(2*iwbmax_small+1)
#endif
    write(*,*)'rank=',mpi_wrank, 'qwstart=', qwstart, 'qwstop=', qwstop
  end subroutine mpi_distribute


  subroutine mpi_close()
    implicit none
#ifdef MPI
      deallocate(rct,disp)
      call MPI_finalize(ierr)
#endif
  end subroutine mpi_close

end module mpi_org
