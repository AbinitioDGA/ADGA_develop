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
  integer :: qwstart,qwstop

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

  subroutine mpi_stop(stopmsg,er,mpi_er)
   use parameters_module, ONLY: ounit
   implicit none
   character(len=*),intent(in)  :: stopmsg      !< message to print
   integer,optional,intent(in)  :: er           !< error flag to print
   integer,optional,intent(in)  :: mpi_er          !< rank of the processor printing this message
#ifdef MPI
   character(len=MPI_MAX_ERROR_STRING) :: str      ! temp string
   integer                             :: strlen      ! length
   integer                             :: mpi_errorflag   ! temp error flag
   character(len=5)                    :: fancynr  ! temp string for the mpi rank
   logical                             :: unitopened ! tells if the file unit is open
#endif
      inquire(unit=ounit,opened=unitopened)
      if (unitopened) call flush(ounit)
#if MPI 
      write(fancynr,'(i5)') mpi_wrank
      if (.not. present(mpi_er)) then
         if (present(er)) then
            if (unitopened) then
               write(ounit,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," er=",i5)') TRIM(ADJUSTL(fancynr)),stopmsg,er
               close(ounit)
            endif
            write(*,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," er=",i5)') TRIM(ADJUSTL(fancynr)),stopmsg,er
         else
            if (unitopened) then
               write(ounit,'(1x,"Rank ",a,", STOP in abinitiodga: ",a)') TRIM(ADJUSTL(fancynr)),stopmsg
               close(ounit)
            endif
            write(*,'(1x,"Rank ",a,", STOP in abinitiodga: ",a)') TRIM(ADJUSTL(fancynr)),stopmsg
         endif
         call sleep(1)
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_errorflag)
      else
         if (unitopened) then
            if (present(er)) then
               write(ounit,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," er=",i25)') TRIM(ADJUSTL(fancynr)),stopmsg,er
            endif
            write(ounit,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," mpi=",i5)') TRIM(ADJUSTL(fancynr)),stopmsg,mpi_er
            close(ounit)
         endif
         if (present(er)) write(*,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," er=",i25)') TRIM(ADJUSTL(fancynr)),stopmsg,er
         write(*,'(1x,"Rank ",a,", STOP in abinitiodga: ",a," mpi=",i5)') TRIM(ADJUSTL(fancynr)),stopmsg,mpi_er
         call MPI_ERROR_STRING(mpi_er, str, strlen,mpi_errorflag)
         if (mpi_errorflag .ne. 0) then
            write(*,*) "Error in MPI_ERROR_STRING:",mpi_errorflag
         else
            write(*,*) "STOP in abinitiodga:" // str(1:strlen)
         endif
         call sleep(1)
         call MPI_Abort(MPI_COMM_WORLD,mpi_er,mpi_errorflag)
      endif
#else
      if (present(er)) then
         write(*,*) "STOP in abinitiodga: " // stopmsg // " er=",er
      else if (.not. present(mpi_er)) then
         write(*,*) "STOP in abinitiodga: " // stopmsg
      endif
      if (present(mpi_er)) then
         write(*,*) "ERROR: Got an mpi error but abinitiodga is compiled without MPI support!"
         write(*,*) "STOP in abinitiodga: " // stopmsg // " mpi=",mpi_er
      endif
#endif

      stop
   end subroutine mpi_stop

end module mpi_org
