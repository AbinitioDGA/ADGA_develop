module susc_module
  implicit none

contains

subroutine calc_chi_qw(chi_qw,interm3,chi0_sum)
   use parameters_module
   implicit none
   complex(kind=8), intent(inout) :: chi_qw(ndim2,ndim2)
   complex(kind=8), intent(in)    :: interm3(ndim2,maxdim)
   complex(kind=8), intent(in)    :: chi0_sum(ndim2,ndim2,-iwfmax_small:iwfmax_small-1)
   complex(kind=8)                :: chi_tmp(ndim2,ndim2)
   integer :: dum,iwf,i3,i1,i2,i
   ! Initialize
   chi_tmp = 0
   do i2=1,ndim2
      do dum=0,2*iwfmax_small-1
         do i3=1,ndim2
            iwf = dum - iwfmax_small
            i = i3 + dum*ndim2 ! = {i3,iwf} 
            chi_tmp(:,i2)=chi_tmp(:,i2)+interm3(:,i)*chi0_sum(i3,i2,iwf)
         end do
      end do
   end do
   ! Add to chi_qw
   chi_qw=chi_qw + chi_tmp/(beta**2)
   return
end subroutine calc_chi_qw

!subroutine calc_bubble(bubble,chi0_sum)
!  use parameters_module
!  implicit none
!  complex(kind=8), intent(in) :: chi0_sum(ndim2,ndim2,iwstart:iwstop)
!  complex(kind=8) :: bubble(ndim2,ndim2)
!  integer :: iwf
!  bubble=0.d0
!  do iwf=iwstart,iwstop
!    bubble(:,:)=bubble(:,:)+chi0_sum(:,:,iwf)
!  end do
!  bubble=bubble/(beta**2)
!end subroutine calc_bubble

! subroutine to output susceptibility
subroutine output_chi_qw(chi_qw,qw,filename_output)
  use parameters_module
  implicit none
  complex(kind=8),intent(in)  :: chi_qw(ndim2,ndim2,nqp*(2*iwbmax_small+1))
  integer,intent(in)          :: qw(2,nqp*(2*iwbmax+1))
  character(len=*),intent(in) :: filename_output
  character(len=100) :: format_str
  real*8 :: chi_qw_1q(2*ndim2**2)
  integer :: iwb,iq,i
  integer :: i1

  if (ounit .ge. 1 .and. (verbose .and. (index(verbstr,"Output") .ne. 0))) then
   write(ounit,*) 'Output chi_qw'
  endif

  ! create a format string that works for various orbital dimensions
  write(format_str,'(A,I5,A)') '(I5,2X,E14.7E2,2X,I5,2X,',2*ndim2**2+3,'(E14.7E2,2X))'

  ! open file and write head line
  open(unit=10,file=trim(output_dir)//filename_output)
  ! the file header should contain important information
  write(10,*) '#iwbmax, nqp, ndim, beta, mu'
  write(10,'(I5,2X,I5,2X,I5,2X,E14.7E2,2X,E14.7E2)') iwbmax_small,nqp,ndim,beta,mu
  write(10,*) '#iwb  ','wb    ','iq  ','      (q)      ','chi_qw'

  ! loop over all entries
  do i1=1,nqp*(2*iwbmax_small+1)
     iq = qw(2,i1)
     iwb = qw(1,i1)
     do i=1,ndim2**2!flatten the band matrix into one output line
        chi_qw_1q(2*i-1) =  real(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1,i1),8)
        chi_qw_1q(2*i)   = dimag(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1,i1))
     end do
     write(10,format_str) iwb,iwb_data(iwb),iq,k_data(:,q_data(iq)),chi_qw_1q

     ! insert an empty line after each omega block. could be useful for plotting.
     if (mod(i1,nqp).eq.0) then
        write(10,*) ' '
     end if

  end do

  close(10)
end subroutine output_chi_qw

! subroutine to output one line of susceptibility
subroutine output_chi_qw_1line(chi_qw,qw,iqw,rank,filename_output)
  use parameters_module
  implicit none
  character(len=*) :: filename_output
  character(len=100) :: format_str
  real*8 :: chi_qw_1q(2*ndim2**2)
  complex(kind=8) :: chi_qw(ndim2,ndim2)
  integer :: iwb,iq,qw(2,nqp*(2*iwbmax+1)),i,iqw,rank,un

  ! create a format string that works for various orbital dimensions
  write(format_str,'(A,I5,A)') '(I5,2X,E14.7E2,2X,I5,2X,',2*ndim2**2+3,'(E14.7E2,2X))'

  ! open file and write head line
  un=10000+rank
  open(unit=un,file=trim(output_dir)//filename_output,access='append')
     iq = qw(2,iqw)
     iwb = qw(1,iqw)
     do i=1,ndim2**2!flatten the band matrix into one output line
        chi_qw_1q(2*i-1) =  real(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1),8)
        chi_qw_1q(2*i)   = dimag(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1))
     end do
     write(un,format_str) iwb,iwb_data(iwb),iq,k_data(:,q_data(iq)),chi_qw_1q
  close(un)
end subroutine output_chi_qw_1line

subroutine output_chi_loc(chi_qw,filename_output)
  use parameters_module
  implicit none
  complex(kind=8),intent(in) :: chi_qw(ndim2,ndim2,2*iwbmax_small+1)
  character(len=*),intent(in) :: filename_output
  character(len=100) :: format_str
  real*8 :: chi_qw_1q(2*ndim2**2)
  integer :: iwb,i,i1

  if (ounit .ge. 1 .and. (verbose .and. (index(verbstr,"Output") .ne. 0))) then
   write(ounit,*) 'output chi_loc'
  endif

  ! create a format string that works for various orbital dimensions
  write(format_str,'(A,I5,A)') '(I5,2X,E14.7E2,2X',2*ndim2**2,'(E14.7E2,2X))'

  ! open file and write head line
  open(unit=10,file=trim(output_dir)//filename_output)
  ! the file header should contain important information
  write(10,*) '#iwbmax, nqp, ndim, beta, mu'
  write(10,'(I5,2X,I5,2X,I5,2X,E14.7E2,2X,E14.7E2)') iwbmax_small,0,ndim,beta,mu !nqp=0 can be used for the postprocessing
  write(10,*) '#iwb  ','wb    ','chi_loc'

  ! loop over all entries
  do i1=1,2*iwbmax_small+1
     iwb = i1-iwbmax_small-1
     do i=1,ndim2**2!flatten the band matrix into one output line
        chi_qw_1q(2*i-1) =  real(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1,i1),8)
        chi_qw_1q(2*i)   = dimag(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1,i1))
     end do
     write(10,format_str) iwb,iwb_data(iwb),chi_qw_1q

  end do
  close(10)
end subroutine output_chi_loc

subroutine output_chi_loc_1line(chi_qw,iwb,rank,filename_output)
  use parameters_module
  implicit none
  character(len=*) :: filename_output
  character(len=100) :: format_str
  real*8 :: chi_qw_1q(2*ndim2**2)
  complex(kind=8) :: chi_qw(ndim2,ndim2)
  integer :: iwb,i,un,rank

  ! create a format string that works for various orbital dimensions
  write(format_str,'(A,I5,A)') '(I5,2X,E14.7E2,2X',2*ndim2**2,'(E14.7E2,2X))'

  ! open file and write head line
  un=20000+rank
  open(unit=un,file=trim(output_dir)//filename_output,access='append')

     do i=1,ndim2**2!flatten the band matrix into one output line
        chi_qw_1q(2*i-1) =  real(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1),8)
        chi_qw_1q(2*i)   = dimag(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1))
     end do
     write(un,format_str) iwb,iwb_data(iwb),chi_qw_1q

  close(un)
end subroutine output_chi_loc_1line

!subroutine calc_bubble_old(bubble,iwb,iq,kq_ind)
!  use parameters_module
!  use one_particle_quant_module
!  implicit none
!  integer,intent(in) :: iwb,iq,kq_ind(nkp,nqp)
!  complex(kind=8) :: bubble(ndim2,ndim2),g1(ndim,ndim),g2(ndim,ndim)
!  integer :: ik,iwf,ikq1,ikq2
!  integer :: i1,i2,i3,i4
!  bubble=0.d0
!  do ik=1,nkp
!    do iwf=-iwfmax_small,iwfmax_small-1
!      call get_gkiw(kq_ind(ik,1),  iwf, 0,   g1)
!      call get_gkiw(kq_ind(ik,iq), iwf, iwb, g2)
!      do i1=1,ndim
!        do i2=1,ndim
!          do i3=1,ndim
!            do i4=1,ndim
!              bubble(i3+ndim*(i1-1),i4+ndim*(i2-1))=bubble(i3+ndim*(i1-1),i4+ndim*(i2-1))+g1(i1,i2)*g2(i3,i4)
!            end do
!          end do
!        end do
!      end do
!    end do
!  end do
!end subroutine calc_bubble_old



end module
