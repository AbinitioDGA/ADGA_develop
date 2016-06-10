module susc_module
  use parameters_module
  implicit none

  contains

    subroutine calc_chi_qw(chi_qw,interm3,chi0_sum)
      use parameters_module
      implicit none
      complex(kind=8), intent(out) :: chi_qw(ndim2,ndim2)
      complex(kind=8), intent(in) :: interm3(ndim2,maxdim)
      complex(kind=8), intent(in) :: chi0_sum(ndim2,ndim2,-iwfmax_small:iwfmax_small-1)

      do i1=1,2*iwfmax_small
         do i2=1,ndim2
            do i3=1,ndim2
               do i4=1,ndim2
                  chi_qw(i2,i3)=chi_qw(i2,i3)+interm3(i2,(i1-1)*ndim2+i4)*chi0_sum(i4,i3,i1-1-iwfmax_small)
               end do
            end do
         end do
      end do

    end subroutine calc_chi_qw


    ! subroutine to output susceptibility
    subroutine output_chi_qw(chi_qw,iwb_data,q_data,qw,filename_output)
      use parameters_module
      implicit none
      character(len=*) :: filename_output
      character(len=100) :: format_str
      real*8 :: iwb_data(-iwbmax:iwbmax), q_data(3,nqp), chi_qw_1q(2*ndim2**2)
      complex(kind=8) :: chi_qw(ndim2,ndim2,nqp*(2*iwbmax_small+1))
      integer :: iwb,iq,qw(2,nqp*(2*iwbmax+1)),i

      ! create a format string that works for various orbital dimensions
      write(format_str,'((A)I3(A))') '(I5,2X,E14.7E2,2X,I5,2X,',2*ndim2**2+3,'(E14.7E2,2X))'

      ! open file and write head line
      open(unit=10,file=trim(output_dir)//filename_output)
      write(10,*) '#iwb  ','wb    ','iq  ','      (q)      ','chi_qw'
  
      ! loop over all entries
      do i1=1,nqp*(2*iwbmax_small+1)
         iq = qw(2,i1)
         iwb = qw(1,i1)
         do i=1,ndim2**2!flatten the band matrix into one output line
            chi_qw_1q(2*i-1) =  real(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1,i1),8)
            chi_qw_1q(2*i)   = dimag(chi_qw((i-1)/ndim2+1,mod(i-1,ndim2)+1,i1))
         end do
         write(10,format_str) iwb,iwb_data(iwb),iq,q_data(:,iq),chi_qw_1q

         ! insert an empty line after each omega block. could be useful for plotting.
         if (mod(i1,nqp).eq.0) then
            write(10,*) ' '
         end if

      end do

      close(10)
    end subroutine output_chi_qw

end module

  
