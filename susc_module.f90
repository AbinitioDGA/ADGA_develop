module susc_module
  implicit none

  contains

    subroutine calc_chi_qw(chi_qw,interm3,chi0_sum)
      use parameters_module
      implicit none
      complex(kind=8), intent(out) :: chi_qw(ndim2,ndim2)
      complex(kind=8), intent(in) :: interm3(ndim2,maxdim)
      complex(kind=8), intent(in) :: chi0_sum(ndim2,ndim2,-iwfmax_small:iwfmax_small-1)
      integer :: iwf,iband,i1,i2
      do i1=1,ndim2
         do i2=1,ndim2
            do iwf=1,2*iwfmax_small
               do iband=1,ndim2
                  chi_qw(i1,i2)=chi_qw(i1,i2)+interm3(i1,(iwf-1)*ndim2+iband)*chi0_sum(iband,i2,iwf-1-iwfmax_small)
               end do
            end do
         end do
      end do
      chi_qw=chi_qw/(beta**2)
    end subroutine calc_chi_qw

    subroutine calc_bubble(bubble,chi0_sum)
      use parameters_module
      implicit none
      complex(kind=8), intent(in) :: chi0_sum(ndim2,ndim2,iwstart:iwstop)
      complex(kind=8) :: bubble(ndim2,ndim2)
      integer :: iwf
      bubble=0.d0
      do iwf=iwstart,iwstop
        bubble(:,:)=bubble(:,:)+chi0_sum(:,:,iwf)
      end do
      bubble=bubble/(beta**2)
    end subroutine calc_bubble

    ! subroutine to output susceptibility
    subroutine output_chi_qw(chi_qw,iwb_data,qw,filename_output)
      use parameters_module
      implicit none
      character(len=*) :: filename_output
      character(len=100) :: format_str
      real*8 :: iwb_data(-iwbmax:iwbmax), chi_qw_1q(2*ndim2**2)
      complex(kind=8) :: chi_qw(ndim2,ndim2,nqp*(2*iwbmax_small+1))
      integer :: iwb,iq,qw(2,nqp*(2*iwbmax+1)),i
      integer :: i1

      write(*,*) 'Output chi_qw'

      ! create a format string that works for various orbital dimensions
      write(format_str,'((A)I3(A))') '(I5,2X,E14.7E2,2X,I5,2X,',2*ndim2**2+3,'(E14.7E2,2X))'

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
    subroutine output_chi_qw_1line(chi_qw,iwb_data,qw,iqw,rank,filename_output)
      use parameters_module
      implicit none
      character(len=*) :: filename_output
      character(len=100) :: format_str
      real*8 :: iwb_data(-iwbmax:iwbmax), chi_qw_1q(2*ndim2**2)
      complex(kind=8) :: chi_qw(ndim2,ndim2)
      integer :: iwb,iq,qw(2,nqp*(2*iwbmax+1)),i,iqw,rank,un

      ! create a format string that works for various orbital dimensions
      write(format_str,'((A)I3(A))') '(I5,2X,E14.7E2,2X,I5,2X,',2*ndim2**2+3,'(E14.7E2,2X))'

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

    subroutine output_chi_loc(chi_qw,iwb_data,filename_output)
      use parameters_module
      implicit none
      character(len=*) :: filename_output
      character(len=100) :: format_str
      real*8 :: iwb_data(-iwbmax:iwbmax), chi_qw_1q(2*ndim2**2)
      complex(kind=8) :: chi_qw(ndim2,ndim2,2*iwbmax_small+1)
      integer :: iwb,i,i1

      write(*,*) 'output chi_loc'

      ! create a format string that works for various orbital dimensions
      write(format_str,'((A)I3(A))') '(I5,2X,E14.7E2,2X',2*ndim2**2,'(E14.7E2,2X))'

      ! open file and write head line
      open(unit=10,file=trim(output_dir)//filename_output)
      ! the file header should contain important information
      write(10,*) '#iwbmax, nqp, ndim, beta, mu'
      write(10,'(I5,2X,I5,2X,I5,2X,E14.7E2,2X,E14.7E2)') iwbmax_small,0,ndim,beta,mu!nqp=0 can be used as a flag for the postprocessing
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

    subroutine output_chi_loc_1line(chi_qw,iwb_data,iwb,rank,filename_output)
      use parameters_module
      implicit none
      character(len=*) :: filename_output
      character(len=100) :: format_str
      real*8 :: iwb_data(-iwbmax:iwbmax), chi_qw_1q(2*ndim2**2)
      complex(kind=8) :: chi_qw(ndim2,ndim2)
      integer :: iwb,i,un,rank

      ! create a format string that works for various orbital dimensions
      write(format_str,'((A)I3(A))') '(I5,2X,E14.7E2,2X',2*ndim2**2,'(E14.7E2,2X))'

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

    subroutine calc_bubble_old(bubble,iwb,iq,kq_ind,iw_data,siw,hk,dc)
      use parameters_module
      use one_particle_quant_module
      implicit none
      integer,intent(in) :: iwb,iq,kq_ind(nkp,nqp)
      double precision,intent(in) :: iw_data(-iwmax:iwmax-1),dc(2,ndim)
      complex(kind=8),intent(in) :: siw(-iwmax:iwmax-1,ndims), hk(ndim,ndim,nkp)
      complex(kind=8) :: bubble(ndim2,ndim2),g1(ndim,ndim),g2(ndim,ndim)
      integer :: ik,iwf,ikq1,ikq2
      integer :: i1,i2,i3,i4
      bubble=0.d0
      do ik=1,nkp
        do iwf=-iwfmax_small,iwfmax_small-1
          call get_gkiw(kq_ind(ik,1),  iwf, 0,   iw_data, siw, hk, dc, g1)
          call get_gkiw(kq_ind(ik,iq), iwf, iwb, iw_data, siw, hk, dc, g2)
          do i1=1,ndim
            do i2=1,ndim
              do i3=1,ndim
                do i4=1,ndim
                  bubble(i3+ndim*(i1-1),i4+ndim*(i2-1))=bubble(i3+ndim*(i1-1),i4+ndim*(i2-1))+g1(i1,i2)*g2(i3,i4)
                end do
              end do
            end do
          end do
        end do
      end do
    end subroutine calc_bubble_old

end module
