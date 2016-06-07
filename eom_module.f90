module eom_module
  use parameters_module
  use one_particle_quant_module
  implicit none

contains

  subroutine calc_eom(interm3_dens,interm3_magn,gamma_dmft_dens,gamma_dmft_magn,gamma_loc_sum_left,sigma,kq_ind,iwb,iq,iw_data,u,v,u_tilde,hk,dc,siw)
  use parameters_module
  use lapack_module

  implicit none
  integer :: dum,i,j,iwf,iwf2,l,k,ik,ikq
  integer,intent(in) :: kq_ind(nkp,nqp),iwb,iq
  real*8,intent(in) :: iw_data(-iwmax:iwmax-1)
  complex(kind=8),intent(in) :: siw(-iwmax:iwmax-1,ndims)
  complex(kind=8),intent(in) :: hk(ndim,ndim,nkp)
  double precision,intent(in) :: dc(2,ndim)
  complex(kind=8),intent(in) :: u(ndim2,ndim2),u_tilde(ndim2,ndim2),v(ndim2,ndim2)
  complex(kind=8),intent(in) :: gamma_dmft_dens(ndim2,maxdim), gamma_dmft_magn(ndim2,maxdim)
  complex(kind=8) :: interm3_dens(ndim2,maxdim),interm3_magn(ndim2,maxdim)
  complex(kind=8) :: gamma_loc_sum_left(ndim2,maxdim)
  complex(kind=8),intent(inout) :: sigma(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp)
  complex(kind=8) :: m_tot_array(ndim,ndim,ndim,ndim,-iwfmax_small:iwfmax_small),m_tot(ndim2,maxdim)
  complex(kind=8) :: u_work(ndim2,ndim2), m_work(ndim2,maxdim)
  complex(kind=8) :: gkiw(ndim,ndim)
  complex(kind=8) :: alpha, delta

  !subtract -1 in the diagonal of the orbital blocks:
  do i1=1,ndim2
     do dum=0,2*iwfmax_small-1
        
        interm3_dens(i1,i1+dum*ndim2) = interm3_dens(i1,i1+dum*ndim2)-1.d0
        interm3_magn(i1,i1+dum*ndim2) = interm3_magn(i1,i1+dum*ndim2)-1.d0

     enddo
  enddo
  

  !call cpu_time(finish)
  !write(*,*)'doing matrix multiplication:', finish-start


  
  !call cpu_time(start)

  !EQUATION OF MOTION:
  m_work = 0.d0
  m_tot = 0.d0
  u_work = 0.d0

  !density part:
  u_work = v + u - 0.5d0*u_tilde

  !subtract dmft part (dmft self energy will be added in the end):
  interm3_dens = interm3_dens - gamma_dmft_dens
  
  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, interm3_dens, ndim2, delta, m_work, ndim2)

  m_tot = m_work

  !magnetic part:
  u_work = -1.5d0*u_tilde
  m_work = 0.d0

  !subtract dmft part:
  interm3_magn = interm3_magn - gamma_dmft_magn

  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, interm3_magn, ndim2, delta, m_work, ndim2)

  m_tot = m_tot + m_work

  !local part:
  u_work = v + u
  m_work = 0.d0

  gamma_loc_sum_left = gamma_loc_sum_left - gamma_dmft_dens
  
  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, gamma_loc_sum_left, ndim2, delta, m_work, ndim2)

  m_tot = m_tot - m_work

  !break up the compound index:
  m_tot_array = 0.d0

  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1 = i1+1
        i2 = 0
        do iwf2=-iwfmax_small,iwfmax_small-1
           do l=1,ndim
              do k=1,ndim
                 i2 = i2+1
                 
                 m_tot_array(i,j,k,l,iwf2) = m_tot(i1,i2)
              enddo
           enddo
        enddo
     enddo
  enddo

  
  !compute k-dependent self energy:
  do ik=1,nkp
     ikq = kq_ind(ik,iq) 
     do iwf=-iwfmax_small,iwfmax_small-1

        call get_gkiw(ikq, iwf, iwb, iw_data, siw, hk, dc, gkiw)
        
        do i=1,ndim
           do l=1,ndim
              do j=1,ndim
                 do k=1,ndim

                    sigma(i,l,iwf,ik) = sigma(i,l,iwf,ik)+m_tot_array(i,j,k,l,iwf)*gkiw(k,j)
                    !sigma(i,l,iwf,ik) = sigma(i,l,iwf,ik)+m_tot_array(i,j,j,l,iwf)*gkiw(j,j) !test
                 enddo
              enddo
           enddo
        enddo

     enddo
  enddo
             


end subroutine calc_eom



subroutine output_eom(iw_data,k_data,sigma_sum,sigma_loc)
  use parameters_module
  implicit none
  real*8 :: iw_data(-iwmax:iwmax-1)
  real*8 :: k_data(3,nkp)
  complex(kind=8) :: sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp)
  complex(kind=8) :: sigma_loc(ndim, ndim, -iwfmax_small:iwfmax_small-1)
  integer :: ik,iwf,i,iband

  open(34, file=trim(output_dir)//"siw_0_0_0_nostraight.dat", status='unknown')
  open(35, file=trim(output_dir)//"siw_0_0_0.5_nostraight.dat", status='unknown')
  open(36, file=trim(output_dir)//"siw_0_0.5_0.5_nostraight.dat", status='unknown')
  open(37, file=trim(output_dir)//"siw_loc_nostraight.dat", status='unknown')
  open(38, file=trim(output_dir)//"siw_iwf_0_nostraight.dat", status='unknown')
  open(39, file=trim(output_dir)//"siw_0.5_0.5_0.5_nostraight.dat", status='unknown')
  open(40, file=trim(output_dir)//"siw_all_nostraight.dat",status='unknown')

  do ik=1,100
     write(38,'(100F12.6)') k_data(2,ik), k_data(3,ik), (real(sigma_sum(i,i,0,ik)), aimag(sigma_sum(i,i,0,ik)), i=1,3)
  enddo 

  do iwf=-iwfmax_small,iwfmax_small-1
     write(34,'(100F12.6)')iw_data(iwf), (real(sigma_sum(i,i,iwf,1)), aimag(sigma_sum(i,i,iwf,1)), i=1,3)
     write(35,'(100F12.6)')iw_data(iwf), (real(sigma_sum(i,i,iwf,6)), aimag(sigma_sum(i,i,iwf,6)),i=1,3)
     write(36,'(100F12.6)')iw_data(iwf), (real(sigma_sum(i,i,iwf,56)), aimag(sigma_sum(i,i,iwf,56)),i=1,3)
     write(37,'(100F12.6)')iw_data(iwf), (real(sigma_loc(i,i,iwf)), aimag(sigma_loc(i,i,iwf)),i=1,3)
     write(39,'(100F12.6)')iw_data(iwf), (real(sigma_sum(i,i,iwf,556)),aimag(sigma_sum(i,i,iwf,556)),i=1,3)
  enddo

  do iwf=-iwfmax_small,iwfmax_small-1
    do ik=1,nkp
      write(40,'(100F12.6)') dble(iwf), iw_data(iwf), dble(ik), k_data(1,ik), k_data(2,ik), k_data(3,ik), (real(sigma_sum(i,i,iwf,ik)), aimag(sigma_sum(i,i,iwf,ik)),i=1,3)
    end do
  end do


  close(34)
  close(35)
  close(36)
  close(37)
  close(38)
  close(39)
  close(40)

end subroutine output_eom


end module

