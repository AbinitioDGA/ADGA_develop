module eom_module

  use parameters_module
  use one_particle_quant_module
  use kq_tools
  implicit none

contains

!================================================================================================
  subroutine calc_eom(interm3_dens, interm3_magn, gamma_dmft_dens, gamma_dmft_magn, gamma_loc_sum_left, sigma, kq_ind, iwb, iq, iw_data, u, v, u_tilde, hk, dc, siw)
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
  
  !EQUATION OF MOTION:
  m_work = 0.d0
  m_tot = 0.d0
  u_work = 0.d0

  !density part:
  u_work = v + u - 0.5d0*u_tilde
  interm3_dens = interm3_dens - gamma_dmft_dens !subtract dmft part (dmft self energy will be added in the end)
  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, interm3_dens, ndim2, delta, m_work, ndim2)
  m_tot = m_work

  !magnetic part:
  u_work = -1.5d0*u_tilde
  m_work = 0.d0
  interm3_magn = interm3_magn - gamma_dmft_magn !subtract dmft part:
  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, interm3_magn, ndim2, delta, m_work, ndim2)
  m_tot = m_tot + m_work

  !local part:
  u_work = u
  m_work = 0.d0
  gamma_loc_sum_left = gamma_loc_sum_left - gamma_dmft_dens
  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, gamma_loc_sum_left, ndim2, delta, m_work, ndim2)
  m_tot = m_tot - m_work

  !v(q) part
  u_work = v
  m_work = 0.d0
  alpha = 1.d0
  delta = 0.d0
  call zgemm('n', 'n', ndim2, maxdim, ndim2, alpha, u_work, ndim2, gamma_dmft_dens, ndim2, delta, m_work, ndim2)
  m_tot = m_tot + m_work
  

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

  
  !compute k-dependent self energy (convolution with Greens function gkiw): 
  do ik=1,nkp
     ikq = kq_ind(ik,iq) 
     do iwf=-iwfmax_small,iwfmax_small-1
        call get_gkiw(ikq, iwf, iwb, iw_data, siw, hk, dc, gkiw)
        
        do i=1,ndim
           do l=1,ndim
              do j=1,ndim
                 do k=1,ndim
                    sigma(i,l,iwf,ik) = sigma(i,l,iwf,ik)+m_tot_array(i,j,k,l,iwf)*gkiw(k,j)
                 enddo
              enddo
           enddo
        enddo

     enddo
  enddo
             
end subroutine calc_eom
!=============================================================================================


!=============================================================================================
subroutine add_siw_dmft(siw, sigma_sum, sigma_loc) 
  implicit none
  complex(kind=8), intent(in) :: siw(-iwmax:iwmax-1,ndims) 
  complex(kind=8) :: sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp)
  complex(kind=8) :: sigma_loc(ndim, ndim, -iwfmax_small:iwfmax_small-1)
  integer :: ik, iwf, iband

 ! local contribution is replaced by the DMFT self energy for better asymptotics
    do ik=1,nkp
       do iwf=-iwfmax_small,iwfmax_small-1
          do iband=1,ndim
             sigma_sum(iband, iband, iwf, ik) = sigma_sum(iband, iband, iwf, ik) + siw(iwf, iband)
          enddo
       enddo
    enddo

    sigma_loc = 0.d0
    do ik=1,nkp
      sigma_loc(:,:,:) = sigma_loc(:,:,:)+sigma_sum(:,:,:,ik)
    enddo
    sigma_loc = sigma_loc/dble(nkp)

end subroutine add_siw_dmft
!===============================================================================================


!==============================================================================================
subroutine output_eom(iw_data, k_data, sigma_sum, sigma_loc)
  implicit none

  real*8, intent(in) :: iw_data(-iwmax:iwmax-1)
  real*8, intent(in) :: k_data(3,nkp)
  complex(kind=8), intent(in) :: sigma_sum(ndim, ndim, -iwfmax_small:iwfmax_small-1, nkp)
  complex(kind=8), intent(in) :: sigma_loc(ndim, ndim, -iwfmax_small:iwfmax_small-1)
  complex(kind=8) :: sigma_tmp(ndim*(ndim+1)/2)
  integer :: ik, iwf, i, j, iband,nkp_eom,ii
  character(len=50) :: eom_format

  !TODO generate the filenames automatically
  open(34, file=trim(output_dir)//"siw_0_0_0.dat", status='unknown')
  open(35, file=trim(output_dir)//"siw_0_0_0.5.dat", status='unknown')
 
  open(44, file=trim(output_dir)//"siw_loc.dat", status='unknown')
  open(45, file=trim(output_dir)//"siw_all_k.dat",status='unknown')

  !TODO read symmetry points from file
  do iwf=-iwfmax_small,iwfmax_small-1
     write(34,'(100F12.6)')iw_data(iwf), ((real(sigma_sum(i,j,iwf,1)), aimag(sigma_sum(i,j,iwf,1)), j=i,ndims), i=1,ndims)
     write(35,'(100F12.6)')iw_data(iwf), ((real(sigma_sum(i,j,iwf,6)), aimag(sigma_sum(i,j,iwf,6)), j=i,ndims), i=1,ndims)
     write(44,'(100F12.6)')iw_data(iwf), ((real(sigma_loc(i,j,iwf)), aimag(sigma_loc(i,j,iwf)), j=i,ndims), i=1,ndims)
  enddo

  do ik=1,nkp
     do iwf=-iwfmax_small,iwfmax_small-1
        write(45,'(100F12.6)') k_data(1,ik), k_data(2,ik), k_data(3,ik), iw_data(iwf), &
            ((real(sigma_sum(i,j,iwf,ik)), aimag(sigma_sum(i,j,iwf,ik)), j=i,ndims), i=1,ndims)
     enddo
  enddo

  if (k_path_eom) then
    nkp_eom=n_segments()*nkp1/2+1
    allocate(k_data_eom(nkp_eom))
    call generate_q_path(nkp1,k_data_eom)
    write(*,*) k_data_eom
    open(unit=46,file='skiw_path.dat')
    write(eom_format,'((A)I2(A))') '(I5,2X,F12.6,2X,',2*ndim*(ndim+1)/2,'(F12.6,2X))'
    do ik=1,nkp_eom
      do iwf=-iwfmax_small,iwfmax_small-1
        ii=0
        do i1=1,ndim
          do i2=1,i1
            ii=ii+1
            sigma_tmp(ii)=sigma_sum(i1,i2,iwf,k_data_eom(ik))
          end do
        end do
        write(46,eom_format) k_data_eom(ik),iw_data(iwf),sigma_tmp
      end do
    end do
    close(46)
  end if

  close(34)
  close(35)
  close(44)
  close(45)

end subroutine output_eom
!================================================================================================

end module

