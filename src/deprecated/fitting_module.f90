module fitting_module

contains
!Fitting routine for sum in fermionic frequency box:
subroutine fit_sum_fermionic_box(apart, w, chi, chi_sum_left)
  use parameters_module
  implicit none

  integer, intent(in) :: apart(:)
  complex(kind=8), intent(in) :: w(:), chi(maxdim,maxdim)
  complex(kind=8), intent(out) :: chi_sum_left(ndim2,maxdim)
  complex(kind=8),allocatable :: chi_sum_left_part(:,:,:)
  integer :: ipart, istart, istop, vstart, vstop, i1, i2, dum
  integer :: ipartmax,nreal,nimag
  double precision :: diff_fit, diff_scale

  ipartmax = size(apart)
  allocate(chi_sum_left_part(ipartmax,ndim2,maxdim))
  chi_sum_left_part = 0.d0
           
  do ipart=1,ipartmax
     istart=iwfmax_small-apart(ipart)
     istop=iwfmax_small-1+apart(ipart)
   
     do i1=1,ndim2
        do dum=istart,istop
           chi_sum_left_part(ipart,i1,:) = chi_sum_left_part(ipart,i1,:)+chi(i1+dum*ndim2,:)
        enddo
     enddo
  enddo

  vstart = iwfmax_small-apart(1)+15
  vstop = iwfmax_small-1+apart(1)-15
  chi_sum_left = 0d0
  do i1=1,ndim2
     do i2=1,vstart*ndim2-1
        chi_sum_left(i1,i2) = chi_sum_left_part(ipartmax,i1,i2)
     enddo

     do i2=vstart*ndim2,vstop*ndim2
        if (ipartmax .eq. 1) then
           !chi_sum_left(i1,i2) = chi_sum_left_part(ipartmax,i1,i2)
           chi_sum_left(i1,i2) = dot_product(w,chi_sum_left_part(:,i1,i2)) ! CHECK
        else
           nreal = count(dble(chi_sum_left_part(2:ipartmax,i1,i2)-chi_sum_left_part(1:ipartmax-1,i1,i2)) .ge. 0)
           nimag = count(dimag(chi_sum_left_part(2:ipartmax,i1,i2)-chi_sum_left_part(1:ipartmax-1,i1,i2)) .ge. 0)
           if ((nreal .eq. 0 .or. nreal .eq. ipartmax-1).and. (nimag .eq. 0 .or. nimag .eq. ipartmax-1)) then
              chi_sum_left(i1,i2) = dot_product(w,chi_sum_left_part(:,i1,i2))
           else if ((nreal .eq. 0 .or. nreal .eq. ipartmax-1)) then
              chi_sum_left(i1,i2) = cmplx(dble(dot_product(w,dble(chi_sum_left_part(:,i1,i2)))),&
                                          dimag(chi_sum_left_part(ipartmax,i1,i2)),KIND=8)
           else if ((nimag .eq. 0 .or. nimag .eq. ipartmax-1)) then
              chi_sum_left(i1,i2) = cmplx(dble(chi_sum_left_part(ipartmax,i1,i2)),&
                                          dble(dot_product(w,dimag(chi_sum_left_part(:,i1,i2)))),KIND=8)
           else
              chi_sum_left(i1,i2) = chi_sum_left_part(ipartmax,i1,i2)
           endif
        endif

        ! Check
        diff_fit = abs(chi_sum_left(i1,i2)-chi_sum_left_part(ipartmax,i1,i2))
        diff_scale = abs(chi_sum_left_part(ipartmax,i1,i2)-chi_sum_left_part(1,i1,i2))
        if(diff_fit .gt. 4.d0*diff_scale)then
          if (ounit .gt. 0) then
           write(ounit,'(a,2i4,2f16.10,10x,2f16.10)') 'failed fit?', i1,&
                   (i2/ndim2)-iwfmax_small,diff_fit,diff_scale,chi_sum_left(i1,i2)
           write(ounit,'("real:",999f16.10)') dble(chi_sum_left_part(:,i1,i2))
           write(ounit,'("imag:",999f16.10)') dimag(chi_sum_left_part(:,i1,i2))
          endif
          chi_sum_left(i1,i2) = chi_sum_left_part(ipartmax,i1,i2)
        endif
     enddo

     do i2=vstop*ndim2+1,maxdim
        chi_sum_left(i1,i2) = chi_sum_left_part(ipartmax,i1,i2)
     enddo
  enddo

  deallocate(chi_sum_left_part)

end subroutine fit_sum_fermionic_box


end module

