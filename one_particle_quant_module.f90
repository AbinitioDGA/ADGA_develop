module one_particle_quant_module
  use lapack_module
  use parameters_module
  use aux
  implicit none

  contains

subroutine get_giw()
  use lapack_module
  use parameters_module
  implicit none
  integer :: ik, iw, i
  complex(kind=8) :: g(ndim,ndim), g2(ndim,ndim)

  giw = 0.d0
  do ik=1,nkp
     g(:,:) = -hk(:,:,ik)
     do iw=0,iwmax-1 !use symmetry of giw(-w)=giw^*(w)
        do i=1,ndim
           g(i,i) = ci*iw_data(iw)+mu-hk(i,i,ik)-dc(1,i)
        enddo
        do i=1,ndim
           g(i,i) = g(i,i)-siw(iw,i) !no spin dependence in single particle Greens function
        enddo
        g2 = g(:,:)
        call inverse_matrix(g2)
        do i=1,ndim
           giw(iw,i) = giw(iw,i)+g2(i,i)
        enddo
     enddo
  enddo

  do iw=0,iwmax-1
     do i=1,ndim
        giw(-iw-1,i) = real(giw(iw,i),kind=8)-ci*aimag(giw(iw,i))
     enddo
  enddo

  giw = giw/dble(nkp)

end subroutine get_giw



subroutine get_gkiw(ikq, iwf, iwb, gkiw)
  use lapack_module
  use parameters_module
  implicit none
  integer :: i
  integer :: iwf, iwb, ikq
  complex(kind=8), intent(out) :: gkiw(ndim,ndim)


  gkiw(:,:) = -hk(:,:,ikq)
  do i=1,ndim
     gkiw(i,i) = ci*iw_data(iwf-iwb)+mu-hk(i,i,ikq)-dc(1,i)
  enddo
  do i=1,ndim
  gkiw(i,i) = gkiw(i,i)-siw(iwf-iwb,i)
  enddo
  call inverse_matrix(gkiw)


end subroutine get_gkiw


subroutine get_gloc(sigma_sum, gloc)
  use lapack_module
  use parameters_module
  implicit none

  complex(kind=8) :: sigma_sum(ndim,ndim,-iwfmax_small:iwfmax_small-1,nkp)
  integer :: ik, iw, i
  complex(kind=8) :: g(ndim,ndim), g2(ndim,ndim)
  complex(kind=8), intent (out) :: gloc(-iwmax:iwmax-1,ndim,ndim)

  gloc = 0.d0
  do ik=1,nkp
     g(:,:) = -hk(:,:,ik)
     do iw=0,iwfmax_small-1
        do i=1,ndim
           g(i,i) = ci*iw_data(iw)+mu-hk(i,i,ik)-dc(1,i)
        enddo
        do i=1,ndim
           g(i,i) = g(i,i)-sigma_sum(i,i,iw,ik)
        enddo
        g2 = g(:,:)
        call inverse_matrix(g2)
        gloc(iw,:,:) = gloc(iw,:,:)+g2(:,:)
     enddo

     do iw=iwfmax_small,iwmax-1
        do i=1,ndim
           g(i,i) = ci*iw_data(iw)+mu-hk(i,i,ik)-dc(1,i)
        enddo
        do i=1,ndim
           g(i,i) = g(i,i)-siw(iw,i)
        enddo
        g2 = g(:,:)
        call inverse_matrix(g2)
        gloc(iw,:,:) = gloc(iw,:,:)+g2(:,:)
     enddo

  enddo

  do iw=0,iwmax-1
     gloc(-iw-1,:,:) = real(gloc(iw,:,:),kind=8)-ci*aimag(gloc(iw,:,:))
  enddo

  gloc = gloc/dble(nkp)
end subroutine get_gloc



subroutine get_chi0_loc(iwf, iwb, chi0_loc)
  use parameters_module
  implicit none
  integer :: i, j, k, l, i1, i2
  integer :: iwf, iwb
  complex(kind=8), intent(out) :: chi0_loc(ndim*ndim,ndim*ndim)

  chi0_loc = 0.d0
  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1=i1+1
        i2=0
        do l=1,ndim
           do k=1,ndim
              i2=i2+1
              if(i==l .and. j==k)then
                if(index2ineq(nineq,ndims,i,j,k,l)) then
                  chi0_loc(i1,i2) = - beta*giw(iwf,i)*giw(iwf-iwb,j)
                endif
              endif
           enddo
        enddo
     enddo
  enddo

end subroutine get_chi0_loc




subroutine get_chi0_loc_inv(iwf, iwb, chi0_loc)
  use parameters_module
  implicit none
  integer :: i, j, k, l, i1, i2
  integer :: iwf, iwb
  complex(kind=8), intent(out) :: chi0_loc(ndim*ndim,ndim*ndim)

  chi0_loc = 0.d0
  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1=i1+1
        i2=0
        do l=1,ndim
           do k=1,ndim
              i2=i2+1
              if(i==l .and. j==k)then
                if(index2ineq(nineq,ndims,i,j,k,l)) then
                  chi0_loc(i1,i2) = -1.d0/(beta*giw(iwf,i)*giw(iwf-iwb,j))
                endif
              endif
           enddo
        enddo
     enddo
  enddo

end subroutine get_chi0_loc_inv


subroutine get_chi0(ik, ikq, iwf, iwb, chi0)
  use lapack_module
  use parameters_module
  implicit none
  integer :: i, j, k, l, i1, i2
  integer :: iwf, iwb, ik, ikq
  complex(kind=8) :: g1(ndim,ndim), g2(ndim,ndim)
  complex(kind=8), intent(out) :: chi0(ndim*ndim,ndim*ndim)

  g1(:,:) = -hk(:,:,ik)
  do i=1,ndim
     g1(i,i) = ci*iw_data(iwf)+mu-hk(i,i,ik)-dc(1,i)
  enddo
  do i=1,ndim
  g1(i,i) = g1(i,i)-siw(iwf,i)
  enddo


  if (ndim .eq. 1) then
    g1(1,1)=1.d0/g1(1,1)
  else if (ndim .eq. 2) then
    call inverse_matrix_2(g1)
  else if (ndim .eq. 3) then
    call inverse_matrix_3(g1)
  else
    call inverse_matrix(g1)
  end if

  g2(:,:) = -hk(:,:,ikq)
  do i=1,ndim
     g2(i,i) = ci*iw_data(iwf-iwb)+mu-hk(i,i,ikq)-dc(1,i)
  enddo
  do i=1,ndim
  g2(i,i) = g2(i,i)-siw(iwf-iwb,i)
  enddo

  if (ndim .eq. 1) then
    g2(1,1)=1.d0/g2(1,1)
  else if (ndim .eq. 2) then
    call inverse_matrix_2(g2)
  else if (ndim .eq. 3) then
    call inverse_matrix_3(g2)
  else
    call inverse_matrix(g2)
  end if

  chi0 = 0.d0
  i1 = 0
  do i=1,ndim
     do j=1,ndim
        i1=i1+1
        i2=0
        do l=1,ndim
           do k=1,ndim
              i2=i2+1
                if(index2ineq(nineq,ndims,i,j,k,l)) then
                  chi0(i1,i2) = - beta*g1(i,l)*g2(k,j)
                endif
           enddo
        enddo
     enddo
  enddo

end subroutine get_chi0

end module
