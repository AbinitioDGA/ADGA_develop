module index_module
  implicit none
  public

contains

subroutine index_kq_search(k_data, q_data, index)
  
      use parameters_module
      implicit none
      integer ikp,jkp,i,n,icheck
      integer :: index(nkp,nqp)
      double precision kq(3),dum,dk
      double precision :: k_data(3,nkp), q_data(3,nqp)

      index = 0

      !determine distance between the first 2 k-points as a measure of how fine the k-mesh is: dk=|k1-k2|
      dum = 0.d0
      do i=1,3
         dum = dum+(k_data(i,1)-k_data(i,2))**2.d0
      enddo
      dk = dsqrt(dum)
      write(*,*)'using internal dk=',dk


      do ikp=1,nkp  !k
         do jkp=1,nqp  !q
            kq(:) = k_data(:,ikp)-q_data(:,jkp)  !k+q

!fold back into 1.BZ (k-points defined on a grid between 0(=0) and 1(=2\pi)
            do i=1,3
               if (kq(i).lt.0.d0) then
                  kq(i) = kq(i)+1.d0
               endif
            enddo
            

!search for its number
            i=0
            icheck=0
            do while ((i.lt.nkp).and.(icheck.eq.0))
               i=i+1
               !checks if kq and bk(:,i) are less than dk*0.9 apart
               call checkk(kq,k_data(:,i),icheck,dk*0.9d0)
            enddo
            index(ikp,jkp)=i

            !write(*,*)'ikp=', ikp, 'jkp=', jkp, 'kqind=', i
         enddo !jkp

         
         if (abs(mod(real(ikp),real(500))).lt.0.1d0) then
            write(*,*)ikp,'/',nkp
         endif
        

      enddo !ikp


! consistency check
      do ikp=1,nkp
         do jkp=1,nqp
            if ((index(ikp,jkp).lt.1).or.(index(ikp,jkp).gt.nkp)) then
               write(*,*)ikp,jkp,nkp,index(ikp,jkp)
               write(*,*)k_data(:,ikp)
               write(*,*)k_data(:,jkp)
               STOP 'k out of bound'
            endif
!! Consistency check of new integer-based method
!            if (index(ikp,jkp).ne.k_minus_q(ikp,jkp)) then
!              write(*,*) 'k-q mismatch at ik',ikp,'iq',jkp
!              write(*,*) 'k-q',k_minus_q(ikp,jkp),'correct',index(ikp,jkp)
!              write(*,*) ikp,k_minus_q(ikp,jkp)
!              STOP 
!            end if
!            if (index(ikp,jkp) .eq. k_minus_q(ikp,jkp)) then
!              cnt=cnt+1
!            end if
         enddo
      enddo
!      write(*,*) cnt, 'correct additions'

! write index of k+q into file
!      open(11,file='HMLT.index.kpq',status='unknown')
!      do ikp=1,nkp
!         write(11,*)(index(ikp,jkp),jkp=1,nqp)
!      enddo
!      close(11)

      !do ikp=1,nkp
      !   do jkp=1,nqp
      !      write(*,*)ikp,k_data(:,ikp)
      !      write(*,*)jkp,q_data(:,jkp)
      !      write(*,*)index(ikp,jkp),k_data(:,index(ikp,jkp))
      !      write(*,*)
      !   enddo
      !enddo

end subroutine index_kq_search


subroutine index_kq(index)
  
      use parameters_module
      implicit none
      integer ikp,jkp,i,n,icheck,cnt
      integer :: index(nkp,nqp)
      double precision kq(3),dum,dk

      index = 0

! New method

      do ikp=1,nkp
        do jkp=1,nqp
          index(ikp,jkp)=k_minus_q(ikp,jkp)
        end do
      end do

end subroutine index_kq 


subroutine checkk(k,q,ic,dk)
      implicit none
      integer ic
      double precision k(3),q(3),p(3),op,dk

      op=real(dk,kind=8)
      ic=0
      p=k-q
!      write(*,*) p
      if ( (abs(p(1)).lt.op).and.(abs(p(2)).lt.op).and.(abs(p(3)).lt.op) ) then
        ic =1
!        write(*,*) p
      end if
return
end subroutine checkk


! The following function calculates the index of \vec{k} - \vec{q}.
! It uses only integers
! \vec{k} is associated to (ix,iy,iz)
! \vec{q} is associated to (lx,ly,lz)
! k-space is assumed to have nkpx*nkpy*nkpz points
! q-space is assumed to have nqpx*nqpy*nqpz points,
! where each element of the q-space has to be an element of the k-space.
! subtractions are done in integers, 
! fold-back to BZ is achieved by modulo division.
function k_minus_q(ik,iq)
  use parameters_module
  implicit none
  integer :: ik,iq,k_minus_q
  integer :: ix,iy,iz,lx,ly,lz

  ix=mod(ik-1,nkpx)
  lx=mod(iq-1,nqpx)
  iy=mod((ik-1)/nkpx,nkpy)
  ly=mod((iq-1)/nqpx,nqpy)
  iz=(ik-1)/(nkpx*nkpy)
  lz=(iq-1)/(nqpx*nqpy)

  k_minus_q=1+mod(nkpx+ix-lx*nkpx/nqpx,nkpx) + &
              mod(nkpy+iy-ly*nkpy/nqpy,nkpy)*nkpz + &
              mod(nkpz+iz-lz*nkpz/nqpz,nkpz)*nkpy*nkpz
end function k_minus_q

end module index_module
