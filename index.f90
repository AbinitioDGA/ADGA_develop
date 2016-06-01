subroutine index_kq(k_data, q_data, index)
  
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

         enddo
      enddo

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

end subroutine index_kq 


subroutine checkk(k,q,ic,dk)
      implicit none
      integer ic
      double precision k(3),q(3),p(3),op,dk

      op=real(dk,kind=8)
      ic=0
      p=k-q
      if ( (abs(p(1)).lt.op).and.(abs(p(2)).lt.op).and.(abs(p(3)).lt.op) ) ic=1

return
end subroutine checkk


