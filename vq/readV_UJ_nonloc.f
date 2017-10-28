c read W(dR,0,R',R'+dR';n1,n2,n3,n4)
c rcut1: cutoff for R'(ir1)
c rcut2: cutoff for dR(ir2) and dR'(ir3)
      program readV
      use hdf5
      use hdf5_module

      implicit real*8(a-h,o-z)

      real(8),allocatable :: freq(:),freq2(:)
      real(8),allocatable :: rwc(:,:,:,:,:,:),
     &                       rv(:,:,:,:,:), cv(:,:,:,:,:),
     &                       rw(:,:,:,:,:,:),cw(:,:,:,:,:,:),
     &                       rws1(:,:),rws2(:,:),rdum(:,:,:,:),
     &                       rws1_2(:,:),rws2_2(:,:)
      integer(4),allocatable :: irws1(:),irws2(:),irws1_2(:),irws2_2(:)
      character :: ch3*3
      integer i1,i2,i3,i4,ir1,ir2,ir3
c computing U(iv)  [R=0, dens-dens, iw-mesh]
      double precision beta,pi,R,Rp,dR,dRp,RpdRp,veclen
      external veclen
      integer niw,iw,i
      double precision, allocatable :: iom(:)
      complex(kind=8), allocatable :: Wiv(:,:,:)
      complex(kind=8) ci
      parameter (pi=3.1415926535898d0)
      parameter (ci=(0.d0,1.d0))

      
      character(len=13), parameter :: filename="V-nonloc.hdf5"
      integer(hid_t) :: dspace_r_id,r_id,axes_id 
      integer, parameter :: rank_r=2
      character(len=20) :: name_buffer
      integer(hid_t) :: grp_id,dset_id
      integer(hid_t) :: dspace_v_id
      integer(hsize_t), dimension(1) :: dim_v
      integer, parameter :: rank_v=1

      integer ind
      double precision, allocatable :: tmp_v(:)
      

c read V 
      write(*,*)'reading VMATU'
      open(91,file="VMATU")
      read(91,*)
      read(91,*)nwf,alat
      read(91,*)
      read(91,*)rcut1,rcut2
      read(91,*)
      read(91,*)nrws1,nrws2,nrws ! nrws = nrws1*nrws2*nrws2
      allocate(rv(nwf,nwf,nwf,nwf,nrws),cv(nwf,nwf,nwf,nwf,nrws),
     &         rws1(3,nrws1),irws1(nrws1),
     &         rws2(3,nrws2),irws2(nrws2))
      read(91,*)
      read(91,*)rws1,irws1
      read(91,*)
      read(91,*)rws2,irws2
      read(91,*)
      write(*,*)'oh my...'
      read(91,*)rv
      read(91,*)
      write(*,*)'oh my...'
      read(91,*)cv
      write(*,*)'Done read VMATU. puuuh!'
      close(91)

c read Wc
      if (1.eq.1) then
      write(*,*)'reading WcMATU'
      open(92,file="WcMATU")
      read(92,*)
      read(92,*)nwf2,nw,alat2
      read(92,*)
      read(92,*)rcut1_2,rcut2_2
      read(92,*)
      read(92,*)nrws1_2,nrws2_2,nrws_2

      if (nwf .ne. nwf2) stop 'nwf error in WcMAT'
      if (nrws .ne. nrws_2) stop 'nrws error in WcMAT'
      if (nrws1 .ne. nrws1_2) stop 'nrws1 error in WcMAT'
      if (nrws2 .ne. nrws2_2) stop 'nrws2 error in WcMAT'
      if (abs(alat-alat2).gt.1.d-4) stop 'alat error in WcMAT'
      if (abs(rcut1-rcut1_2).gt.1.d-4) stop 'rcut1 error in WcMAT'
      if (abs(rcut2-rcut2_2).gt.1.d-4) stop 'rcut2 error in WcMAT'

      allocate(freq(nw),
     &         rwc(nwf,nwf,nwf,nwf,nrws,nw)
     &        ,rws1_2(3,nrws1),irws1_2(nrws1)
     &        ,rws2_2(3,nrws2),irws2_2(nrws2))
!      allocate(cwc(nwf,nwf,nwf,nwf,nrws,nw))
      read(92,*)
      read(92,*)freq(1:nw)
      read(92,*)
      read(92,*)rws1_2,irws1_2
      read(92,*)
      read(92,*)rws2_2,irws2_2
      read(92,*)
      write(*,*)'oh my...'
      read(92,*)rwc
      read(92,*)
!      read(92,*)cwc
      write(*,*)'done WcMATU. puuh!'
      close(92)

c W = v + Wc
      allocate(rw(nwf,nwf,nwf,nwf,nrws,nw))
      do iw = 1,nw
         rw(:,:,:,:,:,iw) = rwc(:,:,:,:,:,iw) + rv(:,:,:,:,:)
      enddo

      deallocate(rwc)
      endif




     

      call h5open_f(err)

      call h5fcreate_f(filename, h5f_acc_trunc_f, new_file_id, err)

      allocate(tmp_v(nrws1))

      !create dataspaces:
      dims = (/3,nrws1/)
      call h5screate_simple_f(rank_r,dims,dspace_r_id,err)

      dim_v = nrws1
      call h5screate_simple_f(rank_v,dim_v,dspace_v_id,err)

      !create R-list:
      call h5gcreate_f(new_file_id,".axes",axes_id,err)
      call h5dcreate_f(axes_id,"R-points",h5t_native_double,
     & dspace_r_id,r_id,err)
      call h5dwrite_f(r_id,h5t_native_double,rws1,dims,err)
      call h5dclose_f(r_id,err)
      call h5gclose_f(axes_id,err)
 
      
       do i1 = 1,nwf
         do i2 = 1,nwf
           do i3 = 1,nwf
              do i4 = 1,nwf

                !assumes Kanamori interactions 
                if (((i1==i2).and.(i3==i4)) .or. ((i1==i3).and.
     &              (i2==i4)) .or. ((i1==i4).and.(i2==i3))) then
                  
                  tmp_v=0.d0

                  call component2index_band(nwf, ind, i1, i2, i3, i4)

                  write(name_buffer, '(I5.5)'), ind
                  call h5gcreate_f(new_file_id,name_buffer,
     &                 grp_id,err)
                  call h5dcreate_f(grp_id,"value",h5t_native_double,
     &                 dspace_v_id, dset_id, err)


                  ! Ordnung in ctqmc ist U^ijkl c^+i c^+j ck cl
                  !hier hingegen JMiyake=UMiyake^ijji
                  ! also U^ijkl = UMiyake^iljk --> i1 i4 i2 i3
                  tmp_v(:) = rw(i1,i4,i2,i3,:,1)
              
                 
                  call h5dwrite_f(dset_id,h5t_native_double,
     &                         tmp_v,dim_v,err)

                 


                  !the following is for cases where nrws2>1 (when the
                  !screened interaction rw has more than one R-index)
                  !has not yet been used/needs modifications with hdf5-output 


                  !do iw = 1,1              ! loop over frequency ... here we only do STATIC (SVO files have only nw=1 anyways)

                  ! loops over centres of W(dR,0,R',R'+dR')  dR'=dRp etc...
                  !  do ir2=1,nrws2
                  !    dR=veclen(rws2(:,ir2))

                  !    do ir1=1,nrws1
                  !      Rp=veclen(rws1(:,ir1))

                  !      do ir3=1,nrws2
                  !        dRp=veclen(rws2(:,ir3))
                  !        RpdRp=veclen(rws1(:,ir1)+rws2(:,ir3)) !this is R'+dR'

         
                          ! laufender Index fuer Koordinate
                   !       ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1


                   !       if ((dR.lt.1.d-5).and.(dRp.lt.1.d-5)) then ! = two-centre density-density between R=0 and R'=Rp

                                  
                          ! Ordnung in ctqmc ist U^ijkl c^+i c^+j ck cl
                          ! d.h. dichte-dichte ist U^ijji, J=U^ijij

                          !hier hingegen JMiyake=UMiyake^ijji
                          ! also U^ijkl = UMiyake^iljk --> i1 i4 i2 i3
               
                          !write(111,'(1F12.6,4I4,2F12.6)')Rp,i1,i2,i3,i4,
                          !  rw(i1,i4,i2,i3,ir,iw),rv(i1,i4,i2,i3,ir) ! screened interaction, bare interaction

                            

                    !      endif        ! centres
     
                 !       enddo           ! ir3
                 !     enddo              ! ir2
                 !   enddo                 ! ir1
                 ! enddo                    ! nw
                  
                  call h5dclose_f(dset_id,err)   
                  call h5gclose_f(grp_id,err)             
   
                endif !Kanamori
              enddo !orbitals
            enddo
          enddo
        enddo

       
       call h5fclose_f(new_file_id,err)
       call h5close_f(err)
       

       end program readV
c-----------------------------------------------------------------------


      double precision function veclenold(vec)
      implicit none
      double precision vec(3),dum
      integer i
      dum=0.d0
      do i=1,3
         dum=dum+vec(i)*vec(i)
      enddo
      veclenold=dsqrt(dum)
      return
      end

      double precision function veclen(vec)
      implicit none
      double precision vec(3)
      integer i
      veclen=0.d0
      do i=1,3
         veclen=veclen+vec(i)*vec(i)
      enddo
      veclen=dsqrt(veclen)
      return
      end


      subroutine component2index_band(Nbands, ind, b1, b2, b3, b4)
      implicit none
      integer,intent(in) :: Nbands
      integer,intent(in) :: b1, b2, b3, b4
      integer,intent(out) :: ind

      ind =  Nbands**3*(b1-1) + Nbands**2*(b2-1) + Nbands*(b3-1) + b4

      end subroutine component2index_band

