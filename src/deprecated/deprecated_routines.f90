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


subroutine index_kq_search(k_data, q_data, index)
  
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
      if (ounit .gt. 0) write(ounit,*)'using internal dk=',dk


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
            if (ounit .gt. 0) write(ounit,*)ikp,'/',nkp
         endif
        

      enddo !ikp


! consistency check
      do ikp=1,nkp
         do jkp=1,nqp
            if ((index(ikp,jkp).lt.1).or.(index(ikp,jkp).gt.nkp)) then
               if (ounit .gt. 0) write(ounit,*)ikp,jkp,nkp,index(ikp,jkp)
               if (ounit .gt. 0) write(ounit,*)k_data(:,ikp)
               if (ounit .gt. 0) write(ounit,*)k_data(:,jkp)
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

subroutine generate_q_path(n1,qdata,er,erstr)
  implicit none
  integer,intent(in) :: n1
  integer,intent(out) :: qdata(n1**3)
  integer,intent(out) :: er
  character(len=*),intent(out) :: erstr
  integer :: i,nsegments,iostatus,start,seg_len,i1,i2
  character(len=1) :: qpoint_name
  character,allocatable :: qpoints_str(:)
!  integer,allocatable :: qpoints(:)

   er = 0
   erstr = ''
   qdata = 0

  if (.not. (nkpx.eq.nkpy.and.nkpy.eq.nkpz)) then
    erstr='Error: condition nkpx==nkpy==nkpz not fulfilled'
    er = 1
    return
  end if

  if (ounit .gt. 0) write(ounit,*) 'generating q path from: ',filename_q_path
  nsegments=n_segments()
  if (ounit .gt. 0) write(ounit,*) nsegments, 'q segments'
  
  open(unit=1,file=filename_q_path)
  allocate(qpoints_str(nsegments+1))
!  allocate(qpoints(nsegments*nkp1/2+1))
  do i=1,nsegments+1
    read(1,*,iostat=iostatus) qpoints_str(i)
  end do
  close(1)

  if (ounit .gt. 0) write(ounit,*) qpoints_str

  if (nsegments.gt.0) then
  do i=1,nsegments
    if (i.eq.1) then
      start=0
    else
      start=1
    end if
    i1=(i-1)*n1/2+1+start
    i2=i*n1/2+1
    if (ounit .gt. 0) write(ounit,*) i1,i2
    
!    call q_path_segment_old(start,i1,i2,qpoints_str(i),qpoints_str(i+1),qdata(i1:i2))
    call q_path_segment(n1,start,i1,i2,qpoints_str(i),qpoints_str(i+1),qdata(i1:i2),er,erstr)
    if (er  .ne. 0) exit
  end do
  else if (nsegments.eq.0) then
    qdata(1)=q_index_from_code(qpoints_str(1))
  else
    er = 2
    erstr='Error: wrong number of q path segments, probably negative.'
  end if

  return
end subroutine generate_q_path

function q_index_from_code(q_code)
  implicit none
  integer :: q_index_from_code
  character,intent(in) :: q_code

  q_index_from_code=0

  select case(q_code)
  case ('G')
    q_index_from_code=k_index(0,0,0)
  case ('X')
    q_index_from_code=k_index(0,0,nkp1/2)
  case ('M')
    q_index_from_code=k_index(0,nkp1/2,nkp1/2)
  case ('R') 
    q_index_from_code=k_index(nkp1/2,nkp1/2,nkp1/2)
  case default
    ! error
    q_index_from_code = -1
  end select

  return
end function q_index_from_code

subroutine q_path_segment(n1,start,istart,istop,qpoint_1,qpoint_2,segment,er,erstr) 
  implicit none
  integer,intent(out) :: er
  character(len=*),intent(out) :: erstr
  character :: qpoint_1,qpoint_2
  integer :: q_ind_1,q_ind_2,q_vec_1(3),q_vec_2(3),distance(3),step(3)
  integer :: start,istart,istop,j,n1
  integer :: segment(istart:istop),i
  integer :: vector(3)

  er = 0
  erstr = ''

  q_ind_1=q_index_from_code(qpoint_1)
  q_ind_2=q_index_from_code(qpoint_2)

  if (q_ind_1 .le. -1 .or. q_ind_2 .le. -1) then
   er = 1
   erstr = 'q_path_segment: Error: unknown symbol(s): '//qpoint_1//' '//qpoint_2
   return
  endif

  call k_vector(q_ind_1,q_vec_1)

  call k_vector(q_ind_2,q_vec_2)

  distance=q_vec_2-q_vec_1
  step=0
  do i=1,3
    if (distance(i).ne.0) then
      step(i)=distance(i)/abs(distance(i))
    end if
  end do

  do i=istart,istop
    j=(i-istart+start)*nkp1/n1
    vector=q_vec_1+j*step
    segment(i)=k_index(vector)
  end do
end subroutine q_path_segment


function n_segments()
  implicit none
  integer :: i,nsegments,n_segments,iostatus
  character(len=1) :: qpoint_name

  iostatus=0
  nsegments=-2
  open(unit=1,file=filename_q_path)
    
  do while (iostatus.eq.0)
    read(1,*,iostat=iostatus) qpoint_name
    nsegments=nsegments+1
  end do
  close(1)
!  write(*,*) nsegments, 'q segments'
  n_segments=nsegments
end function n_segments


! initialize an orthogonal lattice with 
!   * spacing a and nkpx points in x direction
!   *         b     nkpy           y
!   *         c     nkpz           z
subroutine init_k_grid_cubic(k_data,nkpx,nkpy,nkpz,a,b,c)
  implicit none
  integer :: ik,ikx,iky,ikz,iband,nkpx,nkpy,nkpz
  real*8  :: k_data(:,:),a,b,c
  ik=0
  do ikz=1,nkpz
    do iky=1,nkpy
      do ikx=1,nkpx
        ik=ik+1
        k_data(1,ik)=dble(ikx-1)/dble(nkpx)/a
        k_data(2,ik)=dble(iky-1)/dble(nkpy)/b
        k_data(3,ik)=dble(ikz-1)/dble(nkpz)/c
      end do
    end do
  end do
!  open(35, file=trim(output_dir)//"k_data.dat", status='unknown')
!  do ik=1,nkp
!    write(35,*) k_data(:,ik)
!  end do
!  close(35)
end subroutine init_k_grid_cubic

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
t
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

