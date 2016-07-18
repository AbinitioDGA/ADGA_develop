module parameters_module
  implicit none
  public

  complex(kind=8) ci
  parameter (ci=(0.d0,1.d0))
  double precision,parameter :: pi=3.1415926535897932385d0,t=0.25d0
  integer :: nkp, ndim, ndims, ndim2, maxdim,i1,i2,i3,i4
  integer :: nkpx,nkpy,nkpz,nkp1,nqpx,nqpy,nqpz,nqp1
  integer :: iwmax, iwbmax, iwfmax, iwbmax_small, iwfmax_small,nk_frac
  integer :: iwstart,iwstop
  double precision :: mu, beta
  double precision, allocatable :: k_data(:,:)
  integer :: nqp,nkp_eom
  integer,allocatable :: q_data(:),k_data_eom(:)
  character(len=150) :: filename,filename_vertex,filename_umatrix,output_dir,filename_q_path
  logical :: orb_sym,small_freq_box,full_chi0
  logical :: do_eom,do_chi
  logical :: q_path,q_vol
  integer :: vertex_type
  integer,parameter :: full_g4=0,connected_g4=1,chi_g4=2
  interface k_vector
    module procedure k_vector_1, k_vector_3
  end interface k_vector
  interface k_index
    module procedure k_index_1, k_index_3
  end interface k_index

contains
subroutine read_config()
  implicit none
  character(len=100) :: cmd_arg
  character(len=100) :: config_file
  character(len=150) :: str_tmp
  integer :: int_tmp_1,int_tmp_2,int_tmp_3

  call getarg(1,cmd_arg)
  config_file=trim(cmd_arg)
  write(*,*) 'Reading config: ',config_file

  open(unit=1,file=config_file)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_vertex=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  filename_umatrix=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2,int_tmp_3
  orb_sym=int_tmp_1
  small_freq_box=int_tmp_2
  full_chi0=int_tmp_3

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2
  iwfmax_small=int_tmp_1
  iwbmax_small=int_tmp_2

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1
  if (int_tmp_1.eq.1) then
    q_path=.true.
    q_vol=.false.
  else if (int_tmp_1.eq.2) then 
    q_path=.false.
    q_vol=.true.
  end if

  read(1,*)
  read(1,*)
  read(1,*) str_tmp
  filename_q_path=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) nk_frac

  read(1,*)
  read(1,*)
  read(1,'(A)') str_tmp
  output_dir=trim(str_tmp)

  read(1,*)
  read(1,*)
  read(1,*) int_tmp_1,int_tmp_2
  do_chi=int_tmp_1
  do_eom=int_tmp_2
  close(1)
end subroutine read_config

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
  write(*,*) nsegments, 'q segments'
  close(1)
  n_segments=nsegments
end function n_segments

subroutine init_k_grid_cubic()
  implicit none
  integer :: ik,ikx,iky,ikz,iband
  ik=0
  do ikz=1,nkpz
    do iky=1,nkpy
      do ikx=1,nkpx
        ik=ik+1
        k_data(1,ik)=dble(ikx-1)/dble(nkpx)
        k_data(2,ik)=dble(iky-1)/dble(nkpy)
        k_data(3,ik)=dble(ikz-1)/dble(nkpz)
      end do
    end do
  end do
  open(35, file=trim(output_dir)//"k_data.dat", status='unknown')
  do ik=1,nkp
    write(35,*) k_data(:,ik)
  end do
  close(35)
end subroutine init_k_grid_cubic

subroutine init_hubbard_hamiltonian(hk)
  implicit none
  integer :: ik,iband
  complex(kind=8) :: hk(ndim,ndim,nkp)
  
  hk=0.d0
  do ik=1,nkp
    do iband=1,ndim
      hk(iband,iband,ik)=-2.d0*t*(dcos(2.d0*pi*k_data(1,ik))+dcos(2.d0*pi*k_data(2,ik))+dcos(2.d0*pi*k_data(3,ik)))
    end do
  end do
end subroutine init_hubbard_hamiltonian

subroutine generate_q_vol(qdata)
  implicit none
  integer :: qdata(nqp),i,j,k,i1

  i1=0
  do i=0,nqp1-1
    do j=0,nqp1-1
      do k=0,nqp1-1
        i1 = i1+1
        qdata(i1)=k_index(i*nkpx/nqpx,j*nkpy/nqpy,k*nkpz/nqpz)
      enddo
    enddo
  enddo

end subroutine generate_q_vol

subroutine generate_q_path(qdata)
  implicit none
  integer :: i,nsegments,iostatus,start,seg_len,i1,i2
  character(len=1) :: qpoint_name
  character,allocatable :: qpoints_str(:)
  integer :: qdata(nqp)
!  integer,allocatable :: qpoints(:)

  if (.not. (nkpx.eq.nkpy.and.nkpy.eq.nkpz)) then
    write(*,*) 'Error: condition nkpx==nkpy==nkpz not fulfilled'
    stop
  end if

  write(*,*) 'generating q path from: ',filename_q_path
  nsegments=n_segments()
  write(*,*) nsegments, 'q segments'
  
  open(unit=1,file=filename_q_path)
  allocate(qpoints_str(nsegments+1))
!  allocate(qpoints(nsegments*nkp1/2+1))
  do i=1,nsegments+1
    read(1,*,iostat=iostatus) qpoints_str(i)
  end do
  close(1)

  write(*,*) qpoints_str

  if (nsegments.gt.0) then
  do i=1,nsegments
    if (i.eq.1) then
      seg_len=nkp1/2+1
      start=0
    else
      seg_len=nkp1/2
      start=1
    end if
    i1=(i-1)*nkp1/2+1+start
    i2=i*nkp1/2+1
    write(*,*) i1,i2
    write(*,*)
    
!    call q_path_segment_old(start,i1,i2,qpoints_str(i),qpoints_str(i+1),qdata(i1:i2))
    call q_path_segment(start,i1,i2,qpoints_str(i),qpoints_str(i+1),qdata(i1:i2))
  end do
  else if (nsegments.eq.0) then
    qdata(1)=q_index_from_code(qpoints_str(1))
  else
    write(*,*) 'Error: wrong number of q path segments, probably negative.'
    stop
  end if
!  write(*,*) qdata

end subroutine generate_q_path

function q_index_from_code(q_code)
  implicit none
  character(len=1) :: q_code
  integer :: q_index_from_code

  q_index_from_code=0

  if (q_code.eq.'G') then
    q_index_from_code=k_index(0,0,0)
  else if (q_code.eq.'X') then
    q_index_from_code=k_index(0,0,nkp1/2)
  else if (q_code.eq.'M') then
    q_index_from_code=k_index(0,nkp1/2,nkp1/2)
  else if (q_code.eq.'R') then
    q_index_from_code=k_index(nkp1/2,nkp1/2,nkp1/2)
  else
    write(*,*) 'Error: unknown q code'
    stop
  end if
end function q_index_from_code

subroutine q_path_segment_old(start,istart,istop,qpoint_1,qpoint_2,segment) 
  implicit none
  character :: qpoint_1,qpoint_2
  integer :: start,istart,istop,j
  integer :: segment(istart:istop),i

  if (qpoint_1.eq.'G' .and. qpoint_2.eq.'X') then
    do i=istart,istop
      j=(i-istart+start)*nkp1/nqp1
      segment(i)=k_index(0,0,j)
    end do
  else if (qpoint_1.eq.'X' .and. qpoint_2.eq.'M') then
    do i=istart,istop
      j=(i-istart+start)*nkp1/nqp1
      segment(i)=k_index(0,j,nkp1/2)
    end do
  else if (qpoint_1.eq.'M' .and. qpoint_2.eq.'R') then
    do i=istart,istop
      j=(i-istart+start)*nkp1/nqp1
      segment(i)=k_index(j,nkp1/2,nkp1/2)
    end do
  else if (qpoint_1.eq.'R' .and. qpoint_2.eq.'G') then
    do i=istart,istop
      j=(nkp1/2-i+istart-start)*nkp1/nqp1
      segment(i)=k_index(j,j,j)
    end do
  else
  write(*,*) qpoint_1,qpoint_2,'not yet implemented'
  end if

end subroutine q_path_segment_old
subroutine q_path_segment(start,istart,istop,qpoint_1,qpoint_2,segment) 
  implicit none
  character :: qpoint_1,qpoint_2
  integer :: q_ind_1,q_ind_2,q_vec_1(3),q_vec_2(3),distance(3),step(3)
  integer :: start,istart,istop,j
  integer :: segment(istart:istop),i
  integer :: vector(3)

  q_ind_1=q_index_from_code(qpoint_1)
  q_ind_2=q_index_from_code(qpoint_2)

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
    j=(i-istart+start)*nkp1/nqp1
    vector=q_vec_1+j*step
    segment(i)=k_index(vector)
  end do
end subroutine q_path_segment


subroutine init() 
  implicit none
  maxdim = ndim*ndim*2*iwfmax_small
  ndim2 = ndim*ndim
  if (full_chi0) then
    iwstart=-iwmax+iwbmax
    iwstop=iwmax-iwbmax-1
  else
    iwstart=-iwfmax_small
    iwstop=iwfmax_small-1
  end if
!  nkp=nkpx*nkpy*nkpz

  if (q_vol) then
    nqp = nqpx*nqpy*nqpz
  if (mod(nkpx,nqpx).ne.0 .or. mod(nkpy,nqpy).ne.0 .or. mod(nkpz,nqpz).ne.0) then
     stop 'mismatch between k- and q-grid!'
   endif
  end if
  if (q_path) then
    nkp1=nkpx
    nqp1=nqpx
  end if
end subroutine init

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


subroutine index_kq(ind)
      implicit none
      integer ikp,jkp
      integer :: ind(nkp,nqp)
      ind = 0

! New method

      do ikp=1,nkp
        do jkp=1,nqp
          ind(ikp,jkp)=k_minus_q(ikp,q_data(jkp))
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
  implicit none
  integer :: ik,iq,k_minus_q
  integer :: ix,iy,iz,lx,ly,lz

  call k_vector(ik,ix,iy,iz)
  call k_vector(iq,lx,ly,lz)

!  iz=mod(ik-1,nkpz)
!  lz=mod(iq-1,nqpz)
!  iy=mod((ik-1)/nkpz,nkpy)
!  ly=mod((iq-1)/nqpz,nqpy)
!  ix=(ik-1)/(nkpz*nkpy)
!  lx=(iq-1)/(nqpz*nqpy)

  k_minus_q=1+mod(nkpz+iz-lz,nkpz) + &
              mod(nkpy+iy-ly,nkpy)*nkpx + &
              mod(nkpx+ix-lx,nkpx)*nkpy*nkpx
!  k_minus_q=1+mod(nkpz+iz-lz*nkpz/nqpz,nkpz) + &
!              mod(nkpy+iy-ly*nkpy/nqpy,nkpy)*nkpx + &
!              mod(nkpx+ix-lx*nkpx/nqpx,nkpx)*nkpy*nkpx
end function k_minus_q

function k_index_1(k)
  implicit none
  integer,intent(in) :: k(3)
  integer :: k_index_1

  k_index_1 = 1 + k(3) + k(2)*nkpz + k(1)*nkpy*nkpz

end function k_index_1

function k_index_3(kx,ky,kz)
  implicit none
  integer :: k_index_3,kx,ky,kz

  k_index_3 = 1 + kz + ky*nkpz + kx*nkpy*nkpz

end function k_index_3

subroutine k_vector_1(ik,k)
  implicit none
  integer,intent(in) :: ik
  integer,intent(out) :: k(3)

  k(3)=mod(ik-1,nkpz)
  k(2)=mod((ik-1)/nkpz,nkpy)
  k(1)=(ik-1)/(nkpy*nkpz)
end subroutine k_vector_1

subroutine k_vector_3(ik,kx,ky,kz)
  implicit none
  integer :: ik,kx,ky,kz

  kz=mod(ik-1,nkpz)
  ky=mod((ik-1)/nkpz,nkpy)
  kx=(ik-1)/(nkpy*nkpz)
end subroutine k_vector_3

end module parameters_module
