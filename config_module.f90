module config_module
  use parameters_module

contains

subroutine read_config()
  use lookup_module

  implicit none
  integer :: i,j,k,l,ineq,stat
  integer :: int_temp
  character(len=150) :: str_temp, str_ineq


  integer :: search_start, search_end
  integer :: subsearch_start, subsearch_end
  character(len=150) :: config_file, output
  integer :: pst, empty

  ! Config File checks
  if (iargc() .ne. 1) then
    stop 'The program has to be executed with exactly one argument. (Name of config file)'
  end if
  call getarg(1,config_file)

  open(unit=10,file=trim(config_file),action='read',iostat=stat)
    if (stat .ne. 0) then
      close(10)
      stop 'Input file cannot be opened'
    endif

    ! line counting
    lines=0
    read_count: do
      read(10,'(A)',END=200) str_temp ! read whole line as string, doesnt skip empty lines
      lines=lines+1
    enddo read_count

    200 continue
    rewind 10

    allocate(file_temp(lines))

    ! remove empty lines and comment strings
    empty=0
    read_temp: do i=1,lines
      read(10,'(A)') str_temp
        str_temp = trim(adjustl(str_temp))
        pst=scan(str_temp,cmnt) ! find out possible comment symbol
        if (pst .eq. 1) then ! whole line is commented
          file_temp(i) = ''
        elseif (pst .eq. 0) then ! no comment symbol found
          file_temp(i) = str_temp
        else  ! getting left side of comment
          file_temp(i) = trim(str_temp(1:pst-1))
        endif

        if (len_trim(file_temp(i)) .eq. 0) then ! filter out empty lines
          empty=empty+1
        endif
    enddo read_temp

    ! rewrite everything to a new clean string array
    allocate(file_save(lines-empty))
    j=1
    read_save: do i=1,lines
      if(len_trim(file_temp(i)) .ne. 0) then
        file_save(j)=file_temp(i)
        j=j+1
      endif
    enddo read_save
    deallocate(file_temp)

    lines=lines-empty

  close(unit=10)


  ! FREE FORMAT LOOKUP
  ! defining default values
  do_chi=.false.
  do_eom=.true.
  q_vol=.true.
  k_path_eom=.false.
  q_path_susc=.false.
  output_dir='Output/'
  filename_hk=''; filename=''; filename_vertex_sym=''
  filename_vq=''; filename_q_path=''; filename_umatrix=''
  nineq=1
  iwfmax_small=-1; iwbmax_small=-1 ! maximum number of frequencies -- see check_freq_range

  ! search for General stuff + Allocation of values
  call group_find('[General]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    stop 'General Group not found'
  endif
  call bool_find('calc-susc', do_chi, search_start, search_end)
  call bool_find('calc-eom', do_eom, search_start, search_end)
  call int_find('NAt', nineq, search_start, search_end)
  call int_find('N4iwf', iwfmax_small, search_start, search_end)
  call int_find('N4iwb', iwbmax_small, search_start, search_end)
  call string_find('HkFile', filename_hk, search_start, search_end)
  if (trim(adjustl(filename_hk)) .eq. '') then
    read_ext_hk = .false.
  else
    read_ext_hk = .true.
  endif
  call string_find('VqFile', filename_vq, search_start, search_end)
  if (trim(adjustl(filename_vq)) .eq. '') then
    do_vq = .false. 
  else
    do_vq = .true.
  endif
  call string_find('QpathFile', filename_q_path, search_start, search_end)
  call int3_find('k-grid', nkpx, nkpy, nkpz, search_start, search_end)
  call int3_find('q-grid', nqpx, nqpy, nqpz, search_start, search_end)
  call bool_find('qvol', q_vol, search_start, search_end)
  call bool_find('k-path-eom', k_path_eom, search_start, search_end)
  call bool_find('q-path-susc', q_path_susc, search_start, search_end)
  call string_find('Output', output_dir, search_start, search_end)
  str_temp = trim(adjustl(output_dir))
  if (scan(trim(str_temp),'/',.true.) .ne. len_trim(str_temp)) then   ! no / at the end
    output_dir = trim(str_temp) // '/'  ! add it
  endif
  call string_find('UFile', filename_umatrix, search_start, search_end)
  if(trim(adjustl(filename_umatrix)) .eq. '') then
    read_ext_u = .false.
  else
    read_ext_u = .true.
  endif


  allocate(interaction(nineq))
  allocate(interaction_mode(nineq))
  allocate(ndims(nineq,2))
  allocate(Udd(nineq),Vdd(nineq),Jdd(nineq))
  allocate(Udp(nineq),Vdp(nineq),Jdp(nineq))
  allocate(Upp(nineq),Vpp(nineq),Jpp(nineq))
  interaction=''
  interaction_mode=0; ndims=0
  Udd=0.d0; Vdd=0.d0; Jdd=0.d0
  Udp=0.d0; Vdp=0.d0; Jdp=0.d0
  Upp=0.d0; Vpp=0.d0; Jpp=0.d0


  ! search for Atoms (interaction parameters for umatrix)
  call group_find('[Atoms]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    stop 'Atoms Group not found'
  endif
  do ineq=1,nineq
    write(str_ineq,'(A2,I1,A2)') '[[',ineq,']]'
    call subgroup_find(str_ineq, search_start, search_end, subsearch_start, subsearch_end)
    if (subsearch_start .eq. 0) then ! group was not found
      stop 'Atomnumber subgroup not found'
    endif

    call string_find('Interaction',interaction(ineq),subsearch_start,subsearch_end)
    call int_find('Nd',ndims(ineq,1),subsearch_start,subsearch_end)
    call int_find('Np',ndims(ineq,2),subsearch_start,subsearch_end)
    call float_find('Udd',Udd(ineq),subsearch_start,subsearch_end)
    call float_find('Vdd',Vdd(ineq),subsearch_start,subsearch_end)
    call float_find('Jdd',Jdd(ineq),subsearch_start,subsearch_end)
    call float_find('Upp',Upp(ineq),subsearch_start,subsearch_end)
    call float_find('Vpp',Vpp(ineq),subsearch_start,subsearch_end)
    call float_find('Jpp',Jpp(ineq),subsearch_start,subsearch_end)
    call float_find('Udp',Udp(ineq),subsearch_start,subsearch_end)
    call float_find('Vdp',Vdp(ineq),subsearch_start,subsearch_end)
    call float_find('Jdp',Jdp(ineq),subsearch_start,subsearch_end)

    select case (interaction(ineq))
      case ('Kanamori')
        interaction_mode(ineq) = 1
      case default
        interaction_mode(ineq) = 0
    end select
  enddo


  ndim=sum(ndims)

  if (ndim .eq. 0) then
    stop 'Number of bands per atom is required in [Atoms] section'
  endif

  allocate(u(ndim**2,ndim**2), u_tilde(ndim**2,ndim**2))
  deallocate(interaction)


  ! search for 1particle and 2particle files / parameters
  call group_find('[One-Particle]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    stop 'One-Particle Group not found'
  endif
  call string_find('1PFile', filename, search_start, search_end)
  call bool_find('orb-sym', orb_sym, search_start, search_end)

  call group_find('[Two-Particle]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    stop 'Two-Particle Group not found'
  endif
  call string_find('2PFile', filename_vertex_sym, search_start, search_end)
  call int_find('vertex-type', vertex_type, search_start, search_end)

end subroutine read_config


subroutine init()
  implicit none
  integer :: i
  maxdim = ndim*ndim*2*iwfmax_small
  ndim2 = ndim*ndim
  if (full_chi0 .and. do_chi) then
    iwstart=min(-iwmax+iwbmax,-iwfmax_small)
    iwstop=max(iwmax-iwbmax-1,iwfmax_small-1)
  else
    iwstart=-iwfmax_small
    iwstop=iwfmax_small-1
  end if

  ! Since currently the BZ sum has to go over all points of the Hamiltonian, we
  ! do a consistency check here.
  if (nkpx*nkpy*nkpz .ne. nkp) then
    stop 'Wrong number of k points in config file.'
  end if
  nkp=nkpx*nkpy*nkpz


  if (q_vol) then
    nqp = nqpx*nqpy*nqpz
    if (mod(nkpx,nqpx).ne.0 .or. mod(nkpy,nqpy).ne.0 .or. mod(nkpz,nqpz).ne.0) then
      stop 'mismatch between k- and q-grid!'
    endif
  end if

  if (q_path_susc .or. k_path_eom) then
    stop 'q paths currently not stable'
    nkp1=nkpx
    nqp1=nqpx
  end if

  ! create arrays with Matsubara frequencies
  allocate(iw_data(-iwmax:iwmax-1),iwb_data(-iwbmax:iwbmax),iwf_data(-iwfmax:iwfmax-1))
  do i=-iwmax,iwmax-1
    iw_data(i)=pi*(2*i+1)/beta
  end do
  do i=-iwfmax,iwfmax-1
    iwf_data(i)=pi*(2*i+1)/beta
  end do
  do i=-iwbmax,iwbmax
    iwb_data(i)=pi*2*i/beta
  end do
  
end subroutine init

subroutine finalize()
  implicit none
  deallocate(iw_data,iwf_data,iwb_data)
end subroutine finalize


subroutine check_freq_range(mpi_wrank,master)
  implicit none
  integer :: mpi_wrank, master

  if (iwfmax_small .le. 0) then
    iwfmax_small = iwfmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Calculating with maximum number of fermionic frequencies &
                  =', iwfmax
    endif
  endif

  if (iwbmax_small .lt. 0) then
    iwbmax_small = iwbmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Calculating with maximum number of bosonic frequencies &
                  =', iwbmax
    endif
  endif

  if (iwfmax_small .gt. iwfmax) then
    iwfmax_small = iwfmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Error: Wrong input for fermionic frequencies'
      write(*,*) 'Calculating with maximum number of fermionic frequencies &
                  =', iwfmax
    endif
  endif

  if (iwbmax_small .gt. iwbmax) then
    iwbmax_small = iwbmax
    if (mpi_wrank .eq. master) then
      write(*,*) 'Error: Wrong input for bosonic frequencies'
      write(*,*) 'Calculating with maximum number of bosonic frequencies &
                  =', iwbmax
    endif
  endif

end subroutine check_freq_range


subroutine check_config()
  implicit none
  logical :: there
  integer :: ineq

  exist_p = .false.
  do ineq=1,nineq
    if(ndims(ineq,2) .ne. 0) exist_p=.true.
  enddo

  if (read_ext_hk .eqv. .true.) then
    inquire (file=trim(filename_hk), exist=there)
    if (.not. there) then
      stop "Error: Hamiltonian file does not exist"
    endif
  endif

  if (read_ext_u .eqv. .true.) then
    inquire (file=trim(filename_umatrix), exist=there)
    if (.not. there) then
      stop "Error: Umatrix file does not exist"
    endif
  endif

  if (q_vol .eqv. .false.) then
    inquire (file=trim(filename_q_path), exist=there)
    if (.not. there) then
      stop "Error: Q-Path file does not exist"
    endif
  endif

  if (do_vq .eqv. .true.) then
    inquire (file=trim(filename_vq), exist=there)
    if (.not. there) then
      stop "Error: V(Q) file does not exist"
    endif
  endif

  inquire (file=trim(filename), exist=there)
  if (.not. there) then
    stop "Error: One-Particle data file does not exist"
  endif

  ! inquire (file=trim(filename_umatrix), exist=there)
  ! if (.not. there) then
  !   stop "Error: Umatrix file does not exist"
  ! endif

  if (vertex_type .lt. 0 .or. vertex_type .gt. 2) then
    stop "Error: Choose appropriate vertex type"
  endif

  inquire (file=trim(filename_vertex_sym), exist=there)
  if (.not. there) then
    stop "Error: Two-Particle data file does not exist"
  endif

  end subroutine check_config

end module config_module
