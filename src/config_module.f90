module config_module
  use parameters_module

contains

subroutine read_config(er,erstr)
  use lookup_module

  implicit none
  integer,intent(out) :: er
  character(len=*),intent(out) :: erstr
  integer :: i,j,k,l,ineq,stat
  integer :: int_temp
  character(len=150) :: str_temp, str_ineq

  integer :: search_start, search_end
  integer :: subsearch_start, subsearch_end
  character(len=150) :: config_file, output
  integer :: pst, empty

  character(len=50), allocatable  :: general_dict(:)
  character(len=50), allocatable  :: atom_dict(:)
  character(len=50), allocatable  :: oneparticle_dict(:)
  character(len=50), allocatable  :: twoparticle_dict(:)
  character(len=50), allocatable  :: output_dict(:)
! variables for date-time string
  character(20) :: date,time,zone
  integer,dimension(8) :: time_date_values

  er = 0
  erstr = ''

  ! Config File checks
  if (iargc() .ne. 1) then
    er = 1
    erstr = 'The program has to be executed with exactly one argument. (Name of config file)'
    return
  end if
  call getarg(1,config_file)

  open(unit=10,file=trim(config_file),action='read',iostat=stat)
    if (stat .ne. 0) then
      close(10)
      er = 2
      erstr = 'Input file cannot be opened'
      return
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


  !=================================================================================
  ! FREE FORMAT LOOKUP
  !================================================================================
  ! defining sensible default values
  do_chi=.true.                       ! chi-calculation on
  do_eom=.true.                       ! eom-calculation on
  do_cond=.false.                     ! conductivity calculation off
  cond_dmftlegs = .true.              ! calculate conductivity with legs from DMFT
  do_cond_ph = .false.                ! calculate particle-hole vertex contribution of the conductivity
  do_cond_phbar = .true.              ! calculate particle hole bar vertex contribution of the conductivity
  extend_cond_bubble = .false.        ! extend frequency box of conductivity bubble
  q_vol=.true.                        ! homogeneous q-volume on
  q_path_susc=.false.                 ! q-path disabled
  k_path_eom=.false.                  ! k-path disabled
  external_chi_loc=.false.            ! no external local chi
  external_threelegs=.false.          ! no external gamma^wv

  susc_full_output=.false.            ! 2 leg dependencies instead of all 4
  text_output=.false.                 ! text files disabled
  gzip_compression=4                  ! gzip compression of large hdf5 dataset - ranges from 0 to 9
  output_dir='output/'                ! default output folder that gets created

  filename_hk=''; filename_1p=''; filename_vertex_sym=''
  filename_vq=''; filename_qdata=''; filename_umatrix=''
  filename_kdata=''
  filename_chi_loc=''; filename_threelegs=''
  filename_hkder=''
  filename_condlegs = ''
  output_filename=''


  dmft_iter='dmft-last'
  summation_order = -1

  nineq=1
  orb_sym = .false.
  vertex_type = -1
  iwfmax_small=-1 ! default -> calculate in the full vertex frequency box
  iwbmax_small=-1
  iwbcond=0
  nkpx = 0; nkpy = 0; nkpz = 0
  nqpx = 0; nqpy = 0; nqpz = 0

  calc_eigen=.false.
  number_eigenvalues = 1
  number_eigenvectors = 1

  bse_inversion = .true.
  sc_mode = .false.
  !================================================================================


  ! search for General stuff + Allocation of values
  !--------------------------------------------------------------------------------
  call group_find('[General]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    er = 3
    erstr = 'General Group not found'
    return
  endif
  if (search_start .eq. -1) then ! group was not found
    er = 4
    erstr = 'General Group empty'
    return
  endif

  allocate(general_dict(21))
  ! defining dictionary (filling of general_dict)
  general_dict(1)  = 'calc-susc'
  general_dict(2)  = 'calc-eom'
  general_dict(3)  = 'NAt'
  general_dict(4)  = 'N4iwf'
  general_dict(5)  = 'N4iwb'
  general_dict(6)  = 'HkFile'
  general_dict(7)  = 'VqFile'
  general_dict(8)  = 'QDataFile'
  general_dict(9)  = 'KDataFile'
  general_dict(10) = 'k-grid'
  general_dict(11) = 'q-grid'
  general_dict(12) = 'Output'
  general_dict(13) = 'Outfile'
  general_dict(14) = 'UFile'
  general_dict(15) = 'calc-cond'
  general_dict(16) = 'N1iwbc'
  general_dict(17) = 'HkdkFile'
  general_dict(18) = 'cond-legs'
  general_dict(19) = 'extend-cond-bubble'
  general_dict(20) = 'cond-ph'
  general_dict(21) = 'cond-phbar'
  ! spell checking for General group

  call spell_check(search_start, search_end, 'General', general_dict, er, erstr)
  if (er .ne. 0) return
  deallocate(general_dict)

  ! read values
  call bool_find('calc-susc', do_chi, search_start, search_end)
  call bool_find('calc-eom', do_eom, search_start, search_end)
  call bool_find('calc-cond', do_cond, search_start, search_end)
  call int_find('NAt', nineq, search_start, search_end)
  call int_find('N4iwf', iwfmax_small, search_start, search_end)
  call int_find('N4iwb', iwbmax_small, search_start, search_end)
  call int_find('N1iwbc', iwbcond, search_start, search_end)
  call string_find('HkFile', filename_hk, search_start, search_end)
  if (trim(adjustl(filename_hk)) .eq. '') then
    read_ext_hk = .false.
  else
    read_ext_hk = .true.
  endif
  call string_find('HkdkFile', filename_hkder, search_start, search_end)
  if (do_cond .and. trim(adjustl(filename_hkder)) .eq. '') then
    er = 15; erstr='do_cond requires HkdkFile with band derivatives'
    return
  endif
  call string_find('cond-legs', filename_condlegs, search_start, search_end)
  if (trim(adjustl(filename_condlegs)) .eq. '') then ! not found or empty
    cond_dmftlegs = .true.
  else
    cond_dmftlegs = .false. ! use the provided Greens function for the outer legs
  endif

  call bool_find('cond-ph', do_cond_ph, search_start, search_end)
  call bool_find('cond-phbar', do_cond_phbar, search_start, search_end)

  call bool_find('extend-cond-bubble', extend_cond_bubble, search_start, search_end)
  call string_find('VqFile', filename_vq, search_start, search_end)
  if (trim(adjustl(filename_vq)) .eq. '') then
    do_vq = .false.
  else
    do_vq = .true.
  endif
  call string_find('QDataFile', filename_qdata, search_start, search_end)
  if (trim(adjustl(filename_qdata)) .eq. '') then
    q_path_susc = .false.
    q_vol = .true.
  else
    q_path_susc = .true.
    q_vol = .false.
  end if
  call string_find('KDataFile', filename_kdata, search_start, search_end)
  if (trim(adjustl(filename_kdata)) .eq. '') then
    k_path_eom = .false.
  else
    k_path_eom = .true.
  end if
  call int3_find('k-grid', nkpx, nkpy, nkpz, search_start, search_end)
  call int3_find('q-grid', nqpx, nqpy, nqpz, search_start, search_end)
  call string_find('Output', output_dir, search_start, search_end)
  if (len_trim(adjustl(output_dir)) .ge. 1) then
   str_temp = trim(adjustl(output_dir))
   if (scan(trim(str_temp),'/',.true.) .ne. len_trim(str_temp)) then   ! no / at the end
     output_dir = trim(str_temp) // '/'  ! add it
   endif
  else
   output_dir='output/'
  endif

  ! If the field 'Outfile' is present in the config file, use the specified
  ! filename. Otherwise restore the old behaviour (date and time string)
  call string_find('Outfile', output_filename, search_start, search_end)
  if (len_trim(adjustl(output_filename)) .ge. 1) then
    write(*,*) 'using user-specified output name'
   output_filename = trim(output_dir)//trim(adjustl(output_filename))
  else
   write(*,*) 'using default output name' 
   call date_and_time(date,time,zone,time_date_values)
   output_filename=trim(output_dir)//'adga-'//trim(date)//'-'//trim(time)//'-output.hdf5'
   write(*,*) output_filename
  endif
  call string_find('UFile', filename_umatrix, search_start, search_end)
  if(trim(adjustl(filename_umatrix)) .eq. '') then
    read_ext_u = .false.
  else
    read_ext_u = .true.
  endif

  !--------------------------------------------------------------------------------
  verbose = .false.
  verbstr = ''
  call group_find('[Verbose]', search_start, search_end)
  if (search_start .ge. 1) then
     verbose = .true.
     verbstr = file_save(search_start)
  endif

  !--------------------------------------------------------------------------------
  debug = .false.
  dbgstr = ''
  call group_find('[Debug]', search_start, search_end)
  if (search_start .ge. 1) then
     debug = .true.
     dbgstr = file_save(search_start)
  endif

  ! Make sure that we only use 1 q-point when we run with Onlydmft
  if (debug .and. (index(dbgstr,"Onlydmft") .ne. 0)) then
     ! Only local quantities and overwrite everything else
     nqpx = 1
     nqpy = 1
     nqpz = 1
     q_vol = .true.
     q_path_susc =.false.
     k_path_eom = .false.
  endif

  if (debug .and. (index(dbgstr,"Onlyph") .ne. 0)) then
    do_ph = .true.
  else
    do_ph = .false.
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
  !--------------------------------------------------------------------------------
  call group_find('[Atoms]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    er = 5
    erstr = 'Atoms Group not found'
    return
  endif
  if (search_start .eq. -1) then ! group empty
    er = 6
    erstr = 'Atoms Group empty'
    return
  endif

  allocate(atom_dict(12))
  ! defining dictionary (filling of atom_dict)
  atom_dict(1)  = 'Interaction'
  atom_dict(2)  = 'Nd'
  atom_dict(3)  = 'Np'
  atom_dict(4)  = 'Udd'
  atom_dict(5)  = 'Vdd'
  atom_dict(6)  = 'Jdd'
  atom_dict(7)  = 'Upp'
  atom_dict(8)  = 'Vpp'
  atom_dict(9)  = 'Jpp'
  atom_dict(10) = 'Udp'
  atom_dict(11) = 'Vdp'
  atom_dict(12) = 'Jdp'

  do ineq=1,nineq
    write(str_ineq,'(A2,I1,A2)') '[[',ineq,']]'
    call subgroup_find(str_ineq, search_start, search_end, subsearch_start, subsearch_end)
    if (subsearch_start .eq. 0) then ! group was not found
      er = 7
      write(erstr,'("Atomnumber ", I1," subgroup not found")') ineq
      return
    endif
    if (subsearch_start .eq. -1) then ! group empty
      er = 8
      write(erstr,'("Atomnumber ", I1," subgroup empty")') ineq
      return
    endif
    if ((ineq .eq. nineq) .and. (subsearch_end .ne. search_end)) then
      er = 9
      erstr = 'More Atom descriptions than provided in NAt'
      return
    endif

    ! spell checking for General group
    call spell_check(subsearch_start, subsearch_end, 'Atom', atom_dict, er, erstr)
    if (er .ne. 0) return

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
      case ('kanamori') ! because fortran, that's why
        interaction_mode(ineq) = 1
      case default
        interaction_mode(ineq) = 0
    end select
  enddo
  deallocate(atom_dict)


  ndim=sum(ndims)

  if (ndim .eq. 0) then
    er = 10
    erstr = 'Number of bands per atom is required in [Atoms] section'
    return
  endif

  allocate(u(ndim**2,ndim**2), u_tilde(ndim**2,ndim**2))
  deallocate(interaction)


  ! search for 1particle and 2particle files / parameters
  !--------------------------------------------------------------------------------
  call group_find('[One-Particle]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    er = 11
    erstr = 'One-Particle Group not found'
    return
  endif
  if (search_start .eq. -1) then ! group was not found
    er = 12
    erstr = 'One-Particle Group empty'
    return
  endif

  allocate(oneparticle_dict(3))
  oneparticle_dict(1) = '1PFile'
  oneparticle_dict(2) = 'dmft-iter'
  oneparticle_dict(3) = 'orb-sym'
  call spell_check(search_start, search_end, 'One-particle', oneparticle_dict, er, erstr)
  if (er .ne. 0) return
  deallocate(oneparticle_dict)

  call string_find('1PFile', filename_1p, search_start, search_end)
  call string_find('dmft-iter', dmft_iter, search_start, search_end)
  call bool_find('orb-sym', orb_sym, search_start, search_end)

  !--------------------------------------------------------------------------------
  call group_find('[Two-Particle]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    er = 13
    erstr= 'Two-Particle Group not found'
    return
  endif
  if (search_start .eq. -1) then ! group was not found
    er = 14
    erstr= 'Two-Particle Group empty'
    return
  endif
  ! defining dictionary (filling of general_dict)
  allocate(twoparticle_dict(4))
  twoparticle_dict(1) = '2PFile'
  twoparticle_dict(2) = 'vertex-type'
  twoparticle_dict(3) = 'chi-loc-file'
  twoparticle_dict(4) = 'threeleg-file'
  call spell_check(search_start, search_end, 'Two-particle', twoparticle_dict, er, erstr)
  if (er .ne. 0) return
  deallocate(twoparticle_dict)

  call string_find('2PFile', filename_vertex_sym, search_start, search_end)
  call int_find('vertex-type', vertex_type, search_start, search_end)

  call string_find('chi-loc-file',filename_chi_loc, search_start,search_end)
  if (trim(adjustl(filename_chi_loc)) .eq. '') then
    external_chi_loc = .false.
  else
    external_chi_loc = .true.
  endif

  call string_find('threeleg-file',filename_threelegs,search_start,search_end)
  if (trim(adjustl(filename_threelegs)) .eq. '') then
    external_threelegs = .false.
  else
    external_threelegs = .true.
  endif

  !--------------------------------------------------------------------------------
  call group_find('[Output]', search_start, search_end)
  if (search_start .gt. 0) then ! group was found -- this is an optional group
    ! defining dictionary (filling of general_dict)
    allocate(output_dict(3))
    output_dict(1) = 'susc-full-output'
    output_dict(2) = 'gzip-compression'
    output_dict(3) = 'text-output'
    call spell_check(search_start, search_end, 'Output', output_dict, er, erstr)
    if (er .ne. 0) return
    deallocate(output_dict)

    call bool_find('susc-full-output', susc_full_output, search_start, search_end)
    call int_find('gzip-compression', gzip_compression, search_start, search_end)
    call bool_find('text-output', text_output, search_start, search_end)
  endif

  call group_find('[Selfconsistency]',search_start,search_end)
  if (search_start .gt. 0) then 
    call int_find('summation-order',summation_order,search_start,search_end)
    if (summation_order .lt. 0) then
      bse_inversion = .true.
    else
      bse_inversion = .false. ! default value: true
    end if

    ! this can be omitted if it turns out that k-dependent selfenergy is always needed for SC
    call bool_find('se-nonloc',sc_mode,search_start,search_end)

  endif

  call group_find('[Eigenvalues]', search_start, search_end)
  if (search_start .gt. 0) then ! group was found -- this is an optional group
    calc_eigen = .true.
    call int_find('Nvalues', number_eigenvalues, search_start, search_end)
    call int_find('Nvectors', number_eigenvectors, search_start, search_end)
  else
    calc_eigen = .false.
  endif

  ! call group_find('[Conductivity]', search_start, search_end)
  ! if (search_start .gt. 0) then ! group was found -- this is an optional group
  !   do_cond = .true.
  !   call string_find('matrix-elements', filename_hkder, search_start, search_end)
  !   if (trim(adjustl(filename_hkder)) .eq. '') then ! not found or empty
  !     er = 15
  !     erstr = 'matrix-elements have to be provided for a conductivity calculation'
  !     return
  !   endif

  !   call int_find('Nf', iwbcond, search_start, search_end)
  !   if (iwbcond < 0) then ! not found or empty
  !     er = 16
  !     erstr = 'number of bosonic frequencies for conductivity must be equal or greater to 0'
  !     return
  !   endif

  !   call string_find('cond-legs', filename_condlegs, search_start, search_end)
  !   if (trim(adjustl(filename_condlegs)) .eq. '') then ! not found or empty
  !     cond_dmftlegs = .true.
  !   else
  !     cond_dmftlegs = .false. ! use the provided Greens function for the outer legs
  !   endif

  !   call bool_find('extend-bubble', extend_cond_bubble, search_start, search_end)

  ! else
  !   do_cond = .false.
  ! endif


  deallocate(file_save)
  return
end subroutine read_config


subroutine config_init(er,erstr)
  implicit none
  integer,intent(out) :: er
  character(len=*),intent(out) :: erstr
  integer :: i
  real(KIND=8),parameter :: pi = 4d0*atan(1d0)
  logical :: chi0flag

  er = 0
  erstr = ''
  ! Check if we should use the big range for chi0^w and chi0^q
  chi0flag = (debug .and. (index(dbgstr,"Bubble") .ne. 0))

  maxdim = ndim*ndim*2*iwfmax_small
  ndim2 = ndim*ndim
  if (chi0flag .and. do_chi) then
    iwstart=-iwmax+iwbmax_small ! we have a config check for the previous min statement now
    iwstop=iwmax-iwbmax_small-1 ! same here for the previous max statement
  else
    iwstart=-iwfmax_small
    iwstop=iwfmax_small-1
  end if

  ! Since currently the BZ sum has to go over all points of the Hamiltonian, we
  ! do a consistency check here.
  if (nkpx*nkpy*nkpz .ne. nkp) then
    er = 1
    erstr = 'Wrong number of k points in config file.'
    return
  end if
  nkp=nkpx*nkpy*nkpz


  if (q_vol) then
    nqp = nqpx*nqpy*nqpz
    if ((nqpx .le. 0) .or. (nqpy .le. 0) .or. (nqpz .le. 0)) then
      er = 2
      erstr = 'invalid q-grid'
      return
    endif
    if (mod(nkpx,nqpx).ne.0 .or. mod(nkpy,nqpy).ne.0 .or. mod(nkpz,nqpz).ne.0) then
      er = 3
      erstr = 'mismatch between k- and q-grid!'
      return
    endif
  end if

  if (calc_eigen) then
    if (number_eigenvalues .lt. 0) then
      number_eigenvalues = maxdim ! all of them
    else if (number_eigenvalues .eq. 0) then
      calc_eigen = .false. ! deactivate it again .. reasoning 0 EVs -> no diagonalization
    endif
    if (number_eigenvectors .lt. 0) then
      number_eigenvectors = maxdim ! all of them
    endif
  endif

  ! create arrays with Matsubara frequencies
  allocate(iw_data(-iwmax:iwmax-1),iwb_data(-iwbmax_small:iwbmax_small),iwf_data(-iwfmax_small:iwfmax_small-1))
  do i=-iwmax,iwmax-1
    iw_data(i)=pi*(2*i+1)/beta
  end do
  do i=-iwfmax_small,iwfmax_small-1
    iwf_data(i)=pi*(2*i+1)/beta
  end do
  do i=-iwbmax_small,iwbmax_small
    iwb_data(i)=pi*2*i/beta
  end do
  
end subroutine config_init

subroutine finalize()
  implicit none
  deallocate(iw_data,iwf_data,iwb_data)
end subroutine finalize


subroutine check_freq_range(er)
  implicit none
  integer :: mpi_wrank, master, er

  if (ounit .gt. 0) write(ounit,'(1x)')

  if (iwfmax_small .le. 0) then
    iwfmax_small = iwfmax
    if (ounit .gt. 0) then
      write(ounit,*) 'Calculating with maximum number of fermionic frequencies &
                  =', iwfmax_small
    endif
  endif

  if (iwbmax_small .lt. 0) then
    iwbmax_small = iwbmax
    if (ounit .gt. 0) then
      write(ounit,*) 'Calculating with maximum number of bosonic frequencies &
                  =', iwbmax_small
    endif
  endif

  if (iwfmax_small .gt. iwfmax) then
    iwfmax_small = iwfmax
    if (ounit .gt. 0) then
      write(ounit,*) 'Error: Wrong input for fermionic frequencies'
      write(ounit,*) 'Calculating with maximum number of fermionic frequencies &
                  =', iwfmax_small
    endif
  endif

  if (iwbmax_small .gt. iwbmax) then
    iwbmax_small = iwbmax
    if (ounit .gt. 0) then
      write(ounit,*) 'Error: Wrong input for bosonic frequencies'
      write(ounit,*) 'Calculating with maximum number of bosonic frequencies &
                  =', iwbmax_small
    endif
  endif

  if (n3iwb .lt. iwbmax_small .and. external_threelegs) then
    er = 1
    if (ounit .gt. 0) then
      write(ounit,*) 'Error: N3iwb must be greater or equal to N4iwb'
      write(ounit,*) 'N3iwb=',n3iwb,'  N4iwb=',iwbmax_small
    end if
    return
  end if

  if (n3iwf .lt. iwfmax_small .and. external_threelegs) then
    er = 2
    if (ounit .gt. 0) then
      write(ounit,*) 'Error: N3iwf must be greater or equal to N4iwf'
      write(ounit,*) 'N3iwf=',n3iwf,'  N4iwf=',iwfmax_small
    end if
    return
  end if

  if (n2iwb .lt. iwbmax_small .and. external_chi_loc) then
    er = 3
    if (ounit .gt. 0) then
      write(ounit,*) 'Error: N2iwb must be greater or equal to N4iwb'
      write(ounit,*) 'N2iwb=',n2iwb,'  N4iwb=',iwbmax_small
    end if
    return
  end if

  ! in order to calculate the susceptibility (lower leg: v-w)
  ! the DMFT frequency box has to be larger(>=) than
  ! v_vertex + w_vertex, so that v_vertex - ( -w_vertex)
  ! is still contained.
  if (iwmax .lt. (iwfmax_small + iwbmax_small)) then
    er = 4
    if (ounit .gt. 0) then
      write(ounit,*) 'Error: N1iwf must be greater than N4iwf + N4iwb'
      write(ounit,*) 'N1iwf=',iwmax,'  N4iwf=',iwfmax_small,'  N4iwb=',iwbmax_small
    endif
    return
  endif


  if (do_cond) then
    if (iwbcond < 0) then
      er = 5
      write(ounit,*) 'Error: Number of frequencies for conductivity must be >= 0'
      return
    endif

    if (iwbcond > iwfmax_small) then
      er = 6
      if (ounit .gt. 0) then
        write(ounit,*) 'Error: Number of frequencies for conductivity must be smaller than fermionic box'
        write(ounit,*) 'N1bc=',iwbcond,'  N4iwf=',iwfmax_small
      endif
      return
    else
      iwfcond = iwfmax_small-iwbcond
    endif

    if (extend_cond_bubble) then
      iwcstart = -iwmax+iwbcond
      iwcstop  = iwmax-iwbcond-1
    else
      iwcstart = -iwfcond
      iwcstop  = iwfcond-1
    endif
  endif


  if (ounit .gt. 0) write(ounit,'(1x)')

end subroutine check_freq_range


subroutine check_config(er,erstr)
  implicit none
  integer,intent(out) :: er
  character(len=*),intent(out) :: erstr
  logical :: there
  integer :: ineq

  er = 0
  erstr = ''
  exist_p = .false.
  do ineq=1,nineq
    if(ndims(ineq,2) .ne. 0) exist_p=.true.
  enddo

  if (read_ext_hk .eqv. .true.) then
    inquire (file=trim(filename_hk), exist=there)
    if (.not. there) then
      er = 1
      erstr = TRIM(ADJUSTL(erstr))//"Error: Can not find the Ham.hk file: "//trim(filename_hk)
      return
    endif
  endif

  if (read_ext_u .eqv. .true.) then
    inquire (file=trim(filename_umatrix), exist=there)
    if (.not. there) then
      er = 2
      erstr = TRIM(ADJUSTL(erstr))//"Error: Can not find the Umatrix file: "//trim(filename_umatrix)
      return
    endif
  endif

  if (q_vol .eqv. .false.) then
    inquire (file=trim(filename_qdata), exist=there)
    if (.not. there) then
      er = 3
      erstr = TRIM(ADJUSTL(erstr))//"Error: Can not find the Q-Path file: "//trim(filename_qdata)
      return
    endif
  endif

  if (do_vq .eqv. .true.) then
    inquire (file=trim(filename_vq), exist=there)
    if (.not. there) then
      er = 4
      erstr = TRIM(ADJUSTL(erstr))//"Error: Can not find the V(q) file: "//trim(filename_vq)
      return
    endif
  endif

  inquire (file=trim(filename_1p), exist=there)
  if (.not. there) then
      er = 5
      erstr = TRIM(ADJUSTL(erstr))//"Error: Can not find the 1P data file: "//trim(filename_1p)
      return
  endif

  inquire (file=trim(filename_vertex_sym), exist=there)
  if (.not. there) then
      er = 6
      erstr = TRIM(ADJUSTL(erstr))//"Error: Can not find the 2P data file: "//trim(filename_vertex_sym)
      return
  endif

  if (vertex_type .lt. 0 .or. vertex_type .gt. 2) then
    er = 7
    erstr = "Error: Choose appropriate vertex type"
  endif


  if (k_path_eom) then
    inquire (file=trim(filename_kdata), exist=there)
    if (.not. there) then
        er = 8
        erstr = TRIM(ADJUSTL(erstr))//"Error: Can not find the K-Path file: "//trim(filename_kdata)
        return
    endif
  endif

  if (do_vq .and. do_cond) then
    er = 9
    erstr = TRIM(ADJUSTL(erstr))//"Error: Conductivity calculation does not support V(q) at the moment"
    return
  endif

  if (.not. bse_inversion .and. do_cond) then
    er = 9
    erstr = TRIM(ADJUSTL(erstr))//"Error: Conductivity calculation does not support the use of geometric series"
    return
  endif

  return
end subroutine check_config

end module config_module
