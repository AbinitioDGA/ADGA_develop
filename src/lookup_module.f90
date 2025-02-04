module lookup_module
  use parameters_module
  implicit none
  private
  integer :: i,j,pst
  character(len=150) :: str_temp, str_split
  public :: string_find, int_find, int3_find, float_find, bool_find, group_find, subgroup_find, &
            spell_check

  contains

  subroutine string_find(search_string, save_string, search_start, search_end)
    character(*), intent(in)  :: search_string
    character(len=150), intent(inout) :: save_string ! keep default string
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        save_string=trim(adjustl(str_temp(pst+1:)))
      endif
    enddo
  end subroutine string_find

  subroutine int_find(search_string, save_int, search_start, search_end)
    character(*), intent(in)  :: search_string
    integer, intent(inout) :: save_int ! keep default values
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=trim(adjustl(str_temp(pst+1:)))
        read(str_temp,*) save_int
      endif
    enddo
  end subroutine int_find

  subroutine int3_find(search_string, save_int1, save_int2, save_int3, search_start, search_end)
    character(*), intent(in)  :: search_string
    integer, intent(inout) :: save_int1, save_int2, save_int3 ! keep default values
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=trim(adjustl(str_temp(pst+1:)))
        pst=scan(str_temp,multseperator)
        str_split=trim(adjustl(str_temp(:pst-1)))
        read(str_split,*) save_int1
        str_temp=trim(adjustl(str_temp(pst+1:)))
        pst=scan(str_temp,multseperator)
        str_split=trim(adjustl(str_temp(:pst-1)))
        str_temp=trim(adjustl(str_temp(pst+1:)))
        read(str_split,*) save_int2
        read(str_temp,*) save_int3
      endif
    enddo
  end subroutine int3_find

  subroutine float_find(search_string, save_float, search_start, search_end)
    character(*), intent(in)  :: search_string
    real(8), intent(inout) :: save_float ! keep default values
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=trim(adjustl(str_temp(pst+1:)))
        read(str_temp,*) save_float
      endif
    enddo
  end subroutine float_find

  subroutine bool_find(search_string, save_bool, search_start, search_end)
    character(*), intent(in)  :: search_string
    logical, intent(inout) :: save_bool
    integer, intent(in) :: search_start, search_end

    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,seperator)
        str_temp=trim(adjustl(str_temp(pst+1:)))
        read(str_temp,*) save_bool
      endif
    enddo
  end subroutine bool_find

  subroutine group_find(search_string, save_start, save_end)
    character(*), intent(in) :: search_string
    integer, intent(out) :: save_start, save_end
    save_start=0
    save_end=0

    do i=1,lines
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        save_start=i+1
        exit
      endif
    enddo

    if (save_start .ge. 1) then ! group was found
      do i=save_start, lines
        if (index(trim(file_save(i)),'[') .eq. 1) then
          if (index(trim(file_save(i)),'[[') .eq. 1) then ! skip subgroups
            cycle
          endif
          save_end=i-1 ! one above the next session
          exit
        endif
      enddo

      if(save_end .eq. 0) then
        save_end = lines ! if nothing else is found, until the end of the file
      endif

      if (save_start .gt. save_end) then ! group found, but no content
        save_start = -1
      endif
    endif
    return
    ! save_start -> 0: not found; -1: found, but empty
  end subroutine group_find

  subroutine subgroup_find(search_string, search_start, search_end, save_start, save_end)
    character(*), intent(in) :: search_string
    integer, intent(in) :: search_start, search_end
    integer, intent(out) :: save_start, save_end
    save_start=0
    save_end=0

    do i=search_start, search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        save_start=i+1
        exit
      endif
    enddo

    if (save_start .ge. 1) then ! subgroup found
      do i=save_start, search_end
        if (index(trim(file_save(i)),'[') .eq. 1) then
          save_end=i-1 ! one above the next session
          exit
        endif
      enddo

      if(save_end .eq. 0) then
        save_end = search_end ! if nothing else is found, until the end of the group
                              ! whose size was already determined by group_find
      endif

      if (save_start .gt. save_end) then ! subgroup found, but no content
        save_start = -1
      endif
    endif
    return
    ! save_start -> 0: not found; -1: found, but empty
  end subroutine subgroup_find

  subroutine spell_check(search_start, search_end, grname, dictionary, er, erstr)
    character(*), intent(in) :: grname
    character(*), intent(in) :: dictionary(:)
    integer, intent(in) :: search_start, search_end
    integer, intent(out) :: er
    character(len=200), intent(out) :: erstr

    do i=search_start,search_end
      str_temp = file_save(i)
      pst=scan(str_temp,seperator)
      if (pst .eq. 0) then
        er = 100
        erstr = 'Variable in '//trim(grname)//' group without argument: '//str_temp
        return
      endif
      str_temp=trim(adjustl(str_temp(:(pst-1))))
      er = 101
      do j = 1,size(dictionary)
        if (str_temp == dictionary(j)) then
          er = 0
          exit
        endif
      enddo
      if (er .ne. 0) then
        erstr = 'Spelling error or unknown variable in '//trim(grname)//' group: '//str_temp
        return
      endif
    enddo
  end subroutine spell_check

end module lookup_module
