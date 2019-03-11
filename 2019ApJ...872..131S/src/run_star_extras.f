! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_star_extras

  use star_lib
  use star_def
  use const_def
  use crlibm_lib

  implicit none

  integer, parameter :: extra_info_alloc = 1
  integer, parameter :: extra_info_get = 2
  integer, parameter :: extra_info_put = 3

  ! these routines are called by the standard run_star check_model
contains

  subroutine how_many_other_mesh_fcns(id, n)
    integer, intent(in) :: id
    integer, intent(out) :: n
    type (star_info), pointer :: s
    integer :: ierr

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    n = 0
    if (s% x_logical_ctrl(1)) n = n+1
    if (s% x_logical_ctrl(2)) n = n+1
  end subroutine how_many_other_mesh_fcns

  function urca(s,k) result(res)
    type (star_info), pointer :: s
    integer, intent(in) :: k
    real(dp) :: res
    res = s% eta(k) / s% x_ctrl(2)
  end function urca

  subroutine urca_other_mesh_fcn_data( &
       id, nfcns, names, gval_is_xa_function, vals1, ierr)
    integer, intent(in) :: id
    integer, intent(in) :: nfcns
    character (len=*) :: names(:)
    logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
    real(dp), pointer :: vals1(:) ! =(nz, nfcns)
    integer, intent(out) :: ierr
    integer :: nz, k
    real(dp), pointer :: vals(:,:)
    real(dp) :: weight
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    nz = s% nz
    vals(1:nz,1:nfcns) => vals1(1:nz*nfcns)

    ! ensure N zones per decade in radius
    if (s% x_logical_ctrl(1)) then
       names(1) = 'logR'
       gval_is_xa_function(1) = .false.
       do k=1,nz
          vals(k,1) = s% x_ctrl(1) * log10_cr(s% r(k))
       end do
    end if

    ! remesh based on chemical potential
    if (s% x_logical_ctrl(2)) then
       names(2) = 'urca1'
       gval_is_xa_function(1) = .false.
       do k=1,nz
          vals(k,2) = urca(s, k)
       end do
    end if

  end subroutine urca_other_mesh_fcn_data

  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! this is the place to set any procedure pointers you want to change
    ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
    s% how_many_other_mesh_fcns => how_many_other_mesh_fcns
    s% other_mesh_fcn_data => urca_other_mesh_fcn_data

    ! Uncomment these lines if you wish to use the functions in this file,
    ! otherwise we use a null_ version which does nothing.
    s% extras_startup => extras_startup
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    !s% how_many_extra_history_columns => how_many_extra_history_columns
    !s% data_for_extra_history_columns => data_for_extra_history_columns
    !s% how_many_extra_profile_columns => how_many_extra_profile_columns
    !s% data_for_extra_profile_columns => data_for_extra_profile_columns

    ! Once you have set the function pointers you want,
    ! then uncomment this (or set it in your star_job inlist)
    ! to disable the printed warning message,
    s% job% warn_run_star_extras =.false.

  end subroutine extras_controls


  integer function extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_startup = 0
    if (.not. restart) then
       call alloc_extra_info(s)
    else ! it is a restart
       call unpack_extra_info(s)
    end if
  end function extras_startup


  ! returns either keep_going, retry, backup, or terminate.
  integer function extras_check_model(id, id_extra)

    use chem_def

    integer, intent(in) :: id, id_extra
    integer :: i, k, i_burn_max, nz, ierr
    type (star_info), pointer :: s
    real(dp) :: max_diff_eta, orunaway, crunaway
    real(dp) :: dlogmg24, dlogne24, dlogmax
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going

    ! apply varcontrol-esque timestep limit to eta
    if (s% x_logical_ctrl(3)) then

       i = maxloc(abs(s% eta(1:s% nz) - s% eta_start(1:s% nz)), dim=1)
       max_diff_eta = s% eta(i) - s% eta_start(i)

       if (abs(max_diff_eta) .gt. s% x_ctrl(3)) then
          extras_check_model = redo
          s% dt = 0.5d0 * s% dt
          write(*,*) 'redo for delta eta limit'
       endif

    endif

    ! get fine time resolution around central captures
    nz = s% nz
    dlogmg24 = abs(log10_cr(s% xa(s% net_iso(img24), nz)) - &
         log10_cr(s% xa_old(s% net_iso(img24), nz)))
    dlogne24 = abs(log10_cr(s% xa(s% net_iso(ine24), nz)) - &
         log10_cr(s% xa_old(s% net_iso(ine24), nz)))

    dlogmax = 0
    if (min(s% xa(s% net_iso(img24), nz), s% xa(s% net_iso(ine24), nz)) > 1e-6) then
       if (s% xa(s% net_iso(img24), nz) > s% xa(s% net_iso(ine24), nz)) then
          dlogmax = dlogne24
       else
          dlogmax = dlogmg24
       end if
    end if


    ! ! terminate at carbon runaway
    ! k = s% max_eps_nuc_k
    ! crunaway = s% eps_nuc_categories(icc,k) - s% non_nuc_neu(k)
    
    
    ! if ((crunaway .gt. 0) .and. (s%T(k) .gt. 1e8)) then
    !    extras_check_model = terminate
    !    s% termination_code = t_xtra2
    !    termination_code_str(t_xtra2) = 'C+C runaway'
    !    return
    ! end if

    
    ! terminate at oxygen runaway
    k = s% max_eps_nuc_k
    orunaway = s% eps_nuc_categories(ioo,k) - s% non_nuc_neu(k)

    if ((orunaway .gt. 0) .and. (s%T(k) .gt. 8e8)) then
       extras_check_model = terminate
       s% termination_code = t_xtra2
       termination_code_str(t_xtra2) = 'O+O runaway'
       return
    end if


    ! if you want to check multiple conditions, it can be useful
    ! to set a different termination code depending on which
    ! condition was triggered.  MESA provides 9 customizeable
    ! termination codes, named t_xtra1 .. t_xtra9.  You can
    ! customize the messages that will be printed upon exit by
    ! setting the corresponding termination_code_str value.
    ! termination_code_str(t_xtra1) = 'my termination condition'

    ! by default, indicate where (in the code) MESA terminated
    if (extras_check_model == terminate) s% termination_code = t_extras_check_model
  end function extras_check_model


  integer function how_many_extra_history_columns(id, id_extra)
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_columns = 13
  end function how_many_extra_history_columns


  subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)

    use chem_def

    integer, intent(in) :: id, id_extra, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k, nz, col, kbot

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    !note: do NOT add the extras names to history_columns.list
    ! the history_columns.list is only for the built-in log column options.
    ! it must not include the new column names you are adding here.

    nz = s% nz

    names = 'unset'
    vals = 0

    names(1) = 'num_unstable_zones'
    vals(1) = 0
    do k = 1, nz
       if (s% brunt_N2(k) .lt. 0) vals(1) = vals(1) + 1
    enddo

    names(2) = 'center_t_heat'
    vals(2) =  s% cp(nz) * s% T(nz) / s% eps_nuc(nz)

    names(3) = 'center_t_compress'
    vals(3) = 1d0 / s% dlnd_dt(nz)

    names(4) = 'center_fermi_energy'
    vals(4) = s% eta(nz) * (1d-6 * kev * s% T(nz))

    names(5) = 'brunt_B_center'
    vals(5) = s% brunt_B(nz)

    names(6) = 'brunt_nonB_center'
    vals(6) = -s% gradT_sub_grada(nz)

    names(7) = 'unstable_mass'
    vals(7) = 0
    do k = 1, nz
       if (s% brunt_N2(k) .lt. 0) vals(7) = vals(7) + s% dm(k)
    enddo

    names(8) = 'unstable_core_top_r'
    names(9) = 'unstable_core_top_m'
    vals(8) = -1
    vals(9) = -1
    ! first, check if the center is unstable
    if (s% brunt_N2(nz) .lt. 0) then
       write(*,*) 'doing check'
       ! start at center
       ! go out until we find where it stabilizes
       do k = nz, 1, -1
          if (s% brunt_N2(k) .gt. 0) then
             vals(8) = s% r(k)
             vals(9) = s% m(k)/Msun
             exit
          end if
       enddo
    end if

    names(10) = 'unstable_shell_bot_r'
    names(11) = 'unstable_shell_bot_m'
    vals(10) = -1
    vals(11) = -1
    do k = nz-1, 1, -1
       if ((s% brunt_N2(k+1) .gt. 0) .and. (s% brunt_N2(k) .lt. 0)) then
          vals(10) = s% r(k)
          vals(11) = s% m(k)/Msun
          kbot = k
          exit
       end if
    enddo


    names(12) = 'unstable_shell_top_r'
    names(13) = 'unstable_shell_top_m'
    vals(12) = -1
    vals(13) = -1
    do k = kbot-1, 1, -1
       if ((s% brunt_N2(k+1) .lt. 0) .and. (s% brunt_N2(k) .gt. 0)) then
          vals(12) = s% r(k)
          vals(13) = s% m(k)/Msun
          exit
       end if
    enddo

    col = 13
    ! get radii associated with various weak "shells"
    call trace_shell("mgnane24",col, img24, ine24)
    call trace_shell("nane23",  col, ina23, ine23)
    call trace_shell("mgna25",  col, img25, ina25)
    call trace_shell("nane25",  col, ina25, ine25)

    ! names(50) = foo

  contains

    subroutine trace_shell(name, col, iso1, iso2)
      character(len=*), intent(in) :: name
      integer, intent(inout) :: col
      integer, intent(in) :: iso1, iso2

      integer :: jx, jy, idx
      jx = s% net_iso(iso1)
      jy = s% net_iso(iso2)

      do k = 1, nz
         if (s% xa(jy,k) > s% xa(jx, k)) then
            idx = k
            exit
         endif
      enddo

      names(col+1) = trim(name) // '_r'
      names(col+2) = trim(name) // '_m'
      names(col+3) = trim(name) // '_Rho'
      names(col+4) = trim(name) // '_T'
      names(col+5) = trim(name) // '_eta'
      names(col+6) = trim(name) // '_H'

      if ((idx .ge. nz) .or. (idx .le. 1)) then
         vals(col+1:col+6) = -1
      else
         vals(col+1) = s% r(idx)
         vals(col+2) = s% m(idx) / Msun
         vals(col+3) = s% Rho(idx)
         vals(col+4) = s% T(idx)
         vals(col+5) = s% eta(idx)
         vals(col+6) = s% scale_height(idx)
      end if

      col = col + 6

    end subroutine trace_shell



  end subroutine data_for_extra_history_columns


  integer function how_many_extra_profile_columns(id, id_extra)
    use star_def, only: star_info
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 6
  end function how_many_extra_profile_columns


  subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
    use star_def, only: star_info, maxlen_profile_column_name
    use const_def, only: dp
    use chem_def, only: ina24, ine24
    integer, intent(in) :: id, id_extra, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: j, k
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    !note: do NOT add the extra names to profile_columns.list
    ! the profile_columns.list is only for the built-in profile column options.
    ! it must not include the new column names you are adding here.

    ! here is an example for adding a profile column
    !if (n /= 1) stop 'data_for_extra_profile_columns'

    names(1) = 'mesh_urca'
    if (s% x_logical_ctrl(2)) then 
       do k = 1, nz
          vals(k,1) = urca(s, k)
       end do
    else
       vals(:,1) = -1
    end if
    

    names(2) = 'sign_brunt_N2'
    j = s% net_iso(ine24)
    do k = 1, nz
       vals(k,2) = sign(1d0, s% brunt_N2(k))
    end do

    names(3) = 'logq'
    do k = 1, nz
       vals(k,3) = log10_cr(s% q(k))
    end do

    names(4) = 'nz_minus_k'
    do k = 1, nz
       vals(k,4) = (s% nz + 1) - k
    end do

    names(5) = 'log_t_heat_nuc'
    do k = 1, nz
       vals(k,5) =  safe_log10_cr(s% cp(k) * s% T(k) / s% eps_nuc(k))
    end do

    names(6) = 'log_t_heat_cond'
    do k = 1, nz-1
       vals(k,6) = safe_log10_cr(s% dm(k) * s% cp(k) * s% T(k) / abs(s% L(k)-s%L(k+1)))
    end do
    vals(nz, 6) = 1d99

    ! names(2) = 'energy_per_na24_capture'
    ! j = s% net_iso(ine24)
    ! do k = 1, nz
    !   vals(k,2) = 24.0 * mp * s% eps_nuc(k) / abs(s% dxdt_nuc(j, k)) / mev_to_ergs
    ! end do

  end subroutine data_for_extra_profile_columns


  ! returns either keep_going or terminate.
  ! note: cannot request retry or backup; extras_check_model can do that.
  integer function extras_finish_step(id, id_extra)
    type (star_info), pointer :: s
    integer, intent(in) :: id, id_extra
    integer :: ierr
    integer :: f
    extras_finish_step = keep_going

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    call store_extra_info(s)

    ! MESA provides a number of variables that make it easy to get user input.
    ! these are part of the star_info structure and are named
    ! x_character_ctrl, x_integer_ctrl, x_logical_ctrl, and x_ctrl.
    ! by default there are num_x_ctrls, which defaults to 100, of each.
    ! they can be specified in the controls section of your inlist.

    f = s% x_integer_ctrl(1)

    ! MESA also provides a number variables that are useful for implementing
    ! algorithms which require a state. if you just use these variables
    ! restarts, retries, and backups will work without doing anything special.
    ! they are named xtra1 .. xtra30, ixtra1 .. ixtra30, and lxtra1 .. lxtra30.
    ! they are automatically versioned, that is if you set s% xtra1, then
    ! s% xtra1_old will contains the value of s% xtra1 from the previous step
    ! and s% xtra1_older contains the one from two steps ago.

    if (s% log_center_density .lt. 9.0) return
    s% xtra1 = s% log_center_density

    ! this expression will evaluate to true if f times the log center density
    ! has crossed an integer during the last step.  If f = 5, then we will get
    ! output at log center density = {... 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 ... }
    if ((floor(f * s% xtra1_old) - floor(f * s% xtra1) .ne. 0)) then

       ! save a profile & update the history
       s% need_to_update_history_now = .true.
       s% need_to_save_profiles_now = .true.

       ! by default the priority is 1; you can change that if you'd like
       s% save_profiles_model_priority = 3

    endif

    ! see extras_check_model for information about custom termination codes
    ! by default, indicate where (in the code) MESA terminated
    if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

  end function extras_finish_step


  subroutine extras_after_evolve(id, id_extra, ierr)
    integer, intent(in) :: id, id_extra
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine extras_after_evolve


  ! routines for saving and restoring extra data so can do restarts

  subroutine alloc_extra_info(s)
    type (star_info), pointer :: s
    call move_extra_info(s,extra_info_alloc)
  end subroutine alloc_extra_info


  subroutine unpack_extra_info(s)
    type (star_info), pointer :: s
    call move_extra_info(s,extra_info_get)
  end subroutine unpack_extra_info


  subroutine store_extra_info(s)
    type (star_info), pointer :: s
    call move_extra_info(s,extra_info_put)
  end subroutine store_extra_info


  subroutine move_extra_info(s,op)
    type (star_info), pointer :: s
    integer, intent(in) :: op

    integer :: i, j, num_ints, num_dbls, ierr

    i = 0
    ! call move_int or move_flg
    num_ints = i

    i = 0
    ! call move_dbl

    num_dbls = i

    if (op /= extra_info_alloc) return
    if (num_ints == 0 .and. num_dbls == 0) return

    ierr = 0
    call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in star_alloc_extras'
       write(*,*) 'alloc_extras num_ints', num_ints
       write(*,*) 'alloc_extras num_dbls', num_dbls
       stop 1
    end if

  contains

    subroutine move_dbl(dbl)
      real(dp) :: dbl
      i = i+1
      select case (op)
      case (extra_info_get)
         dbl = s% extra_work(i)
      case (extra_info_put)
         s% extra_work(i) = dbl
      end select
    end subroutine move_dbl

    subroutine move_int(int)
      integer :: int
      i = i+1
      select case (op)
      case (extra_info_get)
         int = s% extra_iwork(i)
      case (extra_info_put)
         s% extra_iwork(i) = int
      end select
    end subroutine move_int

    subroutine move_flg(flg)
      logical :: flg
      i = i+1
      select case (op)
      case (extra_info_get)
         flg = (s% extra_iwork(i) /= 0)
      case (extra_info_put)
         if (flg) then
            s% extra_iwork(i) = 1
         else
            s% extra_iwork(i) = 0
         end if
      end select
    end subroutine move_flg

  end subroutine move_extra_info


end module run_star_extras
