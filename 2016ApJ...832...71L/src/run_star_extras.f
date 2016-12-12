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

  ! these routines are called by the standard run_star check_model
contains

  subroutine my_other_D_mix(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr

    integer :: k, k_peak, nz, k_conv_lower, k_conv_upper
    real(dp) :: z, other_D_mix, sig, N2_peak

    logical :: mix_here
    
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    nz = s% nz
    k_peak = maxloc(s% brunt_N2(1:nz),1)
    N2_peak = s% brunt_N2(k_peak)

    write(*,*) "called other D mix"

    ! unified mixing routine
    ! x_ctrl(1)
    !     mix for N2 < x_ctrl(1)**2 * N2peak

    ! x_logical_ctrl(1)
    !    .true.  -> D_mix = D_th / x_ctrl(2)
    !    .false. ->  D_mix = x_ctrl(2)

    do k = nz, 1, -1
       z = (s% R(k) - s% R(k_peak)) / s% scale_height(k_peak)
       if ((z .gt. 0.0) .and. (z .lt. 0.3)) then ! don't let things go too far
          if (s% brunt_N2(k) .lt. s% x_ctrl(1)**2 * N2_peak) then

             if (s% x_logical_ctrl(1)) then
                ! Lewis number mixing
                sig = 4d0*crad*clight*s% T(k)**3 / (s% Rho(k)*s% opacity(k))
                other_D_mix = sig / (s% cp(k) * s% Rho(k) * s% x_ctrl(2))
             else
                ! constant diffusivity mixing
                other_D_mix = s% x_ctrl(2)
             endif

             ! apply mixing
             if (other_D_mix .gt. s% D_mix(k)) then 
                s% D_mix(k) = other_D_mix
                s% mixing_type(k) = anonymous_mixing
             endif

          endif
       endif
    enddo

  end subroutine my_other_D_mix

  subroutine how_many_my_other_mesh_fcns(id, n)
    integer, intent(in) :: id
    integer, intent(out) :: n
    n = 1
  end subroutine how_many_my_other_mesh_fcns

  subroutine my_other_mesh_fcn_data( &
       id, nfcns, names, gval_is_xa_function, vals1, ierr)
    integer, intent(in) :: id
    integer, intent(in) :: nfcns
    character (len=*) :: names(:)
    logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
    real(dp), pointer :: vals1(:) ! =(nz, nfcns)
    integer, intent(out) :: ierr
    integer :: nz, k, k_peak
    real(dp), pointer :: vals(:,:)
    real(dp), parameter :: weight = 300
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    names(1) = 'kap_function'
    gval_is_xa_function(1) = .false.

    nz = s% nz

    k_peak = maxloc(s% brunt_N2(1:nz),1)
    
    vals(1:nz,1:nfcns) => vals1(1:nz*nfcns)
    do k=1,nz
       vals(k,1) = weight*tanh((s% R(k) - s%R(k_peak)) / (0.1 * s% scale_height(k_peak)))
    end do
  end subroutine my_other_mesh_fcn_data
  

  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! this is the place to set any procedure pointers you want to change
    ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

    ! Uncomment these lines if you wish to use the functions in this file,
    ! otherwise we use a null_ version which does nothing.
    s% extras_startup => extras_startup
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns  

    s% other_D_mix => my_other_D_mix
    
    s% how_many_other_mesh_fcns => how_many_my_other_mesh_fcns
    s% other_mesh_fcn_data => my_other_mesh_fcn_data

    ! Once you have set the function pointers you want,
    ! then uncomment this (or set it in your star_job inlist)
    ! to disable the printed warning message,
    s% job% warn_run_star_extras =.false.       

  end subroutine extras_controls

  ! None of the following functions are called unless you set their
  ! function point in extras_control.


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
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going         
    if (.false. .and. s% star_mass_h1 < 0.35d0) then
       ! stop when star hydrogen mass drops to specified level
       extras_check_model = terminate
       write(*, *) 'have reached desired hydrogen mass'
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
    how_many_extra_history_columns = 34
  end function how_many_extra_history_columns


  subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
    integer, intent(in) :: id, id_extra, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr

    integer :: k, k_heat, k_peak, k_burn, k_conv_lower, k_conv_upper

    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return


    !note: do NOT add these names to history_columns.list
    ! the history_columns.list is only for the built-in log column options.
    ! it must not include the new column names you are adding here.

    ierr = 0

    ! calculate all the things
    do k = s% nz, 1, -1
       if (s% dlnT_dt(k) .gt. 0) then 
          k_heat = k          
          exit
       endif
    end do

    k_peak = maxloc(s% brunt_N2(1:s% nz),1)


    do k = s% nz, 1, -1
       if (s% eps_nuc(k) - s% non_nuc_neu(k) .gt. 0) then
          k_burn = k
          exit
       end if
    end do


    do k = s% nz, 1, -1
       if (s% mixing_type(k) .gt. 0) then
          k_conv_lower = k
          exit
       end if
    end do

    do k = k_conv_lower-1, 1, -1
       if (s% mixing_type(k) .le. 0) then
          k_conv_upper = k
          exit
       end if
    end do

    k_heat = max(1, min(k_heat, s% nz))
    k_peak = max(1, min(k_peak, s% nz))
    k_burn = max(1, min(k_burn, s% nz))
    k_conv_lower = max(1, min(k_conv_lower, s% nz))
    k_conv_upper = max(1, min(k_conv_upper, s% nz))

    ! calculate w_c
    s% xtra1_array(1:s% nz) = 2 * pi * s% conv_vel(1:s% nz) &
         / (s% scale_height(1:s% nz))

    names(1) = 'w_conv_max'
    vals(1) = maxval(s% xtra1_array(k_conv_upper:k_conv_lower))

    ! calculate mach
    s% xtra2_array(1:s% nz) = s% conv_vel(1:s% nz) &
         / (s% csound(1:s% nz))

    names(2) = 'mach_conv_max'
    vals(2) = maxval(s% xtra2_array(k_conv_upper:k_conv_lower))

    ! upstream N2
    names(3) = 'brunt_N2_heat'
    vals(3) = s% brunt_N2(k_heat)

    ! peak N2
    names(4) = 'brunt_N2_peak'
    vals(4) = s% brunt_N2(k_peak)

    names(5) = 'k_heat'
    vals(5) = k_heat

    names(6) = 'k_peak'
    vals(6) = k_peak

    names(7) = 'k_burn'
    vals(7) = k_burn

    names(8) = 'k_conv_lower'
    vals(8) = k_conv_lower

    names(9) = 'k_conv_upper'
    vals(9) = k_conv_upper

    names(10) = 'q_heat'
    vals(10) = s% q(k_heat)

    names(11) = 'q_peak'
    vals(11) = s% q(k_peak)

    names(12) = 'q_burn'
    vals(12) = s% q(k_burn)

    names(13) = 'q_conv_lower'
    vals(13) = s% q(k_conv_lower)

    names(14) = 'q_conv_upper'
    vals(14) = s% q(k_conv_upper)

    names(15) = 'R_heat'
    vals(15) = s% R(k_heat)

    names(16) = 'R_peak'
    vals(16) = s% R(k_peak)

    names(17) = 'R_burn'
    vals(17) = s% R(k_burn)

    names(18) = 'R_conv_lower'
    vals(18) = s% R(k_conv_lower)

    names(19) = 'R_conv_upper'
    vals(19) = s% R(k_conv_upper)

    names(20) = 'scale_height_heat'
    vals(20) = s% scale_height(k_heat)

    names(21) = 'scale_height_peak'
    vals(21) = s% scale_height(k_peak)

    names(22) = 'scale_height_burn'
    vals(22) = s% scale_height(k_burn)

    names(23) = 'scale_height_conv_lower'
    vals(23) = s% scale_height(k_conv_lower)

    names(24) = 'scale_height_conv_upper'
    vals(24) = s% scale_height(k_conv_upper)

    names(25) = 'T_heat'
    vals(25) = s% T(k_heat)

    names(26) = 'T_peak'
    vals(26) = s% T(k_peak)

    names(27) = 'T_burn'
    vals(27) = s% T(k_burn)

    names(28) = 'T_conv_lower'
    vals(28) = s% T(k_conv_lower)

    names(29) = 'T_conv_upper'
    vals(29) = s% T(k_conv_upper)

    names(30) = 'Rho_heat'
    vals(30) = s% Rho(k_heat)

    names(31) = 'Rho_peak'
    vals(31) = s% Rho(k_peak)

    names(32) = 'Rho_burn'
    vals(32) = s% Rho(k_burn)

    names(33) = 'Rho_conv_lower'
    vals(33) = s% Rho(k_conv_lower)

    names(34) = 'Rho_conv_upper'
    vals(34) = s% Rho(k_conv_upper)


  end subroutine data_for_extra_history_columns


  integer function how_many_extra_profile_columns(id, id_extra)
    use star_def, only: star_info
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 4
  end function how_many_extra_profile_columns


  subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
    use star_def, only: star_info, maxlen_profile_column_name
    use const_def, only: dp
    integer, intent(in) :: id, id_extra, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k, k_peak

    real(dp) :: sig, Le

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    names(1) = 'Le'
    do k = 1, nz
       sig = 4d0*crad*clight*s% T(k)**3 / (s% Rho(k)*s% opacity(k))
       vals(k,1) = sig / (s% cp(k) * s% Rho(k) * s% D_mix(k))
    end do
    names(2) = 'log_Le'
    do k = 1, nz
       vals(k,2) = log10(vals(k,1))
    end do

    k_peak = maxloc(s% brunt_N2(1:s% nz),1)
    names(3) = 'z_div_H'
    do k = 1, nz
       vals(k,3) = (s% R(k) - s% R(k_peak)) / s% scale_height(k_peak)
    end do

    names(4) = 'log_N_div_Npeak'
    do k = 1, nz
       vals(k,4) = 0.5 * safe_log10_cr(s% brunt_N2(k) / s% brunt_N2(k_peak))
    end do
    
  end subroutine data_for_extra_profile_columns


  ! returns either keep_going or terminate.
  ! note: cannot request retry or backup; extras_check_model can do that.
  integer function extras_finish_step(id, id_extra)
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going
    call store_extra_info(s)

    ! to save a profile, 
    ! s% need_to_save_profiles_now = .true.
    ! to update the star log,
    ! s% need_to_update_history_now = .true.

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

  ! put these defs at the top and delete from the following routines
  !integer, parameter :: extra_info_alloc = 1
  !integer, parameter :: extra_info_get = 2
  !integer, parameter :: extra_info_put = 3


  subroutine alloc_extra_info(s)
    integer, parameter :: extra_info_alloc = 1
    type (star_info), pointer :: s
    call move_extra_info(s,extra_info_alloc)
  end subroutine alloc_extra_info


  subroutine unpack_extra_info(s)
    integer, parameter :: extra_info_get = 2
    type (star_info), pointer :: s
    call move_extra_info(s,extra_info_get)
  end subroutine unpack_extra_info


  subroutine store_extra_info(s)
    integer, parameter :: extra_info_put = 3
    type (star_info), pointer :: s
    call move_extra_info(s,extra_info_put)
  end subroutine store_extra_info


  subroutine move_extra_info(s,op)
    integer, parameter :: extra_info_alloc = 1
    integer, parameter :: extra_info_get = 2
    integer, parameter :: extra_info_put = 3
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

