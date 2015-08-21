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

      real(dp) :: logrhoc, logrhoc_old
      real(dp) :: logyec, logyec_old
      
      ! these routines are called by the standard run_star check_model
      contains

      real(dp) function center_avg_x(s,j)
         type (star_info), pointer :: s
         integer, intent(in) :: j
         real(dp) :: sum_x, sum_dq, dx, dq
         integer :: k
         sum_x = 0
         sum_dq = 0
         do k = s% nz, 1, -1
            dq = s% dq(k)
            dx = s% xa(j,k)*dq
            if (sum_dq+dq >= s% center_avg_value_dq) then
               sum_x = sum_x+ dx*(s% center_avg_value_dq - sum_dq)/dq
               sum_dq = s% center_avg_value_dq
               exit
            end if
            sum_x = sum_x + dx
            sum_dq = sum_dq + dq
         end do
         center_avg_x = sum_x/sum_dq
      end function center_avg_x

      subroutine extras_controls(s, ierr)

         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

      end subroutine extras_controls


      integer function extras_startup(s, id, restart, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         ierr = 0
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(s, id, id_extra)

         use chem_def, only: img24, ina24, ioo
         use chem_def, only : i_burn_ne, category_name
         use rates_def, only : i_rate

         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer :: mg24, na24, k
         integer :: i_burn_max
         real(dp) :: center_mg24, center_na24, x_limit, orunaway
         extras_check_model = keep_going

         ! ! get isotopes
         ! mg24 = s% net_iso(img24)
         ! na24 = s% net_iso(ina24)

         ! center_mg24 = center_avg_x(s,mg24)
         ! center_na24 = center_avg_x(s,na24)

         ! x_limit = s% x_ctrl(1)

         ! if ( (center_na24 < x_limit) .and. (center_mg24 < x_limit)) then

         !   ! stop when isotopes drop too lo
         !   extras_check_model = terminate
         !   s% termination_code = t_xtra1
         !   termination_code_str(t_xtra1) = 'Mg24 & Na24 depletion'
         !   return
         ! end if

         ! if (s% log_center_density .gt. 9.7) s% mass_change = 0

         k = s% max_eps_nuc_k
         orunaway = s% eps_nuc_categories(i_rate,ioo,k) - s% non_nuc_neu(k)

         if (orunaway .gt. 0) then
            extras_check_model = terminate
            s% termination_code = t_xtra2
            termination_code_str(t_xtra2) = 'O+O runaway'
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depenending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination conditon'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         how_many_extra_history_columns = 7
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(s, id, id_extra, n, names, vals, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr

         real(dp) :: center_CpT

         !note: do NOT add these names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         center_CpT = s% Cp(s% nz) * s% T(s% nz)
         
         names(1) = 'tc_accretion'
         vals(1) = safe_log10_cr(s% star_mass / s% star_mdot * secyer)

         logrhoc = log_cr(s% Rho(s% nz))
         names(2) = 'tc_compression'
         vals(2) = safe_log10_cr(s% dt_old / abs(logrhoc - logrhoc_old))
         logrhoc_old = logrhoc
         
         names(3) = 'tc_dynamic'
         vals(3) = 1d0 / sqrt(standard_cgrav * s% Rho(s% nz))

         logyec = log_cr(s% center_ye)
         names(4) = 'tc_weak'
         vals(4) = safe_log10_cr(s% dt_old / abs(logyec - logyec_old))
         logyec_old = logyec
         
         names(5) = 'tc_neutrino'
         vals(5) = safe_log10_cr(center_CpT / s% center_non_nuc_neu)

         names(6) = 'tc_grav'
         vals(6) = safe_log10_cr(center_CpT / s% center_eps_grav)

         names(7) = 'tc_nuc'
         vals(7) = safe_log10_cr(center_CpT / s% center_eps_nuc)

         
         ierr = 0
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 2
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(s, id, id_extra, n, nz, names, vals, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         real(dp) :: radiative_conductivity
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0

         !note: do NOT add these names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'

         names(1) = 'semiconvectionD'
         do k = 1, nz
            radiative_conductivity = (4*crad*clight / 3)*s% T(k)*s% T(k)*s% T(k) / (s% opacity(k)*s% rho(k)) ! erg / (K cm sec)
            vals(k,1) = radiative_conductivity/(6*s% Cp(k) *s% rho(k)) &
                  *(s% gradT(k) - s% grada(k))/(s% gradL(k) - s% gradT(k))
         end do

         names(2) = 'B'
         do k = 1, nz
            vals(k,2) = s% gradL(k) - s% grada(k)
         end do

      end subroutine data_for_extra_profile_columns

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(s, id, id_extra)

         use chem_def, only: img24, ina24, ioo
         use chem_def, only : i_burn_ne, category_name
         use rates_def, only : i_rate

         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer :: ierr
         integer :: f

         integer :: mg24, na24, k
         integer :: i_burn_max
         real(dp) :: center_mg24, center_na24, x_limit, orunaway

         extras_finish_step = keep_going
         call store_extra_info(s)

         ! ! get isotopes
         ! mg24 = s% net_iso(img24)
         ! center_mg24 = center_avg_x(s,mg24)
         s% xtra1 = s% log_center_density

         f = 100
         ! this expression will evaluate to true if f times the log center density
         ! has crossed an integer during the last step.  If f = 5, then we will get
         ! output at log center density = {... 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 ... }
         if ((floor(f * s% xtra1_old) - floor(f * s% xtra1) .ne. 0)) then

            ! save a profile & update the history
            s% need_to_update_history_now = .true.
            s% need_to_save_profiles_now = .true.

            ! by default the priority is 1; you can change that if you'd like
            s% save_profiles_model_priority = 3

            write(*,*) "saving profile at rho_c", s% xtra1

         endif

         ! ! see extras_check_model for information about custom termination codes
         ! ! by default, indicate where (in the code) MESA terminated
         ! if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

         ! ! now, if we're getting close to the critical mass, wind the timestep way down
         ! if (s% star_mass > 1.375) then
         !    s% delta_lgRho_cntr_limit = 0.0001
         ! endif

      end function extras_finish_step


      subroutine extras_after_evolve(s, id, id_extra, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         ierr = 0
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
