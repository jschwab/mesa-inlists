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

      subroutine ZJ_capped_other_wind(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         real(dp) :: P0, vinf, Mdot_max

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         vinf = 300 * 1d5 ! cm/s
         P0 = 0.012 * pow_cr(Rsurf/Rsun, 1.86d0) * pow_cr(Msurf/Msun,-0.73d0)

         w = 4d-13*(Lsurf*Rsurf/Msurf)/(Lsun*Rsun/Msun)
         
         if (P0 .gt. 100) then
            w = 4.83d-9 * pow_cr(Msurf/Msun,-2.1d0) * pow_cr(Lsurf/Lsun,2.7d0) * w
         end if

         w = w * s% Reimers_scaling_factor
         
         if (s% x_logical_ctrl(2)) then
            Mdot_max = Lsurf/(vinf*vinf)/(Msun/secyer)
            w = min(w, Mdot_max)
         end if

         if (s% x_logical_ctrl(3)) then
            Mdot_max = 10**(-s% x_ctrl(3))
            w = min(w, Mdot_max)
         end if
         
         ierr = 0
      end subroutine ZJ_capped_other_wind

      subroutine frac_max_other_wind(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf ! surface values (cgs)
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         real(dp) :: P0, vinf, Mdot_max

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         vinf = 300 * 1d5 ! cm/s
         Mdot_max = 2d0*Lsurf/(vinf*vinf)/(Msun/secyer)
         w = s% x_ctrl(4) * Mdot_max
         
         ierr = 0
      end subroutine frac_max_other_wind
        
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
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
         s% job% warn_run_star_extras = .false.

         !s% other_wind => ZJ_capped_other_wind
         s% other_wind => frac_max_other_wind

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
            ! allow to turn off gold tolerances at the start
            if (s% x_logical_ctrl(1)) s% use_gold_tolerances = .false.
         else ! it is a restart
            call unpack_extra_info(s)
         end if
         
      end function extras_startup
      

      integer function extras_start_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


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
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns

      subroutine how_many_extra_history_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols=0
      end subroutine how_many_extra_history_header_items
      
      subroutine data_for_extra_history_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra history header item
      !set num_cols=1 in how_many_extra_history_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items


      subroutine how_many_extra_profile_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols = 0
      end subroutine how_many_extra_profile_header_items
      
      subroutine data_for_extra_profile_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra profile header item
      !set num_cols=1 in how_many_extra_profile_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         integer :: f
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         
         if (s% x_logical_ctrl(1)) then
            if ((s% L_surf .gt. 1e3) .and. (s% Teff .lt. 1e4)) then 
               s% use_gold_tolerances = .true.
            end if
         end if
         
         ! save profiles regularly in c core mass
         f = s% x_integer_ctrl(1)
         s% xtra1 = s% c_core_mass

         if ((floor(f * s% xtra1_old) - floor(f * s% xtra1) .ne. 0)) then

            ! save a profile & update the history
            s% need_to_update_history_now = .true.
            s% need_to_save_profiles_now = .true.

         endif

         ! time to stop if we've gone back to the blue
         ! we know we've "gone back" if there's not that much He left
         if (s% star_mass_he4 .lt. 0.05) then
            if (s% Teff .gt. 3e4) then
               extras_finish_step = terminate
               termination_code_str(t_xtra1) = 'back to blue'
               s% termination_code = t_xtra1
               return
            endif
         end if
         
         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         integer :: k, nzhi, iounit
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is a hack
         ! when the model stops, we save the cell location of the bottom of the envelope
         ! to an inlist and then we will read this when we start again

         do k = 1, s% nz
            if (s% T(k) .gt. 1e8) exit
         end do

         open(newunit=iounit, file='inlist_nzhi')
         write(iounit, *) '&star_job'
         write(iounit, '(A,I4)') '  set_nzhi = ', k
         write(iounit, *) '/'
         close(iounit)

         
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
      