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
      use chem_def, only : ic12

      implicit none

      real(dp) :: initial_c12_mass, initial_age
      real(dp) :: central_conv_work, total_central_conv_work
      real(dp) :: cumulative_nuc_neu

      ! these routines are called by the standard run_star check_model
      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         s% other_split_mix => SGC06_other_split_mix

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

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
         s% job% warn_run_star_extras =.false.       

      end subroutine extras_controls

      subroutine SGC06_other_split_mix( &
            id, dt_total, species, pass, num_passes, result_code, ierr)
         use const_def, only: dp
         integer, intent(in) :: id
         real(dp), intent(in) :: dt_total
         integer, intent(in) :: species, pass, num_passes
         integer, intent(out) :: result_code, ierr

         integer :: nz, i, j, k, m, ktop, kbot
         real(dp) :: dr, taujk, fjk, mconv
         real(dp), dimension(:), allocatable :: dtau_mix
         real(dp), dimension(:,:), allocatable :: dxa
         type (star_info), pointer :: s
         result_code = keep_going
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         nz = s% nz

         ! this is a quick-and-dirty implementation of the "advective"
         ! mixing algorithm is described in equations (4)-(7) of

         !     Straniero, O., Gallino, R., & Cristallo, S. (2006)
         !     s process in low-mass asymptotic giant branch stars
         !     http://adsabs.harvard.edu/abs/2006NuPhA.777.311S

         ! we want to calculate the change in abundances, so
         ! X_j - X0_j = (1/Mconv) sum_k (X0_k - X0_j) f_j,k dm_k

         ! allocate arrays
         allocate(dtau_mix(nz))
         allocate(dxa(species, nz))
         dxa = 0

         ! calculate dtau_mix(k)
         do k = 1, nz
            if (k == s% nz) then
               dr = s% r(k) - s% R_center
            else
               dr = s% r(k) - s% r(k+1)
            end if
            if (s% conv_vel(k) .gt. 1d-99) then
               dtau_mix(k) = dr / s% conv_vel(k)
            else
               dtau_mix(k) = 0
            end if
         end do

         ! loop over convection zones

         ! integer :: num_conv_boundaries
         !    ! boundaries of regions with mixing_type = convective_mixing
         !    ! boundaries are ordered from center to surface
         ! real(dp), pointer :: conv_bdy_q(:) ! (num_conv_boundaries)
         !    ! subcell location of boundary
         ! logical, pointer :: top_conv_bdy(:) ! (num_conv_boundaries)
         ! integer, pointer :: conv_bdy_loc(:) ! (num_conv_boundaries)
         !    ! if top_conv_bdy, top of region is between loc and loc+1
         !    ! else bottom of region is between loc and loc-1

         do m = 1, s% num_conv_boundaries

            if (s% top_conv_bdy(m)) then
               if (m == 1) then
                  kbot = nz
               else
                  kbot = s% conv_bdy_loc(m-1)
               endif
               ktop = s% conv_bdy_loc(m) + 1
            else
               cycle ! it was a bottom, so go get the top
            endif

            ! get mass in convective zone
            mconv = s% m(ktop) - s% m(kbot)

            do j = ktop, kbot ! loop over zone
               do i = 1, species ! loop over species
                  do k = ktop, kbot ! do k sum

                     if (j .eq. k) cycle ! dxa is 0 for j == k

                     ! calculate tau(j,k)
                     taujk = sum(dtau_mix(min(j,k):max(j,k)))

                     ! convert to f(j,k)
                     if (s% dt .gt. taujk) then
                        fjk = 1d0
                     else
                        fjk = dt_total / taujk
                     end if

                     ! add contribution to dxa (this is the k-sum)
                     dxa(i,j) = dxa(i,j) + &
                          (s% xa(i,k) - s% xa(i,j)) * fjk * s% dm(k)

                  end do

                  ! now divide by total convective zone mass
                  dxa(i,j) = dxa(i,j) / mconv

               end do
            end do

         end do

         ! adjust composition
         s% xa(1:species, 1:nz) = s% xa(1:species, 1:nz) + dxa(1:species, 1:nz)

         ! Bill says:
         ! Your other_split_mix is free to change the s% xa values as it pleases.
         ! Just be careful to conserve species (i.e., don't change sum(s% xa(j,1:nz)) for any j)
         ! and make the sum(s% xa(1:species,k)) = 1 for each k.

         ! in practice, even though the above is not guaranteed to do
         ! this to machine precision, it doesn't seem to be a problem

         deallocate(dxa, dtau_mix)

      end subroutine SGC06_other_split_mix


       function eval_conv_work(s, k) result(conv_work)

         use chem_def, only : ina23

         type (star_info), pointer :: s
         integer :: k
         real(dp) :: conv_work

         integer :: j
         real(dp) :: delta_X23, delta_mu, face_flux, delta_Ye

         j = s% net_iso(ina23)
         
         delta_X23 = s% xa(j,k) - s% xa(j,k-1)
         delta_mu = s% eta(k) * kerg * s%T (k) - s% eta(k-1) * kerg * s%T (k-1)
         face_flux = 4 * pi * s% r(k) * s% r(k) * s% rho_face(k) * s% mlt_vc(k)

         delta_Ye = s% ye(k) - s% ye(k-1)
         
         conv_work = s% sig(k) * (delta_Ye / amu) * delta_mu * s% dt

       end function eval_conv_work

       
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0

         total_central_conv_work = 0
         cumulative_nuc_neu = 0

         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if

         initial_c12_mass = dot_product(s% dm(1:s% nz), &
              s% xa(s% net_iso(ic12), 1:s% nz)) / Msun
         initial_age = s% star_age

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

         if (s% x_logical_ctrl(1)) then
            if (s% mass_conv_core > 0) extras_check_model = terminate
         endif

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
         how_many_extra_history_columns = 5
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: current_C12_mass, current_age, mdot, surface_c12, accreted_C12_mass
         integer :: k

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         names(1) = 'center_neutron_excess'
         vals(1) = 1d0 - 2d0 * s% ye(s% nz)

         current_C12_mass = dot_product(s% dm(1:s% nz), &
              s% xa(s% net_iso(ic12), 1:s% nz)) / Msun
         current_age = s% star_age
         surface_c12 = s% xa(s% net_iso(ic12), s% nz)
         mdot = s% mass_change
         accreted_C12_mass = surface_c12 * mdot * (s% star_age - initial_age)

         names(2) = 'carbon_mass_burned'
         vals(2) = -(current_C12_mass - initial_C12_mass - accreted_C12_mass)
         write(*,*) trim(names(2)), vals(2)

         names(3) = 'central_conv_work'
         vals(3) = central_conv_work
         write(*,*) trim(names(3)), ' = ', vals(3)

         names(4) = 'total_central_conv_work'
         vals(4) = total_central_conv_work
         write(*,*) trim(names(4)), ' = ', vals(4)

         names(5) = 'cumulative_nuc_neu'
         vals(5) = cumulative_nuc_neu
         
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 2
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


         names(1) = 'conv_work'
         vals(1,1) = 0
         do k = 2, nz
            vals(k,1) = eval_conv_work(s, k)
         end do

         names(2) = 'int_conv_work'
         vals(nz,2) = vals(nz,1)
         do k = nz - 1, 1, -1
            vals(k,2) = vals(k+1,2) + vals(k,1)
         end do

      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! central temperature (T8)
         s% xtra1 = s% T(s% nz) / 1e8

         if ((floor(s% xtra1_old) - floor(s% xtra1) .ne. 0)) then

            ! save a profile & update the history
            s% need_to_update_history_now = .true.
            s% need_to_save_profiles_now = .true.

         endif

         ! loop & sum over central convection zone
         central_conv_work = 0
         if (s% top_conv_bdy(1)) then
            do k = s% conv_bdy_loc(1), s% nz
               central_conv_work = central_conv_work + eval_conv_work(s, k)
            end do
         end if

         total_central_conv_work = total_central_conv_work + central_conv_work

         cumulative_nuc_neu = cumulative_nuc_neu + &
              dot_product(s% dm(1:s% nz),s% eps_nuc_neu_total(1:s% nz)) * s% dt
         
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
