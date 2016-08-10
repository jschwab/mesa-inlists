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
      use interp_2d_lib_sg
      implicit none

      logical :: kap_retry_flag
      
      real(dp), parameter :: logT_blend = 4.15 !3.85
      real(dp), parameter :: delta_logT_blend = 0.05
      
      integer, parameter :: kk_nr = 46 ! for (logR=-8.0; logR<=1.0; logR += 0.2)
      integer, parameter :: kk_nt = 43 ! for (logT=3.325; logT<=4.40;logT += 0.025)

      real, dimension(kk_nr) :: kk_r
      real, dimension(kk_nt) :: kk_t
      real, dimension(kk_nr, kk_nt) :: kk_k

      real, pointer :: wk1(:) ! 4*kk_nr*kk_nt
      real, pointer, dimension(:, :, :) :: wk ! 4,kk_nr,kk_nt

      ! maybe I should use not-a-knot?! (=0)
      integer, parameter :: ibcxmin=4, ibcxmax=4, ibcymin=4, ibcymax=4

      real, dimension(kk_nr) :: bcxmin, bcxmax
      real, dimension(kk_nt) :: bcymin, bcymax
      integer :: ilinx, iliny
      
     ! subroutine interp_mkbicub_sg(x,nx,y,ny,f1,nf2,
     ! >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     ! >   ibcymin,bcymin,ibcymax,bcymax,
     ! >   ilinx,iliny,ier)

     !     use bicub_sg

     !     integer, intent(in) :: nx                        ! length of x vector
     !     integer, intent(in) :: ny                        ! length of y vector
     !     real, intent(in) :: x(:) ! (nx)                        ! x vector, strict ascending
     !     real, intent(in) :: y(:) ! (ny)                        ! y vector, strict ascending
     !     integer, intent(in) :: nf2                       ! 2nd dimension of f, nf2.ge.nx
     !     real, intent(inout), pointer :: f1(:) ! =(4,nf2,ny)               ! data & spline coefficients




      ! these routines are called by the standard run_star check_model
       contains

         real(dp) function interp_val_to_pt(v,k,sz,dq,str)
           use interp_1d_lib, only: interp_4_to_1
           integer, intent(in) :: k, sz
           real(dp), pointer :: v(:), dq(:)
           character (len=*), intent(in) :: str
           integer :: ierr
           include 'formats'
           if (k == 1) then
              interp_val_to_pt = v(k)
              return
           end if
           if (k > 2 .and. k < sz) then
              ierr = 0
              call interp_4_to_1( &
                   0.5d0*(dq(k-2)+dq(k-1)), &
                   0.5d0*(dq(k-1)+dq(k)), &
                   0.5d0*(dq(k)+dq(k+1)), &
                   0.5d0*dq(k-2)+dq(k-1), &
                   v(k-2), v(k-1), v(k), v(k+1), &
                   interp_val_to_pt, str, ierr)
              if (ierr == 0) return
              write(*,1) '0.5d0*(dq(k-2)+dq(k-1))', 0.5d0*(dq(k-2)+dq(k-1))
              write(*,1) '0.5d0*(dq(k-1)+dq(k))', 0.5d0*(dq(k-1)+dq(k))
              write(*,1) '0.5d0*(dq(k)+dq(k+1))', 0.5d0*(dq(k)+dq(k+1))
              write(*,2) 'dq(k-2)', k-2, dq(k-2)
              write(*,2) 'dq(k-1)', k-1, dq(k-1)
              write(*,2) 'dq(k)', k, dq(k)
              write(*,2) 'dq(k+1)', k+1, dq(k+1)

              stop 'interp_val_to_pt'
           endif
           interp_val_to_pt = (v(k)*dq(k-1) + v(k-1)*dq(k))/(dq(k-1) + dq(k))
         end function interp_val_to_pt

         real(dp) function get_L_rad(s,k)
           type (star_info), pointer :: s         
           integer, intent(in) :: k
           integer :: j
           real(dp) :: kap_face, del_m, del_T4
           if (k == 1) then
              j = 2
           else
              j = k
           end if
           kap_face = interp_val_to_pt(s% opacity,j,s% nz,s% dq,'get_L_rad')
           del_m = 0.5d0*(s% dm(j-1) + s% dm(j))
           del_T4 = pow4(s% T(j-1)) - pow4(s% T(j))
           get_L_rad = -s% area(j)*s% area(j)*crad*clight/(3*kap_face)*(del_T4/del_m)
         end function get_L_rad

         integer function extras_startup(s, id, restart, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr

         integer :: ii, jj
         real :: logt, logr, logd, logkap
         ierr = 0
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if

         ! read tables
         open(unit = 13, file = "CO-table-with-metals.dat")
         do ii = 1, kk_nt
            do jj = 1, kk_nr
               read(13, *) logt, logr, logd, logkap
               kk_k(jj,ii) = logkap
               kk_r(jj) = logr
            end do
            kk_t(ii) = logt
         end do
         close(13)

         ! set up pointers
         allocate(wk(1:4, 1:kk_nr, 1:kk_nt), wk1(1:4*kk_nr*kk_nt))
         wk(1:4, 1:kk_nr, 1:kk_nt) => wk1(1:4*kk_nr*kk_nt)

         ! put tables in work array
         wk(1, 1:kk_nr, 1:kk_nt) = kk_k(1:kk_nr, 1:kk_nt)
         
         ! construct interpolant
         call interp_mkbicub_sg(kk_r, kk_nr, kk_t, kk_nt, wk1, kk_nr, &
              ibcxmin,bcxmin,ibcxmax,bcxmax, & 
              ibcymin,bcymin,ibcymax,bcymax, & 
              ilinx,iliny,ierr)

         if (ierr .ne. 0) stop "Failed while constructing interpolant"

         if (restart) then

            ! flip opacity on/off
            if (log10_cr(s% T(1)) .lt. logT_blend + 1.5 * delta_logT_blend ) then
               s% use_other_kap = .true.
            endif
            
            if (log10_cr(s% T(1)) .gt. logT_blend + 2.5 * delta_logT_blend ) then
               s% use_other_kap = .false.
            endif
            
            write(*,*) "restart: use_other_kap", s% use_other_kap

         endif

         kap_retry_flag = .false.

      end function extras_startup
      


      subroutine extras_controls(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

         s% other_kap_get_Type2 => cold_CO_kap_get_Type2

      end subroutine extras_controls

      subroutine kasen_CO(log10_rho, log10_T, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
        real(dp), intent(in) :: log10_rho
        real(dp), intent(in) :: log10_T
        real(dp), intent(out) :: kap
        real(dp), intent(out) :: dln_kap_dlnRho
        real(dp), intent(out) :: dln_kap_dlnT

        integer :: ierr
        
        real :: logd, logt, logr

        integer :: ict(6)
        real :: fval(6)

        logd = real(log10_Rho)
        logt = real(log10_T)
        logr = logd - 3*logT + 18

        if ((logr .gt. maxval(kk_r)) .or. (logr .lt. minval(kk_r))) then 
           write(*,*) "Off table in R: ", maxval(kk_r), "logr = ", logr, minval(kk_r)
           ierr = -1
           return
        end if
        if ((logt .gt. maxval(kk_t)) .or. (logt .lt. minval(kk_t))) then
           write(*,*) "Off table in T: ", maxval(kk_t), "logt = ", logt, minval(kk_t)
           ierr = -1
           return
        end if

        ! get f, df/dx, df/dy
        ict(1:3) = 1
        ict(4:6) = 0
        
        call interp_evbicub_sg(logr,logt,kk_r,kk_nr,kk_t,kk_nt,&
             ilinx,iliny,wk1,kk_nr,ict,fval,ierr)

        kap = 10d0 ** fval(1)
        dln_kap_dlnRho = fval(2)
        dln_kap_dlnT = fval(3) - 3.0 * fval(2)
        
      end subroutine kasen_CO

      subroutine cold_CO_kap_get_Type2( &
           id, k, handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, &
           log10_rho, log10_T, species, chem_id, net_iso, xa, &
           lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
           frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)

        use kap_lib, only: kap_get_Type2
        integer, intent(in) :: id
        integer, intent(in) :: k
        integer, intent(in) :: handle
        real(dp), intent(in) :: zbar
        real(dp), intent(in) :: X, Z, Zbase, XC, XN, XO, XNe
        real(dp), intent(in) :: log10_rho
        real(dp), intent(in) :: log10_T
        double precision, intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
        integer, intent(in) :: species
        integer, pointer :: chem_id(:)
        integer, pointer :: net_iso(:)
        real(dp), intent(in) :: xa(:)
        real(dp), intent(out) :: frac_Type2
        real(dp), intent(out) :: kap
        real(dp), intent(out) :: dln_kap_dlnRho
        real(dp), intent(out) :: dln_kap_dlnT
        integer, intent(out) :: ierr

        real(dp) :: log10_R
        
        real(dp) :: w_kap, w_dln_kap_dlnRho, w_dln_kap_dlnT
        real(dp) :: c_kap, c_dln_kap_dlnRho, c_dln_kap_dlnT
        real(dp) :: alfa, beta, logT, lower_bdy, upper_bdy

        logical :: dbg = .false.
        real(dp), parameter :: kap_min = 1e-6
                
        type (star_info), pointer :: s
        call get_star_ptr(id, s, ierr)

        ! aliases
        logT = log10_T
        lower_bdy = logT_blend - delta_logT_blend
        upper_bdy = logT_blend + delta_logT_blend

        alfa = (logT - lower_bdy) / (upper_bdy - lower_bdy)
        alfa = max(min(alfa,1d0),0d0) ! force in [0,1]
        beta = 1d0 - alfa

        if (logT >= lower_bdy) then
           call kap_get_Type2( &
                handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, &
                log10_rho, log10_T, &
                lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                frac_Type2, w_kap, w_dln_kap_dlnRho, w_dln_kap_dlnT, ierr)
        else
           w_kap = 0
           w_dln_kap_dlnRho = 0
           w_dln_kap_dlnT = 0
        end if

        if (logT <= upper_bdy) then

           call kasen_CO(log10_rho, log10_T, &
                c_kap, c_dln_kap_dlnRho, c_dln_kap_dlnT, ierr)

           if (ierr .lt. 0) then 
              write(*,*) "k =", k
              write(*,*) "X =", X
              write(*,*) "Z =", Z
              write(*,*) "Zbase =", Zbase
              write(*,*) "XC =", XC
              write(*,*) "XN =", XN
              write(*,*) "XO =", XO
              write(*,*) "XNe =", XNe
              write(*,*) "log10_rho =", log10_rho
              write(*,*) "log10_T =", log10_T

              write(*,*) "c_kap =", c_kap
              write(*,*) "c_dln_kap_dlnRho=", c_dln_kap_dlnRho
              write(*,*) "c_dln_kap_dlnT=", c_dln_kap_dlnT

              if (dbg) then 
                 stop
              else
                 kap_retry_flag = .true.

                 c_kap = kap_min
                 c_dln_kap_dlnRho = 0
                 c_dln_kap_dlnT = 0

              end if
           endif

        else
           c_kap = 0
           c_dln_kap_dlnRho = 0
           c_dln_kap_dlnT = 0
        end if


        kap = alfa*w_kap + beta*c_kap
        dln_kap_dlnRho = (alfa*w_kap*w_dln_kap_dlnRho + beta*c_kap*c_dln_kap_dlnRho)/kap
        dln_kap_dlnT = (alfa*w_kap*w_dln_kap_dlnT + beta*c_kap*c_dln_kap_dlnT)/kap

      end subroutine cold_CO_kap_get_Type2
      
      
      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(s, id, id_extra)

        use chem_def, only : i_burn_ne, i_burn_si
        use rates_def, only : i_rate

        integer :: k
        real(dp) :: runaway
        
        type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         extras_check_model = keep_going         

         ! retry bad steps where we went off the kap table
         if (kap_retry_flag) then
            extras_check_model = retry
            kap_retry_flag = .false.
         endif

         ! stop at ne burning
         if (s% x_logical_ctrl(1)) then 
            ! there are two conditions where we want to stop
            termination_code_str(t_xtra1) = 'started ne burning'

            k = s% max_eps_nuc_k
            runaway = s% eps_nuc(k) - s% non_nuc_neu(k)

            if (runaway .gt. 0) then

               ! check if max energy is neon burning
               if (maxloc(s% eps_nuc_categories(i_rate,:,k),1) == i_burn_ne) then
                  extras_check_model = terminate
                  s% termination_code = t_xtra1
               end if
            end if
         end if

         ! stop at si burning
         if (s% x_logical_ctrl(2)) then 
            ! there are two conditions where we want to stop
            termination_code_str(t_xtra2) = 'started si burning'

            k = s% max_eps_nuc_k
            runaway = s% eps_nuc(k) - s% non_nuc_neu(k)
            
            if (runaway .gt. 0) then

               ! check if max energy is si burning
               if (maxloc(s% eps_nuc_categories(i_rate,:,k),1) == i_burn_si) then
                  extras_check_model = terminate
                  s% termination_code = t_xtra2
               end if
            end if
         end if

         
         ! write(*,*) "logT: ", kk_t
         ! write(*,*) "logR: ", kk_r
         ! stop

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


      integer function how_many_extra_history_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         how_many_extra_history_columns = 39
      end function how_many_extra_history_columns
      
      subroutine data_for_extra_history_columns(s, id, id_extra, n, names, vals, ierr)
        type (star_info), pointer :: s
        integer, intent(in) :: id, id_extra, n
        character (len=maxlen_history_column_name) :: names(n)
        real(dp) :: vals(n)
        integer, intent(out) :: ierr

        integer :: k, k_heat, k_peak, k_burn, k_conv_lower, k_conv_upper
        real(dp) :: maxT, Lrad, Lconv, L_conv_max, L_conv_avg, n1, d1

        logical, dimension(:), allocatable :: my_mask

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

        ! practice safe indexing...badly
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


        names(1+34) = 'first_max_T_lgT'
        names(2+34) = 'first_max_T_lgRho'
        names(3+34) = 'first_max_T_m'

        maxT = 0
        do k = s% nz, 1, -1
           if (s% T(k) .lt. maxT) then
              vals(1+34) = log10_cr(s% T(k))
              vals(2+34) = log10_cr(s% Rho(k))
              vals(3+34) = s% m(k) / Msun
              exit
           end if
           maxT = max(maxT, s% T(k)) 
        enddo

        names(4+34) = 'L_conv_avg'
        names(5+34) = 'L_conv_max'

        L_conv_max = 0
        n1 = 0
        d1 = 0
        do k = 1, s% nz
           Lrad = get_L_rad(s,k)
           Lconv = (s% L(k) - Lrad)/Lsun
           L_conv_max = max(L_conv_max, Lconv)
           if (Lconv .gt. 0) then
              n1 = n1 + s% dq(k) * Lconv
              d1 = d1 + s% dq(k)
           endif
        enddo
        L_conv_avg = n1!/d1

        vals(4+34) = L_conv_avg
        vals(5+34) = L_conv_max

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(s, id, id_extra, n, nz, names, vals, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
         
         !note: do NOT add these names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer :: ierr
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! flip opacity on/off
         if (log10_cr(s% T(1)) .lt. logT_blend + 1.5 * delta_logT_blend ) then
            s% use_other_kap = .true.
         endif

         if (log10_cr(s% T(1)) .gt. logT_blend + 2.5 * delta_logT_blend ) then
            s% use_other_kap = .false.
         endif

         write(*,*) "use_other_kap", s% use_other_kap
         
         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
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
      
