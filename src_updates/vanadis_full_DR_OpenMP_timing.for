!
! VANADIS.FOR - VER. 1.07_6  07.05.2026 - full DR
!
! autor : Marek Chodorski
! e-mail: marek_ac@wp.pl
!
!
! Vanadis 3D - A three-dimensional dynamic Eulerian finite-element model for atmospheric pollutant transport
!
!
! In the program I used the following procedure:
!
!            https://www.netlib.org/slatec/lin/dbcg.f
!
!            Preconditioned BiConjugate Gradient Sparse Ax = b Solver.
!            Routine to solve a Non-Symmetric linear system  Ax = b
!            using the Preconditioned BiConjugate Gradient method.
!            library SLATEC (SLAP)
!            author Greenbaum, Anne, (Courant Institute)
!            Seager, Mark K., (LLNL)
!
      program vanadis

         parameter (nax = 7, nay = 20, naz = 20)
         parameter (max_time = 2000)
         parameter (i_time_step = 100)

         implicit real*8 (a - h, o - z)
         dimension dx(nax), dy(nay), dz(naz)
         dimension xx(nax * nay * naz), yy(nax * nay * naz),
     &             zz(nax * nay * naz)
         dimension t(6)
         dimension iss((nax - 1) * (nay - 1) * (naz - 1), 8)
         dimension kod((nax - 1) * (nay - 1) * (naz - 1), 6)
         dimension bb((nax - 1) * (nay - 1) * (naz - 1) * 64)
         dimension tnw2(nax * nay * naz)
		 
		 dimension kod2((nax - 1) * (nay - 1) * (naz - 1), 6)
		 
		 real*8 mass

         common /d_time/ i_d_czas

         il_iteracji = 1000
         i_d_czas = i_time_step
         ttol = 1d-8

         call copyright
         call read_cfg
         call par_3d
         call wsp_na_o (nax, nay, naz, dx, dy, dz)
         call siatka (xx, yy, zz, dx, dy, dz, nax, nay, naz)
         call w_el (iss, nax, nay, naz)

         call kod_el_dr(nax, nay, naz, kod)
		 
         i_czas = 0

         do while (i_czas .lt. max_time)

            call przest3d(i_czas, xx, yy, zz, iss, kod, t, nax, nay,
     &                    naz, bb, il_iteracji, ttol, tnw2)
	 
c------------------------------------------------------------
c  LICZENIE MASY PO KAŻDYM KROKU CZASOWYM
c------------------------------------------------------------
            nelem = (nax-1)*(nay-1)*(naz-1)
		    
            call integrate_mass_hex8(nelem, iss, xx, yy, zz,
     &                               tnw2, mass)
		    
            write(6,*) 'czas=', i_czas, ' masa=', mass	 

         enddo

         call vanadis_files(nax, nay, naz, max_time, i_time_step)

      end

      subroutine wsp_na_o(nax, nay, naz, dx, dy, dz)

         implicit real*8 (a - h, o - z)
         dimension dx(*), dy(*), dz(*)

         common /cal_el/ ilosc_el, ilosc_wezlow, i_ot_elem(8)
!$omp threadprivate(/cal_el/)

         ilosc_el = (nax - 1) * (nay - 1) * (naz - 1)
         ilosc_wezlow = nax * nay * naz

         dx(1) = 0.0d+00
         dx(2) = 20.0d+00
         dx(3) = 35.0d+00
         dx(4) = 85.0d+00
         dx(5) = 135.0d+00
         dx(6) = 230.0d+00
         dx(7) = 420.0d+00
   
         dy(1) = -4000.0d+00
         dy(2) = -1500.0d+00
         dy(3) = -900.0d+00
         dy(4) = -700.0d+00
         dy(5) = -500.0d+00
         dy(6) = -400.0d+00
         dy(7) = -270.0d+00
         dy(8) = -180.0d+00
         dy(9) = -90.0d+00
         dy(10) = -25.0d+00
         dy(11) = 25.0d+00
         dy(12) = 90.0d+00
         dy(13) = 180.0d+00
         dy(14) = 270.0d+00
         dy(15) = 400.0d+00
         dy(16) = 500.0d+00
         dy(17) = 700.0d+00
         dy(18) = 900.0d+00
         dy(19) = 1500.0d+00
         dy(20) = 4000.0d+00
   
   	     dz(1) = -4000.0d+00
         dz(2) = -1500.0d+00
         dz(3) = -900.0d+00
         dz(4) = -700.0d+00
         dz(5) = -500.0d+00
         dz(6) = -400.0d+00
         dz(7) = -270.0d+00
         dz(8) = -180.0d+00
         dz(9) = -90.0d+00
         dz(10) = -25.0d+00
         dz(11) = 25.0d+00
         dz(12) = 90.0d+00
         dz(13) = 180.0d+00
         dz(14) = 270.0d+00
         dz(15) = 400.0d+00
         dz(16) = 500.0d+00
         dz(17) = 700.0d+00
         dz(18) = 900.0d+00
         dz(19) = 1500.0d+00
         dz(20) = 4000.0d+00

         return
      end

      subroutine licz_alfa (i_czas, nr_el, kod, alfa, qs, nax, nay,
     &                      naz)

         implicit real*8 (a - h, o - z)
         dimension kod((nax - 1) * (nay - 1) * (naz - 1), 6)
         real * 8 q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg, alpha_cfg
         common /cfg/ q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg,
     &                alpha_cfg
!
!      QS=-ALFA*TF    ! TF-TEMP.PLYNU OMYWAJACEGO POWIERZCHNIE
!
!      ALFA=0.001D+00 ! WSPOLCZYNNIK WNIKANIA - PARAMETR
!
         alfa = alpha_cfg
         qs = 0.0d+00     ! STRUMIEN CIEPLA - 0 DLA WARUNKU BRZEGOWEGO III RODZAJU (W NASZYM PRZYPADKU)

         return
      end

      subroutine atmosfera (i_czas, nr_el, lambda, pzanik, v, t, q_in)

         implicit real*8 (a - h, o - z)
         real*8 lambda, obj, dlx, dly, dlz
         real*8 q_in
         dimension t(6), v(3), lambda(3)
         common /wsp_lok/ xxe(8), yye(8), zze(8)
!$omp threadprivate(/wsp_lok/)
         real * 8 q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg, alpha_cfg
         common /cfg/ q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg,
     &                alpha_cfg

         t(1) = 0.0d0
         t(2) = 0.0d0
         t(3) = 0.0d0
         t(4) = 0.0d0
         t(5) = 0.0d0
         t(6) = 0.0d0

         lambda(1) = k_cfg
         lambda(2) = lambda(1)
         lambda(3) = lambda(1)

         pzanik = p_cfg

! here is a place for pollution source in element nr_el

         q_in = 0.0d+00                    ! kg/s

         if (nr_el .eq. 1083) q_in = q_cfg ! kg/s

         v(1) = v1_cfg
         v(2) = v2_cfg
         v(3) = v3_cfg

         return
      end

      subroutine siatka (xx, yy, zz, dx, dy, dz, nax, nay, naz)

         implicit real*8 (a - h, o - z)
         dimension xx(*), yy(*), zz(*), dx(*), dy(*), dz(*)

         i = 1
         do jz = 1, naz
            do jy = 1, nay
               do jx = 1, nax
                  xx(i) = dx(jx)
                  yy(i) = dy(jy)
                  zz(i) = dz(jz)
                  i = i + 1
               enddo
            enddo
         enddo

         return
      end

      subroutine w_el (iss, nax, nay, naz)
         dimension iss((nax - 1) * (nay - 1) * (naz - 1), 8)

         k = 1
         ik = 0
         ik2 = nay * nax

         do j = 1, (nax - 1) * (nay - 1) * (naz - 1)

            iss(j, 1) = k + ik2 * ik
            iss(j, 2) = k + 1 + ik2 * ik
            iss(j, 3) = k + 1 + nax + ik2 * ik
            iss(j, 4) = k + nax + ik2 * ik
            iss(j, 5) = k + nax * nay + ik2 * ik
            iss(j, 6) = k + 1 + nax * nay + ik2 * ik
            iss(j, 7) = k + 1 + nax * nay + nax + ik2 * ik
            iss(j, 8) = k + nax * nay + nax + ik2 * ik

            k = k + 1

            if (mod(j, (nax - 1)) .eq. 0) k = k + 1

            if (mod(j, (nax - 1) * (nay - 1)) .eq. 0) then
               k = 1
               ik = ik + 1
            endif

         enddo

         return
      end
c
      subroutine kod_el_dr(nax, nay, naz, kod)
         dimension kod((nax - 1) * (nay - 1) * (naz - 1), 6)

         il_el = (nax - 1) * (nay - 1) * (naz - 1)

         do i = 1, il_el
            do j = 1, 6
               kod(i, j) = 0
            enddo
         enddo

         do i = 1, ((nax - 1) * (nay - 1))
            kod(i, 1) = 1
         enddo

         do i = il_el - (nax - 1) * (nay - 1) + 1, il_el
            kod(i, 2) = 2
         enddo

         do i = (nax - 1), il_el, (nax - 1)
            kod(i, 4) = 4
         enddo

         do i = 1, il_el, nax - 1
            kod(i, 3) = 3
         enddo

         do k = 0, naz - 2
            do i = k * (nax - 1) * (nay - 1) + 1,
     &             k * (nax - 1) * (nay - 1) + nax - 1
               kod(i, 6) = 6
            enddo
         enddo

         do k = 1, naz - 1
            do i = k * (nax - 1) * (nay - 1) - nax + 1 + 1,
     &             k * (nax - 1) * (nay - 1)
               kod(i, 5) = 5
            enddo
         enddo

         return
      end
c
      subroutine par_3d

         implicit real*8 (a - h, o - z)
         real*8 ksi, eta, zeta
         real*8 ksi_g, eta_g, zeta_g

         common /lok/ ksi(8), eta(8), zeta(8)
         common /lok_g/ ksi_g(8), eta_g(8), zeta_g(8)
         common /pow_el/ i_pow(6, 4)

         ksi(1) = -1.0d+00
         ksi(2) = 1.0d+00
         ksi(3) = 1.0d+00
         ksi(4) = -1.0d+00
         ksi(5) = -1.0d+00
         ksi(6) = 1.0d+00
         ksi(7) = 1.0d+00
         ksi(8) = -1.0d+00

         eta(1) = -1.0d+00
         eta(2) = -1.0d+00
         eta(3) = 1.0d+00
         eta(4) = 1.0d+00
         eta(5) = -1.0d+00
         eta(6) = -1.0d+00
         eta(7) = 1.0d+00
         eta(8) = 1.0d+00

         zeta(1) = -1.0d+00
         zeta(2) = -1.0d+00
         zeta(3) = -1.0d+00
         zeta(4) = -1.0d+00
         zeta(5) = 1.0d+00
         zeta(6) = 1.0d+00
         zeta(7) = 1.0d+00
         zeta(8) = 1.0d+00

         do j = 1, 8
            ksi_g(j) = 0.577350269189626d0 * ksi(j)
            eta_g(j) = 0.577350269189626d0 * eta(j)
            zeta_g(j) = 0.577350269189626d0 * zeta(j)
         enddo

         i_pow(1, 1) = 1
         i_pow(1, 2) = 2
         i_pow(1, 3) = 3
         i_pow(1, 4) = 4

         i_pow(2, 1) = 5
         i_pow(2, 2) = 6
         i_pow(2, 3) = 7
         i_pow(2, 4) = 8

         i_pow(3, 1) = 1
         i_pow(3, 2) = 4
         i_pow(3, 3) = 8
         i_pow(3, 4) = 5

         i_pow(4, 1) = 2
         i_pow(4, 2) = 3
         i_pow(4, 3) = 7
         i_pow(4, 4) = 6

         i_pow(5, 1) = 4
         i_pow(5, 2) = 3
         i_pow(5, 3) = 7
         i_pow(5, 4) = 8

         i_pow(6, 1) = 1
         i_pow(6, 2) = 2
         i_pow(6, 3) = 6
         i_pow(6, 4) = 5

         return
      end

      subroutine data3d (kel, iss, xx, yy, zz, nax, nay, naz)

         implicit real*8 (a - h, o - z)
         dimension iss((nax - 1) * (nay - 1) * (naz - 1), 8)
         dimension xx(*), yy(*), zz(*)

!    KEL      - NUMER ELEMENTU

         common /wsp_lok/ xxe(8), yye(8), zze(8)
!$omp threadprivate(/wsp_lok/)
         common /cal_el/ ilosc_el, ilosc_wezlow, i_ot_elem(8)
!$omp threadprivate(/cal_el/)

         do j = 1, 8
            i_ot_elem(j) = iss(kel, j)
            xxe(j) = xx(i_ot_elem(j))
            yye(j) = yy(i_ot_elem(j))
            zze(j) = zz(i_ot_elem(j))
         enddo

         return
      end

      subroutine m_m(a, b, c, n, m)

         real*8 a(*), b(*), c(n, m)
         integer*4 i, j, n, m

         do i = 1, n
            do j = 1, m
               c(i, j) = 0.0d0
            enddo
         enddo

         do i = 1, n
            do j = 1, m
               c(i, j) = c(i, j) + a(i) * b(j)
            enddo
         enddo

         return
      end

      subroutine ilo_m_v(a, b, n, c)

         implicit real*8 (a - h, o - z)
         real*8 a(n, n), b(n), c(n)

         do i = 1, n
            c(i) = 0.0d0
            do k = 1, n
               c(i) = c(i) + a(i, k) * b(k)
            enddo
         enddo

         return
      end

      subroutine el3d (i_czas, kel1, t, kod, nax, nay, naz)

         implicit real*8 (a - h, o - z)

!    KEL1 -NUMER ELEMENTU

         real*8 lambda
         dimension t(*)
         dimension kod((nax - 1) * (nay - 1) * (naz - 1), 6)
         dimension zx1(4), zx2(4), zx3(4)
         dimension lambda(3), v(3)
         real*8 n1(8), dnx1(8), dny1(8), dnz1(8), djac
         real * 8 detjj
         real*8 dnxw(8), dnyw(8), dnzw(8)
         real*8 n1b(8), w(8)
         real*8 ksi, eta, zeta
         real*8 ksi_g, eta_g, zeta_g
         real * 8 detj(4), ze2_g(4), ze3_g(4), detjs, w1(8)

         real*8 dnx2(8), dny2(8), dnz2(8), n2(8)

         real*8 dnj(24),jak(9),detjq
		 
         real*8 dn_dxi(3, 8),dw_dxi(3, 8)
		 
		 real*8 vol

         real*8 gp2(2), gw2(2)
         data gp2 / -0.577350269189626d0, 0.577350269189626d0 /
         data gw2 / 1.d0, 1.d0 /

         common /lok_g/ ksi_g(8), eta_g(8), zeta_g(8)
         common /wsp_lok/ xxe(8), yye(8), zze(8)
!$omp threadprivate(/wsp_lok/)
         common /pow_el/ i_pow(6, 4)
         common /skl_el/ hma(8, 8), cma(8, 8), fmv(8), cg(8, 8),
     &                   hx(8,8), hy(8, 8), hz(8, 8), hvx(8, 8),
     &                   hvy(8, 8), hvz(8, 8), hmav(8, 8), cma_n(8, 8)
!$omp threadprivate(/skl_el/)

         keltmp = kel1

         call atmosfera (i_czas, keltmp, lambda, pzanik, v, t, q_in)

         do ji = 1, 8
            fmv(ji) = 0.0d+00
            do jj = 1, 8
               hma(ji, jj) = 0.0d+00
               cma(ji, jj) = 0.0d+00
               hmav(ji, jj) = 0.0d+00
               cma_n(ji, jj) = 0.0d+00
            enddo
         enddo

         vol = 0.0d0
	     
         do jl = 1, 8
	     
            zss = ksi_g(jl)
            ztt = eta_g(jl)
            zuu = zeta_g(jl)
	     
            call gauss_point_hex8DR(zss, ztt, zuu, xxe, yye, zze, djac,
     &           lambda, v, n1, w, dnx1, dny1, dnz1, dnxw, dnyw, dnzw)
	     
            vol = vol + djac

         enddo
		 
		 qv = q_in / vol

         do 1000 jl = 1, 8

            zss = ksi_g(jl)
            ztt = eta_g(jl)
            zuu = zeta_g(jl)

            call gauss_point_hex8DR(zss, ztt, zuu, xxe, yye, zze, djac,
     &           lambda,v, n1, w, dnx1, dny1, dnz1, dnxw, dnyw, dnzw)

            call m_m(dnxw, dnx1, hx, 8, 8)
            call m_m(dnyw, dny1, hy, 8, 8)
            call m_m(dnzw, dnz1, hz, 8, 8)

            call m_m(w, n1, cg, 8, 8)
            call m_m(w, dnx1, hvx, 8, 8)
            call m_m(w, dny1, hvy, 8, 8)
            call m_m(w, dnz1, hvz, 8, 8)

            z_nic1 = djac * lambda(1)
            z_nic2 = djac * lambda(2)
            z_nic3 = djac * lambda(3)

! ZRODLO CIEPLA  W ELEMENCIE

            do i = 1, 8
               fmv(i) = fmv(i) + qv * w(i) * djac
            enddo

            do i = 1, 8
               do j = 1, 8

                  hma(i, j) = hma(i, j) + hx(i, j) * z_nic1 + hy(i, j)
     &                        * z_nic2 + hz(i, j) * z_nic3
                  hmav(i, j) = hmav(i, j) + (hvx(i, j) * v(1)
     &                       + hvy(i,j) * v(2) +hvz(i, j) * v(3)) * djac
                  cma(i, j) = cma(i, j) + cg(i, j) * djac * pzanik

! DO STANU NIESTAC
                  cma_n(i, j) = cma_n(i, j) + cg(i, j) * djac

               enddo
            enddo

1000     continue

!
! DO WARUNKU BRZEGOWEGO III RODZAJU (TYLKO DLA POWIERZCHNI-"I_WB")
!
! ---------------------------------------------------------------
!  Warunek brzegowy III rodzaju na dolnej ścianie (xi = -1)
! ---------------------------------------------------------------

         i_wb = 3     ! dolna ściana - tak jak ja widzi uzytkownik
	 
         if (kod(kel1, i_wb) .eq. i_wb) then

            ! --- współrzędne węzłów ściany (4 węzły) ---
            do i2 = 1, 4
               zx1(i2) = xxe(i_pow(i_wb, i2))
               zx2(i2) = yye(i_pow(i_wb, i2))
               zx3(i2) = zze(i_pow(i_wb, i2))
            enddo
         
            do ig1 = 1, 2
               do ig2 = 1, 2
         
                  ze1 = -1.0d0          ! dolna ściana HEX8
                  ze2 = gp2(ig1)
                  ze3 = gp2(ig2)

                  call quad4_3d(zx1, zx2, zx3, ze2, ze3, detjs)

                  call gauss_point_hex8DR(ze1, ze2, ze3,
     &                   xxe, yye, zze,
     &                   djac, lambda, v, n1, w, dnx1, dny1, dnz1,
     &                   dnxw, dnyw, dnzw)

                  call m_m(w, n1, cg, 8, 8)

!                  call shape_hex8dr(ze1, ze2, ze3, n1b, dn_dxi, w,
!     &                              dw_dxi, lambda, v)
!         
!                  call m_m(w, n1b, cg, 8, 8)

                  call licz_alfa(i_czas, keltmp, kod, alfaa, qs,
     &                           nax, nay, naz)
         
                  z_nic = detjs * gw2(ig1) * gw2(ig2) * alfaa

                  do i5 = 1, 8
                     fmv(i5) = fmv(i5) - qs * w(i5) * detjs
                     do j5 = 1, 8
                        cma(i5, j5) = cma(i5, j5) + cg(i5, j5) * z_nic
                     enddo
                  enddo
         
               enddo
            enddo
         
         endif

         do k = 1, 8
            do l = 1, 8
               hma(k, l) = hma(k, l) + cma(k, l) + hmav(k, l)
            enddo
         enddo

         return
      end

      subroutine przest3d (i_czas, xx, yy, zz, iss, kod, t, nax, nay,
     &                     naz, bb, il_iteracji, ttol, tnw2)

         implicit real*8 (a - h, o - z)

         dimension xx(*), yy(*), zz(*), t(*), bb(*)
         dimension iss((nax - 1) * (nay - 1) * (naz - 1), 8)
         dimension kod((nax - 1) * (nay - 1) * (naz - 1), 6)

         dimension tnw2(nax * nay * naz), hma_t(8, 8), fmv_t(8),
     &                  wekt_t(8)

         character*20 wynik_t

         common /pow_el/ i_pow(6, 4)
         common /skl_el/ hma(8, 8), cma(8, 8), fmv(8), cg(8, 8),
     &                   hx(8,8), hy(8, 8), hz(8, 8), hvx(8, 8),
     &                   hvy(8, 8), hvz(8, 8), hmav(8, 8), cma_n(8, 8)
!$omp threadprivate(/skl_el/)
         common /cal_el/ ilosc_el, ilosc_wezlow, i_ot_elem(8)
!$omp threadprivate(/cal_el/)

         real * 8 hma2d(8, 8), fmv2d(8)
         real*8 t_loc(6)
         real*8 tmont1, tmont2, tsolv1, tsolv2
         real*8 omp_get_wtime

         real*8 fmv3(nax * nay * naz)
         real*8 diag(nax * nay * naz)
         integer*4 k, l
		 
		 logical is_diagonally_dominant
		 
		 integer*4 bc_dirichlet(nax * nay * naz)

         common /d_time/ i_d_czas

         if (i_czas .eq. 0) then
            do ii2 = 1, ilosc_wezlow
               tnw2(ii2) = 0.d0
            enddo
         endif

         i_czas = i_czas + i_d_czas

         do i = 1, ilosc_wezlow
            diag(i) = 0.0d+00
         enddo

         do i = 1, (64 * ilosc_el)
            bb(i) = 0.0d+00
         enddo

         do i = 1, ilosc_wezlow
            fmv3(i) = 0.0d+00
         enddo

         print *
         write(*, 199) ilosc_el, ilosc_wezlow, i_czas
 199     format('matrix/vector creation...',' nr of elements: ',i6,
     &          '  nodes: ',i7, ' time: ', i6,' sec')

         tmont1 = omp_get_wtime()

!$omp parallel do default(shared)
!$omp& private(jk,i,j,i4,k,l,iin,hma_t,wekt_t,fmv_t,
!$omp& fmv2d,hma2d,t_loc)
!$omp& reduction(+:fmv3,diag)
!$omp& copyin(/cal_el/)
         do 100 jk = 1, ilosc_el

!     USUWANIE 3 ELEMENTOW - PRZESZKODA
!        IF (JK.EQ.739) GOTO 100
!        IF (JK.EQ.740) GOTO 100
!        IF (JK.EQ.741) GOTO 100

            call data3d (jk, iss, xx, yy, zz, nax, nay, naz)
            call el3d (i_czas, jk, t_loc, kod, nax, nay, naz)

            do i = 1, 8
               do j = 1, 8
                  hma_t(i, j) = hma(i, j) - 3.0d+00 / i_d_czas *
     &            cma_n(i, j)
               enddo
            enddo

            do i4 = 1, 8
               wekt_t(i4) = tnw2(i_ot_elem(i4))
            enddo

            call ilo_m_v(hma_t, wekt_t, 8, fmv_t)

!--------------------------------------------------------------------
! wypelnianie macierzy i wektora dla kroku czasowego
!
            do i = 1, 8
               fmv2d(i) = 3.0d+00 * fmv(i) - fmv_t(i)
            enddo

            do i = 1, 8
               do j = 1, 8
               hma2d(i, j) = 2.0d0 * hma(i, j) + 3.0d0 / i_d_czas *
     &         cma_n(i, j)
               enddo
            enddo
!
! ----------------------------------------
! podstawienie do GS

            iin = 64 * (jk - 1) + 1
            call sklad_wd (bb(iin), hma2d)

            do j = 1, 8
               fmv3(i_ot_elem(j)) = fmv3(i_ot_elem(j)) + fmv2d(j)
            enddo
! ---------------------------------------

            do i = 1, 8
               do j = 1, 8
                  k = i_ot_elem(i)
                  l = i_ot_elem(j)
                  if (k .eq. l) then
                     diag(k) = diag(k) + hma2d(i, j)
                  endif
               enddo
            enddo

 100     continue
!$omp end parallel do

         tmont2 = omp_get_wtime()
         print *, 'Czas montowania macierzy/wektora = ',
     &            tmont2 - tmont1, ' sekund'

         if (i_czas.eq. 1000)
     &       call ascii_heatmap_ppm (ilosc_el, ilosc_wezlow, iss, bb)
	 
	 	 print*,'Diagonally dominant matrix [T/F] ', 
     &       is_diagonally_dominant(ilosc_el, ilosc_wezlow, iss, bb)
	 
         print*,'wspolczynnik dominacji ',
     &    diagonal_dominance_factor(ilosc_el, ilosc_wezlow, iss, bb) 

         print *,('solving a system of equations...')

         tol = ttol
         il_i = il_iteracji
         non = ilosc_wezlow
		 
		 call build_bc_nodes(iss, kod, bc_dirichlet, nax, nay, naz)
		 
c		 
c--- narzucenie Dirichleta na poziomie układu równań
c
         do i = 1, ilosc_wezlow
            if (bc_dirichlet(i).eq.1) then
               fmv3(i)  = 0.0d0      ! prawa strona = 0
               diag(i)  = 1.0d0      ! preconditioner: M_ii = 1
               tnw2(i)  = 0.0d0      ! rozwiązanie na brzegu = 0
            endif
         enddo

         tsolv1 = omp_get_wtime()

         call dbcg_simple(non, fmv3, tnw2, tol, il_i, iter, err, bb,
     &                    iss, diag, nax, nay, naz, bc_dirichlet)

         tsolv2 = omp_get_wtime()
         print *, 'Czas solvera = ', tsolv2 - tsolv1, ' sekund'

!
! for CUDA change dbcg_simple --> kernel_dbcg_simple and "output\" --> "output/" ( Windows / Linux )
!
!        CALL kernel_dbcg_simple(NON,FMV3,TNW2,TOL,IL_I,ITER,ERR,BB,
!    &                           ISS,DIAG,NAX,NAY,NAZ)
!
         if (err.gt.ttol) then
		    print*, 'ERROR - stopping criterion not satisfied: ',iter, ttol, err
         else			
           write(*, 1111) iter, ttol, err
         endif
1111     format(' OK - number of iterations: ', i5,' --- ',
     &          ' desired convergence tolerance: ', e10.4,
     &          ' residual norm: ', e10.4, ' --- ')

         write(wynik_t, 105) i_czas
105      format('output\T_', i5.5, '.TXT')

         nout10 = 10
         open(nout10, file = wynik_t, status = 'UNKNOWN')
         do i = 1, ilosc_wezlow
            write(nout10, 122) xx(i), yy(i), zz(i), tnw2(i) * 1d6 !!! RESULTS [mg/m^3] !!!
         enddo
 122     format (e12.5, 2x, e12.5, 2x, e12.5, 2x, e12.5)
         close(nout10)

         return
      end

      subroutine matvec(xin, xout, iss, bb, nax, nay, naz) !matrix vector multiply operation  Y = A*X  given A and X

         implicit real*8 (a - h, o - z)
         dimension xin(*), xout(*)
         dimension iss((nax - 1) * (nay - 1) * (naz - 1), 8), bb(*)

         common /cal_el/ ilosc_el, ilosc_wezlow, i_ot_elem(8)
!$omp threadprivate(/cal_el/)

         do j = 1, ilosc_wezlow
            xout(j) = 0.0d0
         enddo

         ik = 1
         do je = 1, ilosc_el
            do j = 1, 8
               i_ot_elem(j) = iss(je, j)
            enddo

            do i = 1, 8
               do j = 1, 8
                  xout(i_ot_elem(i)) = xout(i_ot_elem(i))
     &            +bb(ik) * xin(i_ot_elem(j))
                  ik = ik + 1
               enddo
            enddo
         enddo

         return
      end

      subroutine mttvec(xin, xout, iss, bb, nax, nay, naz) !matrix transpose vector multiply y = A'*X given A and X

         implicit real*8 (a - h, o - z)
         dimension xin(*), xout(*), b3(8, 8)
         dimension iss((nax - 1) * (nay - 1) * (naz - 1), 8), bb(*)

         common /cal_el/ ilosc_el, ilosc_wezlow, i_ot_elem(8)
!$omp threadprivate(/cal_el/)

         do j = 1, ilosc_wezlow
            xout(j) = 0.0d0
         enddo

         do je = 1, ilosc_el
            do j = 1, 8
               i_ot_elem(j) = iss(je, j)
            enddo

            ik3 = 0
            do i = 1, 8
               ik3 = ik3 + 1
               do j = 1, 8
                  b3(i, j) = bb(j + 8 * (ik3 - 1) + 64 * (je - 1))
               enddo
            enddo

            do i = 1, 8
               do j = 1, 8
                  xout(i_ot_elem(i)) = xout(i_ot_elem(i)) + 
     &            b3(j, i) * xin(i_ot_elem(j))

               enddo
            enddo
         enddo

         return
      end

      subroutine sklad_wd (bb, hmad)

         implicit real*8 (a - h, o - z)
         dimension bb(*)
         real * 8 hmad(8, 8)

         common /skl_el/ hma(8, 8), cma(8, 8), fmv(8), cg(8, 8),
     &                   hx(8,8), hy(8, 8), hz(8, 8), hvx(8, 8),
     &                   hvy(8, 8), hvz(8, 8), hmav(8, 8), cma_n(8, 8)
!$omp threadprivate(/skl_el/)

         common / duza_t / hma2(64 * 3, 43 * 3), fmv2(64 * 3)

         common /d_time/ i_d_czas

         k = 1
         do i = 1, 8
            do j = 1, 8
               bb(k) = hmad(i, j)
               k = k + 1
            enddo
         enddo

         return
      end

      subroutine dbcg_simple(n, b, x, tol, itmax, iter, err, bb, iss,
     &                       diag,nax, nay, naz, bc_dirichlet)

         integer*4 n, iter, itmax, nax, nay, naz, i
         real*8 tol, err, b(*), bb(*), dnrm2
         real*8 diag(*), x(*)
         integer*4 iss((nax - 1) * (nay - 1) * (naz - 1), 8)
         real*8 ak, akden, bk, bkden, bknum, bnrm
         real*8 p(n), pp(n), r(n), rr(n), z(n), zz(n)
		 
		 integer*4 bc_dirichlet(*)

c         call matvec(x, r, iss, bb, nax, nay, naz)

         call matvec_0(x, r, iss, bb, nax, nay, naz, bc_dirichlet)

         do i = 1, n
            r(i) = b(i) - r(i)
            rr(i) = r(i)
         enddo
         call msolve(n, r, z, diag)
         call mtsolve(n, rr, zz, diag)
         bnrm = dnrm2(n, b)                   !ISDBCG FUNCTION
         err = dnrm2(n, r) / bnrm
         if (err .le. tol) goto 200
!
!         ***** iteration loop *****
!
         do 100 k = 1, itmax
            iter = k
            bknum = 0.d0
            do i = 1, n
               bknum = bknum + z(i) * rr(i)  !DDOT FUNCTION
            enddo
            if (iter .eq. 1) then            !DCOPY SUBROUTINE
               do i = 1, n
                  p(i) = z(i)
                  pp(i) = zz(i)
               enddo
            else
               bk = bknum / bkden
               do i = 1, n
                  p(i) = z(i) + bk * p(i)
                  pp(i) = zz(i) + bk * pp(i)
               enddo
            endif
            bkden = bknum
			
c            call matvec(p, z, iss, bb, nax, nay, naz)

             call matvec_0(p, z, iss, bb, nax, nay, naz, bc_dirichlet)

            akden = 0.d0
            do i = 1, n
               akden = akden + z(i) * pp(i)   !DAXPY SUBROUTINE
            enddo
            ak = bknum / akden
            do i = 1, n
               x(i) = x(i) + ak * p(i)
               r(i) = r(i) - ak * z(i)
            enddo

c            call mttvec(pp, zz, iss, bb, nax, nay, naz)

            call mttvec_0(pp, zz, iss, bb, nax, nay, naz, bc_dirichlet)		

            do i = 1, n
               rr(i) = rr(i) - ak * zz(i)
            enddo
            call msolve(n, r, z, diag)
            call mtsolve(n, rr, zz, diag)
            err = dnrm2(n, r) / bnrm
!CC        WRITE (*,*) ' ITER=',ITER,' ERR=',ERR
            if (err .le. tol) goto 200
100      continue

         print*, 'ERROR - stopping criterion not satisfied'
200      return
      end

      subroutine msolve(n, r, z, diag) !solves a linear system MZ = R with the preconditioning matrix diag
         integer n, i
         real*8 z(n), r(n)
         real*8 diag(*)
         do i = 1, n
            z(i) = r(i) / diag(i)
         enddo
         return
      end

      subroutine mtsolve(n, rr, zz, diag)
         integer n, i
         real*8 zz(n), rr(n)
         real*8 diag(*)
         do i = 1, n
            zz(i) = rr(i) / diag(i)
         enddo
         return
      end

      function dnrm2(n, dx)
         integer n, i
         real*8 dx(n), dnrm2
         dnrm2 = 0.0d0
         do i = 1, n
            dnrm2 = dnrm2 + dx(i)**2
         enddo
         dnrm2 = sqrt(dnrm2)
         return
      end

      subroutine read_cfg

         real*8 q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg, alpha_cfg
         integer*4 i
         character*80 chrtmp
         character*1 ans

         common /cfg/ q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg,
     &   alpha_cfg

         open(1, file = 'vanadis.cfg',status='OLD',form='formatted')
!
         do i = 1, 5
            read(1, 1000) chrtmp
         enddo

         read(1, 2000) q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg,
     &   alpha_cfg

         close(1)
!
         print *, ' '
         print *, ' Configuration parameters:'
         print *, ' '
         write(*, 101) q_cfg
         write(*, 102) k_cfg
         write(*, 103) p_cfg
         write(*, 104) v1_cfg
         write(*, 105) v2_cfg
         write(*, 106) v3_cfg
         write(*, 107) alpha_cfg
         print *,' '
         print *,' Results in T_ttttt.TXT X[m],Y[m],Z[m],
     &		       S(x,y,z)[mg/m^3]'

	     print *,' '
		 print *,' Create calalog OUTPUT in current DIR for results'
		 print *,' '
         print *, ' IS IT CORRECT?  [Y/N]'
         read(*, *) ans
!
         if (.not. (ans .eq. 'Y' .or. ans .eq. 'y')) stop
!
1000     format(a80)
2000     format(d11.2, 1x, d11.2, 5(1x, d10.2))

101      format(' Q[kg/s]    = ', d11.2)
102      format(' K[m^2/s]   = ', d11.2)
103      format(' P[1/s]     =  ', d10.2)
104      format(' V1[m/s]    =  ', d10.2)
105      format(' V2[m/s]    =  ', d10.2)
106      format(' V3[m/s]    =  ', d10.2)
107      format(' ALPHA[m/s] =  ', d10.2)

         return
      end

      subroutine copyright

       print *, '******************************************************'
       print *, ' '
       print *, ' Vanadis 3D - A three-dimensional FEM model '
       print *, '              for atmospheric pollutant transport '
       print *, ' '
       print *, ' author:   Marek Chodorski'
	   print *, ' e-mail:   marek_ac@wp.pl'
       print *, ' '
       print *, '******************************************************'

      return
      end

      subroutine vanadis_files(nx1, nx2, nx3, iczas, istep)

         implicit real*8 (a - h, o - z)
         character*25 dane(iczas)
         character*25 dane_x(nx1, iczas), dane_y(nx2, iczas),
     &   dane_z(nx3, iczas)
         integer x, y, z

         dimension x1((nx1) * (nx2) * (nx3))
     &   , x2((nx1) * (nx2) * (nx3)), x3((nx1) * (nx2) * (nx3))
     &   , s((nx1) * (nx2) * (nx3))

         non = (nx1) * (nx2) * (nx3)

         x = nx1
         y = nx2
         z = nx3

         print*
		 print*,'creating 2D data files...'

         do j = istep, iczas, istep
            do i = 1, x
               write(dane_x(i, j), 105) j, i
            enddo
         enddo
105      format('output\X_', i5.5, '_', i3.3, '.TXT')

         do j = istep, iczas, istep
            do i = 1, y
               write(dane_y(i, j), 106) j, i
            enddo
         enddo
106      format('output\Y_', i5.5, '_', i3.3, '.TXT')

         do j = istep, iczas, istep
            do i = 1, z
               write(dane_z(i, j), 107) j, i
            enddo
         enddo
107      format('output\Z_', i5.5, '_', i3.3, '.TXT')


125      format (e12.5, 2x, e12.5, 2x, e12.5)

         do j = 1, iczas
            write(dane(j), 108) j
         enddo
108      format('output\T_', i5.5, '.TXT')

         do j_iczas = istep, iczas, istep

            open (9, file = dane(j_iczas), status = 'OLD')
            do i = 1, non
               read(9, 222) x1(i), x2(i), x3(i), s(i)
            enddo
            close(9)
222         format (e12.5, 2x, e12.5, 2x, e12.5, 2x, e12.5)

            do j = 1, x
               open (j, file = dane_x(j, j_iczas), status = 'unknown')
               do i = 1, (y * z)
                  write(j, 125) x2(x * (i - 1) + 1 + j - 1),
     &            x3(x * (i - 1) + 1 + j - 1), s(x * (i - 1) + 1 + j -
     &            1)
               enddo
               close(j)
            enddo

            do j = 1, y
               k = 1
               open (j, file = dane_y(j, j_iczas), status = 'unknown')
               do i = 1, (x * z)
                  if (mod(i - 1, x) .eq. 0 .and. (i - 1) .ne. 0) k = k
     &            + 1
                  write(j, 125) x3(i + x * (j - 1) +x * (y - 1) * (k -
     &            1))
     &            , x1(i + x * (j - 1) +x * (y - 1) * (k - 1))
     &            , s(i + x * (j - 1) +x * (y - 1) * (k - 1))
               enddo
               close(j)
            enddo

            do j = 1, z
               open (j, file = dane_z(j, j_iczas), status = 'unknown')
               do i = 1, (x * y)
                  write(j, 125) x2(i + x * y * (j - 1))
     &            , x1(i + x * y * (j - 1)), s(i + x * y * (j - 1))
               enddo
               close(j)
            enddo

         enddo

         return
      end

!--------------------------------------------------------------------
! 4-node quad in 3D
! Input:
!   x(4), y(4), z(4)  - nodal coordinates
!   r, s              - natural coordinates (Gaussian point)
!
! Output:
!   N(4)              - shape functions
!   dNdr(4), dNds(4)  - derivatives wrt r, s
!   J(3,2)            - Jacobian matrix [dx/dr dx/ds; dy/dr dy/ds; dz/dr dz/ds]
!   detJs             - surface Jacobian
!--------------------------------------------------------------------

      subroutine quad4_3d(x, y, z, r, s, detjs)

         real*8 x(4), y(4), z(4)
         real*8 r, s
         real*8 n(4), dndr(4), dnds(4)
         real*8 j(3, 2), detjs
         real*8 xr, xs, yr, ys, zr, zs
         real*8 cx, cy, cz
         integer*4 i

! Shape functions
         n(1) = 0.25d0 * (1.0d0 - r) * (1.0d0 - s)
         n(2) = 0.25d0 * (1.0d0 + r) * (1.0d0 - s)
         n(3) = 0.25d0 * (1.0d0 + r) * (1.0d0 + s)
         n(4) = 0.25d0 * (1.0d0 - r) * (1.0d0 + s)

! Derivatives wrt r
         dndr(1) = -0.25d0 * (1.0d0 - s)
         dndr(2) = 0.25d0 * (1.0d0 - s)
         dndr(3) = 0.25d0 * (1.0d0 + s)
         dndr(4) = -0.25d0 * (1.0d0 + s)

! Derivatives wrt s
         dnds(1) = -0.25d0 * (1.0d0 - r)
         dnds(2) = -0.25d0 * (1.0d0 + r)
         dnds(3) = 0.25d0 * (1.0d0 + r)
         dnds(4) = 0.25d0 * (1.0d0 - r)

! Jacobian components
         xr = 0.0d0
         xs = 0.0d0
         yr = 0.0d0
         ys = 0.0d0
         zr = 0.0d0
         zs = 0.0d0

         do i = 1, 4
            xr = xr + dndr(i) * x(i)
            xs = xs + dnds(i) * x(i)

            yr = yr + dndr(i) * y(i)
            ys = ys + dnds(i) * y(i)

            zr = zr + dndr(i) * z(i)
            zs = zs + dnds(i) * z(i)
         enddo

         j(1, 1) = xr
         j(1, 2) = xs
         j(2, 1) = yr
         j(2, 2) = ys
         j(3, 1) = zr
         j(3, 2) = zs

! Surface Jacobian
         cx = yr * zs - zr * ys
         cy = zr * xs - xr * zs
         cz = xr * ys - yr * xs

         detjs = sqrt(cx * cx + cy * cy + cz * cz)

         return
      end

      subroutine shape_hex8dr_base(xi2, eta2, zeta2, n, dn_dxi, w,
     &                      dw_dxi, lambda, v)
!     8-node brick shape functions and derivatives wrt (xi,eta,zeta)

         real*8 xi2, eta2, zeta2
         real*8 n(8), dn_dxi(3, 8)
         real*8 xm, xp, ym, yp, zm, zp
         real*8 w(8), dw_dxi(3, 8)
         real*8 lambda(3), v(3)

         real*8 xxe(8), yye(8), zze(8)
         real*8 alfa1, alfa2, alfa3
         real*8 gamma1, gamma2, gamma3
         real*8 waga1, waga2, waga3
         real*8 dlug_el_z, dlug_el_y, dlug_el_x

         integer*4 j

         real*8 ksi(8), eta(8), zeta(8)

         common /lok/ ksi, eta, zeta
         common /wsp_lok/ xxe, yye, zze
!$omp threadprivate(/wsp_lok/)

!     Local coordinates
         xm = 1.0d0 - xi2
         xp = 1.0d0 + xi2
         ym = 1.0d0 - eta2
         yp = 1.0d0 + eta2
         zm = 1.0d0 - zeta2
         zp = 1.0d0 + zeta2

!     Shape functions
         n(1) = 0.125d0 * xm * ym * zm
         n(2) = 0.125d0 * xp * ym * zm
         n(3) = 0.125d0 * xp * yp * zm
         n(4) = 0.125d0 * xm * yp * zm
         n(5) = 0.125d0 * xm * ym * zp
         n(6) = 0.125d0 * xp * ym * zp
         n(7) = 0.125d0 * xp * yp * zp
         n(8) = 0.125d0 * xm * yp * zp

!     dN/dxi
         dn_dxi(1, 1) = -0.125d0 * ym * zm
         dn_dxi(1, 2) = 0.125d0 * ym * zm
         dn_dxi(1, 3) = 0.125d0 * yp * zm
         dn_dxi(1, 4) = -0.125d0 * yp * zm
         dn_dxi(1, 5) = -0.125d0 * ym * zp
         dn_dxi(1, 6) = 0.125d0 * ym * zp
         dn_dxi(1, 7) = 0.125d0 * yp * zp
         dn_dxi(1, 8) = -0.125d0 * yp * zp

!     dN/deta
         dn_dxi(2, 1) = -0.125d0 * xm * zm
         dn_dxi(2, 2) = -0.125d0 * xp * zm
         dn_dxi(2, 3) = 0.125d0 * xp * zm
         dn_dxi(2, 4) = 0.125d0 * xm * zm
         dn_dxi(2, 5) = -0.125d0 * xm * zp
         dn_dxi(2, 6) = -0.125d0 * xp * zp
         dn_dxi(2, 7) = 0.125d0 * xp * zp
         dn_dxi(2, 8) = 0.125d0 * xm * zp

!     dN/dzeta
         dn_dxi(3, 1) = -0.125d0 * xm * ym
         dn_dxi(3, 2) = -0.125d0 * xp * ym
         dn_dxi(3, 3) = -0.125d0 * xp * yp
         dn_dxi(3, 4) = -0.125d0 * xm * yp
         dn_dxi(3, 5) = 0.125d0 * xm * ym
         dn_dxi(3, 6) = 0.125d0 * xp * ym
         dn_dxi(3, 7) = 0.125d0 * xp * yp
         dn_dxi(3, 8) = 0.125d0 * xm * yp

!     Element lengths (orthogonal brick)
         dlug_el_z = zze(5) - zze(1)
         dlug_el_y = yye(4) - yye(1)
         dlug_el_x = xxe(2) - xxe(1)

!     Peclet numbers
         gamma3 = dabs(v(3) * dlug_el_z / lambda(3))
         gamma2 = dabs(v(2) * dlug_el_y / lambda(2))
         gamma1 = dabs(v(1) * dlug_el_x / lambda(1))

         alfa3 = 0.0d0
         if (dabs(v(3)) .gt. 0.0d0) then
            alfa3 = (dexp(0.5d0 * gamma3) + dexp(-0.5d0 * gamma3))
     &       / (dexp(0.5d0 * gamma3) - dexp(-0.5d0 * gamma3))
     &      - 2.0d0 / gamma3
         endif

         alfa2 = 0.0d0
         if (dabs(v(2)) .gt. 0.0d0) then
            alfa2 = (dexp(0.5d0 * gamma2) + dexp(-0.5d0 * gamma2))
     &       / (dexp(0.5d0 * gamma2) - dexp(-0.5d0 * gamma2))
     &      - 2.0d0 / gamma2
         endif

         alfa1 = 0.0d0
         if (dabs(v(1)) .gt. 0.0d0) then
            alfa1 = (dexp(0.5d0 * gamma1) + dexp(-0.5d0 * gamma1))
     &       / (dexp(0.5d0 * gamma1) - dexp(-0.5d0 * gamma1))
     &      - 2.0d0 / gamma1
         endif

!     Weights W(J) and their derivatives
         do 10 j = 1, 8
            waga3 = ((1.0d0 + zeta(j) * zeta2) * (1.0d0 - zeta(j) *
     &      zeta2))
     &      * (-0.75d0) * alfa3
            waga2 = ((1.0d0 + eta(j) * eta2) * (1.0d0 - eta(j) * eta2))
     &      * (-0.75d0) * alfa2
            waga1 = ((1.0d0 + ksi(j) * xi2) * (1.0d0 - ksi(j) * xi2))
     &      * (-0.75d0) * alfa1

            w(j) = ((1.0d0 + ksi(j) * xi2) * 0.5d0 + waga1)
     &      * ((1.0d0 + eta(j) * eta2) * 0.5d0 + waga2)
     &      * ((1.0d0 + zeta(j) * zeta2) * 0.5d0 + waga3)

            dw_dxi(1, j) = (ksi(j) * 0.5d0
     &      + 1.5d0 * alfa1 * ksi(j)**2 * xi2)
     &      * ((1.0d0 + zeta(j) * zeta2) * 0.5d0 + waga3)
     &      * ((1.0d0 + eta(j) * eta2) * 0.5d0 + waga2)

            dw_dxi(2, j) = ((1.0d0 + ksi(j) * xi2) * 0.5d0 + waga1)
     &      * ((1.0d0 + zeta(j) * zeta2) * 0.5d0 + waga3)
     &      * (eta(j) * 0.5d0
     &      + 1.5d0 * alfa2 * eta(j)**2 * eta2)

            dw_dxi(3, j) = (zeta(j) * 0.5d0
     &      + 1.5d0 * alfa3 * zeta(j)**2 * zeta2)
     &      * ((1.0d0 + ksi(j) * xi2) * 0.5d0 + waga1)
     &      * ((1.0d0 + eta(j) * eta2) * 0.5d0 + waga2)
 10      continue

      end subroutine shape_hex8dr_base
	  
      subroutine shape_hex8DR(xi2, eta2, zeta2, n, dn_dxi, w,
     &                      dw_dxi, lambda, v)
!  8-node brick shape functions and derivatives wrt (xi,eta,zeta)
!  W tej wersji DR: w i dw_dxi zostaną nadpisane w gauss_point_hex8DR,
!  tutaj liczymy tylko N i dN/dxi, dN/deta, dN/dzeta.

         real*8 xi2, eta2, zeta2
         real*8 n(8), dn_dxi(3, 8)
         real*8 xm, xp, ym, yp, zm, zp
         real*8 w(8), dw_dxi(3, 8)
         real*8 lambda(3), v(3)

         real*8 xxe(8), yye(8), zze(8)
         integer*4 j
         real*8 ksi(8), eta(8), zeta(8)

         common /lok/ ksi, eta, zeta
         common /wsp_lok/ xxe, yye, zze
!$omp threadprivate(/wsp_lok/)

! --- Local coordinates helpers
         xm = 1.0d0 - xi2
         xp = 1.0d0 + xi2
         ym = 1.0d0 - eta2
         yp = 1.0d0 + eta2
         zm = 1.0d0 - zeta2
         zp = 1.0d0 + zeta2

! --- Shape functions N(i)
         n(1) = 0.125d0 * xm * ym * zm
         n(2) = 0.125d0 * xp * ym * zm
         n(3) = 0.125d0 * xp * yp * zm
         n(4) = 0.125d0 * xm * yp * zm
         n(5) = 0.125d0 * xm * ym * zp
         n(6) = 0.125d0 * xp * ym * zp
         n(7) = 0.125d0 * xp * yp * zp
         n(8) = 0.125d0 * xm * yp * zp

! --- dN/dxi
         dn_dxi(1,1) = -0.125d0 * ym * zm
         dn_dxi(1,2) =  0.125d0 * ym * zm
         dn_dxi(1,3) =  0.125d0 * yp * zm
         dn_dxi(1,4) = -0.125d0 * yp * zm
         dn_dxi(1,5) = -0.125d0 * ym * zp
         dn_dxi(1,6) =  0.125d0 * ym * zp
         dn_dxi(1,7) =  0.125d0 * yp * zp
         dn_dxi(1,8) = -0.125d0 * yp * zp

! --- dN/deta
         dn_dxi(2,1) = -0.125d0 * xm * zm
         dn_dxi(2,2) = -0.125d0 * xp * zm
         dn_dxi(2,3) =  0.125d0 * xp * zm
         dn_dxi(2,4) =  0.125d0 * xm * zm
         dn_dxi(2,5) = -0.125d0 * xm * zp
         dn_dxi(2,6) = -0.125d0 * xp * zp
         dn_dxi(2,7) =  0.125d0 * xp * zp
         dn_dxi(2,8) =  0.125d0 * xm * zp

! --- dN/dzeta
         dn_dxi(3,1) = -0.125d0 * xm * ym
         dn_dxi(3,2) = -0.125d0 * xp * ym
         dn_dxi(3,3) = -0.125d0 * xp * yp
         dn_dxi(3,4) = -0.125d0 * xm * yp
         dn_dxi(3,5) =  0.125d0 * xm * ym
         dn_dxi(3,6) =  0.125d0 * xp * ym
         dn_dxi(3,7) =  0.125d0 * xp * yp
         dn_dxi(3,8) =  0.125d0 * xm * yp

! W i dW/dxi zostaną ustawione w gauss_point_hex8DR (DR)
         do j = 1, 8
            w(j) = n(j)
            dw_dxi(1,j) = dn_dxi(1,j)
            dw_dxi(2,j) = dn_dxi(2,j)
            dw_dxi(3,j) = dn_dxi(3,j)
         enddo

      end subroutine shape_hex8DR
	  
	  

      subroutine gauss_point_hex8DR_base(xi, eta, zeta, xnod, ynod, znod,
     &detj, lambda, v, n, w, dnx1, dny1, dnz1, dnxw, dnyw, dnzw)

! Evaluate N, dN/dx, dN/dy, dN/dz and det(J) at Gauss point

         real*8 xi, eta, zeta
         real*8 xnod(8), ynod(8), znod(8)
         real*8 n(8), detj
         real*8 dn_dxi(3, 8)
         real*8 j(3, 3), invj(3, 3)
         real*8 tmp(3)
         integer*4 i, jj, k
		 real*8 detji

         real*8 lambda(3), v(3)
         real*8 w(8), dw_dxi(3, 8), dwdx(3, 8)

         real*8 dnx1(8), dny1(8), dnz1(8)
         real*8 dnxw(8), dnyw(8), dnzw(8)

         call shape_hex8dr_base(xi, eta, zeta, n, dn_dxi, w, dw_dxi,
     &   lambda, v)

! Build Jacobian J = d(x,y,z)/d(xi,eta,zeta)

         do i=1,3
           do jj=1,3
             j(i,jj) = 0.0d0
		   enddo
         enddo

         do i = 1, 8
            j(1, 1) = j(1, 1) + xnod(i) * dn_dxi(1, i)
            j(1, 2) = j(1, 2) + xnod(i) * dn_dxi(2, i)
            j(1, 3) = j(1, 3) + xnod(i) * dn_dxi(3, i)

            j(2, 1) = j(2, 1) + ynod(i) * dn_dxi(1, i)
            j(2, 2) = j(2, 2) + ynod(i) * dn_dxi(2, i)
            j(2, 3) = j(2, 3) + ynod(i) * dn_dxi(3, i)

            j(3, 1) = j(3, 1) + znod(i) * dn_dxi(1, i)
            j(3, 2) = j(3, 2) + znod(i) * dn_dxi(2, i)
            j(3, 3) = j(3, 3) + znod(i) * dn_dxi(3, i)
         end do

! Determinant of J
         detj = j(1, 1) * (j(2, 2) * j(3, 3) - j(2, 3) * j(3, 2))
     &        - j(1, 2) * (j(2, 1) * j(3, 3) - j(2, 3) * j(3, 1))
     &        + j(1, 3) * (j(2, 1) * j(3, 2) - j(2, 2) * j(3, 1))

! Inverse of J (explicit 3x3 inverse)
         detji=1.0d0/detj

         invj(1, 1) = detji * (j(2, 2) * j(3, 3) - j(2, 3) * j(3, 2))
         invj(1, 2) = detji * (j(1, 3) * j(3, 2) - j(1, 2) * j(3, 3))
         invj(1, 3) = detji * (j(1, 2) * j(2, 3) - j(1, 3) * j(2, 2))

         invj(2, 1) = detji * (j(2, 3) * j(3, 1) - j(2, 1) * j(3, 3))
         invj(2, 2) = detji * (j(1, 1) * j(3, 3) - j(1, 3) * j(3, 1))
         invj(2, 3) = detji * (j(1, 3) * j(2, 1) - j(1, 1) * j(2, 3))

         invj(3, 1) = detji * (j(2, 1) * j(3, 2) - j(2, 2) * j(3, 1))
         invj(3, 2) = detji * (j(1, 2) * j(3, 1) - j(1, 1) * j(3, 2))
         invj(3, 3) = detji * (j(1, 1) * j(2, 2) - j(1, 2) * j(2, 1))

! dN/dx = invJ * dN/dxi

         do i = 1, 8
            dnx1(i) = invj(1, 1) * dn_dxi(1, i) + invj(2, 1) * 
     &                dn_dxi(2, i) + invj(3, 1) * dn_dxi(3, i)
            dny1(i) = invj(1, 2) * dn_dxi(1, i) + invj(2, 2) *
     &                dn_dxi(2, i) + invj(3, 2) * dn_dxi(3, i)
            dnz1(i) = invj(1, 3) * dn_dxi(1, i) + invj(2, 3) *
     &                dn_dxi(2, i) + invj(3, 3) * dn_dxi(3, i)
         enddo

! dW/dx = invJ * dW/dxi

         do i = 1, 8
            dnxw(i) = invj(1, 1) * dw_dxi(1, i) + invj(2, 1) * 
     &                dw_dxi(2, i) + invj(3, 1) * dw_dxi(3, i)
            dnyw(i) = invj(1, 2) * dw_dxi(1, i) + invj(2, 2) *
     &                dw_dxi(2, i) + invj(3, 2) * dw_dxi(3, i)
            dnzw(i) = invj(1, 3) * dw_dxi(1, i) + invj(2, 3) *
     &                dw_dxi(2, i) + invj(3, 3) * dw_dxi(3, i)
         enddo

      end subroutine gauss_point_hex8DR_base	  
c
      subroutine gauss_point_hex8DR(xi, eta, zeta, xnod, ynod, znod,
     & detj, lambda, v, n, w, dnx1, dny1, dnz1, dnxw, dnyw, dnzw)

! Evaluate N, dN/dx, dN/dy, dN/dz, W_DR, and det(J) at Gauss point
! W_DR = N + tau_DR * (v · grad N)   (Directional Residual)
! tau_DR = tau(Pe_stream), Pe wzdłuż kierunku prędkości

         implicit real*8 (a-h,o-z)

         real*8 xi, eta, zeta
         real*8 xnod(8), ynod(8), znod(8)
         real*8 n(8), detj
         real*8 dn_dxi(3, 8)
         real*8 j(3, 3), invj(3, 3)
         real*8 lambda(3), v(3)
         real*8 w(8), dw_dxi(3, 8)

         real*8 dnx1(8), dny1(8), dnz1(8)
         real*8 dnxw(8), dnyw(8), dnzw(8)
		 
		 real*8 dlx, dly, dlz
		 real*8 vmag, eps, h_stream, lambda_s, Pe, tau_DR
         real*8 ex, ey, ez, vdotg
         real*8 epos, eneg, coshPe, sinhPe, cothPe

         real*8 detji
         integer*4 i, jj

! --- Najpierw klasyczne N i dN/dxi (shape_hex8DR)
         call shape_hex8DR(xi, eta, zeta, n, dn_dxi, w, dw_dxi,
     &                     lambda, v)

! --- Jacobian J = d(x,y,z)/d(xi,eta,zeta)
         do i = 1, 3
            do jj = 1, 3
               j(i,jj) = 0.0d0
            enddo
         enddo

         do i = 1, 8
            j(1,1) = j(1,1) + xnod(i) * dn_dxi(1,i)
            j(1,2) = j(1,2) + xnod(i) * dn_dxi(2,i)
            j(1,3) = j(1,3) + xnod(i) * dn_dxi(3,i)

            j(2,1) = j(2,1) + ynod(i) * dn_dxi(1,i)
            j(2,2) = j(2,2) + ynod(i) * dn_dxi(2,i)
            j(2,3) = j(2,3) + ynod(i) * dn_dxi(3,i)

            j(3,1) = j(3,1) + znod(i) * dn_dxi(1,i)
            j(3,2) = j(3,2) + znod(i) * dn_dxi(2,i)
            j(3,3) = j(3,3) + znod(i) * dn_dxi(3,i)
         enddo

! --- det(J)
         detj = j(1,1) * (j(2,2)*j(3,3) - j(2,3)*j(3,2))
     &        - j(1,2) * (j(2,1)*j(3,3) - j(2,3)*j(3,1))
     &        + j(1,3) * (j(2,1)*j(3,2) - j(2,2)*j(3,1))

         detji = 1.0d0 / detj

! --- J^{-1}
         invj(1,1) = detji * (j(2,2)*j(3,3) - j(2,3)*j(3,2))
         invj(1,2) = detji * (j(1,3)*j(3,2) - j(1,2)*j(3,3))
         invj(1,3) = detji * (j(1,2)*j(2,3) - j(1,3)*j(2,2))

         invj(2,1) = detji * (j(2,3)*j(3,1) - j(2,1)*j(3,3))
         invj(2,2) = detji * (j(1,1)*j(3,3) - j(1,3)*j(3,1))
         invj(2,3) = detji * (j(1,3)*j(2,1) - j(1,1)*j(2,3))

         invj(3,1) = detji * (j(2,1)*j(3,2) - j(2,2)*j(3,1))
         invj(3,2) = detji * (j(1,2)*j(3,1) - j(1,1)*j(3,2))
         invj(3,3) = detji * (j(1,1)*j(2,2) - j(1,2)*j(2,1))

! --- dN/dx = invJ * dN/dxi
         do i = 1, 8
            dnx1(i) = invj(1,1)*dn_dxi(1,i) + invj(2,1)*dn_dxi(2,i)
     &               + invj(3,1)*dn_dxi(3,i)
            dny1(i) = invj(1,2)*dn_dxi(1,i) + invj(2,2)*dn_dxi(2,i)
     &               + invj(3,2)*dn_dxi(3,i)
            dnz1(i) = invj(1,3)*dn_dxi(1,i) + invj(2,3)*dn_dxi(2,i)
     &               + invj(3,3)*dn_dxi(3,i)
         enddo

! --- Długości elementu w osiach (zakładamy prostopadłościan)

         dlx = dabs(xnod(2) - xnod(1))
         dly = dabs(ynod(4) - ynod(1))
         dlz = dabs(znod(5) - znod(1))

! --- Prędkość i Peclet wzdłuż kierunku v
         eps   = 1.0d-14
         vmag  = dsqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))

! efektywna długość wzdłuż strugi (klasyczna definicja 1/h = sum |v_i|/h_i / |v|)
         if (vmag .gt. eps) then
            h_stream = 2.0d0 / ( dabs(v(1))/dlx + dabs(v(2))/dly
     &                          + dabs(v(3))/dlz + eps )
         else
            h_stream = 0.0d0
         endif

! dyfuzja izotropowa
         lambda_s = lambda(1)

         if (vmag .gt. eps .and. lambda_s .gt. eps) then
            Pe = vmag * h_stream / (2.0d0 * lambda_s)
         else
            Pe = 0.0d0
         endif

! --- tau_DR(Pe) – klasyczna funkcja SUPG/DR:
!     tau = h/(2|v|) * (coth(Pe) - 1/Pe)
         tau_DR = 0.0d0
         if (Pe .gt. 1.0d-12 .and. vmag .gt. eps) then
            ! cosh/sinh z exp
            epos   = dexp(Pe)
            eneg   = dexp(-Pe)
            coshPe = 0.5d0 * (epos + eneg)
            sinhPe = 0.5d0 * (epos - eneg)
            cothPe = coshPe / sinhPe
            tau_DR = 0.5d0 * h_stream / vmag * (cothPe - 1.0d0/Pe)
         endif

! --- Kierunek jednostkowy prędkości
         if (vmag .gt. eps) then
            ex = v(1) / vmag
            ey = v(2) / vmag
            ez = v(3) / vmag
         else
            ex = 0.0d0
            ey = 0.0d0
            ez = 0.0d0
         endif

! --- Budowa funkcji testowej DR:
!     w_j = N_j + tau_DR * (v · grad N_j)
!     v · grad N_j = vmag * (e · grad N_j)
         do i = 1, 8
            vdotg = v(1)*dnx1(i) + v(2)*dny1(i) + v(3)*dnz1(i)
            w(i)  = n(i) + tau_DR * vdotg
         enddo

! --- Dla dyfuzji przyjmujemy grad W ≈ grad N (czysto DR w konwekcji)
         do i = 1, 8
            dnxw(i) = dnx1(i)
            dnyw(i) = dny1(i)
            dnzw(i) = dnz1(i)
         enddo

      end subroutine gauss_point_hex8DR
c
c
      subroutine build_bc_nodes(iss, kod, bc_dirichlet, nax, nay, naz)
      implicit real*8 (a-h,o-z)

      integer*4 bc_dirichlet(*)
      integer*4 iss((nax-1)*(nay-1)*(naz-1),8)
      integer*4 kod((nax-1)*(nay-1)*(naz-1),6)

      common /cal_el/ ilosc_el, ilosc_wezlow, i_ot_elem(8)
!$omp threadprivate(/cal_el/)

c--- wyzeruj maskę
      do i = 1, ilosc_wezlow
         bc_dirichlet(i) = 0
      enddo

c--- przejście po wszystkich elementach
      do kel = 1, ilosc_el

c--- pobierz węzły elementu
         do j = 1, 8
            i_ot_elem(j) = iss(kel, j)
         enddo

c--- Dirichlet na ścianach 1,2,4,5,6 (wszystkie oprócz 3)
         do iwb = 1, 6
            if (iwb .ne. 3) then
               if (kod(kel,iwb) .eq. iwb) then
                  call mark_face_nodes(iwb, i_ot_elem, bc_dirichlet)
               endif
            endif
         enddo

      enddo

      return
      end
c
      subroutine mark_face_nodes(iwb, nodes, bc_dirichlet)
      implicit real*8 (a-h,o-z)

      integer*4 nodes(8), bc_dirichlet(*)

c--- Ściana 1: eta = -1
      if (iwb .eq. 1) then
         bc_dirichlet(nodes(1)) = 1
         bc_dirichlet(nodes(2)) = 1
         bc_dirichlet(nodes(3)) = 1
         bc_dirichlet(nodes(4)) = 1
      endif

c--- Ściana 2: eta = +1
      if (iwb .eq. 2) then
         bc_dirichlet(nodes(5)) = 1
         bc_dirichlet(nodes(6)) = 1
         bc_dirichlet(nodes(7)) = 1
         bc_dirichlet(nodes(8)) = 1
      endif

c--- Ściana 3: ksi = -1 (WYŁĄCZONA, ale zostawiamy dla kompletności)
      if (iwb .eq. 3) then
         bc_dirichlet(nodes(1)) = 1
         bc_dirichlet(nodes(4)) = 1
         bc_dirichlet(nodes(8)) = 1
         bc_dirichlet(nodes(5)) = 1
      endif

c--- Ściana 4: ksi = +1
      if (iwb .eq. 4) then
         bc_dirichlet(nodes(2)) = 1
         bc_dirichlet(nodes(3)) = 1
         bc_dirichlet(nodes(7)) = 1
         bc_dirichlet(nodes(6)) = 1
      endif

c--- Ściana 5: zeta = -1
      if (iwb .eq. 5) then
         bc_dirichlet(nodes(1)) = 1
         bc_dirichlet(nodes(2)) = 1
         bc_dirichlet(nodes(6)) = 1
         bc_dirichlet(nodes(5)) = 1
      endif

c--- Ściana 6: zeta = +1
      if (iwb .eq. 6) then
         bc_dirichlet(nodes(4)) = 1
         bc_dirichlet(nodes(3)) = 1
         bc_dirichlet(nodes(7)) = 1
         bc_dirichlet(nodes(8)) = 1
      endif

      return
      end
c
c======================================================================
c  Zaznaczanie węzłów ściany na podstawie numeru ściany iwb
c  Konwencja hexa8 (lokalne numery węzłów):
c
c      z
c      ^
c      |
c      4--------3
c     /|       /|
c    8--------7 |
c    | |      | |
c    | 1------|-2 --> y
c    |/       |/
c    5--------6
c
c  Ściany:
c    iwb = 1 : eta = -1  -> (1,2,3,4)
c    iwb = 2 : eta = +1  -> (5,6,7,8)
c    iwb = 3 : ksi = -1  -> (1,4,8,5)
c    iwb = 4 : ksi = +1  -> (2,3,7,6)
c    iwb = 5 : zeta = -1 -> (1,2,6,5)
c    iwb = 6 : zeta = +1 -> (4,3,7,8)
c======================================================================

c
      subroutine matvec_0(xin, xout, iss, bb, nax, nay, naz,
     &                    bc_dirichlet)

      implicit real*8 (a - h, o - z)
      dimension xin(*), xout(*)
      dimension iss((nax - 1) * (nay - 1) * (naz - 1), 8), bb(*)
	  integer*4 bc_dirichlet
      dimension bc_dirichlet(*)

      common /cal_el/ ilosc_el, ilosc_wezlow, i_ot_elem(8)
!$omp threadprivate(/cal_el/)

c----- wyzeruj wynik
      do j = 1, ilosc_wezlow
         xout(j) = 0.0d0
      enddo

c----- wymuś Dirichlet: traktuj xin=0 na węzłach brzegowych
      do j = 1, ilosc_wezlow
         if (bc_dirichlet(j).eq.1) xin(j)=0.0d0
      enddo

      ik = 1
      do je = 1, ilosc_el
         do j = 1, 8
            i_ot_elem(j) = iss(je, j)
         enddo

         do i = 1, 8
            ii = i_ot_elem(i)
            do j = 1, 8
               jj = i_ot_elem(j)
               xout(ii) = xout(ii) + bb(ik) * xin(jj)
               ik = ik + 1
            enddo
         enddo
      enddo

c----- wymuś Dirichlet na wyniku
      do j = 1, ilosc_wezlow
         if (bc_dirichlet(j).eq.1) xout(j)=0.0d0
      enddo

      return
      end
c
      subroutine mttvec_0(xin, xout, iss, bb, nax, nay, naz,
     &                    bc_dirichlet)

      implicit real*8 (a - h, o - z)
      dimension xin(*), xout(*), b3(8, 8)
      dimension iss((nax - 1) * (nay - 1) * (naz - 1), 8), bb(*)
	  integer*4 bc_dirichlet
      dimension bc_dirichlet(*)

      common /cal_el/ ilosc_el, ilosc_wezlow, i_ot_elem(8)
!$omp threadprivate(/cal_el/)

c----- wyzeruj wynik
      do j = 1, ilosc_wezlow
         xout(j) = 0.0d0
      enddo

c----- wymuś Dirichlet: xin=0 na węzłach brzegowych
      do j = 1, ilosc_wezlow
         if (bc_dirichlet(j).eq.1) xin(j)=0.0d0
      enddo

      do je = 1, ilosc_el
         do j = 1, 8
            i_ot_elem(j) = iss(je, j)
         enddo

         ik3 = 0
         do i = 1, 8
            ik3 = ik3 + 1
            do j = 1, 8
               b3(i, j) = bb(j + 8 * (ik3 - 1) + 64 * (je - 1))
            enddo
         enddo

         do i = 1, 8
            ii = i_ot_elem(i)
            do j = 1, 8
               jj = i_ot_elem(j)
               xout(ii) = xout(ii) + b3(j, i) * xin(jj)
            enddo
         enddo
      enddo

c----- wymuś Dirichlet na wyniku
      do j = 1, ilosc_wezlow
         if (bc_dirichlet(j).eq.1) xout(j)=0.0d0
      enddo

      return
      end
	  
	  
      logical function is_diagonally_dominant(ilosc_el, ilosc_wezlow,
     &                                        iss, bb)
          implicit none
          integer ilosc_el, ilosc_wezlow
          integer iss(ilosc_el, 8)
          double precision bb(64*ilosc_el)
      
          double precision, allocatable :: diag(:), sumrow(:)
          integer je, i, j, gi, gj, ik3
          double precision b3(8,8)
      
          allocate(diag(ilosc_wezlow))
          allocate(sumrow(ilosc_wezlow))
      
          diag   = 0.0d0
          sumrow = 0.0d0
      
          ! --- przejście po elementach ---
          do je = 1, ilosc_el
      
              ! odczyt bloku 8x8
              ik3 = 0
              do i = 1, 8
                  ik3 = ik3 + 1
                  do j = 1, 8
                      b3(i,j) = bb(j + 8*(ik3-1) + 64*(je-1))
                  enddo
              enddo
      
              ! wkład do sum diagonalnych i pozostałych
              do i = 1, 8
                  gi = iss(je,i)
                  do j = 1, 8
                      gj = iss(je,j)
      
                      if (gi == gj) then
                          diag(gi) = diag(gi) + abs(b3(i,j))
                      else
                          sumrow(gi) = sumrow(gi) + abs(b3(i,j))
                      endif
      
                  enddo
              enddo
      
          enddo
      
          ! --- sprawdzenie warunku dominacji ---
          do i = 1, ilosc_wezlow
              if (diag(i) <= sumrow(i)) then
                  is_diagonally_dominant = .false.
                  return
              endif
          enddo
      
          is_diagonally_dominant = .true.
      end function is_diagonally_dominant	   
	  
	  
      subroutine ascii_heatmap_ppm(ilosc_el, ilosc_wezlow, iss, bb)
      implicit none
      integer ilosc_el, ilosc_wezlow
      integer iss(ilosc_el, 8)
      double precision bb(64*ilosc_el)

      integer NX, NY
      parameter (NX=300, NY=300)

      double precision H(NY, NX)
      double precision b3(8,8)
      double precision sx, sy
      double precision val, maxH
      integer je, i, j, gi, gj
      integer ix, iy
      integer ii, jj
      integer dx, dy

      integer R, G, B
      integer unit

C --- zerowanie heatmapy ---
      do 10 ii = 1, NY
         do 20 jj = 1, NX
            H(ii,jj) = 0.0d0
20       continue
10    continue

C --- skala indeksów ---
      sx = dble(ilosc_wezlow) / dble(NX)
      sy = dble(ilosc_wezlow) / dble(NY)

C --- przejście po elementach ---
      do 100 je = 1, ilosc_el

C --- odczyt bloku 8x8 ---
         do 30 i = 1, 8
            do 40 j = 1, 8
               b3(i,j) = dabs(bb(j + 8*(i-1) + 64*(je-1)))
40          continue
30       continue

C --- oversampling 4x4 ---
         do 50 i = 1, 8
            gi = iss(je,i)
            iy = int((gi-1)/sy) + 1
            if (iy .gt. NY) iy = NY

            do 60 j = 1, 8
               gj = iss(je,j)
               ix = int((gj-1)/sx) + 1
               if (ix .gt. NX) ix = NX

               val = b3(i,j) * 0.0625d0

               do 70 dy = 0, 3
                  if (iy+dy .gt. NY) goto 70
                  do 80 dx = 0, 3
                     if (ix+dx .gt. NX) goto 80
                     H(iy+dy, ix+dx) = H(iy+dy, ix+dx) + val
80                continue
70             continue

60          continue
50       continue

100   continue

C --- normalizacja ---
      maxH = 0.0d0
      do 110 ii = 1, NY
         do 120 jj = 1, NX
            if (H(ii,jj) .gt. maxH) maxH = H(ii,jj)
120      continue
110   continue
      if (maxH .eq. 0.0d0) maxH = 1.0d0

C --- otwarcie pliku ---
      unit = 99
      open(unit, file='heatmap.ppm', status='unknown')

C --- nagłówek PPM ---
      write(unit,*) 'P3'
      write(unit,*) NX, NY
      write(unit,*) 255

C --- zapis pikseli ---
      do 200 iy = 1, NY
         do 210 ix = 1, NX
            val = H(iy,ix) / maxH

C --- paleta A: niebieski → zielony → czerwony ---
C --- R, G, B ∈ [0,255]

            if (val .le. 0.25d0) then
C --- niebieski → cyjan
               R = 0
               G = int(1020.0d0 * val)
               B = 255

            else if (val .le. 0.50d0) then
C --- cyjan → zielony
               R = 0
               G = 255
               B = int(255.0d0 - 1020.0d0*(val-0.25d0))

            else if (val .le. 0.75d0) then
C --- zielony → żółty
               R = int(1020.0d0*(val-0.50d0))
               G = 255
               B = 0

            else
C --- żółty → czerwony
               R = 255
               G = int(255.0d0 - 1020.0d0*(val-0.75d0))
               B = 0
            endif

            write(unit,*) R, G, B

210      continue
200   continue

      close(unit)
      return
      end
	  
           double precision function diagonal_dominance_factor(
     &                         ilosc_el, ilosc_wezlow, iss, bb)
           implicit none
           integer ilosc_el, ilosc_wezlow
           integer iss(ilosc_el, 8)
           double precision bb(64*ilosc_el)
       
           double precision, allocatable :: diag(:), sumrow(:)
           integer je, i, j, gi, gj
           integer idx
           double precision aij, ti, M
       
           allocate(diag(ilosc_wezlow))
           allocate(sumrow(ilosc_wezlow))
       
           diag   = 0.0d0
           sumrow = 0.0d0
       
           ! --- przejście po elementach ---
           do je = 1, ilosc_el
       
               do i = 1, 8
                   gi = iss(je, i)
       
                   do j = 1, 8
                       gj = iss(je, j)
       
                       ! indeks w bb
                       idx = (j + 8*(i-1)) + 64*(je-1)
                       aij = bb(idx)
       
                       if (gi == gj) then
                           diag(gi) = diag(gi) + abs(aij)
                       else
                           sumrow(gi) = sumrow(gi) + abs(aij)
                       endif
       
                   enddo
               enddo
       
           enddo
       
           ! --- obliczenie M(A) ---
           M = 0.0d0
           do i = 1, ilosc_wezlow
               if (diag(i) > 0.0d0) then
                   ti = sumrow(i) / diag(i)
                   if (ti > M) M = ti
               endif
           enddo
       
           diagonal_dominance_factor = M
       
           deallocate(diag, sumrow)
       
       end function diagonal_dominance_factor
	   
	   
	   
	  subroutine integrate_mass_hex8(nelem, conn, xx, yy, zz, c, mass)
c-----------------------------------------------------------------
c  Oblicza masę M = ∫ c(x) dΩ dla elementów HEX8 (2×2×2 Gauss)
c-----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer nelem, conn(nelem,8)
      real*8 xx(*), yy(*), zz(*), c(*)
      real*8 mass

      real*8 xi_g(2), w_g(2)
      real*8 xi, eta, zeta, wgt
      real*8 n(8), dN_dxi(3,8)
      real*8 J(3,3), detJ
      real*8 xc(8), yc(8), zc(8), cc(8)
      real*8 c_val
      integer e, a, ig, jg, kg

c--- Gauss 2×2×2
      xi_g(1) = -1.0d0/dsqrt(3.0d0)
      xi_g(2) =  1.0d0/dsqrt(3.0d0)
      w_g(1)  =  1.0d0
      w_g(2)  =  1.0d0

      mass = 0.0d0

c--- pętla po elementach
      do 100 e = 1, nelem

c-------- pobierz współrzędne i stężenia węzłów elementu
         do 10 a = 1, 8
            xc(a) = xx(conn(e,a))
            yc(a) = yy(conn(e,a))
            zc(a) = zz(conn(e,a))
            cc(a) = c(conn(e,a))
10       continue

c-------- Gauss 2×2×2
         do 200 ig = 1, 2
            xi = xi_g(ig)
            do 210 jg = 1, 2
               eta = xi_g(jg)
               do 220 kg = 1, 2
                  zeta = xi_g(kg)
                  wgt  = w_g(ig)*w_g(jg)*w_g(kg)

c-------------- funkcje kształtu
                  call shape_hex8_only_n(xi, eta, zeta, n, dN_dxi)

c-------------- Jacobian
                  call jacobian_hex8(dN_dxi, xc, yc, zc, J, detJ)

c-------------- interpolacja stężenia
                  c_val = 0.0d0
                  do 30 a = 1, 8
                     c_val = c_val + n(a)*cc(a)
30                continue

c-------------- wkład do masy
                  mass = mass + c_val * detJ * wgt

220            continue
210         continue
200      continue

100   continue

      return
      end


      subroutine jacobian_hex8(dN_dxi, x, y, z, J, detJ)
      implicit real*8 (a-h,o-z)
      real*8 dN_dxi(3,8), x(8), y(8), z(8)
      real*8 J(3,3), detJ
      integer a

c--- zerowanie
      do 5 a = 1, 3
         J(a,1) = 0.0d0
         J(a,2) = 0.0d0
         J(a,3) = 0.0d0
5     continue

c--- składanie Jacobianu
      do 10 a = 1, 8
         J(1,1) = J(1,1) + dN_dxi(1,a)*x(a)
         J(1,2) = J(1,2) + dN_dxi(2,a)*x(a)
         J(1,3) = J(1,3) + dN_dxi(3,a)*x(a)

         J(2,1) = J(2,1) + dN_dxi(1,a)*y(a)
         J(2,2) = J(2,2) + dN_dxi(2,a)*y(a)
         J(2,3) = J(2,3) + dN_dxi(3,a)*y(a)

         J(3,1) = J(3,1) + dN_dxi(1,a)*z(a)
         J(3,2) = J(3,2) + dN_dxi(2,a)*z(a)
         J(3,3) = J(3,3) + dN_dxi(3,a)*z(a)
10    continue

c--- wyznacznik
      detJ = J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))
     &     - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1))
     &     + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))

      return
      end
	  
	  
c only N
c
      subroutine shape_hex8_only_n(xi2, eta2, zeta2, n, dn_dxi)
!     8-node brick shape functions and derivatives wrt (xi,eta,zeta)

         real*8 xi2, eta2, zeta2
         real*8 n(8), dn_dxi(3, 8)
         real*8 xm, xp, ym, yp, zm, zp
         real*8 w(8), dw_dxi(3, 8)
         real*8 lambda(3), v(3)

         real*8 xxe(8), yye(8), zze(8)
         real*8 alfa1, alfa2, alfa3
         real*8 gamma1, gamma2, gamma3
         real*8 waga1, waga2, waga3
         real*8 dlug_el_z, dlug_el_y, dlug_el_x

         integer*4 j

         real*8 ksi(8), eta(8), zeta(8)

         common /lok/ ksi, eta, zeta
         common /wsp_lok/ xxe, yye, zze
!$omp threadprivate(/wsp_lok/)

!     Local coordinates
         xm = 1.0d0 - xi2
         xp = 1.0d0 + xi2
         ym = 1.0d0 - eta2
         yp = 1.0d0 + eta2
         zm = 1.0d0 - zeta2
         zp = 1.0d0 + zeta2

!     Shape functions
         n(1) = 0.125d0 * xm * ym * zm
         n(2) = 0.125d0 * xp * ym * zm
         n(3) = 0.125d0 * xp * yp * zm
         n(4) = 0.125d0 * xm * yp * zm
         n(5) = 0.125d0 * xm * ym * zp
         n(6) = 0.125d0 * xp * ym * zp
         n(7) = 0.125d0 * xp * yp * zp
         n(8) = 0.125d0 * xm * yp * zp

!     dN/dxi
         dn_dxi(1, 1) = -0.125d0 * ym * zm
         dn_dxi(1, 2) = 0.125d0 * ym * zm
         dn_dxi(1, 3) = 0.125d0 * yp * zm
         dn_dxi(1, 4) = -0.125d0 * yp * zm
         dn_dxi(1, 5) = -0.125d0 * ym * zp
         dn_dxi(1, 6) = 0.125d0 * ym * zp
         dn_dxi(1, 7) = 0.125d0 * yp * zp
         dn_dxi(1, 8) = -0.125d0 * yp * zp

!     dN/deta
         dn_dxi(2, 1) = -0.125d0 * xm * zm
         dn_dxi(2, 2) = -0.125d0 * xp * zm
         dn_dxi(2, 3) = 0.125d0 * xp * zm
         dn_dxi(2, 4) = 0.125d0 * xm * zm
         dn_dxi(2, 5) = -0.125d0 * xm * zp
         dn_dxi(2, 6) = -0.125d0 * xp * zp
         dn_dxi(2, 7) = 0.125d0 * xp * zp
         dn_dxi(2, 8) = 0.125d0 * xm * zp

!     dN/dzeta
         dn_dxi(3, 1) = -0.125d0 * xm * ym
         dn_dxi(3, 2) = -0.125d0 * xp * ym
         dn_dxi(3, 3) = -0.125d0 * xp * yp
         dn_dxi(3, 4) = -0.125d0 * xm * yp
         dn_dxi(3, 5) = 0.125d0 * xm * ym
         dn_dxi(3, 6) = 0.125d0 * xp * ym
         dn_dxi(3, 7) = 0.125d0 * xp * yp
         dn_dxi(3, 8) = 0.125d0 * xm * yp

      end subroutine
c	  