! Vanadis Solver / DR Stabilization
! Copyright (c) 1995-2026 Marek Chodorski
! Dual licensing:
!   - GPL-3.0 (open-source)
!   - Commercial license (contact author)
! See LICENSE.md for details.
!
! VANADIS 3D v3.1.0 development line - 12.09.2026
! Terrain-following mesh with terrain-consistent velocity field
! Quadrature: 3x3x3 HEX8 volume, 3x3 Q4 boundary
! DR: Jacobian metric h_stream, safe Pe, pure-advection K=0 limit
! Fields: v(x,y,z,t), K(x,y,z,t), alpha(x,y,z,t) evaluated locally
!
! autor : Marek Chodorski
! e-mail: marek_ac@wp.pl
!
! Vanadis 3D - A three-dimensional dynamic Eulerian finite-element model
! for atmospheric pollutant transport
!
      program vanadis

         parameter (nax = 7, nay = 20, naz = 20)
         parameter (max_time = 2000)
         parameter (i_time_step = 100)

         implicit real*8 (a - h, o - z)

c--- Full mass-balance QA diagnostic.  Set to 0 for production runs.
         integer*4 imass_balance
         parameter (imass_balance = 1)
         dimension dx(nax), dy(nay), dz(naz)
         dimension xx(nax * nay * naz), yy(nax * nay * naz),
     &             zz(nax * nay * naz)
         dimension t(6)
         dimension iss((nax - 1) * (nay - 1) * (naz - 1), 8)
         dimension kod((nax - 1) * (nay - 1) * (naz - 1), 6)
         dimension bb((nax - 1) * (nay - 1) * (naz - 1) * 64)
         dimension tnw2(nax * nay * naz)
         dimension tnw_prev(nax * nay * naz)

         real*8 mass, mass_old, dm
         real*8 rq_n, rq_np1, rin_n, rin_np1, rout_n, rout_np1
         real*8 rp_n, rp_np1, ra_n, ra_np1
         real*8 rdiv_n, rdiv_np1
         real*8 mb_q, mb_in, mb_out, mb_div, mb_p, mb_alpha
         real*8 mb_dir, eps_num
         real*8 cum_q, cum_in, cum_out, cum_div, cum_p
         real*8 cum_alpha, cum_dir
         real*8 cum_eps, dtb

         real*8 tsplit1, tsplit2, tall1, tall2
         real*8 omp_get_wtime

         common /d_time/ i_d_czas

         il_iteracji = 1000
         i_d_czas = i_time_step
         ttol = 1d-8

         call copyright
         call read_cfg
         call make_output_dir
         call par_3d
         call wsp_na_o (nax, nay, naz, dx, dy, dz)
         call terrain_init(dy, dz, nay, naz)
         call siatka (xx, yy, zz, dx, dy, dz, nax, nay, naz)
         call w_el (iss, nax, nay, naz)
         call kod_el_dr(nax, nay, naz, kod)

         i_czas = 0
         nelem = (nax-1)*(nay-1)*(naz-1)

         do i = 1, nax*nay*naz
            tnw_prev(i) = 0.0d0
         enddo

         if (imass_balance .eq. 1) then
            cum_q = 0.0d0
            cum_in = 0.0d0
            cum_out = 0.0d0
            cum_div = 0.0d0
            cum_p = 0.0d0
            cum_alpha = 0.0d0
            cum_dir = 0.0d0
            cum_eps = 0.0d0

            open(21, file='output/MASS_BALANCE.TXT', status='REPLACE')
            write(21,'(a)')
     &      '# t M_n M_np1 dM M_Q M_out M_in M_div M_P'
     &      //' M_alpha M_dir eps_num Cum_Q Cum_out Cum_in Cum_div'
     &      //' Cum_P Cum_alpha Cum_dir Cum_eps'
         endif

         tall1 = omp_get_wtime()

         do while (i_czas .lt. max_time)

            call przest3d(i_czas, xx, yy, zz, iss, kod, t, nax, nay,
     &                    naz, bb, il_iteracji, ttol, tnw2,
     &                    imass_balance, mb_dir)

            if (imass_balance .eq. 1) then

c------------------------------------------------------------
c  FULL MASS BALANCE AFTER EACH TIME STEP
c
c  dM = M_Q - M_out + M_in + M_div - M_P - M_alpha
c    - M_dir + eps_num
c
c  All rate terms are integrated in time with theta=2/3:
c       M_X = dt/3 * (X^n + 2*X^(n+1))
c
c  M_alpha is the signed net Robin removal rate integral
c       integral (alpha*S - q_s) dGamma.
c
c  M_dir is the algebraic reaction of the homogeneous Dirichlet
c  constraint, evaluated from the full unmasked EBE operator.
c------------------------------------------------------------
               call integrate_mass_hex8(nelem, iss, xx, yy, zz,
     &                                  tnw_prev, mass_old)
               call integrate_mass_hex8(nelem, iss, xx, yy, zz,
     &                                  tnw2, mass)

               call mass_balance_rates(i_czas-i_d_czas, nelem, iss,
     &              kod, xx, yy, zz, tnw_prev, nax, nay, naz, rq_n,
     &              rin_n, rout_n, rdiv_n, rp_n, ra_n)

               call mass_balance_rates(i_czas, nelem, iss, kod,
     &              xx, yy, zz, tnw2, nax, nay, naz, rq_np1,
     &              rin_np1, rout_np1, rdiv_np1, rp_np1,
     &              ra_np1)

               dtb = dble(i_d_czas)
               mb_q     = dtb/3.0d0*(rq_n + 2.0d0*rq_np1)
               mb_in    = dtb/3.0d0*(rin_n + 2.0d0*rin_np1)
               mb_out   = dtb/3.0d0*(rout_n + 2.0d0*rout_np1)
               mb_div   = dtb/3.0d0*(rdiv_n + 2.0d0*rdiv_np1)
               mb_p     = dtb/3.0d0*(rp_n + 2.0d0*rp_np1)
               mb_alpha = dtb/3.0d0*(ra_n + 2.0d0*ra_np1)

               dm = mass - mass_old
               eps_num = dm - (mb_q - mb_out + mb_in + mb_div
     &                         - mb_p - mb_alpha - mb_dir)

               cum_q = cum_q + mb_q
               cum_in = cum_in + mb_in
               cum_out = cum_out + mb_out
               cum_div = cum_div + mb_div
               cum_p = cum_p + mb_p
               cum_alpha = cum_alpha + mb_alpha
               cum_dir = cum_dir + mb_dir
               cum_eps = cum_eps + eps_num

               write(6,*) 'czas=', i_czas, ' masa=', mass
               write(6,*) 'MASS BALANCE [kg]: dM=', dm,
     &                    ' Q=', mb_q, ' OUT=', mb_out, ' IN=', mb_in
               write(6,*) '                     DIV=', mb_div,
     &                    ' P=', mb_p, ' ALPHA=', mb_alpha,
     &                    ' DIR=', mb_dir, ' EPS=', eps_num
               write(6,*) 'CUMULATIVE [kg]: Q=', cum_q,
     &                    ' OUT=', cum_out, ' IN=', cum_in,
     &                    ' DIV=', cum_div, ' P=', cum_p,
     &                    ' ALPHA=', cum_alpha,
     &                    ' DIR=', cum_dir, ' EPS=', cum_eps

               write(21,2100) i_czas, mass_old, mass, dm, mb_q,
     &                        mb_out, mb_in, mb_div, mb_p, mb_alpha,
     &                        mb_dir, eps_num, cum_q, cum_out, cum_in,
     &                        cum_div, cum_p, cum_alpha, cum_dir,
     &                        cum_eps
2100           format(i8,1x,19(e18.10,1x))

               do i = 1, nax*nay*naz
                  tnw_prev(i) = tnw2(i)
               enddo

            else

c--- Lightweight production diagnostic: mass in the domain only.
               call integrate_mass_hex8(nelem, iss, xx, yy, zz,
     &                                  tnw2, mass)
               write(6,*) 'czas=', i_czas, ' masa=', mass

            endif

         enddo

         if (imass_balance .eq. 1) close(21)

         tall2 = omp_get_wtime()
         print *
         print *, 'Czas calkowity = ', tall2 - tall1, ' sekund'
         print *

         tsplit1 = omp_get_wtime()
         call vanadis_files(nax, nay, naz, max_time, i_time_step)
         tsplit2 = omp_get_wtime()
         print *, 'Czas podzialu plikow = ', tsplit2 - tsplit1,
     &            ' sekund'

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
C
      subroutine velocity_field(i_czas, x, y, z, v)
C=======================================================================
C  Velocity field v(x,y,z,t) [m/s].
C
C  x is the vertical coordinate in Vanadis.  For terrain x=xg(y,z),
C  the vertical component is adjusted so that v.n=0 at the ground:
C
C       v1 = v2*dxg/dy + v3*dxg/dz .
C
C  With constant v2,v3 this field is divergence-free inside each
C  terrain cell.  Replace v2,v3 below by measured/CFD/meteo data when
C  required; terrain_geometry supplies the local terrain slopes.
C=======================================================================
         implicit real*8 (a-h,o-z)
         integer*4 i_czas
         real*8 x, y, z, v(3)
         real*8 xg, dxgdy, dxgdz
         real*8 q_cfg, k_cfg, p_cfg
         real*8 v1_cfg, v2_cfg, v3_cfg, alpha_cfg

         common /cfg/ q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg,
     &                v3_cfg, alpha_cfg

C        Horizontal velocity components.
         v(2) = v2_cfg
         v(3) = v3_cfg

C        Ground elevation and local slopes.
         call terrain_geometry(y, z, xg, dxgdy, dxgdz)

C        Terrain-consistent vertical velocity component.
         v(1) = v(2)*dxgdy + v(3)*dxgdz

         return
      end

      subroutine velocity_divergence_numeric(i_czas,x,y,z,divv)
c=======================================================================
c  Numerical divergence of v(x,y,z,t) for the optional physical mass
c  balance.  Vanadis solves the advective form v.grad(S).  Therefore,
c  for a spatially variable velocity field:
c
c    integral v.grad(S) dOmega
c      = integral_boundary S v.n dGamma - integral S div(v) dOmega.
c
c  This centered difference is used only when imass_balance=1.
c=======================================================================
         implicit real*8 (a-h,o-z)
         integer*4 i_czas
         real*8 x,y,z,divv,hx,hy,hz
         real*8 vp(3),vm(3)

         hx=dmax1(1.0d-5,1.0d-5*(1.0d0+dabs(x)))
         hy=dmax1(1.0d-5,1.0d-5*(1.0d0+dabs(y)))
         hz=dmax1(1.0d-5,1.0d-5*(1.0d0+dabs(z)))

         call velocity_field(i_czas,x+hx,y,z,vp)
         call velocity_field(i_czas,x-hx,y,z,vm)
         divv=(vp(1)-vm(1))/(2.0d0*hx)

         call velocity_field(i_czas,x,y+hy,z,vp)
         call velocity_field(i_czas,x,y-hy,z,vm)
         divv=divv+(vp(2)-vm(2))/(2.0d0*hy)

         call velocity_field(i_czas,x,y,z+hz,vp)
         call velocity_field(i_czas,x,y,z-hz,vm)
         divv=divv+(vp(3)-vm(3))/(2.0d0*hz)

         return
      end

      subroutine diffusion_field(i_czas, x, y, z, lambda)
c=======================================================================
c  Diagonal effective diffusion field K(x,y,z,t) [m2/s].
c
c  lambda(1:3) may be made anisotropic.  The default is the historical
c  isotropic K_CFG.  Values are evaluated locally at each Gauss point.
c=======================================================================
         implicit real*8 (a-h,o-z)
         integer*4 i_czas
         real*8 x, y, z, lambda(3)
         real*8 q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg, alpha_cfg
         common /cfg/ q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg,
     &                alpha_cfg

         lambda(1) = k_cfg
         lambda(2) = k_cfg
         lambda(3) = k_cfg

         return
      end

      subroutine alpha_field(i_czas, x, y, z, alfa, qs)
c=======================================================================
c  Robin/deposition field alpha(x,y,z,t) [m/s] and optional q_s.
c
c  Evaluated at every 3x3 surface Gauss point.  Default alpha is the
c  configuration value and q_s=0.
c=======================================================================
         implicit real*8 (a-h,o-z)
         integer*4 i_czas
         real*8 x, y, z, alfa, qs
         real*8 q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg, alpha_cfg
         common /cfg/ q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg,
     &                alpha_cfg

         alfa = alpha_cfg
         qs = 0.0d0

         return
      end

      subroutine transport_fields(i_czas, x, y, z, lambda, v)
         implicit real*8 (a-h,o-z)
         integer*4 i_czas
         real*8 x, y, z, lambda(3), v(3)

         call velocity_field(i_czas, x, y, z, v)
         call diffusion_field(i_czas, x, y, z, lambda)

         return
      end

c--------------------------------------------------------------------
c  Zrodlo zanieczyszczenia Q(t)c--------------------------------------------------------------------
c  Zrodlo zanieczyszczenia Q(t) [kg/s]
c
c  Wersja testowa: zrodlo stale, aktywne od t=0.
c  Q(0)=Q_CFG i Q(t)=Q_CFG dla t>0 w elemencie 1083.
c  Ta procedura jest jedynym miejscem definiujacym zaleznosc Q(t).
c--------------------------------------------------------------------
      subroutine source_term (i_czas, nr_el, q_in)

         implicit real*8 (a - h, o - z)
         real*8 q_in
         real*8 q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg, alpha_cfg
         common /cfg/ q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg,
     &                alpha_cfg

         q_in = 0.0d0

c--- stale zrodlo aktywne od t=0: Q(0)=Q_CFG
         if (nr_el .eq. 1083 .and. i_czas .ge. 0) q_in = q_cfg

         return
      end

c--------------------------------------------------------------------
c  Nieliniowy wspolczynnik zaniku P(S)
c
c  P(S) = P0 + K2*S
c
c  c_el  - srednie stezenie w elemencie S [kg/m3]
c  P0    - bazowy wspolczynnik zaniku    [1/s]
c  K2    - wspolczynnik drugiego rzedu   [m3/(kg*s)]
c  s_ref - stezenie, dla ktorego K2*S = P0 [kg/m3]
c--------------------------------------------------------------------
      subroutine licz_p (c_el, pzanik)

         implicit real*8 (a - h, o - z)
         real*8 c_el, pzanik, s_abs, s_ref, k2
         real*8 q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg, alpha_cfg
         common /cfg/ q_cfg, k_cfg, p_cfg, v1_cfg, v2_cfg, v3_cfg,
     &                alpha_cfg

c--- Wszystkie obliczenia prowadzone sa w jednostkach SI.
c--- S [kg/m3], P0 [1/s], K2 [m3/(kg*s)].
c--- Dla s_ref = 5e-6 kg/m3 (= 5 mg/m3): P(s_ref) = 2*P0.
         s_ref = 5.0d-6
         k2 = p_cfg / s_ref

c--- dabs zachowuje dodatni wspolczynnik takze dla malych
c--- numerycznych podstrzelen ponizej zera.
         s_abs = dabs(c_el)

c--- Dla fizycznego S>=0: P(S)*S = P0*S + K2*S^2
c--- P(S)*S odpowiada kinetyce P0*S + K2*S^2.
         pzanik = p_cfg + k2 * s_abs

         return
      end
C
      subroutine terrain_init(dy, dz, nay, naz)
C=======================================================================
C  Initialize terrain elevations at the horizontal mesh nodes.
C
C  xground(jy,jz) is the ground elevation xg [m] at
C       y = dy(jy), z = dz(jz).
C
C  The current block reproduces the verified 0,1,4,7 m warped-patch
C  test.  For another terrain, change only xground(:,:) here (or load
C  it once from a DEM/data file).  The FEM mesh and velocity field then
C  use exactly the same terrain through terrain_geometry.
C=======================================================================
         implicit real*8 (a-h,o-z)

         integer*4 max_ty, max_tz
         parameter (max_ty = 200, max_tz = 200)

         dimension dy(*), dz(*)
         real*8 xground(max_ty,max_tz)
         real*8 ygrid(max_ty), zgrid(max_tz)
         integer*4 nyterr, nzterr

         common /terrain_grid/ xground, ygrid, zgrid
         common /terrain_size/ nyterr, nzterr

         if (nay .gt. max_ty .or. naz .gt. max_tz) then
            write(*,*) 'ERROR terrain_init: terrain grid too large'
            stop
         endif

         nyterr = nay
         nzterr = naz

         do jy = 1, nay
            ygrid(jy) = dy(jy)
         enddo

         do jz = 1, naz
            zgrid(jz) = dz(jz)
         enddo

C        Require strictly increasing horizontal coordinates.
         do jy = 1, nay-1
            if (ygrid(jy+1) .le. ygrid(jy)) then
               write(*,*) 'ERROR terrain_init: y grid not increasing'
               stop
            endif
         enddo

         do jz = 1, naz-1
            if (zgrid(jz+1) .le. zgrid(jz)) then
               write(*,*) 'ERROR terrain_init: z grid not increasing'
               stop
            endif
         enddo

C        Default terrain: flat xg = 0 m.
         do jz = 1, naz
            do jy = 1, nay
               xground(jy,jz) = 0.0d0
            enddo
         enddo

C        Verified warped-patch test on the central 50 x 50 m cell:
C          y=-25,z=-25 : xg=0 m
C          y= 25,z=-25 : xg=1 m
C          y=-25,z= 25 : xg=4 m
C          y= 25,z= 25 : xg=7 m
C
C        These indices correspond to the current dy,dz arrays.
         if (nay .ge. 11 .and. naz .ge. 11) then
            if (dabs(ygrid(10)+25.0d0) .lt. 1.0d-12 .and.
     &          dabs(ygrid(11)-25.0d0) .lt. 1.0d-12 .and.
     &          dabs(zgrid(10)+25.0d0) .lt. 1.0d-12 .and.
     &          dabs(zgrid(11)-25.0d0) .lt. 1.0d-12) then
               xground(10,10) = 0.0d0
               xground(11,10) = 1.0d0
               xground(10,11) = 4.0d0
               xground(11,11) = 7.0d0
            else
               write(*,*) 'ERROR terrain_init: test grid mismatch'
               stop
            endif
         else
            write(*,*) 'ERROR terrain_init: test grid too small'
            stop
         endif

         return
      end
C
      subroutine terrain_geometry(y, z, xg, dxgdy, dxgdz)
C=======================================================================
C  Bilinear terrain interpolation and analytical terrain slopes.
C
C  The same nodal terrain data are used by the terrain-following mesh
C  and by velocity_field, so the ground geometry and v.n=0 condition
C  remain consistent.
C=======================================================================
         implicit real*8 (a-h,o-z)

         integer*4 max_ty, max_tz
         parameter (max_ty = 200, max_tz = 200)

         real*8 y, z, xg, dxgdy, dxgdz
         real*8 xground(max_ty,max_tz)
         real*8 ygrid(max_ty), zgrid(max_tz)
         real*8 yq, zq, y1, y2, z1, z2, r, s
         real*8 h11, h21, h12, h22
         integer*4 nyterr, nzterr

         common /terrain_grid/ xground, ygrid, zgrid
         common /terrain_size/ nyterr, nzterr

         if (nyterr .lt. 2 .or. nzterr .lt. 2) then
            write(*,*) 'ERROR terrain_geometry: terrain not initialized'
            stop
         endif

C        Clamp only for robustness at the outer domain boundary.
         yq = dmax1(ygrid(1), dmin1(y, ygrid(nyterr)))
         zq = dmax1(zgrid(1), dmin1(z, zgrid(nzterr)))

C        Locate horizontal terrain cell.
         jy = nyterr-1
         do 10 j = 1, nyterr-1
            if (yq .le. ygrid(j+1)) then
               jy = j
               goto 20
            endif
 10      continue
 20      continue

         jz = nzterr-1
         do 30 k = 1, nzterr-1
            if (zq .le. zgrid(k+1)) then
               jz = k
               goto 40
            endif
 30      continue
 40      continue

         y1 = ygrid(jy)
         y2 = ygrid(jy+1)
         z1 = zgrid(jz)
         z2 = zgrid(jz+1)

         r = (yq-y1)/(y2-y1)
         s = (zq-z1)/(z2-z1)

         h11 = xground(jy  ,jz  )
         h21 = xground(jy+1,jz  )
         h12 = xground(jy  ,jz+1)
         h22 = xground(jy+1,jz+1)

C        Bilinear ground elevation xg(y,z).
         xg = (1.0d0-r)*(1.0d0-s)*h11
     &       +r*(1.0d0-s)*h21
     &       +(1.0d0-r)*s*h12
     &       +r*s*h22

C        Analytical slopes of the bilinear surface.
         dxgdy = ((1.0d0-s)*(h21-h11)
     &            +s*(h22-h12))/(y2-y1)

         dxgdz = ((1.0d0-r)*(h12-h11)
     &            +r*(h22-h21))/(z2-z1)

         return
      end
C
      subroutine terrain_height(y, z, xg)

         implicit real*8 (a-h,o-z)
         real*8 y, z, xg, dxgdy, dxgdz

         call terrain_geometry(y, z, xg, dxgdy, dxgdz)

         return
      end
C
      subroutine siatka (xx, yy, zz, dx, dy, dz,
     &                   nax, nay, naz)
C=======================================================================
C  Terrain-following mesh.  dx(:) remains the reference vertical grid.
C  For each horizontal node, all vertical layers are mapped between the
C  local ground xg(y,z) and the unchanged top boundary xtop.
C=======================================================================
         implicit real*8 (a-h,o-z)
         dimension xx(*), yy(*), zz(*)
         dimension dx(*), dy(*), dz(*)

         real*8 xbase, xtop, sigma, xg

         xbase = dx(1)
         xtop  = dx(nax)

         if (xtop .le. xbase) then
            write(*,*) 'ERROR in siatka: invalid vertical range'
            stop
         endif

         i = 1

         do jz = 1, naz
            do jy = 1, nay

               call terrain_height(dy(jy), dz(jz), xg)

               if (xg .ge. xtop) then
                  write(*,*) 'ERROR in siatka: terrain reaches top'
                  stop
               endif

               do jx = 1, nax

                  sigma = (dx(jx)-xbase)/(xtop-xbase)

                  xx(i) = dx(jx)
     &                  + (1.0d0-sigma)*(xg-xbase)
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
         real*8 ksi_g, eta_g, zeta_g, w_g
         real*8 gp3(3), gw3(3)

         common /lok/ ksi(8), eta(8), zeta(8)
         common /lok_g/ ksi_g(27), eta_g(27), zeta_g(27), w_g(27)
         common /pow_el/ i_pow(6, 4)

         data gp3 / -0.774596669241483d0, 0.0d0,
     &               0.774596669241483d0 /
         data gw3 /  0.555555555555556d0, 0.888888888888889d0,
     &               0.555555555555556d0 /

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

c--- Gauss 3x3x3: 27 punktow i iloczynowe wagi
         j = 0
         do ig = 1, 3
            do jg = 1, 3
               do kg = 1, 3
                  j = j + 1
                  ksi_g(j) = gp3(ig)
                  eta_g(j) = gp3(jg)
                  zeta_g(j) = gp3(kg)
                  w_g(j) = gw3(ig) * gw3(jg) * gw3(kg)
               enddo
            enddo
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

      subroutine el3d(i_czas, kel1, t, kod, nax, nay, naz,
     &                c_el, c_el_n, pzanik, pzanik_n, fmv_n,
     &                hma_old, cmass_old)

         implicit real*8 (a-h,o-z)

         real*8 c_el, c_el_n, pzanik, pzanik_n
         real*8 lambda(3), lambda_old(3), v(3), v_old(3)
         real*8 n1(8), n_old(8), w(8), w_old(8)
         real*8 dnx1(8), dny1(8), dnz1(8)
         real*8 dnxw(8), dnyw(8), dnzw(8)
         real*8 dnxo(8), dnyo(8), dnzo(8)
         real*8 dnxwo(8), dnywo(8), dnzwo(8)
         real*8 djac, djac_old, detjs
         real*8 xgp, ygp, zgp, xgo, ygo, zgo
         real*8 vol, q_in, q_in_n, qv, qv_n
         real*8 fmv_n(8), hma_old(8,8), cmass_old(8,8)
         real*8 zx1(4), zx2(4), zx3(4)
         real*8 gp3(3), gw3(3), wgt, swgt
         real*8 alfaa, alfa_old, qs, qs_old
         real*8 advj, advj_old, diffij, diffij_old
         integer*4 i_time_old
         dimension t(*)
         dimension kod((nax-1)*(nay-1)*(naz-1),6)

         data gp3 / -0.774596669241483d0, 0.0d0,
     &               0.774596669241483d0 /
         data gw3 /  0.555555555555556d0, 0.888888888888889d0,
     &               0.555555555555556d0 /

         real*8 ksi_g, eta_g, zeta_g, w_g
         common /lok_g/ ksi_g(27), eta_g(27), zeta_g(27), w_g(27)
         common /wsp_lok/ xxe(8), yye(8), zze(8)
!$omp threadprivate(/wsp_lok/)
         common /pow_el/ i_pow(6,4)
         common /skl_el/ hma(8,8), cma(8,8), fmv(8), cg(8,8),
     &                   hx(8,8), hy(8,8), hz(8,8), hvx(8,8),
     &                   hvy(8,8), hvz(8,8), hmav(8,8), cma_n(8,8)
!$omp threadprivate(/skl_el/)
         common /d_time/ i_d_czas

         i_time_old = i_czas - i_d_czas

c--- Homogeneous Dirichlet values used by the current model.
         do i = 1, 6
            t(i) = 0.0d0
         enddo

         call source_term(i_czas, kel1, q_in)
         call source_term(i_time_old, kel1, q_in_n)
         call licz_p(c_el, pzanik)
         call licz_p(c_el_n, pzanik_n)

         do i = 1, 8
            fmv(i) = 0.0d0
            fmv_n(i) = 0.0d0
            do j = 1, 8
               hma(i,j) = 0.0d0
               cma(i,j) = 0.0d0
               cma_n(i,j) = 0.0d0
               hmav(i,j) = 0.0d0
               hma_old(i,j) = 0.0d0
               cmass_old(i,j) = 0.0d0
            enddo
         enddo

c--- Geometric element volume.  The current-time call also validates J.
         vol = 0.0d0
         do jl = 1, 27
            call gauss_point_hex8DR(i_czas, ksi_g(jl), eta_g(jl),
     &           zeta_g(jl), xxe, yye, zze, djac, lambda, v,
     &           n1, w, dnx1, dny1, dnz1, dnxw, dnyw, dnzw,
     &           xgp, ygp, zgp)
            vol = vol + djac*w_g(jl)
         enddo

         if (vol .le. 0.0d0) then
            write(*,*) 'ERROR in el3d: non-positive element volume.'
            write(*,*) 'element=',kel1,' volume=',vol
            stop
         endif

         qv = q_in/vol
         qv_n = q_in_n/vol

c=======================================================================
c  Assemble both time levels.
c
c  This is required when v(x,y,z,t), K(x,y,z,t), alpha(x,y,z,t)
c  change in time, because both H and the DR-weighted mass matrix M
c  depend on the local fields through W.
c=======================================================================
         do jl = 1, 27

            zss = ksi_g(jl)
            ztt = eta_g(jl)
            zuu = zeta_g(jl)
            wgt = w_g(jl)

c----------- n+1 fields and DR test function
            call gauss_point_hex8DR(i_czas, zss, ztt, zuu,
     &           xxe, yye, zze, djac, lambda, v, n1, w,
     &           dnx1, dny1, dnz1, dnxw, dnyw, dnzw,
     &           xgp, ygp, zgp)

c----------- n fields and DR test function
            call gauss_point_hex8DR(i_time_old, zss, ztt, zuu,
     &           xxe, yye, zze, djac_old, lambda_old, v_old,
     &           n_old, w_old, dnxo, dnyo, dnzo,
     &           dnxwo, dnywo, dnzwo, xgo, ygo, zgo)

            do i = 1, 8
               fmv(i) = fmv(i) + qv*w(i)*djac*wgt
               fmv_n(i) = fmv_n(i) + qv_n*w_old(i)*djac_old*wgt

               do j = 1, 8

                  diffij = dnx1(i)*lambda(1)*dnx1(j)
     &                   + dny1(i)*lambda(2)*dny1(j)
     &                   + dnz1(i)*lambda(3)*dnz1(j)
                  advj = v(1)*dnx1(j) + v(2)*dny1(j)
     &                 + v(3)*dnz1(j)

                  hma(i,j) = hma(i,j)
     &              + (diffij + w(i)*advj
     &              +  w(i)*n1(j)*pzanik)*djac*wgt

                  cma_n(i,j) = cma_n(i,j)
     &              + w(i)*n1(j)*djac*wgt

                  diffij_old = dnxo(i)*lambda_old(1)*dnxo(j)
     &                       + dnyo(i)*lambda_old(2)*dnyo(j)
     &                       + dnzo(i)*lambda_old(3)*dnzo(j)
                  advj_old = v_old(1)*dnxo(j)
     &                     + v_old(2)*dnyo(j)
     &                     + v_old(3)*dnzo(j)

                  hma_old(i,j) = hma_old(i,j)
     &              + (diffij_old + w_old(i)*advj_old
     &              +  w_old(i)*n_old(j)*pzanik_n)*djac_old*wgt

                  cmass_old(i,j) = cmass_old(i,j)
     &              + w_old(i)*n_old(j)*djac_old*wgt

               enddo
            enddo
         enddo

c=======================================================================
c  Robin boundary on face 3.  alpha and the DR weight are evaluated
c  independently at n and n+1 at every 3x3 surface Gauss point.
c=======================================================================
         i_wb = 3

         if (kod(kel1,i_wb) .eq. i_wb) then

            do i2 = 1, 4
               zx1(i2) = xxe(i_pow(i_wb,i2))
               zx2(i2) = yye(i_pow(i_wb,i2))
               zx3(i2) = zze(i_pow(i_wb,i2))
            enddo

            do ig1 = 1, 3
               do ig2 = 1, 3

                  ze1 = -1.0d0
                  ze2 = gp3(ig1)
                  ze3 = gp3(ig2)
                  swgt = gw3(ig1)*gw3(ig2)

                  call quad4_3d(zx1,zx2,zx3,ze2,ze3,detjs)

                  call gauss_point_hex8DR(i_czas, ze1, ze2, ze3,
     &                 xxe, yye, zze, djac, lambda, v, n1, w,
     &                 dnx1, dny1, dnz1, dnxw, dnyw, dnzw,
     &                 xgp, ygp, zgp)
                  call alpha_field(i_czas,xgp,ygp,zgp,alfaa,qs)

                  call gauss_point_hex8DR(i_time_old, ze1, ze2, ze3,
     &                 xxe, yye, zze, djac_old, lambda_old, v_old,
     &                 n_old, w_old, dnxo, dnyo, dnzo,
     &                 dnxwo, dnywo, dnzwo, xgo, ygo, zgo)
                  call alpha_field(i_time_old,xgo,ygo,zgo,
     &                             alfa_old,qs_old)

                  do i = 1, 8
                     fmv(i) = fmv(i) - qs*w(i)*detjs*swgt
                     fmv_n(i) = fmv_n(i)
     &                        - qs_old*w_old(i)*detjs*swgt
                     do j = 1, 8
                        hma(i,j) = hma(i,j)
     &                    + alfaa*w(i)*n1(j)*detjs*swgt
                        hma_old(i,j) = hma_old(i,j)
     &                    + alfa_old*w_old(i)*n_old(j)*detjs*swgt
                     enddo
                  enddo

               enddo
            enddo
         endif

         return
      end

      subroutine przest3d(i_czas, xx, yy, zz, iss, kod, t,
     &                     nax, nay, naz, bb, il_iteracji, ttol,
     &                     tnw2, imass_balance, mb_dir)

         implicit real*8 (a-h,o-z)

         dimension xx(*), yy(*), zz(*), t(*), bb(*)
         dimension iss((nax-1)*(nay-1)*(naz-1),8)
         dimension kod((nax-1)*(nay-1)*(naz-1),6)
         dimension tnw2(nax*nay*naz), hma_t(8,8), fmv_t(8),
     &             wekt_t(8)
         dimension tnw_n(nax*nay*naz), tnw_nl(nax*nay*naz)
         dimension hma_old(8,8), cmass_old(8,8)

         character*20 wynik_t
         integer*4 ibuf, jbuf, nbuf, lrec
         parameter (nbuf=100)
         parameter (lrec=54)
         character*54 bufor(nbuf)

c--- Expensive EBE diagonal-dominance diagnostics: 0=off, 1=on.
         integer*4 idiag_check
         parameter (idiag_check=0)

         common /pow_el/ i_pow(6,4)
         common /skl_el/ hma(8,8), cma(8,8), fmv(8), cg(8,8),
     &                   hx(8,8), hy(8,8), hz(8,8), hvx(8,8),
     &                   hvy(8,8), hvz(8,8), hmav(8,8), cma_n(8,8)
!$omp threadprivate(/skl_el/)
         common /cal_el/ ilosc_el, ilosc_wezlow, i_ot_elem(8)
!$omp threadprivate(/cal_el/)
         common /d_time/ i_d_czas

         real*8 hma2d(8,8), fmv2d(8), fmv_n(8)
         real*8 t_loc(6), c_el, c_el_n, pzanik, pzanik_n
         real*8 nl_num, nl_den, nl_err, nl_tol
         real*8 p_min, p_max
         integer*4 nl_iter, nl_max
         real*8 tmont1, tmont2, tsolv1, tsolv2, twrite1, twrite2
         real*8 omp_get_wtime
         real*8 fmv3(nax*nay*naz), diag(nax*nay*naz)
         integer*4 bc_dirichlet(nax*nay*naz)
         integer*4 imass_balance, istat
         real*8 mb_dir, rhs_dir_sum, ax_dir_sum
         logical is_diagonally_dominant, finite8

         if (i_czas .eq. 0) then
            do i = 1, ilosc_wezlow
               tnw2(i) = 0.0d0
            enddo
         endif

c--- C^n remains fixed during Picard iterations for step n+1.
         do i = 1, ilosc_wezlow
            tnw_n(i) = tnw2(i)
            tnw_nl(i) = tnw2(i)
         enddo

         nl_tol = 1.0d-8
         nl_max = 20
         nl_iter = 0
         mb_dir = 0.0d0

         i_czas = i_czas + i_d_czas

2000     continue
         nl_iter = nl_iter + 1

         do i = 1, ilosc_wezlow
            diag(i) = 0.0d0
            fmv3(i) = 0.0d0
         enddo
         do i = 1, 64*ilosc_el
            bb(i) = 0.0d0
         enddo

         print *
         write(*,199) ilosc_el, ilosc_wezlow, i_czas
199      format('matrix/vector creation...',' nr of elements: ',i6,
     &          '  nodes: ',i7,' time: ',i6,' sec')

         tmont1 = omp_get_wtime()
         p_min = huge(1.0d0)
         p_max = -huge(1.0d0)

!$omp parallel do default(shared)
!$omp& private(jk,i,j,i4,k,l,iin,hma_t,wekt_t,fmv_t,fmv2d,fmv_n,
!$omp& hma2d,hma_old,cmass_old,t_loc,c_el,c_el_n,pzanik,pzanik_n)
!$omp& reduction(+:fmv3,diag) reduction(min:p_min) reduction(max:p_max)
!$omp& copyin(/cal_el/)
         do 100 jk = 1, ilosc_el

            call data3d(jk,iss,xx,yy,zz,nax,nay,naz)

            c_el = 0.0d0
            c_el_n = 0.0d0
            do i4 = 1, 8
               c_el = c_el + tnw_nl(i_ot_elem(i4))/8.0d0
               c_el_n = c_el_n + tnw_n(i_ot_elem(i4))/8.0d0
            enddo

            call el3d(i_czas,jk,t_loc,kod,nax,nay,naz,
     &                 c_el,c_el_n,pzanik,pzanik_n,fmv_n,
     &                 hma_old,cmass_old)

            p_min = dmin1(p_min,pzanik)
            p_max = dmax1(p_max,pzanik)

c-----------------------------------------------------------------------
c theta=2/3 with time-dependent operator and DR-weighted mass matrix:
c
c (2 H^(n+1) + 3/dt M^(n+1)) C^(n+1)
c    = (3/dt M^n - H^n) C^n + 2 F^(n+1) + F^n
c-----------------------------------------------------------------------
            do i = 1, 8
               do j = 1, 8
                  hma_t(i,j) = hma_old(i,j)
     &                        - 3.0d0/i_d_czas*cmass_old(i,j)
                  hma2d(i,j) = 2.0d0*hma(i,j)
     &                        + 3.0d0/i_d_czas*cma_n(i,j)
               enddo
            enddo

            do i4 = 1, 8
               wekt_t(i4) = tnw_n(i_ot_elem(i4))
            enddo
            call ilo_m_v(hma_t,wekt_t,8,fmv_t)

            do i = 1, 8
               fmv2d(i) = 2.0d0*fmv(i) + fmv_n(i) - fmv_t(i)
            enddo

            iin = 64*(jk-1) + 1
            call sklad_wd(bb(iin),hma2d)

            do j = 1, 8
               fmv3(i_ot_elem(j)) = fmv3(i_ot_elem(j)) + fmv2d(j)
            enddo

            do i = 1, 8
               do j = 1, 8
                  k = i_ot_elem(i)
                  l = i_ot_elem(j)
                  if (k .eq. l) diag(k) = diag(k) + hma2d(i,j)
               enddo
            enddo

100      continue
!$omp end parallel do

         write(*,1113) p_min,p_max
1113     format(' P(S) range [1/s]: min=',e12.5,'  max=',e12.5)

         tmont2 = omp_get_wtime()
         print *, 'Czas montowania macierzy/wektora = ',
     &            tmont2-tmont1,' sekund'

         if (i_czas .eq. 1000)
     &      call ascii_heatmap_ppm(ilosc_el,ilosc_wezlow,iss,bb)

         if (idiag_check .eq. 1) then
            print *,'Diagonally dominant matrix [T/F] ',
     &         is_diagonally_dominant(ilosc_el,ilosc_wezlow,iss,bb)
            print *,'wspolczynnik dominacji ',
     &         diagonal_dominance_factor(ilosc_el,ilosc_wezlow,iss,bb)
         endif

         print *,'solving a system of equations...'

         tol = ttol
         il_i = il_iteracji
         non = ilosc_wezlow

         call build_bc_nodes(iss,kod,bc_dirichlet,nax,nay,naz)

c--- Preserve full RHS on constrained rows for optional mass balance.
         rhs_dir_sum = 0.0d0
         if (imass_balance .eq. 1) then
            do i = 1, ilosc_wezlow
               if (bc_dirichlet(i) .eq. 1)
     &            rhs_dir_sum = rhs_dir_sum + fmv3(i)
            enddo
         endif

c--- Homogeneous Dirichlet system mask.
         do i = 1, ilosc_wezlow
            if (bc_dirichlet(i) .eq. 1) then
               fmv3(i) = 0.0d0
               diag(i) = 1.0d0
               tnw2(i) = 0.0d0
            endif
         enddo

c--- Fail fast before entering BiCG.
         do i = 1, ilosc_wezlow
            if (.not.finite8(fmv3(i))) then
               write(*,*) 'ERROR: non-finite RHS at node ',i
               stop
            endif
            if (.not.finite8(tnw2(i))) then
               write(*,*) 'ERROR: non-finite initial X at node ',i
               stop
            endif
            if (.not.finite8(diag(i)) .or.
     &          dabs(diag(i)) .le. 1.0d-30) then
               write(*,*) 'ERROR: invalid Jacobi diagonal at node ',i
               write(*,*) 'diag=',diag(i)
               stop
            endif
         enddo

         iter = 0
         err = 0.0d0
         istat = 0
         tsolv1 = omp_get_wtime()

         call kernel_dbcg_simple(non,fmv3,tnw2,tol,il_i,iter,err,bb,
     &        iss,diag,nax,nay,naz,bc_dirichlet,istat)

c--- CPU fallback uses the same status convention:
c        call dbcg_simple(non,fmv3,tnw2,tol,il_i,iter,err,bb,
c    &        iss,diag,nax,nay,naz,bc_dirichlet,istat)

         tsolv2 = omp_get_wtime()
         print *, 'Czas solvera = ',tsolv2-tsolv1,' sekund'

         if (istat .ne. 0) then
            write(*,*) 'ERROR: DBCG failed. status=',istat,
     &                 ' iter=',iter,' err=',err
            stop
         endif
         if (.not.finite8(err)) then
            write(*,*) 'ERROR: DBCG returned non-finite residual.'
            stop
         endif
         if (err .gt. ttol) then
            write(*,*) 'ERROR: DBCG tolerance not satisfied: ',
     &                 iter,ttol,err
            stop
         endif

         write(*,1111) iter,ttol,err
1111     format(' OK - number of iterations: ',i5,' --- ',
     &          ' desired convergence tolerance: ',e10.4,
     &          ' residual norm: ',e10.4,' --- ')

c--- Picard convergence for P(S).
         nl_num = 0.0d0
         nl_den = 0.0d0
         do i = 1, ilosc_wezlow
            nl_num = nl_num + (tnw2(i)-tnw_nl(i))**2
            nl_den = nl_den + tnw2(i)**2
         enddo

         if (nl_den .gt. 1.0d-30) then
            nl_err = dsqrt(nl_num/nl_den)
         else
            nl_err = dsqrt(nl_num)
         endif

         write(*,1112) nl_iter,nl_err
1112     format(' nonlinear P(S) iteration: ',i3,
     &          '  relative change: ',e10.4)

         if (nl_err .gt. nl_tol .and. nl_iter .lt. nl_max) then
            do i = 1, ilosc_wezlow
               tnw_nl(i) = tnw2(i)
            enddo
            goto 2000
         endif

         if (nl_err .gt. nl_tol) then
            write(*,*) 'ERROR: nonlinear P(S) iteration not converged.'
            write(*,*) 'iterations=',nl_iter,' relative change=',nl_err
            stop
         endif

c--- Algebraic Dirichlet reaction for full mass-balance QA.
         mb_dir = 0.0d0
         if (imass_balance .eq. 1) then
            call dirichlet_reaction_sum(ilosc_el,iss,bb,tnw2,
     &           bc_dirichlet,ax_dir_sum)
            mb_dir = dble(i_d_czas)/3.0d0*
     &               (rhs_dir_sum-ax_dir_sum)
         endif

         write(wynik_t,105) i_czas
105      format('output/T_',i5.5,'.TXT')

         twrite1 = omp_get_wtime()
         nout10 = 10
         open(nout10,file=wynik_t,status='UNKNOWN')
         ibuf = 0

         do i = 1, ilosc_wezlow
            ibuf = ibuf + 1
            write(bufor(ibuf),122) xx(i),yy(i),zz(i),tnw2(i)*1d6
            if (ibuf .eq. nbuf) then
               write(nout10,123) (bufor(jbuf),jbuf=1,ibuf)
               ibuf = 0
            endif
         enddo

         if (ibuf .gt. 0) then
            write(nout10,123) (bufor(jbuf),jbuf=1,ibuf)
         endif
122      format(e12.5,2x,e12.5,2x,e12.5,2x,e12.5)
123      format(a54)
         close(nout10)

         twrite2 = omp_get_wtime()
         print *, 'Czas zapisu na dysk = ',twrite2-twrite1,' sekund'

         return
      end

      subroutine matvec(xin, xout, iss, bb, nax, nay, naz)
! matrix vector multiply operation Y = A*X given A and X

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

      subroutine mttvec(xin, xout, iss, bb, nax, nay, naz)
! matrix transpose vector multiply y = A'*X given A and X

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

      subroutine dbcg_simple(n,b,x,tol,itmax,iter,err,bb,iss,
     &                       diag,nax,nay,naz,bc_dirichlet,istat)
c=======================================================================
c  Robust CPU BiCG/Jacobi fallback.
c
c  ISTAT = 0  converged
c          1  maximum iterations reached
c          2  invalid/zero Jacobi diagonal
c          3  BiCG breakdown denominator
c          4  non-finite arithmetic detected
c=======================================================================
         integer*4 n,iter,itmax,nax,nay,naz,i,istat
         integer*4 bc_dirichlet(*)
         integer*4 iss((nax-1)*(nay-1)*(naz-1),8)
         real*8 tol,err,b(*),bb(*),diag(*),x(*)
         real*8 dnrm2,ak,akden,bk,bkden,bknum,bnrm,rnrm
         real*8 p(n),pp(n),r(n),rr(n),z(n),zz(n)
         real*8 break_tol
         logical finite8

         break_tol = 1.0d-30
         iter = 0
         err = 0.0d0
         istat = 0

         do i = 1, n
            if (.not.finite8(diag(i)) .or.
     &          dabs(diag(i)) .le. break_tol) then
               istat = 2
               err = huge(1.0d0)
               return
            endif
         enddo

         call matvec_0(x,r,iss,bb,nax,nay,naz,bc_dirichlet)
         do i = 1, n
            r(i) = b(i)-r(i)
            rr(i) = r(i)
         enddo

         call msolve(n,r,z,diag)
         call mtsolve(n,rr,zz,diag)

         bnrm = dnrm2(n,b)
         rnrm = dnrm2(n,r)
         if (bnrm .gt. 1.0d-30) then
            err = rnrm/bnrm
         else
            err = rnrm
         endif

         if (.not.finite8(err)) then
            istat = 4
            err = huge(1.0d0)
            return
         endif
         if (err .le. tol) return

         bkden = 0.0d0

         do 100 k = 1, itmax
            iter = k
            bknum = 0.0d0
            do i = 1, n
               bknum = bknum + z(i)*rr(i)
            enddo

            if (.not.finite8(bknum)) then
               istat = 4
               err = huge(1.0d0)
               return
            endif

            if (k .eq. 1) then
               do i = 1, n
                  p(i) = z(i)
                  pp(i) = zz(i)
               enddo
            else
               if (dabs(bkden) .le.
     &             break_tol*dmax1(1.0d0,dabs(bknum))) then
                  istat = 3
                  return
               endif
               bk = bknum/bkden
               if (.not.finite8(bk)) then
                  istat = 4
                  err = huge(1.0d0)
                  return
               endif
               do i = 1, n
                  p(i) = z(i) + bk*p(i)
                  pp(i) = zz(i) + bk*pp(i)
               enddo
            endif
            bkden = bknum

            call matvec_0(p,z,iss,bb,nax,nay,naz,bc_dirichlet)
            akden = 0.0d0
            do i = 1, n
               akden = akden + z(i)*pp(i)
            enddo

            if (.not.finite8(akden)) then
               istat = 4
               err = huge(1.0d0)
               return
            endif
            if (dabs(akden) .le.
     &          break_tol*dmax1(1.0d0,dabs(bknum))) then
               istat = 3
               return
            endif

            ak = bknum/akden
            if (.not.finite8(ak)) then
               istat = 4
               err = huge(1.0d0)
               return
            endif

            do i = 1, n
               x(i) = x(i) + ak*p(i)
               r(i) = r(i) - ak*z(i)
            enddo

            call mttvec_0(pp,zz,iss,bb,nax,nay,naz,bc_dirichlet)
            do i = 1, n
               rr(i) = rr(i) - ak*zz(i)
            enddo

            call msolve(n,r,z,diag)
            call mtsolve(n,rr,zz,diag)

            rnrm = dnrm2(n,r)
            if (bnrm .gt. 1.0d-30) then
               err = rnrm/bnrm
            else
               err = rnrm
            endif

            if (.not.finite8(err)) then
               istat = 4
               err = huge(1.0d0)
               return
            endif
            if (err .le. tol) return
100      continue

         istat = 1
         return
      end

      subroutine msolve(n,r,z,diag)
         integer*4 n,i
         real*8 z(n),r(n),diag(*)
         do i = 1, n
            z(i) = r(i)/diag(i)
         enddo
         return
      end

      subroutine mtsolve(n,rr,zz,diag)
         integer*4 n,i
         real*8 zz(n),rr(n),diag(*)
         do i = 1, n
            zz(i) = rr(i)/diag(i)
         enddo
         return
      end

      logical function finite8(x)
         implicit real*8 (a-h,o-z)
         real*8 x
         finite8 = .false.
         if (x .ne. x) return
         if (dabs(x) .gt. huge(1.0d0)) return
         finite8 = .true.
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
         print *, ' v,K,alpha are evaluated locally at Gauss points.'
         print *,' '
         print *,' Results in T_ttttt.TXT X[m],Y[m],Z[m],
     &		       S(x,y,z)[mg/m^3]'

	     print *,' '
		 print *,' Results are written to output directory'
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

      subroutine make_output_dir

         logical output_exists

         inquire(file = 'output', exist = output_exists)
         if (.not. output_exists) call system('mkdir output')

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

         character*25 dane
         character*25 dane_x, dane_y, dane_z

         integer x, y, z
         integer*4 ibuf, jbuf, ipos, ind
         integer*4 ix, iy, iz
         integer*4 nbuf, lrec
         integer*4 nout

         parameter (nbuf = 100)
         parameter (lrec = 40)

         character*40 bufor(nbuf)
         character*4200 blok

         dimension x1((nx1) * (nx2) * (nx3))
     &   , x2((nx1) * (nx2) * (nx3)), x3((nx1) * (nx2) * (nx3))
     &   , s((nx1) * (nx2) * (nx3))

         non = (nx1) * (nx2) * (nx3)

         x = nx1
         y = nx2
         z = nx3

         nout = 10

         print*
         print*,'creating 2D data files...'


125      format (e12.5, 2x, e12.5, 2x, e12.5)
126      format (a40)

222      format (e12.5, 2x, e12.5, 2x, e12.5, 2x, e12.5)


         do j_iczas = istep, iczas, istep


! --------------------------------------------------------------------
! Odczyt calego pliku 3D dla danego kroku czasowego
! --------------------------------------------------------------------

            write(dane, 108) j_iczas
108         format('output/T_', i5.5, '.TXT')

            open (9, file = dane, status = 'OLD')

            do i = 1, non
               read(9, 222) x1(i), x2(i), x3(i), s(i)
            enddo

            close(9)


! --------------------------------------------------------------------
! Plaszczyzny X = const
! --------------------------------------------------------------------

            do j = 1, x

               write(dane_x, 105) j_iczas, j
105            format('output/X_', i5.5, '_', i3.3, '.TXT')

               open (nout, file = dane_x, status = 'unknown')

               ibuf = 0

               do iz = 1, z
                  do iy = 1, y

                     ind = j + x * (iy - 1)
     &                       + x * y * (iz - 1)

                     ibuf = ibuf + 1

                     write(bufor(ibuf), 125)
     &                    x2(ind), x3(ind), s(ind)

                     if (ibuf .eq. nbuf) then

                        write(nout, 126)
     &                       (bufor(jbuf), jbuf = 1, ibuf)

                        ibuf = 0

                     endif

                  enddo
               enddo


! --- zapis pozostalosci bufora

               if (ibuf .gt. 0) then

                  write(nout, 126)
     &                 (bufor(jbuf), jbuf = 1, ibuf)

               endif

               close(nout)

            enddo


! --------------------------------------------------------------------
! Plaszczyzny Y = const
! --------------------------------------------------------------------

            do j = 1, y

               write(dane_y, 106) j_iczas, j
106            format('output/Y_', i5.5, '_', i3.3, '.TXT')

               open (nout, file = dane_y, status = 'unknown')

               ibuf = 0

               do iz = 1, z
                  do ix = 1, x

                     ind = ix + x * (j - 1)
     &                        + x * y * (iz - 1)

                     ibuf = ibuf + 1

                     write(bufor(ibuf), 125)
     &                    x3(ind), x1(ind), s(ind)

                     if (ibuf .eq. nbuf) then

                        write(nout, 126)
     &                       (bufor(jbuf), jbuf = 1, ibuf)

                        ibuf = 0

                     endif

                  enddo
               enddo


! --- zapis pozostalosci bufora

               if (ibuf .gt. 0) then

                  write(nout, 126)
     &                 (bufor(jbuf), jbuf = 1, ibuf)

               endif

               close(nout)

            enddo


! --------------------------------------------------------------------
! Plaszczyzny Z = const
! --------------------------------------------------------------------

            do j = 1, z

               write(dane_z, 107) j_iczas, j
107            format('output/Z_', i5.5, '_', i3.3, '.TXT')

               open (nout, file = dane_z, status = 'unknown')

               ibuf = 0

               do iy = 1, y
                  do ix = 1, x

                     ind = ix + x * (iy - 1)
     &                        + x * y * (j - 1)

                     ibuf = ibuf + 1

                     write(bufor(ibuf), 125)
     &                    x2(ind), x1(ind), s(ind)

                     if (ibuf .eq. nbuf) then

                        write(nout, 126)
     &                       (bufor(jbuf), jbuf = 1, ibuf)

                        ibuf = 0

                     endif

                  enddo
               enddo


! --- zapis pozostalosci bufora

               if (ibuf .gt. 0) then

                  write(nout, 126)
     &                 (bufor(jbuf), jbuf = 1, ibuf)

               endif

               close(nout)

            enddo


         enddo

         return
      end
	  
      subroutine vanadis_files_old(nx1, nx2, nx3, iczas, istep)

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
105      format('output/X_', i5.5, '_', i3.3, '.TXT')

         do j = istep, iczas, istep
            do i = 1, y
               write(dane_y(i, j), 106) j, i
            enddo
         enddo
106      format('output/Y_', i5.5, '_', i3.3, '.TXT')

         do j = istep, iczas, istep
            do i = 1, z
               write(dane_z(i, j), 107) j, i
            enddo
         enddo
107      format('output/Z_', i5.5, '_', i3.3, '.TXT')


125      format (e12.5, 2x, e12.5, 2x, e12.5)

         do j = 1, iczas
            write(dane(j), 108) j
         enddo
108      format('output/T_', i5.5, '.TXT')

         do j_iczas = istep, iczas, istep

            open (9, file = dane(j_iczas), status = 'OLD')
            do i = 1, non
               read(9, 222) x1(i), x2(i), x3(i), s(i)
            enddo
            close(9)
222         format (e12.5, 2x, e12.5, 2x, e12.5, 2x, e12.5)

            do j = 1, x
               open (30, file = dane_x(j, j_iczas), status = 'unknown')
               do i = 1, (y * z)
                  write(30, 125) x2(x * (i - 1) + 1 + j - 1),
     &            x3(x * (i - 1) + 1 + j - 1), s(x * (i - 1) + 1 + j -
     &            1)
               enddo
               close(30)
            enddo

            do j = 1, y
               k = 1
               open (30, file = dane_y(j, j_iczas), status = 'unknown')
               do i = 1, (x * z)
                  if (mod(i - 1, x) .eq. 0 .and. (i - 1) .ne. 0) k = k
     &            + 1
                  write(30, 125) x3(i + x * (j - 1) +x * (y - 1) * (k -
     &            1))
     &            , x1(i + x * (j - 1) +x * (y - 1) * (k - 1))
     &            , s(i + x * (j - 1) +x * (y - 1) * (k - 1))
               enddo
               close(30)
            enddo

            do j = 1, z
               open (30, file = dane_z(j, j_iczas), status = 'unknown')
               do i = 1, (x * y)
                  write(30, 125) x2(i + x * y * (j - 1))
     &            , x1(i + x * y * (j - 1)), s(i + x * y * (j - 1))
               enddo
               close(30)
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
     &                      dw_dxi)
!  8-node brick shape functions and derivatives wrt (xi,eta,zeta)
!  W tej wersji DR: w i dw_dxi zostaną nadpisane w gauss_point_hex8DR,
!  tutaj liczymy tylko N i dN/dxi, dN/deta, dN/dzeta.

         real*8 xi2, eta2, zeta2
         real*8 n(8), dn_dxi(3, 8)
         real*8 xm, xp, ym, yp, zm, zp
         real*8 w(8), dw_dxi(3, 8)

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
	  
	  

      subroutine gauss_point_hex8DR_base(xi, eta, zeta, xnod, ynod,
     &                                    znod,
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
      subroutine gauss_point_hex8DR(i_czas,xi,eta,zeta,
     & xnod,ynod,znod,detj,lambda,v,n,w,dnx1,dny1,dnz1,
     & dnxw,dnyw,dnzw,xgp,ygp,zgp)
c=======================================================================
c  Local HEX8/DR evaluation with pointwise fields.
c
c  v(x,y,z,t) and K(x,y,z,t) are evaluated at the physical Gauss point.
c  DR:
c       W_i = N_i + tau_DR * (v . grad N_i)
c       h_stream = 2 / ||J^{-1} e||, e=v/|v|
c       lambda_s = e^T Lambda e
c       tau_DR = h/(2|v|) * (coth(Pe)-1/Pe)
c  and for K=0: tau_DR -> h/(2|v|).
c=======================================================================
         implicit real*8 (a-h,o-z)
         integer*4 i_czas
         real*8 xi,eta,zeta,xnod(8),ynod(8),znod(8)
         real*8 n(8),detj,dn_dxi(3,8),j(3,3),invj(3,3)
         real*8 lambda(3),v(3),w(8),dw_dxi(3,8)
         real*8 dnx1(8),dny1(8),dnz1(8)
         real*8 dnxw(8),dnyw(8),dnzw(8)
         real*8 vmag,h_stream,lambda_s,pe,tau_dr,betape
         real*8 ex,ey,ez,gx,gy,gz,gnorm,vdotg,detji,advscale
         real*8 jnorm1,jnorm2,jnorm3,jscale,jac_quality,jac_tol
         real*8 jnormf,invnormf,jac_cond,jac_cond_max
         real*8 eps_v,eps_g,xgp,ygp,zgp
         integer*4 i,jj

         call shape_hex8DR(xi,eta,zeta,n,dn_dxi,w,dw_dxi)

c--- Physical location of this Gauss point.
         xgp = 0.0d0
         ygp = 0.0d0
         zgp = 0.0d0
         do i = 1, 8
            xgp = xgp + n(i)*xnod(i)
            ygp = ygp + n(i)*ynod(i)
            zgp = zgp + n(i)*znod(i)
         enddo

         call transport_fields(i_czas,xgp,ygp,zgp,lambda,v)

c--- Every diagonal diffusion component must be non-negative.
         do i = 1, 3
            if (lambda(i) .lt. 0.0d0) then
               write(*,*) 'ERROR in gauss_point_hex8DR:'
               write(*,*) 'Negative diffusion component ',i,lambda(i)
               write(*,*) 'x,y,z,t=',xgp,ygp,zgp,i_czas
               stop
            endif
         enddo

c--- J = d(x,y,z)/d(xi,eta,zeta)
         do i = 1, 3
            do jj = 1, 3
               j(i,jj) = 0.0d0
            enddo
         enddo
         do i = 1, 8
            j(1,1)=j(1,1)+xnod(i)*dn_dxi(1,i)
            j(1,2)=j(1,2)+xnod(i)*dn_dxi(2,i)
            j(1,3)=j(1,3)+xnod(i)*dn_dxi(3,i)
            j(2,1)=j(2,1)+ynod(i)*dn_dxi(1,i)
            j(2,2)=j(2,2)+ynod(i)*dn_dxi(2,i)
            j(2,3)=j(2,3)+ynod(i)*dn_dxi(3,i)
            j(3,1)=j(3,1)+znod(i)*dn_dxi(1,i)
            j(3,2)=j(3,2)+znod(i)*dn_dxi(2,i)
            j(3,3)=j(3,3)+znod(i)*dn_dxi(3,i)
         enddo

         detj = j(1,1)*(j(2,2)*j(3,3)-j(2,3)*j(3,2))
     &        - j(1,2)*(j(2,1)*j(3,3)-j(2,3)*j(3,1))
     &        + j(1,3)*(j(2,1)*j(3,2)-j(2,2)*j(3,1))

         jnorm1 = dsqrt(j(1,1)**2+j(1,2)**2+j(1,3)**2)
         jnorm2 = dsqrt(j(2,1)**2+j(2,2)**2+j(2,3)**2)
         jnorm3 = dsqrt(j(3,1)**2+j(3,2)**2+j(3,3)**2)
         jscale = jnorm1*jnorm2*jnorm3
         jac_tol = 1.0d-12

         if (jscale .le. 0.0d0) then
            write(*,*) 'ERROR in gauss_point_hex8DR: degenerate J.'
            write(*,*) 'xi,eta,zeta=',xi,eta,zeta,' detJ=',detj
            stop
         endif

         jac_quality = dabs(detj)/jscale
         if (detj .le. 0.0d0) then
            write(*,*) 'ERROR in gauss_point_hex8DR: detJ <= 0.'
            write(*,*) 'xi,eta,zeta=',xi,eta,zeta,' detJ=',detj
            write(*,*) 'J quality=',jac_quality
            stop
         endif
         if (jac_quality .lt. jac_tol) then
            write(*,*) 'ERROR in gauss_point_hex8DR: nearly singular J.'
            write(*,*) 'J quality=',jac_quality
            stop
         endif

c--- Inverse Jacobian.
         detji = 1.0d0/detj
         invj(1,1)=detji*(j(2,2)*j(3,3)-j(2,3)*j(3,2))
         invj(1,2)=detji*(j(1,3)*j(3,2)-j(1,2)*j(3,3))
         invj(1,3)=detji*(j(1,2)*j(2,3)-j(1,3)*j(2,2))
         invj(2,1)=detji*(j(2,3)*j(3,1)-j(2,1)*j(3,3))
         invj(2,2)=detji*(j(1,1)*j(3,3)-j(1,3)*j(3,1))
         invj(2,3)=detji*(j(1,3)*j(2,1)-j(1,1)*j(2,3))
         invj(3,1)=detji*(j(2,1)*j(3,2)-j(2,2)*j(3,1))
         invj(3,2)=detji*(j(1,2)*j(3,1)-j(1,1)*j(3,2))
         invj(3,3)=detji*(j(1,1)*j(2,2)-j(1,2)*j(2,1))

c--- Conservative condition-number guard (Frobenius estimate).
         jnormf = 0.0d0
         invnormf = 0.0d0
         do i = 1, 3
            do jj = 1, 3
               jnormf = jnormf + j(i,jj)*j(i,jj)
               invnormf = invnormf + invj(i,jj)*invj(i,jj)
            enddo
         enddo
         jac_cond = dsqrt(jnormf)*dsqrt(invnormf)
         jac_cond_max = 1.0d12
         if (jac_cond .gt. jac_cond_max) then
            write(*,*) 'ERROR in gauss_point_hex8DR: ill-conditioned J.'
            write(*,*) 'cond_F(J)=',jac_cond,' limit=',jac_cond_max
            write(*,*) 'x,y,z,t=',xgp,ygp,zgp,i_czas
            stop
         endif

c--- Physical gradients: grad_x(N)=J^{-T} grad_xi(N).
         do i = 1, 8
            dnx1(i)=invj(1,1)*dn_dxi(1,i)
     &             +invj(2,1)*dn_dxi(2,i)
     &             +invj(3,1)*dn_dxi(3,i)
            dny1(i)=invj(1,2)*dn_dxi(1,i)
     &             +invj(2,2)*dn_dxi(2,i)
     &             +invj(3,2)*dn_dxi(3,i)
            dnz1(i)=invj(1,3)*dn_dxi(1,i)
     &             +invj(2,3)*dn_dxi(2,i)
     &             +invj(3,3)*dn_dxi(3,i)
         enddo

         eps_v = 1.0d-14
         eps_g = 1.0d-14
         vmag = dsqrt(v(1)**2+v(2)**2+v(3)**2)
         if (vmag .gt. eps_v) then
            ex=v(1)/vmag
            ey=v(2)/vmag
            ez=v(3)/vmag
         else
            ex=0.0d0
            ey=0.0d0
            ez=0.0d0
         endif

         gx=invj(1,1)*ex+invj(1,2)*ey+invj(1,3)*ez
         gy=invj(2,1)*ex+invj(2,2)*ey+invj(2,3)*ez
         gz=invj(3,1)*ex+invj(3,2)*ey+invj(3,3)*ez
         gnorm=dsqrt(gx*gx+gy*gy+gz*gz)

         if (vmag .gt. eps_v .and. gnorm .gt. eps_g) then
            h_stream=2.0d0/gnorm
         else
            h_stream=0.0d0
         endif

         lambda_s=lambda(1)*ex*ex+lambda(2)*ey*ey+lambda(3)*ez*ez

         tau_dr=0.0d0
         pe=0.0d0
         betape=0.0d0
         if (vmag .gt. eps_v .and. h_stream .gt. 0.0d0) then
            advscale=0.5d0*vmag*h_stream
            if (lambda_s .eq. 0.0d0) then
               betape=1.0d0
               tau_dr=0.5d0*h_stream/vmag
            else if (lambda_s .lt. advscale/50.0d0) then
               betape=1.0d0-lambda_s/advscale
               tau_dr=0.5d0*h_stream/vmag*betape
            else
               pe=advscale/lambda_s
               if (pe .lt. 1.0d-3) then
                  betape=pe/3.0d0-pe**3/45.0d0
     &                  +2.0d0*pe**5/945.0d0
               else
                  betape=1.0d0/dtanh(pe)-1.0d0/pe
               endif
               tau_dr=0.5d0*h_stream/vmag*betape
            endif
         endif

         do i = 1, 8
            vdotg=v(1)*dnx1(i)+v(2)*dny1(i)+v(3)*dnz1(i)
            w(i)=n(i)+tau_dr*vdotg
            dnxw(i)=dnx1(i)
            dnyw(i)=dny1(i)
            dnzw(i)=dnz1(i)
         enddo

         return
      end
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

c--- Ściana 1: zeta = -1
      if (iwb .eq. 1) then
         bc_dirichlet(nodes(1)) = 1
         bc_dirichlet(nodes(2)) = 1
         bc_dirichlet(nodes(3)) = 1
         bc_dirichlet(nodes(4)) = 1
      endif

c--- Ściana 2: zeta = +1
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

c--- Ściana 5: eta = +1
      if (iwb .eq. 5) then
         bc_dirichlet(nodes(4)) = 1
         bc_dirichlet(nodes(3)) = 1
         bc_dirichlet(nodes(7)) = 1
         bc_dirichlet(nodes(8)) = 1
      endif

c--- Ściana 6: eta = -1
      if (iwb .eq. 6) then
         bc_dirichlet(nodes(1)) = 1
         bc_dirichlet(nodes(2)) = 1
         bc_dirichlet(nodes(6)) = 1
         bc_dirichlet(nodes(5)) = 1
      endif

      return
      end
c
c======================================================================
c  Zaznaczanie wezlow sciany na podstawie numeru sciany iwb
c  Konwencja hexa8 (lokalne numery wezlow):
c
c      eta
c       ^
c       |
c       4--------3
c      /|       /|
c     8--------7 |
c     | |      | |
c     | 1------|-2 --> ksi
c     |/       |/
c     5--------6
c
c  Kierunki lokalne:
c    1 -> 2 : ksi  = -1 -> +1
c    1 -> 4 : eta  = -1 -> +1
c    1 -> 5 : zeta = -1 -> +1
c
c  Sciany:
c    iwb = 1 : zeta = -1 -> (1,2,3,4)
c    iwb = 2 : zeta = +1 -> (5,6,7,8)
c    iwb = 3 : ksi  = -1 -> (1,4,8,5)
c    iwb = 4 : ksi  = +1 -> (2,3,7,6)
c    iwb = 5 : eta  = +1 -> (4,3,7,8)
c    iwb = 6 : eta  = -1 -> (1,2,6,5)
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
	   
	   
	   
c=======================================================================
c  Sum A*x on Dirichlet rows of the full unmasked EBE operator.
c=======================================================================
      subroutine dirichlet_reaction_sum(ilosc_el,iss,bb,x,
     &                                  bc_dirichlet,ax_dir_sum)
         implicit real*8 (a-h,o-z)
         integer*4 ilosc_el,iss(ilosc_el,8),bc_dirichlet(*)
         real*8 bb(*),x(*),ax_dir_sum
         integer*4 je,i,j,ii,jj,ik

         ax_dir_sum=0.0d0
         ik=1
         do je=1,ilosc_el
            do i=1,8
               ii=iss(je,i)
               if (bc_dirichlet(ii) .eq. 1) then
                  do j=1,8
                     jj=iss(je,j)
                     ax_dir_sum=ax_dir_sum+bb(ik)*x(jj)
                     ik=ik+1
                  enddo
               else
                  ik=ik+8
               endif
            enddo
         enddo
         return
      end

c=======================================================================
c  Instantaneous physical mass-balance rates.
c
c  3x3x3 volume integration and 3x3 surface integration are used.
c  v and alpha are evaluated pointwise as v(x,y,z,t), alpha(x,y,z,t).
c  M_div = integral S*div(v) dOmega makes the diagnostic valid also
c  for a non-divergence-free spatial velocity field.
c=======================================================================
      subroutine mass_balance_rates(i_czas,nelem,conn,kod,
     &           xx,yy,zz,c,nax,nay,naz,rq_rate,
     &           rin_rate,rout_rate,rdiv_rate,rp_rate,ra_rate)
         implicit real*8 (a-h,o-z)
         integer*4 i_czas,nelem,nax,nay,naz
         integer*4 conn(nelem,8)
         integer*4 kod((nax-1)*(nay-1)*(naz-1),6)
         real*8 xx(*),yy(*),zz(*),c(*)
         real*8 rq_rate,rin_rate,rout_rate,rdiv_rate
         real*8 rp_rate,ra_rate
         real*8 gp(3),gw(3),n(8),dn_dxi(3,8),jmat(3,3),detj
         real*8 xc(8),yc(8),zc(8),cc(8),v(3)
         real*8 q_in,pzanik,c_el,c_val,xi,eta,zeta,r,s
         real*8 sx,sy,sz,area,vdotda,flux,alfa,qs,wgt
         real*8 xgp,ygp,zgp,divv
         integer*4 e,a,ig,jg,kg,iwb

         data gp / -0.774596669241483d0,0.0d0,
     &              0.774596669241483d0 /
         data gw /  0.555555555555556d0,0.888888888888889d0,
     &              0.555555555555556d0 /

         rq_rate=0.0d0
         rin_rate=0.0d0
         rout_rate=0.0d0
         rdiv_rate=0.0d0
         rp_rate=0.0d0
         ra_rate=0.0d0

         do 100 e=1,nelem
            c_el=0.0d0
            do a=1,8
               xc(a)=xx(conn(e,a))
               yc(a)=yy(conn(e,a))
               zc(a)=zz(conn(e,a))
               cc(a)=c(conn(e,a))
               c_el=c_el+cc(a)/8.0d0
            enddo

            call source_term(i_czas,e,q_in)
            rq_rate=rq_rate+q_in
            call licz_p(c_el,pzanik)

c----------- Volume reaction integral P(S_e)*S.
            do ig=1,3
               xi=gp(ig)
               do jg=1,3
                  eta=gp(jg)
                  do kg=1,3
                     zeta=gp(kg)
                     wgt=gw(ig)*gw(jg)*gw(kg)
                     call shape_hex8_only_n(xi,eta,zeta,n,dn_dxi)
                     call jacobian_hex8(dn_dxi,xc,yc,zc,jmat,detj)
                     c_val=0.0d0
                     xgp=0.0d0
                     ygp=0.0d0
                     zgp=0.0d0
                     do a=1,8
                        c_val=c_val+n(a)*cc(a)
                        xgp=xgp+n(a)*xc(a)
                        ygp=ygp+n(a)*yc(a)
                        zgp=zgp+n(a)*zc(a)
                     enddo
                     call velocity_divergence_numeric(i_czas,xgp,ygp,
     &                                                zgp,divv)
                     rdiv_rate=rdiv_rate+c_val*divv*detj*wgt
                     rp_rate=rp_rate+pzanik*c_val*detj*wgt
                  enddo
               enddo
            enddo

c----------- External faces: advective IN/OUT and Robin net removal.
            do iwb=1,6
               if (kod(e,iwb) .eq. iwb) then
                  do ig=1,3
                     r=gp(ig)
                     do jg=1,3
                        s=gp(jg)
                        wgt=gw(ig)*gw(jg)

                        if (iwb .eq. 1) then
                           xi=r
                           eta=s
                           zeta=-1.0d0
                        else if (iwb .eq. 2) then
                           xi=r
                           eta=s
                           zeta=1.0d0
                        else if (iwb .eq. 3) then
                           xi=-1.0d0
                           eta=r
                           zeta=s
                        else if (iwb .eq. 4) then
                           xi=1.0d0
                           eta=r
                           zeta=s
                        else if (iwb .eq. 5) then
                           xi=r
                           eta=1.0d0
                           zeta=s
                        else
                           xi=r
                           eta=-1.0d0
                           zeta=s
                        endif

                        call shape_hex8_only_n(xi,eta,zeta,n,dn_dxi)
                        call jacobian_hex8(dn_dxi,xc,yc,zc,jmat,detj)

                        c_val=0.0d0
                        xgp=0.0d0
                        ygp=0.0d0
                        zgp=0.0d0
                        do a=1,8
                           c_val=c_val+n(a)*cc(a)
                           xgp=xgp+n(a)*xc(a)
                           ygp=ygp+n(a)*yc(a)
                           zgp=zgp+n(a)*zc(a)
                        enddo

                        call velocity_field(i_czas,xgp,ygp,zgp,v)

c---------------------- Outward oriented dA/d(local coordinates).
                        if (iwb .eq. 1 .or. iwb .eq. 2) then
                           sx=jmat(2,1)*jmat(3,2)
     &                       -jmat(3,1)*jmat(2,2)
                           sy=jmat(3,1)*jmat(1,2)
     &                       -jmat(1,1)*jmat(3,2)
                           sz=jmat(1,1)*jmat(2,2)
     &                       -jmat(2,1)*jmat(1,2)
                           if (iwb .eq. 1) then
                              sx=-sx
                              sy=-sy
                              sz=-sz
                           endif
                        else if (iwb .eq. 3 .or. iwb .eq. 4) then
                           sx=jmat(2,2)*jmat(3,3)
     &                       -jmat(3,2)*jmat(2,3)
                           sy=jmat(3,2)*jmat(1,3)
     &                       -jmat(1,2)*jmat(3,3)
                           sz=jmat(1,2)*jmat(2,3)
     &                       -jmat(2,2)*jmat(1,3)
                           if (iwb .eq. 3) then
                              sx=-sx
                              sy=-sy
                              sz=-sz
                           endif
                        else
                           sx=jmat(2,3)*jmat(3,1)
     &                       -jmat(3,3)*jmat(2,1)
                           sy=jmat(3,3)*jmat(1,1)
     &                       -jmat(1,3)*jmat(3,1)
                           sz=jmat(1,3)*jmat(2,1)
     &                       -jmat(2,3)*jmat(1,1)
                           if (iwb .eq. 6) then
                              sx=-sx
                              sy=-sy
                              sz=-sz
                           endif
                        endif

                        vdotda=v(1)*sx+v(2)*sy+v(3)*sz
                        flux=c_val*vdotda*wgt
                        if (flux .ge. 0.0d0) then
                           rout_rate=rout_rate+flux
                        else
                           rin_rate=rin_rate-flux
                        endif

                        if (iwb .eq. 3) then
                           area=dsqrt(sx*sx+sy*sy+sz*sz)
                           call alpha_field(i_czas,xgp,ygp,zgp,
     &                                      alfa,qs)
                           ra_rate=ra_rate+
     &                       (alfa*c_val-qs)*area*wgt
                        endif

                     enddo
                  enddo
               endif
            enddo
100      continue

         return
      end

      subroutine integrate_mass_hex8(nelem, conn, xx, yy, zz, c,
     &                                mass)
c-----------------------------------------------------------------
c  Oblicza masę M = ∫ c(x) dΩ dla elementów HEX8 (3x3x3 Gauss)
c-----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer nelem, conn(nelem,8)
      real*8 xx(*), yy(*), zz(*), c(*)
      real*8 mass

      real*8 xi_g(3), w_g(3)
      real*8 xi, eta, zeta, wgt
      real*8 n(8), dN_dxi(3,8)
      real*8 J(3,3), detJ
      real*8 xc(8), yc(8), zc(8), cc(8)
      real*8 c_val
      integer e, a, ig, jg, kg

c--- Gauss 3x3x3
      xi_g(1) = -0.774596669241483d0
      xi_g(2) =  0.0d0
      xi_g(3) =  0.774596669241483d0
      w_g(1)  =  0.555555555555556d0
      w_g(2)  =  0.888888888888889d0
      w_g(3)  =  0.555555555555556d0

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

c-------- Gauss 3x3x3
         do 200 ig = 1, 3
            xi = xi_g(ig)
            do 210 jg = 1, 3
               eta = xi_g(jg)
               do 220 kg = 1, 3
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
C
      SUBROUTINE  WB_1(T,II,HMA,FMV)

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION HMA(8,8),FMV(8)

      I=II
      DO 10 J=1,8
        HMA(I,I)=1
        FMV(I)=T
        IF (I.EQ.J) GO TO 10
        FMV(J)=FMV(J)-HMA(J,I)*T
        HMA(I,J)=0
        HMA(J,I)=0
 10   CONTINUE

      RETURN
      END
C	  