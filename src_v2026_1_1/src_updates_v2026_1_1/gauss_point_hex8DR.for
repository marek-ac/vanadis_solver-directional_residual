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
		 
		 real*8 gx, gy, gz, gnorm

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

C! --- Długości elementu w osiach (zakładamy prostopadłościan)
C
C         dlx = dabs(xnod(2) - xnod(1))
C         dly = dabs(ynod(4) - ynod(1))
C         dlz = dabs(znod(5) - znod(1))
C
C! --- Prędkość i Peclet wzdłuż kierunku v
C         eps   = 1.0d-14
C         vmag  = dsqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
C
C! efektywna długość wzdłuż strugi (klasyczna definicja 1/h = sum |v_i|/h_i / |v|)
C         if (vmag .gt. eps) then
C            h_stream = 2.0d0 / ( dabs(v(1))/dlx + dabs(v(2))/dly
C     &                          + dabs(v(3))/dlz + eps )
C         else
C            h_stream = 0.0d0
C         endif


! --- Kierunek jednostkowy prędkości
        eps   = 1.0d-14
        vmag  = dsqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
        
        if (vmag .gt. eps) then
           ex = v(1) / vmag
           ey = v(2) / vmag
           ez = v(3) / vmag
        else
           ex = 0.0d0
           ey = 0.0d0
           ez = 0.0d0
        endif
        
        ! --- Kierunek prędkości w lokalnych współrzędnych (xi,eta,zeta):
        !     g = J^{-1} * v_hat
!        real*8 gx, gy, gz, gnorm
        
        gx = invj(1,1)*ex + invj(1,2)*ey + invj(1,3)*ez
        gy = invj(2,1)*ex + invj(2,2)*ey + invj(2,3)*ez
        gz = invj(3,1)*ex + invj(3,2)*ey + invj(3,3)*ez
        
        gnorm = dsqrt(gx*gx + gy*gy + gz*gz)
        
        if (gnorm .gt. eps) then
           h_stream = 2.0d0 / gnorm
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
