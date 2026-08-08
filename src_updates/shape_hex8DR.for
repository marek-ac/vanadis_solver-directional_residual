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
	  