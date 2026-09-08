      subroutine gauss_point_hex8DR(xi, eta, zeta, xnod, ynod, znod,
     & detj, lambda, v, n, w, dnx1, dny1, dnz1,
     & dnxw, dnyw, dnzw)

c=======================================================================
c  Evaluate:
c    N, dN/dx, dN/dy, dN/dz,
c    directional-residual test function W,
c    det(J)
c
c  DR test function:
c
c       W_i = N_i + tau_DR * (v . grad N_i)
c
c  Streamwise element length:
c
c       h_stream = 2 / || J^{-1} e ||
c
c       e = v / |v|
c
c  Directional diffusion:
c
c       lambda_s = e^T Lambda e
c
c  Streamwise Peclet number:
c
c       Pe = |v| h_stream / (2 lambda_s)
c
c  Stabilization parameter:
c
c       tau_DR = h_stream/(2|v|)
c    &           * (coth(Pe) - 1/Pe)
c
c  Pure-advection limit:
c
c       lambda_s -> 0
c       Pe       -> infinity
c       tau_DR   -> h_stream/(2|v|)
c
c  Physical diffusion remains Galerkin:
c       grad(N) is used in the diffusion term.
c       The DR perturbation is not differentiated there.
c=======================================================================

         implicit real*8 (a-h,o-z)

         real*8 xi, eta, zeta
         real*8 xnod(8), ynod(8), znod(8)

         real*8 n(8), detj
         real*8 dn_dxi(3,8)

         real*8 j(3,3), invj(3,3)

         real*8 lambda(3), v(3)

         real*8 w(8), dw_dxi(3,8)

         real*8 dnx1(8), dny1(8), dnz1(8)
         real*8 dnxw(8), dnyw(8), dnzw(8)

         real*8 vmag
         real*8 h_stream
         real*8 lambda_s
         real*8 pe, tau_DR, betaPe

         real*8 ex, ey, ez
         real*8 gx, gy, gz, gnorm
         real*8 vdotg

         real*8 detji

c----- Jacobian quality control
         real*8 jnorm1, jnorm2, jnorm3
         real*8 jscale, jac_quality
         real*8 jac_tol

c----- Advection scale = |v| h_stream / 2
         real*8 advscale

c----- Numerical threshold only for velocity / metric tests
         real*8 eps_v, eps_g

         integer*4 i, jj


c=======================================================================
c  Classical HEX8 shape functions and local derivatives
c=======================================================================

         call shape_hex8DR(xi, eta, zeta, n, dn_dxi, w, dw_dxi,
     &                     lambda, v)


c=======================================================================
c  Jacobian
c
c       J = d(x,y,z) / d(xi,eta,zeta)
c=======================================================================

         do i = 1, 3
            do jj = 1, 3
               j(i,jj) = 0.0d0
            enddo
         enddo

         do i = 1, 8

            j(1,1) = j(1,1) + xnod(i)*dn_dxi(1,i)
            j(1,2) = j(1,2) + xnod(i)*dn_dxi(2,i)
            j(1,3) = j(1,3) + xnod(i)*dn_dxi(3,i)

            j(2,1) = j(2,1) + ynod(i)*dn_dxi(1,i)
            j(2,2) = j(2,2) + ynod(i)*dn_dxi(2,i)
            j(2,3) = j(2,3) + ynod(i)*dn_dxi(3,i)

            j(3,1) = j(3,1) + znod(i)*dn_dxi(1,i)
            j(3,2) = j(3,2) + znod(i)*dn_dxi(2,i)
            j(3,3) = j(3,3) + znod(i)*dn_dxi(3,i)

         enddo


c=======================================================================
c  det(J)
c=======================================================================

         detj = j(1,1)*(j(2,2)*j(3,3) - j(2,3)*j(3,2))
     &        - j(1,2)*(j(2,1)*j(3,3) - j(2,3)*j(3,1))
     &        + j(1,3)*(j(2,1)*j(3,2) - j(2,2)*j(3,1))


c=======================================================================
c  Jacobian validity check
c
c  Instead of testing det(J) against an arbitrary dimensional number,
c  use a dimensionless quality measure:
c
c       q_J = |det(J)| /
c             ( ||J_row1|| ||J_row2|| ||J_row3|| )
c
c  q_J -> 0 means nearly singular geometry.
c=======================================================================

         jnorm1 = dsqrt(j(1,1)*j(1,1)
     &                + j(1,2)*j(1,2)
     &                + j(1,3)*j(1,3))

         jnorm2 = dsqrt(j(2,1)*j(2,1)
     &                + j(2,2)*j(2,2)
     &                + j(2,3)*j(2,3))

         jnorm3 = dsqrt(j(3,1)*j(3,1)
     &                + j(3,2)*j(3,2)
     &                + j(3,3)*j(3,3))

         jscale = jnorm1*jnorm2*jnorm3

         jac_tol = 1.0d-12

         if (jscale .le. 0.0d0) then

            write(*,*) 'ERROR in gauss_point_hex8DR:'
            write(*,*) 'Degenerate Jacobian.'
            write(*,*) 'xi,eta,zeta = ', xi, eta, zeta
            write(*,*) 'detJ        = ', detj
            stop

         endif

         jac_quality = dabs(detj)/jscale


c----- Negative det(J) means inverted element

         if (detj .le. 0.0d0) then

            write(*,*) 'ERROR in gauss_point_hex8DR:'
            write(*,*) 'Negative or zero Jacobian determinant.'
            write(*,*) 'Element is inverted or degenerate.'
            write(*,*) 'xi,eta,zeta = ', xi, eta, zeta
            write(*,*) 'detJ        = ', detj
            write(*,*) 'J quality   = ', jac_quality
            stop

         endif


c----- Nearly singular Jacobian

         if (jac_quality .lt. jac_tol) then

            write(*,*) 'ERROR in gauss_point_hex8DR:'
            write(*,*) 'Nearly singular Jacobian.'
            write(*,*) 'xi,eta,zeta = ', xi, eta, zeta
            write(*,*) 'detJ        = ', detj
            write(*,*) 'J quality   = ', jac_quality
            stop

         endif


c=======================================================================
c  Inverse Jacobian
c=======================================================================

         detji = 1.0d0/detj

         invj(1,1) =
     &      detji*(j(2,2)*j(3,3) - j(2,3)*j(3,2))

         invj(1,2) =
     &      detji*(j(1,3)*j(3,2) - j(1,2)*j(3,3))

         invj(1,3) =
     &      detji*(j(1,2)*j(2,3) - j(1,3)*j(2,2))


         invj(2,1) =
     &      detji*(j(2,3)*j(3,1) - j(2,1)*j(3,3))

         invj(2,2) =
     &      detji*(j(1,1)*j(3,3) - j(1,3)*j(3,1))

         invj(2,3) =
     &      detji*(j(1,3)*j(2,1) - j(1,1)*j(2,3))


         invj(3,1) =
     &      detji*(j(2,1)*j(3,2) - j(2,2)*j(3,1))

         invj(3,2) =
     &      detji*(j(1,2)*j(3,1) - j(1,1)*j(3,2))

         invj(3,3) =
     &      detji*(j(1,1)*j(2,2) - j(1,2)*j(2,1))


c=======================================================================
c  Physical gradients
c
c       grad_x(N) = J^{-T} grad_xi(N)
c=======================================================================

         do i = 1, 8

            dnx1(i) =
     &         invj(1,1)*dn_dxi(1,i)
     &       + invj(2,1)*dn_dxi(2,i)
     &       + invj(3,1)*dn_dxi(3,i)

            dny1(i) =
     &         invj(1,2)*dn_dxi(1,i)
     &       + invj(2,2)*dn_dxi(2,i)
     &       + invj(3,2)*dn_dxi(3,i)

            dnz1(i) =
     &         invj(1,3)*dn_dxi(1,i)
     &       + invj(2,3)*dn_dxi(2,i)
     &       + invj(3,3)*dn_dxi(3,i)

         enddo


c=======================================================================
c  Unit velocity vector
c=======================================================================

         eps_v = 1.0d-14
         eps_g = 1.0d-14

         vmag = dsqrt(v(1)*v(1)
     &              + v(2)*v(2)
     &              + v(3)*v(3))

         if (vmag .gt. eps_v) then

            ex = v(1)/vmag
            ey = v(2)/vmag
            ez = v(3)/vmag

         else

            ex = 0.0d0
            ey = 0.0d0
            ez = 0.0d0

         endif


c=======================================================================
c  Stream direction expressed in reference coordinates
c
c       g = J^{-1} e
c=======================================================================

         gx = invj(1,1)*ex
     &      + invj(1,2)*ey
     &      + invj(1,3)*ez

         gy = invj(2,1)*ex
     &      + invj(2,2)*ey
     &      + invj(2,3)*ez

         gz = invj(3,1)*ex
     &      + invj(3,2)*ey
     &      + invj(3,3)*ez

         gnorm = dsqrt(gx*gx + gy*gy + gz*gz)


c=======================================================================
c  Directional element length
c=======================================================================

         if (vmag .gt. eps_v .and. gnorm .gt. eps_g) then

            h_stream = 2.0d0/gnorm

         else

            h_stream = 0.0d0

         endif


c=======================================================================
c  Effective diffusion along the streamline
c
c       lambda_s = e^T Lambda e
c
c  Lambda is diagonal here.
c=======================================================================

         lambda_s =
     &        lambda(1)*ex*ex
     &      + lambda(2)*ey*ey
     &      + lambda(3)*ez*ez


c----- Negative physical diffusion is invalid

         if (lambda_s .lt. 0.0d0) then

            write(*,*) 'ERROR in gauss_point_hex8DR:'
            write(*,*) 'Negative streamwise diffusion.'
            write(*,*) 'lambda_s = ', lambda_s
            stop

         endif


c=======================================================================
c  tau_DR
c=======================================================================

         tau_DR = 0.0d0
         pe     = 0.0d0
         betaPe = 0.0d0

         if (vmag .gt. eps_v .and.
     &       h_stream .gt. 0.0d0) then


c-------- advscale has units of diffusion coefficient:
c
c             advscale = |v| h_stream / 2
c
c        therefore Pe = advscale/lambda_s

            advscale = 0.5d0*vmag*h_stream


c=======================================================================
c  PURE ADVECTION
c
c       lambda_s = 0
c       Pe -> infinity
c       beta(Pe) -> 1
c
c       tau_DR = h_stream/(2|v|)
c=======================================================================

            if (lambda_s .eq. 0.0d0) then

               betaPe = 1.0d0

               tau_DR =
     &            0.5d0*h_stream/vmag


c=======================================================================
c  Very large Peclet number
c
c  Avoid forming a potentially huge Pe explicitly.
c
c       if Pe > 50:
c
c       beta = 1 - 1/Pe
c            = 1 - lambda_s/advscale
c=======================================================================

            else if (lambda_s .lt. advscale/50.0d0) then

               betaPe =
     &            1.0d0 - lambda_s/advscale

               tau_DR =
     &            0.5d0*h_stream/vmag*betaPe


c=======================================================================
c  Finite Peclet number
c=======================================================================

            else

               pe = advscale/lambda_s


c----------- Small Pe: series expansion

               if (pe .lt. 1.0d-3) then

                  betaPe =
     &                 pe/3.0d0
     &               - pe**3/45.0d0
     &               + 2.0d0*pe**5/945.0d0


c----------- Moderate Pe

               else

                  betaPe =
     &               1.0d0/dtanh(pe)
     &               - 1.0d0/pe

               endif


               tau_DR =
     &            0.5d0*h_stream/vmag*betaPe

            endif

         endif


c=======================================================================
c  Directional Residual test function
c
c       W_i = N_i + tau_DR (v . grad N_i)
c=======================================================================

         do i = 1, 8

            vdotg =
     &           v(1)*dnx1(i)
     &         + v(2)*dny1(i)
     &         + v(3)*dnz1(i)

            w(i) =
     &         n(i) + tau_DR*vdotg

         enddo


c=======================================================================
c  Physical diffusion remains Galerkin.
c
c  The diffusion term uses grad(N).
c  The DR perturbation of the test function is NOT differentiated
c  in the physical diffusion term.
c=======================================================================

         do i = 1, 8

            dnxw(i) = dnx1(i)
            dnyw(i) = dny1(i)
            dnzw(i) = dnz1(i)

         enddo


         return
      end