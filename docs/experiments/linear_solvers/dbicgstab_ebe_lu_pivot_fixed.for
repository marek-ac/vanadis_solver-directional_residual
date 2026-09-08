      subroutine lu_pp8_bicg(a, lu, ipiv, nreg, izero)
      implicit real*8 (a-h,o-z)

      real*8 a(8,8), lu(8,8)
      real*8 anorm, pmax, pivmin, tmp
      integer ipiv(8), nreg, izero
      integer i, j, k, p

c====================================================================
c     In-place LU factorization with partial pivoting for one HEX8
c     element matrix.  The result satisfies approximately
c
c         P*A = L*U
c
c     in the packed array LU.  ipiv(k) stores the row interchanged
c     with row k at elimination step k.
c====================================================================

      anorm = 0.0d0
      do 20 i = 1, 8
         do 10 j = 1, 8
            lu(i,j) = a(i,j)
            if (dabs(a(i,j)) .gt. anorm) anorm = dabs(a(i,j))
   10    continue
         ipiv(i) = i
   20 continue

      nreg  = 0
      izero = 0

c----- Empty local block: use identity as a harmless local inverse.
      if (anorm .eq. 0.0d0) then
         izero = 1
         do 40 i = 1, 8
            do 30 j = 1, 8
               lu(i,j) = 0.0d0
   30       continue
            lu(i,i) = 1.0d0
            ipiv(i) = i
   40    continue
         return
      endif

c----- Relative pivot floor.  It is used only after pivoting failed
c----- to find a sufficiently large pivot in the current column.
      pivmin = 1.0d-14 * anorm
      if (pivmin .lt. 1.0d-300) pivmin = 1.0d-300

      do 100 k = 1, 8

c-------- Partial pivoting.
         p = k
         pmax = dabs(lu(k,k))
         do 50 i = k+1, 8
            if (dabs(lu(i,k)) .gt. pmax) then
               pmax = dabs(lu(i,k))
               p = i
            endif
   50    continue

         ipiv(k) = p

         if (p .ne. k) then
            do 60 j = 1, 8
               tmp     = lu(k,j)
               lu(k,j) = lu(p,j)
               lu(p,j) = tmp
   60       continue
         endif

c-------- Last-resort regularization of a tiny pivot.
         if (dabs(lu(k,k)) .lt. pivmin) then
            if (lu(k,k) .lt. 0.0d0) then
               lu(k,k) = -pivmin
            else
               lu(k,k) =  pivmin
            endif
            nreg = nreg + 1
         endif

c-------- Gaussian elimination, packed L below diagonal and U on/above.
         if (k .lt. 8) then
            do 80 i = k+1, 8
               lu(i,k) = lu(i,k) / lu(k,k)
               do 70 j = k+1, 8
                  lu(i,j) = lu(i,j) - lu(i,k)*lu(k,j)
   70          continue
   80       continue
         endif

  100 continue

      return
      end


      subroutine element_lu_all_pp_bicg(bb, nelem, lu_el,
     &                                  ipiv_el)
      implicit real*8 (a-h,o-z)

      integer nelem, ipiv_el(8,nelem)
      real*8 bb(*), lu_el(8,8,nelem)
      real*8 ke(8,8)
      integer e, i, j, ik, nreg, izero
      integer nreg_total, nzero_total

      nreg_total  = 0
      nzero_total = 0

      do 100 e = 1, nelem
         ik = 64*(e-1)
         do 20 i = 1, 8
            do 10 j = 1, 8
               ke(i,j) = bb(ik + (i-1)*8 + j)
   10       continue
   20    continue

         call lu_pp8_bicg(ke, lu_el(1,1,e), ipiv_el(1,e),
     &                    nreg, izero)

         nreg_total  = nreg_total  + nreg
         nzero_total = nzero_total + izero
  100 continue

      if (nreg_total .gt. 0) then
         write(*,*) 'BiCGSTAB LU: regularized pivots = ', nreg_total
      endif
      if (nzero_total .gt. 0) then
         write(*,*) 'BiCGSTAB LU: empty element blocks = ', nzero_total
      endif

      return
      end


      subroutine msolve_el_pp_bicg(n, r, z, iss, nelem,
     &                             lu_el, ipiv_el)
      implicit real*8 (a-h,o-z)

      integer n, nelem, iss(nelem,8), ipiv_el(8,nelem)
      real*8 r(n), z(n), lu_el(8,8,nelem)
      real*8 rhs(8), y(8), zloc(8), sum, tmp
      integer e, a, i, j, g, p

c====================================================================
c     Additive EBE preconditioner:
c
c       z = sum_e R_e^T * K_e^{-1} * R_e * r
c
c     Local solves use the pivoted LU factors prepared once before
c     the Krylov iteration.
c====================================================================

      do 10 i = 1, n
         z(i) = 0.0d0
   10 continue

      do 100 e = 1, nelem

c-------- Gather local residual.
         do 20 a = 1, 8
            g = iss(e,a)
            rhs(a)  = r(g)
            y(a)    = 0.0d0
            zloc(a) = 0.0d0
   20    continue

c-------- Apply the same row permutations used during factorization.
         do 30 i = 1, 8
            p = ipiv_el(i,e)
            if (p .ne. i) then
               tmp    = rhs(i)
               rhs(i) = rhs(p)
               rhs(p) = tmp
            endif
   30    continue

c-------- Forward substitution: L*y = P*r_e, diag(L)=1.
         do 50 i = 1, 8
            sum = 0.0d0
            do 40 j = 1, i-1
               sum = sum + lu_el(i,j,e)*y(j)
   40       continue
            y(i) = rhs(i) - sum
   50    continue

c-------- Back substitution: U*zloc = y.
         do 70 i = 8, 1, -1
            sum = 0.0d0
            do 60 j = i+1, 8
               sum = sum + lu_el(i,j,e)*zloc(j)
   60       continue
            zloc(i) = (y(i)-sum) / lu_el(i,i,e)
   70    continue

c-------- Additive scatter to global vector.
         do 80 a = 1, 8
            g = iss(e,a)
            z(g) = z(g) + zloc(a)
   80    continue

  100 continue

      return
      end


      subroutine dbicgstab_full(n, b, x, tol, itmax, max_restart,
     &       iter, err, bb, iss, nax, nay, naz, ilosc_el)

      implicit real*8 (a-h,o-z)

      integer n, itmax, max_restart, iter
      integer nax, nay, naz, ilosc_el
      real*8 tol, err
      real*8 b(n), x(n), bb(*)
      integer iss((nax-1)*(nay-1)*(naz-1),8)

      real*8, allocatable :: r(:), r0(:), p(:), v(:)
      real*8, allocatable :: phat(:), shat(:), tt(:)
      real*8, allocatable :: lu_el(:,:,:)
      integer, allocatable :: ipiv_el(:,:)

      real*8 rho, rho_old, alpha, omega, beta
      real*8 bnrm, resid, denom, numer
      real*8 nr0, nr, nv, scale, eps_break, tiny_abs
      real*8 last_resid, growth_limit
      integer stagnation_count
      integer i, k, restart

      allocate(r(n), r0(n), p(n), v(n))
      allocate(phat(n), shat(n), tt(n))
      allocate(lu_el(8,8,ilosc_el), ipiv_el(8,ilosc_el))

c====================================================================
c     PRECONDITIONER: PIVOTED ELEMENT LU (EBE)
c====================================================================
      call element_lu_all_pp_bicg(bb, ilosc_el, lu_el, ipiv_el)

      eps_break = 1.0d-14
      tiny_abs  = 1.0d-300

      bnrm = dnrm2(n, b)
      if (bnrm .eq. 0.0d0) bnrm = 1.0d0

      iter = 0
      err  = 1.0d0

c====================================================================
c     RESTART LOOP
c
c     Restarts are used only as breakdown/stagnation recovery.
c     At every restart the true residual b-A*x is recomputed and used
c     as a new shadow residual r0.
c====================================================================
      do 1000 restart = 0, max_restart

c-------- True residual r = b - A*x.
         call matvec(x, r, iss, bb, nax, nay, naz)
         do 10 i = 1, n
            r(i)  = b(i) - r(i)
            r0(i) = r(i)
            p(i)  = 0.0d0
            v(i)  = 0.0d0
   10    continue

         resid = dnrm2(n, r) / bnrm
         err   = resid
         if (resid .le. tol) goto 9000

         rho_old = 1.0d0
         alpha   = 1.0d0
         omega   = 1.0d0

         nr0 = dnrm2(n, r0)
         last_resid = resid
         growth_limit = 2.0d0
         stagnation_count = 0

c====================================================================
c     MAIN RIGHT-PRECONDITIONED BICGSTAB LOOP
c====================================================================
         do 500 k = 1, itmax

            iter = iter + 1

c----------- rho = (r0,r)
            rho = 0.0d0
            do 20 i = 1, n
               rho = rho + r0(i)*r(i)
   20       continue

            nr = dnrm2(n, r)
            scale = nr0*nr
            if (scale .lt. tiny_abs) scale = tiny_abs

            if (dabs(rho) .le. eps_break*scale) then
               write(*,*) 'Restart: rho breakdown'
               goto 600
            endif

c----------- Search direction.  At every restart k=1 must start p=r.
            if (k .eq. 1) then
               do 30 i = 1, n
                  p(i) = r(i)
   30          continue
            else
               if (dabs(omega) .le. tiny_abs) then
                  write(*,*) 'Restart: omega breakdown'
                  goto 600
               endif

               beta = (rho/rho_old)*(alpha/omega)
               do 40 i = 1, n
                  p(i) = r(i) + beta*(p(i)-omega*v(i))
   40          continue
            endif

c----------- phat = M^{-1} p
            call msolve_el_pp_bicg(n, p, phat, iss, ilosc_el,
     &                             lu_el, ipiv_el)

c----------- v = A*phat.  IMPORTANT: v must not be overwritten by
c----------- the second matrix-vector product later in the iteration.
            call matvec(phat, v, iss, bb, nax, nay, naz)

c----------- alpha = rho / (r0,v)
            denom = 0.0d0
            do 50 i = 1, n
               denom = denom + r0(i)*v(i)
   50       continue

            nv = dnrm2(n, v)
            scale = nr0*nv
            if (scale .lt. tiny_abs) scale = tiny_abs

            if (dabs(denom) .le. eps_break*scale) then
               write(*,*) 'Restart: alpha breakdown'
               goto 600
            endif

            alpha = rho / denom

c----------- s = r - alpha*v.  Reuse r as the BiCGSTAB s vector.
            do 60 i = 1, n
               r(i) = r(i) - alpha*v(i)
   60       continue

            resid = dnrm2(n, r) / bnrm
            err   = resid

c----------- Early convergence after alpha step.  The old routine
c----------- returned here without updating x.  That is incorrect.
            if (resid .le. tol) then
               do 70 i = 1, n
                  x(i) = x(i) + alpha*phat(i)
   70          continue

c-------------- Verify with the true residual b-A*x.
               call matvec(x, tt, iss, bb, nax, nay, naz)
               do 80 i = 1, n
                  r(i) = b(i) - tt(i)
   80          continue
               err = dnrm2(n, r) / bnrm
               if (err .le. tol) goto 9000

               write(*,*) 'Restart: recursive residual drift'
               goto 600
            endif

c----------- shat = M^{-1} s
            call msolve_el_pp_bicg(n, r, shat, iss, ilosc_el,
     &                             lu_el, ipiv_el)

c----------- tt = A*shat.  This MUST be a separate vector from v.
            call matvec(shat, tt, iss, bb, nax, nay, naz)

c----------- omega = (tt,s)/(tt,tt)
            numer = 0.0d0
            denom = 0.0d0
            do 90 i = 1, n
               numer = numer + tt(i)*r(i)
               denom = denom + tt(i)*tt(i)
   90       continue

            if (denom .le. tiny_abs) then
               write(*,*) 'Restart: omega denominator breakdown'
               goto 600
            endif

            omega = numer / denom
            if (dabs(omega) .le. tiny_abs) then
               write(*,*) 'Restart: omega breakdown'
               goto 600
            endif

c----------- x = x + alpha*phat + omega*shat
c----------- r = s - omega*tt
            do 100 i = 1, n
               x(i) = x(i) + alpha*phat(i) + omega*shat(i)
               r(i) = r(i) - omega*tt(i)
  100       continue

            resid = dnrm2(n, r) / bnrm
            err   = resid

            write(*,'(A,I6,A,E12.5)') 'BiCGSTAB iter=',
     &           iter,'  resid=',resid

c----------- Never report convergence from the recursive residual only.
            if (resid .le. tol) then
               call matvec(x, tt, iss, bb, nax, nay, naz)
               do 110 i = 1, n
                  r(i) = b(i) - tt(i)
  110          continue
               err = dnrm2(n, r) / bnrm
               if (err .le. tol) goto 9000

               write(*,*) 'Restart: recursive residual drift'
               goto 600
            endif

c====================================================================
c     ADAPTIVE RESTART CONDITIONS
c====================================================================

c----------- Strong residual growth.
            if (resid .gt. growth_limit*last_resid) then
               write(*,*) 'Restart: residual blow-up'
               goto 600
            endif

c----------- Stagnation: less than about five percent improvement
c----------- for eight consecutive iterations.
            if (resid .gt. 0.95d0*last_resid) then
               stagnation_count = stagnation_count + 1
            else
               stagnation_count = 0
            endif

            if (stagnation_count .ge. 8) then
               write(*,*) 'Restart: stagnation'
               goto 600
            endif

            last_resid = resid
            rho_old = rho

  500    continue

c-------- itmax iterations completed in this restart cycle.
         write(*,*) 'Restart: iteration limit in current cycle'

  600    continue

c-------- If no further restart is allowed, return the TRUE residual.
         if (restart .eq. max_restart) then
            call matvec(x, tt, iss, bb, nax, nay, naz)
            do 610 i = 1, n
               r(i) = b(i) - tt(i)
  610       continue
            err = dnrm2(n, r) / bnrm
            write(*,*) 'BiCGSTAB: reached max_restart, resid=', err
            goto 9000
         endif

 1000 continue

 9000 continue
      deallocate(r, r0, p, v, phat, shat, tt)
      deallocate(lu_el, ipiv_el)

      return
      end
