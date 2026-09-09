          program gaaf_cal1
          PARAMETER (Nx1=6,Nx2=13,nx3=22,iczas=1980)

          integer (kind=1) irgb(3,256)
          dimension x1(322),x2(322),s(322),i_t(286,4),stezenia(511,157)
          dimension wsp_w(2,4)
          dimension f_n(4),s_w(4)
          character*12 dane_x(nx1+1,iczas),dane_pcx(iczas)


          non=(nx1+1)*(nx2+1)*(nx3+1)
          non2=(nx2+1)*(nx3+1)

        do j=1,iczas
        Kt=INT(j/10)
        Lt=INT(j/100)
        Mt=INT(j/1000)
        do i=1,nx1+1
        K=INT(I/10)
        L=INT(I/100)
        dane_x(i,j)='X'//'_'//CHAR(48+Mt)//CHAR(48+Lt-10*Mt)
     &   //CHAR(48+Kt-10*Lt)//CHAR(48+j-10*Kt)//'.'//CHAR(48+L)//
     &   CHAR(48+K-10*L)//CHAR(48+I-10*K)
        enddo
        enddo

        do jt=1,iczas
        Kt=INT(jt/10)
        Lt=INT(jt/100)
        Mt=INT(jt/1000)
        dane_pcx(jt)='X1_'//CHAR(48+Mt)//CHAR(48+Lt-10*Mt)//
     &   CHAR(48+Kt-10*Lt)//CHAR(48+jt-10*Kt)//'.pcx'
        enddo


            j=0
            k=1
         do i=1,nx2*nx3
            j=j+1
            i_t(j,1)=k
            i_t(j,4)=k+nx2+1
            i_t(j,3)=k+nx2+2
            i_t(j,2)=k+1
            k=k+1
            if (mod(k,(nx2+1)).eq.0) k=k+1
         enddo




         call vga@

         do j=1,256
         do i=1,3

         ir=42-j*2+5
         ig=21-j*2+5
         ib=0

         if (j.eq.8) then
         ir=0
         ig=50
         ib=53
         endif

         if (i.eq.1) irgb(i,j)=ir
         if (i.eq.2) irgb(i,j)=ig
         if (i.eq.3) irgb(i,j)=ib

         enddo
         enddo

         call set_video_dac_block@(0,256,irgb)
         call set_video_dac@(0,0,0,0)



        do i_iczas=5,iczas,5

        open (9,FILE=dane_x(1,i_iczas),STATUS='OLD')
            do i=1,non2
            read(9,3)x1(i),x2(i),s(i)
            enddo
        close(9)
3       format (e10.4,2x,e10.4,2x,e10.4)


         do j=1+9,157+9
         do i=1+9,511+9


c       wspolrzedne piksela
        x=-650+1300./156.*(j-1-9)
        y=-250+4250./510.*(i-1-9)

        do iii=1,nx2*nx3
        if ( x.ge.x1(i_t(iii,1)).and.x.lt.x1(i_t(iii,2)).and.
     &   y.ge.x2(i_t(iii,2)).and.y.lt.x2(i_t(iii,3)) ) goto 4
        enddo
4       ii=iii

        wsp_w(1,1)=x1(i_t(ii,1))
        wsp_w(1,2)=x1(i_t(ii,2))
        wsp_w(1,3)=x1(i_t(ii,3))
        wsp_w(1,4)=x1(i_t(ii,4))
        wsp_w(2,1)=x2(i_t(ii,1))
        wsp_w(2,2)=x2(i_t(ii,2))
        wsp_w(2,3)=x2(i_t(ii,3))
        wsp_w(2,4)=x2(i_t(ii,4))

        s_w(1)=s(i_t(ii,1))
        s_w(2)=s(i_t(ii,2))
        s_w(3)=s(i_t(ii,3))
        s_w(4)=s(i_t(ii,4))

        call wsp_ksi_eta(wsp_w,x,y,p_ksi,p_eta)
        call F_k(f_n,p_ksi,p_eta)
        call stezenie(f_n,s_w,ss)
        stezenia(i-9,j-9)=ss

         if (stezenia(i-9,j-9).lt.0.001) ij=23
         if (stezenia(i-9,j-9).ge.0.001) ij=24
         if (stezenia(i-9,j-9).ge.0.003) ij=25
         if (stezenia(i-9,j-9).ge.0.004) ij=26
         if (stezenia(i-9,j-9).ge.0.005) ij=27
         if (stezenia(i-9,j-9).ge.0.006) ij=28
         if (stezenia(i-9,j-9).ge.0.007) ij=29
         if (stezenia(i-9,j-9).ge.0.008) ij=30
         if (stezenia(i-9,j-9).ge.0.010) ij=31
         if (stezenia(i-9,j-9).ge.0.012) ij=33
         if (stezenia(i-9,j-9).ge.0.016) ij=34
         if (stezenia(i-9,j-9).ge.0.020) ij=35
         if (stezenia(i-9,j-9).ge.0.024) ij=36
         if (stezenia(i-9,j-9).ge.0.028) ij=37

         call set_pixel@(i+1,j+1,ij)

         enddo
         enddo

      call beep@

      call get_screen_block@(10,10,520,166,buffer)
      if(buffer.eq.-1) then
      error_code=1
      goto 10
      endif


      call screen_block_to_pcx@(dane_pcx(i_iczas),buffer,error_code)
      if (error_code.lt.0) then
      erro_code=1
      goto 10
      else
      call doserr@(error_code)
      endif

      call clear_screen@
      call return_storage@(buffer)

      enddo


      call text_mode@
      stop '****** OK'

 10   call text_mode@
      if(error_code.eq.1) stop'error: out of memory'

      end

        Subroutine F_k(f_n,p_ksi,p_eta)

        dimension f_n(4)

        f_n(1)=0.25*(1-p_ksi)*(1-p_eta)
        f_n(2)=0.25*(1+p_ksi)*(1-p_eta)
        f_n(3)=0.25*(1+p_ksi)*(1+p_eta)
        f_n(4)=0.25*(1-p_ksi)*(1+p_eta)

        return
        end

        Subroutine wsp_xy(wsp_w,f_n,x,y)

        dimension f_n(4),wsp_w(2,4)
        x=0
        y=0
        do i=1,4
        x=x+f_n(i)*wsp_w(1,i)
        y=y+f_n(i)*wsp_w(2,i)
        enddo

        return
        end

        Subroutine wsp_ksi_eta(wsp_w,x,y,p_ksi,p_eta)

        dimension wsp_w(2,4)

        Sx1=wsp_w(1,1)+wsp_w(1,2)+wsp_w(1,3)+wsp_w(1,4)
        Sx3=-wsp_w(1,1)+wsp_w(1,2)+wsp_w(1,3)-wsp_w(1,4)
        Sy1=wsp_w(2,1)+wsp_w(2,2)+wsp_w(2,3)+wsp_w(2,4)
        Sy2=-wsp_w(2,1)-wsp_w(2,2)+wsp_w(2,3)+wsp_w(2,4)

        p_ksi=(4*x-Sx1)/Sx3
        p_eta=(4*y-Sy1)/Sy2


        return
        end

        Subroutine stezenie(f_n,s_w,s)

        dimension f_n(4),s_w(4)
        s=0
        do i=1,4
        s=f_n(i)*s_w(i)+s
        enddo

        return
        end

