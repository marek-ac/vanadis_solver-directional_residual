          program analiza_1
          external mouse_trap

          integer (kind=1) irgb(3,256)
          dimension x1(322),x2(322),s(322),i_t(286,4),stezenia(511,157)
          dimension wsp_w(2,4)
          dimension f_n(4),s_w(4)
          character*12 dane
          character*12 dane_x,dane_y,dane_s
          integer (kind=3) q
          integer (kind=2) cursor_h,cursor_v,button_status
          common cursor_h,cursor_v,button_status,dane_x,dane_y,dane_s,
     &     stezenia,buffer

          nx1=6
          nx2=13
          nx3=22

          non=(nx1+1)*(nx2+1)*(nx3+1)
          non2=(nx2+1)*(nx3+1)

  1     format(a12)
  3     format (e10.4,2x,e10.4,2x,e10.4)
        call soua@('Podaj nazwe zbioru wejsciowego ZBIOR --->   ')
        read(*,1)dane
        open (9,FILE=dane,STATUS='OLD')
            do i=1,non2
            read(9,3)x1(i),x2(i),s(i)
            enddo
        close(9)

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

         if (stezenia(i-9,j-9).le.0.001) ij=23
         if (stezenia(i-9,j-9).ge.0.002) ij=24
         if (stezenia(i-9,j-9).ge.0.003) ij=25
         if (stezenia(i-9,j-9).ge.0.004) ij=26
         if (stezenia(i-9,j-9).ge.0.005) ij=27
         if (stezenia(i-9,j-9).ge.0.006) ij=28
         if (stezenia(i-9,j-9).ge.0.007) ij=29
         if (stezenia(i-9,j-9).ge.0.008) ij=30
         if (stezenia(i-9,j-9).ge.0.010) ij=31
         if (stezenia(i-9,j-9).ge.0.012) ij=33
         if (stezenia(i-9,j-9).ge.0.014) ij=34
         if (stezenia(i-9,j-9).ge.0.016) ij=35
         if (stezenia(i-9,j-9).ge.0.018) ij=36
         if (stezenia(i-9,j-9).ge.0.020) ij=37

         call set_pixel@(i+1,j+1,ij)

         enddo
         enddo

        i=10
        j=200
        is=50
        is2=5
        ig3=-38
        ig=2
        call fill_rectangle@(i,j       ,i+is,j+10*ig,37)
        call fill_rectangle@(i,j+10*ig ,i+is,j+20*ig,36)
        call fill_rectangle@(i,j+20*ig ,i+is,j+30*ig,35)
        call fill_rectangle@(i,j+30*ig ,i+is,j+40*ig,34)
        call fill_rectangle@(i,j+40*ig ,i+is,j+50*ig,33)
        call fill_rectangle@(i,j+50*ig ,i+is,j+60*ig,31)
        call fill_rectangle@(i,j+60*ig ,i+is,j+70*ig,30)
        call fill_rectangle@(i,j+70*ig ,i+is,j+80*ig,29)
        call fill_rectangle@(i,j+80*ig ,i+is,j+90*ig,28)
        call fill_rectangle@(i,j+90*ig ,i+is,j+100*ig,27)
        call fill_rectangle@(i,j+100*ig,i+is,j+110*ig,26)
        call fill_rectangle@(i,j+110*ig,i+is,j+120*ig,25)
        call fill_rectangle@(i,j+120*ig,i+is,j+130*ig,24)
        call fill_rectangle@(i,j+130*ig,i+is,j+140*ig,23)


        call draw_text@('mg/m^3',i+is2-4  ,j+20*ig+ig3-18  ,38)

        call draw_text@('0.020',i+is2  ,j+20*ig+ig3  ,32)
        call draw_text@('0.018',i+is2  ,j+30*ig+ig3  ,32)
        call draw_text@('0.016',i+is2  ,j+40*ig+ig3  ,32)
        call draw_text@('0.014',i+is2  ,j+50*ig+ig3  ,32)
        call draw_text@('0.012',i+is2  ,j+60*ig+ig3  ,32)
        call draw_text@('0.010',i+is2  ,j+70*ig+ig3  ,32)
        call draw_text@('0.008',i+is2  ,j+80*ig+ig3  ,32)
        call draw_text@('0.007',i+is2  ,j+90*ig+ig3  ,32)
        call draw_text@('0.006',i+is2  ,j+100*ig+ig3 ,32)
        call draw_text@('0.005',i+is2  ,j+110*ig+ig3 ,32)
        call draw_text@('0.004',i+is2  ,j+120*ig+ig3  ,32)
        call draw_text@('0.003',i+is2  ,j+130*ig+ig3  ,32)
        call draw_text@('0.002',i+is2  ,j+140*ig+ig3  ,32)
        call draw_text@('0.000',i+is2  ,j+150*ig+ig3  ,32)

      call beep@

c mouse
        call initialise_mouse@
        call set_trap@(mouse_trap,q,4)
        call set_mouse_interrupt_mask@(11)
        call set_mouse_bounds@(10,10,10+510,10+156)
        call display_mouse_cursor@

c         call get_key@(k)

         call draw_text@('X=',570,11,7)
         call draw_text@('Y=',570,31,7)
         call draw_text@('S=',570,51,7)

         do
         enddo


        call get_key@(k)
        call return_storage@(buffer)

        call text_mode@
       stop '****** OK'


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

        subroutine mouse_trap
        integer (kind=2) cursor_h,cursor_v,button_status
        character*12 dane_x,dane_y,dane_s
        dimension stezenia(511,157)

        common cursor_h,cursor_v,button_status,dane_x,dane_y,dane_s,
     &   stezenia,buffer

        call get_mouse_position@(cursor_h,cursor_v,button_status)

         if (button_status.eq.1) then

         ii=cursor_h
         ii2=cursor_v

          y=-650+1300./156.*(ii2-1-9)
          x=-250+4250./510.*(ii-1-9)

         m=iabs(int(x))
         K=iabs(INT(x/10))
         L=iabs(INT(x/100))
         n=iabs(INT(x/1000))
         if (x.ge.0.) dane_x=' '//CHAR(48+n)//CHAR(48+L-10*n)
     &   //CHAR(48+K-10*L)//CHAR(48+m-10*K)

         if (x.lt.0.) dane_x='-'//CHAR(48+n)//CHAR(48+L-10*n)
     &   //CHAR(48+K-10*L)//CHAR(48+m-10*K)

         m=iabs(int(y))
         K=iabs(INT(y/10))
         L=iabs(INT(y/100))
         n=iabs(INT(y/1000))
         if (y.ge.0.) dane_y=' '//CHAR(48+n)//CHAR(48+L-10*n)
     &   //CHAR(48+K-10*L)//CHAR(48+m-10*K)
         if (y.lt.0.) dane_y='-'//CHAR(48+n)//CHAR(48+L-10*n)
     &   //CHAR(48+K-10*L)//CHAR(48+m-10*K)

         si=(stezenia(ii-9,ii2-9))*1e5

         m=iabs(int(si))
         K=iabs(INT(si/10))
         L=iabs(INT(si/100))
         n=iabs(INT(si/1000))
         dane_s=' '//CHAR(48+n)//CHAR(48+L-10*n)
     &   //CHAR(48+K-10*L)//CHAR(48+m-10*K)


         call fill_rectangle@(600,0,639,70,0)

         call draw_text@(dane_x,600,10,7)
         call draw_text@(dane_y,600,30,7)
         call draw_text@(dane_s,600,50,7)


         endif

        if (button_status.eq.2) then

        call text_mode@
        call beep@
        stop '****** OK'
        endif

        return
        end
