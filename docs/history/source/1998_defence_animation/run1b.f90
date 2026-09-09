          program odczyt_pcx
          PARAMETER (Nx1=6,Nx2=13,nx3=22,iczas=1980)
          integer (kind=2) error_code(iczas)
          integer (kind=3) buffer(iczas)
          integer (kind=1) irgb(3,256)
          character*12 dane_pcx(iczas),czas(iczas)
          character*17 palette

        do jt=1,iczas
        Kt=INT(jt/10)
        Lt=INT(jt/100)
        Mt=INT(jt/1000)
        dane_pcx(jt)='X1_'//CHAR(48+Mt)//CHAR(48+Lt-10*Mt)//
     &   CHAR(48+Kt-10*Lt)//CHAR(48+jt-10*Kt)//'.pcx'
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


      call new_page@

      do i_czas=5,iczas,5
      call pcx_to_screen_block@(dane_pcx(i_czas),buffer(i_czas)
     &,palette,error_code(i_czas))
      call doserr@(error_code(i_czas))
      enddo


        call fill_rectangle@(10,10,520,166,23)

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

        call draw_text@('0.028',i+is2  ,j+20*ig+ig3  ,32)
        call draw_text@('0.024',i+is2  ,j+30*ig+ig3  ,32)
        call draw_text@('0.020',i+is2  ,j+40*ig+ig3  ,32)
        call draw_text@('0.016',i+is2  ,j+50*ig+ig3  ,32)
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

c opis osi

        call draw_text@('-650',510+10+5,10+5,23)
        call draw_text@(' 650',510+10+5,156+10-15,23)
        call draw_text@('-250',10+5,166+5,23)
        call draw_text@('4000',510+10-35,166+5,23)



        call set_text_attribute@(102,4.,0.,0.)

        call draw_text@('czas: 0000s',220,250,23)


        call get_key@(k)


      do i_czas=5,iczas,5
      call restore_screen_block@(10,10,buffer(i_czas),0,
     & error_code(i_czas))

      if(buffer(i_czas).eq.-1) then
      error_code(i_czas)=1
      goto 10
      endif


        jt=i_czas
        Kt=INT(jt/10)
        Lt=INT(jt/100)
        Mt=INT(jt/1000)
        czas(jt)='czas: '//CHAR(48+Mt)//CHAR(48+Lt-10*Mt)//
     &   CHAR(48+Kt-10*Lt)//CHAR(48+jt-10*Kt)//'s'

        call fill_rectangle@(220,210,560,250,0)

        call draw_text@(czas(i_czas),220,250,23)


      call return_storage@(buffer(i_czas))

      call sleep@(0.1)
      enddo

      call get_key@(k)
      call text_mode@


10    call text_mode@
      if(error_code(i_czas).eq.1) stop'*******error: out of memory'

      end



