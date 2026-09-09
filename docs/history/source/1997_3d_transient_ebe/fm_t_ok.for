c uklad stacjonarny gradienty sp. i eliminacja Gaussa + czas
c spr. dokl. jest poprawny

      program f_m_gr_sp
      parameter (nax=7,nay=14,naz=23)
      implicit real*8 (a-h,o-z)
      dimension dx(nax),dy(nay),dz(naz)
      dimension XX(nax*nay*naz),YY(nax*nay*naz),ZZ(nax*nay*naz)
      dimension t(6)
      dimension iSS((nax-1)*(nay-1)*(naz-1),8)
      dimension kod((nax-1)*(nay-1)*(naz-1),6)
      dimension bb((nax-1)*(nay-1)*(naz-1)*64)

c      common /duza_t/ hma2(650,213),fmv2(650)

      common /gs_ps/ fmv3(2254)
      common /d_time/ i_d_czas

           eps1B=1.0d-20
           znorma=1.0d-15
           IL_ITERACJI=1000
           i_d_czas=5


      call par_3d
      call wsp_na_o (nax,nay,naz,dx,dy,dz)
      call siatka   (XX,YY,ZZ,dx,dy,dz,nax,nay,naz)
      call w_el     (iSS,nax,nay,naz)
      call kod_el   (nax,nay,naz,kod)

      i_czas=0
1010  continue
      call przest3d(i_czas,XX,YY,ZZ,iSS,kod,t,nax,nay,naz,
     &              bb,EPS1B,ZNORMA,IL_ITERACJI)
      if (i_czas.le.2000) goto 1010


c      open (10,FILE='gauss.txt',STATUS='unknown')
c      do i=1,nax*nay*naz
c      write(10,125)xx(i),yy(i),zz(i),fmv2(i)*1e6
c      enddo
c      close(10)
c 125  format (e10.4,2x,e10.4,2x,e10.4,2x,e10.4)

      end

      subroutine wsp_na_o(nax,nay,naz,dx,dy,dz)
      implicit real*8 (a-h,o-z)
      dimension dx(*),dy(*),dz(*)
      common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband

      ilosc_el=(nax-1)*(nay-1)*(naz-1)
      ilosc_wezlow=nax*nay*naz
      mband=nax*nay+nax+2

      dx(1)=0.0d+00
      dx(2)=20.0d+00
      dx(3)=35.0d+00
      dx(4)=85.0d+00
      dx(5)=135.0d+00
      dx(6)=230.0d+00
      dx(7)=420.0d+00


      dy(1)=-650.0d+00
      dy(2)=-520.0d+00
      dy(3)=-390.0d+00
      dy(4)=-270.0d+00
      dy(5)=-180.0d+00
      dy(6)=-90.0d+00
      dy(7)=-25.0d+00
      dy(8)=25.0d+00
      dy(9)=90.0d+00
      dy(10)=180.0d+00
      dy(11)=270.0d+00
      dy(12)=390.0d+00
      dy(13)=520.0d+00
      dy(14)=650.0d+00

      dz(1)=-250.0d+00
      dz(2)=-100.0d+00
      dz(3)=-50.0d+00
      dz(4)=-25.0d+00
      dz(5)=25.0d+00
      dz(6)=80.0d+00
      dz(7)=150.0d+00
      dz(8)=250.0d+00
      dz(9)=450.0d+00
      dz(10)=600.0d+00
      dz(11)=750.0d+00
      dz(12)=900.0d+00
      dz(13)=1100.0d+00
      dz(14)=1300.0d+00
      dz(15)=1500.0d+00
      dz(16)=1800.0d+00
      dz(17)=2100.0d+00
      dz(18)=2400.0d+00
      dz(19)=2700.0d+00
      dz(20)=3000.0d+00
      dz(21)=3300.0d+00
      dz(22)=3600.0d+00
      dz(23)=4000.0d+00

      return
      end


      subroutine licz_alfa (i_czas,nr_el,kod,ALFA,qs,NAX,NAY,NAZ)
      implicit real*8 (a-h,o-z)
      dimension kod((nax-1)*(nay-1)*(naz-1),6)
c     qs=-alfa*Tf    Tf-temp.plynu omywajacego powierzchnie
      alfa=0.0d+00
      qs=0.0d+00
      return
      end


      subroutine atmosfera  (i_czas,nr_el,lambda,pzanik,v,t,qv,ro,cp)
      implicit real*8 (a-h,o-z)
      real*8 lambda
      dimension t(6),v(3),lambda(3)
      ro=1.0d+00
      cp=1.0d+00

      IF (I_CZAS.LE.1000) THEN
      ab_pasquil=1.43d+00
      lambda(1)=8.15d+00
      lambda(2)=lambda(1)*ab_pasquil
      lambda(3)=lambda(1)*ab_pasquil
      pzanik=0.000693147
      qv=0.0d+00
      if (nr_el.eq.273) qv=8.0d+00*1d-09
      v(1)=0.
      v(2)=0.
      v(3)=1.554128
      t(1)=0.
      t(2)=0.
      t(3)=0.
      t(4)=0.
      t(5)=0.
      t(6)=0.
      ENDIF

      IF (I_CZAS.GT.1000) THEN
      ab_pasquil=1.430d+00
      lambda(1)=8.150d+00
      lambda(2)=lambda(1)*ab_pasquil
      lambda(3)=lambda(1)*ab_pasquil
      pzanik=0.000693147
      qv=0.0d+00
      if (nr_el.eq.273) qv=8.0d+00*1d-09
C 1g/S PRZY 50*50*50

       if (nr_el.eq.687) qv=1.111*1d-09
C 1g/S PRZY 150*50*120

      v(1)=0.0d+00
      v(2)=-0.50d+00
      v(3)=1.554128
      t(1)=0.0d+00
      t(2)=0.0d+00
      t(3)=0.0d+00
      t(4)=0.0d+00
      t(5)=0.0d+00
      t(6)=0.0d+00
      ENDIF

c     wstrzymanie konwekcji przed przeszkoda
      if (nr_el.eq.661) v(3)=0.0d+00
      if (nr_el.eq.662) v(3)=0.0d+00
      if (nr_el.eq.663) v(3)=0.0d+00

      return
      end



      subroutine siatka (XX,YY,ZZ,dx,dy,dz,nax,nay,naz)
      implicit real*8 (a-h,o-z)
      dimension XX(*),YY(*),ZZ(*),dx(*),dy(*),dz(*)

          i=1
          do jz=1,naz
          do jy=1,nay
          do jx=1,nax
          XX(i)=dx(jx)
          YY(i)=dy(jy)
          ZZ(i)=dz(jz)
          i=i+1
          enddo
          enddo
          enddo

          return
          end

        subroutine w_el (iSS,nax,nay,naz)
        dimension iSS((nax-1)*(nay-1)*(naz-1),8)

        k=1
        ik=0
        ik2=nay*nax

        do j=1,(nax-1)*(nay-1)*(naz-1)

        iSS(j,1)=k                    +  ik2*ik
        iSS(j,2)=k+1                  +  ik2*ik
        iSS(j,3)=k+1    + nax         +  ik2*ik
        iSS(j,4)=k      + nax         +  ik2*ik
        iSS(j,5)=k      + nax*nay     +  ik2*ik
        iSS(j,6)=k+1    + nax*nay     +  ik2*ik
        iSS(j,7)=k+1    + nax*nay+nax +  ik2*ik
        iSS(j,8)=k      + nax*nay+nax +  ik2*ik

        k=k+1

        if( mod(j,(nax-1)).eq.0 ) k=k+1

        if( mod(j,(nax-1)*(nay-1)).eq.0 ) then
        k=1
        ik=ik+1
        endif

        enddo
        return
        end

        subroutine kod_el(nax,nay,naz,kod)
        dimension kod((nax-1)*(nay-1)*(naz-1),6)

        il_el=(nax-1)*(nay-1)*(naz-1)
        do i=1,il_el
        do j=1,6
        kod(i,j)=0
        enddo
        enddo

        do i=1,((nax-1)*(nay-1))
        kod(i,1)=1
        enddo

        do i=il_el-(nax-1)*(nay-1)+1,il_el
        kod(i,2)=2
        enddo

        do i=(nax-1),il_el,(nax-1)
        kod(i,4)=4
        enddo

        do i=1,il_el,nax-1
        kod(i,3)=3
        enddo

        do k=0,naz-2
        do i=k*(nax-1)*(nay-1)+1,k*(nax-1)*(nay-1)+nax-1
        kod(i,6)=6
        enddo
        enddo

        do k=1,naz-1
        do i=k*(nax-1)*(nay-1)-nax+1+1,k*(nax-1)*(nay-1)
        kod(i,5)=5
        enddo
        enddo

        return
        end


      subroutine par_3d
      implicit real*8 (a-h,o-z)
      real*8 ksi,eta,zeta
      real*8 KSI_G,ETA_G,ZETA_G

      common /LOK/ KSI(8),ETA(8),ZETA(8)
      common /LOK_G/ KSI_G(8),ETA_G(8),ZETA_G(8)
      common /POW_EL/ i_pow(6,4)

            KSI(1)=-1.0d+00
            KSI(2)=1.0d+00
            KSI(3)=1.0d+00
            KSI(4)=-1.0d+00
            KSI(5)=-1.0d+00
            KSI(6)=1.0d+00
            KSI(7)=1.0d+00
            KSI(8)=-1.0d+00

            ETA(1)=-1.0d+00
            ETA(2)=-1.0d+00
            ETA(3)=1.0d+00
            ETA(4)=1.0d+00
            ETA(5)=-1.0d+00
            ETA(6)=-1.0d+00
            ETA(7)=1.0d+00
            ETA(8)=1.0d+00

            ZETA(1)=-1.0d+00
            ZETA(2)=-1.0d+00
            ZETA(3)=-1.0d+00
            ZETA(4)=-1.0d+00
            ZETA(5)=1.0d+00
            ZETA(6)=1.0d+00
            ZETA(7)=1.0d+00
            ZETA(8)=1.0d+00

         do j=1,8
            KSI_G(j)  =0.577350269189626D0*KSI(j)
            ETA_G(j)  =0.577350269189626D0*ETA(j)
            ZETA_G(j) =0.577350269189626D0*ZETA(j)
         enddo

            i_pow(1,1)=4
            i_pow(1,2)=3
            i_pow(1,3)=2
            i_pow(1,4)=1

            i_pow(2,1)=5
            i_pow(2,2)=6
            i_pow(2,3)=7
            i_pow(2,4)=8

            i_pow(3,1)=4
            i_pow(3,2)=1
            i_pow(3,3)=5
            i_pow(3,4)=8

            i_pow(4,1)=2
            i_pow(4,2)=3
            i_pow(4,3)=7
            i_pow(4,4)=6

            i_pow(5,1)=4
            i_pow(5,2)=3
            i_pow(5,3)=7
            i_pow(5,4)=8

            i_pow(6,1)=1
            i_pow(6,2)=2
            i_pow(6,3)=6
            i_pow(6,4)=5

      return
      end



      subroutine data3d (kel,Iss,XX,YY,ZZ,NAX,NAY,NAZ)
      implicit real*8 (a-h,o-z)
      dimension iss((nax-1)*(nay-1)*(naz-1),8)
      dimension XX(*),YY(*),ZZ(*)

c    kel      - numer elementu

       common /WSP_LOK/ xxe(8),yye(8),zze(8)
       common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband


         do j=1,8
            I_OT_ELEM(j)    =iSS(kel,j)
            xxe(j)   =XX(I_OT_ELEM(j))
            yye(j)   =YY(I_OT_ELEM(j))
            zze(j)   =ZZ(I_OT_ELEM(j))
         enddo

      return
      end

        subroutine m_m(a,b,c,n,mm,m)
        implicit real*8 (a-h,o-z)
        real*8 a(n,mm),b(mm,m),c(n,m)

           do i=1,n
           do j=1,m
             c(i,j)=0
           do k=1,mm
             c(i,j)=c(i,j)+a(i,k)*b(k,j)
           enddo
           enddo
           enddo

        return
        end

      subroutine ilo_m_v(a,b,n,c)
      implicit real*8 (a-h,o-z)
      real*8 a(n,n),b(n),c(n)
      do 103 i=1,n
            c(i)=0
            do 103 k=1,n
               c(i)=c(i)+a(i,k)*b(k)
103    continue
      return
      end

      subroutine el3d (i_czas,kel1,T,kod,nax,nay,naz)
      implicit real*8 (a-h,o-z)

c    kel1 -numer elementu

      real*8 lambda
      dimension T(*)
      dimension kod((nax-1)*(nay-1)*(naz-1),6)
      dimension zx1(4),zx2(4),zx3(4)
      dimension lambda(3),v(3)
      real*8 N1(8),DNX1(8),DNY1(8),DNZ1(8),DJAC
      real*8 DNXw(8),DNYw(8),DNZw(8)
      real*8  n1b(8),w(8)
      real*8 ksi,eta,zeta
      real*8 KSI_G,ETA_G,ZETA_G



      common /LOK/ KSI(8),ETA(8),ZETA(8)
      common /LOK_G/ KSI_G(8),ETA_G(8),ZETA_G(8)
      common /WSP_LOK/ xxe(8),yye(8),zze(8)
      common /POW_EL/ i_pow(6,4)
      common /SKL_EL/ hma(8,8),cma(8,8),fmv(8),cg(8,8),hx(8,8),hy(8,8),
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8),CMA_N(8,8),HMAT(8,8)
      common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband


         call atmosfera (i_czas,kel1,lambda,pzanik,v,T,qv,ro,cp)


         do ji=1,8
            fmv(ji)=0.0d+00
         do jj=1,8
            hma(ji,jj)=0.0d+00
            cma(ji,jj)=0.0d+00
            hmav(ji,jj)=0.0d+00
            CMA_N(ji,jj)=0.0d+00

         enddo
         enddo

      do 40 jl=1,8
            zss=KSI_G(jl)
            ztt=ETA_G(jl)
            zuu=ZETA_G(jl)

            call baza (xxe,yye,zze,zss,ztt,zuu,N1,DNX1,DNY1,DNZ1,DJAC
     &                ,w,lambda,v,dnxw,dnyw,dnzw)

            call m_m(dnxw,dnx1,hx,8,1,8)
            call m_m(dnyw,dny1,hy,8,1,8)
            call m_m(dnzw,dnz1,hz,8,1,8)

            call m_m(w,n1,cg,8,1,8)
            call m_m(w,dnx1,hvx,8,1,8)
            call m_m(w,dny1,hvy,8,1,8)
            call m_m(w,dnz1,hvz,8,1,8)


          z_nic1=DJAC*LAMBDA(1)
          z_nic2=DJAC*LAMBDA(2)
          z_nic3=DJAC*LAMBDA(3)


c zrodlo ciepla  w elemencie

               do i=1,8
               fmv(i)=fmv(i)+qv*w(i)*djac
               enddo

         do i=1,8
         do j=1,8

          hma(i,j)=hma(i,j)+hx(i,j)*z_nic1+hy(i,j)*z_nic2+hz(i,j)*z_nic3
         hmav(i,j)=hmav(i,j)+(hvx(i,j)*v(1)+hvy(i,j)*v(2)
     &             +hvz(i,j)*v(3))*djac
         cma(i,j)=cma(i,j)+cg(i,j)*djac*pzanik

C DO STANU NIESTAC
         CMA_N(i,j)=cma_N(i,j)+ro*cp*cg(i,j)*DJac

         enddo
         enddo

 40      continue

c do warunku brzegowego III rodzaju (tylko dla powierzchni-"i_wb")
c 3               ZE1=-1.
c 1               ZE3=-1.
c 4               ZE1=1.

         i_wb=3
         if ((kod(kel1,i_wb)).eq.i_wb) then
               do i2=1,4
               zx1(i2)=xxe(i_pow(i_wb,i2))
               zx2(i2)=yyE(i_pow(i_wb,i2))
               zx3(i2)=zzE(i_pow(i_wb,i2))
               ENDDO

         do 403 i=1,4
               ZE1=KSI_G(i_pow(i_wb,i))
               ze2=ETA_G(i_pow(i_wb,i))
               ze3=ZETA_G(i_pow(i_wb,i))
               ZE1=-1.0d+00

               call jak2d(zx1,zx2,zx3,detj)

               call  FUN_K(ze1,ze2,ze3,N1b,W,v,lambda)

               call m_m(w,n1b,cg,8,1,8)

               call licz_alfa (i_czas,kel1,kod,alfaa,qs,NAX,NAY,NAZ)

c fmv(*) - war. brzegowy II rodzaju  na pow. i_wb
                     z_nic=detj*alfaa
                     do i5=1,8
                        fmv(i5)=fmv(i5)-qs*w(i5)*detj
                     do  j5=1,8
                        cma(i5,j5)=cma(i5,j5)+cg(i5,j5)*z_nic
                     enddo
                     enddo
 403     continue
         endif



               do k=1,8
c               fmv(K)=fmv(K)
               do l=1,8
               HMAT(K,L)=hma(k,l)
               hma(k,l)=hma(k,l)+cma(k,l)+hmav(k,l)
               ENDDO
               ENDDO


      return
      end


      subroutine przest3d (i_czas,XX,YY,ZZ,iSS,kod,t,nax,nay,naz,
     &              bb,EPS1B,ZNORMA,IL_ITERACJI)
      implicit real*8 (a-h,o-z)

      dimension XX(*),YY(*),ZZ(*),t(*),bb(*)
      dimension Iss((nax-1)*(nay-1)*(naz-1),8)
      dimension kod((nax-1)*(nay-1)*(naz-1),6)

      dimension  TNW2(2254),hma_t(8,8),fmv_t(8),wekt_t(8)

      dimension HMA3(2254),HMA4(2254)
      character*12 wynik_t(100)

      common /WSP_LOK/ xxe(8),yye(8),zze(8)
      common /POW_EL/ i_pow(6,4)
      common /SKL_EL/ hma(8,8),cma(8,8),fmv(8),cg(8,8),hx(8,8),hy(8,8),
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8),CMA_N(8,8),HMAT(8,8)
      common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband

c      common /duza_t/ hma2(650,213),fmv2(650)

      common /gs_ps/ fmv3(2254)
      common /d_time/ i_d_czas



c      DO I=1,ilosc_wezlow
c      FMV2(I)=0.
c      DO J=1,(2*mband-1)
c      HMA2(i,J)=0.
c      enddo
c      enddo

c odczyt z pliku wektora rozwiazan
      nout17=17
      open(nout17,file='wyn.bin',form='unformatted',status='old')
      do 1141 ii2=1,ilosc_wezlow
      read(nout17)tnw2(ii2)
 1141 continue
      close(nout17)

        i_czas=i_czas+i_d_czas
        print*,'czas= ',i_czas


        do  i=1,ilosc_wezlow
        HMA3(i)=0.0d+00
        HMA4(i)=0.0d+00
        enddo


        do  i=1,(64*ilosc_el)
        bb(i)=0.0d+00
        enddo
        do  i=1,ilosc_wezlow
        fmv3(i)=0.0d+00
        enddo

      print*
      write(1,199)ilosc_el,ilosc_wezlow
 199  format('Wypelnianie macierzy...',i5,' elementow przestrzeni & '
     + ,i5,'  wezlow')

      do 100 jk=1,ilosc_el
c      print 1999
c 1999 format('±',$)
        print*,'--',jk

C     USUWANIE 3 ELEMENTOW - przeszkoda

      IF (JK.EQ.739) GOTO 100
      IF (JK.EQ.740) GOTO 100
      IF (JK.EQ.741) GOTO 100


         call data3d (JK,iSS,XX,YY,ZZ,NAX,NAY,NAZ)
         call el3d (i_czas,JK,T,kod,nax,nay,naz)

c warunki brzegowe I rodzaju
      do i=1,4
         if (kod(jk,1).eq.1) call wb_1(t(1),i_pow(1,i),hma,fmv)
         if (kod(jk,2).eq.2) call wb_1(t(2),i_pow(2,i),hma,fmv)
c         if (kod(jk,3).eq.3) call wb_1(t(3),i_pow(3,i),hma,fmv)
         if (kod(jk,4).eq.4) call wb_1(t(4),i_pow(4,i),hma,fmv)
         if (kod(jk,5).eq.5) call wb_1(t(5),i_pow(5,i),hma,fmv)
         if (kod(jk,6).eq.6) call wb_1(t(6),i_pow(6,i),hma,fmv)
      enddo



c do gradientow sprzezonych

         iin=64*(jk-1)+1
         call sklad_w (bb(iin))


        do 125 i=1,8
        do 125 j=1,8
        hma_t(i,j)=hma(i,j)-3.0d+00/i_d_czas*CMA_N(i,j)
125     continue


        do i4=1,8
        wekt_t(i4)=tnw2(I_OT_ELEM(i4))
        enddo

        call ilo_m_v(hma_t,wekt_t,8,fmv_t)


c gs
        do j=1,8
        fmv3(I_OT_ELEM(j))=fmv3(I_OT_ELEM(j))+3.0d+00*fmv(j)-fmv_t(j)
        enddo


c podstawienie dla macierzy w metodzie Gaussa
c        do i=1,8
c        fmv2(I_OT_ELEM(i))=fmv2(I_OT_ELEM(i))+fmv(i)
c        enddo

c        do I=1,8
c        DO j=1,8
c        k=I_OT_ELEM(i)
c        l=I_OT_ELEM(j)+mband-k
c        hma2(k,l)=hma2(k,l)+hma(i,j)
c        enddo
c        enddo

        do 34 I=1,8
        k=I_OT_ELEM(i)
        hma3(K)=hma3(K)+hmaT(I,I)
        hma4(K)=hma4(K)+CMA_N(I,I)

  34    continue

 100  continue

C SZUKANIE NAJKROTSZEGO CZASU
        CZAS_test=1000000
        CZAS_test1=1000000

        DO I=1,ilosc_wezlow
        IF (HMA3(I).GT.0.) CZAS_test1=HMA4(I)/HMA3(I)
        IF (CZAS_test.GT.ABS(CZAS_test1)) CZAS_test=CZAS_test1
        ENDDO
        print*
        PRINT*,'max. krok czasowy ',CZAS_test

      print*
      print 999
 999  format('Rozwiazywanie ukladu rownan...')


      print*,('Gradienty sprz.')

      call sparse(fmv3,ilosc_wezlow,tnw2,rsq,iSS,bb,EPS1B
     &            ,ZNORMA,IL_ITERACJI,nax,nay,naz)

c dane po przejsciu GS -bin
      nout14=14
      open(nout14,file='wyn.bin',form='unformatted',status='unknown')
      do 1113 ii=1,ilosc_wezlow
      write(nout14)tnw2(ii)
 1113 continue
      close(nout14)


        K_t=INT(i_czas/10)
        L_t=INT(i_czas/100)
        M_t=INT(i_czas/1000)
        wynik_t(i_czas)='t_'//CHAR(48+M_t)//CHAR(48+L_t-10*M_t)
     &   //CHAR(48+K_t-10*L_t)//CHAR(48+I_czas-10*K_t)//'.TXT'


      nout10=10
      open(nout10,file=wynik_t(i_czas),status='unknown')
      do 10 i=1,ilosc_wezlow
      write(nout10,122)Xx(I),yy(I),zz(I),(tnw2(i)*1e6)
 10   continue
 122  format (e10.4,2x,e10.4,2x,e10.4,2x,e10.4)
      close(nout10)

c      call GAUSSASY(nax,nay,naz)


      return
      end

      SUBROUTINE  wb_1(t,ii,hma,fmv)
      implicit real*8 (a-h,o-z)
      dimension hma(8,8),fmv(8)
      i=ii
      do 10 j=1,8
      hma(i,i)=1
      fmv(i)=t
      if (i.eq.j) go to 10
      fmv(j)=fmv(j)-hma(j,i)*t
      hma(i,j)=0
      hma(j,i)=0
 10   continue
      return
      end


c      SUBROUTINE GAUSSASY(nax,nay,naz)
c      IMPLICIT REAL*8 (A-H,O-Z)
c      common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband
c      common /duza_t/ hma2(650,213),fmv2(650)
c
c      NN=ilosc_wezlow
c      MM=2*MBAND-1
c      N=0
c   71 N=N+1
c      IF(N.EQ.NN) GO TO 74
c      IF(hma2(N,MBAND).EQ.0.) GO TO 71
c      MMM=MIN0(MM,NN-N+MBAND)
c      I=N
c      LL=MBAND+1
c      D=fmv2(N)/hma2(N,MBAND)
c      DO 73 L=LL,MMM
c      I=I+1
c      C=hma2(N,L)/hma2(N,MBAND)
c      IF(C.EQ.0.) GO TO 73
c      J=1
c      DO 72 K=LL,MMM
c      NJ=N+J
c      hma2(NJ,L-J)=hma2(NJ,L-J)-C*hma2(NJ,MBAND-J)
c   72 J=J+1
c      hma2(N,L)=C
c   73 fmv2(I)=fmv2(I)-D*hma2(I,MBAND-I+N)
c      fmv2(N)=D
c      GO TO 71
c   74 fmv2(N)=fmv2(N)/hma2(N,MBAND)
c   75 I=N
c      N=N-1
c      IF(N.EQ.0) RETURN
c      NNN=MIN0(NN,N+MBAND-1)
c      DO 76 L=I,NNN
c      J=L-N+MBAND
c   76 fmv2(N)=fmv2(N)-fmv2(L)*hma2(N,J)
c      GO TO 75
c      END


      SUBROUTINE   baza(XQ,YQ,ZQ,SS,TT,UU,N,DNX,DNY,DNZ,DJAC
     &                 ,w,lambda,v,dnxw,dnyw,dnzw)
C
C     + + + PURPOSE + + +
C     To evaluate the base functions and their derivatives with respect
C     to x, y, and z, and the determinant of the Jacobian at a Gaussian
C     point
C
C     + + + DUMMY ARGUMENTS + + +
      REAL*8  XQ(8),YQ(8),ZQ(8),SS,TT,UU,N(8),DNX(8),DNY(8),DNZ(8),DJAC
      real*8  DNXw(8),DNYw(8),DNZw(8)
C
C     + + + ARGUMENT DEFINITIONS + + +
C     XQ    - X-coordinate at eight nodes of the element
C     YQ    - Y-coordinate at eight nodes of the element
C     ZQ    - Z-coordinate at eight nodes of the element
C     SS    - Xsi-coordinate of the Gaussian point
C     TT    - Eta-coordinate of the Gaussian point
C     UU    - Zeta-coordinate of the Gaussian point
C     N     - Base functions associated with eight nodes of the element
C     DNX   - Partial derivative of the base function with respect to x
C     DNY   - Partial derivative of the base function with respect to y
C     DNZ   - Partial derivative of the base function with respect to z
C     DJAC  - Determinant of the Jacobian
C
C     + + + LOCAL VARIABLES + + +
      INTEGER            I
      DOUBLE PRECISION   SM,SP,TM,TP,UM,UP,DJACI,
     >                   SUM1,SUM2,SUM3,SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,
     >                   SUMI1,SUMI2,SUMI3,SUMI4,SUMI5,SUMI6,SUMI7,
     >                   SUMI8,SUMI9,DNSS(8),DNTT(8),DNUU(8)

      real*8 DNSSw(8),DNTTw(8),DNUUw(8)
      real*8 w(8),lambda(3),v(3)
      real*8 ksi,eta,zeta
      real*8 xxe,yye,zze
      real*8 alfa1,alfa2,alfa3,gamma1,gamma2,gamma3,waga1,waga2,waga3
     & ,dlug_el_z,dlug_el_y,dlug_el_x

      common /LOK/ KSI(8),ETA(8),ZETA(8)
      common /WSP_LOK/ xxe(8),yye(8),zze(8)

C
C     + + + END SPECIFICATIONS + + +
C
C     compute some grouped variables
      SM  = 1.0D0 - SS
      SP  = 1.0D0 + SS
      TM  = 1.0D0 - TT
      TP  = 1.0D0 + TT
      UM  = 1.0D0 - UU
      UP  = 1.0D0 + UU
C
C     compute base functions
      N(1) = .125D0*SM*TM*UM
      N(2) = .125D0*SP*TM*UM
      N(3) = .125D0*SP*TP*UM
      N(4) = .125D0*SM*TP*UM
      N(5) = .125D0*SM*TM*UP
      N(6) = .125D0*SP*TM*UP
      N(7) = .125D0*SP*TP*UP
      N(8) = .125D0*SM*TP*UP

      dlug_el_z=zze(5)-zze(1)
      dlug_el_y=yye(4)-yye(1)
      dlug_el_x=xxe(2)-xxe(1)

      gamma3=v(3)*dlug_el_z/lambda(3)
      gamma2=v(2)*dlug_el_y/lambda(2)
      gamma1=v(1)*dlug_el_x/lambda(1)

      if (abs(v(3)).gt.0.0d+00)
     &alfa3=(exp(gamma3*0.50d+00)+exp(-0.50d+00*gamma3))
     //(exp(gamma3*0.50d+00)-exp(-0.50d+00*gamma3))-2.0d+00/gamma3
      if (abs(v(3)).lt.(0.00000001)) alfa3=0.0d+00
      alfa3=abs(alfa3)

      if (abs(v(2)).gt.0.0d+00)
     &alfa2=(exp(gamma2*0.50d+00)+exp(-0.50d+00*gamma2))
     //(exp(gamma2*0.50d+00)-exp(-0.50d+00*gamma2))-2.0d+00/gamma2
      if (abs(v(2)).lt.(0.00000001)) alfa2=0.0d+00
      alfa2=abs(alfa2)

      if (abs(v(1)).gt.0.0d+00)
     &alfa1=(exp(gamma1*0.50d+00)+exp(-0.50d+00*gamma1))
     //(exp(gamma1*0.50d+00)-exp(-0.50d+00*gamma1))-2.0d+00/gamma1
      if (abs(v(1)).lt.(0.00000001)) alfa1=0.0d+00
      alfa1=abs(alfa1)

      do 1 j=1,8
      waga3=((1.0d+00+zeta(j)*uu)*(1.0d+00-zeta(j)*uu))
     & *(-0.750d+00)*alfa3
      waga2=((1.0d+00+eta(j)*tt)*(1.0d+00-eta(j)*tt))
     & *(-0.750d+00)*alfa2
      waga1=((1.0d+00+ksi(j)*ss)*(1.0d+00-ksi(j)*ss))
     & *(-0.750d+00)*alfa1

      w(j)=((1.0d+00+ksi(j)*ss)*0.50d+00+waga1)*((1.0d+00+
     & eta(j)*tt)*0.5d+00+waga2)*( (1.0d+00+zeta(j)*uu)*0.5d+00+waga3)
 1    continue

C
C     compute the partial derivatives of base functions with respect to
C     local coordinates xsi
      DNSS(1) = -.125D0*TM*UM
      DNSS(2) =  .125D0*TM*UM
      DNSS(3) =  .125D0*TP*UM
      DNSS(4) = -.125D0*TP*UM
      DNSS(5) = -.125D0*TM*UP
      DNSS(6) =  .125D0*TM*UP
      DNSS(7) =  .125D0*TP*UP
      DNSS(8) = -.125D0*TP*UP
C
C     compute the partial derivatives of base functions with respect to
C     local coordinates eta
      DNTT(1) = -.125D0*SM*UM
      DNTT(2) = -.125D0*SP*UM
      DNTT(3) =  .125D0*SP*UM
      DNTT(4) =  .125D0*SM*UM
      DNTT(5) = -.125D0*SM*UP
      DNTT(6) = -.125D0*SP*UP
      DNTT(7) =  .125D0*SP*UP
      DNTT(8) =  .125D0*SM*UP
C
C     compute the partial derivatives of base functions with respect to
C     local coordinates zeta
      DNUU(1) = -.125D0*SM*TM
      DNUU(2) = -.125D0*SP*TM
      DNUU(3) = -.125D0*SP*TP
      DNUU(4) = -.125D0*SM*TP
      DNUU(5) =  .125D0*SM*TM
      DNUU(6) =  .125D0*SP*TM
      DNUU(7) =  .125D0*SP*TP
      DNUU(8) =  .125D0*SM*TP


      do 2 j=1,8
      waga3=((1.0d+00+zeta(j)*uu)*(1.0d+00-zeta(j)*uu))
     & *(-0.750d+00)*alfa3
      waga2=((1.0d+00+eta(j)*tt)*(1.0d+00-eta(j)*tt))
     & *(-0.750d+00)*alfa2
      waga1=((1.0d+00+ksi(j)*ss)*(1.0d+00-ksi(j)*ss))
     & *(-0.750d+00)*alfa1

       dnssw(j)=(ksi(j)*0.5d+00+1.5d+00*ALFA1*ksi(J)**2*ss)*((1.d+00
     &  +zeta(j)*uu)*0.5d+00+waga3)*((1.d+00+eta(j)*tt)*0.5d+00+waga2)

      dnttw(j)=((1.+ksi(j)*ss)*0.5d+00+waga1)*((1.+zeta(j)*uu)
     & *0.5d+00+waga3)*
     &(eta(j)*0.5d+00+1.5d+00*ALFA2*eta(J)**2*tt)

      dnuuw(j)=(zeta(j)*0.5d+00+1.5d+00*ALFA3*zeta(J)**2*uu)*
     &((1.d+00+ksi(j)*ss)*0.5d+00+waga1)*((1.d+00+eta(j)
     & *tt)*0.5d+00+waga2)

 2    continue

C
C     initiate the nine entries of the Jacobian matrix
      SUM1 = 0.0
      SUM2 = 0.0
      SUM3 = 0.0
      SUM4 = 0.0
      SUM5 = 0.0
      SUM6 = 0.0
      SUM7 = 0.0
      SUM8 = 0.0
      SUM9 = 0.0
C
C     compute the nine entries of the Jacobian matrix
      DO 290 I = 1,8
        SUM1 = SUM1 + XQ(I)*DNSS(I)
        SUM2 = SUM2 + YQ(I)*DNSS(I)
        SUM3 = SUM3 + ZQ(I)*DNSS(I)
        SUM4 = SUM4 + XQ(I)*DNTT(I)
        SUM5 = SUM5 + YQ(I)*DNTT(I)
        SUM6 = SUM6 + ZQ(I)*DNTT(I)
        SUM7 = SUM7 + XQ(I)*DNUU(I)
        SUM8 = SUM8 + YQ(I)*DNUU(I)
        SUM9 = SUM9 + ZQ(I)*DNUU(I)
  290 CONTINUE
C
C     compute the determinant of the Jacobian matrix
      DJAC = SUM1*(SUM5*SUM9-SUM6*SUM8) + SUM2*(SUM6*SUM7-SUM4*SUM9) +
     >       SUM3*(SUM4*SUM8-SUM5*SUM7)
C
C     compute the inverse of the determinant of the Jacobian matrix
      DJACI = 1.0D0/DJAC
C
C     compute the nine entries of the inverse Jacobian matrix
      SUMI1 = DJACI*(SUM5*SUM9 - SUM6*SUM8)
      SUMI2 = DJACI*(SUM3*SUM8 - SUM2*SUM9)
      SUMI3 = DJACI*(SUM2*SUM6 - SUM3*SUM5)
      SUMI4 = DJACI*(SUM6*SUM7 - SUM4*SUM9)
      SUMI5 = DJACI*(SUM1*SUM9 - SUM3*SUM7)
      SUMI6 = DJACI*(SUM3*SUM4 - SUM1*SUM6)
      SUMI7 = DJACI*(SUM4*SUM8 - SUM5*SUM7)
      SUMI8 = DJACI*(SUM2*SUM7 - SUM1*SUM8)
      SUMI9 = DJACI*(SUM1*SUM5 - SUM2*SUM4)
C
C     compute the partial derivatives of base functions with respect to
C     global coordinate x, y, and z.
      DO 390 I = 1,8
        DNX(I) = SUMI1*DNSS(I) + SUMI2*DNTT(I) + SUMI3*DNUU(I)
        DNY(I) = SUMI4*DNSS(I) + SUMI5*DNTT(I) + SUMI6*DNUU(I)
        DNZ(I) = SUMI7*DNSS(I) + SUMI8*DNTT(I) + SUMI9*DNUU(I)
  390 CONTINUE

      DO 392 I = 1,8
        DNXw(I) = SUMI1*DNSSw(I) + SUMI2*DNTTw(I) + SUMI3*DNUUw(I)
        DNYw(I) = SUMI4*DNSSw(I) + SUMI5*DNTTw(I) + SUMI6*DNUUw(I)
        DNZw(I) = SUMI7*DNSSw(I) + SUMI8*DNTTw(I) + SUMI9*DNUUw(I)
  392 CONTINUE
C
      RETURN
      END

      SUBROUTINE   FUN_K(SS,TT,UU,N,W,v,lambda)

      REAL*8  SS,TT,UU,N(8)
C     SS    - Xsi-coordinate of the Gaussian point
C     TT    - Eta-coordinate of the Gaussian point
C     UU    - Zeta-coordinate of the Gaussian point
C     N     - Base functions associated with eight nodes of the element

      REAL*8   SM,SP,TM,TP,UM,UP
      real*8   w(8),lambda(3),v(3)
      real*8 ksi,eta,zeta
      real*8 xxe,yye,zze
      real*8 alfa1,alfa2,alfa3,gamma1,gamma2,gamma3,waga1,waga2,waga3
     & ,dlug_el_z,dlug_el_y,dlug_el_x
      common /LOK/ KSI(8),ETA(8),ZETA(8)
      common /WSP_LOK/ xxe(8),yye(8),zze(8)

C     compute some grouped variables
      SM  = 1.0D0 - SS
      SP  = 1.0D0 + SS
      TM  = 1.0D0 - TT
      TP  = 1.0D0 + TT
      UM  = 1.0D0 - UU
      UP  = 1.0D0 + UU

C     compute base functions
      N(1) = .125D0*SM*TM*UM
      N(2) = .125D0*SP*TM*UM
      N(3) = .125D0*SP*TP*UM
      N(4) = .125D0*SM*TP*UM
      N(5) = .125D0*SM*TM*UP
      N(6) = .125D0*SP*TM*UP
      N(7) = .125D0*SP*TP*UP
      N(8) = .125D0*SM*TP*UP

c
      dlug_el_z=zze(5)-zze(1)
      dlug_el_y=yye(4)-yye(1)
      dlug_el_x=xxe(2)-xxe(1)

      gamma3=v(3)*dlug_el_z/lambda(3)
      gamma2=v(2)*dlug_el_y/lambda(2)
      gamma1=v(1)*dlug_el_x/lambda(1)

      if (abs(v(3)).gt.0.0d+00)
     &alfa3=(exp(gamma3*0.50d+00)+exp(-0.50d+00*gamma3))
     //(exp(gamma3*0.50d+00)-exp(-0.50d+00*gamma3))-2.0d+00/gamma3
      if (abs(v(3)).lt.(0.00000001)) alfa3=0.0d+00
      alfa3=abs(alfa3)

      if (abs(v(2)).gt.0.0d+00)
     &alfa2=(exp(gamma2*0.50d+00)+exp(-0.50d+00*gamma2))
     //(exp(gamma2*0.50d+00)-exp(-0.50d+00*gamma2))-2.0d+00/gamma2
      if (abs(v(2)).lt.(0.00000001)) alfa2=0.0d+00
      alfa2=abs(alfa2)

      if (abs(v(1)).gt.0.0d+00)
     &alfa1=(exp(gamma1*0.50d+00)+exp(-0.50d+00*gamma1))
     //(exp(gamma1*0.50d+00)-exp(-0.50d+00*gamma1))-2.0d+00/gamma1
      if (abs(v(1)).lt.(0.00000001)) alfa1=0.0d+00
      alfa1=abs(alfa1)

      do 1 j=1,8
      waga3=((1.0d+00+zeta(j)*uu)*(1.0d+00-zeta(j)*uu))
     & *(-0.750d+00)*alfa3
      waga2=((1.0d+00+eta(j)*tt)*(1.0d+00-eta(j)*tt))
     & *(-0.750d+00)*alfa2
      waga1=((1.0d+00+ksi(j)*ss)*(1.0d+00-ksi(j)*ss))
     & *(-0.750d+00)*alfa1

      w(j)=((1.0d+00+ksi(j)*ss)*0.50d+00+waga1)*((1.0d+00+
     & eta(j)*tt)*0.5d+00+waga2)*( (1.0d+00+zeta(j)*uu)*0.5d+00+waga3)
 1    continue

      RETURN
      END

      SUBROUTINE   jak2d(XQ,YQ,ZQ,detj)

      real*8  XQ(4),YQ(4),ZQ(4)
C     XQ    - x-coordinate at four nodes of the surface segment
C     YQ    - y-coordinate at four nodes of the surface segment
C     ZQ    - z-coordinate at four nodes of the surface segment

C     + + + LOCAL VARIABLES + + +
      INTEGER            IQ,KG
      real*8             P,SS,TT,SM,SP,TM,TP,
     >                   DXDSS,DYDSS,DZDSS,DXDTT,DYDTT,DZDTT,
     >                   DETZ,DETY,DETX,DETj,
     >                   N(4),S(4),T(4),DNSS(4),DNTT(4)
C
C     + + + INTRINSICS + + +
      INTRINSIC DSQRT
C
C     + + + DATA INITIALIZATIONS + + +
      DATA P/ 0.577350269189626D0/
      DATA S/-1.0D+00, 1.0D+00, 1.0D+00,-1.0D+00/
      DATA T/-1.0D+00,-1.0D+00, 1.0D+00, 1.0D+00/
C
C     + + + END SPECIFICATIONS + + +

C     *** Perform integration with Gaussian quadrature
C
      DO KG = 1,4
C
C       determine local coordinate at the Gaussian point KG
        SS = P*S(KG)
        TT = P*T(KG)
      enddo
C
C       compute some grouped variables
        SM  = 1.0D0 - SS
        SP  = 1.0D0 + SS
        TM  = 1.0D0 - TT
        TP  = 1.0D0 + TT
C
C       compute base functions
        N(1) = 0.25D0*SM*TM
        N(2) = 0.25D0*SP*TM
        N(3) = 0.25D0*SP*TP
        N(4) = 0.25D0*SM*TP
C
C       compute partial derivatives of base functions with respect to
C       local coordinate xi
        DNSS(1) = -0.25D0*TM
        DNSS(2) =  0.25D0*TM
        DNSS(3) =  0.25D0*TP
        DNSS(4) = -0.25D0*TP
C
C       compute partial derivatives of base functions with respect to
C       local coordinate eta
        DNTT(1) = -0.25D0*SM
        DNTT(2) = -0.25D0*SP
        DNTT(3) =  0.25D0*SP
        DNTT(4) =  0.25D0*SM
C
C       initiate six entries of the
C       (partial r/partial xsi) X (partial r/partial eta).
        DXDSS = 0.0D0
        DYDSS = 0.0D0
        DZDSS = 0.0D0
        DXDTT = 0.0D0
        DYDTT = 0.0D0
        DZDTT = 0.0D0
C
C       compute six entries of the
C       (partial r/partial xsi) X (partial r/partial eta).
        DO 290 IQ = 1,4
          DXDSS = DXDSS + XQ(IQ)*DNSS(IQ)
          DYDSS = DYDSS + YQ(IQ)*DNSS(IQ)
          DZDSS = DZDSS + ZQ(IQ)*DNSS(IQ)
          DXDTT = DXDTT + XQ(IQ)*DNTT(IQ)
          DYDTT = DYDTT + YQ(IQ)*DNTT(IQ)
          DZDTT = DZDTT + ZQ(IQ)*DNTT(IQ)
  290   CONTINUE
C
C       compute the determinant of the Jacobian matrix
        DETZ =  DXDSS*DYDTT - DYDSS*DXDTT
        DETY = -DXDSS*DZDTT + DZDSS*DXDTT
        DETX =  DYDSS*DZDTT - DZDSS*DYDTT
        DETJ =  DSQRT(DETX*DETX + DETY*DETY + DETZ*DETZ)

      RETURN
      END


        subroutine sparse(b,n2,x,rsq,iSS,bb,EPS,ZNORMA,IL_ITERACJI
     &   ,nax,nay,naz)
           implicit real*8 (a-h,o-z)
           dimension g(2254),h(2254),xi(2254),xj(2254),x(*),b(*)
           dimension bb(*),iSS((nax-1)*(nay-1)*(naz-1),8)
           common /gs_ps/ fmv3(2254)

           eps2 = n2*(eps)* eps
           irst = 0

1       continue
        irst = irst+1
        call  asub(x,xi,n2,iSS,bb,nax,nay,naz)
        rp = 0.0d0
        bsq = 0.0d0

        do 10 j=1,n2
        bsq = bsq+b(j)**2
        xi(j) = xi(j)-b(j)
        rp = rp+xi(j)**2
10      continue


        call  atsub(xi,g,n2,iSS,bb,nax,nay,naz)


      do 20 j=1,n2
      g(j) = -g(j)
      h(j) = g(j)
20    continue

c      ilosc iteracji mozna zmienic
c      do 70 iter=1,(10*n2)
      do 70 iter=1,IL_ITERACJI

      call asub(h,xi,n2,iSS,bb,nax,nay,naz)
      anum = 0.0d0
      aden = 0.0d0

      do 30 j=1,n2
         anum = anum+g(j)*h(j)
         aden = aden+xi(j)**2
30       continue
      IF (aden.eq.0.0) THEN
         print*,('pause in routine SPARSE , very singular matrix')
         endif
      anum = anum/aden

      do 40 j=1,n2
         xi(j) = x(j)
         x(j) = x(j)+anum*h(j)
40       continue
      call asub(x,xj,n2,iSS,bb,nax,nay,naz)
      rsq = 0.0d0

      do 50 j=1,n2
         xj(j) = xj(j)-b(j)
         rsq = rsq+xj(j)**2
50       continue

      If   (rsq.eq.rp)   THEN
      print*,('promien poszukiwan mozna zmniejszyc-dobra dokladnosc')
      GOTO 99
      endif

      If  (rsq.le.(bsq*eps2))  THEN
      print*,('promien poszukiwan mozna zmniejszyc-dobra dokladnosc')
      GOTO 99
      endif



      IF (rsq.gt.rp) THEN
        do 60 j=1,n2
            x(j)= xi(j)
60          continue

         IF (irst.ge.3) THEN
         print*,('miesza to samo- mozna zwiekszyc irst')
         GOTO 99

         endif

         GOTO 1
         endif
c 70    continue

      rp = rsq
      call atsub(xj,xi,n2,iSS,bb,nax,nay,naz)
      gg = 0.0d0
      dgg = 0.0d0

      do 80 j=1,n2
         gg = gg+g(j)**2
         dgg = dgg+(xi(j)+g(j))*xi(j)
80       continue

      IF (gg.eq.0.0) THEN

      print*,('rozwiazanie idealne')
      GOTO 99

      endif

      gam = dgg/gg

         do 90 j=1,n2
         g(j) = -xi(j)
         h(j) = g(j)+gam*h(j)
90       continue

902    format (2x,'iter.= ',i5,' norma= ',2x,e10.4)
       print 902,iter,rsq

      if (znorma.gt.rsq) then
      print*,('norma spelniona')
      GOTO 99
      endif

70    continue

      print*,('pause in routine SPARSE')
      print*,('too many iterations')


99    return
      END



      subroutine asub(xin,xout,n,iSS,bb,nax,nay,naz)
      implicit real*8 (a-h,o-z)
      dimension xin(*),xout(*)
c      ,b3(8,8)
      dimension iSS((nax-1)*(nay-1)*(naz-1),8),bb(*)

      common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband


         do j=1,ilosc_wezlow
         xout(j)=0.0d0
         enddo

         ik=1

         do 20 je=1,ilosc_el
         do j=1,8
         I_OT_ELEM(j)=iSS(je,j)
         enddo


c         IK3=0
c         DO 5 I=1,8
c         IK3=IK3+1
c         DO 5 J=1,8
c5        B3(I,J)=Bb(J+8*(IK3-1)+64*(JE-1) )

         do 20 i=1,8
         do 20 j=1,8
         xout(I_OT_ELEM(i))=xout(I_OT_ELEM(i))
     &    +bb(ik)*xin(I_OT_ELEM(j))
 20      ik=ik+1

      return
      end

      subroutine atsub(xin,xout,n,iSS,bb,nax,nay,naz)
      implicit real*8 (a-h,o-z)
      dimension xin(*),xout(*),B3(8,8)
      dimension iSS((nax-1)*(nay-1)*(naz-1),8),bb(*)

      common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband


         do  j=1,ilosc_WEZLOW
         xout(j)=0.0d0
         enddo

         do 20 je=1,ilosc_el
         do j=1,8
         I_OT_ELEM(j)=iSS(je,j)
         enddo


         IK3=0
         DO 5 I=1,8
         IK3=IK3+1
         DO 5 J=1,8
5        B3(I,J)=Bb(J+8*(IK3-1)+64*(JE-1) )

         do 20 i=1,8
            do 20 j=1,8
            xout(I_OT_ELEM(i))=xout(I_OT_ELEM(i))+
     &       b3(J,I)*xin(I_OT_ELEM(j))

20     continue
      return
      end

      subroutine sklad_w (bb)
      implicit real*8 (a-h,o-z)
      dimension bb(*)
      common /SKL_EL/ hma(8,8),cma(8,8),fmv(8),cg(8,8),hx(8,8),hy(8,8),
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8),CMA_N(8,8),HMAT(8,8)

      common /d_time/ i_d_czas

         k=1
         do i=1,8
         do j=1,8
         bb(k)=2.0d0*hma(i,j)+3.0d0/i_d_czas*CMA_N(i,j)
         k=k+1
         enddo
         enddo
      return
      end


