c uklad NIEST. LU TEST - sprawdzony liczy caly czas macierze
c wszystko OK!!!

      program f_m_LU
      parameter (nax=5,nay=10,naz=13)
      implicit real*8 (a-h,o-z)
      dimension dx(nax),dy(nay),dz(naz)
      dimension XX(nax*nay*naz),YY(nax*nay*naz),ZZ(nax*nay*naz)
      dimension t(6)
      dimension iSS((nax-1)*(nay-1)*(naz-1),8)
      dimension kod((nax-1)*(nay-1)*(naz-1),6)
      dimension bb((nax-1)*(nay-1)*(naz-1)*64)

C      dimension W_W((nax-1)*(nay-1)*(naz-1)*8)

      common /duza_tB/ hma2B(169,650),fmv2B(650)

      COMMON /TNW2_OK/ TNW2(650),W_W(3456)

      common /d_time/ i_d_czas

      i_d_czas=1


      call par_3d
      call wsp_na_o (nax,nay,naz,dx,dy,dz)
      call siatka   (XX,YY,ZZ,dx,dy,dz,nax,nay,naz)
      call w_el     (iSS,nax,nay,naz)
      call kod_el   (nax,nay,naz,kod)

      i_czas=0
1010  continue
      call przest3d(i_czas,XX,YY,ZZ,iSS,kod,t,nax,nay,naz,bb)
      if (i_czas.le.200) goto 1010


      end

      subroutine wsp_na_o(nax,nay,naz,dx,dy,dz)
      implicit real*8 (a-h,o-z)
      dimension dx(*),dy(*),dz(*)
      common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband

      ilosc_el=(nax-1)*(nay-1)*(naz-1)
      ilosc_wezlow=nax*nay*naz
      mband=nax*nay+nax+2

      dx(1)=0.
      dx(2)=10.
      dx(3)=20.
      dx(4)=30.
      dx(5)=40.


      dy(1)=0.
      dy(2)=10.
      dy(3)=20.
      dy(4)=30.
      dy(5)=40.
      dy(6)=50.
      dy(7)=60.
      dy(8)=70.
      dy(9)=80.
      dy(10)=90.

      dz(1)=0.
      dz(2)=10.
      dz(3)=20.
      dz(4)=30.
      dz(5)=40.
      dz(6)=50.
      dz(7)=60.
      dz(8)=70.
      dz(9)=80.
      dz(10)=90.
      dz(11)=100.
      dz(12)=110.
      dz(13)=120.

      return
      end


      subroutine licz_alfa (i_czas,nr_el,kod,ALFA,qs,NAX,NAY,NAZ)
      implicit real*8 (a-h,o-z)
      dimension kod((nax-1)*(nay-1)*(naz-1),6)
c     qs=-alfa*Tf    Tf-temp.plynu omywajacego powierzchnie
      alfa=0.
      qs=0.
      return
      end


      subroutine atmosfera  (i_czas,nr_el,lambda,pzanik,v,t,qv,ro,cp)
      implicit real*8 (a-h,o-z)
      real*8 lambda
      dimension t(6),v(3),lambda(3)
      ro=1.
      cp=1.

      ab_pasquil=1.43
      lambda(1)=8.15
      lambda(2)=lambda(1)*ab_pasquil
      lambda(3)=lambda(1)*ab_pasquil
      pzanik=0.000693147

         v(1)=0.
         v(2)=0.
         v(3)=0.5


         t(1)=0.
         t(2)=0.
         t(3)=0.
         t(4)=0.
         t(5)=0.
         t(6)=0.

         qv=0.

         if (nr_el.eq.126) qv=8*1e-9

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

            KSI(1)=-1.
            KSI(2)=1.
            KSI(3)=1.
            KSI(4)=-1.
            KSI(5)=-1.
            KSI(6)=1.
            KSI(7)=1.
            KSI(8)=-1.

            ETA(1)=-1.
            ETA(2)=-1.
            ETA(3)=1.
            ETA(4)=1.
            ETA(5)=-1.
            ETA(6)=-1.
            ETA(7)=1.
            ETA(8)=1.

            ZETA(1)=-1.
            ZETA(2)=-1.
            ZETA(3)=-1.
            ZETA(4)=-1.
            ZETA(5)=1.
            ZETA(6)=1.
            ZETA(7)=1.
            ZETA(8)=1.

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
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8),CMA_N(8,8)
      common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband


         call atmosfera (i_czas,kel1,lambda,pzanik,v,T,qv,ro,cp)


         do ji=1,8
            fmv(ji)=0.
         do jj=1,8
            hma(ji,jj)=0.
            cma(ji,jj)=0.
            hmav(ji,jj)=0.
            CMA_N(ji,jj)=0.

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

         i_wb=4
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
               ZE1=1.

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
               hma(k,l)=hma(k,l)+cma(k,l)+hmav(k,l)
               ENDDO
               ENDDO


      return
      end


      subroutine przest3d (i_czas,XX,YY,ZZ,iSS,kod,t,nax,nay,naz
     &                    ,bb)
      implicit real*8 (a-h,o-z)

      dimension XX(*),YY(*),ZZ(*),t(*),bb(*)
      dimension Iss((nax-1)*(nay-1)*(naz-1),8)
      dimension kod((nax-1)*(nay-1)*(naz-1),6)

      dimension  hma_t(8,8),fmv_t(8),wekt_t(8)

      dimension HMA3(650),HMA4(650)
      character*12 wynik_t(100)

      dimension ipvt(650),z(650)
      COMMON /TNW2_OK/ TNW2(650),W_W(3456)
      common /WSP_LOK/ xxe(8),yye(8),zze(8)
      common /POW_EL/ i_pow(6,4)
      common /SKL_EL/ hma(8,8),cma(8,8),fmv(8),cg(8,8),hx(8,8),hy(8,8),
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8),CMA_N(8,8)
      common /CAL_EL/ ilosc_el,ilosc_wezlow,I_OT_ELEM(8),mband

      common /duza_tB/ hma2B(169,650),fmv2B(650)

      common /d_time/ i_d_czas


           i_czas=i_czas+i_d_czas
           print*,'czas= ',i_czas

C          3*MBAND-2

           DO I=1,ilosc_wezlow
           FMV2B(I)=0.
           ENDDO

           IF (I_CZAS.EQ.1) THEN
           DO I=1,ilosc_wezlow
           DO J=1,3*MBAND-2
           HMA2B(J,I)=0.
           enddo
           enddo
           ENDIF

c odczyt z pliku wektora rozwiazan
c      nout17=17
c      open(nout17,file='wyn.bin',form='unformatted',status='old')
c      do 1141 i=1,ilosc_wezlow
c      read(nout17)tnw2(i)
c 1141 continue
c      close(nout17)

       IF (I_CZAS.EQ.1) THEN
       do i=1,ilosc_wezlow
       tnw2(i)=0.
       enddo
       endif

       IF (I_CZAS.EQ.1) THEN
        do  i=1,ilosc_wezlow
        HMA3(i)=0.
        HMA4(i)=0.
        enddo
        do  i=1,(64*ilosc_el)
        bb(i)=0.
        enddo
        do  i=1,(8*ilosc_el)
        W_W(i)=0.
        enddo
       endif


      IF (I_CZAS.EQ.1) THEN
      print*
      write(1,9)ilosc_el,ilosc_wezlow
 9    format('Wypelnianie macierzy...',i5,' elementow przestrzeni i '
     + ,i5,'  wezlow')

      do 100 jk=1,ilosc_el
      print 1999
 1999 format('±',$)

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



c do LU

         iin=64*(jk-1)+1
         call sklad_w (bb(iin))

         call sklad_w2(JK)


        do 125 i=1,8
        do 125 j=1,8
        hma_t(i,j)=hma(i,j)-3./i_d_czas*CMA_N(i,j)
125     continue

c z wekt rozw
        do i4=1,8
        wekt_t(i4)=tnw2(I_OT_ELEM(i4))
        enddo

        call ilo_m_v(hma_t,wekt_t,8,fmv_t)


c podstawienie dla macierzy w metodzie Gaussa LU
        do i=1,8
        fmv2B(I_OT_ELEM(i))=fmv2B(I_OT_ELEM(i))+3.*fmv(I)-fmv_t(I)
        enddo

        IF (I_CZAS.EQ.1)THEN
        M=2*MBAND-1
        do I=1,8
        DO j=1,8
        hma2B(I_OT_ELEM(I)-I_OT_ELEM(j)+M,I_OT_ELEM(j))
     &   =hma2B(I_OT_ELEM(I)-I_OT_ELEM(j)+M,I_OT_ELEM(j))
     &   +2.*hma(i,j)+3./i_d_czas*CMA_N(i,j)
        enddo
        enddo
        ENDIF

        do I=1,8
        k=I_OT_ELEM(i)
        hma3(K)=hma3(K)+hma(I,I)
        hma4(K)=hma4(K)+CMA_N(I,I)
        ENDDO

 100  continue

c z poczatku petli 100
         endif

c na probe

        if (i_czas.NE.1) then

         do i=1,ilosc_wezlow
         fmv2B(I)=0.
         enddo

         do 102 iel=1,ilosc_el

         do j=1,8
         I_OT_ELEM(j)=iSS(iel,j)
         enddo

         IK3=0
         DO 5 I=1,8
         IK3=IK3+1
         DO 5 J=1,8
5        hma_t(I,J)=Bb(J+8*(IK3-1)+64*(iel-1) )

         DO 6 i=1,8
6        fmv(I)=w_w(i+8*(iel-1))

         do i4=1,8
         wekt_t(i4)=tnw2(I_OT_ELEM(i4))
         enddo

         call ilo_m_v(hma_t,wekt_t,8,fmv_t)


c podstawienie dla macierzy w metodzie Gaussa LU
        do i=1,8
        fmv2B(I_OT_ELEM(i))=fmv2B(I_OT_ELEM(i))+3.*fmv(I)-fmv_t(I)
        enddo
 102    continue
        endif

c dla czasu roznego od 1 sek.

C SZUKANIE NAJKROTSZEGO CZASU
c        CZAS_test=1000000
c        CZAS_test1=1000000
c
c        DO I=1,ilosc_el
c        IF (HMA3(I).GT.0.) CZAS_test1=HMA4(I)/HMA3(I)
c        IF (CZAS_test.GT.ABS(CZAS_test1)) CZAS_test=CZAS_test1
c        ENDDO
c        print*
c        PRINT*,'max. krok czasowy ',CZAS_test

      print*
      print 999
 999  format('Rozwiazywanie ukladu rownan LU ')

      lda=3*MBAND-2
      n=ilosc_wezlow
      ml=MBAND-1
      mu=MBAND-1
      job=0

      IF (I_CZAS.EQ.1)
     &      call SGBCO (HMA2B, LDA, N, ML, MU, IPVT, RCOND, Z)

      call dgbsl(HMA2B,lda,n,ml,mu,ipvt,FMV2B,job)

c dane po przejsciu GS -bin
c      nout14=14
c      open(nout14,file='wyn.bin',form='unformatted',status='unknown')
c      do 1113 ii=1,ilosc_wezlow
c      write(nout14)FMV2B(ii)
c 1113 continue
c      close(nout14)

        do i=1,ilosc_wezlow
        tnw2(i)=FMV2B(i)
        enddo

        if (i_czas.eq.100) then
        K_t=INT(i_czas/10)
        L_t=INT(i_czas/100)
        M_t=INT(i_czas/1000)
        wynik_t(i_czas)='t_'//CHAR(48+M_t)//CHAR(48+L_t-10*M_t)
     &   //CHAR(48+K_t-10*L_t)//CHAR(48+I_czas-10*K_t)//'.TXT'

      nout10=10
      open(nout10,file=wynik_t(i_czas),status='unknown')
      do 10 i=1,ilosc_wezlow
      write(nout10,122)Xx(I),yy(I),zz(I),(FMV2B(i)*1e6)
 10   continue
 122  format (e10.4,2x,e10.4,2x,e10.4,2x,e10.4)
      close(nout10)
      stop
      endif

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
           dimension g(650),h(650),xi(650),xj(650),x(*),b(*)
           dimension bb(*),iSS((nax-1)*(nay-1)*(naz-1),8)
           common /gs_ps/ fmv3(650)

           eps2 = n2*(eps)* eps
           irst = 0

1       continue
        irst = irst+1
        call  asub(x,xi,n2,iSS,bb,nax,nay,naz)
        rp = 0.0
        bsq = 0.0

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
      anum = 0.0
      aden = 0.0

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
      rsq = 0.0

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
      gg = 0.0
      dgg = 0.0

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
         xout(j)=0.
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
         xout(j)=0.
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
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8),CMA_N(8,8)

      common /d_time/ i_d_czas

         k=1
         do i=1,8
         do j=1,8
         bb(k)=hma(i,j)-3./i_d_czas*CMA_N(i,j)
         k=k+1
         enddo
         enddo

      return
      end

      subroutine sklad_w2(JK)
      implicit real*8 (a-h,o-z)
      COMMON /TNW2_OK/ TNW2(650),W_W(3456)
      common /SKL_EL/ hma(8,8),cma(8,8),fmv(8),cg(8,8),hx(8,8),hy(8,8),
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8),CMA_N(8,8)

      common /d_time/ i_d_czas

        INN2=8*(jk-1)+1

         K=1
         do I=1,8
         W_W(8*(jk-1)+K)=FMV(I)
         k=k+1
         enddo

      return
      end


      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)
c
c     dgbsl solves the double precision band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgbco or dgbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran min0
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end


*DECK SGBCO
      SUBROUTINE SGBCO (ABD, LDA, N, ML, MU, IPVT, RCOND, Z)

C***PURPOSE  Factor a band matrix by Gaussian elimination and
C            estimate the condition number of the matrix.
C
C     SBGCO factors a real band matrix by Gaussian
C     elimination and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, SGBFA is slightly faster.
C     To solve  A*X = B , follow SBGCO by SGBSL.
C     To compute  INVERSE(A)*C , follow SBGCO by SGBSL.
C     To compute  DETERMINANT(A) , follow SBGCO by SGBDI.
C
C     On Entry
C
C        ABD     REAL(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT. N .
C                More efficient if  ML .LE. MU .
C
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        RCOND   REAL
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       REAL(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX(1, J-MU)
C                      I2 = MIN(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C     Example:  If the original matrix is
C
C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66
C
C      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain
C
C            *  *  *  +  +  +  , * = not used
C            *  * 13 24 35 46  , + = used for pivoting
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  SASUM, SAXPY, SDOT, SGBFA, SSCAL
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SGBCO
      INTEGER LDA,N,ML,MU,IPVT(*)
      REAL*8 ABD(LDA,*),Z(*)
      REAL*8 RCOND
C
      REAL*8 SDOT,EK,T,WK,WKM
      REAL*8 ANORM,S,SASUM,SM,YNORM
      INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  SGBCO
      ANORM = 0.0E0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = MAX(ANORM,SASUM(L,ABD(IS,J),1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C
C     FACTOR
C
      CALL SGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(ABD(M,K))) GO TO 30
            S = ABS(ABD(M,K))/ABS(EK-Z(K))
            CALL SSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (ABD(M,K) .EQ. 0.0E0) GO TO 40
            WK = WK/ABD(M,K)
            WKM = WKM/ABD(M,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0E0
            WKM = 1.0E0
   50    CONTINUE
         KP1 = K + 1
         JU = MIN(MAX(JU,MU+IPVT(K)),N)
         MM = M
         IF (KP1 .GT. JU) GO TO 90
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + ABS(Z(J)+WKM*ABD(MM,J))
               Z(J) = Z(J) + WK*ABD(MM,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN(ML,N-K)
         IF (K .LT. N) Z(K) = Z(K) + SDOT(LM,ABD(M+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0) GO TO 110
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN(ML,N-K)
         IF (K .LT. N) CALL SAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0) GO TO 130
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = W
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .LE. ABS(ABD(M,K))) GO TO 150
            S = ABS(ABD(M,K))/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (ABD(M,K) .NE. 0.0E0) Z(K) = Z(K)/ABD(M,K)
         IF (ABD(M,K) .EQ. 0.0E0) Z(K) = 1.0E0
         LM = MIN(K,M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL SAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END


      SUBROUTINE DAXPY( N, DA, DX, INCX, DY, INCY )
*
*     constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*
*     .. Scalar Arguments ..
      INTEGER           INCX, INCY, N
      DOUBLE PRECISION  DA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  DX( 1 ), DY( 1 )
*     ..
*     .. Local Scalars ..
      INTEGER           I, IX, IY, M, MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC         MOD
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
      IF( DA.EQ.0.0D0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 )
     $   GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF( INCX.LT.0 )
     $   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 )
     $   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         DY( IY ) = DY( IY ) + DA*DX( IX )
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD( N, 4 )
      IF( M.EQ.0 )
     $   GO TO 40
      DO 30 I = 1, M
         DY( I ) = DY( I ) + DA*DX( I )
   30 CONTINUE
      IF( N.LT.4 )
     $   RETURN
   40 MP1 = M + 1
      DO 50 I = MP1, N, 4
         DY( I ) = DY( I ) + DA*DX( I )
         DY( I+1 ) = DY( I+1 ) + DA*DX( I+1 )
         DY( I+2 ) = DY( I+2 ) + DA*DX( I+2 )
         DY( I+3 ) = DY( I+3 ) + DA*DX( I+3 )
   50 CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION DDOT( N, DX, INCX, DY, INCY )
*
*     forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*
*     .. Scalar Arguments ..
      INTEGER                         INCX, INCY, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION                DX( 1 ), DY( 1 )
*     ..
*     .. Local Scalars ..
      INTEGER                         I, IX, IY, M, MP1
      DOUBLE PRECISION                DTEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC                       MOD
*     ..
*     .. Executable Statements ..
*
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 )
     $   GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF( INCX.LT.0 )
     $   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 )
     $   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         DTEMP = DTEMP + DX( IX )*DY( IY )
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD( N, 5 )
      IF( M.EQ.0 )
     $   GO TO 40
      DO 30 I = 1, M
         DTEMP = DTEMP + DX( I )*DY( I )
   30 CONTINUE
      IF( N.LT.5 )
     $   GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1, N, 5
         DTEMP = DTEMP + DX( I )*DY( I ) + DX( I+1 )*DY( I+1 ) +
     $           DX( I+2 )*DY( I+2 ) + DX( I+3 )*DY( I+3 ) +
     $           DX( I+4 )*DY( I+4 )
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END


      subroutine sgbfa(abd,lda,n,ml,mu,ipvt,info)
      integer lda,n,ml,mu,ipvt(1),info
      real*8 abd(lda,1)
c
c     sgbfa factors a real band matrix by elimination.
c
c     sgbfa is usually called by sgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     real(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgbsl will divide by zero if
c                     called.  use  rcond  in sgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c     fortran max0,min0
c
c     internal variables
c
      real*8 t
      integer i,isamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0e0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0e0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = isamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0e0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0e0/abd(m,k)
            call sscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call saxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0e0) info = n
      return
      end

C                                                                  ************
      SUBROUTINE SAXPY(N,DA,DX,INCX,DY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
      END
C                                                                     SSCAL
C                                                                  ************
      SUBROUTINE  SSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,IX
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END

C                                                                     ISAMAX
C                                                                  ************
      INTEGER FUNCTION ISAMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N
C
      ISAMAX = 0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      ISAMAX = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
         IF (ABS(DX(IX)) .LE. DMAX) GO TO 5
         ISAMAX = I
         DMAX = ABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
         IF (ABS(DX(I)) .LE. DMAX) GO TO 30
         ISAMAX = I
         DMAX = ABS(DX(I))
   30 CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION SASUM(N,DX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DTEMP
      INTEGER I,INCX,M,MP1,N,IX
C
      SASUM = 0.0D0
      DTEMP = 0.0D0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + ABS(DX(IX))
        IX = IX + INCX
   10 CONTINUE
      SASUM = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + ABS(DX(I))
   30 CONTINUE
      IF (N .LT. 6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + ABS(DX(I)) + ABS(DX(I+1)) + ABS(DX(I+2))
     *  + ABS(DX(I+3)) + ABS(DX(I+4)) + ABS(DX(I+5))
   50 CONTINUE
   60 SASUM = DTEMP
      RETURN
      END

C                                                                      SDOT
C                                                                  ************
      DOUBLE PRECISION FUNCTION SDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      SDOT = 0.0D0
      DTEMP = 0.0D0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      SDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     *   DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 SDOT = DTEMP
      RETURN
      END
