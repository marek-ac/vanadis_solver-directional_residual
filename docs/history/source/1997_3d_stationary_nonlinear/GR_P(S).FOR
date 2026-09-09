c    kwe(*)   - tablica numerow wezlow elementow, noe*8 elementow
c    x1(*)    - wspolrzedne x1 wezlow elementow , non
c    kel(*) - kod elementu (spr. czy element jest na brzegu),noe
c    d1(*) -wsp. wezlow wzdluz kierunku 1 , nx1+1


      program test_met_Gaus_i_grad_nowe
      implicit real*8 (a-h,o-z)
      dimension d1(7),d2(14),d3(23)
      dimension kel(1716),kwe(13728),kkod(8),t(6),bb(109824)
c      common /duza_t/ hma2(2254,213),fmv2(2254)
      common /conoe/ noe,non,lnd(8),mband
      common /cokel/ lkod(6)
      COMMON /XYZ/x1(2254),x2(2254),x3(2254)
      common /gs_ps/ fmv3(2254)


           eps1B=1.0e-30
           znorma=1.0e-10
           IL_ITERACJI=2000

      call tset3d
      lkod(1)=1
      lkod(2)=2
      lkod(3)=3
      lkod(4)=4
      lkod(5)=5
      lkod(6)=6

      call aux3db (nx1,nx2,nx3,d1,d2,d3)
      call grd3db (d1,d2,d3,nx1,nx2,nx3)
      call nwe3db (kwe,nx1,nx2,nx3)
      call kod3db (x1,x2,x3,kwe,kel,nx1,nx2,nx3)
      call temp3d (x1,x2,x3,kwe,kel,kkod,t,bb,EPS1B,ZNORMA,IL_ITERACJI)

c      do Gaussa
c      nout10=10
c      open(nout10,file='wyn_ga.tmp',status='unknown')
c      do 10 i=1,non
c      write(nout10,122)X1(I),X2(I),X3(I),fmv2(i)
c 10   continue
c 122  format (e10.4,2x,e10.4,2x,e10.4,2x,e10.4)
c      close(nout10)
      end

      subroutine aux3db(nx1,nx2,nx3,d1,d2,d3)
      implicit real*8 (a-h,o-z)
      dimension d1(*),d2(*),d3(*)
      common /conoe/ noe,non,lnd(8),mband
c
c    kwe(*) - tablica numerow wezlow elementow,nel*8 elementow
c    kel(*) - kod elementu (spr. czy element jest na brzegu)
c    nx1,nx2,nx3 - liczba elementow w kier. 1,2,3
c    d1(*),d2(*),d3(*) -wsp. wezlow wzdluz kierunkow 1,2,3
c    x1(*)    - wspolrzedne x1 wezlow elementow
c    x2(*)    - wspolrzedne x2 wezlow elementow
c    x3(*)    - wspolrzedne x3 wezlow elementow
c
      nx1=6
      nx2=13
      nx3=22
      noe=nx1*nx2*nx3
      non=(nx1+1)*(nx2+1)*(nx3+1)
      mband=(nx1+1)*(nx2+1)+nx1+3
c
      d1(1)=0.
      d1(2)=20.
      d1(3)=35
      d1(4)=85
      d1(5)=135.
      d1(6)=230.
      d1(7)=420.


      d2(1)=-650.
      d2(2)=-520.
      d2(3)=-390.
      d2(4)=-270.
      d2(5)=-180.
      d2(6)=-90.
      d2(7)=-25.
      d2(8)=25.
      d2(9)=90.
      d2(10)=180.
      d2(11)=270.
      d2(12)=390.
      d2(13)=520.
      d2(14)=650.

      d3(1)=-250.
      d3(2)=-100.
      d3(3)=-50.
      d3(4)=-25.
      d3(5)=25.
      d3(6)=80.
      d3(7)=150.
      d3(8)=250.
      d3(9)=450.
      d3(10)=600.
      d3(11)=750.
      d3(12)=900.
      d3(13)=1100.
      d3(14)=1300.
      d3(15)=1500.
      d3(16)=1800.
      d3(17)=2100.
      d3(18)=2400.
      d3(19)=2700.
      d3(20)=3000.
      d3(21)=3300.
      d3(22)=3600.
      d3(23)=4000.


      return
      end


      subroutine alfint (kel,kkod,ALFA,qs)
      implicit real*8 (a-h,o-z)
c     qs=-alfa*Tf    Tf-temp.plynu omywajacego powierzchnie
c zskala - przeskalpwanie wymiarow
      alfa=0.0
      qs=0.
      return
      end
c
c
      subroutine alrola(kel,lambdaX,lambdaY,lambdaZ
     & ,pzanik,vx,vy,vz,t,qv)
      implicit real*8 (a-h,o-z)
      dimension t(6)
      real*8 lambdaX,lambdaY,lambdaZ,vx,vy,vz
      common /coeld/ x123(24)
      COMMON /XYZ/x1(2254),x2(2254),x3(2254)
      common /kostka/tnw3(1716)

      ab_pasquil=1.43
      lambdaX=8.15
      lambdaY=lambdaX*ab_pasquil
      lambdaZ=lambdaX*ab_pasquil


c      pzanik=0.000693147


      pzanik=4e-6   +dabs(tnw3(kel))*2e-4+( dabs(tnw3(kel)) )**2*0.1


      qv=0.
c      if ( kel.eq. 135.or.kel.eq.1500) qv=1e-6

c      if ( kel.eq.202) qv=1e-9
c      1G/S PRZY 100*100*100

      if ( kel.eq.273) qv=8*1e-9
      vx=0.
      vy=0.
      vz=1.554128
      t(1)=0.
      t(2)=0.
      t(3)=0.
      t(4)=0.
      t(5)=0.
      t(6)=0.

      return
      end

c
c
      subroutine grd3db (d1,d2,d3,nx1,nx2,nx3)
      implicit real*8 (a-h,o-z)
      dimension d1(*),d2(*),d3(*)
      COMMON /XYZ/x1(2254),x2(2254),x3(2254)
c     generacja siatki elementow
      ix1=nx1+1
      ix2=nx2+1
      ix3=nx3+1
      i=1
c  tworzenie siatki rombowej
c       al=d2(1)
c       a_kat=atan(d2(1)/al)
c       print*,2*a_kat*180/3.141592636
c      nout1=1
c      open(nout1,file='a.tmp',status='unknown')
      do 10 j3=1,ix3
         do 10 j2=1,ix2
            do 10 j1=1,ix1
            x1(i)=d1(j1)
c            x2(i)=d2(j2)+d2(j2)/al*d3(j3)
             x2(i)=d2(j2)
            x3(i)=d3(j3)
c      write(nout1,101)i,x1(i),x2(i),x3(i)
 10   i=i+1
c 101  format (2x,' el. ',i5,2x,3f10.2)
c      close(nout1)
      return
      end
c
c
        subroutine nwe3db (kwe,nx1,nx2,nx3)
        dimension kwe(*)
c
c podprogram podstawia do tablicy kwe(*)  numery wezlow
c otaczajacych dany element
c nx1 ilosc elementow w kierunku x1
c nx2 ilosc elementow w kierunku x2
c nx3 ilosc elementow w kierunku x3
c
        i=1
        iel=nx1*nx2*nx3
        iw1=nx1+1
        ipl=iw1*(nx2+1)
        do 1 j3=1,nx3
           ipo=ipl*(j3-1)
           do 1 j2=1,nx2
           do 1 j1=1,nx1
              il=ipo+j2*iw1+j1
              kwe(i)=il-iw1
              kwe(i+iel)=il-iw1+1
              kwe(i+2*iel)=il+1
              kwe(i+3*iel)=il
              kwe(i+4*iel)=kwe(i)+ipl
              kwe(i+5*iel)=kwe(i+iel)+ipl
              kwe(i+6*iel)=kwe(i+2*iel)+ipl
              kwe(i+7*iel)=kwe(i+3*iel)+ipl
    1   i=i+1
        return
        end
c
c
      subroutine kod3db (x1,x2,x3,kwe,kel,nx1,nx2,nx3)
      implicit real*8 (a-h,o-z)
c
c
      dimension x1(*),x2(*),x3(*),kwe(*),kel(*)
c
c
      common /cokel/ lkod(6)
      common /conoe/ noe,non,lnd(8),mband
c
c
      ikod1=lkod(1)
      ikod2=lkod(2)
      ikod3=lkod(3)
      ikod4=lkod(4)
      ikod5=lkod(5)
      ikod6=lkod(6)
c
      do 70 i=1,noe
         kel(i)=0
 70    continue
      ipl=nx1*nx2
      do 80 j=1,ipl
         kel(j)=ikod1
 80    continue
      i=ipl*(nx3-1)
      do 90 j=1,ipl
         kel(i+j)=10*kel(i+j)+ikod2
 90    continue
      i=0
      do 100 j3=1,nx3
         do 99 j=1,ipl,nx1
 99       kel(i+j)=10*kel(i+j)+ikod3
 100   i=i+ipl
       i=0
       do 110 j3=1,nx3
         do 109 j=nx1,ipl,nx1
 109      kel(i+j)=10*kel(i+j)+ikod4
 110   i=i+ipl
      i=ipl-nx1
      do 120 j3=1,nx3
         do 119 j=1,nx1
 119      kel(j+i)=10*kel(j+i)+ikod5
 120   i=i+ipl
       i=0
       do 130 j3=1,nx3
         do 129 j=1,nx1
 129      kel(j+i)=10*kel(j+i)+ikod6
 130   i=i+ipl
c
c
c      print 190, (j,j=1,8)
C       nout2=2
C      open(nout2,file='b.tmp',status='unknown')
C      write(nout2,190) (j,j=1,8)
      do 150 je=1,noe
C         write(nout2,160) je,(kwe(je+(j-1)*noe),j=1,8),kel(je)
         i1=kwe(je)
         i2=kwe(je+noe)
         i3=kwe(je+2*noe)
         i4=kwe(je+3*noe)
         i5=kwe(je+4*noe)
         i6=kwe(je+5*noe)
         i7=kwe(je+6*noe)
         i8=kwe(je+7*noe)
C         write(nout2,180) x1(i1),x1(i2),x1(i3),x1(i4),x1(i5),x1(i6)
C     +   ,x1(i7),x1(i8)
C         write(nout2,180) x2(i1),x2(i2),x2(i3),x2(i4),x2(i5),x2(i6)
C     +   ,x2(i7),x2(i8)
C         write(nout2,180) x3(i1),x3(i2),x3(i3),x3(i4),x3(i5),x3(i6)
C     +   ,x3(i7),x3(i8)
 150   continue
C     close(nout2)
C        DO 777 L=1,20
C        print 334,L,(kwe(l+(j-1)*noe),j=1,8)
C 777    CONTINue
c 334    format('kwe(',I4,')=',2X,8i6)

      return

C 160   format(1x,i5,8x, 8i6,5x,i8)
C 180   format(5x,8f9.2)
C 190   format(///1x,'el.',3x,'wezel',8i6,10x,'kod'/)

       end


      subroutine tset3d
      implicit real*8 (a-h,o-z)

      common /coeta/ eta1(8),eta2(8),eta3(8),psi1(4),psi2(4)
      common /coetg/ etg1(8),etg2(8),etg3(8)
      common /cokod/ nelb(24)
      common /conoe/ noe,non,lnd(8),mband
      do 1 j=1,8
         eta1(j)=-1.
         eta2(j)=-1.
 1    eta3(j)=-1.
      eta1(2)=1.
      eta1(3)=1.
      eta1(6)=1.
      eta1(7)=1.
c
      eta2(3)=1.
      eta2(4)=1.
      eta2(7)=1.
      eta2(8)=1.
c
      eta3(5)=1.
      eta3(6)=1.
      eta3(7)=1.
      eta3(8)=1.
c
      psi1(1)=-1.
      psi1(2)= 1.
      psi1(3)= 1.
      psi1(4)=-1.
c
      psi2(1)=-1.
      psi2(2)=-1.
      psi2(3)=1.
      psi2(4)=1.
c
c   common /coetg/ --------
c
      zmult=0.577350269189626D0
      do 2 j=1,8
         etg1(j)=zmult*eta1(j)
         etg2(j)=zmult*eta2(j)
 2    etg3(j)=zmult*eta3(j)
c
c   common /cokod/ --------
c
      do 3 j=1,8
 3    nelb(j)=j
      nelb(9) =1
      nelb(10)=4
      nelb(11)=8
      nelb(12)=5
      nelb(13)=2
      nelb(14)=3
      nelb(15)=7
      nelb(16)=6
      nelb(17)=4
      nelb(18)=3
      nelb(19)=7
      nelb(20)=8
      nelb(21)=1
      nelb(22)=2
      nelb(23)=6
      nelb(24)=5

      return
      end



      subroutine data3d (kel,kwe,x1,x2,x3)
      implicit real*8 (a-h,o-z)
      dimension kwe(*)
      dimension x1(*),x2(*),x3(*)
c    kel      - numer elementu
c    kwe(*)   - tablica numerow wezlow elementow,nel*8 elementow
c    x1(*)    - wspolrzedne x1 wezlow elementow
c    x2(*)    - wspolrzedne x2 wezlow elementow
c    x3(*)    - wspolrzedne x3 wezlow elementow
c
       common /coeld/ x123(24)
       common /conoe/ noe,non,lnd(8),mband
c

      do 1 j=1,8
         iwel=kwe(kel+(j-1)*noe)
         lnd(j)    =iwel
         x123(j)   =x1(iwel)
         x123(j+8) =x2(iwel)
         x123(j+16)=x3(iwel)
1     continue
      return
      end
c
c
      subroutine ilovec (a,b,n,c)
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      c=0.
      do 1 j=1,n
 1    c=c+a(j)*b(j)
      return
      end
c
c
      subroutine ilo(a,b,n,nm,m,c)
      implicit real*8 (a-h,o-z)
      real*8 a(n,nm),b(nm,m),c(n,m)
      do 10 i=1,n
         do 10 j=1,m
            c(i,j)=0
            do 10 k=1,nm
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
10    continue
      return
      end
c
c
      subroutine dkod (kod,ikod,kkod)
      implicit real*8 (a-h,o-z)
      dimension kkod(8)

      do 56 i=1,8
 56   kkod(i)=0

      k0=kod
      ikod=0
10    k1=k0/10
      ikod=ikod+1
      kkod(ikod)=k0-10*k1
      if(kkod(ikod).eq.0)     ikod=ikod-1
      k0=k1
      k0=k1
      if(k0.ne.0) go to 10

      if (ikod.eq.2) then
      kk1=kkod(1)
      kkod(1)=kkod(2)
      kkod(2)=kk1
      endif
      if (ikod.eq.3) then
      kk1=kkod(1)
      kkod(1)=kkod(3)
      kkod(3)=kk1
      endif
      do 66 i=1,6
      j=kkod(i)
      if (j.eq.i) go to 66
      if (j.eq.0) go to 66
      kkod(j+2)=kkod(j+1)
      kkod(j+1)=kkod(j)
      kkod(j)=j
      kkod(i)=0
 66   continue
      kkod(7)=0
      kkod(8)=0
      return
      end
c
c
      subroutine fj2(x1,x2,et1,et2,detj)
      implicit real*8 (a-h,o-z)
c
      real*8 jak(4),n(4)
      dimension x1(*),x2(*),dnj(8)
      common /coeta/ eta1(8),eta2(8),eta3(8),psi1(4),psi2(4)
c
c     x1(*)  - tablica 4 elementowa, wspolrzedne wezlow
c     x2(*)  - tablica 4 elementowa, wspolrzedne wezlow
c     et1,et2- wspolrzedne lokalne punktu w elemencie
c     detj   - wyznacznik jakobianu przeksztalcenia
c
      do 1 j=1,4
         dnj(j)   =psi1(j)*(1.+psi2(j)*et2)/4.
         dnj(j+4) =psi2(j)*(1.+psi1(j)*et1)/4.
 1    n(j)=(1.+psi1(j)*et1)*(1.+psi2(j)*et2)/4.
      i=1
      do 4 ji=1,2
         do 4 jj=1,2
            jak(i)=0.
            do 3 jl=1,4
               if (ji.eq.1) then
                  jak(i)=jak(i)+dnj(jl+4*(jj-1))*x1(jl)
                  else
                     jak(i)=jak(i)+dnj(jl+4*(jj-1))*x2(jl)
                endif
 3      continue
 4    i=i+1
      detj=jak(1)*jak(4)-jak(2)*jak(3)
      return
      end

      subroutine fn3(x123,et1,et2,et3,n)
      implicit real*8 (a-h,o-z)
      real*8 n(*),x123(*)
      common /coeta/ eta1(8),eta2(8),eta3(8),psi1(4),psi2(4)

      do 10 jj=1,8
      n(jj)=(1.+eta1(jj)*et1)*(1.+eta2(jj)*et2)*(1.+eta3(jj)*et3)/8.
 10   continue
      return
      end

      subroutine fw3(x123,et1,et2,et3,w,vz,vy,vx,
     &lambdaX,lambdaY,lambdaZ)
      implicit real*8 (a-h,o-z)
      real*8 w(*),x123(*),lambdaX,lambdaY,lambdaZ,gamma3,gamma2,gamma1
     &,vx,vy,vz
      common /coeta/ eta1(8),eta2(8),eta3(8),psi1(4),psi2(4)
c
      dlug_el_z=x123(21)-x123(17)
      dlug_el_y=x123(12)-x123(9)
      dlug_el_x=x123(2)-x123(1)

      gamma3=vz*dlug_el_z/lambdaZ
      gamma2=vy*dlug_el_y/lambdaY
      gamma1=vx*dlug_el_x/lambdaX

      if (dabs(vz).gt.0.)
     &alfa3=(exp(gamma3*0.5)+exp(-0.5*gamma3))
     //(exp(gamma3*0.5)-exp(-0.5*gamma3))-2/gamma3
      if (dabs(vz).lt.((1e-12))) alfa3=0.
      alfa3=dabs(alfa3)

      if (dabs(vy).gt.0.)
     &alfa2=(exp(gamma2*0.5)+exp(-0.5*gamma2))
     //(exp(gamma2*0.5)-exp(-0.5*gamma2))-2/gamma2
      if (dabs(vy).lt.((1e-12))) alfa2=0.
      alfa2=dabs(alfa2)

      if (dabs(vx).gt.0.)
     &alfa1=(exp(gamma1*0.5)+exp(-0.5*gamma1))
     //(exp(gamma1*0.5)-exp(-0.5*gamma1))-2/gamma1
      if (dabs(vx).lt.((1e-12))) alfa1=0.
      alfa1=dabs(alfa1)

      do 1 j=1,8
      waga3=((1.+eta3(j)*et3)*(1.-eta3(j)*et3))*(-3)/4*alfa3
      waga2=((1.+eta2(j)*et2)*(1.-eta2(j)*et2))*(-3)/4*alfa2
      waga1=((1.+eta1(j)*et1)*(1.-eta1(j)*et1))*(-3)/4*alfa1

      w(j)=((1.+eta1(j)*et1)/2.+waga1)*((1.+eta2(j)*et2)/2.+waga2)*
     & ( (1.+eta3(j)*et3)/2.+waga3)
 1    continue
      return
      end
c
c
      subroutine fj3(x123,et1,et2,et3,dnj,DNJ2,jak,
     &detj,vz,vy,vx,lambdaX,lambdaY,lambdaZ)

      implicit real*8 (a-h,o-z)
c
      real*8 jak(*),x123(*),dnj(*),DNJ2(*),lambdaX,lambdaY,lambdaZ,
     &gamma3,gamma2,gamma1
      common /coeta/ eta1(8),eta2(8),eta3(8),psi1(4),psi2(4)
c
c     x123(*)- tablica 24 elementowa, wspolrzedne wezlow
c     et1,et2,et3 - wspolrzedne lokalne punktu w elemencie
c     dnj(*) - tablica 24 elementowa, pochodne funkcji ksztaltu
c     jak(*) - tablica 9  elementowa, jakobian przeksztalcenia
c     detj   - wyznacznik jakobianu przeksztalcenia
      dlug_el_z=x123(21)-x123(17)
      dlug_el_y=x123(12)-x123(9)
      dlug_el_x=x123(2)-x123(1)

      gamma3=vz*dlug_el_z/lambdaZ
      gamma2=vy*dlug_el_y/lambdaY
      gamma1=vx*dlug_el_x/lambdaX

      if (abs(vz).gt.0.)
     &alfa3=(exp(gamma3*0.5)+exp(-0.5*gamma3))
     //(exp(gamma3*0.5)-exp(-0.5*gamma3))-2/gamma3
      if (abs(vz).lt.((1e-12))) alfa3=0.
      alfa3=abs(alfa3)
      if (abs(vy).gt.0.)
     &alfa2=(exp(gamma2*0.5)+exp(-0.5*gamma2))
     //(exp(gamma2*0.5)-exp(-0.5*gamma2))-2/gamma2
      if (abs(vy).lt.((1e-12))) alfa2=0.
      alfa2=abs(alfa2)

      if (abs(vx).gt.0.)
     &alfa1=(exp(gamma1*0.5)+exp(-0.5*gamma1))
     //(exp(gamma1*0.5)-exp(-0.5*gamma1))-2/gamma1
      if (abs(vx).lt.((1e-12))) alfa1=0.
      alfa1=abs(alfa1)

         do 2 j=1,8
      waga3=((1.+eta3(j)*et3)*(1.-eta3(j)*et3))*(-3)/4*alfa3
      waga2=((1.+eta2(j)*et2)*(1.-eta2(j)*et2))*(-3)/4*alfa2
      waga1=((1.+eta1(j)*et1)*(1.-eta1(j)*et1))*(-3)/4*alfa1

       dnj2(j)=(eta1(j)/2+3/2*ALFA1*ETA1(J)**2*et1)*
     &((1.+eta3(j)*et3)/2+waga3)*((1.+eta2(j)*et2)/2+waga2)

      dnj2(j+8)=((1.+eta1(j)*et1)/2+waga1)*((1.+eta3(j)*et3)/2+waga3)*
     &(eta2(j)/2+3/2*ALFA2*ETA2(J)**2*et2)

       dnj2(j+16)=(eta3(j)/2+3/2*ALFA3*ETA3(J)**2*et3)*
     &((1.+eta1(j)*et1)/2+waga1)*((1.+eta2(j)*et2)/2+waga2)

 2    continue

         do 1 j=1,8

         dnj(j)   =eta1(j)*(1.+eta2(j)*et2)*(1.+eta3(j)*et3)/8.
         dnj(j+8) =eta2(j)*(1.+eta1(j)*et1)*(1.+eta3(j)*et3)/8.
         dnj(j+16)=eta3(j)*(1.+eta1(j)*et1)*(1.+eta2(j)*et2)/8.

 1    continue

      i=1
      do 4 ji=1,3
         do 4 jj=1,3
            jak(i)=0.
            do 3 jl=1,8
 3          jak(i)=jak(i)+dnj(jl+8*(jj-1))*x123(jl+8*(ji-1))
 4    i=i+1
c
      detj=jak(1)*jak(5)*jak(9)+jak(2)*jak(6)*jak(7)
      detj=detj+jak(3)*jak(4)*jak(8)-jak(3)*jak(5)*jak(7)
      detj=detj-jak(1)*jak(6)*jak(8)-jak(2)*jak(4)*jak(9)

      return
      end

c
      subroutine frw3(dnj,rij,dnx)
      implicit real*8 (a-h,o-z)
c
      real*8 dnj(*),rij(*),dnx(*)
      do 1 jl=1,8
         do 1 ji =1,3
            i=jl+8*(ji-1)
            dnx(i)=0.
            do 1 jj=1,3
 1    dnx(i)=dnx(i)+rij(jj+3*(ji-1))*dnj(jl+8*(jj-1))
      return
      end

      subroutine fr3(dnj,jak,detj,rij,dnx,dnj2)
      implicit real*8 (a-h,o-z)
c
      common /coeta/ eta1(8),eta2(8),eta3(8),psi1(4),psi2(4)
      real*8 jak(*),dnj(*),rij(*),dnx(*),dnj2(*)
c     dnj(*) - tablica 24 elementowa, pochodne funkcji ksztaltu
c     jak(*) - tablica 9  elementowa, jakobian przeksztalcenia
c     detj   - wyznacznik jakobianu przeksztalcenia
c     dnx(*) - tablica 24 elementowa, pochodne f.k. wzgledem x1,x2,x3
c     rij(*) - odwrotnosc jakobianu przeksztalcenia
c
      rij(1)= (jak(5)*jak(9)-jak(6)*jak(8))/detj
      rij(2)=-(jak(4)*jak(9)-jak(6)*jak(7))/detj
      rij(3)= (jak(4)*jak(8)-jak(5)*jak(7))/detj
      rij(4)=-(jak(2)*jak(9)-jak(3)*jak(8))/detj
      rij(5)= (jak(1)*jak(9)-jak(3)*jak(7))/detj
      rij(6)=-(jak(1)*jak(8)-jak(2)*jak(7))/detj
      rij(7)= (jak(2)*jak(6)-jak(3)*jak(5))/detj
      rij(8)=-(jak(1)*jak(6)-jak(3)*jak(4))/detj
      rij(9)= (jak(1)*jak(5)-jak(2)*jak(4))/detj

      do 1 jl=1,8
         do 1 ji =1,3
            i=jl+8*(ji-1)
            dnx(i)=0.
            do 1 jj=1,3
 1    dnx(i)=dnx(i)+rij(jj+3*(ji-1))*dnj(jl+8*(jj-1))
      return
      end
c

      subroutine el3d (kel1,T,kel)
c      ,kod,thxyz,tcg,tfmv)
      implicit real*8 (a-h,o-z)
c    kod(*)   - tablica zawierajaca kod elementow
c    kel1 -numer elementu

      real*8 lambdaX,lambdaY,lambdaZ
      dimension T(*),kel(*)
      dimension kkod(8),zx1(4),zx2(4)
c
      common /coetg/ etg1(8),etg2(8),etg3(8)
      common /coeld/ x123(24)
      common /cofrm/ dnj(24),DNJ2(24),dnx(24),fnkw(8),fjak(9),rij(9)
     +,fnk2(4),dwx(24),fnk(8)
      common /cokod/ nelb(24)
      common /comac/ hma(8,8),cma(8,8),fmv(8),cg(8,8),hx(8,8),hy(8,8),
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8)
      common /conoe/ noe,non,lnd(8),mband

c
      inf=8
      ii1=1
c      zvol=0.
      do 1 ji=1,inf

         fmv(ji)=0.
         do 1 jj=1,inf
            hma(ji,jj)=0.
            cma(ji,jj)=0.
            hmav(ji,jj)=0.
 1    continue
      call alrola (kel1,lambdaX,lambdaY,lambdaZ,pzanik,vx,vy,vz,T,qv)

      do 40 jl=1,inf
            zet1=etg1(jl)
            zet2=etg2(jl)
            zet3=etg3(jl)
c           call alrola (kel,lambdaX,lambdaY,lambdaZ,pzanik,vx,vy,vz,T,qv)
           call fn3 (x123,zet1,zet2,zet3,fnk)
           call fw3(x123,zet1,zet2,zet3,fnkw,vz,vy,vx,
     &      lambdaX,lambdaY,lambdaZ)
           call fj3 (x123,zet1,zet2,zet3,dnj,DNJ2,fjak,detj
     +       ,vz,vy,vx,lambdaX,lambdaY,lambdaZ)
            call fr3 (dnj,fjak,detj,rij,dnx,dnj2)
            call frw3(dnj2,rij,dwx)
            call ilo(dwx(1),dnx(1),inf,ii1,inf,hx)
            call ilo(dwx(9),dnx(9),inf,ii1,inf,hy)
            call ilo(dwx(17),dnx(17),inf,ii1,inf,hz)
            call ilo(fnkw(1),fnk(1),inf,ii1,inf,cg)
            call ilo(fnkw(1),dnx(1),inf,ii1,inf,hvx)
            call ilo(fnkw(1),dnx(9),inf,ii1,inf,hvy)
            call ilo(fnkw(1),dnx(17),inf,ii1,inf,hvz)

            zdtl=detj
            zdro=detj*pzanik
c            zvol=zvol+detj
c zrodlo ciepla  w elemencie
c               if ( kel1.eq. 135.or.kel1.eq.1500) then
               do 70 i=1,inf
 70            fmv(i)=fmv(i)+qv*fnkw(i)*detj
c               endif
                 do 37 i=1,inf
                 do 37 j=1,inf
          hma(i,j)=hma(i,j)+(hx(i,j)*LAMBDAX+hy(i,j)*LAMBDAY
     &    +hz(i,j)*LAMBDAZ)*zdtl
          hmav(i,j)=hmav(i,j)+(hvx(i,j)*vx+hvy(i,j)*vy+hvz(i,j)*vz)*detj
 37       cma(i,j)=cma(i,j)+cg(i,j)*zdro
 40      continue
c
c do warunku brzegowego III rodzaju (tylko dla dolnej powierzchni-"3")

               call dkod(kel(kel1),iii,kkod)

               if ((kkod(3)).eq.3) then
               do 50 i2=1,4
               zx1(i2)=x123(8+nelb(8+i2))
               zx2(i2)=x123(16+nelb(8+i2))
 50           continue

               do 400 i=1,4
               ZE1=-1.
               ze2=etg2(nelb(8+i))
               ze3=etg3(nelb(8+i))
               call fj2 (zx1,zx2,ze2,ze3,detj)
               call fn3 (x123,ze1,ze2,ze3,fnk)
               call fw3(x123,ze1,ze2,ze3,fnkw,vz,vy,vx,
     &          lambdaX,LAMBDAY,LAMBDAZ)
               call ilo(fnkw(1),fnk(1),INF,ii1,INF,cg)
               call alfint (kel1,kkod,alfaa,qs)
c fmv(*) - war. brzegowy II rodzaju  na trzeciej powierzchni
                  zdro=detj*alfaa
                  do 800 i5=1,inf
                     fmv(i5)=fmv(i5)-qs*fnkw(i5)*detj
                     do 800 j5=1,inf
                        cma(i5,j5)=cma(i5,j5)+cg(i5,j5)*zdro
 800               continue

 400  continue
               endif

               do 20 k=1,inf
c               fmv(i)=fmv(i)/zvol
               do 20 l=1,inf
               hma(k,l)=hma(k,l)+(cma(k,l)+hmav(k,l))
 20            continue


      return
      end


      subroutine temp3d (x1,x2,x3,kwe,kel,kkod,t,bb
     &,EPS1B,ZNORMA,IL_ITERACJI)
      implicit real*8 (a-h,o-z)
c      REAL*8 LAMBDAX,LAMBDAY,LAMBDAZ
      dimension x1(*),x2(*),x3(*),kwe(*),kel(*),kkod(*),t(*),bb(*)
      dimension  TNW2(2254),tnw4(2254),lnd2(8)

c      common /duza_t/ hma2(2254,213),fmv2(2254)
      common /coeld/ x123(24)
      common /cokod/ nelb(24)
      common /comac/ hma(8,8),cma(8,8),fmv(8),cg(8,8),hx(8,8),hy(8,8),
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8)
      common /conoe/ noe,non,lnd(8),mband
      common /gs_ps/ fmv3(2254)
      common /kostka/tnw3(1716)
      inf=8
c      call alrola (kel,lambdaX,LAMBDAY,LAMBDAZ,pzanik,vx,vy,vz,t,qv)

c      DO 90 I=1,non
c      fmv2(I)=0.
c      DO 90 J=1,(2*mband-1)
c      HMA2(i,J)=0.
c 90   CONTINUE

        do  i=1,Non
        tnw4(i)=0.
        enddo
        do  i=1,Noe
        tnw3(i)=0
        enddo



 2323   continue


        do  i=1,(64*NOE)
        bb(i)=0
        enddo
        do  i=1,NON
        fmv3(i)=0
        enddo


      do 100 jk=1,noe
C      print 1999,jk
 1999 format('--',i5)

      iel=jk
         call data3d (iel,kwe,x1,x2,x3)
         call el3d (iel,T,kel)
c         ,kod,thxyz,tcg,tfmv)

c        warunki brzegowe I rodzaju
      call dkod(kel(jk),iii,kkod)
      do 80 i=1,4
        if (kkod(5).eq.5) call wb_1(t(5),nelb(4*(kkod(5)-1)+i),hma,fmv)
        if (kkod(6).eq.6) call wb_1(t(6),nelb(4*(kkod(6)-1)+i),hma,fmv)

        if (kkod(1).eq.1) call wb_1(t(1),nelb(4*(kkod(1)-1)+i),hma,fmv)
        if (kkod(2).eq.2) call wb_1(t(2),nelb(4*(kkod(2)-1)+i),hma,fmv)
C       if (kkod(3).eq.3) call wb_1(t(3),nelb(4*(kkod(3)-1)+i),hma,fmv)
        if (kkod(4).eq.4) call wb_1(t(4),nelb(4*(kkod(4)-1)+i),hma,fmv)
  80  continue




c do gradientow sprzezonych

         ib1=inf*inf*(jk-1)+1
         call bmdv (bb(ib1),inf)


c podstawienie dla macierzy w metodzie Gaussa
c        do j=1,inf
c        fmv2(lnd(j))=fmv2(lnd(j))+fmv(j)*1e6
c        enddo
c        do 34 I=1,8
c        DO 34 j=1,8
c        k=lnd(i)
c        l=lnd(j)+mband-k
c        hma2(k,l)=hma2(k,l)+hma(i,j)
c  34    continue

c gs
        do j=1,inf
        fmv3(lnd(j))=fmv3(lnd(j))+fmv(j)*1e6
        enddo


 100  continue

c       nowe gs
c        do  i=1,(64*NOE)
c        b2(i)=b(i+non)
c        enddo

c        DO I=1,NON
c        fmv3(i)=fmv2(i)
c        ENDDO

c     odremowac przy pierwszym uruchomieniu
c      DO I=1,NON
c      TNW2(I)=0
c      ENDDO

C     jezeli mam wektor  bin to wyl
      nout14=14
      open(nout14,file='wyn.bin',form='unformatted',status='unknown')
      do 1133 ii=1,non
      write(nout14)tnw2(ii)
 1133 continue
      close(nout14)


      print*,('gradienty nowe')

c odczyt z pliku wektora rozwiazan
      nout17=17
      open(nout17,file='wyn.bin',form='unformatted',status='old')
      do 1141 ii2=1,non
      read(nout17)tnw2(ii2)
 1141 continue
      close(nout17)

      call sparse(fmv3,non,tnw2,rsq,kwe,bb,EPS1B,ZNORMA,IL_ITERACJI)

c dane po przejsciu GS -bin
      nout14=14
      open(nout14,file='wyn.bin',form='unformatted',status='unknown')
      do 1113 ii=1,non
      write(nout14)tnw2(ii)
 1113 continue
      close(nout14)

c      nout10=10
c      open(nout10,file='wyn_gs_n.tmp',status='unknown')
c      do 10 i=1,non
c      write(nout10,122)X1(I),X2(I),X3(I),tnw2(i)
c 10   continue
c 122  format (e10.4,2x,e10.4,2x,e10.4,2x,e10.4)
c      close(nout10)

        do  i=1,Noe
        tnw3(i)=0
        enddo


      do 220 je=1,noe
        do 101 j=1,inf
        lnd2(j)=kwe(je+(j-1)*noe)
 101    tnw3(je)=tnw2(lnd2(j))/8+tnw3(je)
 220  continue

      z_kontrola=0.
      do i=1,non
      z_kontrola=(dabs(tnw4(i))-dabs(tnw2(i)))**2+z_kontrola
      enddo
      print*, 'kontrola zgodnosci ',z_kontrola

c do poczatku
      do je=1,non
      tnw4(je)=tnw2(je)
      enddo

      if (z_kontrola.gt.0.00000000000001) goto 2323



      nout10=10
      open(nout10,file='wyn_gs_n.tmp',status='unknown')
      do 10 i=1,non
      write(nout10,122)X1(I),X2(I),X3(I),tnw2(i)
 10   continue
 122  format (e10.4,2x,e10.4,2x,e10.4,2x,e10.4)
      close(nout10)


c      nout11=11
c      open(nout11,file='kostki.tmp',status='unknown')
c      do 1027 i=1,noe
c      write(nout11,1222)tnw3(i)
c 1027 continue
c 1222 format (e10.4)
c      close(nout11)



c      call GAUSSASY

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


c      SUBROUTINE GAUSSASY
c      IMPLICIT REAL*8 (A-H,O-Z)
c      common /conoe/ noe,non,lnd(8),mband
cc      common /duza_t/ hma2(2254,213),fmv2(2254)

c      NN=NON
c      MM=2*MBAND-1
c      N=0
c   71 N=N+1
c      print*,n
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
c
        subroutine sparse(b,n2,x,rsq,k,bb,EPS,ZNORMA,IL_ITERACJI)
           implicit real*8 (a-h,o-z)
           dimension g(2254),h(2254),xi(2254),xj(2254),x(*),b(*)
           dimension k(*),bb(*)
           common /gs_ps/ fmv3(2254)

           eps2 = n2*eps**2
           irst = 0

1       continue
        irst = irst+1
        call  asub(x,xi,n2,k,bb)
        rp = 0.0
        bsq = 0.0

        do 10 j=1,n2
        bsq = bsq+b(j)**2
        xi(j) = xi(j)-b(j)
        rp = rp+xi(j)**2
10      continue


        call  atsub(xi,g,n2,k,bb)


      do 20 j=1,n2
      g(j) = -g(j)
      h(j) = g(j)
20    continue

c      ilosc iteracji mozna zmienic
c      do 70 iter=1,(10*n2)
      do 70 iter=1,IL_ITERACJI

      call asub(h,xi,n2,k,bb)
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
      call asub(x,xj,n2,k,bb)
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
      call atsub(xj,xi,n2,k,bb)
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
CCC       print 902,iter,rsq

      if (znorma.gt.rsq) then
      print*,('norma spelniona')
      GOTO 99
      endif

70    continue

      print*,('pause in routine SPARSE')
      print*,('too many iterations')


99    return
      END



      subroutine asub(xin,xout,n,k,bb)
      implicit real*8 (a-h,o-z)
      dimension xin(*),xout(*)
c      ,b3(8,8)
      dimension k(*),bb(*)

      common /conoe/ noe,non,lnd(8),mband
c
c      b2 macierz z elementami 64 dla kazdej kostki
c      k tablica zawierajaca numery wezlow elementow

C        do  i=1,non
C        print*,i
C        print*,K(i)
C        enddo


      inf=8

      do 1 j=1,non
 1    xout(j)=0.

      ik=1
      do 20 je=1,noe
         do 10 j=1,inf
 10      lnd(j)=k(je+(j-1)*noe)

c         IK3=0
c         DO 5 I=1,8
c         IK3=IK3+1
c         DO 5 J=1,8
c5        B3(I,J)=Bb(J+8*(IK3-1)+64*(JE-1) )

         do 20 i=1,inf
         do 20 j=1,inf
         xout(lnd(i))=xout(lnd(i))+bb(ik)*xin(lnd(j))
20       ik=ik+1

      return
      end

      subroutine atsub(xin,xout,n,k,bb)
      implicit real*8 (a-h,o-z)
      dimension xin(*),xout(*),B3(8,8)
      dimension k(*),bb(*)

      common /conoe/ noe,non,lnd(8),mband
C      common /gs_ps/ fmv3(2254)
c
c      b2 macierz z elementami 64 dla kazdej kostki
c      k tablica zawierajaca numery wezlow elementow

      inf=8
c
      do 1 j=1,non
 1    xout(j)=0.

c      ik=1

      do 20 je=1,noe
         do 10 j=1,inf
 10      lnd(j)=k(je+(j-1)*noe)

         IK3=0
         DO 5 I=1,8
         IK3=IK3+1
         DO 5 J=1,8
5        B3(I,J)=Bb(J+8*(IK3-1)+64*(JE-1) )

         do 20 i=1,inf
            do 20 j=1,inf
            xout(lnd(i))=xout(lnd(i))+b3(J,I)*xin(lnd(j))
c 20   ik=ik+1
20     continue
      return
      end

      subroutine bmdv (bb,knf)
      implicit real*8 (a-h,o-z)
      dimension bb(*)

      common /coeld/ x123(24)
      common /comac/ hma(8,8),cma(8,8),fmv(8),cg(8,8),hx(8,8),hy(8,8),
     +hz(8,8),hvx(8,8),hvy(8,8),hvz(8,8),hmav(8,8)
c      common /gs_ps/ fmv3(600)

      k=1
      do 1 i=1,knf
      do 1 j=1,knf
            bb(k)=hma(i,j)
 1    k=k+1
      return
      end
c

