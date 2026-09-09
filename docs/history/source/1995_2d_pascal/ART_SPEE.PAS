{$A+,B-,D+,E-,F-,G+,I+,L+,N+,O-,P-,Q-,R+,S+,T-,V+,X+}
{$M 39000,0,655360}

uses crt,dos;

const
dzi=2 ;
ob=0;   { polozenie zrodla emisji  0 - srodek}
ab=1;      {1 - stosowany upwinding , 0 - niestosowany upwinding}
tp1=0;
tp2=0;    {temp. na brzebu prawym - tp1 i lewym - tp2}
kx=10    ;   {wspolczynniki dyfuzji}
ky=10    ;
qv=2e-7{/dzi};
alfa=0.00  ;   {do war. III-go rodzaju}
vx=0    ;   {predkosci unoszenia}
vy=0     ;  {vy ze znakiem "-"}
pzanik=0.0000;   {wsp. zaniku}
qs=0  ;   {do war. II-go rodzaju}
             {qs=-alfa*Tf   Tf:temp. plynu omywajacego powierzchnie}
             {war.brzeg. III-go rodzaju}

n3=12{10}   ;   { N3/N4 - KSZTALT SIATKI (* N3+1=WYMIAR MACIERZY 'h' *)}
n4=504{680}   ;   { ilosc rownan . w programie to jest takze "n" }
              {n4 musi byc podzielne przez n3,a calosc przez 2}


type

  m       = array [1..n4]                             of extended;
  macierz = array [1..N4,1..N3+1]                     of DOUBLE{EXTENDED};
  m1 = array [1..2*(n3-1)*round(n4/n3-1),1..3]        of integer;
  m3 = array [1..n4]                                  of EXTENDED;
  m4 = array [1..3,1..3]                              of extended;
  m5 = array [1..5]                                   of extended;
  m6 = array [1..3]                                   of extended;
  m7 = array [1..n4,1..n3]                            of DOUBLE{EXTENDED};

var
 k2,k3,kmac,ii,
 znakali         :   integer;
 troj            :   m1;
 sx,sy           :   m3;
 h2              :   m4;
 znak            :   char;


 k,n,i,j,w1,w2,mband,nn                          : longInt;
 wyniki                                       : Char;
 p,aaa,bbb,p2                                 : m;
 txt,txt2,txt3                                : text;
 s,xmax,NIC,ww,www,pp,
 xb1,yb1,xb2,yb2,xb3,yb3,sp1,sp2,
 pole,xs,ys,dN,dN2,gamma,gamma3,okr,
 gamma2,alfa3,alfa2,dlugelx,dlugely,kq,KSR    : extended;
 h                                            : ^macierz ;
 HS                                           :^M7;
 x,y                                          : m5;
 xxx,yyy                                      : m6;



    function pole_troj(x1,y1,x2,y2,x3,y3 :extended) :extended;

             begin
             pole_troj:=( (x1*y2-x2*y1) + (x2*y3-x3*y2) + (x3*y1-x1*y3) )/2;
             end;

function min0 (a,b : integer ):integer;
begin
if a>=b then
min0:=b else min0:=a;
end;

procedure gaussasy2;

label 71,72,73,74,75,76;
var nnn,k,mm,mmm,n,nj,ll,j,i,l  : integer;
    d,c                         : {REAL}extended;
begin
mm:=2*mband-1;
n:=0;
71: inc(n);
if n=nn then goto 74;
If (h^[n,1]=0) then goto 71;

mmm:=min0(mm,nn-n+mband);
i:=n;
ll:=mband+1;

d:=p[n]/h^[n,1] ;

for l:=ll to mmm  do
begin
inc(i);

if L<=N3 THEN
c:=hS^[n,l]/h^[n,1] ELSE c:=h^[n,l-N3]/h^[n,1];

if c=0 then continue {goto 73};
j:=1;
for k:=ll to mmm do
begin
nj:=n+j;
if L-J<=N3 THEN
hS^[nj,l-j]:=hS^[nj,l-j]-c*hs^[nj,mband-j] else
h^[nj,l-j-N3]:=h^[nj,l-j-N3]-c*hs^[nj,mband-j] ;
inc(j);
72: end;

if L<=N3 THEN
hS^[n,l]:=c ELSE h^[n,l-N3]:=c ;

if MBAND-I+N<=N3 THEN
p[i]:=p[i]-d*hS^[i,mband-i+n]  ELSE p[i]:=p[i]-d*h^[i,mband-i+n-N3] ;
73:end;

p[n]:=d;
goto 71;

74: p[n]:=p[n]/h^[n,1];
75:i:=n;
dec(n);
if n=0 then exit;
nnn:=min0(nn,n+mband-1);
for l:=i to nnn do
begin
j:=l-n+mband;

if J<=N3 THEN
p[n]:=p[n]-p[l]*hS^[n,j] ELSE p[n]:=p[n]-p[l]*h^[n,j-N3];
76:end;
goto 75;
end;  {gaussasy2}


procedure SetCursor(StartLine,EndLine:integer);
{ -Set Cursor Size}
var r : registers;
begin
    with r do
     begin
       CH:=StartLine;
       CL:=EndLine;
       Ah:=1;
       intr($10,r);
    end;
end;


procedure trojkat( x1,y1,x2,y2,x3,y3 :extended;
                            var  xb1,yb1,xb2,yb2,xb3,yb3 :extended);

        begin
        xb1:=(x1+x2)/2;
        yb1:=(y1+y2)/2;
        xb2:=(x2+x3)/2;
        yb2:=(y2+y3)/2;
        xb3:=(x3+x1)/2;
        yb3:=(y3+y1)/2;

        end; {liczy srodki bokow}


procedure trojkat_sr( x1,y1,x2,y2        :extended;
                              var  xa,ya :extended);


        begin
        xa:=(x1+x2)/2;
        ya:=(y1+y2)/2;

        end;{liczy srodek boku}

        procedure trojkatsc( x1,y1,x2,y2,x3,y3 :extended;
                            var  x,y :extended);
        var xa,ya :extended;
        begin
        xa:=(x1+x2)/2;
        ya:=(y1+y2)/2;

        x:=(x3-xa)/3+xa;
        y:=(y3-ya)/3+ya; {to sa wsp. sr. ciez.}

        end;{liczy srodki ciezkosci}


begin

if round(n4/n3) mod 2 <> 0 then
begin
clrscr;
writeln('Zmien wymiary siatki');
halt;
end;


{zakladanie siatki
*******************}
SetCursor($20, 0);
ClrScr;
{wsp X}
(*
for i:= round(n4/2)+1 to round(n4/2)+n3 do
sx[i]:=5;
for i:= round(n4/2)+n3+1 to round(n4/2)+2*n3 do
sx[i]:=10;
for i:= round(n4/2)+2*n3+1 to round(n4/2)+3*n3 do
sx[i]:=15;
for i:= round(n4/2)+3*n3+1 to round(n4/2)+4*n3 do
sx[i]:=20;
for i:= round(n4/2)+4*n3+1 to round(n4/2)+5*n3 do
sx[i]:=25;
for i:= round(n4/2)+5*n3+1 to round(n4/2)+6*n3 do
sx[i]:=35;
for i:= round(n4/2)+6*n3+1 to round(n4/2)+7*n3 do
sx[i]:=45;
for i:= round(n4/2)+7*n3+1 to round(n4/2)+8*n3 do
sx[i]:=60;
for i:= round(n4/2)+8*n3+1 to round(n4/2)+9*n3 do
sx[i]:=75;
for i:= round(n4/2)+9*n3+1 to round(n4/2)+10*n3 do
sx[i]:=90;
for i:= round(n4/2)+10*n3+1 to round(n4/2)+11*n3 do
sx[i]:=120;
for i:= round(n4/2)+11*n3+1 to round(n4/2)+12*n3 do
sx[i]:=140;
for i:= round(n4/2)+12*n3+1 to round(n4/2)+13*n3 do
sx[i]:=160;
for i:= round(n4/2)+13*n3+1 to round(n4/2)+14*n3 do
sx[i]:=190;
for i:= round(n4/2)+14*n3+1 to round(n4/2)+15*n3 do
sx[i]:=220;
for i:= round(n4/2)+15*n3+1 to round(n4/2)+16*n3 do
sx[i]:=250;
for i:= round(n4/2)+16*n3+1 to round(n4/2)+17*n3 do
sx[i]:=300;
for i:= round(n4/2)+17*n3+1 to round(n4/2)+18*n3 do
sx[i]:=350;
for i:= round(n4/2)+18*n3+1 to round(n4/2)+19*n3 do
sx[i]:=400;
for i:= round(n4/2)+19*n3+1 to round(n4/2)+20*n3 do
sx[i]:=450;
for i:= round(n4/2)+20*n3+1 to round(n4/2)+21*n3 do
sx[i]:=500;
for i:= round(n4/2)+21*n3+1 to round(n4/2)+22*n3 do
sx[i]:=600;
for i:= round(n4/2)+22*n3+1 to round(n4/2)+23*n3 do
sx[i]:=700;
for i:= round(n4/2)+23*n3+1 to round(n4/2)+24*n3 do
sx[i]:=800;
for i:= round(n4/2)+24*n3+1 to round(n4/2)+25*n3 do
sx[i]:=900;
for i:= round(n4/2)+25*n3+1 to round(n4/2)+26*n3 do
sx[i]:=1000;
for i:= round(n4/2)+26*n3+1 to round(n4/2)+27*n3 do
sx[i]:=1100;
for i:= round(n4/2)+27*n3+1 to round(n4/2)+28*n3 do
sx[i]:=1200;
for i:= round(n4/2)+28*n3+1 to round(n4/2)+29*n3 do
sx[i]:=1500;
for i:= round(n4/2)+29*n3+1 to round(n4/2)+30*n3 do
sx[i]:=1800;
for i:= round(n4/2)+30*n3+1 to round(n4/2)+31*n3 do
sx[i]:=2100;
for i:= round(n4/2)+31*n3+1 to round(n4/2)+32*n3 do
sx[i]:=2500;
for i:= round(n4/2)+32*n3+1 to round(n4/2)+33*n3 do
sx[i]:=3200;
for i:= round(n4/2)+33*n3+1 to round(n4/2)+34*n3 do
sx[i]:=4000;


for j:=0 to round((n4/n3)/2)-1 do
for i:=1+n3*j to n3*(j+1) do
sx[i]:=-sx[n4-i+1]; {ta linia decyduje o wartosci wsp. X}

for j:=0 to round(n4/n3-1) do
for i:=1+n3*j to n3*(j+1) do
begin
sy[i]:=(i-1-(n3*j)+1)*10;{ta linia decyduje o wartosci wsp. Y}
if  sy[i]=100 then sy[i]:=1000;
if  sy[i]=90 then sy[i]:=700;
if  sy[i]=80 then sy[i]:=500;
if  sy[i]=70 then sy[i]:=300;
if  sy[i]=60 then sy[i]:=100;
end;
while keypressed do znak:=readkey;

writeln('wsp.X');
for i:=1 to round(n4/n3) do
write(sx[i*n3]:0:0,',');
repeat until keypressed;while keypressed do znak:=readkey;
clrscr;
writeln;
writeln('wsp.Y');
for i:=1 to n3 do
write(sy[i]:0:0,',');
repeat until keypressed;
*)


kq:=1;
writeln;
sx[1]:=0;
for j:=round((n4/n3)/2) to round(n4/n3-1) do
begin
kq:=kq*1.15;
for i:=1+n3*j to n3*(j+1) do
begin
sx[i]:=j*{j*}sqrt(j)*kq-105;
end;
end;

for i:= round(n4/2)+1 to round(n4/2)+n3 do
sx[i]:=10*dzi;

for j:=0 to round((n4/n3)/2)-1 do
for i:=1+n3*j to n3*(j+1) do
sx[i]:=-sx[n4-i+1]; {ta linia decyduje o wartosci wsp. X}

(*
writeln;
sx[1]:=0;
for j:=0 to round(n4/n3-1) do
for i:=1+n3*j to n3*(j+1) do
sx[i]:=j+1; {ta linia decyduje o wartosci wsp. X}
*)

for j:=0 to round(n4/n3-1) do
for i:=1+n3*j to n3*(j+1) do
begin
sy[i]:=(i-1-(n3*j)+1)*50;{ta linia decyduje o wartosci wsp. Y}

{
if sy[i]=10 then sy[i]:=5;
if sy[i]=20 then sy[i]:=sy[i]-355;
if sy[i]=245 then sy[i]:=400;
if sy[i]=195 then sy[i]:=250;
}

if sy[i]<=350 then sy[i]:=sy[i]/10;
if sy[i]>350 then sy[i]:=sy[i]-355;
if sy[i]=245 then sy[i]:=400;
if sy[i]=195 then sy[i]:=250;

end;
while keypressed do znak:=readkey;

writeln('wsp.X');
for i:=1 to round(n4/n3) do
write(sx[i*n3]:0:0,',');
repeat until keypressed;while keypressed do znak:=readkey;
clrscr;
writeln;
writeln('wsp.Y');
for i:=1 to n3 do
write(sy[i]:0:0,',');
repeat until keypressed;
clrscr;
{wypelnianie macierzy z nr. wierzcholkow trojkatow
**************************************************}
j:=0;
k:=1;k2:=1;
repeat

j:=j+1;
if j mod 2 <>0 then begin
troj[j,1]:=k;
troj[j,2]:=k+n3;
troj[j,3]:=k+1;
k:=k+1;
if k mod n3=0 then k:=k+1;end;

if j mod 2= 0 then begin

troj[j,1]:=k2+N3;
troj[j,2]:=k2+n3+1;
troj[j,3]:=k2+1;
k2:=k2+1;
if k2 mod n3=0 then k2:=k2+1;end;

until j=2*(n3-1)*round(n4/n3-1);
writeln;

n:=n4;
nn:=n4;
mband:=n3+1;
{wprowadzanie macierzy pasmowej 'h' oraz wektora wyrazow wolnych 'p'}
{**********************************************************************}

for i:=1 to N4 do
p[i]:=0;

new (h);
for i:=1 to N4 do
for j:=1 to N3+1 do
h^[i,j]:=0;

new (hs);
for i:=1 to n4 do
for j:=1 to n3 do
hs^[i,j]:=0;

k3:=1;
w1:=0;w2:=N3;


(*h2 - macierz 3*3 *)
writeln('Wypelnianie macierzy dyfuzji i wektora wyrazow wolnych ...');
for kmac:=1 to 2*(n3-1)*round(n4/n3-1) do



                   begin{poczatek}
                   write('°');
                   {wsp. x,y wierzcholkow trojkatow}
                   x[1]:=sx[troj[kmac,1]];
                   x[2]:=sx[troj[kmac,2]];
                   x[3]:=sx[troj[kmac,3]];
                   x[4]:=x[1];x[5]:=x[2];
                   y[1]:=sy[troj[kmac,1]];
                   y[2]:=sy[troj[kmac,2]];
                   y[3]:=sy[troj[kmac,3]];
                   y[4]:=y[1];y[5]:=y[2];

                   pole:=pole_troj(x[1],y[1],x[2],y[2],x[3],y[3]);

                   trojkat( x[1],y[1],x[2],y[2],x[3],y[3],
                    xxx[1],yyy[1],xxx[2],yyy[2],xxx[3],yyy[3]);


                   for i:=1 to 3 do
                   for j:=1 to 3 do
                   h2[i,j]:=0;


for i:=1 to 3 do
for j:=1 to 3 do

begin

                   if kmac mod 2 <> 0 then
                   begin
                   if i=1 then znakali:=-1 ;
                   if i=2 then znakali:=1;
                   if i=3 then znakali:=-1 ;
                   end;

                   if kmac mod 2 = 0 then
                   begin
                   if i=1 then znakali:=1 ;
                   if i=2 then znakali:=1;
                   if i=3 then znakali:=-1 ;
                   end;


                   dlugelx:=x[2]-x[3];
                   gamma:=vx*dlugelx/kx;
                   if vx=0 then alfa2 :=0 else
                   alfa2:=(exp(gamma*0.5)+exp(-0.5*gamma))
                   /(exp(gamma*0.5)-exp(-0.5*gamma))-2/gamma;

                   alfa2:=alfa2*ab;


for ii:=1 to 3 do
begin

if xxx[ii] = 0 then
             BEGIN
                okr:=0;
                KSR:=2;
             END

             else

              BEGIN
                  okr:=kx/(xxx[ii]);
                  KSR:=1;
              END;


IF ABS(XXX[II]) < 100 THEN
KQ:=1 ELSE KQ:=1;

{ksr:=1;okr:=0;kq:=0;}
h2[i,j]:=
         (
{KSR *}
kx* ( 0.5*(y[i+1]-y[i+2])/pole
     +znakali*abs ( (-3)/sqr(dlugelx)*alfa2*(2*abs(-xxx[ii]+x[i])-dlugelx) )
    )*
    0.5*(y[j+1]-y[j+2])/pole +

ky*  0.5*(x[i+2]-x[i+1])/pole*
     0.5*(x[j+2]-x[j+1])/pole+

(vx {-  okr} )*0.5*(y[j+1]-y[j+2])/pole*
{Ni} (   0.5*( xxx[ii]*y[i+1]-x[i+1]*yyy[ii]+ x[i+1]*y[i+2]-x[i+2]*y[i+1]
+x[i+2]*yyy[ii]-xxx[ii]*y[i+2])/pole
+znakali*abs ( (-3)/sqr(dlugelx)*alfa2*abs(-abs(xxx[ii])+abs(x[i]))
*(abs(-abs(xxx[ii])+abs(x[i]))-dlugelx))
     )+

vy*0.5*(x[j+2]-x[j+1])/pole*
{Ni} (   0.5*( xxx[ii]*y[i+1]-x[i+1]*yyy[ii]+ x[i+1]*y[i+2]-x[i+2]*y[i+1]
+x[i+2]*yyy[ii]-xxx[ii]*y[i+2])/pole
+znakali*abs ( (-3)/sqr(dlugelx)*alfa2*
abs(-abs(xxx[ii])+abs(x[i]))*(abs(-abs(xxx[ii])+abs(x[i]))-dlugelx))

     )+

pzanik*kq*
(
0.5*( xxx[ii]*y[i+1]-x[i+1]*yyy[ii] + x[i+1]*y[i+2]-x[i+2]*y[i+1]
+x[i+2]*yyy[ii]-xxx[ii]*y[i+2])/pole
+znakali*abs ( (-3)/sqr(dlugelx)*alfa2*
abs(-abs(xxx[ii])+abs(x[i]))*(abs(-abs(xxx[ii])+abs(x[i]))-dlugelx))

        )*
(
0.5*( xxx[ii]*y[j+1]-x[j+1]*yyy[ii] + x[j+1]*y[j+2]-x[j+2]*y[j+1]
+x[j+2]*yyy[ii]-xxx[ii]*y[j+2])/pole
)
         )*pole/3

+h2[i,j] ;
end;
end;
                   if (kmac-1) mod (2*n3-2) =0 then

                    begin
{warunek brzegowy 3-go rodzaju}
trojkat_sr( x[1],y[1],x[2],y[2],xs,ys);

                   dlugelx:=x[2]-x[1];
                   gamma:=vx*dlugelx/kx;
                   if vx=0 then alfa2 :=0 else
                   alfa2:=(exp(gamma*0.5)+exp(-0.5*gamma))
                   /(exp(gamma*0.5)-exp(-0.5*gamma))-2/gamma;

                   alfa2:=alfa2*ab;

znakali:=-1;
h2[1,1]:=h2[1,1]+alfa*(x[2]-x[1])/6/cos(arctan((y[2]-y[1])/(x[2]-x[1])))*
       (  1*1  + 0*0 +
 4 *   ( 0.5*(xs*y[2]-x[2]*ys+x[2]*y[3]-x[3]*y[2]+x[3]*ys-xs*y[3] )/pole
+znakali*abs ( (-3)/sqr(dlugelx)*alfa2*
abs(-abs(xs)+abs(x[1]))*(abs(-abs(xs)+abs(x[1]))-dlugelx))
 )*
       ( 0.5*(xs*y[2]-x[2]*ys+x[2]*y[3]-x[3]*y[2]+x[3]*ys-xs*y[3] )/pole)
 );

znakali:=-1;
h2[1,2]:=h2[1,2]+alfa*(x[2]-x[1])/6/cos(arctan((y[2]-y[1])/(x[2]-x[1])))*
4*( 0.5*(xs*y[2]-x[2]*ys+x[2]*y[3]-x[3]*y[2]+x[3]*ys-xs*y[3] )/pole
+znakali*abs ( (-3)/sqr(dlugelx)*alfa2*
abs(-abs(xs)+abs(x[1]))*(abs(-abs(xs)+abs(x[1]))-dlugelx))
)*
( 0.5*(x[1]*ys-xs*y[1]+xs*y[3]-x[3]*ys+x[3]*y[1]-x[1]*y[3] )/pole);

znakali:=1;
h2[2,1]:=h2[2,1]+alfa*(x[2]-x[1])/6/cos(arctan((y[2]-y[1])/(x[2]-x[1])))*
4*( 0.5*(xs*y[2]-x[2]*ys+x[2]*y[3]-x[3]*y[2]+x[3]*ys-xs*y[3] )/pole
+znakali*abs ( (-3)/sqr(dlugelx)*alfa2*
abs(-abs(xs)+abs(x[2]))*(abs(-abs(xs)+abs(x[2]))-dlugelx))
)*
( 0.5*(x[1]*ys-xs*y[1]+xs*y[3]-x[3]*ys+x[3]*y[1]-x[1]*y[3] )/pole);

znakali:=1;
h2[2,2]:=h2[2,2]+alfa*(x[2]-x[1])/6/cos(arctan((y[2]-y[1])/(x[2]-x[1])))*
(1+0+ 4*
( 0.5*(x[1]*ys-xs*y[1]+xs*y[3]-x[3]*ys+x[3]*y[1]-x[1]*y[3] )/pole
+znakali*abs ( (-3)/sqr(dlugelx)*alfa2*
abs(-abs(xs)+abs(x[2]))*(abs(-abs(xs)+abs(x[2]))-dlugelx))
)*
( 0.5*(x[1]*ys-xs*y[1]+xs*y[3]-x[3]*ys+x[3]*y[1]-x[1]*y[3] )/pole));
{wektor wyrazow wolnych}
for i:=1 to 2 do
begin
if i=1 then znakali:=-1 ;
if i=2 then znakali:=1;

p2[i]:=-0.5*qs/cos(arctan((y[2]-y[1])/(x[2]-x[1])))*(x[2]-x[1])*
(
(xs*y[i+1]-x[i+1]*ys+x[i+1]*y[i+2]-x[i+2]*y[i+1]+x[i+2]*ys-xs*y[i+2])/pole
+znakali*abs ( (-3)/sqr(dlugelx)*alfa2*
abs(-abs(xs)+abs(x[i]))*(abs(-abs(xs)+abs(x[i]))-dlugelx))
);
end;
p[k3]:=p[k3]+p2[1];
p[k3+n3]:=p2[2];
k3:=k3+n3;

                    end;



{WYPELNIANIE MACIERZY WSPOLCZYNNIKOW H[I,J]}

if kmac mod 2 <>0 then
begin
                  w1:=w1+1;
                  if w1 mod N3 = 0 then w1:=w1+1;
                  h^[w1,1]:=h2[1,1]+h^[w1,1];
                  h^[w1,N3+1]:=h2[1,2]+h^[w1,N3+1];
                  h^[w1,2]:=h2[1,3]+h^[w1,2];
                  h^[w1+N3,1]:=h2[2,2]+h^[w1+N3,1];
                  h^[w1+1,N3]:=h2[3,2]+h^[w1+1,N3];
                  h^[w1+1,1]:=h2[3,3]+h^[w1+1,1];
                  hs^[w1+n3,1]:=hs^[w1+n3,1]+h2[2,1];
                  hs^[w1+1,n3]:=hs^[w1+1,n3]+h2[3,1];
                  hs^[w1+n3,2]:=hs^[w1+n3,2]+h2[2,3];
end;

if kmac mod 2=0 then
begin
                   w2:=w2+1;
                   if w2 mod N3 = 0 then w2:=w2+1;

                   h^[w2,1]:=h2[1,1]+h^[w2,1];
                   h^[w2,2]:=h2[1,2]+h^[w2,2];
                   h^[w2-N3+1,N3]:=h2[3,1]+h^[w2-N3+1,N3];
                   h^[w2+1,1]:=h2[2,2]+h^[w2+1,1];
                   h^[w2-N3+1,N3+1]:=h2[3,2]+h^[w2-N3+1,N3+1];
                   h^[w2-N3+1,1]:=h2[3,3]+h^[w2-N3+1,1];
                   hs^[w2+1,n3]:=hs^[w2+1,n3]+h2[2,1];
                   hs^[w2,2]:=hs^[w2,2]+h2[1,3];
                   hs^[w2+1,1]:=hs^[w2+1,1]+h2[2,3];

end;
                   end;{poczatek}
{
for i:=1 to n4 do
for j:=1 to n3+1 do
writeln('h[',i,',',j,']=',h^[i,j]:0:3);
}

{******************************************************************
*******************************************************************}

                   dlugelx:=abs(sx[round(n4/2)+trunc(n3/2)-ob]
                   -sx[round(n4/2)-trunc(n3/2)-ob]);
                   gamma:=vx*dlugelx/kx;
                   if vx=0 then alfa2 :=0 else
                   alfa2:=(exp(gamma*0.5)+exp(-0.5*gamma))
                   /(exp(gamma*0.5)-exp(-0.5*gamma))-2/gamma;

                   alfa2:=alfa2*ab;

pole:=pole_troj(sx[round(n4/2)-trunc(n3/2)-ob],sy[round(n4/2)-trunc(n3/2)-ob],
           sx[round(n4/2)+trunc(n3/2)-ob],sy[round(n4/2)+trunc(n3/2)-ob],
           sx[round(n4/2)-trunc(n3/2)+1-ob],sy[round(n4/2)-trunc(n3/2)+1-ob]);

trojkatsc( sx[round(n4/2)-trunc(n3/2)-ob],sy[round(n4/2)-trunc(n3/2)-ob],
           sx[round(n4/2)+trunc(n3/2)-ob],sy[round(n4/2)+trunc(n3/2)-ob],
           sx[round(n4/2)-trunc(n3/2)+1-ob],sy[round(n4/2)-trunc(n3/2)+1-ob],
           xs,ys);

i:=round(n4/2)-trunc(n3/2)-ob ;
znakali:=-1 ;
p[round(round(n4/2)-trunc(n3/2))-ob]:=p[round(round(n4/2)-trunc(n3/2))-ob]+
qv*pole*
(1/3+znakali*
abs ( (-3)/sqr(dlugelx)*alfa2*abs(-abs(xs)+abs(sx[i]))
*(abs(-abs(xs)+abs(sx[i]))-dlugelx))
);

i:=round(n4/2)+trunc(n3/2)-ob ;
znakali:=1 ;
p[round(n4/2)+trunc(n3/2)-ob]:=p[round(n4/2)+trunc(n3/2)-ob]+
qv*pole*
(1/3+znakali*
abs ( (-3)/sqr(dlugelx)*alfa2*abs(-abs(xs)+abs(sx[i]))
*(abs(-abs(xs)+abs(sx[i]))-dlugelx))
);

i:=round(n4/2)-trunc(n3/2)+1-ob ;
znakali:=-1 ;
p[round(n4/2)-trunc(n3/2)+1-ob]:=p[round(n4/2)-trunc(n3/2)+1-ob]+
qv*pole*
(1/3+znakali*
abs ( (-3)/sqr(dlugelx)*alfa2*abs(-abs(xs)+abs(sx[i]))
*(abs(-abs(xs)+abs(sx[i]))-dlugelx))
);

trojkatsc( sx[round(n4/2)+trunc(n3/2)+1-ob],sy[round(n4/2)+trunc(n3/2)+1-ob],
           sx[round(n4/2)+trunc(n3/2)-ob],sy[round(n4/2)+trunc(n3/2)-ob],
           sx[round(n4/2)-trunc(n3/2)+1-ob],sy[round(n4/2)-trunc(n3/2)+1-ob],
           xs,ys);

i:=round(n4/2)+trunc(n3/2)-ob;
znakali:=1 ;
p[round(n4/2)+trunc(n3/2)-ob]:=p[round(n4/2)+trunc(n3/2)-ob]+
qv*pole*
(1/3+znakali*
abs ( (-3)/sqr(dlugelx)*alfa2*abs(-abs(xs)+abs(sx[i]))
*(abs(-abs(xs)+abs(sx[i]))-dlugelx))
);

i:=round(n4/2)+trunc(n3/2)+1-ob;
znakali:=1 ;
p[i]:=p[i]+
qv*pole*
(1/3+znakali*
abs ( (-3)/sqr(dlugelx)*alfa2*abs(-abs(xs)+abs(sx[i]))
*(abs(-abs(xs)+abs(sx[i]))-dlugelx))
);

i:=round(n4/2)-trunc(n3/2)+1-ob;
znakali:=-1 ;
p[round(n4/2)-trunc(n3/2)+1-ob]:=p[round(n4/2)-trunc(n3/2)+1-ob]+
qv*pole*
(1/3+znakali*
abs ( (-3)/sqr(dlugelx)*alfa2*abs(-abs(xs)+abs(sx[i]))
*(abs(-abs(xs)+abs(sx[i]))-dlugelx))
);

{*************************************************************************
**************************************************************************}

(*
{brzeg prawy i lewy}
{tutaj wprowadzam warunki brzegowe}
{Tb temp. na brzegu ciala}
for i:=1 to n3 do
begin
h^[i,1]:=h^[i,1]+(1E15);
p[i]:=p[i]+(1E15)*tp1 {Tb. st. C};
end;
{**************************}
{Tb2 temp. na brzegu ciala}

for i:=n4-n3+1 to n4 do
begin
h^[i,1]:=h^[i,1]+(1E15);
p[i]:=p[i]+(1E15)*tp2 {Tb2. st. C};
end;
{*************************}
*)

(*
{warunki brzegowe - brzeg gorny i dolny}
for i:=1 to round(n4/n3) do
begin
p[i*n3]:=p[i*n3]+tp1*1e15;
h^[i*n3,1]:=h^[i*n3,1]+1e15;
end;

for i:=0 to round(n4/n3)-1 do
begin
p[1+i*n3]:=p[i*n3+1]+tp2*1e15;
h^[i*n3+1,1]:=h^[i*n3+1,1]+1e15;
end;
{*********************************}
*)

{brzeg prawy i lewy}
{tutaj wprowadzam warunki brzegowe  I -go ;inny sp.}
{Tb temp. na brzegu ciala}

{*************************************************************
**************************************************************}
for i:=1 to n3 do
for j:=1 to n4 do
if  i<>j then
   begin
   if i>j then
          begin
          p[j]:=p[j]-h^[j,i-j+1]* tp1  {Tp1};
          end;

   if j>i then
          begin
          if n3+1-j+i>0 then
                        begin
                        p[j]:=p[j]-hs^[j,n3+1-j+i]* tp1  {Tp1};
                        hs^[j,n3+1-j+i]:=0;
                        end;
          end;

    end;


for i:=1 to n3 do
for j:=1 to n3+1 do
h^[i,j]:=0;

for i:=1 to n3 do
begin
h^[i,1]:=1;
p[i]:=tp1;  {Tp1}
end;



for i:=n4-n3+1 to n4 do
for j:=1 to n4 do
if  i<>j then
   begin
   if i>j then
          begin
          if i-j+1<n3+2 then
          begin
          p[j]:=p[j]-h^[j,i-j+1]* tp2  {Tp2};
          h^[j,i-j+1]:=0;
          end;
          end;

   if j>i then
          begin
          if n3+1-j+i>0 then
                        begin
                        p[j]:=p[j]-hs^[j,n3+1-j+i]* tp2  {Tp2};
                        end;
          end;

    end;


for i:=n4-n3+1 to n4 do
for j:=1 to n3+1 do
h^[i,j]:=0;

for i:=n4-n3+1 to n4 do
for j:=1 to n3 do
hs^[i,j]:=0;

for i:=n4-n3+1 to n4 do
begin
h^[i,1]:=1;
p[i]:=tp2;  {Tp2}
end;
{*********************************************************
**********************************************************}

(*
FOR I:=1 TO N4 DO
P[I]:=0;

for i:=1 to round(n4/n3) do
begin
h^[i*n3,1]:=h^[i*n3,1]+(1E20);
p[i*n3]:=p[i*n3]+(1E20)*TP1;
end;
*)

writeln;
writeln('... skonczone');
while keypressed do znak:=readkey;
repeat until keypressed;

{granice obszaru wewnetrznego zrodla ciepla}
clrscr;
writeln('Granice obszaru zrodla zanieczyszczen');
writeln('X [',round(n4/2)-trunc(n3/2)-ob,'] = ',sx[round(n4/2)-trunc(n3/2)-ob]:5:1);
writeln('X [',round(n4/2)-trunc(n3/2)+1-ob,'] = ',sx[round(n4/2)-trunc(n3/2)+1-ob]:5:1);
writeln('X [',round(n4/2)+trunc(n3/2)-ob,'] = ',sx[round(n4/2)+trunc(n3/2)-ob]:5:1);
writeln('X [',round(n4/2)+trunc(n3/2)+1-ob,'] = ',sx[round(n4/2)+trunc(n3/2)+1-ob]:5:1);
writeln;
writeln('Y [',round(n4/2)-trunc(n3/2)-ob,'] = ',sy[round(n4/2)-trunc(n3/2)-ob]:5:1);
writeln('Y [',round(n4/2)-trunc(n3/2)+1-ob,'] = ',sy[round(n4/2)-trunc(n3/2)+1-ob]:5:1);
writeln('Y [',round(n4/2)+trunc(n3/2)-ob,'] = ',sy[round(n4/2)+trunc(n3/2)-ob]:5:1);
writeln('Y [',round(n4/2)+trunc(n3/2)+1-ob,'] = ',sy[round(n4/2)+trunc(n3/2)+1-ob]:5:1);
while keypressed do znak:=readkey;
repeat until keypressed;
{macierze AAA oraz BBB do sprawdzenia
poprawnosci rozwiazania ukladu rownan}
SP1:=P[1];
SP2:=P[2];

FOR I:=1 TO N3+1 DO
AAA[I]:=H^[1,I];
for i:=n3+2 to n4 do
aaa[i]:=0;

FOR I:=1 TO N3+1 DO
bbb[i+1]:=H^[2,I];

for i:=n3+3 to n4 do
bbb[i]:=0;
bbb[1]:=hs^[2,n3];

 begin
  clrscr;
  writeln;
  writeln ('Rozwiazywanie ukladu rownan ...');
  gaussasy2 ;
  dispose(h);dispose(hs);
  clrscr;
  writeln('Wyniki : ');
        i:=1;
        k:=0;
        while i<=n do

        begin

        WrIteln('x[',i:3,'] = ',P[I]:0:12);
        i:=i+1;
        k:=k+1;
        if k=20 then
                begin
                k:=0;
                while KeyPressed do
                wyniki:=Readkey;
                Writeln;
                write ('kontynuacja - ');
                writeln ('dowolny klawisz');
                repeat
                until KeyPressed;
                clrscr;
                end
        end
       end;
          ww:=0;www:=0;
          for i:=1 to n4 do
          ww:=AAA[I]*P[i]+ww;
          www:=ww-Sp1;
          writeln;
          writeln('spr.wiersza NR1[=0!!] = ', www{:0:10});

          ww:=0;www:=0;
          for i:=1 to n4 do
          ww:=BBB[I]*P[i]+ww;
          www:=ww-SP2;
          writeln;
          writeln('spr.wiersza NR2[=0!!] = ', www{:0:10});

      { Wydruk wynikow do pliku }

      assign(txt,'c:\txt\S_1.txt');
      rewrite(txt);
      FOR I:=1 TO N DO
      writeln(txt,P[I]*1e6);
      close(txt);
      assign(txt2,'c:\txt\x_1.txt');
      rewrite(txt2);
      FOR I:=1 TO N DO
      writeln(txt2,sx[i]);
      close(txt2);

      assign(txt3,'c:\txt\y_1.txt');
      rewrite(txt3);
      FOR I:=1 TO N DO
      writeln(txt3,sy[i]);
      close(txt3);

{SetCursor(6,7);}
while keypressed do znak:=readkey;   repeat until keypressed;
end.
