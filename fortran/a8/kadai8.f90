program kadai8
implicit none

double precision,parameter::h=2.0d-1
double precision,parameter::a=1.2d0
double precision,parameter::b=(a+1.0d0)/(a-1.0d0)
integer,parameter::im=50
integer,parameter::jm=50
double precision,parameter::r=2.871d2       !ÉKÉXíËêî
double precision,parameter::g=1.4d0     !î‰îMî‰
double precision,parameter::m=2.9d0         !É}ÉbÉnêî
double precision,parameter::ii=1.0d0
double precision,parameter::eps=1.0d-6
double precision,parameter::pr=0.72d0
double precision,parameter::sz=110.0d0
double precision,parameter::t0=293.0d0
double precision,parameter::mu0=1.0d-3
double precision,parameter::prt=0.9d0
double precision,parameter::ret=500.0d0
double precision,parameter::cmu=0.09d0
double precision,parameter::c1=1.44d0
double precision,parameter::c2=1.92d0
double precision,parameter::sigmak=1.0d0
double precision,parameter::sigmaipu=1.3d0
double precision,parameter::cp=g*r/(g-1.0d0)   !íËà≥î‰îM
double precision,dimension(0:im,0:jm)::x
double precision,dimension(0:im,0:jm)::y
double precision,dimension(0:im,0:jm)::rho=1.2d0,p=0.0d0,u=0.0d0,v=0.0d0,&
energy=0.0d0,eb,t,c,vn1,vn2,vt,v2,un,vn,uxi,ueta,vxi,veta,txi,teta
double precision,dimension(0:im,0:jm)::xix,xiy,etax,etay,xxi,yxi,xeta,yeta,&
jj,k,ipu,mu,mut,kappat,dt=1.0d-8,nu,kwall=0.0d0,ipuwall=0.0d0
double precision,dimension(1:6,0:im,0:jm)::q,ee,ff,rv,sv,pp
double precision,dimension(1:6,0:im,0:jm,0:5)::qq
integer i,j,n,mm,kk
double precision dx,dy,jmm,p1,rho1,slope,veps,dd,&
txx,txy,tyx,tyy,r4,s4,r5,r6,s5,s6,pro,sxx,sxy,syx,syy


!main 
print*,"*****START CALCULATION*****"
call mesh
call initialsetting
call metrix

do n=0,100000  !åvéZâÒêî
do j=0,jm
   do i=0,im
      do kk=1,6
QQ(kk,i,j,0)=Q(kk,i,j)
      end do
   end do
end do

do mm=1,4
call TVDxi
call TVDeta
call calpara
call visxi
call viseta
call prod
call wall
call runge_kutta
call calculate
call ltsm
end do

print*,n,dd
if(dd<eps)exit
end do

print*,"*****END CALCULATION*****"

call yuplus
call output


contains

double precision function sland(t)   !ÉTÉUÅ[ÉâÉìÉhÇÃéÆÇ…ÇÊÇÈÉ ÇÃì±èo
double precision t
sland=mu0*((t/t0)**1.5)*((t0+sz)/(t+sz))
end function


subroutine mesh !ïsìôä‘äuäiéqê∂ê¨

do j=0,jm
   do i=0,im
      if(i<=5)then
      dx=(50.0d0*h)/dble(im)
      x(i,j)=dble(i)*dx
      else
      jmm=dble(i-5)/(dble(im-5))
      x(i,j)=x(5,j)+45.0d0*h*(b**jmm-1.0d0)/(b-1.0d0)
      end if
      jmm=dble(j)/(dble(jm))
      y(i,j)=h*(b**jmm-1.0d0)/(b-1.0d0)
   end do
end do


end subroutine mesh

subroutine initialsetting !èâä˙ê›íË
do j=0,jm
   do i=0,im
t(i,j)=293.15d0
c(i,j)=dsqrt(g*r*t(i,j)) !âπë¨
u(i,j)=m*c(i,j)
v(i,j)=0.0d0
p(i,j)=rho(i,j)*r*t(i,j)
energy(i,j)=rho(i,j)*((r*t(i,j)/(g-1.0d0))+((u(i,j)**2.0)+v(i,j)**2.0)/2.0d0)
mu(i,j)=1.82d-5
k(i,j)=1.5d0*(1.0d-2*u(i,j))**2.0
nu(i,j)=mu(i,j)/rho(i,j)
ipu(i,j)=(k(i,j)**2.0)/(nu(i,j)*ret)
   end do
end do

call wall

do i=6,im
   u(i,0)=0.0d0
   v(i,0)=0.0d0
   energy(i,0)=rho(i,0)*(r*t(i,0)/(g-1.0d0)+(u(i,0)**2.0+v(i,0)**2.0)/2.0d0)
   nu(i,0)=mu(i,0)/rho(i,0)
   k(i,0)=1.5d0*(1.0d-2*u(i,0))**2.0
   ipu(i,0)=(k(i,0)**2.0)/(nu(i,0)*ret)
end do

do j=0,jm
   do i=0,im
mut(i,j)=cmu*rho(i,j)*(k(i,j)**2.0)/ipu(i,j)
kappat(i,j)=mut(i,j)*cp/prt
   end do
end do

end subroutine initialsetting

subroutine metrix !ç¿ïWïœä∑

do j=1,jm-1
   do i=0,im
      yeta(i,j)=((y(i,j+1)-y(i,j-1))/2.0d0)
      xeta(i,j)=((x(i,j+1)-x(i,j-1))/2.0d0)
   end do
end do
do j=0,jm
   do i=1,im-1
      yxi(i,j)=((y(i+1,j)-y(i-1,j))/2.0d0)
      xxi(i,j)=((x(i+1,j)-x(i-1,j))/2.0d0)
   end do
end do
do i=0,im
   yeta(i,0)=(-3.0d0*y(i,0)+4.0d0*y(i,1)-y(i,2))/(2.0d0*ii)
   yeta(i,jm)=(3.0d0*y(i,jm)-4.0d0*y(i,jm-1)+y(i,jm-2))/(2.0d0*ii)
   xeta(i,0)=(-3.0d0*x(i,0)+4.0d0*x(i,1)-x(i,2))/(2.0d0*ii)
   xeta(i,jm)=(3.0d0*x(i,jm)-4.0d0*x(i,jm-1)+x(i,jm-2))/(2.0d0*ii)
end do
do j=0,jm
   yxi(0,j)=(-3.0d0*y(0,j)+4.0d0*y(1,j)-y(2,j))/(2.0d0*ii)
   yxi(im,j)=(3.0d0*y(im,j)-4.0d0*y(im-1,j)+y(im-2,j))/(2.0d0*ii)
   xxi(0,j)=(-3.0d0*x(0,j)+4.0d0*x(1,j)-x(2,j))/(2.0d0*ii)
   xxi(im,j)=(3.0d0*x(im,j)-4.0d0*x(im-1,j)+x(im-2,j))/(2.0d0*ii)
end do

do j=0,jm
   do i=0,im
      jj(i,j)=1.0d0/(xxi(i,j)*yeta(i,j)-yxi(i,j)*xeta(i,j))
      xix(i,j)=yeta(i,j)*jj(i,j)
      xiy(i,j)=-xeta(i,j)*jj(i,j)
      etax(i,j)=-yxi(i,j)*jj(i,j)
      etay(i,j)=xxi(i,j)*jj(i,j)
   end do
end do

do j=0,jm
   do i=0,im
      q(1,i,j)=rho(i,j)/jj(i,j)
      q(2,i,j)=(rho(i,j)*u(i,j))/jj(i,j)
      q(3,i,j)=(rho(i,j)*v(i,j))/jj(i,j)
      q(4,i,j)=energy(i,j)/jj(i,j)
      q(5,i,j)=rho(i,j)*k(i,j)/jj(i,j)
      q(6,i,j)=rho(i,j)*ipu(i,j)/jj(i,j)
   end do
end do

end subroutine metrix

subroutine calpara  !äeéÌÉpÉâÉÅÅ[É^åvéZ
do j=0,jm
   do i=1,im
mu(i,j)=sland(t(i,j))
nu(i,j)=mu(i,j)/rho(i,j)
mut(i,j)=cmu*rho(i,j)*(k(i,j)**2.0)/ipu(i,j)kappat(i,j)=mut(i,j)*cp/prt
   end do
end do
end subroutine calpara

subroutine visxi   !ÉÃï˚å¸îSê´çÄ
double precision,dimension(0:im,0:jm)::xixi,xiyi,etaxi,etayi,jji,ui,vi,ti,&
uxii,vxii,txii,uetai,vetai,tetai,rhoi,muti,kxii,ketai,ipuxii,ipuetai,ki,ipui,&
kappati

do j=1,jm-1
   do i=0,im-1
xixi(i,j)=(xix(i,j)+xix(i+1,j))/2.0d0
xiyi(i,j)=(xiy(i,j)+xiy(i+1,j))/2.0d0
etaxi(i,j)=(etax(i,j)+etax(i+1,j))/2.0d0
etayi(i,j)=(etay(i,j)+etay(i+1,j))/2.0d0
jji(i,j)=(jj(i,j)+jj(i+1,j))/2.0d0
ui(i,j)=(u(i,j)+u(i+1,j))/2.0d0
vi(i,j)=(v(i,j)+v(i+1,j))/2.0d0
ti(i,j)=(t(i,j)+t(i+1,j))/2.0d0
ki(i,j)=(k(i,j)+k(i+1,j))/2.0d0
ipui(i,j)=(ipu(i,j)+ipu(i+1,j))/2.0d0
rhoi(i,j)=(q(1,i,j)*jj(i,j)+q(1,i+1,j)*jj(i+1,j))/(jj(i,j)+jj(i+1,j))*jji(i,j)
muti(i,j)=cmu*rhoi(i,j)*(ki(i,j)**2.0)/ipui(i,j)
uxii(i,j)=u(i+1,j)-u(i,j)
vxii(i,j)=v(i+1,j)-v(i,j)
txii(i,j)=t(i+1,j)-t(i,j)
kxii(i,j)=k(i+1,j)-k(i,j)
ipuxii(i,j)=ipu(i+1,j)-ipu(i,j)
uetai(i,j)=(((u(i,j+1)+u(i+1,j+1))/2.0d0)-((u(i,j-1)+u(i+1,j-1))/2.0d0))/2.0d0
vetai(i,j)=(((v(i,j+1)+v(i+1,j+1))/2.0d0)-((v(i,j-1)+v(i+1,j-1))/2.0d0))/2.0d0
tetai(i,j)=(((t(i,j+1)+t(i+1,j+1))/2.0d0)-((t(i,j-1)+t(i+1,j-1))/2.0d0))/2.0d0
ketai(i,j)=(((k(i,j+1)+k(i+1,j+1))/2.0d0)-((k(i,j-1)+k(i+1,j-1))/2.0d0))/2.0d0
ipuetai(i,j)=(((ipu(i,j+1)+ipu(i+1,j+1))/2.0d0)-&
             ((ipu(i,j-1)+ipu(i+1,j-1))/2.0d0))/2.0d0
kappati(i,j)=muti(i,j)*cp/prt
   end do
end do

do j=1,jm-1
   do i=0,im-1
sxx=2.0d0/3.0d0*(2.0d0*(uxii(i,j)*xixi(i,j)+uetai(i,j)*etaxi(i,j))-&
                   (vxii(i,j)*xiyi(i,j)+vetai(i,j)*etayi(i,j)))
syy=2.0d0/3.0d0*(2.0d0*(vxii(i,j)*xiyi(i,j)+vetai(i,j)*etayi(i,j))-&
                   (uxii(i,j)*xixi(i,j)+uetai(i,j)*etaxi(i,j)))
sxy=((uxii(i,j)*xiyi(i,j)+uetai(i,j)*etayi(i,j))+&
        (vxii(i,j)*xixi(i,j)+vetai(i,j)*etaxi(i,j)))
syx=sxy
txx=muti(i,j)*sxx-2.0d0*rhoi(i,j)*ki(i,j)/3.0d0
txy=muti(i,j)*sxy
tyx=muti(i,j)*syx
tyy=muti(i,j)*syy-2.0d0*rhoi(i,j)*ki(i,j)/3.0d0
r4=txx*ui(i,j)+txy*vi(i,j)+(kappati(i,j)*&
   (txii(i,j)*xixi(i,j)+tetai(i,j)*etaxi(i,j)))
s4=tyx*ui(i,j)+tyy*vi(i,j)+(kappati(i,j)*&
   (txii(i,j)*xiyi(i,j)+tetai(i,j)*etayi(i,j)))
r5=muti(i,j)*(kxii(i,j)*xixi(i,j)+ketai(i,j)*etaxi(i,j))/sigmak
s5=muti(i,j)*(kxii(i,j)*xiyi(i,j)+ketai(i,j)*etayi(i,j))/sigmak
r6=muti(i,j)*(ipuxii(i,j)*xixi(i,j)+ipuetai(i,j)*etaxi(i,j))/sigmaipu
s6=muti(i,j)*(ipuxii(i,j)*xiyi(i,j)+ipuetai(i,j)*etayi(i,j))/sigmaipu
rv(1,i,j)=0.0d0
rv(2,i,j)=(xixi(i,j)*txx+xiyi(i,j)*txy)/jji(i,j)
rv(3,i,j)=(xixi(i,j)*tyx+xiyi(i,j)*tyy)/jji(i,j)
rv(4,i,j)=(xixi(i,j)*r4+xiyi(i,j)*s4)/jji(i,j)
rv(5,i,j)=(xixi(i,j)*r5+xiyi(i,j)*s5)/jji(i,j)
rv(6,i,j)=(xixi(i,j)*r6+xiyi(i,j)*s6)/jji(i,j)
   end do
end do

end subroutine visxi

subroutine viseta  !É≈ï˚å¸îSê´çÄ
double precision,dimension(0:im-1,0:jm-1)::xixe,xiye,etaxe,etaye,jje,ue,ve,te,&
uxie,vxie,txie,uetae,vetae,tetae,rhoe,mute,kxie,ketae,ipuxie,ipuetae,ke,ipue,&
kappate

do j=0,jm-1
   do i=1,im-1
xixe(i,j)=(xix(i,j)+xix(i,j+1))/2.0d0
xiye(i,j)=(xiy(i,j)+xiy(i,j+1))/2.0d0
etaxe(i,j)=(etax(i,j)+etax(i,j+1))/2.0d0
etaye(i,j)=(etay(i,j)+etay(i,j+1))/2.0d0
jje(i,j)=(jj(i,j)+jj(i,j+1))/2.0d0
ue(i,j)=(u(i,j)+u(i,j+1))/2.0d0
ve(i,j)=(v(i,j)+v(i,j+1))/2.0d0
te(i,j)=(t(i,j)+t(i,j+1))/2.0d0
ke(i,j)=(k(i,j)+k(i,j+1))/2.0d0
ipue(i,j)=(ipu(i,j)+ipu(i,j+1))/2.0d0
rhoe(i,j)=(q(1,i,j)*jj(i,j)+q(1,i,j+1)*jj(i,j+1))/(jj(i,j)+jj(i,j+1))*jje(i,j)
mute(i,j)=cmu*rhoe(i,j)*(ke(i,j)**2.0)/ipue(i,j)
uetae(i,j)=u(i,j+1)-u(i,j)
vetae(i,j)=v(i,j+1)-v(i,j)
tetae(i,j)=t(i,j+1)-t(i,j)
ketae(i,j)=k(i,j+1)-k(i,j)
ipuetae(i,j)=ipu(i,j+1)-ipu(i,j)
uxie(i,j)=(((u(i+1,j)+u(i+1,j+1))/2.0d0)-((u(i-1,j)+u(i-1,j+1))/2.0d0))/2.0d0
vxie(i,j)=(((v(i+1,j)+v(i+1,j+1))/2.0d0)-((v(i-1,j)+v(i-1,j+1))/2.0d0))/2.0d0
txie(i,j)=(((t(i+1,j)+t(i+1,j+1))/2.0d0)-((t(i-1,j)+t(i-1,j+1))/2.0d0))/2.0d0
kxie(i,j)=(((k(i+1,j)+k(i+1,j+1))/2.0d0)-((k(i-1,j)+k(i-1,j+1))/2.0d0))/2.0d0
ipuxie(i,j)=(((ipu(i+1,j)+ipu(i+1,j+1))/2.0d0)-&
            ((ipu(i-1,j)+ipu(i-1,j+1))/2.0d0))/2.0d0
kappate(i,j)=mute(i,j)*cp/prt
   end do
end do

do j=0,jm-1
   do i=1,im-1
sxx=2.0d0/3.0d0*(2.0d0*(uxie(i,j)*xixe(i,j)+uetae(i,j)*etaxe(i,j))-&
                   (vxie(i,j)*xiye(i,j)+vetae(i,j)*etaye(i,j)))
syy=2.0d0/3.0d0*(2.0d0*(vxie(i,j)*xiye(i,j)+vetae(i,j)*etaye(i,j))-&
                   (uxie(i,j)*xixe(i,j)+uetae(i,j)*etaxe(i,j)))
sxy=((uxie(i,j)*xiye(i,j)+uetae(i,j)*etaye(i,j))+&
        (vxie(i,j)*xixe(i,j)+vetae(i,j)*etaxe(i,j)))
syx=sxy
txx=mute(i,j)*sxx-2.0d0*rhoe(i,j)*ke(i,j)/3.0d0
txy=mute(i,j)*sxy
tyx=mute(i,j)*syx
tyy=mute(i,j)*syy-2.0d0*rhoe(i,j)*ke(i,j)/3.0d0
r4=txx*ue(i,j)+txy*ve(i,j)+(kappate(i,j)*&
        (txie(i,j)*xixe(i,j)+tetae(i,j)*etaxe(i,j)))
s4=tyx*ue(i,j)+tyy*ve(i,j)+(kappate(i,j)*&
        (txie(i,j)*xiye(i,j)+tetae(i,j)*etaye(i,j)))
r5=mute(i,j)*(kxie(i,j)*xixe(i,j)+ketae(i,j)*etaxe(i,j))/sigmak
s5=mute(i,j)*(kxie(i,j)*xiye(i,j)+ketae(i,j)*etaye(i,j))/sigmak
r6=mute(i,j)*(ipuxie(i,j)*xixe(i,j)+ipuetae(i,j)*etaxe(i,j))/sigmaipu
s6=mute(i,j)*(ipuxie(i,j)*xiye(i,j)+ipuetae(i,j)*etaye(i,j))/sigmaipu
sv(1,i,j)=0.0d0
sv(2,i,j)=(etaxe(i,j)*txx+etaye(i,j)*txy)/jje(i,j)
sv(3,i,j)=(etaxe(i,j)*tyx+etaye(i,j)*tyy)/jje(i,j)
sv(4,i,j)=(etaxe(i,j)*r4+etaye(i,j)*s4)/jje(i,j)
sv(5,i,j)=(etaxe(i,j)*r5+etaye(i,j)*s5)/jje(i,j)
sv(6,i,j)=(etaxe(i,j)*r6+etaye(i,j)*s6)/jje(i,j)
   end do
end do

end subroutine viseta


subroutine prod   !ê∂ê¨çÄåvéZ
do j=1,jm-1
   do i=1,im-1
uxi(i,j)=(u(i+1,j)-u(i-1,j))/2.0d0
vxi(i,j)=(v(i+1,j)-v(i-1,j))/2.0d0
ueta(i,j)=(u(i,j+1)-u(i,j-1))/2.0d0
veta(i,j)=(v(i,j+1)-v(i,j-1))/2.0d0
sxx=2.0d0/3.0d0*(2.0d0*(uxi(i,j)*xix(i,j)+ueta(i,j)*etax(i,j))-&
                   (vxi(i,j)*xiy(i,j)+veta(i,j)*etay(i,j)))
syy=2.0d0/3.0d0*(2.0d0*(vxi(i,j)*xiy(i,j)+veta(i,j)*etay(i,j))-&
                   (uxi(i,j)*xix(i,j)+ueta(i,j)*etax(i,j)))
sxy=((uxi(i,j)*xiy(i,j)+ueta(i,j)*etay(i,j))+&
        (vxi(i,j)*xix(i,j)+veta(i,j)*etax(i,j)))
syx=sxy
txx=mut(i,j)*sxx-2.0d0/3.0d0*rho(i,j)*k(i,j)
txy=mut(i,j)*sxy
tyx=mut(i,j)*syx
tyy=mut(i,j)*sxx-2.0d0/3.0d0*rho(i,j)*k(i,j)
pro=txx*(uxi(i,j)*xix(i,j)+ueta(i,j)*etax(i,j))&
+txy*((vxi(i,j)*xix(i,j)+veta(i,j)*etax(i,j))&
+(uxi(i,j)*xiy(i,j)+ueta(i,j)*etay(i,j)))&
+tyy*(vxi(i,j)*xiy(i,j)+veta(i,j)*etay(i,j))

pp(1,i,j)=0.0d0
pp(2,i,j)=0.0d0
pp(3,i,j)=0.0d0
pp(4,i,j)=0.0d0
pp(5,i,j)=(pro-rho(i,j)*ipu(i,j))/jj(i,j)
pp(6,i,j)=(c1*ipu(i,j)*pro/k(i,j)-c2*rho(i,j)*(ipu(i,j)**2.0)/k(i,j))/jj(i,j)
   end do
end do

end subroutine prod

subroutine wall   !ï«ñ@ë•
double precision::kp,ipup,ypp,utau,yp
double precision,parameter::kappar=0.42d0
double precision,parameter::largee=9.7d0
do i=5,im
yp=y(i,1)
kp=k(i,1)
ipup=ipu(i,1)
utau=(cmu*(kp**2.0))/(kappar*yp*ipup)
ypp=yp*utau/nu(i,1)
if(ypp>11.63d0)then
utau=kappar*u(i,1)/dlog(ypp*largee)
else
utau=dsqrt(nu(i,1)*u(i,1)/yp)
end if

k(i,1)=(utau**2.0)/dsqrt(cmu)
ipu(i,1)=(utau**3.0)/(kappar*yp)
kwall(i,1)=k(i,1)
ipuwall(i,1)=ipu(i,1)
end do

end subroutine wall


!---------------------------------------------------------------------
!*********************************************************************
!***** 				TVD				******
!*********************************************************************
subroutine TVDxi

double precision::RWL,RWR,AKL,AKX,AKY,AJACM,UM,VM,DUM,CM,PHI,BETA,&
                  AKXT,AKYT,THIT,AKM,SI
double precision,dimension(1:4,1:4,0:im)::RR
double precision,dimension(1:4,1:4,0:im)::RI
double precision,dimension(1:6,0:im)::EIGM
double precision,dimension(1:6)::D
double precision,dimension(1:6,-1:im)::ALPHA
double precision,dimension(1:6,0:im)::EIG
double precision,dimension(1:6,0:im)::GG
double precision,dimension(1:6,0:im)::E
double precision,dimension(1:6,0:im)::PHIM
double precision::DELTA,GAMMA,UCONT
 DO J=1,jm-1
  DO I=0,im-1
   RWL=dSQRT(Q(1,I,J)*jj(I,J))/(dSQRT(Q(1,I+1,J)*jj(I+1,J))+dSQRT(Q(1,I,J)*jj(I,J)))
   RWR=dSQRT(Q(1,I+1,J)*jj(I+1,J))/(dSQRT(Q(1,I+1,J)*jj(I+1,J))+dSQRT(Q(1,I,J)*jj(I,J)))
   AKX=0.5d0*(xix(I,J)+xix(I+1,J))
   AKY=0.5d0*(xiy(I,J)+xiy(I+1,J))
   AJACM=0.5d0*(jj(I,J)+jj(I+1,J))

   UM=RWL*U(I,J)+RWR*U(I+1,J)
   VM=RWL*V(I,J)+RWR*V(I+1,J)
   DUM=RWL*Q(4,I,J)/Q(1,I,J)+RWR*Q(4,I+1,J)/Q(1,I+1,J)
   CM=dSQRT(g*(g-1.0d0)*dabs(DUM-0.5d0*(UM**2.0+VM**2.0)))

!   NOMENCLATURE
   PHI=0.5d0*(g-1.0d0)*(UM**2.0+VM**2.0)
   BETA=1.0d0/(2.0d0*CM**2.0)
   AKXT=AKX/dSQRT(AKX**2.0+AKY**2.0)
   AKYT=AKY/dSQRT(AKX**2.0+AKY**2.0)
   THIT=AKXT*UM+AKYT*VM

!   RIGHT ENGEN-VECTORS
   RR(1,1,I)=1.0d0
   RR(1,2,I)=0.0d0
   RR(1,3,I)=1.0d0
   RR(1,4,I)=1.0d0
   RR(2,1,I)=UM
   RR(2,2,I)=AKYT
   RR(2,3,I)=UM+AKXT*CM
   RR(2,4,I)=UM-AKXT*CM
   RR(3,1,I)=VM
   RR(3,2,I)=-AKXT
   RR(3,3,I)=VM+AKYT*CM
   RR(3,4,I)=VM-AKYT*CM
   RR(4,1,I)=PHI/(g-1.0d0)
   RR(4,2,I)=AKYT*UM-AKXT*VM
   RR(4,3,I)=(PHI+CM**2.0)/(g-1.0d0)+CM*THIT
   RR(4,4,I)=(PHI+CM**2.0)/(g-1.0d0)-CM*THIT

!   INVERS OR RIGHT EIGEN-VECTORS
   RI(1,1,I)=1.0d0-PHI/CM**2.0
   RI(1,2,I)=(g-1.0d0)*UM/CM**2.0
   RI(1,3,I)=(g-1.0d0)*VM/CM**2.0
   RI(1,4,I)=-(g-1.0d0)/CM**2.0
   RI(2,1,I)=-AKYT*UM+AKXT*VM
   RI(2,2,I)=AKYT
   RI(2,3,I)=-AKXT
   RI(2,4,I)=0.0d0
   RI(3,1,I)=BETA*(PHI-CM*THIT)
   RI(3,2,I)=BETA*(AKXT*CM-(g-1.0d0)*UM)
   RI(3,3,I)=BETA*(AKYT*CM-(g-1.0d0)*VM)
   RI(3,4,I)=BETA*(g-1.0d0)
   RI(4,1,I)=BETA*(PHI+CM*THIT)
   RI(4,2,I)=-BETA*(AKXT*CM+(g-1.0d0)*UM)
   RI(4,3,I)=-BETA*(AKYT*CM+(g-1.0d0)*VM)
   RI(4,4,I)=BETA*(g-1.0d0)

!   ENGEN-VALUES AT INTERMEDIATE
   EIGM(1,I)=AKX*UM+AKY*VM
   EIGM(2,I)=EIGM(1,I)
   EIGM(3,I)=EIGM(1,I)+CM*dSQRT(AKX**2.0+AKY**2.0)
   EIGM(4,I)=EIGM(1,I)-CM*dSQRT(AKX**2.0+AKY**2.0)
   EIGM(5,I)=EIGM(1,I)
   EIGM(6,I)=EIGM(1,I)
!   ALPHA
   DO Kk=1,6
    D(Kk)=(Q(Kk,I+1,J)*jj(I+1,J)-Q(Kk,I,J)*jj(I,J))/AJACM
   enddo

   DO Kk=1,6
    IF(Kk.LE.4) THEN
    ALPHA(Kk,I)=RI(Kk,1,I)*D(1)+RI(Kk,2,I)*D(2)+RI(Kk,3,I)*D(3)+RI(Kk,4,I)*D(4)
    ELSE
    ALPHA(Kk,I)=D(Kk)
    END IF
   enddo
  enddo

  DO Kk=1,6
   ALPHA(Kk,-1)=ALPHA(Kk,0)
   ALPHA(Kk,im)=ALPHA(Kk,im-1)
  enddo

  DO I=0,im
   DUM=dSQRT(g*R*T(I,J)*(xix(I,J)**2.0+xiy(I,J)**2.0))
   EIG(1,I)=xix(I,J)*U(I,J)+xiy(I,J)*V(I,J)
   EIG(2,I)=EIG(1,I)
   EIG(3,I)=EIG(1,I)+DUM
   EIG(4,I)=EIG(1,I)-DUM
   EIG(5,I)=EIG(1,I)
   EIG(6,I)=EIG(1,I)
  enddo

!   CAL LIMITER FUNCTION G
  DO Kk=1,6
   DO I=0,im
!     MINMOD LIMITER
  !  if(ALPHA(Kk,I)>=0.0)b=1.0
  !  if(ALPHA(Kk,I)<0.0)b=-1.0

         SI=dSIGN(1.0d0,ALPHA(Kk,I))

    GG(Kk,I)=SI*dMAX1(0.0,dMIN1(dabs(ALPHA(Kk,I)),SI*ALPHA(Kk,I-1)))
   enddo
  enddo

  DO Kk=1,6
   DO I=0,im-1
!    GAMMA
    DELTA=1.0E-10
  !   DELTA=AMAX1(0.0,EIGM(Kk,I)-EIG(Kk,I),EIG(Kk,I+1)-EIGM(Kk,I))
    IF(ALPHA(Kk,I).NE.0.0)THEN
     GAMMA=0.5d0*FPSI(EIGM(Kk,I),DELTA)*(GG(Kk,I+1)-GG(Kk,I))/ALPHA(Kk,I)
    ELSE
     GAMMA=0.0
    ENDIF

!    PHI
    PHIM(Kk,I)=0.5d0*FPSI(EIGM(Kk,I),DELTA)*(GG(Kk,I+1)+GG(Kk,I))-FPSI(EIGM(Kk,I)+GAMMA,DELTA)*ALPHA(Kk,I)
   enddo
  enddo


!   CONVECTION COMPORNENTS
  DO I=0,im
   UCONT=xix(I,J)*U(I,J)+xiy(I,J)*V(I,J)
   E(1,I)=Q(1,I,J)*UCONT
   E(2,I)=Q(2,I,J)*UCONT+P(I,J)/jj(i,j)*xix(I,J)
   E(3,I)=Q(3,I,J)*UCONT+P(I,J)/jj(i,j)*xiy(I,J)
   E(4,I)=Q(4,I,J)*UCONT+P(I,J)/jj(i,j)*UCONT
   E(5,I)=Q(5,I,J)*UCONT
   E(6,I)=Q(6,I,J)*UCONT	
  enddo

!   XI-DIRECTION CONVECTION FLUX
  DO  Kk=1,6
   DO  I=0,im-1
    IF(Kk.LE.4) THEN
     EE(Kk,I,J)=0.5*(E(Kk,I)+E(Kk,I+1)+RR(Kk,1,I)*PHIM(1,I)&
         &    +RR(Kk,2,I)*PHIM(2,I)+RR(Kk,3,I)*PHIM(3,I)+RR(Kk,4,I)*PHIM(4,I))
    else
     EE(Kk,I,J)=0.5*(E(Kk,I)+E(Kk,I+1)+PHIM(Kk,I))
    END IF
   enddo
  enddo
 enddo
return
END subroutine TVDxi

subroutine TVDeta

double precision,dimension(1:4,1:4,0:JM)::RRI=0.0d0
double precision,dimension(1:4,1:4,0:JM)::RII=0.0d0
double precision,dimension(1:6,0:JM)::FIGM=0.0d0
double precision,dimension(1:6)::DI=0.0d0
double precision,dimension(1:6,-1:JM)::ALPHAI=0.0d0
double precision,dimension(1:6,0:JM)::EIGI=0.0d0
double precision,dimension(1:6,0:JM)::GGI=0.0d0
double precision,dimension(1:6,0:JM)::PHIMI=0.0d0
double precision,dimension(1:6,0:JM)::F=0.0d0
double precision RWLI,RWRI,AKXI,AKYI,AJACMI,UMI,VMI,DUMI,CMI,PHII,BETAI,AKXTI
double precision AKYTI,THITI,SI,DELTAI,GAMMAI,UCONTI


 do I=1,IM-1
    do J=0,JM-1
 
       RWLI=dSQRT(Q(1,I,J)*JJ(I,J))&
    &      /(dSQRT(Q(1,I,J+1)*JJ(I,J+1))+dSQRT(Q(1,I,J)*JJ(I,J)))
       RWRI=dSQRT(Q(1,I,J+1)*JJ(I,J+1))&
    &      /(dSQRT(Q(1,I,J+1)*JJ(I,J+1))+dSQRT(Q(1,I,J)*JJ(I,J)))
       AKXI=0.5d0*(etax(I,J)+etax(I,J+1))
       AKYI=0.5d0*(etay(I,J)+etay(I,J+1))
       AJACMI=0.5d0*(JJ(I,J)+JJ(I,J+1))
       UMI=RWLI*U(I,J)+RWRI*U(I,J+1)
       VMI=RWLI*V(I,J)+RWRI*V(I,J+1)
       DUMI=RWLI*Q(4,I,J)/Q(1,I,J)+RWRI*Q(4,I,J+1)/Q(1,I,J+1)
       CMI=dSQRT(g*(g-1.0d0)*dabs(DUMI-0.5d0*(UMI**2+VMI**2)))

       !NOMENCLATURE
       PHII=0.5d0*(g-1.0d0)*(UMI**2+VMI**2)
       BETAI=1.0d0/(2.0d0*CMI**2)
       AKXTI=AKXI/dSQRT(AKXI**2+AKYI**2)
       AKYTI=AKYI/dSQRT(AKXI**2+AKYI**2)
       THITI=AKXTI*UMI+AKYTI*VMI

       !RIGHT ENGEN-VECTORS
       RRI(1,1,J)=1.0d0
       RRI(1,2,J)=0.0d0
       RRI(1,3,J)=1.0d0
       RRI(1,4,J)=1.0d0
       RRI(2,1,J)=UMI
       RRI(2,2,J)=AKYTI
       RRI(2,3,J)=UMI+AKXTI*CMI
       RRI(2,4,J)=UMI-AKXTI*CMI
       RRI(3,1,J)=VMI
       RRI(3,2,J)=-AKXTI
       RRI(3,3,J)=VMI+AKYTI*CMI
       RRI(3,4,J)=VMI-AKYTI*CMI
       RRI(4,1,J)=PHII/(g-1.0d0)
       RRI(4,2,J)=AKYTI*UMI-AKXTI*VMI
       RRI(4,3,J)=(PHII+CMI**2)/(g-1.0d0)+CMI*THITI
       RRI(4,4,J)=(PHII+CMI**2)/(g-1.0d0)-CMI*THITI

       !INVERS OR RIGHT EIGEN-VECTORS
       RII(1,1,J)=1.0d0-PHII/CMI**2
       RII(1,2,J)=(g-1.0d0)*UMI/CMI**2
       RII(1,3,J)=(g-1.0d0)*VMI/CMI**2
       RII(1,4,J)=-(g-1.0d0)/CMI**2
       RII(2,1,J)=-AKYTI*UMI+AKXTI*VMI
       RII(2,2,J)=AKYTI
       RII(2,3,J)=-AKXTI
       RII(2,4,J)=0.0d0
       RII(3,1,J)=BETAI*(PHII-CMI*THITI)
       RII(3,2,J)=BETAI*(AKXTI*CMI-(g-1.0d0)*UMI)
       RII(3,3,J)=BETAI*(AKYTI*CMI-(g-1.0d0)*VMI)
       RII(3,4,J)=BETAI*(g-1.0d0)
       RII(4,1,J)=BETAI*(PHII+CMI*THITI)
       RII(4,2,J)=-BETAI*(AKXTI*CMI+(g-1.0d0)*UMI)
       RII(4,3,J)=-BETAI*(AKYTI*CMI+(g-1.0d0)*VMI)
       RII(4,4,J)=BETAI*(g-1.0d0)

       !ENGEN-VALUES AT INTERMEDIATE
       FIGM(1,J)=AKXI*UMI+AKYI*VMI
       FIGM(2,J)=FIGM(1,J)
       FIGM(3,J)=FIGM(1,J)+CMI*dSQRT(AKXI**2+AKYI**2)
       FIGM(4,J)=FIGM(1,J)-CMI*dSQRT(AKXI**2+AKYI**2)
       FIGM(5,J)=FIGM(1,J)
       FIGM(6,J)=FIGM(1,J)

       !ALPHAI
    do Kk=1,6
       DI(Kk)=(Q(Kk,I,J+1)*JJ(I,J+1)-Q(Kk,I,J)*JJ(I,J))/AJACMI
    end do

    do Kk=1,6
       if(kk <= 4)then
       ALPHAI(kK,J)=RII(Kk,1,J)*DI(1)+RII(Kk,2,J)*DI(2)+RII(Kk,3,J)*DI(3)&
                   &+RII(Kk,4,J)*DI(4)
       else
       ALPHAI(Kk,J)=DI(kk)
       end if
    end do
 end do

 do Kk=1,6
       ALPHAI(Kk,-1)=ALPHAI(Kk,0)
       ALPHAI(kK,JM)=ALPHAI(Kk,JM-1)
 end do

 do J=0,JM
       DUMI=dSQRT(g*R*T(I,J)*(etax(I,J)**2+etay(I,J)**2))
       EIGI(1,J)=etax(I,J)*U(I,J)+etay(I,J)*V(I,J)
       EIGI(2,J)=EIGI(1,J)
       EIGI(3,J)=EIGI(1,J)+DUMI
       EIGI(4,J)=EIGI(1,J)-DUMI
       EIGI(5,J)=EIGI(1,J)
       EIGI(6,J)=EIGI(1,J)
 end do

       !CAL LIMITER FUNCTION G
 do Kk=1,6
    do J=0,JM

       !MINMOD LIMITER
       SI=dSIGN(1.0d0,ALPHAI(Kk,J))
       GGI(Kk,J)=SI*dMAX1(0.0d0,dMIN1(dabs(ALPHAI(Kk,J)),SI*ALPHAI(Kk,J-1)))
    end do
 end do

 do Kk=1,6
    do J=0,JM-1

       !GAMMAI
       DELTAI=1.0E-10
       !DELTAI=dMAX1(0.0d0,FIGM(Kk,J)-EIGI(Kk,J),EIGI(Kk,J+1)-FIGM(Kk,J))
       IF(ALPHAI(Kk,J).NE.0.0d0)THEN
       GAMMAI=0.5d0*FPSI(FIGM(Kk,J),DELTAI)&
    &         *(GGI(Kk,J+1)-GGI(Kk,J))/ALPHAI(Kk,J)
       ELSE
       GAMMAI=0.0d0
       END IF

       !PHII
       PHIMI(Kk,J)=0.5d0*FPSI(FIGM(Kk,J),DELTAI)*(GGI(Kk,J+1)+GGI(Kk,J))&
    &             -FPSI(FIGM(Kk,J)+GAMMAI,DELTAI)*ALPHAI(Kk,J)
    end do
 end do


       !CONVECTION COMPORNENTS
 do J=0,JM
       UCONTI=etax(I,J)*U(I,J)+etay(I,J)*V(I,J)
       F(1,J)=Q(1,I,J)*UCONTI
       F(2,J)=Q(2,I,J)*UCONTI+P(I,J)/JJ(i,j)*etax(I,J)
       F(3,J)=Q(3,I,J)*UCONTI+P(I,J)/JJ(i,j)*etay(I,J)
       F(4,J)=Q(4,I,J)*UCONTI+P(I,J)/JJ(i,j)*UCONTI
       F(5,J)=Q(5,I,J)*UCONTI
       F(6,J)=Q(6,I,J)*UCONTI
 end do

       !XI-DIRECTION CONVECTION FLUX
 do Kk=1,6
    do J=0,JM-1
     if(kk <= 4)then
       FF(KK,I,J)=0.5d0*(F(Kk,J)+F(Kk,J+1)&
    &            +RRI(Kk,1,J)*PHIMI(1,J)+RRI(Kk,2,J)*PHIMI(2,J)&
    &            +RRI(Kk,3,J)*PHIMI(3,J)+RRI(Kk,4,J)*PHIMI(4,J))
     else
       FF(KK,I,J)=0.5d0*(F(Kk,J)+F(Kk,J+1)+PHIMI(Kk,J))
     end if
    end do
 end do

end do
return

end subroutine TVDeta



!*****FPSIä÷êîÇÃíËã`*****           Éµ(z)ÇÃåvéZ
double precision FUNCTION FPSI(Z,DELTA)
   double precision Z,DELTA
      IF(dabs(Z).GE.DELTA)THEN       !.ge. >=
        FPSI=dabs(Z)
      ELSE
        FPSI=0.5d0*(Z**2.0+DELTA**2.0)/DELTA
      END IF
      RETURN
      END function


subroutine runge_kutta

do j=1,jm-1
   do i=1,im-1
      do kk=1,6
QQ(kk,i,j,mm)=QQ(kk,i,j,0)-(1.0d0/dble(5-mm))*dt(i,j)*&
(-EE(kk,i-1,j)+EE(kk,i,j)-FF(kk,i,j-1)+FF(kk,i,j)+(rv(kk,i-1,j)-rv(kk,i,j))+&
(sv(kk,i,j-1)-sv(kk,i,j))-pp(kk,i,j))
q(kk,i,j)=qq(kk,i,j,mm)
      end do
   end do
end do


end subroutine runge_kutta

subroutine calculate    !èîó ÇÃåvéZ
un=u
vn=v

do j=1,jm-1
   do i=1,im-1
      rho(i,j)=jj(i,j)*q(1,i,j)
      u(i,j)=q(2,i,j)/q(1,i,j)
      v(i,j)=q(3,i,j)/q(1,i,j)
      energy(i,j)=jj(i,j)*q(4,i,j)
      k(i,j)=q(5,i,j)/q(1,i,j)
      ipu(i,j)=q(6,i,j)/q(1,i,j)
      eb(i,j)=energy(i,j)/rho(i,j)-(((u(i,j)**2.0)+v(i,j)**2.0)/2.0d0)
      t(i,j)=eb(i,j)*(g-1.0d0)/r
      p(i,j)=(g-1.0d0)*rho(i,j)*eb(i,j)
   end do
end do

!ó¨èo(âE)
do j=0,jm
 slope=(x(im,j)-x(im-2,j))/(x(im-2,j)-x(im-1,j))
   rho(im,j)=rho(im-2,j)+slope*(rho(im-2,j)-rho(im-1,j))
   p(im,j)=p(im-2,j)+slope*(p(im-2,j)-p(im-1,j))
   u(im,j)=u(im-2,j)+slope*(u(im-2,j)-u(im-1,j))
   v(im,j)=v(im-2,j)+slope*(v(im-2,j)-v(im-1,j))
   k(im,j)=k(im-2,j)+slope*(k(im-2,j)-k(im-1,j))
   ipu(im,j)=ipu(im-2,j)+slope*(ipu(im-2,j)-ipu(im-1,j))
   t(im,j)=t(im-2,j)+slope*(t(im-2,j)-t(im-1,j))
   energy(im,j)=energy(im-2,j)+slope*(energy(im-2,j)-energy(im-1,j))

   q(1,im,j)=rho(im,j)/jj(im,j)
   q(2,im,j)=rho(im,j)*u(im,j)/jj(im,j)
   q(3,im,j)=rho(im,j)*v(im,j)/jj(im,j)
   q(4,im,j)=energy(im,j)/jj(im,j)
   q(5,im,j)=rho(im,j)*k(im,j)/jj(im,j)
   q(6,im,j)=rho(im,j)*ipu(im,j)/jj(im,j)
end do

!è„
do i=0,im
   rho(i,jm)=rho(i,jm-1)
   u(i,jm)=u(i,jm-1)
   v(i,jm)=v(i,jm-1)
   p(i,jm)=p(i,jm-1)
   t(i,jm)=t(i,jm-1)
   k(i,jm)=k(i,jm-1)
   ipu(i,jm)=ipu(i,jm-1)
   t(i,jm)=t(i,jm-1)
energy(i,jm)=energy(i,jm-1)

   q(1,i,jm)=rho(i,jm)/jj(i,jm)
   q(2,i,jm)=rho(i,jm)*u(i,jm)/jj(i,jm)
   q(3,i,jm)=rho(i,jm)*v(i,jm)/jj(i,jm)
   q(4,i,jm)=energy(i,jm)/jj(i,jm)
   q(5,i,jm)=rho(i,jm)*k(i,jm)/jj(i,jm)
   q(6,i,jm)=rho(i,jm)*ipu(i,jm)/jj(i,jm)
end do

do i=0,im
!ï«(â∫)
   if(i<6)then
   rho(i,0)=rho(i,1)
   p(i,0)=p(i,1)
   u(i,0)=u(i,1)
   v(i,0)=0.0d0
   k(i,0)=k(i,1)
   ipu(i,0)=ipu(i,1)
   else
   rho(i,0)=rho(i,1)
   p(i,0)=p(i,1)
   u(i,0)=0.0d0
   v(i,0)=0.0d0
   k(i,1)=kwall(i,1)
   k(i,0)=k(i,1)
   ipu(i,1)=ipuwall(i,1)
   ipu(i,0)=ipu(i,1)
   end if
   
   t(i,0)=t(i,1)
   energy(i,0)=p(i,0)/(g-1.0d0)+rho(i,0)*((u(i,0)**2.0+v(i,0)**2.0))/2.0d0

   q(1,i,0)=rho(i,0)/jj(i,0)
   q(2,i,0)=rho(i,0)*u(i,0)/jj(i,0)
   q(3,i,0)=rho(i,0)*v(i,0)/jj(i,0)
   q(4,i,0)=energy(i,0)/jj(i,0)
   q(5,i,0)=rho(i,0)*k(i,0)/jj(i,0)
   q(6,i,0)=rho(i,0)*ipu(i,0)/jj(i,0)
end do



dd=0.0d0
do j=1,jm
   do i=0,im
      veps=(abs(u(i,j)-un(i,j)))/abs(un(i,j))
      if(dd<veps)then
         dd=veps
      end if
      veps=(abs(v(i,j)-vn(i,j)))/abs(vn(i,j))
      if(dd<veps)then
         dd=veps
      end if
   end do
end do

end subroutine calculate

subroutine ltsm    !ã«èäéûä‘çèÇ›ñ@(local time step method)
double precision,parameter::cn=0.2             !ÉNÅ[ÉâÉìêî
double precision::lcxi,lceta,ldxi,ldeta,lipuxi,lipueta,lxi,leta,uu,vv

do j=0,jm
   do i=0,im
!ëŒó¨çÄ
!mu(i,j)=sland(t(i,j))
uu=xix(i,j)*u(i,j)+xiy(i,j)*v(i,j)
vv=etax(i,j)*u(i,j)+etay(i,j)*v(i,j)
c(i,j)=dsqrt(g*r*t(i,j))

lcxi=dabs(uu)+c(i,j)*dsqrt(xix(i,j)**2.0+xiy(i,j)**2.0)
lceta=dabs(vv)+c(i,j)*dsqrt(etax(i,j)**2.0+etay(i,j)**2.0)

!ägéUçÄ
ldxi=2.0d0*(mu(i,j)+mut(i,j))*(xix(i,j)**2.0+xiy(i,j)**2.0)/rho(i,j)
ldeta=2.0d0*(mu(i,j)+mut(i,j))*(etax(i,j)**2.0+etay(i,j)**2.0)/rho(i,j)

!óêÇÍÇÃéUàÌçÄ
lipuxi=ipu(i,j)/k(i,j)
lipueta=lipuxi
!ã«èäéûä‘çèÇ›ïù
lxi=lcxi+ldxi+lipuxi
leta=lceta+ldeta+lipueta
dt(i,j)=cn*dmin1(1.0d0/lxi,1.0d0/leta)

   end do
end do
end subroutine ltsm

subroutine yuplus
double precision tauwall,utau,yp,kp,ipup,ypp
double precision,parameter::kappar=0.42d0
double precision,parameter::largee=9.7d0
double precision,dimension(0:im,0:jm)::yplus=0.0d0,uplus=0.0d0


do j=0,jm
   do i=5,im
tauwall=mu(i,0)*(u(i,1)-u(i,0))/(y(i,1)-y(i,0))
utau=dsqrt(tauwall/rho(i,0))
uplus(i,j)=u(i,j)/utau
yplus(i,j)=utau*y(i,j)/nu(i,j)
   end do
end do

open(19,file='ypup.dat',status='replace')

do j=0,jm
      write(19,'(f13.6,$)') yplus(45,j)
      write(19,'(f13.6)') uplus(45,j)
   end do
close(19)

end subroutine yuplus

subroutine output
open(17,file='kadai8.dat',status='replace')
   do j=0,jm
      do i=0,im
         write(17,'(f13.6,$)') x(i,j)
         write(17,'(f13.6,$)') y(i,j)
         write(17,'(f13.6,$)') u(i,j)
         write(17,'(f13.6,$)') v(i,j)
         write(17,'(f13.3,$)') p(i,j)
         write(17,'(f13.6)') t(i,j)
      end do
   end do
close(17)

open(18,file='kadai8.fld',status='replace')
write(18,'(a)')'# AVS field file'
write(18,'(a)')'ndim=2'
write(18,'(a,i0)')'dim1=',im+1
write(18,'(a,i0)')'dim2=',jm+1
write(18,'(a)')'nspace=2'
write(18,'(a)')'veclen=4'
write(18,'(a)')'data=float'
write(18,'(a)')'field=irregular'
write(18,'(a)')'label=u,v,p,t'
write(18,'(a)')'variable 1 file=kadai8.dat filetype=ascii skip=0 offset=2 stride=6'
write(18,'(a)')'variable 2 file=kadai8.dat filetype=ascii skip=0 offset=3 stride=6'
write(18,'(a)')'variable 3 file=kadai8.dat filetype=ascii skip=0 offset=4 stride=6'
write(18,'(a)')'variable 4 file=kadai8.dat filetype=ascii skip=0 offset=5 stride=6'
write(18,'(a)')'coord 1 file=kadai8.dat filetype=ascii skip=0 offset=0 stride=6'
write(18,'(a)')'coord 2 file=kadai8.dat filetype=ascii skip=0 offset=1 stride=6'
close(18)
end subroutine output

end program kadai8
