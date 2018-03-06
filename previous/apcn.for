	program  hydrogen bonds inside tube
	integer frame,toto,totfra,totc,nacar,kcar,r
	integer i,j,k,totcar,totion,step,binz,cn
	double precision rnamin,deltz,rkmin,zmin,zmax
	parameter (totfra=100001,toto=2310,totc=255,totgri=1284,
     $	rnamin=0.32,totcar=5,totion=30,deltz=0.01,rkmin=0.34,step=600,
     $	zmin=2.262,zmax=3.738)	
	double precision histnacn(step),histnacar(step),histnan(step)
	double precision histkcn(step),histkcar(step),histkn(step)
	double precision hena(step),heana(step),henacar(step)
	double precision perhena(step)
	double precision hek(step),heak(step),hekcar(step)
	double precision perhek(step)
	double precision xcar(totcar),ycar(totcar),zcar(totcar)
	double precision xna(totion),yna(totion),zna(totion)
	double precision xk(totion),yk(totion),zk(totion)
	double precision perznacn(step),perznacar(step)
	double precision perzkcn(step),perzkcar(step)
	double precision rxo(toto),ryo(toto),rzo(toto)
  double precision xh1(toto),yh1(toto),zh1(toto)
  double precision xh2(toto),yh2(toto),zh2(toto)
	double precision rxc(totc),ryc(totc),rzc(totc)
	double precision dist,distcar,difnacar,difkcar
	double precision difhenacar,difhekcar
	double precision factor,epsionOW,sigmaOW,cow
	double precision epsionocar,sigmaocar,cocar,rh1,rh2
	double precision epsionHW,sigmaHW,chw
	double precision epsionna,sigmana,cna,epsionk,sigmak,ck
	double precision epsona,sigmaona,epsok,sigmaok
	double precision epsnacar,sigmanacar,epskcar,sigmakcar

	character*80 line


	do 111 i=1,step
	histnan(i)=0.0
	histnacn(i)=0.0
	histnacar(i)=0.0
	perznacn(i)=0.0
	perznacar(i)=0.0
	histkn(i)=0.0
	histkcn(i)=0.0
	histkcar(i)=0.0
	perzkcn(i)=0.0
	perzkcar(i)=0.0
	hena(i)=0.0
	hek(i)=0.0
	heana(i)=0.0
	heak(i)=0.0
	henacar(i)=0.0
	hekcar(i)=0.0
	perhena(i)=0.0
	perhek(i)=0.0
111	continue



	open(unit=10,file='cntsol33.gro',status='unknown')
	open(unit=11,file='ap_of_coordination_number.txt',
     &	status='unknown')

	difnacar=0.0
	difkcar=0.0
	difhenacar=0.0
	difhekcar=0.0
	nacar=0
	kcar=0

	factor=138.935

	epsionocar=8.78640e-01
	sigmaocar=2.96000e-01
	cocar=-0.37

	epsionOW=6.50194e-01
	sigmaOW=3.16557e-01
	cow=-0.8476

	epsionHW=0.00000e+00
	sigmaHW=0.00000e+00
	chw=0.4238

	epsionna=1.15980e-2
	sigmana=3.33045e-01
	cna=1.00

	epsona=4*sqrt(epsionOW*epsionna)
	sigmaona=(sigmaOW+sigmana)/2
	sigmaona=sigmaona**6

	epsnacar=4*sqrt(epsionocar*epsionna)
	sigmanacar=(sigmaocar+sigmana)/2
	sigmanacar=sigmanacar**6

	epsionk=1.37235e-03
	sigmak=4.93463e-01
	ck=1.00

	epsok=4*sqrt(epsionOW*epsionk)
	sigmaok=(sigmaOW+sigmak)/2
	sigmaok=sigmaok**6

	epskcar=4*sqrt(epsionocar*epsionk)
	sigmakcar=(sigmaocar+sigmak)/2
	sigmakcar=sigmakcar**6

c	loop over all frames
	do 13 frame=1,totfra
	
      read(10,*)
	read(10,*)

	do 10 i=1,totgri
      read(10,*)
10	continue

	do 900 i=1,totc
	read(10,*)
900   continue

	do 11 i=1,totcar
      read(10,*)
	read(10,201)line,xcar(i),ycar(i),zcar(i),line
11	continue

      do 101 i=1,toto
	read(10,201)line,rxo(i),ryo(i),rzo(i),line
	read(10,201)line,xh1(i),yh1(i),zh1(i),line
      read(10,201)line,xh2(i),yh2(i),zh2(i),line
101   continue

	do 103 i=1,totion
	read(10,201)line,xk(i),yk(i),zk(i),line
103	continue

	do 102 i=1,totion
	read(10,201)line,xna(i),yna(i),zna(i),line
102	continue

	do 104 i=1,2*totion
      read(10,*)
104	continue

      read(10,*)

	do 12 i=1,totion
	binz=int(zna(i)/deltz)+1.0
	histnan(binz)=histnan(binz)+1.0
	cn=0
	hena(binz)=0.0
	do 15 j=1,toto
	dist=(xna(i)-rxo(j))*(xna(i)-rxo(j))+
     $	(yna(i)-ryo(j))*(yna(i)-ryo(j))+
     $	(zna(i)-rzo(j))*(zna(i)-rzo(j))
	if (dist.le.(rnamin*rnamin)) then
	cn=cn+1
	histnacn(binz)=histnacn(binz)+1.0
	r=sigmaona/(dist**3)
	hena(binz)=hena(binz)+epsona*(r**2-r)
	hena(binz)=hena(binz)+factor*cow*cna/sqrt(dist)
	rh1=(xna(i)-xh1(j))*(xna(i)-xh1(j))+
     $	(yna(i)-yh1(j))*(yna(i)-yh1(j))+
     $	(zna(i)-zh1(j))*(zna(i)-zh1(j))
	rh2=(xna(i)-xh2(j))*(xna(i)-xh2(j))+
     $	(yna(i)-yh2(j))*(yna(i)-yh2(j))+
     $	(zna(i)-zh2(j))*(zna(i)-zh2(j))
	hena(binz)=hena(binz)+factor*chw*cna/sqrt(rh1)+
     $	factor*chw*cna/sqrt(rh2)
	endif
15	continue
	if (cn.NE.0) then
	hena(binz)=hena(binz)/cn
	heana(binz)=heana(binz)+hena(binz)
	endif
	if (zna(i).ge.zmin.and.zna(i).le.zmax) then
	do 16 k=1,totcar
	distcar=(xna(i)-xcar(k))*(xna(i)-xcar(k))+
     $	(yna(i)-ycar(k))*(yna(i)-ycar(k))+
     $	(zna(i)-zcar(k))*(zna(i)-zcar(k))
	if (distcar.le.(rnamin*rnamin)) then
	histnacar(binz)=histnacar(binz)+1.0
	r=sigmanacar/(distcar**3)
	henacar(binz)=henacar(binz)+epsnacar*(r**2-r)
	henacar(binz)=henacar(binz)+factor*cocar*cna/sqrt(distcar)
	endif
16	continue
	endif
12	continue

	do 22 i=1,totion
	binz=int(zk(i)/deltz)+1.0
	histkn(binz)=histkn(binz)+1.0
	cn=0
	hek(binz)=0.0
	do 25 j=1,toto
	dist=(xk(i)-rxo(j))*(xk(i)-rxo(j))+
     $	(yk(i)-ryo(j))*(yk(i)-ryo(j))+
     $	(zk(i)-rzo(j))*(zk(i)-rzo(j))
	if (dist.le.(rkmin*rkmin)) then
	cn=cn+1
	histkcn(binz)=histkcn(binz)+1.0
	r=sigmaok/(dist**3)
	hek(binz)=hek(binz)+epsok*(r**2-r)
	hek(binz)=hek(binz)+factor*cow*ck/sqrt(dist)
	rh1=(xk(i)-xh1(j))*(xk(i)-xh1(j))+
     $	(yk(i)-yh1(j))*(yk(i)-yh1(j))+
     $	(zk(i)-zh1(j))*(zk(i)-zh1(j))
	rh2=(xk(i)-xh2(j))*(xk(i)-xh2(j))+
     $	(yk(i)-yh2(j))*(yk(i)-yh2(j))+
     $	(zk(i)-zh2(j))*(zk(i)-zh2(j))
	hek(binz)=hek(binz)+factor*chw*ck/sqrt(rh1)+
     $	factor*chw*ck/sqrt(rh2)
	endif
25	continue
	if (cn.NE.0) then
	hek(binz)=hek(binz)/cn
	heak(binz)=heak(binz)+hek(binz)
	endif
	if (zk(i).ge.zmin.and.zk(i).le.zmax) then
	do 26 k=1,totcar
	distcar=(xk(i)-xcar(k))*(xk(i)-xcar(k))+
     $	(yk(i)-ycar(k))*(yk(i)-ycar(k))+
     $	(zk(i)-zcar(k))*(zk(i)-zcar(k))
	if (distcar.le.(rkmin*rkmin)) then
	histkcar(binz)=histkcar(binz)+1.0
	r=sigmakcar/(distcar**3)
	hekcar(binz)=hekcar(binz)+epskcar*(r**2-r)
	hekcar(binz)=hekcar(binz)+factor*cocar*ck/sqrt(distcar)
	endif
26	continue
	endif
22	continue

13	continue

	do 14 binz=1,step
	if (histnan(binz).NE.0) then
	perznacn(binz)=histnacn(binz)/histnan(binz)
	perhena(binz)=heana(binz)/histnan(binz)
	perznacar(binz)=(histnacn(binz)+histnacar(binz))/histnan(binz)
	if (histnacar(binz).NE.0) then
	difhenacar=difhenacar+henacar(binz)/histnan(binz)
	difnacar=difnacar+histnacar(binz)/histnan(binz)
	nacar=nacar+1
	endif
	endif
	write(11,202)binz*deltz,perznacn(binz),perznacar(binz),
     $	perhena(binz)
14	continue

	write(11,*)
	write(11,*)

	do 24 binz=1,step
	if (histkn(binz).NE.0) then
	perzkcn(binz)=histkcn(binz)/histkn(binz)
	perhek(binz)=heak(binz)/histkn(binz)
	perzkcar(binz)=(histkcn(binz)+histkcar(binz))/histkn(binz)
	if (histkcar(binz).NE.0) then
	difhekcar=difhekcar+hekcar(binz)/histkn(binz)
	difkcar=difkcar+histkcar(binz)/histkn(binz)
	kcar=kcar+1
	endif
	endif
	write(11,202)binz*deltz,perzkcn(binz),perzkcar(binz),
     $	perhek(binz)
24	continue

	write(11,203)'difnacar=   ',difnacar/nacar
	write(11,203)'difkcar=    ',difkcar/kcar
	write(11,203)'difhenacar= ',difhenacar/nacar
	write(11,203)'difhekcar=  ',difhekcar/kcar

201   format(A20,3f8.3,A24)
202   format(4f12.3)
203	  format(A12,f12.8)
	 
	close(10)
	close(11)
	end
