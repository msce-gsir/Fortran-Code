	program Dipole_profile_friction
C**********************************************************************C
C                                                                      C
C       Dipole_profile:For friction case                           C
C                                                                      C
C                                                                      C
C       modified by zhuyudan 12/27/2012
C                                                                      C
C----------------------------------------------------------------------C
	integer toto,nin,totfra
	integer i,k,bin,j,binh2o   

	parameter (totfra=10000,toto=493)

	integer maxbin,index
	parameter (maxbin=100) 
	double precision histo(maxbin),histdi(maxbin),histan(maxbin)
	double precision denpdi(maxbin),denpan(maxbin),rod(maxbin)	
	double precision oc1l,c1c20,oc1,oc2,c1c2,diameter,t_delta
	double precision xmax,xmin,ymax,ymin,zmax,zmin
	parameter (ymax=28.8,ymin=0,xmax=59.4,xmin=20.6,zmax=40.08,
     &zmin=23.92)
	double precision oxmax,oxmin,oymax,oymin,ozmax,ozmin	

	double precision x_o(toto),y_o(toto),z_o(toto)
	double precision arrdip(toto,totfra),arrang(toto,totfra)
	integer arrind(toto,totfra)	
	double precision absdip(maxbin),absang(maxbin),absind(maxbin)
	double precision x_h1(toto),y_h1(toto),z_h1(toto)
	double precision x_h2(toto),y_h2(toto),z_h2(toto)
	double precision boxa,boxb,boxc
	double precision oz,h1z,h2z
	double precision maxz,delz,nideal
	double precision vol,dens,const,pi
	double precision c1x,c1y,c1z,c2x,c2y,c2z
	character*2 ohead	

	double precision middlex,middley,middlez
	double precision verticalz
	double precision om,omx,omy,omz,ov,voz,cosmov,mov,dipole,angle


	data denpdi/maxbin*0.0/
	data denpan/maxbin*0.0/
	data histo/maxbin*0.0/
	data histdi/maxbin*0.0/

	parameter (pi=3.1415927)	
	
	do 822 k=1,maxbin
	absdip(k)=0
	absang(k)=0
822	continue


	
	open(unit=10,file='water_conf_5.arc',status='unknown')
	open(unit=11,file='sd_Dipole_5.txt',status='unknown')


	do 13 j=1,totfra

	read(10,*) num

	boxa=xmax-xmin
	boxb=ymax-ymin
	boxc=zmax-zmin
	maxz=boxc
	delz=maxz/maxbin
	nin=0



      do 101 i=1,toto
		read(10,200)index,ohead,x_h1(i),y_h1(i),z_h1(i)
		read(10,200)index,ohead,x_o(i),y_o(i),z_o(i)
		read(10,200)index,ohead,x_h2(i),y_h2(i),z_h2(i)
101   continue


200	FORMAT (I5,2x,A2,3F10.5)

      


	do 14 k=1, toto
c	if(abs(x_o(k)).lt.12.and.abs(y_o(k)).lt.10)then	
! no water near the Wall; 
	middlex=(x_h1(k)+x_h2(k))/2
	middley=(y_h1(k)+y_h2(k))/2
	middlez=(z_h1(k)+z_h2(k))/2

	verticalz=middlez

	omx=x_o(k)-middlex
	omy=y_o(k)-middley
	omz=z_o(k)-middlez
	om=sqrt(omx**2+omy**2+omz**2)

	voz=verticalz-z_o(k)
	ov=abs(voz)
	cosmov=voz/om
	mov=acos(cosmov)/pi*180.0
	dipole=cosmov
	angle=mov

	oz=z_o(k)-zmin
	h1z=z_h1(k)-zmin
	h2z=z_h2(k)-zmin
	bin=int(oz/delz)+1
	arrdip(k,j)=dipole
	arrang(k,j)=angle
	arrind(k,j)=bin
	histo(bin)=histo(bin)+1.0
	histdi(bin)=histdi(bin)+dipole
	histan(bin)=histan(bin)+angle
c	endif
14	continue 	
13	continue

	write(11,'(A3,1x,3A12)')"NUM","Z-Direction","Dipole","Angle"

	do 21 bin=1,maxbin
	rod(bin)=bin*delz/10.0
	denpdi(bin)=histdi(bin)/histo(bin)
	denpan(bin)=histan(bin)/histo(bin)
	write(11,'(I3,3f12.6)')bin,rod(bin),denpdi(bin),denpan(bin)
21	continue
	write(11,*)
	
	do 24 j=1,totfra
	do 22 k=1,toto
	absdip(arrind(k,j))=(arrdip(k,j)-denpdi(arrind(k,j)))**2
     $	+absdip(arrind(k,j))
	absang(arrind(k,j))=(arrang(k,j)-denpan(arrind(k,j)))**2
     $	+absang(arrind(k,j))
22	continue
24	continue

	do 23 bin=1,maxbin
	write(11,'(I3,3f12.6)')bin,rod(bin),
     $	absdip(bin)/histo(bin),absang(bin)/histo(bin)
23	continue



	close(10)
	close(11)
	end
