	program rdf
	integer toto,nin,totfra,totc
	integer i,k,bin,j,binh2o
	parameter (totfra=5001,toto=7823,totc=6528,totgri=1332)
	integer maxbin,maxh2o,index,ctype,otype,r1,r2,r3
	parameter (maxbin=300,maxh2o=1000) 
	double precision hist(maxbin)
	double precision histh2o(0:maxh2o),sumh2o,nic
	double precision denpro(maxbin),rod(maxbin)	
	double precision oc1l,c1c20,oc1,oc2,c1c2,diameter,t_delta
C    t_delta=0.02 ps
	parameter (diameter=1.085564036,t_delta=0.1)
	double precision xmax,xmin,ymax,ymin,zmax,zmin
	double precision oxmax,oxmin,oymax,oymin,ozmax,ozmin	
	double precision rxc(totc),ryc(totc),rzc(totc)
	double precision rxo(toto),ryo(toto),rzo(toto)
	double precision boxa,boxb,boxc
	double precision c1c2sq0,c1c2sq,oc1sq,oc2sq
	double precision oc1x,oc1y,oc1z,oc2x,oc2y,oc2z
	double precision sm,soc1c2,sh,cosoc1c2,cosoc2c1
	double precision maxr,delr,rlower,rupper,nideal
	double precision bxmax,bxmin,bymax,bymin,bzmax,bzmin
	double precision vol,dens,const,pi
	parameter (pi=3.1415927)
	double precision c1x,c1y,c1z,c2x,c2y,c2z
	character*30 line	


	data denpro/maxbin*0.0/
	data hist/maxbin*0.0/


	do 33 binh2o=0,maxh2o
	histh2o(binh2o)=0.0
33    continue
	
	open(unit=10,file='cntsol33.gro',status='unknown')
	open(unit=11,file='number_in_CNT.txt',status='unknown')

	do 13 j=1,totfra
		
      read(10,*)
	read(10,*)

	do 10 i=1,totgri
      read(10,*)
10	continue

	do 900 i=1,totc
	read(10,201)line,rxc(i),ryc(i),rzc(i),line
900   continue

      xmax=rxc(1)
	xmin=rxc(1)
	ymax=ryc(1)
	ymin=ryc(1)
	zmax=rzc(1)
	zmin=rzc(1)

      
      do 908 i=2,totc
	if(rxc(i).gt.xmax) xmax=rxc(i)
		   
	if(rxc(i).lt.xmin) xmin=rxc(i) 
	
	if(ryc(i).gt.ymax) ymax=ryc(i)
	   
	if(ryc(i).lt.ymin) ymin=ryc(i)
	
	if(rzc(i).gt.zmax) zmax=rzc(i)
	   
	if(rzc(i).lt.zmin) zmin=rzc(i)
	
908   continue

      c1x=(xmax+xmin)/2.0
	c1y=(ymax+ymin)/2.0
	c1z=zmin

	c2x=(xmax+xmin)/2.0
      c2y=(ymax+ymin)/2.0
	c2z=zmax

      do 101 i=1,toto
	read(10,201)line,rxo(i),ryo(i),rzo(i),line
	read(10,*)
      read(10,*)
101   continue

	c1c2=c2z-c1z
	c1c2sq=c1c2**2
      
      nin=0


	do 14 k=1, toto

      oc1x=rxo(k)-c1x
	oc1y=ryo(k)-c1y
	oc1z=rzo(k)-c1z
      oc1sq=oc1x**2+oc1y**2+oc1z**2
	oc1=sqrt(oc1sq)

	oc2x=rxo(k)-c2x
	oc2y=ryo(k)-c2y
	oc2z=rzo(k)-c2z
      oc2sq=oc2x**2+oc2y**2+oc2z**2
	oc2=sqrt(oc2sq)

      cosoc1c2=(oc1sq+c1c2sq-oc2sq)/(2.0*oc1*c1c2)
	cosoc2c1=(oc2sq+c1c2sq-oc1sq)/(2.0*oc2*c1c2) 

      sm=(oc1+oc2+c1c2)/2.0
	soc1c2=sqrt(abs(sm*(sm-c1c2)*(sm-oc1)*(sm-oc2)))
	sh=2.0*soc1c2/(c1c2)
	
      if((cosoc1c2.gt.0.0).and.(cosoc1c2.lt.1.0)) then
	if((cosoc2c1.gt.0.0).and.(cosoc2c1.lt.1.0)) then

      if(sh.lt.diameter/2.0) then 
	nin=nin+1
	endif 
	endif
	endif 
14    continue 	

	histh2o(nin)=histh2o(nin)+1.0
	read(10,*)
	read(10,*)
13	continue
	

	sumh2o=0.0
	do 30 binh2o=0,maxh2o
	sumh2o=sumh2o+histh2o(binh2o)
30    continue

      do 31 binh2o=0,maxh2o
	histh2o(binh2o)=histh2o(binh2o)/sumh2o*100
	write(11,'(i6,f12.6)')binh2o,histh2o(binh2o)
31	continue

	write(11,*)

	nic=0.0
	do 32 binh2o=0,maxh2o
	nic=nic+binh2o*histh2o(binh2o)
32	continue

	write(11,1002)'nic=      ',nic
201   format(A20,3f8.3,A24)
1002  format(A10,f10.3)
	 
      close(10)
	close(11)
	end
