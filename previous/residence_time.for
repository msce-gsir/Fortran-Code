	program  dipole distribution
	integer frame,toto,totfra,totc
	integer i,k,cn
	parameter (totfra=10001,totgri=1332,delt=0.1,vale=0.37,
     $	resh2o=7823,totc=6528)
	
	integer bin,hist(180),thita(resh2o),fi(resh2o)
	integer flag(resh2o)
	double precision sum,percent(180),Nax,Nay,Naz,Rt(totfra)
      double precision xo(resh2o),yo(resh2o),zo(resh2o)
	double precision rxh2o1(totfra),ryh2o1(totfra),rzh2o1(totfra)
	double precision rxc(totc),ryc(totc),rzc(totc)
	double precision resit,t
	double precision om,omx,omy,omz,ov,voz,cosmov,mov,dipole

	double precision oc1,oc2,diameter,diameter0
	parameter (diameter=1.085564036)		
	double precision c1c2,c1c2sq,oc1sq,oc2sq
	double precision oc1x,oc1y,oc1z,oc2x,oc2y,oc2z
	double precision sm,soc1c2,sh,cosoc1c2,cosoc2c1
	double precision pi
	parameter (pi=3.1415927)
	character*80 line,line1	
	
	
	
	double precision c1x,c1y,c1z,c2x,c2y,c2z
	double precision c1x0,c1y0,c1z0,c2x0,c2y0,c2z0
C	parameter (c1x0=0.0,c1y0=0.0,c1z0=6.74025)
C	parameter (c2x0=0.0,c2y0=0.0,c2z0=-6.74025)
C	parameter (c1c2sq0=c1c20**2)

	open(unit=10,file='cntsol33.gro')
	open(unit=11,file='residence_time.txt')

c	loop over all frame

	do 12 i=1,resh2o
	flag(i)=0
12	continue

	do 13 frame=1,totfra
	read(10,*)
      read(10,*)

	t=(frame-1)*delt

	do 10 i=1,totgri
	read(10,*)
10	continue

	do 11 i=1,totc
	read(10,203)line,rxc(i),ryc(i),rzc(i),line1
11    continue

	do 20 i=1,resh2o
      read(10,203)line,xo(i),yo(i),zo(i),line1
	read(10,*)
	read(10,*)
20	continue

	read(10,203)line,Nax,Nay,Naz,line1

	if(frame.eq.1)then
	Rt(frame)=1
 	do 22 i=1,resh2o
	dist=(Nax-xo(i))*(Nax-xo(i))+(Nay-yo(i))*(Nay-yo(i))+
     $	(Naz-zo(i))*(Naz-zo(i))
	if (dist.le.(vale*vale)) then
	fi(i)=1
	flag(i)=flag(i)+1
	else fi=0
	endif
22	continue
	endif

	if(.not.(frame.eq.1))then
	Rt(frame)=0
	cn=0
	do 21 i=1,resh2o
	dist=(Nax-xo(i))*(Nax-xo(i))+(Nay-yo(i))*(Nay-yo(i))+
     $	(Naz-zo(i))*(Naz-zo(i))
	if (dist.le.(vale*vale)) then
	flag(i)=flag(i)+1
	cn=cn+1
	thita(i)=1
	if (flag(i).eq.frame) then
	Rt(frame)=Rt(frame)+fi(i)*thita(i)
	endif
	else 
	thita(i)=0
	endif
21	continue
	Rt(frame)=Rt(frame)/cn
	endif

	read(10,*)

	write(11,1004)t,Rt(frame)

13	continue

	resit=0
	do 14 frame=1,totfra
	resit=resit+0.1*Rt(frame)
14	continue

	write(11,1003)'resit=    ',resit

201	format(A36,f8.3,A24)
202   format(A10,I6)
203	format(A20,3f8.3,A24)
1004  format(f9.3,f8.4)
1003  format(A10,f20.3)
1002  format(f10.3,f8.3)
	 
	close(11)
	close(10)
	end
