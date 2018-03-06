	program energy
c	this code calculate the non-bond interactions between molecules
c	they were divided into the CNT-sub and Sub-Sub one.
	integer nmol,nframe,nc
	parameter (nmol=1000,nframe=1001, nc=2560)
	double precision xo,yo,zo
	double precision xh1, yh1, zh1
	double precision xh2, yh2, zh2
	double precision factor
	integer index, index1, frame, n1
	double precision xc(nc), yc(nc), zc(nc)
	character*5 head
	double precision xo1(nmol,nframe),yo1(nmol,nframe),
     &	zo1(nmol,nframe)
	double precision xh11(nmol,nframe),yh11(nmol,nframe),
     &	zh11(nmol,nframe)
	double precision xh21(nmol,nframe),yh21(nmol,nframe),
     &	zh21(nmol,nframe)
	integer flag(nmol,nframe)
c
	double precision disxOWOW(nmol,nmol),disyOWOW(nmol,nmol),
     &	diszOWOW(nmol,nframe), disOWOW(nmol,nframe)
	double precision disxOWHW1(nmol,nmol),disyOWHW1(nmol,nmol),
     &	diszOWHW1(nmol,nmol), disOWHW1(nmol,nmol)
	double precision disxOWHW2(nmol,nmol),disyOWHW2(nmol,nmol),
     &	diszOWHW2(nmol,nmol), disOWHW2(nmol,nmol)
c
	double precision disxHW1OW(nmol,nmol),disyHW1OW(nmol,nmol),
     &	diszHW1OW(nmol,nmol), disHW1OW(nmol,nmol)
	double precision disxHW1HW1(nmol,nmol),disyHW1HW1(nmol,nmol),
     &	diszHW1HW1(nmol,nmol), disHW1HW1(nmol,nmol)
	double precision disxHW1HW2(nmol,nmol),disyHW1HW2(nmol,nmol),
     &	diszHW1HW2(nmol,nmol), disHW1HW2(nmol,nmol)
c
	double precision disxHW2OW(nmol,nmol),disyHW2OW(nmol,nmol),
     &	diszHW2OW(nmol,nmol), disHW2OW(nmol,nmol)
	double precision disxHW2HW1(nmol,nmol),disyHW2HW1(nmol,nmol),
     &	diszHW2HW1(nmol,nmol), disHW2HW1(nmol,nmol)
	double precision disxHW2HW2(nmol,nmol),disyHW2HW2(nmol,nmol),
     &	diszHW2HW2(nmol,nmol), disHW2HW2(nmol,nmol)
c
	double precision epsionC, sigmaC, CC
	double precision epsionOW, sigmaOW, COW
	double precision epsionHW, sigmaHW, CHW
	double precision epison1,epison2
	double precision sigma16,sigma26
	double precision nw
c
	double precision ECNTSUBLJ(nframe), ECNTSUBC(nframe)
	double precision ESUBSUBLJ(nframe), ESUBSUBC(nframe)
c
	double precision ECNTSUBLJT, ECNTSUBCT
	double precision ESUBSUBLJT, ESUBSUBCT
c
	open (unit=10, file='cnt.gro', status='unknown')
	open (unit=11, file='SOLinCNT.in.xvg', status='unknown') 
	open (unit=12, file='NumberinCNT.in.xvg', status='unknown')
	open (unit=14, file='EnergyinCNT.txt', status='unknown')
c	read the coordination for nanotube
	read(10,*)
	read (10,*) 
	do i=1,nc
	read (10,*) n1, head, head, n1, xc(i),yc(i),zc(i)
	end do
	
c	read the number of SOL molecules
	do i=1,12
	read (12,*) 
	end do 
	read (12,*) nw
c	several forcefield parameters
	factor=138.935
c	CNT
	epsionC=2.29288e-01
	sigmaC=3.55000e-01
	cc=0.0
c	OW
	epsionOW=6.50194e-01
	sigmaOW=3.16557e-01
	cow=-0.8476
c	HW
	epsionHW=0.00000e+00
	sigmaHW=0.00000e+00
	chw=0.4238
c	read the coordination in frame
	do i=1,12
	read (11,*)
	end do
c     main loop
	do 100 i=1,nw
	read (11,*) frame,index,xo,yo,zo
	read (11,*) frame,index,xh1,yh1,zh1
	read (11,*) frame,index1,xh2,yh2,zh2
c	calculate the dipole 
	if (frame.ne.frame1) then
	j=0
	frame1=frame
	endif
	j=j+1
c	OW
	xo1(j,frame+1)=xo
	yo1(j,frame+1)=yo
	zo1(j,frame+1)=zo
c	HW1
	xh11(j,frame+1)=xh1
	yh11(j,frame+1)=yh1
	zh11(j,frame+1)=zh1
c	HW2
	xh21(j,frame+1)=xh2
	yh21(j,frame+1)=yh2
	zh21(j,frame+1)=zh2
c
	flag(j,frame+1)=1

c			
100	continue  
c	calculate the energy
c	CNT-Sub C-OW
	sigma1=(sigmaC+sigmaOW)/2.0
	epsion1=4*sqrt(epsionC*epsionOW)
	sigma16=sigma1**6
c	Sub-Sub OW-OW
	sigma2=(sigmaOW+sigmaOW)/2.0
	epsion2=4*sqrt(epsionOW*epsionOW)
	sigma26=sigma2**6
c
	do 200 i=1,nframe
c	CNT-Sub
	do 201 j1=1,nc
	do 202 j2=1,nmol
c	SOL-OW
	if (flag(j2,i).eq.1.0) then
	x1=xo1(j2,i)-xc(j1)
	y1=yo1(j2,i)-yc(j1)
	z1=zo1(j2,i)-zc(j1)
	r1=x1*x1+y1*y1+z1*z1
c	L-J-12-6
	r3=sigma16/r1**3
	ECNTSUBLJ(i)=ECNTSUBLJ(i)+epsion1*(r3**2-r3)
	endif
202	continue
201	continue
c	Sub-Sub
c	calculate the distance between every particles
	do 210 j1=1,nmol-1
	do 220 j2=j1+1,nmol
c	write (*,*) flag(j2,i), flag(j1,i)
	if (flag(j2,i).eq.1.and.flag(j1,i).eq.1) then 
c	Ow-Ow
	disxOwOw(j1,j2)=xo1(j1,i)-xo1(j2,i)
	disyOwOw(j1,j2)=yo1(j1,i)-yo1(j2,i)
	diszOwOw(j1,j2)=zo1(j1,i)-zo1(j2,i)
	disOWOW(j1,j2)=disxOwOw(j1,j2)**2
     &	+disyOwOw(j1,j2)**2
     &    +diszOwOw(j1,j2)**2
c	Ow-HW1
	disxOwHw1(j1,j2)=xo1(j1,i)-xh11(j2,i)
	disyOwHw1(j1,j2)=yo1(j1,i)-yh11(j2,i)
	diszOwHw1(j1,j2)=zo1(j1,i)-zh11(j2,i)
	disOWHW1(j1,j2)=disxOwHW1(j1,j2)**2
     &	+disyOwHW1(j1,j2)**2
     &    +diszOwHW1(j1,j2)**2
c	Ow-HW2
	disxOwHw2(j1,j2)=xo1(j1,i)-xh21(j2,i)
	disyOwHw2(j1,j2)=yo1(j1,i)-yh21(j2,i)
	diszOwHw2(j1,j2)=zo1(j1,i)-zh21(j2,i)
	disOwHW2(j1,j2)=disxOwHW2(j1,j2)**2
     &	+disyOwHW2(j1,j2)**2
     &    +diszOwHW2(j1,j2)**2
c	Hw1-OW
	disxHw1OW(j1,j2)=xh11(j1,i)-xo1(j2,i)
	disyHw1OW(j1,j2)=yh11(j1,i)-yo1(j2,i)
	diszHw1OW(j1,j2)=zh11(j1,i)-zo1(j2,i)
	disHW1OW(j1,j2)=disxHW1Ow(j1,j2)**2
     &	+disyHW1Ow(j1,j2)**2
     &    +diszHW1Ow(j1,j2)**2
c	Hw1-HW1
	disxHw1HW1(j1,j2)=xh11(j1,i)-xh11(j2,i)
	disyHw1HW1(j1,j2)=yh11(j1,i)-yh11(j2,i)
	diszHw1HW1(j1,j2)=zh11(j1,i)-zh11(j2,i)
	disHW1HW1(j1,j2)=disxHW1HW1(j1,j2)**2
     &	+disyHW1HW1(j1,j2)**2
     &    +diszHW1HW1(j1,j2)**2
c	Hw1-HW2
	disxHw1HW2(j1,j2)=xh11(j1,i)-xh21(j2,i)
	disyHw1HW2(j1,j2)=yh11(j1,i)-yh21(j2,i)
	diszHw1HW2(j1,j2)=zh11(j1,i)-zh21(j2,i)
	disHW1HW2(j1,j2)=disxHW1HW2(j1,j2)**2
     &	+disyHW1HW2(j1,j2)**2
     &    +diszHW1HW2(j1,j2)**2
c
c	Hw2-OW
	disxHw2OW(j1,j2)=xh21(j1,i)-xo1(j2,i)
	disyHw2OW(j1,j2)=yh21(j1,i)-yo1(j2,i)
	diszHw2OW(j1,j2)=zh21(j1,i)-zo1(j2,i)
	disHW2OW(j1,j2)=disxHW2OW(j1,j2)**2
     &	+disyHW2OW(j1,j2)**2
     &    +diszHW2OW(j1,j2)**2
c	Hw1-HW1
	disxHw2HW1(j1,j2)=xh21(j1,i)-xh11(j2,i)
	disyHw2HW1(j1,j2)=yh21(j1,i)-yh11(j2,i)
	diszHw2HW1(j1,j2)=zh21(j1,i)-zh11(j2,i)
	disHW2HW1(j1,j2)=disxHW2HW1(j1,j2)**2
     &	+disyHW2HW1(j1,j2)**2
     &    +diszHW2HW1(j1,j2)**2
c	Hw2-HW2
	disxHw2HW2(j1,j2)=xh21(j1,i)-xh21(j2,i)
	disyHw2HW2(j1,j2)=yh21(j1,i)-yh21(j2,i)
	diszHw2HW2(j1,j2)=zh21(j1,i)-zh21(j2,i)
	disHW2HW2(j1,j2)=disxHW2HW2(j1,j2)**2
     &	+disyHW2HW2(j1,j2)**2
     &    +diszHW2HW2(j1,j2)**2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	L-J 
c	L-J Ow-Ow
	r1=disOWOW(j1,j2)
	r3=sigma26/r1**3
	ESUBSUBLJ(i)=ESUBSUBLJ(i)+epsion2*(r3**2-r3)		
c	Columb
c	Columb OW-OW
	r1=disOWOW(j1,j2)
	ESUBSUBC(i)=ESUBSUBC(i)+cow*cow/(sqrt(r1))
c	Columb HW-OW
	r1=disOWHW1(j1,j2)
	ESUBSUBC(i)=ESUBSUBC(i)+cow*chw/(sqrt(r1))
	r1=disOWHW2(j1,j2)
	ESUBSUBC(i)=ESUBSUBC(i)+cow*chw/(sqrt(r1))
	r1=disHW1OW(j1,j2)
	ESUBSUBC(i)=ESUBSUBC(i)+cow*chw/(sqrt(r1))
	r1=disHW2OW(j1,j2)
	ESUBSUBC(i)=ESUBSUBC(i)+cow*chw/(sqrt(r1))
c	Columb HW-HW
	r1=disHW1HW1(j1,j2)
	ESUBSUBC(i)=ESUBSUBC(i)+chw*chw/(sqrt(r1))
	r1=disHW1HW2(j1,j2)
	ESUBSUBC(i)=ESUBSUBC(i)+chw*chw/(sqrt(r1))
	r1=disHW2HW1(j1,j2)
	ESUBSUBC(i)=ESUBSUBC(i)+chw*chw/(sqrt(r1))
	r1=disHW2HW2(j1,j2)
	ESUBSUBC(i)=ESUBSUBC(i)+chw*chw/(sqrt(r1))
	endif
220	continue
210	continue
c 
200	continue
c	calculate the average energy
	do 400 i=1, nframe
	ECNTSUBLJT=ECNTSUBLJT+ECNTSUBLJ(i)
	ECNTSUBCT=ECNTSUBCT+ECNTSUBC(i)
	ESUBSUBLJT=ESUBSUBLJT+ESUBSUBLJ(i)
	ESUBSUBCT=ESUBSUBCT+ESUBSUBC(i)
400	continue
c
	ECNTSUBLJT=ECNTSUBLJT/nw
	ECNTSUBCT=ECNTSUBCT/nw
	ESUBSUBLJT=ESUBSUBLJT/nw
	ESUBSUBCT=ESUBSUBCT/nw	
c	write out the energy
	write (14,98)
	do 300 i=1,nframe
	write (14,99) i, ECNTSUBLJ(i), factor*ECNTSUBC(i),
     &	 ESUBSUBLJ(i), factor*ESUBSUBC(i) 
300	continue
	write (14,99) 000, ECNTSUBLJT, factor*ECNTSUBCT,
     &	 ESUBSUBLJT, factor*ESUBSUBCT 
98	format ('frame, CNT-SUB-LJ CNT-SUB-C SUB-SUB-LJ SUB-SUB-C')
99	format (I5, 4f20.4)
	end