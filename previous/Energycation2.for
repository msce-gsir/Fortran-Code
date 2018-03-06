	program energy
!	this code calculate the non-bond interactions between molecules
!	they were divided into the CNT-ion, ion-water(in shell), and CNT-water(shell).
	integer nsol,nframe,nc,nion,nio
	parameter (nsol=328,nframe=1001, nio=1, nc=2320)
	double precision xo,yo,zo
	double precision xh1, yh1, zh1
	double precision xh2, yh2, zh2
	double precision factor
	integer index, index1, frame, n1,now,nhw1,nhw2
	double precision xc(nc), yc(nc), zc(nc)
	character*8 head,head1,head2
	double precision xion,yion, zion



	double precision sigmaOw,sigmaHw
	double precision epsionC, sigmaC 
	double precision epison1,epison2
	double precision episonOw,episonHw
	double precision sigma1,sigma2
	double precision cIon,cow,chw,rmax,l12,l22

	double precision ECNTWater, EIonWater

	open (unit=10, file='CNT1010.gro', status='unknown')
	open (unit=11,file='cation.gro',status='unknown')
	open (unit=12,file='sol.gro',status='unknown')
	open (unit=13, file='EnergyinCNT.cation.dat', status='unknown')
	open (unit=14, file='r_cation.txt', status='unknown')
!	read the coordination for nanotube
	read(10,*)
	read (10,*) 
	do i=1,nc
	read (10,*) head, head, n1, xc(i),yc(i),zc(i)
	end do
	read (14,*) rmax
	
!	several forcefield parameters
	factor=138.935
	
!	CNT cc represent the charge of carbon,cow,chw,etc have the same mean

	epsionC=2.29288e-01
	sigmaC=3.55000e-01
	cc=0.0
!	OW
	epsionOW=6.50194e-01
	sigmaOW=3.16557e-01
	cow=-0.8476
!	HW
	epsionHW=0.00000e+00
	sigmaHW=0.00000e+00
	chw=0.4238
!	Ion
!	Li
!	epsionIon=0.0765
!	sigmaIon=0.213
!	cIon=1.00
!	K
	epsionIon=1.37235e-03
	sigmaIon=4.93463e-01
	cIon=1.00
!	Na
!	epsionIon=1.15980e-02
!	sigmaIon=3.33045e-01
!	cIon=1.00
!	F
!	epsionIon=3.01248e+00
!	sigmaIon=2.73295e-01
!	cIon=-1.00
!	Cl
!	epsionIon=4.92833e-01
!	sigmaIon=4.41724e-01
!	cIon=-1.00
!	K0
!	epsionIon=1.37235e-03
!	sigmaIon=4.93463e-01
!	cIon=0.00
!	F0
!	epsionIon=3.01248e+00
!	sigmaIon=2.73295e-01
!	cIon=0.00
!	Mg
!	epsionIon=3.66118e+00
!	sigmaIon=1.64447e-01
!	cIon=2.00
!	Ca
!	epsionIon=1.88136e+00
!	sigmaIon=2.41203e-01
!	cIon=2.00
!	K-
!	epsionIon=1.37235e-03
!	sigmaIon=4.93463e-01
!	cIon=-1.00
!	Cl+
!	epsionIon=4.92833e-01
!	sigmaIon=4.41724e-01
!	cIon=1.00
!c    ///////////////////
!	CNT-OW
	epsion1=4*sqrt(epsionOW*epsionC)
	sigma1=(sigmaOW+sigmaC)/2
	sigma1=sigma1**6
!	ion-OW
	epsion2=4*sqrt(epsionOW*epsionIon)
	sigma2=(sigmaOW+sigmaIon)/2
	sigma2=sigma2**6
!	ions
	
!	end of force field parameter
!	read the coordination in frame

!     main loop
!	read the coordination of CNT
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
	do 5000 k=1,nframe
!	read the ion coord.
	read (11,*)
	read (11,*)
!	do 100 i=1,1
	read (11,*) head1,head2, nion, xion,yion,zion
!	write (*,*) xion,yion,zion
!100	continue	
	read (11,*)
!	read the water coord.
	read (12,*)
	read (12,*)
	do 101 i=1,nsol
	read (12,*) head1,head2, now, xo,yo,zo
	read (12,*) head1,head2, nhw1, xh1,yh1,zh1
	read (12,*) head1,head2, nhw2, xh2,yh2,zh2
!	distance
	l22=(xion-xo)**2+(yion-yo)**2+(zion-zo)**2 ! ion-ow
!	write (*,*) l22
	if (l22.lt.rmax**2.and.l22.gt.0.01) then
	k1=k1+1
!	calculate the energy
!	CNT-hydrated water
	do 200 j=1,nc
!	vdw
	l12=(xc(j)-xo)**2+(yc(j)-yo)**2+(zc(j)-zo)**2
	if (l12.lt.2.0.and.l12.gt.0.01) then
	r3=sigma1/(l12**3)
	ECNTWater=ECNTWater+epsion1*(r3**2-r3)
	endif	
200	continue
	
!	Ion-hydrated water
!	vdw
	r3=sigma2/(l22**3)
	EIonWater=EIonWater+epsion2*(r3**2-r3)
!	charge
!	Ion-Ow
	r1=l22	
	EIonWater=EIonWater+factor*cow*cIon/sqrt(r1)
!	Ion-Hw
	r1=(xion-xh1)**2+(yion-yh1)**2+(zion-zh1)**2 ! ion-hw1
	EIonWater=EIonWater+factor*chw*cIon/sqrt(r1)
	r1=(xion-xh2)**2+(yion-yh2)**2+(zion-zh2)**2 ! ion-hw2
	EIonWater=EIonWater+factor*chw*cIon/sqrt(r1)
	endif
101	continue 
	read (12,*)	
5000	continue
!	write (*,*) k1, ECNTWater, EIonWater
!
	ECNTWater=ECNTWater/k1
	EIonWater=EIonWater/k1	
!	write out the energy

	write (13,98)
	write (13,99) ECNTWater, EIonWater 
98	format ('CNT-Water Ion-Water')
99	format (2f20.4)
	end