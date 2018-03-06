	program hydrationeffect
c	the codes to calculate the hydration paramter
	parameter (nframe=10000, nion=1, nsol=328)
	parameter (maxbin=100,delt=0.02)
	double precision xow(nsol),yow(nsol),zow(nsol)
	double precision xhw1(nsol),yhw1(nsol),zhw1(nsol)
	double precision xhw2(nsol),yhw2(nsol),zhw2(nsol)
	double precision xp(nsol),yp(nsol),zp(nsol)
	double precision zdo2(nsol)
	double precision xion(nion),yion(nion),zion(nion)
	double precision l12,l32,l22,theta1
	integer nio, now, nhw1,nhw2
	double precision nslice(maxbin)
	integer bin
	character*6 head1,head2
c	open file
	open (unit=10,file='cation.gro',status='unknown')
	open (unit=11,file='sol.gro',status='unknown')
	open (unit=12,file='orient_cation.txt',status='unknown')
	open (unit=13,file='r_cation.txt',status='unknown')
c
	read (13,*) rmax

c	main loop
	do 5000 k=1,nframe
c	read the ion coord.
	read (10,*)
	read (10,*)
	do 100 i=1,nion
	read (10,*) head1,head2, nio, xion(i),yion(i),zion(i)
100	continue	
	read (10,*)
c	read the water coord.
	read (11,*)
	read (11,*)
	do 101 i=1,nsol
	read (11,*) head1,head2, now, xow(i),yow(i),zow(i)
	read (11,*) head1,head2, nhw1, xhw1(i),yhw1(i),zhw1(i)
	read (11,*) head1,head2, nhw2, xhw2(i),yhw2(i),zhw2(i)
	xp(i)=(xhw1(i)+xhw2(i))*0.5
	yp(i)=(yhw1(i)+yhw2(i))*0.5
	zp(i)=(zhw1(i)+zhw2(i))*0.5
	zdo2(i)=(xp(i)-xow(i))**2+(yp(i)-yow(i))**2+(zp(i)-zow(i))**2
101	continue	
	read (11,*)
c
c	calculate the paramter
	do 200 i=1,nion
	do 201 j=1,nsol
	l22=(xion(i)-xow(j))**2+(yion(i)-yow(j))**2+(zion(i)-zow(j))**2 ! ion-ow
		if (l22.lt.rmax**2) then
	l12=(xp(j)-xion(i))**2+(yp(j)-yion(i))**2+(zp(j)-zion(i))**2 ! ion-dipole
		l32=zdo2(j) ! ow-dipole
		theta1=(l22+l32-l12)/sqrt(l22*l32)/2.0
		bin=int((theta1+1.0)/delt)+1
		nslice(bin)=nslice(bin)+1
		endif
201	continue
200	continue
5000	continue
c	normalize
	do 400 i=1,maxbin
	sum=sum+nslice(i)
400	continue	
	do 500 i=1,maxbin
	nslice(i)=nslice(i)/sum
500	continue
c	write out
	do 600 i=1,maxbin
	write (12,99) (i-1)*delt-1.0, nslice(i)
600	continue
99	format (2f12.6)
	end