 program carbon nanotube array
	character tube_type
	double precision dis

	write(*,*)"Please choose the type of carbon nanotube"
	write(*,*)"a for type armchair"
	write(*,*)"z for type zigzag"

	read(*,'(A1)')tube_type
C# Assume one supercell has two ctcles of C atoms 
	write(*,*)"Input the value of m,n,and n_supercell"
	read(*,'(3I5)')m,n,n_supercell
	
	if(tube_type.eq.'a') then
	
	call armchair(m,n,n_supercell,dis,numberx,numbery)
	
	else
	  if(tube_type.eq.'z') then
	 call zigzag(m,n,n_supercell,dis,numberx,numbery)	
	endif
	endif

	end program carbon nanotube array

C# 子程序 armchair
	subroutine armchair(m,n,nn,dis,numberx,numbery)	
	double precision  pi,sqrt3
	parameter(pi=3.1415927,sqrt3=1.7320508)
	double precision cc,d,r,beta1,beta2
	parameter (cc=1.415)
	double precision ax(2*m),ay(2*m),az(2*m)
	double precision bx(2*m),by(2*m),bz(2*m)
	double precision axnn(nn-1,2*m),aynn(nn-1,2*m),aznn(nn-1,2*m)
	double precision bxnn(nn-1,2*m),bynn(nn-1,2*m),bznn(nn-1,2*m)
      double precision dis

	write(*,*)"Input the displacement"
	read(*,*)dis
	write(*,*)"Input the number of CNTs in x axis"
	read(*,*)numberx
	write(*,*)'Input the number of CNTs in y axis'
	read(*,*)numbery
	
	d=sqrt3/pi*cc*SQRT(REAL(m*m+m*n+n*n))
	r=d/2.0		
	beta1=ASIN(cc/2.0/r) 
	beta2=ASIN(cc/r) 

	do 12 i=1,2*m
	az(i)=0.0
12	continue  
	
	do 13 i=1,2*m
	ax(i)=r*COS(beta1+2*beta2*INT(i/2)+2*beta1*INT((i-1)/2))
	AY(I)=r*SIN(beta1+2*beta2*INT(I/2)+2*beta1*INT((i-1)/2))	
13	continue   
	
	do 14 i=1,2*m
	bz(i)=cc/2.0*sqrt3
14	continue

	
	do 15 I=1,2*m	
	BX(I)=R*COS(BETA2+2*BETA1*INT(I/2)+2*BETA2*INT((I-1)/2))	
	BY(I)=R*SIN(BETA2+2*BETA1*INT(I/2)+2*BETA2*INT((I-1)/2))	
15    continue   

	do 16 i=1,nn-1   
	do 17 j=1,2*m
	axnn(i,j)=ax(j)
	aynn(i,j)=ay(j)
	aznn(i,j)=az(j)+cc*sqrt3*i

	bxnn(i,j)=bx(j)
	bynn(i,j)=by(j)
	bznn(i,j)=bz(j)+cc*sqrt3*i
17	continue
16	continue	

	open(unit=11,file='armchair array.xyz',status='unknown')
	write(11,*)2*2*m+2*(nn-1)*2*m,d	

	do 30 k=1,numberx
	do 18 i=1,2*m
	write(11,'(A5,3F12.6)')'C',ax(i)+(k-1)*dis,ay(i),az(i)
18	continue
	do 19 i=1,2*m
	write(11,'(A5,3F12.6)')'C',bx(i)+(k-1)*dis,by(i),bz(i)
19	continue
	do 20 i=1,nn-1
	do 21 j=1,2*m
	write(11,'(A5,3F12.6)')'C',axnn(i,j)+(k-1)*dis,aynn(i,j),aznn(i,j)
21	continue
	do 22 j=1,2*m
	write(11,'(A5,3F12.6)')'C',bxnn(i,j)+(k-1)*dis,bynn(i,j),bznn(i,j)
22	continue
20	continue
30    continue

      
	do 100 k=1,numbery-1
      do 58 I=1,2*m
	write(11,'(A5,3F12.6)')'C',ax(i),ay(i)+k*dis,az(i)
58	continue
	do 59 i=1,2*m
	write(11,'(A5,3F12.6)')'C',bx(i),by(i)+k*dis,bz(i)
59	continue
	do 60 i=1,nn-1
	do 61 j=1,2*m
	write(11,'(A5,3F12.6)')'C',axnn(i,j),aynn(i,j)+k*dis,aznn(i,j)
61	continue
	do 62 J=1,2*m
	write(11,'(A5,3F12.6)')'C',bxnn(i,j),bynn(i,j)+k*dis,bznn(i,j)
62	continue
60	continue
100   continue

	close(11)
	return
	end



C# 子程序 zigzag
      subroutine zigzag(m,n,nn,dis,numberx,numbery)
	integer m,n,nn
	double precision pi,sqrt3
	parameter(pi=3.415927,sqrt3=1.7320508)
	double precision cc,d, r,z_beta,dis
	parameter(cc=1.415)
	double precision ax(m),ay(m),az(m)
	double precision bx(m),by(m),bz(m)
      double precision axnn(nn-1,m),aynn(nn-1,m),aznn(nn-1,m)
	double precision bxnn(nn-1,m),bynn(nn-1,m),bznn(nn-1,m)
      

      x=dis
      d=sqrt3/pi*cc*sqrt(real(m*m+m*n+n*n))
	r=d/2.0
	z_beta=asin(sqrt3*cc/2.0/r)


	do 122 i=1,m
	az(i)=0.0
122    continue
      
	do 133 i=1,m
	ax(i)=r*cos(z_beta+(i-1)*2*z_beta)
	ay(i)=r*sin(z_beta+(i-1)*2*z_beta)
133    continue


      do 144 i=1,m
	bz(i)=cc/2.0
144    continue
      
	do 155 i=1,m
	bx(i)=r*cos((i-1)*2*z_beta)
	by(i)=r*sin((i-1)*2*z_beta)
155    continue

	
      do 166 i=1,nn-1
	do 177 j=1,m
      axnn(i,j)=ax(j)
	aynn(i,j)=ay(j)
	aznn(i,j)=az(j)+cc*2*int((i+1)/2)+cc*int(i/2)

	bxnn(i,j)=bx(j)
	bynn(i,j)=by(j)
	bznn(i,j)=bz(j)+cc*int((i+1)/2)+cc*2*int(i/2)
177   continue
166   continue

      open(unit=111,file='zigzag.xyz',status='unknown')
	write(111,'(I6)')2*m+2*(nn-1)*m

	do 188 i=1,m	
	write(111,'(A5,3f12.6)')'C',ax(i),ay(i),az(i)
188   continue
      
	do 199 i=1,m
	write(111,'(A5,3f12.6)')'C',bx(i),by(i),bz(i)
199    continue

      do 200 i=1,nn-1
	do 211 j=1,m
      write(111,'(A5,3f12.6)')'C',axnn(i,j),aynn(i,j),aznn(i,j)
211   continue
      
	do 222 j=1,m
	write(111,'(A5,3f12.6)')'C',bxnn(i,j),bynn(i,j),bznn(i,j)
222   continue
200   continue

  
      do 288 i=1,m	
	write(111,'(A5,3f12.6)')'C',ax(i)+dis,ay(i),az(i)
288   continue
      
	do 299 i=1,m
	write(111,'(A5,3f12.6)')'C',bx(i)+dis,by(i),bz(i)
299    continue

      do 300 i=1,nn-1
	do 311 j=1,m
      write(111,'(A5,3f12.6)')'C',axnn(i,j)+dis,aynn(i,j),aznn(i,j)
311   continue
      
	do 322 j=1,m
	write(111,'(A5,3f12.6)')'C',bxnn(i,j)+dis,bynn(i,j),bznn(i,j)
322   continue
300   continue
      
      close(111)

	return
	end

 	
	
	