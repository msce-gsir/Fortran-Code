        program XYdensity
		 integer totmg,totli,totcl,totframe,frame
		 integer totow
		 parameter(totmg=13,totli=13,totcl=39)
		 parameter(totow=1088,totframe=60000)
		 integer binx,biny,maxbin
		 parameter(maxbin=1000)
		 double precision delt,pi
		 parameter(delt=0.015,pi=3.1415927)
		 double precision xp,yp,vol,temp
		 double precision deltx,delty,temp1,temp2
C ///////////////////////////////////////////////////////////////////////
		 double precision cnmg(totframe),cnli(totframe)
		 double precision cncl(totframe),cnow(totframe)
		 double precision cnhw(totframe)
C ///////////////////////////////////////////////////////////////////////
		 double precision histli(maxbin,maxbin),zli(maxbin,maxbin)
		 double precision histcl(maxbin,maxbin),zcl(maxbin,maxbin)
		 double precision histmg(maxbin,maxbin),zmg(maxbin,maxbin)
		 double precision histow(maxbin,maxbin),zow(maxbin,maxbin)
		 double precision histhw(maxbin,maxbin),zhw(maxbin,maxbin)
C ///////////////////////////////////////////////////////////////////////
		 double precision xlowerlimit,xupperlimit
		 double precision ylowerlimit,yupperlimit
		 double precision zlowerlimit,zupperlimit
		 double precision radius,cenx,ceny,cenz
C ///////////////////////////////////////////////////////////////////////
		 double precision lix(totli),liy(totli),liz(totli)
		 double precision mgx(totmg),mgy(totmg),mgz(totmg)
		 double precision clx(totcl),cly(totcl),clz(totcl)
		 double precision owx(totow),owy(totow),owz(totow)
		 double precision hw1x(totow),hw1y(totow),hw1z(totow)
		 double precision hw2x(totow),hw2y(totow),hw2z(totow)
C ///////////////////////////////////////////////////////////////////////
         character*20 lin(totli),cln(totcl),mgn(totmg)
		 character*20 own(totow),hw1n(totow),hw2n(totow)
C ///////////////////////////////////////////////////////////////////////
         open(unit=10,file='radiusandcenter.txt',status='unknown')	
         open(unit=11,file='mg_dist.gro',status='unknown')		
         open(unit=12,file='li_dist.gro',status='unknown')
         open(unit=13,file='cl_dist.gro',status='unknown')
		 open(unit=14,file='sol_dist.gro',status='unknown')
C		 open(unit=20,file='hw_dist.gro',status='unknown')
		 open(unit=15,file='2D-MG.txt',status='unknown')
		 open(unit=16,file='2D-LI.txt',status='unknown')
		 open(unit=17,file='2D-CL.txt',status='unknown')
		 open(unit=18,file='2D-OW.txt',status='unknown')
		 open(unit=19,file='2D-HW.txt',status='unknown')
C ///////////////////////////////////////////////////////////////////////
         do 100 binx=1,maxbin
          do 101 biny=1,maxbin
             histmg(binx,biny)=0.0	
             histli(binx,biny)=0.0	
             histcl(binx,biny)=0.0	
             histow(binx,biny)=0.0	
             histhw(binx,biny)=0.0			 
101       continue
100     continue	

        do 102 frame=1,totframe
		   cnmg(frame)=0.0
		   cnli(frame)=0.0
		   cncl(frame)=0.0
		   cnow(frame)=0.0
		   cnhw(frame)=0.0
102     continue
C ///////////////////////////////////////////////////////////////////////////////
		  read(10,1005)radius,cenx,ceny,cenz
		  xlowerlimit=cenx-radius
	      xupperlimit=cenx+radius
		  ylowerlimit=ceny-radius
		  yupperlimit=ceny+radius
          zlowerlimit=cenz-0.20
		  zupperlimit=cenz+0.20
C ///////////////////////////////////////////////////////////////////////////////	 
        do 200 frame=1,totframe

         read(11,*)
		 read(11,*)
		 do 300 i=1,totmg
		  read(11,1001)mgn(i),mgx(i),mgy(i),mgz(i)
300      continue
         read(11,*)
		 
		 read(12,*)
		 read(12,*)
		 do 301 i=1,totli
		  read(12,1001)lin(i),lix(i),liy(i),liz(i)
301      continue
         read(12,*)
		 
		 read(13,*)
		 read(13,*)
		 do 302 i=1,totcl
		  read(13,1001)cln(i),clx(i),cly(i),clz(i)
302      continue
         read(13,*)
		 
		 read(14,*)
		 read(14,*)
		 do 303 i=1,totow
		  read(14,1001)own(i),owx(i),owy(i),owz(i)
		  read(14,1001)hw1n(i),hw1x(i),hw1y(i),hw1z(i)
		  read(14,1001)hw2n(i),hw2x(i),hw2y(i),hw2z(i)
303      continue
         read(14,*)
		 
C		 read(20,*)
C		 read(20,*)
C		 do 304 i=1,tothw
C		  read(20,1001)hwn(i),hwx(i),hwy(i),hwz(i)
C04      continue
C        read(20,*)
C ///////////////////////////////////////////////////////////////////////////////
        do 400 j=1,totmg
		if(mgz(j).le.zupperlimit.and.mgz(j).ge.zlowerlimit)then
		if(mgx(j).le.xupperlimit.and.mgx(j).ge.xlowerlimit)then
		if(mgy(j).le.yupperlimit.and.mgy(j).ge.ylowerlimit)then
		    deltx=mgx(j)-cenx
			delty=mgy(j)-ceny
			temp1=sqrt(deltx**2+delty**2)
		if(temp1.le.radius)then
		    cnmg(frame)=cnmg(frame)+1.0
		    binx=int(mgx(j)/delt)+1
			biny=int(mgy(j)/delt)+1
			histmg(binx,biny)=histmg(binx,biny)+1.0
		end if
		end if
		end if
		end if
400     continue
C ///////////////////////////////////////////////////////////////////////////////
        do 401 j=1,totli
		if(liz(j).le.zupperlimit.and.liz(j).ge.zlowerlimit)then
		if(lix(j).le.xupperlimit.and.lix(j).ge.xlowerlimit)then
		if(liy(j).le.yupperlimit.and.liy(j).ge.ylowerlimit)then
		    deltx=lix(j)-cenx
			delty=liy(j)-ceny
			temp1=sqrt(deltx**2+delty**2)
		if(temp1.le.radius)then
		    cnli(frame)=cnli(frame)+1.0
		    binx=int(lix(j)/delt)+1
			biny=int(liy(j)/delt)+1
			histli(binx,biny)=histli(binx,biny)+1.0
		end if
		end if
		end if
		end if
401     continue
C ///////////////////////////////////////////////////////////////////////////////
        do 402 j=1,totcl
		if(clz(j).le.zupperlimit.and.clz(j).ge.zlowerlimit)then
		if(clx(j).le.xupperlimit.and.clx(j).ge.xlowerlimit)then
		if(cly(j).le.yupperlimit.and.cly(j).ge.ylowerlimit)then
		    deltx=clx(j)-cenx
			delty=cly(j)-ceny
			temp1=sqrt(deltx**2+delty**2)
		if(temp1.le.radius)then
		    cncl(frame)=cncl(frame)+1.0
		    binx=int(clx(j)/delt)+1
			biny=int(cly(j)/delt)+1
			histcl(binx,biny)=histcl(binx,biny)+1.0
		end if
		end if
		end if
		end if
402     continue
C ///////////////////////////////////////////////////////////////////////////////
        do 404 j=1,totow
		if(owz(j).le.zupperlimit.and.owz(j).ge.zlowerlimit)then
		if(owx(j).le.xupperlimit.and.owx(j).ge.xlowerlimit)then
		if(owy(j).le.yupperlimit.and.owy(j).ge.ylowerlimit)then
		    deltx=owx(j)-cenx
			delty=owy(j)-ceny
			temp1=sqrt(deltx**2+delty**2)
		if(temp1.le.radius)then
		    cnow(frame)=cnow(frame)+1.0
		    binx=int(owx(j)/delt)+1
			biny=int(owy(j)/delt)+1
			histow(binx,biny)=histow(binx,biny)+1.0
			
			binx=int(hw1x(j)/delt)+1
			biny=int(hw1y(j)/delt)+1
			histhw(binx,biny)=histhw(binx,biny)+1.0
			
			binx=int(hw2x(j)/delt)+1
			biny=int(hw2y(j)/delt)+1
			histhw(binx,biny)=histhw(binx,biny)+1.0
			
		end if	
		end if
		end if
		end if
404     continue 

C        do 405 j=1,tothw
C		if(hwz(j).le.zupperlimit.and.hwz(j).ge.zlowerlimit)then
C		if(hwx(j).le.xupperlimit.and.hwx(j).ge.xlowerlimit)then
C		if(hwy(j).le.yupperlimit.and.hwy(j).ge.ylowerlimit)then
C		    deltx=hwx(j)-cenx(frame)
C			delty=hwy(j)-ceny(frame)
C			temp1=sqrt(deltx**2+delty**2)
C		if(temp1.le.radius(frame))then
C		    cnhw(frame)=cnhw(frame)+1.0
C		    binx=int(hwx(j)/delt)+1
C			biny=int(hwy(j)/delt)+1
C			histhw(binx,biny)=histhw(binx,biny)+1.0
C		end if	
C		end if
C		end if
C		end if
C405     continue 
C ///////////////////////////////////////////////////////////////////////////////   
200     continue
C /////////////////////////////////////////////////////////////////////////////// 
        vol=pi*(radius**2)*delt
C /////////////////////////////////////////////////////////////////////////////// 
        xlowerlimit=cenx-radius
		xupperlimit=cenx+radius
		ylowerlimit=ceny-radius
		yupperlimit=ceny+radius
		
        do 501 binx=1,maxbin
		 do 502 biny=1,maxbin
		   xp=binx*delt
		   yp=biny*delt
		if(xp.le.xupperlimit.and.xp.ge.xlowerlimit)then
		if(yp.le.yupperlimit.and.yp.ge.ylowerlimit)then
		   xp=xp-cenx
		   yp=yp-ceny  

           histmg(binx,biny)=histmg(binx,biny)/totframe/vol
           write(15,1008)xp,yp,histmg(binx,biny)

           histli(binx,biny)=histli(binx,biny)/totframe/vol
           write(16,1008)xp,yp,histli(binx,biny)

           histcl(binx,biny)=histcl(binx,biny)/totframe/vol
           write(17,1008)xp,yp,histcl(binx,biny)

           histow(binx,biny)=histow(binx,biny)/totframe/vol
           write(18,1008)xp,yp,histow(binx,biny)

           histhw(binx,biny)=histhw(binx,biny)/totframe/vol
           write(19,1008)xp,yp,histhw(binx,biny)

        end if
        end if
502     continue
501     continue
C ///////////////////////////////////////////////////////////////////////////////
1001    format(A20,3f8.3)
1005    format(4f8.3)
1008    format(3f12.3)
 
        close(10) 
		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)
		close(17)
		close(18)
		close(19)
		close(20)
		
		end
		
