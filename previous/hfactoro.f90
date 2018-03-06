        program hfactor
		 use xdr, only: trrfile
         implicit none
         type(trrfile) :: trr
         real(kind=8), parameter :: pi = 3.141592653589793238462643383
         real(kind=8), parameter :: delt = 0.04 !nm
         real(kind=8), allocatable :: pos(:,:)
		 double precision hydmg,hydmgs,hydli,hydlis
		 double precision htmgf(0:180),htmgs(0:180),htlif(0:180),htlis(0:180)
		 real(kind=8) :: radius,cenx,ceny,cenz
		 real(kind=8) :: xlowerlimit,xupperlimit
		 real(kind=8) :: ylowerlimit,yupperlimit
		 real(kind=8) :: zlowerlimit,zupperlimit
		 real(kind=8) :: deltx,delty,temp
		 integer :: maxbin,bin,ng,i,j,k,npos,cn,cm,na,co1,co2
		 parameter(maxbin=1000,na=3)
		 integer img,fmg,ili,fli,iw,fw,io,fo
		 double precision hw3x,hw3y,hw3z
		 double precision cosa,cosl,cosu,thita ,t                         ! cosl=-0.72,cosu=-1.0
		 parameter(cosl=136.0,cosu=180.0)
		 double precision mgowxsq,mgowysq,mgowzsq,mgow
		 double precision mgoxsq,mgoysq,mgozsq,mgo
		 double precision hw3owxsq,hw3owysq,hw3owzsq,hw3ow
		 double precision hw3mgxsq,hw3mgysq,hw3mgzsq,hw3mg
		 double precision liowxsq,liowysq,liowzsq,liow
		 double precision lioxsq,lioysq,liozsq,lio
		 double precision hw3lixsq,hw3liysq,hw3lizsq,hw3li
		 double precision summgf,summgs,sumlif,sumlis
		 double precision ptmgf,ptmgs,ptlif,ptlis
		 double precision hfmgf,hfmgs,hflif,hflis
		 
		 open(unit=10,file='hyd.txt',status='unknown')
		 open(unit=11,file='serial.txt',status='unknown')
		 open(unit=12,file='radiusandcenter.txt',status='unknown')
		 open(unit=13,file='HF-MG-F.txt',status='unknown')
		 open(unit=14,file='HF-LI-F.txt',status='unknown')
		 open(unit=15,file='HF-MG-S.txt',status='unknown')
		 open(unit=16,file='HF-LI-S.txt',status='unknown')
		 open(unit=17,file='HF-MG-LI-F-S.txt',status='unknown')
		 open(unit=18,file='hfactor.txt',status='unknown')
		 
		 read(10,*)hydmg,hydmgs
		 read(10,*)hydli,hydlis
		 
		 read(11,*)img,fmg
		 read(11,*)ili,fli
		 read(11,*)iw,fw
		 read(11,*)io,fo
		 
		 read(12,1002)radius,cenx,ceny,cenz
		 	 
		 xlowerlimit=cenx-radius
		 xupperlimit=cenx+radius
		 ylowerlimit=ceny-radius
		 yupperlimit=ceny+radius
		 zlowerlimit=cenz-0.25
		 zupperlimit=cenz+0.25
		 
		 htmgf=0.0
		 htmgs=0.0
		 htlif=0.0
		 htlis=0.0
		 
		 ng=0
		 
		 call trr % init("../nvt_0.2.trr")
		 npos = trr % NATOMS
		 allocate(pos(3, npos))
		 
		 call trr % read	

         do while ( trr % STAT == 0 )
            ng = ng + 1                          ! count total frame
            pos = trr % pos(:, 1:trr % NATOMS)   !get the position of atoms sort in npt_em.gro
			if(ng.gt.5000)then                   !avoid first 5ns 

! MG	
            do i=img,fmg
        if(pos(1,i).le.xupperlimit.and.pos(1,i).ge.xlowerlimit)then
		if(pos(2,i).le.yupperlimit.and.pos(2,i).ge.ylowerlimit)then
		if(pos(3,i).le.zupperlimit.and.pos(3,i).ge.zlowerlimit)then
		        deltx=((pos(1,i)-cenx))**2
				delty=((pos(2,i)-ceny))**2
				temp=sqrt(deltx+delty)

                if(temp.le.radius)then
                   
				    cn=0
					cm=0
					co1=0
					co2=0
					
					do j=io,fo
					
					    mgoxsq=(pos(1,i)-pos(1,j))**2
						mgoysq=(pos(2,i)-pos(2,j))**2
						mgozsq=(pos(3,i)-pos(3,j))**2
						mgo=sqrt(mgoxsq+mgoysq+mgozsq)
					    if(mgo.le.hydmg)then
						    co1=co1+1
						end if
						if(mgo.le.hydmgs.and.mgo.gt.hydmg)then
						    co2=co2+1
						end if
					end do
					
						
                    do j=iw,fw,3

                    	mgowxsq=(pos(1,i)-pos(1,j))**2	
                        mgowysq=(pos(2,i)-pos(2,j))**2	
                        mgowzsq=(pos(3,i)-pos(3,j))**2	
                        mgow=sqrt(mgowxsq+mgowysq+mgowzsq)	
                     
                    if(mgow.le.hydmg)then
					if(co1.eq.0)then

					    hw3x=(pos(1,j+1)+pos(1,j+2))/2.0
                        hw3y=(pos(2,j+1)+pos(2,j+2))/2.0
                        hw3z=(pos(3,j+1)+pos(3,j+2))/2.0	

                        hw3owxsq=(hw3x-pos(1,j))**2	
                        hw3owysq=(hw3y-pos(2,j))**2
                        hw3owzsq=(hw3z-pos(3,j))**2		
                        hw3ow=sqrt(hw3owxsq+hw3owysq+hw3owzsq)	
						
						hw3mgxsq=(hw3x-pos(1,i))**2
                        hw3mgysq=(hw3y-pos(2,i))**2
						hw3mgzsq=(hw3z-pos(3,i))**2
						hw3mg=sqrt(hw3mgxsq+hw3mgysq+hw3mgzsq)
						
				        cosa=(mgow**2+hw3ow**2-hw3mg**2)/(2.0*mgow*hw3ow)
				        thita=acos((mgow**2+hw3ow**2-hw3mg**2)/(2.0*mgow*hw3ow))/pi*180.0
				        bin=int(thita)+1.0
				        htmgf(bin)=htmgf(bin)+1.0
				    end if
					end if
					
					if(mgow.le.hydmgs.and.mgow.gt.hydmg)then
					if(co2.ne.0)then
					
					    hw3x=(pos(1,j+1)+pos(1,j+2))/2.0
                        hw3y=(pos(2,j+1)+pos(2,j+2))/2.0
                        hw3z=(pos(3,j+1)+pos(3,j+2))/2.0	

                        hw3owxsq=(hw3x-pos(1,j))**2	
                        hw3owysq=(hw3y-pos(2,j))**2
                        hw3owzsq=(hw3z-pos(3,j))**2		
                        hw3ow=sqrt(hw3owxsq+hw3owysq+hw3owzsq)	
						
						hw3mgxsq=(hw3x-pos(1,i))**2
                        hw3mgysq=(hw3y-pos(2,i))**2
						hw3mgzsq=(hw3z-pos(3,i))**2
						hw3mg=sqrt(hw3mgxsq+hw3mgysq+hw3mgzsq)
						
				        cosa=(mgow**2+hw3ow**2-hw3mg**2)/(2.0*mgow*hw3ow)
				        thita=acos((mgow**2+hw3ow**2-hw3mg**2)/(2.0*mgow*hw3ow))/pi*180.0
				        bin=int(thita)+1.0
				        htmgs(bin)=htmgs(bin)+1.0    
					end if
					end if
					end do
				end if
			end if
			end if
			end if
			    end do
				
				
! Li

            do i=ili,fli
		if(pos(1,i).le.xupperlimit.and.pos(1,i).ge.xlowerlimit)then
		if(pos(2,i).le.yupperlimit.and.pos(2,i).ge.ylowerlimit)then
		if(pos(3,i).le.zupperlimit.and.pos(3,i).ge.zlowerlimit)then
		        deltx=((pos(1,i)-cenx))**2
				delty=((pos(2,i)-ceny))**2
				temp=sqrt(deltx+delty)

                if(temp.le.radius)then
                   
				    cn=0
					cm=0
					co1=0
					co2=0
					
					do j=io,fo
					
					    lioxsq=(pos(1,i)-pos(1,j))**2
						lioysq=(pos(2,i)-pos(2,j))**2
						liozsq=(pos(3,i)-pos(3,j))**2
						lio=sqrt(lioxsq+lioysq+liozsq)
					    if(lio.le.hydli)then
						    co1=co1+1
						end if
						if(lio.le.hydlis.and.lio.gt.hydli)then
						    co2=co2+1
						end if
					end do
						
                    do j=iw,fw,3
					
				        liowxsq=(pos(1,i)-pos(1,j))**2	
                        liowysq=(pos(2,i)-pos(2,j))**2	
                        liowzsq=(pos(3,i)-pos(3,j))**2	
                        liow=sqrt(liowxsq+liowysq+liowzsq)	
				
					if(liow.le.hydli)then
					if(co1.ne.0)then

					    hw3x=(pos(1,j+1)+pos(1,j+2))/2.0
                        hw3y=(pos(2,j+1)+pos(2,j+2))/2.0
                        hw3z=(pos(3,j+1)+pos(3,j+2))/2.0	

                        hw3owxsq=(hw3x-pos(1,j))**2	
                        hw3owysq=(hw3y-pos(2,j))**2
                        hw3owzsq=(hw3z-pos(3,j))**2		
                        hw3ow=sqrt(hw3owxsq+hw3owysq+hw3owzsq)	
						
						hw3lixsq=(hw3x-pos(1,i))**2
                        hw3liysq=(hw3y-pos(2,i))**2
						hw3lizsq=(hw3z-pos(3,i))**2
						hw3li=sqrt(hw3lixsq+hw3liysq+hw3lizsq)
						
				        cosa=(liow**2+hw3ow**2-hw3li**2)/(2.0*liow*hw3ow)
				        thita=acos((liow**2+hw3ow**2-hw3li**2)/(2.0*liow*hw3ow))/pi*180.0
				        bin=int(thita)+1.0
				        htlif(bin)=htlif(bin)+1.0
				    end if	
					end if
						
                    if(liow.le.hydlis.and.liow.gt.hydli)then
					if(co2.ne.0)then
					
					    hw3x=(pos(1,j+1)+pos(1,j+2))/2.0
                        hw3y=(pos(2,j+1)+pos(2,j+2))/2.0
                        hw3z=(pos(3,j+1)+pos(3,j+2))/2.0	

                        hw3owxsq=(hw3x-pos(1,j))**2	
                        hw3owysq=(hw3y-pos(2,j))**2
                        hw3owzsq=(hw3z-pos(3,j))**2		
                        hw3ow=sqrt(hw3owxsq+hw3owysq+hw3owzsq)	
						
						hw3lixsq=(hw3x-pos(1,i))**2
                        hw3liysq=(hw3y-pos(2,i))**2
						hw3lizsq=(hw3z-pos(3,i))**2
						hw3li=sqrt(hw3lixsq+hw3liysq+hw3lizsq)
						
				        cosa=(liow**2+hw3ow**2-hw3li**2)/(2.0*liow*hw3ow)
				        thita=acos((liow**2+hw3ow**2-hw3li**2)/(2.0*liow*hw3ow))/pi*180.0
				        bin=int(thita)+1.0
				        htlis(bin)=htlis(bin)+1.0    
					end if    	
                    end if					
		            end do
				end if
			end if
			end if
			end if
			    end do
				
		    end if
		
		call trr % read
			
        end do		

        ! 5. Close the file
        call trr % close
		
		write(*,1008)ng
		
		summgf=0.0
		summgs=0.0
		sumlif=0.0
		sumlis=0.0
		
		do bin=0,180
		
		    summgf=summgf+htmgf(bin)
			summgs=summgs+htmgs(bin)
			sumlif=sumlif+htlif(bin)
			sumlis=sumlis+htlis(bin)
		end do
		
		hfmgf=0.0
		hfmgs=0.0
		hflif=0.0
		hflis=0.0
		
		do bin=0,180
		
		    ptmgf=(htmgf(bin)/summgf)*100.0
		    ptmgs=(htmgs(bin)/summgs)*100.0
			ptlif=(htlif(bin)/sumlif)*100.0
			ptlis=(htlis(bin)/sumlis)*100.0
			
!			t=cos(bin*pi/180.0)    ! make more attention about this, reference:http://blog.sina.com.cn/s/blog_4b7a74bb0100flm2.html
			
			write(13,1003)bin,ptmgf,htmgf(bin)
			write(14,1003)bin,ptlif,htlif(bin)
			write(15,1003)bin,ptmgs,htmgs(bin)
			write(16,1003)bin,ptlis,htlis(bin)
			write(17,1004)bin,ptmgf,ptlif,ptmgs,ptlis
			
			if(bin.le.cosu.and.bin.ge.cosl)then
			    hfmgf=hfmgf+htmgf(bin)
				hfmgs=hfmgs+htmgs(bin)
				hflif=hflif+htlif(bin)
				hflis=hflis+htlis(bin)
			end if
			
		end do
		
		write(18,1005)'mg first factor = ', hfmgf/summgf
		write(18,1005)'mg second factor = ', hfmgs/summgs
		write(18,1005)'li first factor = ', hflif/sumlif
		write(18,1005)'li second factor = ', hflis/sumlis
		
1002    format(4f8.3)
1003    format(I5,2f8.3)
!1003    format(3f8.3)
1004    format(I5,4f8.3)
1005    format(A20,f8.3)
1008    format(I8)

        close(10)
		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)
		close(17)
		close(18)
		
		end
        			