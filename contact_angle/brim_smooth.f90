    program outline_output
    implicit none
    
    integer id, mol_id, n_atom, n_mol, n_conf, i, j, k, id_new
    parameter (n_atom= 3, n_mol= 15625, n_conf= 601)
    real*8 x, y, z, occupancy, bfactor, lx, ly, lz, bin, den, den_total,&
&   den_ave, vol
    parameter (lx= 413.6, ly= 413.6, lz= 700, bin= 2.82)
    integer n_x, n_y, n_z, ix, iy, iz, num
    real*8,allocatable :: n_bin(:, :, :)
    character(len=12) label, atom_name, resname, label_new, a_type
    
    open(1, file='smoothp3.pdb')
    open(2, file='brim_smoothp3.txt')
    open(3, file='brim_smoothp3.gro')
    open(4, file='brim_smoothp3_2.txt')
    
    n_x= int(lx/bin) + 1 
    n_y= int(ly/bin) + 1 
    n_z= int(lz/bin) + 1
    allocate(n_bin(n_x, n_y, n_z))
    
    write(*,*) n_x, n_y, n_z

    do i= 1, n_x
        
       do j= 1,n_y
            
            do k= 1, n_z
            
            n_bin(i,j,k)= 0
            
            enddo
        
        enddo
        
    enddo
    
    do i= 1, n_conf
        
        do j= 1, 5
            
            read(1,*)
            
        enddo
        
!!!!!!!!! only use the O atom coordinate as dealing !!!!!!!!!!!       
        
        do j= 1, n_mol
            
            read(1,*) label, id, atom_name, resname, mol_id, x, y, z, occupancy, bfactor

            do k= 1, 2

               read(1,*)

            enddo

            ix= int(x/bin) + 1
            iy= int(y/bin) + 1
            iz= int(z/bin) + 1
            n_bin(ix,iy,iz)= n_bin(ix,iy,iz) + 1
            
        enddo

        do j= 1, 2
          
           read(1,*)

        enddo
        
    enddo

    vol= (bin/10.0)**3
    id_new= 0
    label_new= "BAS"
    a_type= "A"
    
    do i= 1, n_x
   
       do j= 1, n_y

          do k= 1, n_z

             n_bin(i, j, k)= n_bin(i, j, k)/n_conf
 
             den= (n_bin(i,j,k)/(6.02E23))*72.0/(vol/(1.0E21))

             if (den > 0.55 .and. den < 0.65)        then

                    write(2,100) i*bin/10.0, j*bin/10.0, k*bin/10.0
     
                    id_new= id_new + 1
  
                    write(3,101) id_new, label_new, a_type, id_new, i*bin/10.0,&
&                                j*bin/10.0, k*bin/10.0

             endif

             if (den > 0.55 .and. den < 0.65 .and. k*bin > 31.0)    then
 
                   write(4,100) i*bin/10.0, j*bin/10.0, k*bin/10.0

             endif


          enddo
  
       enddo

    enddo
            
    close(1)  
    
100 format(F6.2,1X,F6.2,1X,F6.2)
101 format(I4,A3,6X,A1,1X,I4,2X,F6.3,2X,F6.3,2X,F6.3)
    
    close(2)
    close(3)
    close(4)
    
    endprogram outline_output
        
        
            
            
            
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
