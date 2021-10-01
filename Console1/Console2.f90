    program Console1

    implicit none
    real ,allocatable :: xl(:),xr(:),ut(:),G_i(:),ro0_i(:),Y0_i(:),gu_i(:),a_i(:),la_i(:)
    real ,allocatable :: x_new(:),x_old(:),u_new(:),u_old(:),Ec_old(:),Ec_new(:), Pc_old(:), Pc_new(:), P_new(:), syg(:)
	real ,allocatable :: Ef_old(:), Ef_new(:), P_old(:), ro_new(:), e_xx(:), s_xx(:), s_yy(:), s_zz(:),mass(:),mk(:),dx(:)
	real ,allocatable :: q(:),V_new(:),ro_old(:),imp(:),impuls(:),u_con(:)
    integer ,allocatable :: Nb(:),Ne(:)
  

       
    
    integer :: N_body,Nk,i,number,tt,Nt,k,dtk,tk
    real :: ro0,G,Y0,gu,a,la,cf,c21,dt,dt1,t,times,u_n,dv,q0,r5,c_i,Uw,Dw,mindt,newdt,zy0,coef13,zj2,j
    real :: dEg,sumimp

	open(10,file='m.txt')
    open(11,file='km.txt')
	open(12,file='x-u-p.txt')
    open(13,file='imp.txt')
    
    N_body=3
    Nk=100
    Nt=1000
    times=30.0
    tk=0
    dtk=times/Nt
    k=1
    
    allocate (xl(N_body),xr(N_body),ut(N_body),Nb(N_body),Ne(N_body))

    xl(1)=0.0
    xr(1)=20.0
    xl(2)=20.1
    xr(2)=40.1
    xl(3)=40.2
    xr(3)=60.2
    
	ut(1)=0.05
	ut(2)=0.00  
	ut(3)=-0.05
    
    do i=1,N_body
         
        if (i == 1) then 
           Nb(i)=1
           Ne(i)=Nk
        else
            Nb(i)=Ne(i-1)+2
            Ne(i)=Nk+Nb(i)
        endif
    enddo


	allocate (G_i(Ne(N_body)+1),ro0_i(Ne(N_body)+1),Y0_i(Ne(N_body)+1),gu_i(Ne(N_body)+1),a_i(Ne(N_body)+1),la_i(Ne(N_body)+1))
	allocate (x_new(Ne(N_body)+1),x_old(Ne(N_body)+1),u_new(Ne(N_body)+1),u_old(Ne(N_body)+1))
	allocate (Ec_old(Ne(N_body)+1),Ec_new(Ne(N_body)+1), Pc_old(Ne(N_body)+1), Pc_new(Ne(N_body)+1), P_new(Ne(N_body)+1), syg(0:Ne(N_body)+1))
	allocate (Ef_old(Ne(N_body)+1), Ef_new(Ne(N_body)+1), P_old(Ne(N_body)+1), ro_new(Ne(N_body)+1), e_xx(Ne(N_body)+1), s_xx(Ne(N_body)+1), s_yy(Ne(N_body)+1), s_zz(Ne(N_body)+1))
	allocate (mass(Ne(N_body)+1),mk(Ne(N_body)+1),dx(N_body),q(0:Ne(N_body)+1),V_new(Ne(N_body)+1),ro_old(Ne(N_body)+1),imp(N_body),impuls(Nt),u_con(Nt))

	do i=Nb(1),Ne(1)
	number=1
	call material(number,ro0,G,Y0,gu,a,la)
	G_i(i)=G; ro0_i(i)=ro0; ro_old(i)=ro0; Y0_i(i)=Y0; gu_i(i)=gu; a_i(i)=a; la_i(i)=la
!print*, ro0_i(i)
	enddo
	
	do i=Nb(2),Ne(2)
	number=4
	call material(number,ro0,G,Y0,gu,a,la)
	G_i(i)=G; ro0_i(i)=ro0; Y0_i(i)=Y0; ro_old(i)=ro0; gu_i(i)=gu; a_i(i)=a; la_i(i)=la
    enddo
    
    do i=Nb(3),Ne(3)
	number=1
	call material(number,ro0,G,Y0,gu,a,la)
	G_i(i)=G; ro0_i(i)=ro0; Y0_i(i)=Y0; ro_old(i)=ro0; gu_i(i)=gu; a_i(i)=a; la_i(i)=la
    enddo
	
     
	 !do tt=1,N_body
		!x_old(Nb(tt)+1)=xl(tt)
		!
		!u_old(Ne(tt)+1)=ut(tt)
		!dx(tt)=(xr(tt)-xl(tt))/Nk        	
		!mass(Ne(tt)+1)=ro0_i(Ne(tt))*dx(tt) 
		!do i=Nb(tt),Ne(tt)
		!	x_old(i+1)=x_old(i)+dx(tt)
		!	u_old(i)=ut(tt)
		!	mass(i)=ro0_i(i)*dx(tt) 
		!enddo
	!enddo
	

    do tt=1,N_body
		x_old(Nb(tt))=xl(tt)
		x_old(Ne(tt)+1)=xr(tt)    
		dx(tt)=abs(xr(tt)-xl(tt))/Nk
		do i=Nb(tt),Ne(tt)
			x_old(i+1)=x_old(i)+dx(tt)
		enddo
		
		do i=Nb(tt),Ne(tt)+1
		   u_old(i)=ut(tt)
		enddo
		
		do i=Nb(tt),Ne(tt)
		   mass(i)=ro0_i(i)*dx(tt)
		enddo
		
    enddo
!print *,mass
    do i=1,Ne(N_body)+1
        x_new(i)=0.0; u_new(i)=0.0; Ec_old(i)=0.0; Ec_new(i)=0.0; Pc_old(i)=0.0; Pc_new(i)=0.0; P_new(i)=0.0; syg(i)=0.0;
        Ef_old(i)=0.0; Ef_new(i)=0.0; P_old(i)=0.0; ro_new(i)=0.0; e_xx(i)=0.0; s_xx(i)=0.0; s_yy(i)=0.0; s_zz(i)=0.0; mk(i)=0.0
    enddo
	   syg(0)=0.0; q(0)=0.0
    do tt=1,N_body
        do i=Nb(tt),Ne(tt)
            mk(i)=mk(i)+0.5*mass(i)
            mk(i+1)=mk(i+1)+0.5*mass(i)
        enddo
    enddo
!print*,mk

	do i=1,Ne(N_body)+1
        write(10,*)i,mass(i)   !mass.txt
        write(11,*)i,mk(i)     !kmas.txt
    enddo

	cf=0.25
	dt=100.0
	do tt=1,N_body
        do i=Nb(tt),Ne(tt)
		c21=a_i(i)!sqrt(-(a_i(i)**2*ro0_i(i)**2*(la_i(i)*(ro0_i(i)-ro0_i(i))+ro0_i(i)))/(la_i(i)*(ro0_i(i)-ro0_i(i))-ro0_i(i))**3); 
		 dt1=cf*dx(tt)/c21
           if (dt1 < dt) then
               dt=dt1 
           endif
		enddo		
	enddo
!print*,dt
    t=0.0
	q0=0.65
    !ops=1
    tk=1.0
    !dfzz=times/30.0

!print*,x_old
	do while (t < times)
        
        do tt=1,N_body
           do i=Nb(tt),Ne(tt)+1
               x_new(i)=x_old(i)+u_old(i)*dt
           enddo
		enddo
!print*,x_new
		do tt=2,N_body
            if (x_new(Ne(tt-1)+1) > x_new(Nb(tt))) then
!print*,'kkk'
               u_n=(mk(Ne(tt-1)+1)*u_old(Ne(tt-1)+1)+mk(Nb(tt))*u_old(Nb(tt)))/(mk(Ne(tt-1)+1)+mk(Nb(tt)));
               u_old(Nb(tt))=u_n
               u_old(Ne(tt-1)+1)=u_n
               x_new(Nb(tt))=x_old(Nb(tt))+dt*u_n
               x_new(Ne(tt-1)+1)=x_old(Ne(tt-1)+1)+dt*u_n
 !print*,u_n          
           endif
        enddo
		mindt=100.0

		do tt=1,N_body
           do i=Nb(tt),Ne(tt)
				dv=abs(x_new(i+1)-x_new(i))
                V_new(i)=dv
                ro_new(i)=mass(i)/V_new(i)
			   !ro_new(i)=r5
				Uw=a_i(i)*(ro0_i(i)-ro_new(i))/(la_i(i)*(ro_new(i)-ro0_i(i))-ro_new(i))
				Dw=a_i(i)+la_i(i)*Uw
				P_new(i)=ro0_i(i)*Dw*Uw
				c_i=sqrt(-(a_i(i)**2*ro0_i(i)**2*(la_i(i)*(ro_new(i)-ro0_i(i))+ro_new(i)))/(la_i(i)*(ro_new(i)-ro0_i(i))-ro_new(i))**3)
                q(i)=0.0
			   if (u_old(i+1)-u_old(i) < 0) then
                   q(i)=q0*ro_new(i)*c_i*ABS(u_old(i+1)-u_old(i))   
               endif
				newdt=cf*(V_new(i)/c_i)
               if (newdt < mindt) then
                   mindt=newdt
               endif
           enddo
        enddo
!print*,Pc_new        
		 dt=mindt
		 coef13=1./3.
		do tt=1,N_body
			do j=Nb(tt),Ne(tt)

                e_xx(j)=(u_old(j+1)-u_old(j))/V_new(j)
  
               s_xx(j)=s_xx(j)+dt*(4.0*G_i(j)*e_xx(j)*coef13)
               s_yy(j)=s_yy(j)-dt*(2.0*G_i(j)*e_xx(j)*coef13)
               s_zz(j)=s_zz(j)-dt*(2.0*G_i(j)*e_xx(j)*coef13)
	  
               zj2=s_xx(j)*s_xx(j)+s_yy(j)*s_yy(j)+s_zz(j)*s_zz(j)
               
               zy0=2.0*Y0_i(j)*Y0_i(j)*coef13
               if (zj2 > zy0) then
                   s_xx(j)=s_xx(j)*sqrt(zy0/zj2)
                   s_yy(j)=s_yy(j)*sqrt(zy0/zj2)
                   s_zz(j)=s_zz(j)*sqrt(zy0/zj2)
               endif
           enddo
			do j=Nb(tt),Ne(tt)
				Ec_new(j)=Ec_old(j)-0.5*(Pc_new(j)+Pc_old(j))*(1.0/ro_new(j)-1.0/ro_old(j))
			enddo

			do j=Nb(tt),Ne(tt)
				dEg=V_new(j)*s_xx(j)*e_xx(j)*dt
				Ef_new(j)=(Ef_old(j)-(0.5*(Pc_new(j)+P_old(j))+q(j))*(1.0/ro_new(j)-1.0/ro_old(j))+dEg)/(1+0.5*gu_i(j)*(1.0/ro_new(j)-1.0/ro_old(j)))
			enddo	
            do i=Nb(tt),Ne(tt)         
               Pc_new(i)=P_new(i)-gu_i(i)*(Ef_new(i)-Ec_new(i))*ro_new(i)			 
               syg(i)=-P_new(i)+s_xx(i)
              
            enddo
		  imp(tt)=0.0
        enddo       
!print*,Ef_new
          do j=1,Ne(N_body)+1
               u_new(j)=u_old(j)+dt*(syg(j)-q(j)-syg(j-1)+q(j-1))/mk(j)
          enddo
          
! Impuls
        do tt=1,N_body
                imp(tt)=0.0
            do i=Nb(tt),Ne(tt)+1
                imp(tt)=imp(tt)+mk(i)*u_new(i)
                
            enddo
        enddo
                sumimp=0.0
                impuls(k)=0.0
        If (t > tk) then
            do tt=1,N_body
                sumimp=sumimp+imp(tt)
 print*, sumimp
            enddo
                impuls(k)=sumimp
                k=k+1
                tk=tk+dtk
                u_con(k)=u_new(Ne(1)+1)
        endif 
        
        
        
 !print*,k       
      !enddo
          
			do i=1,Ne(N_body)+1 
				x_old(i)=x_new(i) 
				ro_old(i)=ro_new(i) 
				Pc_old(i)=Pc_new(i) 
				P_old(i)=P_new(i) 
				Ec_old(i)=Ec_new(i) 
				Ef_old(i)=Ef_new(i) 
				u_old(i)=u_new(i)

			enddo       
        
        t=t+dt	
!print*,t
	    enddo
	do i=1,Ne(N_body)+1
        write(12,*)i,x_new(i),u_new(i),P_new(i)*100,Pc_new(i)*100
    enddo
    !do i=1,N_body
    
    do i=1,k
        write(13,*)i,impuls(i),u_con(i)
    enddo
    

pause
    end program Console1


subroutine material(number,ro0,G,Y0,gu,a,la)
integer :: number
real :: ro0,G,Y0,gu,a,la

SELECT CASE (number)

   case (1) ! 1 = Al (Aluminium)
    ro0=2.7
    G=0.28; Y0=0.0023; gu=4.29553
    a=0.533; la=1.356   

    case (2) ! 4 = Fe
    z_r50=7.85
    z_c1=1.053854; z_c2=3.7575; z_c3=4.0431
    z_mu=0.793; z_y0=0.003; z_gu=1.9 !z_y0=0.0045
    ag=0.3664; lag=1.790
    
    case (3) ! 6 = Ti (TiTan)
    ro0=4.506
    G=0.335; Y0=0.003; gu=1.84
    a=0.4842; la=1.135

    case (4) ! 8 = Steel(12X18H10T)
    ro0=7.874
    G=0.793; Y0=0.003; gu=1.9 !z_y0=0.0045
    a=0.4599; la=1.388

   END SELECT
end
