
program nastya
implicit none

 integer, allocatable :: Nb(:), Ne(:) 
    real, allocatable :: xl(:), xr(:), us(:), dx(:),imp(:),impuls(:,:),ufs(:),pfs(:)
    real, allocatable :: Gmat(:), Ymat(:),gumat(:),amat(:),lamat(:),m(:),mk(:),q(:),syg(:),ro0m(:),time_exp(:)
    real, allocatable :: ro_old(:),ro_new(:),x_new(:),x_old(:),u_new(:),u_old(:),Pc_new(:),Pc_old(:),Pf_new(:),Pf_old(:)
    real, allocatable :: e_xx(:),s_xx(:),s_yy(:),s_zz(:),Ec_new(:),Ec_old(:),Ef_new(:),Ef_old(:),V_new(:),V_old(:)
    
    integer :: n_body, k, N, i, N1, tt, number,j,ops,l
    real    :: time, coef13, ro0,G0,Y0,gu0,a0,la0,dt,dt1,t,q0,coef,u_n,mindt
    real    :: Uw, Dw,kkk1,kkk2,kkk3,c_i,newdt,zj2,zy0,dEg,dti,t1
    
    open(10,file='mat.txt')
    open(11,file='mass.txt')
    open(12,file='xu.txt')
    open(14,file='xupe.txt')
    open(15,file='imp.txt')
    open(16,file='ufs.txt')
    
    n_body=3
    time = 40.0
    coef = 0.25
    
    q0 = 0.75
    coef13=1.0/3.0
    ops = 1
    
    
    allocate (Nb(n_body), Ne(n_body))
    
    k=0
    N=100
    do i=1, n_body
        Nb(i)=1+k
        Ne(i)=N+k
    k=N+k+1;
    enddo
    
    allocate (xl(n_body), xr(n_body), us(n_body), dx(n_body),imp(n_body))
    
    xl(1)=0.0
    xr(1)=20.0
    xl(2)=20.1
    xr(2)=40.1
    xl(3)=40.2
    xr(3)=60.2
    us(1)=0.05
    us(2)=0.0
    us(3)=-0.05    
    
    do tt=1,n_body
        dx(tt)=abs(xr(tt)-xl(tt))/(Ne(tt)-Nb(tt)+1)
    enddo
    N1=Ne(n_body)+1
    
    allocate (ro_old(N1),Gmat(N1),Ymat(N1),gumat(N1),amat(N1),lamat(N1),ro0m(N1))
    
    !---material---!
    number = 1
    call material(number,ro0,G0,Y0,gu0,a0,la0)
    do i=Nb(1),Ne(1)
        ro_old(i)   = ro0
        ro0m(i)     = ro_old(i)
        Gmat(i)     = G0
        Ymat(i)     = Y0
        gumat(i)    = gu0
        amat(i)     = a0
        lamat(i)    = la0
    enddo
    
    number = 1
    call material(number,ro0,G0,Y0,gu0,a0,la0)
    do i=Nb(2),Ne(2)
        ro_old(i)   = ro0
        ro0m(i)     = ro_old(i)
        Gmat(i)     = G0
        Ymat(i)     = Y0
        gumat(i)    = gu0
        amat(i)     = a0
        lamat(i)    = la0
    enddo
    
    number = 1
    call material(number,ro0,G0,Y0,gu0,a0,la0)
    do i=Nb(3),Ne(3)
        ro_old(i)   = ro0
        ro0m(i)     = ro_old(i)
        Gmat(i)     = G0
        Ymat(i)     = Y0
        gumat(i)    = gu0
        amat(i)     = a0
        lamat(i)    = la0
    enddo
    
    do i=1,N1
        write(10,*) ro_old(i),Gmat(i), Ymat(i),gumat(i),amat(i),lamat(i)
    enddo
    !---material---!
    
    !---mass---!
    allocate (m(N1),mk(N1))
    do tt=1,n_body
        do i=Nb(tt),Ne(tt)
            m(i)=ro_old(i)*dx(tt)
            mk(i)=0.0
        enddo
        m(Ne(tt)+1)=0.0
        do i=Nb(tt),Ne(tt)
            mk(i)=mk(i)+0.5*m(i)
            mk(i+1)=mk(i+1)+0.5*m(i)
        enddo
    enddo
    
    do i=1,N1
        write(11,*) m(i),mk(i)
    enddo
    !---mass---!
    
    !---NULL---!
    allocate (x_new(N1),x_old(N1),u_new(N1),u_old(N1),V_new(N1),V_old(N1))
    allocate (Pc_new(N1),Pc_old(N1),Pf_new(N1),Pf_old(N1))
    allocate (Ec_new(N1),Ec_old(N1),Ef_new(N1),Ef_old(N1))!,impuls(ops))
    allocate (ro_new(N1),e_xx(N1),s_xx(0:N1),s_yy(0:N1),s_zz(0:N1),q(0:N1),syg(0:N1))
    do i=1,N1
    x_new(i)=0.0; x_old(i)=0.0
    u_new(i)=0.0; u_old(i)=0.0
    Pc_old(i)=0.0; Pc_new(i)=0.0
    Pf_old(i)=0.0; Pf_new(i)=0.0
    Ef_old(i)=0.0; Ef_new(i)=0.0
    Ec_old(i)=0.0; Ec_new(i)=0.0
    V_old(i)=0.0; V_new(i)=0.0
    ro_new(i)=0.0
    e_xx(i)=0.0
    s_xx(i)=0.0; s_yy(i)=0.0; s_zz(i)=0.0
    q(i)=0.0; syg(i)=0.0
    enddo
    s_xx(0)=0.0; s_yy(0)=0.0; s_zz(0)=0.0
    q(0)=0.0; syg(0)=0.0
    !---NULL---!
    
    !---SETKA---!
    do tt=1,n_body
        x_old(Nb(tt))=xl(tt)
        x_old(Ne(tt)+1)=xr(tt)
    enddo
    
    do tt=1,n_body
        do i=Nb(tt),Ne(tt)
            x_old(i+1)=x_old(i)+dx(tt)
        enddo
        do i=Nb(tt),Ne(tt)+1
            u_old(i)=us(tt)
        enddo
    enddo
    do i=1,N1
        write(12,*) x_old(i),u_old(i)
    enddo
    !---SETKA---!
    
    !---DT---!
    dt = 1.0
    do tt=1,n_body
        dt1=0.25*dx(tt)/amat(Nb(tt)) !!!!!!
        if (dt1 < dt) then
            dt = dt1
        endif
    enddo
    !write(*,*)dt        
    !---DT---!
    
    !---IMP---!
    dti = dt
    l = 0
    t1 = 0.0
    do while (t1 < time)
        t1 = t1 + dti
        l = l + 1
    enddo
    t1 = dti
    allocate (impuls(n_body+1,l),ufs(l),pfs(l),time_exp(l))
    !---------!
    
    
    !---WAVE---!
    t = 0.0
    do while (t < time)
        !--- coordinate ---!
        do tt=1,n_body
            do i=Nb(tt),Ne(tt)+1
               x_new(i)=x_old(i)+u_old(i)*dt
            enddo
        enddo
        !--- coordinate ---!
        
        !--- contact ---!
        do tt=2,n_body
            if (x_new(Ne(tt-1)+1) > x_new(Nb(tt))) then
                u_n=(mk(Ne(tt-1)+1)*u_old(Ne(tt-1)+1)+mk(Nb(tt))*u_old(Nb(tt)))/((mk(Ne(tt-1))+1)+mk(Nb(tt))) !!!!
                u_old(Nb(tt))=u_n
                u_old(Ne(tt-1)+1)=u_n
                x_new(Nb(tt))=x_old(Nb(tt))+dt*u_n
                x_new(Ne(tt-1)+1)=x_old(Ne(tt-1)+1)+dt*u_n;
            endif
        enddo          
        !--- contact ---!
        
        mindt = 100.0
        
        !--- Pressure ---!
        do tt=1,n_body
            do i=Nb(tt),Ne(tt)
                V_new(i)=x_new(i+1)-x_new(i)
                ro_new(i)=m(i)/ V_new(i)
                Uw=amat(i)*(ro0m(i)-ro_new(i))/(lamat(i)*(ro_new(i)-ro0m(i))-ro_new(i)) 
                Dw=amat(i)+lamat(i)*Uw
                Pf_new(i)=ro0m(i)*Dw*Uw 
                
                kkk1=-ro0m(i)*ro0m(i)*amat(i)*amat(i)
                kkk2=lamat(i)*ro0m(i)-lamat(i)*ro_new(i)-ro_new(i)
                kkk3=(-lamat(i)*ro_new(i)+lamat(i)*ro0m(i)+ro_new(i))**3
                c_i = sqrt(kkk1*kkk2/kkk3)
                newdt=coef*(V_new(i)/c_i)
                q(i)=0.0
                if (u_old(i+1)-u_old(i) < 0.0) then 
                    q(i)=q0*ro_new(i)*c_i*abs(u_old(i+1)-u_old(i))
                endif
                if ( newdt < mindt ) then 
                    mindt=newdt
                endif   
            enddo
        enddo
        !--- Pressure ---!
        dt = mindt
        
        !--- Deformation ---!
        do tt=1,n_body
            do j=Nb(tt),Ne(tt)
                e_xx(j)=(u_old(j+1)-u_old(j))/V_new(j)
                s_xx(j)=s_xx(j)+dt*(4.0*Gmat(j)*e_xx(j)*coef13)
                s_yy(j)=s_yy(j)-dt*(2.0*Gmat(j)*e_xx(j)*coef13)
                s_zz(j)=s_zz(j)-dt*(2.0*Gmat(j)*e_xx(j)*coef13)
                
                zj2=s_xx(j)*s_xx(j)+s_yy(j)*s_yy(j)+s_zz(j)*s_zz(j);
                zy0=2.0*Ymat(j)*Ymat(j)*coef13
                if (zj2 > zy0) then
                    s_xx(j)=s_xx(j)*sqrt(zy0/zj2)
                    s_yy(j)=s_yy(j)*sqrt(zy0/zj2)
                    s_zz(j)=s_zz(j)*sqrt(zy0/zj2)
                endif
            enddo
            !--- Deformation ---!
            
            !--- Energy C ---!
            do j=Nb(tt),Ne(tt)
                Ec_new(j)=Ec_old(j)-0.5*(Pc_new(j)+Pc_old(j))*(1.0/ro_new(j)-1.0/ro_old(j))
            enddo
            !--- Energy C ---!
            
            !--- Energy F ---!
            do j=Nb(tt),Ne(tt)
                dEg=V_new(j)*s_xx(j)*e_xx(j)*dt
                Ef_new(j)=(Ef_old(j)-(0.5*(Pc_new(j)+Pf_old(j))+q(j))*(1.0/ro_new(j)-1.0/ro_old(j))+dEg)/(1+0.5*gumat(j)*(1.0/ro_new(j)-1.0/ro_old(j)))
            enddo
            !--- Energy F ---!
            
            !--- Pressure F ---!
            do i=Nb(tt),Ne(tt)
                Pc_new(i)=Pf_new(i)-gumat(i)*(Ef_new(i)-Ec_new(i))
                syg(i)=-Pf_new(i)+s_xx(i)
            enddo
            !--- Pressure F ---!
            
            !--- Speed ---!
            do j=Nb(tt),Ne(tt)+1
                u_new(j)=u_old(j)+dt*(syg(j)-q(j)-syg(j-1)+q(j-1))/mk(j)
            enddo
            !--- Speed ---!    
        enddo
        
        !--- IMP ---!
        if (t1 < t) then
            do tt=1,n_body
                imp(tt)=0.0
                do j=Nb(tt),Ne(tt)
                    imp(tt)=imp(tt)+mk(j)*u_new(j)
                enddo
                impuls(tt,ops)=impuls(tt,ops)+imp(tt)
                impuls(n_body+1,ops)=impuls(n_body+1,ops)+impuls(tt,ops)
                ufs(ops)=u_new(Ne(n_body)+1)
                pfs(ops)=Pf_new(Ne(n_body))
                time_exp(ops)=t
            enddo
            t1 = t1 + dti
            ops = ops +1
        endif
        !--- IMP ---!
        
        !---old time massive---!
        do i=1,N1
            x_old(i)=x_new(i)
            ro_old(i)=ro_new(i)
            Pc_old(i)=Pc_new(i)
            Pf_old(i)=Pf_new(i)
            Ec_old(i)=Ec_new(i)
            Ef_old(i)=Ef_new(i)
            u_old(i)=u_new(i)
        enddo
        

    t = t + dt
    write(*,*)t
    enddo
    !---WAVE---!
    do i=1,N1
        write(14,*) x_new(i),u_new(i),Pf_new(i),Ef_new(i) ! plot "xupe.txt" u 1:2 with li
    enddo
    
    do i=1,l-1
        write(15,*) time_exp(i),impuls(1,i),impuls(2,i),impuls(3,i),impuls(4,i) ! plot "imp.txt" u 1:2 with li, "imp.txt" u 1:3 with li, "imp.txt" u 1:4 with li, "imp.txt" u 1:5 with li
        write(16,*) time_exp(i),ufs(i),pfs(i)
    enddo
pause
end program nastya

    
subroutine material(number,ro0,G0,Y0,gu0,a0,la0)
   integer :: number
   real :: ro0,G0,Y0,gu0,a0,la0
   
   SELECT CASE (number)
   case (1) ! 1 = Al (Aluminium)
       ro0=2.7
       G0=0.28; Y0=0.005; gu0=4.29553
       a0=0.533; la0=1.356
    
   case (2) ! 4 = Fe
       ro0=7.85
       G0=0.793; Y0=0.003; gu0=1.9
       a0=0.3664; la0=1.790
   case (3) ! 9 = WC0
       ro0=15.66
       G0=0.00301; Y0=0.00953; gu0=1.5 
       a0=0.493; la0=1.309
   END SELECT
end 
