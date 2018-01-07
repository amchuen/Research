% project 3--upwind-Zha-Bilgen flux splitting
clc;
clear
close all;

nx=256; ny=64; xl=4.0;yl=1; tfinal=0.2; time=0; gg=1.4;
p_left=116.5;p_right=1;
r_left=8;r_right=1.4;
u_left=8.25*cos(pi/6); u_right=0.0;
v_left=-8.25*sin(pi/6); v_right=0.0;

r=zeros(nx,ny);p=zeros(nx,ny);rE=zeros(nx,ny); c=zeros(nx,ny);
ru=zeros(nx,ny);u=zeros(nx,ny);mx=zeros(nx,ny);
rv=zeros(nx,ny);v=zeros(nx,ny);my=zeros(nx,ny);
F1=zeros(nx,ny);F2=zeros(nx,ny);F3=zeros(nx,ny);F4=zeros(nx,ny);
G1=zeros(nx,ny);G2=zeros(nx,ny);G3=zeros(nx,ny);G4=zeros(nx,ny);
h=xl/(nx-1);

for i=1:nx
    for j=1:ny
        x(i,j)=h*(i-1);y(i,j)=h*(j-1);
    end
end

for i=1:nx
    for j=1:ny
        r(i,j)=r_right;
        r(i,j)=r_right;
        ru(i,j)=0.0; 
        rv(i,j)=0; 
        rE(i,j)=p_right/(gg-1); 
    end
end

xo=1/6;
for i=1:nx
    for j=1:ny
        if x(i,j) < xo+(y(i,j)/sqrt(3))
            r(i,j)=r_left; rE(i,j)=p_left/(gg-1)+0.5*r_left*u_left^2;
            ru(i,j)=r_left*u_left; rv(i,j)=r_left*v_left;
        end
    end
end

cmax=sqrt( max(gg*p_right/r_right,gg*p_left/r_left) );
dt=0.145*h/cmax; maxstep=tfinal/dt;

for istep=1:maxstep
    
    for i=1:nx
        for j=1:ny
            p(i,j)=(gg-1)*(rE(i,j)-0.5*((ru(i,j)*ru(i,j)+rv(i,j)*rv(i,j))/r(i,j)));
        end
    end
    
    for i=1:nx
        for j=1:ny,c(i,j)=sqrt( gg*p(i,j)/r(i,j) );
        end
    end
    
    for i=1:nx
        for j=1:ny
            u(i,j)=ru(i,j)/r(i,j);v(i,j)=rv(i,j)/r(i,j);
        end
    end
    
    for i=1:nx
        for j=1:ny
            mx(i,j)=u(i,j)/c(i,j);my(i,j)=v(i,j)/c(i,j);
        end
    end
    
    for i=1:nx-1 
        for j=1:ny-1 % Find fluxes
            F1(i,j)=0.5*(ru(i+1,j)+ru(i,j))-0.5*(abs(ru(i+1,j))-abs(ru(i,j)));
            F2(i,j)=0.5*(u(i+1,j)*ru(i+1,j)+p(i+1,j)+u(i,j)*ru(i,j)+p(i,j))...
                -0.5*(abs(u(i+1,j))*ru(i+1,j)-abs(u(i,j))*ru(i,j))...
                -0.5*(p(i+1,j)*mx(i+1,j)-p(i,j)*mx(i,j));
            F3(i,j)=0.5*(u(i+1,j)*rv(i+1,j)+u(i,j)*rv(i,j))...
                -0.5*(abs(u(i+1,j))*rv(i+1,j)-abs(u(i,j))*rv(i,j));
            F4(i,j)=0.5*(u(i+1,j)*(rE(i+1,j)+p(i+1,j))+u(i,j)*(rE(i,j)+p(i,j)))...
                -0.5*(abs(u(i+1,j))*rE(i+1,j)-abs(u(i,j))*rE(i,j))...
                -0.5*(p(i+1,j)*c(i+1,j)-p(i,j)*c(i,j));
            
            if mx(i,j) > 1
                F2(i,j)=ru(i,j)*u(i,j)+p(i,j);
                F4(i,j)=(rE(i,j)+p(i,j))*u(i,j);
            end
            
            if mx(i,j) < -1
                F2(i,j)=ru(i+1,j)*u(i+1,j)+p(i+1,j);
                F4(i,j)=(rE(i+1,j)+p(i+1,j))*u(i+1,j);
            end
            
            G1(i,j)=0.5*(rv(i,j+1)+rv(i,j))-0.5*(abs(rv(i,j+1))-abs(rv(i,j)));
            G2(i,j)=0.5*(v(i,j+1)*ru(i,j+1)+v(i,j)*ru(i,j))...
                -0.5*(abs(v(i,j+1))*ru(i,j+1)-abs(v(i,j))*ru(i,j));
            G3(i,j)=0.5*(v(i,j+1)*rv(i,j+1)+p(i,j+1)+v(i,j)*rv(i,j)+p(i,j))...
                -0.5*(abs(v(i,j+1))*rv(i,j+1)-abs(v(i,j))*rv(i,j))...
                -0.5*(p(i,j+1)*my(i,j+1)-p(i,j)*my(i,j));
            G4(i,j)=0.5*(v(i,j+1)*(rE(i,j+1)+p(i,j+1))+v(i,j)*(rE(i,j)+p(i,j)))...
                -0.5*(abs(v(i,j+1))*rE(i,j+1)-abs(v(i,j))*rE(i,j))...
                -0.5*(p(i,j+1)*c(i,j+1)-p(i,j)*c(i,j));
            
            if my(i,j) > 1
                G3(i,j)=rv(i,j)*v(i,j)+p(i,j);
                G4(i,j)=(rE(i,j)+p(i,j))*v(i,j);
            end
            
            if my(i,j) < -1
                G3(i,j)=rv(i,j+1)*v(i,j+1)+p(i,j+1);
                G4(i,j)=(rE(i,j+1)+p(i,j+1))*v(i,j+1);
            end
        end
    end
    
    for i=2:nx-1 
        for j=2:ny-1 % Update solution
            r(i,j) =r(i,j) -(dt/h)*(F1(i,j)-F1(i-1,j)+G1(i,j)-G1(i,j-1));
            ru(i,j)=ru(i,j)-(dt/h)*(F2(i,j)-F2(i-1,j)+G2(i,j)-G2(i,j-1));
            rv(i,j)=rv(i,j)-(dt/h)*(F3(i,j)-F3(i-1,j)+G3(i,j)-G3(i,j-1));
            rE(i,j)=rE(i,j)-(dt/h)*(F4(i,j)-F4(i-1,j)+G4(i,j)-G4(i,j-1));
        end
    end
    
    % Bottom
    no=floor(nx*xo/xl)+2;
    r(no:nx-2,1)=r(no:nx-2,2);ru(no:nx-2,1)=ru(no:nx-2,2);rE(no:nx-2,1)=rE(no:nx-2,2);
    
    %Top
    xs=xo+((1+20*time))/sqrt(3); ns=floor(nx*xs/xl)+1;
    
    for i=2:ns
        r(i,ny)=r_left; ru(i,ny)=r_left*u_left; rv(i,ny)=r_left*v_left;
        rE(i,ny)=p_left/(gg-1)+0.5*(ru(i,ny)*ru(i,ny)/r(i,ny));
    end
    
    % r(1,2:ny-2)=r(2,2:ny-2);rv(1,2:ny-2)=rv(2,2:ny-2);rE(1,2:ny-2)=rE(2,2:ny-2);
    r(nx,2:ny-2)=r(nx-1,2:ny-2);rv(nx,2:ny-2)=rv(nx-1,2:ny-2);rE(nx,2:ny-2)=rE(nx-1,2:ny-2);
    time=time+dt; 
    disp(istep);
    contour(x,y,r,40); 
    axis equal; 
    axis([0, xl, 0, yl]);
    drawnow;
    pause(0.001)
end

% plot(x,r,'k','linewidth',2);hold on
% set(gca,'Box','on'); set(gca,'Fontsize',24, 'LineWidth',2)
% text(5,0.9,'Density','Fontsize',24)
% print -depsc ZBResults1 mesh(x,y,r) quiver(x,y,u,v)

axis([0, xl*3/4, 0, yl])

%Euler equations--one-dimensional Lax-Friedrich with fluxes

nx=256; maxstep=150;
gg=1.4;p_left=100000;p_right=10000;r_left=1;r_right=0.125;u_left=0;
xl=10.0;h=xl/(nx-1);time=0;

h=xl/(nx-1);time=0;
for i=1:nx
    x(i)=h*(i-1);
end

r=zeros(1,nx);ru=zeros(1,nx);rE=zeros(1,nx);
rn=zeros(1,nx);run=zeros(1,nx);rEn=zeros(1,nx);
c=zeros(1,nx);u=zeros(1,nx);m=zeros(1,nx);

for i=1:nx
    r(i)=r_right;ru(i)=0.0;rE(i)=p_right/(gg-1);
end

for i=1:nx/2
    r(i)=r_left;
    rE(i)=p_left/(gg-1)+0.5*r_left*u_left^2;
    ru(i)=r_left*u_left;
end

cmax=sqrt( max(gg*p_right/r_right,gg*p_left/r_left) ); dt=0.45*h/cmax;
for istep=1:maxstep
    for j=1:nx-1
        rL(j)=r(j);ruL(j)=ru(j);rEL(j)=rE(j); rR(j)=r(j+1);ruR(j)=ru(j+1);rER(j)=rE(j+1);
    end
    rR(nx)=r(nx);ruR(nx)=ru(nx);rER(nx)=rE(nx); % do endpoints
    
    for j=1:nx-1 % find the fluxes
        pL=(gg-1)*(rEL(j)-0.5*(ruL(j)*ruL(j)/rL(j)));
        pR=(gg-1)*(rER(j)-0.5*(ruR(j)*ruR(j)/rR(j)));
        F1(j)=0.5*(ruR(j)+ruL(j))-0.5*(h/dt)*(rR(j)-rL(j));
        F2(j)=0.5*((ruR(j)^2/rR(j))+pR+(ruL(j)^2/rL(j))+pL)-0.5*(h/dt)*(ruR(j)-ruL(j));
        F3(j)=0.5*((ruR(j)/rR(j))*(rER(j)+pR)+(ruL(j)/rL(j))*(rEL(j)+pL))-...
        0.5*(h/dt)*(rER(j)-rEL(j));
    end
    
    for j=2:nx-1 % update the solution
        r(j)=r(j)-(dt/h)*(F1(j)-F1(j-1));
        ru(j)=ru(j)-(dt/h)*(F2(j)-F2(j-1));
        rE(j)=rE(j)-(dt/h)*(F3(j)-F3(j-1));
    end
    
        % do endpoints
        r(1)=r(2); r(nx)=r(nx-4);r(nx-1)=r(nx-4);r(nx-2)=r(nx-4);r(nx-3)=r(nx-4);
        ru(1)=ru(2); ru(nx)=ru(nx-4);ru(nx-1)=ru(nx-4);ru(nx-2)=ru(nx-4);ru(nx-3)=ru(nx-4);
        rE(1)=rE(2); rE(nx)=rE(nx-4);rE(nx-1)=rE(nx-4);rE(nx-2)=rE(nx-4);rE(nx-3)=rE(nx-4);
        % plot(r,'r','linewidth',2); pause
        time=time+dt; %istep
end

plot(r,'g','linewidth',2);% one-dimensional linear advection by WENO

% fifth order scheme. Comparison with first order upwind
%------------------------------------------------------------
n=101; nstep=260; length=2.0;h=length/(n-1);
dt=0.3*h; 
for i=1:n+1
    x(i)=h*(i-1)
end

time=0.0;
gam1=1/10; gam2=3/5; gam3=3/10; eps=0.000001;
f=zeros(n+3,1); fo=zeros(n+3,1); f1=zeros(n+3,1); f2=zeros(n+3,1);
flux=zeros(n+1,1);
fu=zeros(n+3,1);fuo=zeros(n+3,1);

%for i=1:n, f(i)=1+0.5*sin(2*pi*h*(i-1)); end; %initial conditions
for i=10:20
    f(i)=1;
end
for i=1:n+1
    if (x(i)>=0.0 & x(i)<=0.6)
        f(i)=exp(-100*(x(i)-0.3)^2); 
    end
end
for i=1:n+1
    if (x(i)>=0.6 & x(i)<=0.8)
        f(i)=1.0;
    end
end %initial conditions
fu=f;

for m=1:nstep
%         m
%         time
        hold off;plot(f,'linewidt',2); %axis([1 n 2.0, 4.5]); %plot solution
        hold on;plot(fu,'r','linewidt',2);
        set(gcf,'Position',[343 385 500 350]);axis([0 n 0 2]);
        pause(0.1)

    %-----------------------
        fo=f;
    for j=3:n+1
        beta1=(13/12)*(f(j-2)^2-2*f(j-1)+f(j))^2 +(1/4)*(f(j-2)-4*f(j-1)+3*f(j))^2;
        beta2=(13/12)*(f(j-1)^2-2*f(j) +f(j+1))^2+(1/4)*(f(j-1)-f(j+1))^2;
        beta3=(13/12)*(f(j)^2 -2*f(j+1)+f(j+2))^2+(1/4)*(3*f(j)-4*f(j+1)+f(j+2))^2;
        omt1=gam1/(eps+beta1)^2; omt2=gam2/(eps+beta2)^2;omt3=gam3/(eps+beta3)^2;
        omtsum=omt1+omt2+omt3; om1=omt1/omtsum; om2=omt2/omtsum; om3=omt3/omtsum;
        flux(j)=om1*((1/3)*f(j-2)-(7/6)*f(j-1)+(11/6)*f(j))+...
                om2*(-(1/6)*f(j-1)+(5/6)*f(j)+(1/3)*f(j+1))+...
                om3*((1/3)*f(j)+(5/6)*f(j+1)-(1/6)*f(j+2));
    end
    flux(1)=flux(n); flux(2)=flux(n+1);
    for j=2:n
        f(j)=fo(j)-(dt/h)*(flux(j)-flux(j-1));
    end;
    f(1)=f(n);f(n+1)=f(2);

    % Second Step
    f1=f;
    for j=3:n+1
        beta1=(13/12)*(f(j-2)^2-2*f(j-1)+f(j))^2 +(1/4)*(f(j-2)-4*f(j-1)+3*f(j))^2;
        beta2=(13/12)*(f(j-1)^2-2*f(j) +f(j+1))^2+(1/4)*(f(j-1)-f(j+1))^2;
        beta3=(13/12)*(f(j)^2 -2*f(j+1)+f(j+2))^2+(1/4)*(3*f(j)-4*f(j+1)+f(j+2))^2;
        omt1=(gam1)/(eps+beta1)^2; omt2=(gam2)/(eps+beta2)^2;omt3=(gam3)/(eps+beta3)^2;
        omtsum=omt1+omt2+omt3; om1=omt1/omtsum; om2=omt2/omtsum; om3=omt3/omtsum;
        flux(j)=om1*((1/3)*f(j-2)-(7/6)*f(j-1)+(11/6)*f(j))+...
        om2*(-(1/6)*f(j-1)+(5/6)*f(j)+(1/3)*f(j+1))+...
        om3*((1/3)*f(j)+(5/6)*f(j+1)-(1/6)*f(j+2));
    end
    flux(1)=flux(n); flux(2)=flux(n+1);
    for j=2:n
        f(j)=0.75*fo(j)+0.25*f1(j)-0.25*(dt/h)*(flux(j)-flux(j-1));
    end
    f(1)=f(n);f(n+1)=f(2);
    % Third Step
    f2=f;
    for j=3:n+1
        beta1=(13/12)*(f(j-2)^2-2*f(j-1)+f(j))^2 +(1/4)*(f(j-2)-4*f(j-1)+3*f(j))^2;
        beta2=(13/12)*(f(j-1)^2-2*f(j) +f(j+1))^2+(1/4)*(f(j-1)-f(j+1))^2;
        beta3=(13/12)*(f(j)^2 -2*f(j+1)+f(j+2))^2+(1/4)*(3*f(j)-4*f(j+1)+f(j+2))^2;
        omt1=(gam1)/(eps+beta1)^2; omt2=(gam2)/(eps+beta2)^2;omt3=(gam3)/(eps+beta3)^2;
        omtsum=omt1+omt2+omt3; om1=omt1/omtsum; om2=omt2/omtsum; om3=omt3/omtsum;
        flux(j)=om1*((1/3)*f(j-2)-(7/6)*f(j-1)+(11/6)*f(j))+...
            om2*(-(1/6)*f(j-1)+(5/6)*f(j)+(1/3)*f(j+1))+...
            om3*((1/3)*f(j)+(5/6)*f(j+1)-(1/6)*f(j+2));
    end
    flux(1)=flux(n); flux(2)=flux(n+1);
    for j=2:n
        f(j)=(1/3)*fo(j)+(2/3)*f2(j)-(2/3)*(dt/h)*(flux(j)-flux(j-1));
    end;
    f(1)=f(n);f(n+1)=f(2);
    %-----------------------
    % upwind
    fuo=fu; 
    for i=2:n-1
        fu(i)=fuo(i)-(dt/h)*(fuo(i)-fuo(i-1)); 
    end
    fu(1)=fuo(1)-(dt/h)*(fuo(1)-fuo(n-1));fu(n)=fu(1);
    time=time+dt;
end