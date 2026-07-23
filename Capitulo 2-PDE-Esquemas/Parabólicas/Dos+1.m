%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [u]=DosR(J)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Esquema explicito 
%Ecuacion (u)xx+(u)yy=ut, x=[0,1],y=[0,1], u(1,k,n)=f(x,y), u(J+1,k,n)=g(x,y),u(j,1,n)=f(x,y), u(J+1,1,n)=f(x,y) y condicion inicial
%u(x,y,0)=sin(pixy)+sin(2pixy);
%estencil =
%u_{j,k}^{n+1}=u_{j,k}^{n}+r(u_{j+1,k}^{n}-2u_{j,k}^{n}+u_{j-1,k}^{n})+r(u_{j,k+1}^{n}-2u_{j,k}^{n}+u_{j,k-1}^{n})
% para llamar la funcion, guardarla como: Nolineales
% para llamarla: Nolineales(J), J=malla espacial
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=0; x2=1; %J=10;
y1=0; y2=1;
Nt = 50;
dx = (x2-x1) / J;
dy = (y2-y1) / J; 
tf = 0.1;
dt = tf/Nt;
r = dt/(dx)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x1:dx:x2;
y=y1:dy:y2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=sin(pi.*x.*y)+sin(2*pi*x.*y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(J+1,J+1,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:Nt 
   t=n*dt;
   g1 =sin(pi*x1*y1)*exp(-pi*pi*t)+sin(2*pi*x1*y1)*exp(-4*pi*pi*t); 
   g2 = sin(pi*x2*y2)*exp(-pi*pi*t)+sin(2*pi*x2*y2)*exp(-4*pi*pi*t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if n==1
   for k=2:J
     for j=2:J
       u(j,k,n) = f(j) + r*(f(j+1)-2*f(j)+f(j-1));  
     end
   end 
 else 
     for k=2:J
    for j=2:J
  u(j,k,n)=u(j,k,n-1)+r*(u(j+1,k,n-1)-2*u(j,k,n-1)+u(j-1,k,n-1))...
      +r*(u(j,k+1,n-1)-2*u(j,k,n-1)+u(j,k-1,n-1));       
    end 
     end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for k=2:J
    for j=2:J
    u(1,k,n)=g1;
    u(J+1,k,n)=g2;
    u(j,1,n)=g1;
    u(j,J+1,n)=g2;
    end 
     end   
 end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%end 
tt = dt : dt : Nt*dt;
subplot(2,2,1)
contour(u(:,:,1))

subplot(2,2,2)
contour(u(:,:,10))

subplot(2,2,3)
contour(u(:,:,20))

subplot(2,2,4)
contourf(u(:,:,40))
