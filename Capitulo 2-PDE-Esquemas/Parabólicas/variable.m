%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u]=variable(J)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Esquema explicito 
%Ecuacion uxx=ut, x=[0,1], u(1,n)=f(x), u(J+1,n)=g(x) y condicion inicial
%u(x,0)=sin(pix)+sin(2pix);
%estencil = u_{j}^{n}=u_{j}^{n}+r(a(j*1/2))^{n})(u_{j+1}^{n}+u_{j}^{n})+r(a(j-1/2)(u_{j}^{n}+u_{j-1}^{n}))
% para llamar la funcion, guardarla como: ParabolicaL
% para llamarla: Parabolico(J), J=malla espacial
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=0; x2=1; %J=10;
Nt = 50;
dx = (x2-x1) / J; tf = 0.1;
dt = tf/Nt;
r = dt/(dx)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x1:dx:x2;
a=@(x,t) sin(x)*cos(t);
f=sin(pi.*x)+sin(2*pi*x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(J+1,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:Nt 
   t=n*dt;
   g1 =sin(pi*x1)*exp(-pi*pi*t)+sin(2*pi*x1)*exp(-4*pi*pi*t); 
   g2 = sin(pi*x2)*exp(-pi*pi*t)+sin(2*pi*x2)*exp(-4*pi*pi*t);
 if n==1
  for j=2:J
    u(j,n) = f(j) + r*(f(j+1)-2*f(j)+f(j-1));  
    end
 else 
     for j=2:J
     u(j,n)=u(j,n-1)+r*(feval(a,j+1/2,n-1)*(u(j+1,n-1)-u(j,n-1)))...
         +r*(feval(a,j-1/2,n)*(u(j,n-1)-u(j-1,n-1)));
 end 
    u(1,n)=g1;
    u(J+1,n)=g2;
 end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = dt : dt : Nt*dt;
figure(1)
colormap(jet);
surf(x,tt, u');
xlabel('x')
ylabel('t')
zlabel('u')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 