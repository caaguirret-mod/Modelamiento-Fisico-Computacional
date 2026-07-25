function [r2]=Hiperbolic2(J,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Esquema explicito Hiperbolicas
%Ecuacion ux+aut=0, esquema implicito
%Esquema=  u(j,n+1)=u(j,n)+a*r(u_{j-1,n}-u_{j}^{n}));
% para llamar la funcion, guardarla como: Hiperbolic2
% para llamarla: Hiperbolic2(J), J=malla espacial
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=0; x2=1; %J=40; 
Nt = 50;
dx = (x2-x1) / J;
tf = 0.5;
dt = tf/Nt;
r = dt/(dx);
r2=a*r;
c=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x1:dx:x2;
f=exp(-c*(x-0.2).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(J+1,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:Nt 
   t=n*dt;
   g1 =exp(-c*(x1-t-0.2).^2); 
   g2 =exp(-c*(x2-t-0.2).^2);
 if n==1
  for j=2:J
    u(j,n) = f(j)-a*r*(f(j-1)-f(j));  
    end
 else 
     for j=2:J
     u(j,n)=u(j,n-1)-a*r*(u(j-1,n-1)-u(j,n-1));
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
