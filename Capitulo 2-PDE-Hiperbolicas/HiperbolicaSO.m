function [r2]=HiperbolicaSO(J,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Esquema explicito Hiperbolicas
%Ecuacion uxx=autt, esquema implicito
%Esquema=  u(j,n+1)=2(1-r^2)u(j,n)-u(j,n-1)+(r^2)(u(j+1,n)-u(j-1,n));
% para llamar la funcion, guardarla como: Hiperbolic2
% para llamarla: Hiperbolic2(J), J=malla espacial
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=0; x2=1; %J=40; 
Nt = 50;
dx = (x2-x1) / J;
tf = 0.5;
dt = tf/Nt;
r = dt/(dx);
%a=1.0;
r2=a*r;
c=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x1:dx:x2;
f=exp(-c*(x-0.2).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(J+1,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:Nt %ok
   t=n*dt;
   g1 =exp(-c*(x1-t-0.2).^2); 
   g2 =exp(-c*(x2-t-0.2).^2);
 if n==1
  for j=2:J
    u(j,n) = 2*(1-r.^2)*f(j)-f(j)+(r.^2)*(f(j+1)+f(j-1));  
  end
 elseif n==2
       for j=2:J
      u(j,n)=2*(1-r.^2)*u(j,n-1)+(r.^2)*(u(j+1,n-1)+u(j-1,n-1));   
   end
 else 
    for j=2:J
    u(j,n)=2*(1-r.^2)*u(j,n-1)-u(j,n-2)+(r.^2)*(u(j+1,n-1)+u(j-1,n-1));
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