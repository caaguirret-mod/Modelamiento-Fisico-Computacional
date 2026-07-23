
%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=RK4(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ecuacion y´=y-x^2+1, x=[a,b] y(a)=0.5, condicion inicial
% para llamar la funcion, guardarla como: RK4
% para llamarla: RK4(N), N=cantidad de iteraciones
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.0; b=2.0;
f=@(x,y) y-x.^2+1;
%N=10;
h=(b-a)/N;
x(1)=a;
y(1) = 0.5;  %initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N
k1=h.*feval(f,x(k),y(k));
k2 = h.*feval(f,x(k)+0.5.*h,y(k)+0.5.*h.*k1);
k3 = h.*feval(f,x(k)+0.5.*h,y(k)+0.5.*h.*k2);
k4 = h.*feval(f,x(k)+h,y(k)+k3.*h);
y(k+1) = y(k) + (1/6)*(k1+2*k2+2*k3+k4);
x(k+1)=x(k)+h;
T(k)=abs(feval(f,x(k),y(k))-y(k));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=a:h:b;
y2=x2.^2+2.*x2+1-0.5.*exp(x2); %solucion exacta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1);
plot(y,'LineWidth',2); hold on; plot(y2,'LineWidth',2); grid
subplot(2,2,2);
plot(T,'LineWidth',2); grid
end 