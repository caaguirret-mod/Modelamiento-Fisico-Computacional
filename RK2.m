%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=RK2(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ecuacion y´=y-x^2+1, x=[a,b] y(a)=0.5, condicion inicial
% para llamar la funcion, guardarla como: RK2 
% para llamarla: RK2(N), N=cantidad de iteraciones
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.0; b=2.0;
f=@(x,y) y-x.^2+1;
h=(b-a)/N;
x(1)=a;
y(1)=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N
k1=h.*feval(f,x(k),y(k));
k2 = h*f(x(k)+0.5*h,y(k)+0.5*h*k1);
y(k+1) = y(k) + k2;
x(k+1)=x(k)+h;
T(k)=abs(feval(f,x(k),y(k))-y(k));
end
%%%%%%%%%%%%%%%%%%%%%%%%%e
x2=a:h:b;
y2=x2.^2+2.*x2+1-0.5.*exp(x2); %solucion exacta
subplot(2,2,1);
plot(y,'--','LineWidth',2);hold on; plot(y2,'r--','LineWidth',2); grid
subplot(2,2,2);
plot(T,'--o','LineWidth',2); grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%e
end 