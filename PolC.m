%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=PolC(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solucion con interpolacion lineal
%Ecuacion y´=y-x^2+1, x=[a,b] y(a)=0.5, condicion inicial
% para llamar la funcion, guardarla como: PolC
% para llamarla: PolC(N), N=cantidad de iteraciones
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.0; b=2.0;
f=@(x,y) y-x.^2+1;
%N=100;
h=(b-a)/N;
x(1)=a;x(2)=x(1)+h;
y(1) = 0.5; y(2) = 0.2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1 = feval(f,x(1),y(1));
k2 = feval(f,x(1)+h, y(1)+h.*k1);
y(2) = y(1)+ (h./2).*(k1 + k2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:N
y(k+1) = y(k)+(1./(2.*h)).*((x(k-1)-x(k)).^2).*(feval(f,x(k-1),y(k-1))+feval(f,x(k),y(k)));
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