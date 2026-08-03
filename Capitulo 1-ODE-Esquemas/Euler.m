%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u]=Euler(N)
%Esquema Euler
%Ecuacion dy/dx=f(x,y)
% para llamar la funcion, guardarla como: Euler
% para llamarla: Euler(N), N=Numero de pasos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.0; b=2.0;
f=@(x,y) y-x.^2+1;
%N=10;
x(1)=0; y(1)=0.5;
h=(b-a)/N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N
y(k+1)=y(k)+h*feval(f,x(k),y(k));
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