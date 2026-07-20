%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=ADB(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Metodo de Adam-Belford
%Ecuacion y´=y-x^2+1, x=[a,b] y(a)=0.5, condicion inicial
% para llamar la funcion, guardarla como: ADB
% para llamarla: ADB(N), N=cantidad de iteraciones
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.0; b=2.0;
f=@(x,y) y-x.^2+1;
%N=10;
h=(b-a)/N;
x(1)=a;x(2)=x(1)+h;
y(1) = 0.5; y(2) = 0.2; y(3)=0.25; y(4)=0.1; %alphas

k1 = feval(f,x(1),y(1));
k2 = feval(f,x(1)+h, y(1)+h.*k1);
y(2) = y(1)+ (h./2).*(k1 + k2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:N
   f_n = feval(f,x(k), y(k));          
    f_nm1 = feval(f,x(k-1), y(k-1));    
    y(k+1) = y(k) + (h./2).*(3.*f_n-f_nm1);
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