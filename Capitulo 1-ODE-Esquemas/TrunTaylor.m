function [T]=TrunTaylor(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Errores de Truncamiento Taylor
%Ecuacion y´=y-x^2+1, x=[a,b] y(a)=0.5, condicion inicial
% para llamar la funcion, guardarla como: TrunTaylor 
% para llamarla: TrunTaylor(N), N=[0,5], hasta quinto orden.
%%%%%%%%%%%%%%%%%%%%%%%%%
a=0.0; b=3.0;
x(1)=a;
y(1)=0.5;
f=@(x,y) y-x^2+1;
n = 1000;  T(1)=0.0;
h = (b - a) / 1000;
%%%%%%%%%%%%%%%%%%%%%
%orden 2
%%%%%%%%%%%%%%%%%%%%%
if N==2
for k=1:n
y(k+1)=y(k)+(h).*((1+h./2).*(y(k)-x(k).^2+1)-h.*x(k));
x(k+1) = x(k) + h;
T(k)=abs(feval(f,x(k),y(k))-y(k+1));
end
%%%%%%%%%%%%%%%%%%%%%
%orden 3
%%%%%%%%%%%%%%%%%%%%%
elseif N==3
for k=1:n
y(k+1)=y(k)+h.*(y(k)-x(k).^2+1)+0.5.*(h.^2)*(y(k)-x(k).^2-2.*x(k)+1);
x(k+1) = x(k) + h;
T(k)=abs(feval(f,x(k),y(k))-y(k));
end
%%%%%%%%%%%%%%%%%%%%%
%orden 4
%%%%%%%%%%%%%%%%%%%%%
elseif N==4
for k=1:n
y(k+1)=y(k)+h.*(y(k)-x(k).^2+1)+0.5.*(h.^2)*(y(k)-x(k).^2-2.*x(k)+1)+...
    ((h.^3)/6).*(y(k)-x(k).^2-2.*x(k)+1);
x(k+1) = x(k) + h;
T(k)=abs(feval(f,x(k),y(k))-y(k));
end
%%%%%%%%%%%%%%%%%%%%%
%%Orden 5
%%%%%%%%%%%%%%%%%%%%%
elseif N==5
for k=1:n
y(k+1)=y(k)+h.*(y(k)-x(k).^2+1)+0.5.*(h.^2)*(y(k)-x(k).^2-2.*x(k)+1)+...
    ((h.^3)/6).*(y(k)-x(k).^2-2.*x(k)+1)+((h.^4)/24).*(y(k)-x(k).^2-2.*x(k)-1);
x(k+1) = x(k) + h;
T(k+1)=abs(feval(f,x(k),y(k))-y(k));
end 
%%%%%%%%%%%%%%%%%%%%%
    else 
    y=0;
%%%%%%%%%%%%%%%%%%%%%
end 
x2=a:h:b;
y2=x2.^2+2.*x2+1-0.5.*exp(x2); %solucion exacta
subplot(2,2,1);
plot(x,y,'b--o','LineWidth',2),hold on;plot(x2,y2,'r--','LineWidth',1);grid
subplot(2,2,2);
plot(T,'--','LineWidth',2); grid
end 