%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [integral]=Mon(N)
% integracion de Monte carlo
% se genera estadisticamente la funcion f(x)=1./(1+x^2)
% calcula sigma, tiempo y error
% Guardar, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;     
a = 0;            
b = 1;    
ancho = b - a;
real=pi./4;
N=1000;
for n=1:N
 x=rand(1,N);
  fx=1./(1+x.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integral = ancho * mean(fx);
fx2=fx.*fx;
sigma(n)=mean(fx2)-mean(fx).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;
integralM(n)=integral;
E(n)=real-integral;
Tiempo(n)=toc;
fx=0.0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 
subplot(2,2,1);
plot(Tiempo,'b--','LineWidth',2); grid 
subplot(2,2,2);
plot(sigma,'r--','LineWidth',2); grid
subplot(2,2,3);
plot(sigma,'g--','LineWidth',2); grid
subplot(2,2,4);
plot(integralM,'g--','LineWidth',2); grid
