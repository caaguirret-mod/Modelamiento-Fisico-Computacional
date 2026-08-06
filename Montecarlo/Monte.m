%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tiempo,E]=Monte(N)
%Esquema de montecarlo
%integracion de x.^2
%N=numero de elementos aleatorios
%Guardalo como Monte
%Llamarlo Monte(N),  N=Numero de datos estadisticos.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;     
a = 0;            
b = 1;            
real=pi./4;
for n=1:N
 x=rand(1,N);
  fx=1./(1+x.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ancho = b - a;
integral_mc = ancho * mean(fx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;
E(n)=real-integral_mc;
Tiempo(n)=toc;
end
subplot(2,2,1);
plot(Tiempo,'b--','LineWidth',2); grid 
subplot(2,2,2);
plot(E,'r--','LineWidth',2); grid
end