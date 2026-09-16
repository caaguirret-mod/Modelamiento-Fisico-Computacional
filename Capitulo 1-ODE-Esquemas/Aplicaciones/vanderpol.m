function [mu]=vanderpol(mu)
%%%%%%%%%%%%%%%%%%%%%%%%
% Solucion de pendulo, via 2 esquemas RK2
%Pendulo de Van Der pol
% d^2x/dt^2+mu(1-x^2)dx/dt+x=0
%%%%%%%%%%%%%%%%%%%%%%%%
a=-40.0; b=20.0; N=400;
%mu=1.0;
f1=@(y) y;
f2=@(x,y) -x+mu*(1-x.^2)*y;

h=(b-a)/N;
x(1)=1;
y(1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N
%%%%%%%%%%%%%%%%%%%%%%%%%
k1=h.*feval(f1,y(k));
k2 = h*f1(y(k)+0.5*h*k1);
x(k+1) = x(k) + k2;
%%%%%%%%%%%%%%%%%%%%%%%%%
k3=h.*feval(f2,x(k),y(k));
k4 = h*f2(x(k)+0.5*h,y(k)+0.5*h*k3);
y(k+1) = y(k) + k4;
end
figure;
plot(x, y, 'b', 'LineWidth', 3); grid on
view(2)
end