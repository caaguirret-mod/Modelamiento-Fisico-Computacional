close all; clear all; clc;
a=0.0; b=200.0; N=400;
alpha=1;
delta=2;
beta=8/3;
f1=@(x,y,z) -alpha.*(x-y);
f2=@(x,y,z) (delta+1).*x-y-x.*z;
f3=@(x,y,z) -beta*z+x.*y;

h=(b-a)/N;
x(1)=1;
y(1)=1;
z(1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N
%%%%%%%%%%%%%%%%%%%%%%%%%
k1=h.*feval(f1,x(k),y(k));
k2 = h*f1(x(k)+0.5*h,y(k)+0.5*h*k1,z(k)+0.5*h*k1);
x(k+1) = x(k) + k2;
%%%%%%%%%%%%%%%%%%%%%%%%%
k3=h.*feval(f2,x(k),y(k),z(k));
k4 = h*f2(x(k)+0.5*h,y(k)+0.5*h*k3,z(k)+0.5*h*k3);
y(k+1) = y(k) + k4;
%%%%%%%%%%%%%%%%%%%%%%%%%
k5=h.*feval(f3,x(k),y(k),z(k));
k6 = h*f3(x(k)+0.5*h,y(k)+0.5*h*k5,z(k)+0.5*h*k5);
z(k+1) = z(k) + k6;
end
figure;
plot3(x, y, z, 'b', 'LineWidth', 3); grid on
view(3)