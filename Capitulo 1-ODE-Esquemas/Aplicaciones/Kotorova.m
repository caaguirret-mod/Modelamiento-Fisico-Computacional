close all; clear all; clc

x(1)=0.3;
y(1)=0.1;
imax=100;
K=10;

X = linspace(0, pi, imax);
V=K./(4*pi.^2)*cos(X);
Vprime=-K./(2*pi)*sin(X);
f=@(x) -K./(2*pi)*sin(x);

for k=1:imax
y(k+1)=y(k)+feval(f,x(k));
x(k+1)=x(k)+y(k)+feval(f,x(k));
end 


%[X, Y] = meshgrid(x,y);
%Z = X + Y;  % Función z = f(x, y)
%contour(X, Y, Z);
%shading interp;              % Suaviza los colores de la superficie
%colorbar;     
%figure;
%X=
plot(x, y, 'b', 'LineWidth', 2); grid on
view(2)
