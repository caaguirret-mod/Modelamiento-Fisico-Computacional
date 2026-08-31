clear; clc; close all;
Nx = 60;         % Número de puntos en X
Ny = 60;  
L=10; % Número de puntos en Y
h = L / (Nx-1);  % Paso espacial (dx = dy = h)
x = linspace(0, L, Nx);
y = linspace(0, L, Ny);


[X, Y] = meshgrid(x, y);

% --- 2. PARÁMETROS TEMPORALES Y ESTABILIDAD (CFL) ---
cfl = 0.5;                % Factor de seguridad CFL (< 1/sqrt(2) approx 0.707)
dt = cfl * h / sqrt(2);   % Paso del tiempo estable
t_final = 8;             % Tiempo total de simulación
Nt = ceil(t_final/dt);    % Número de pasos de tiempo

r2 = (dt / h)^2; 
cx = L; cy = L/2;       % Centro de la perturbación
sigma = 0.8;% Constante del esquema numérico
f=@(X,Y) exp(-((X-cx).^2+(Y-cy).^2)/ (2*sigma^2));

u=zeros(Nx,Ny,Nt);

for k=2:Ny-1
  for j=2:Nx-1
      u(j,k,1)=2*(1 - 4*r2)*feval(f,j,k)-feval(f,j,k)+ ...
                      r2*(feval(f,j+1,k) + feval(f,j-1,k) + feval(f,j,k+1) + feval(f,j,k-1));
  end
end 
 u(1:Nx,1,1)=0;
 u(1:Nx,Ny,1)=0;
 u(1,1:Ny,1)=0;
 u(Nx,1:Ny,1)=0;

 for n=2:Nt
 for k=2:Ny-1
  for j=2:Nx-1
      u(j,k,n)=2*(1 - 4*r2)*u(j,k,n-1) - u(j,k,n-1) + ...
                      r2*(u(j+1,k,n-1) + u(j-1,k,n-1) + u(j,k+1,n-1)+u(j,k-1,n-1)); 
  end
 end 
 u(1:Nx,1,n)=0;
 u(1:Nx,Ny,n)=0;
 u(1,1:Ny,n)=0;
 u(Nx,1:Ny,n)=0;

                % Vista 3D

 end

 subplot(2,2,1);
 h_surf = surf(X, Y, u(:,:,2));
shading interp;           % Suavizar colores
colormap jet;
colorbar;
zlim([-0.5, 1]);          % Límites del eje Z
title('Simulación de Onda 2D');
xlabel('X'); ylabel('Y'); zlabel('Amplitud (u)');
view(3);  

 subplot(2,2,2);
 h_surf = surf(X, Y, u(:,:,Nx/3));
shading interp;           % Suavizar colores
colormap jet;
colorbar;
zlim([-0.5, 1]);          % Límites del eje Z
title('Simulación de Onda 2D');
xlabel('X'); ylabel('Y'); zlabel('Amplitud (u)');
view(3);  

 subplot(2,2,3);
 h_surf = surf(X, Y, u(:,:,2*Nx/3));
shading interp;           % Suavizar colores
colormap jet;
colorbar;
zlim([-0.5, 1]);          % Límites del eje Z
title('Simulación de Onda 2D');
xlabel('X'); ylabel('Y'); zlabel('Amplitud (u)');
view(3);  

 subplot(2,2,3);
 h_surf = surf(X, Y, u(:,:,Nx));
shading interp;           % Suavizar colores
colormap jet;
colorbar;
zlim([-0.5, 1]);          % Límites del eje Z
title('Simulación de Onda 2D');
xlabel('X'); ylabel('Y'); zlabel('Amplitud (u)');
view(3);  