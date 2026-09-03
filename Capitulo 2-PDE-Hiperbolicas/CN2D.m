clear; clc; close all;

% --- 1. PARÁMETROS DEL DOMINIO Y MALLA ---
L = 8;          % Longitud del dominio cuadrado
Nx = 60;         % Número de puntos en X
Ny = 60;         % Número de puntos en Y
h = L / (Nx-1);  % Paso espacial (dx = dy = h)
x = linspace(0, L, Nx);
y = linspace(0, L, Ny);

[X, Y] = meshgrid(x, y);

% --- 2. PARÁMETROS TEMPORALES Y ESTABILIDAD (CFL) ---
cfl = 0.5;                % Factor de seguridad CFL (< 1/sqrt(2) approx 0.707)
dt = cfl * h / sqrt(2);   % Paso del tiempo estable
t_final = 2;             % Tiempo total de simulación
Nt = ceil(t_final/dt);    % Número de pasos de tiempo

r2 = (dt / h)^2;          % Constante del esquema numérico

% --- 3. ASIGNACIÓN DE MEMORIA ---
u0 = zeros(Nx, Ny);       % Tiempo pasado (n-1)
u1 = zeros(Nx, Ny);       % Tiempo actual (n)
u2 = zeros(Nx, Ny);       % Tiempo futuro (n+1)

% --- 4. CONDICIÓN INICIAL (Perturbación Gaussiana) ---
cx = L/4; cy = L/4;       % Centro de la perturbación
sigma = 0.8;              % Ancho del pulso
u1=exp(-((X-0.2).^2+(Y-cy).^2)/ (2*sigma^2));

% --- 5. CONFIGURACIÓN DE LA FIGURA ---
figure('Position', [100, 100, 800, 600]);
h_surf = surf(X, Y, u2);
shading interp;           % Suavizar colores
colormap jet;
colorbar;
zlim([-0.5, 1]);          % Límites del eje Z
title('Simulación de Onda 2D');
xlabel('X'); ylabel('Y'); zlabel('Amplitud (u)');
view(3);                  % Vista 3D

% --- 6. BUCLE TEMPORAL (Evolución de la onda) ---
for n = 1:Nt
    % Esquema de diferencias finitas para puntos interiores
    for i = 2:Nx-1
        for j = 2:Ny-1
            % Fórmula explícita centrada
            u2(i,j) = 2*u1(i,j)-u0(i,j)+r2*(u2(i+1,j)-2*u2(i,j)+u2(i-1,j)+...
                u0(i-1,j)-2*u0(i,j)+u0(i-1,j)+u2(i,j+1)-2*u2(i,j)+u2(i,j-1)+...
                u0(i,j+1)-2*u0(i,j)+u0(i,j-1));
        end
    end

 
    % Condiciones de Frontera Fijas (Dirichlet: u = 0 en los bordes)
    u2(1, :) = 0; u2(end, :) = 0;
    u2(:, 1) = 0; u2(:, end) = 0;
    
    % Actualizar matrices para el siguiente paso
    u0 = u1;
    u1 = u2;
    
   

    set(h_surf, 'ZData', u1);
    drawnow;
end
