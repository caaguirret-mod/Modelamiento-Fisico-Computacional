%% porous_media_equation.m
% Solucion numerica de la Porous Media Equation en GNU Octave / MATLAB

1; 

%% =================================================================
%  DEFINICIÓN DE FUNCIONES
%% =================================================================

function [U, tiempo_calculo] = resolver_pme(Nx, Nt, dx, dt, r, x, t)
    U = zeros(Nt + 1, Nx + 1);

    % Condicion inicial
    U(1, :) = x;

    tic;

    for n = 1:Nt
        % Condiciones de frontera en t_n
        U(n, 1) = t(n);
        U(n, Nx + 1) = 1 + t(n);

        % Puntos interiores
        for j = 2:Nx
            term_ux = (r / 4) * (U(n, j + 1) - U(n, j - 1))^2;
            term_uxx = r * U(n, j) * (U(n, j + 1) - 2 * U(n, j) + U(n, j - 1));

            U(n + 1, j) = U(n, j) + term_ux + term_uxx;
        end

        % Condiciones de frontera en t_(n+1)
        U(n + 1, 1) = t(n + 1);
        U(n + 1, Nx + 1) = 1 + t(n + 1);
    end

    tiempo_calculo = toc;
end

function U_exacta = solucion_exacta(x, t)
    [X, Tm] = meshgrid(x, t);
    U_exacta = X + Tm;
end


%% ================================================================
%  PROGRAMA PRINCIPAL
%% ================================================================

clear;
clc;
close all;

% Parametros del problema
x_min = 0.0;
x_max = 1.0;
T = 10.0;
Nx = 10;
dx = (x_max - x_min)/Nx;
dt = 0.0001;
Nt = round(T/dt);
r = dt/(dx^2);

x = linspace(x_min, x_max, Nx + 1);
t = linspace(0, T, Nt + 1);

% Informacion del calculo
fprintf('=============================================\n');
fprintf(' POROUS MEDIA EQUATION\n');
fprintf('=============================================\n');
fprintf('dx = %.6f\n', dx);
fprintf('dt = %.6f\n', dt);
fprintf('r  = %.6f\n', r);
fprintf('Nx = %d\n', Nx);
fprintf('Nt = %d\n', Nt);
fprintf('Tiempo final = %.2f\n', T);
fprintf('=============================================\n\n');

% Solucion numerica
tic;
[U, tiempo_calculo] = resolver_pme(Nx, Nt, dx, dt, r, x, t);
toc;

fprintf('Tiempo de procesamiento: %.6f segundos\n\n', tiempo_calculo);

% Solucion exacta
U_exacta = solucion_exacta(x, t);

% Error absoluto
Error = abs(U - U_exacta);
error_maximo = max(Error(:));
error_rms = sqrt(mean(Error(:).^2));

fprintf('Error maximo global = %.12e\n', error_maximo);
fprintf('Error RMS global    = %.12e\n\n', error_rms);

% GRAFICA 1: SUPERFICIE 3D
[X, Tm] = meshgrid(x, t);
figure;
surf(X, Tm, U);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
title('Propagacion temporal de la solucion numerica');
colorbar; shading interp; grid on; view(45, 30);

% GRAFICA 2: NUMERICA VS EXACTA (TIEMPO FINAL)
u_numerica_final = U(end, :);
u_exacta_final = U_exacta(end, :);

figure;
plot(x, u_numerica_final, 'o-', 'LineWidth', 1.2);
hold on;
plot(x, u_exacta_final, '--', 'LineWidth', 1.5);
xlabel('x'); ylabel('u(x,T)');
title('Solucion numerica vs. solucion exacta');
legend('Numerica', 'Exacta', 'Location', 'northwest');
grid on; hold off;

% GRAFICA 3: ERROR ABSOLUTO
figure;
surf(X, Tm, Error);
xlabel('x'); ylabel('t'); zlabel('|u_{numerica}-u_{exacta}|');
title('Error absoluto de la solucion numerica');
colorbar; shading interp; grid on; view(45, 30);

% GRAFICA 4: PERFILES PARA VARIOS TIEMPOS
tiempos_grafica = [0, 2, 5, 8, 10];
figure; hold on;
for k = 1:length(tiempos_grafica)
    indice = round(tiempos_grafica(k)/dt) + 1;
    plot(x, U(indice, :), 'LineWidth', 1.3);
end
xlabel('x'); ylabel('u(x,t)');
title('Evolucion espacial de la solucion para varios tiempos');
legend('t = 0', 't = 2', 't = 5', 't = 8', 't = 10', 'Location', 'northwest');
grid on; hold off;

% GRAFICA 5: TIEMPO VS ITERACIONES
iteraciones = [1000, 5000, 10000, 25000, 50000, 100000];
tiempos_ejecucion = zeros(size(iteraciones));

for k = 1:length(iteraciones)
    Nt_prueba = iteraciones(k);
    T_prueba = Nt_prueba * dt;
    t_prueba = linspace(0, T_prueba, Nt_prueba + 1);

    tic;
    resolver_pme(Nx, Nt_prueba, dx, dt, r, x, t_prueba);
    tiempos_ejecucion(k) = toc;
end

figure;
plot(iteraciones, tiempos_ejecucion, 'o-', 'LineWidth', 1.3);
xlabel('Numero de iteraciones'); ylabel('Tiempo de procesamiento (s)');
title('Tiempo real de procesamiento vs. numero de iteraciones');
grid on;

% TABLA DE RESULTADOS
fprintf('=============================================\n');
fprintf(' TIEMPO VS. ITERACIONES\n');
fprintf('=============================================\n');
fprintf('Iteraciones        Tiempo (s)\n');
fprintf('---------------------------------------------\n');
for k = 1:length(iteraciones)
    fprintf('%10d        %.8f\n', iteraciones(k), tiempos_ejecucion(k));
end
fprintf('=============================================\n');
