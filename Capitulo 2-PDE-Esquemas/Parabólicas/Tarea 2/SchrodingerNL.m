%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u, error_L2_final] = SchrodingerNL(J)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esquema IMEX (Crank-Nicolson + Freezing con Matrices Dispersas)
% Ecuacion i*u_t + alpha*u_xx + kappa*|u|^2*u = 0, x=[-20,20]
% Incluye las 7 visualizaciones completas de análisis físico y computacional.
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Parámetros de la malla espacial
x1 = -20;
x2 = 20;
dx = (x2 - x1) / J;
N_nodes = J + 1;
x = (x1 : dx : x2)';

% 2. Parámetros temporales
tf = 10.0;
Nt = 1000;
dt = tf / (Nt - 1);
tt = linspace(0, tf, Nt);

% 3. Parámetros Físicos de la Ecuación NLSE
alpha = 0.5;
kappa = 1.0;

% 4. Inicialización de matrices
u = zeros(N_nodes, Nt);
u_exacta = zeros(N_nodes, Nt);

% Condición Inicial: Solitón viajero en t=0
A = 1.0;
v = 2.0;
x0 = -10.0;
u(:, 1) = A * sech(A * (x - x0)) .* exp(1i * v * x);
u_exacta(:, 1) = u(:, 1);

% 5. Constantes del esquema Crank-Nicolson
r1 = (1i * alpha * dt) / (2 * dx^2);
r2 = (1i * kappa * dt) / 2;

e = ones(N_nodes, 1);
off_diag_A = -r1 * e;
off_diag_B = r1 * e;

% Vector para medir el tiempo de CPU por iteración
tiempos_por_iter = zeros(Nt-1, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUCLE TEMPORAL PRINCIPAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Ejecutando simulación numérica y calculando solución analítica...\n');

for n = 2:Nt
    t_inicio_iter = tic;

    % --- Paso Numérico IMEX ---
    u_old = u(:, n-1);
    V_n = abs(u_old).^2;

    diag_A = 1 + 2*r1 - r2 * V_n;
    diag_B = 1 - 2*r1 + r2 * V_n;

    A_mat = spdiags([off_diag_A, diag_A, off_diag_A], [-1, 0, 1], N_nodes, N_nodes);
    B_mat = spdiags([off_diag_B, diag_B, off_diag_B], [-1, 0, 1], N_nodes, N_nodes);

    A_mat(1, N_nodes) = -r1;  A_mat(N_nodes, 1) = -r1;
    B_mat(1, N_nodes) = r1;   B_mat(N_nodes, 1) = r1;

    b_vec = B_mat * u_old;
    u(:, n) = A_mat \ b_vec;

    tiempos_por_iter(n-1) = toc(t_inicio_iter);

    % --- Cálculo Analítico Exacto ---
    fase_exacta = v * x - ((v^2 - A^2) / 2) * tt(n);
    u_exacta(:, n) = A * sech(A * (x - x0 - v * tt(n))) .* exp(1i * fase_exacta);
end

tiempo_acumulado = cumsum(tiempos_por_iter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANÁLISIS DE ERROR GLOBAL Y ESPACIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_L2_tiempo = zeros(1, Nt);
for n = 1:Nt
    error_L2_tiempo(n) = norm(u(:, n) - u_exacta(:, n)) / norm(u_exacta(:, n));
end

error_L2_final = error_L2_tiempo(end);
intensidad_num_final = abs(u(:, end)).^2;
intensidad_exacta_final = abs(u_exacta(:, end)).^2;
error_absoluto_espacial = abs(intensidad_num_final - intensidad_exacta_final);

fprintf('\n--- RESULTADOS FINALES ---\n');
fprintf('Error Relativo L2 al tiempo final (t = %.1f s): %.6e\n', tf, error_L2_final);
fprintf('Tiempo total de CPU: %.3f segundos.\n\n', sum(tiempos_por_iter));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZACIÓN 1: SUPERFICIE 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
colormap(jet);
surf(x, tt, abs(u').^2, 'EdgeColor', 'none');
xlabel('Espacio (x)'); ylabel('Tiempo (t)'); zlabel('Intensidad |u|^2');
title('Evolución Espaciotemporal del Solitón (NLSE)');
view([30 45]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZACIÓN 2: PERFILES TRANSVERSALES ESPACIALES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
idx_1 = 1; idx_2 = round(Nt / 3); idx_3 = round(2 * Nt / 3); idx_4 = Nt;

plot(x, abs(u(:, idx_1)).^2, 'k--', 'LineWidth', 1.8, 'DisplayName', sprintf('t = %.2f s', tt(idx_1))); hold on;
plot(x, abs(u(:, idx_2)).^2, 'b-',  'LineWidth', 1.8, 'DisplayName', sprintf('t = %.2f s', tt(idx_2)));
plot(x, abs(u(:, idx_3)).^2, 'g-',  'LineWidth', 1.8, 'DisplayName', sprintf('t = %.2f s', tt(idx_3)));
plot(x, abs(u(:, idx_4)).^2, 'r-',  'LineWidth', 1.8, 'DisplayName', sprintf('t = %.2f s', tt(idx_4))); hold off;
grid on; xlabel('Espacio (x)'); ylabel('Intensidad |u|^2');
title('Perfiles Transversales del Solitón (Fijos en el Tiempo)');
legend('Location', 'northeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZACIÓN 3: RENDIMIENTO (CPU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
subplot(2, 1, 1);
plot(2:Nt, tiempos_por_iter * 1000, 'b-', 'LineWidth', 1.2);
grid on; ylabel('Tiempo / Iteración (ms)'); title('Rendimiento Computacional');

subplot(2, 1, 2);
plot(2:Nt, tiempo_acumulado, 'r-', 'LineWidth', 1.5);
grid on; xlabel('Paso de Tiempo (n)'); ylabel('Tiempo Acumulado (s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZACIÓN 4: NUMÉRICA VS ANALÍTICA AL FINAL
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
plot(x, intensidad_exacta_final, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Solución Analítica'); hold on;
plot(x, intensidad_num_final, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('Solución Numérica (J=%d)', J)); hold off;
grid on; xlabel('Espacio (x)'); ylabel('Intensidad |u|^2');
title(sprintf('Comparación Espacial a t = %.1f s', tf));
legend('Location', 'northeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZACIÓN 5: ERROR ABSOLUTO ESPACIAL (AL FINAL)
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5);
plot(x, error_absoluto_espacial, 'k-', 'LineWidth', 1.5);
grid on; xlabel('Espacio (x)'); ylabel('Error Absoluto | |u_{num}|^2 - |u_{exacta}|^2 |');
title(sprintf('Error Absoluto Espacial al Final (Error L_2 = %.2e)', error_L2_final));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZACIÓN 6: PERFIL TEMPORAL (OSCILOSCOPIO EN x = 0)
%%%%%%%%%%%%%%%%%%%%%%%%%
[~, idx_sonda] = min(abs(x - 0.0));
intensidad_tiempo_num = abs(u(idx_sonda, :)).^2;
intensidad_tiempo_exa = abs(u_exacta(idx_sonda, :)).^2;

figure(6);
plot(tt, intensidad_tiempo_exa, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Señal Teórica'); hold on;
plot(tt, intensidad_tiempo_num, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Señal Simulada'); hold off;
grid on; xlabel('Tiempo (t)'); ylabel('Intensidad |u|^2');
title(sprintf('Corte Temporal en la posición estática x = %.1f', x(idx_sonda)));
legend('Location', 'northeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZACIÓN 7: EVOLUCIÓN DEL ERROR L2 EN EL TIEMPO
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);
plot(tt, error_L2_tiempo, 'k-', 'LineWidth', 1.8);
grid on; xlabel('Tiempo (t)'); ylabel('Error Relativo Norma L_2');
title('Crecimiento del Error Numérico Global en el Tiempo');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZACIÓN 8: ERROR ABSOLUTO TEMPORAL (EN x = 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_absoluto_temporal = abs(intensidad_tiempo_num - intensidad_tiempo_exa);

figure(8);
plot(tt, error_absoluto_temporal, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Tiempo (t)');
ylabel('Error Absoluto | |u_{num}|^2 - |u_{exacta}|^2 |');
title(sprintf('Error Absoluto Temporal en la sonda (x = %.1f)', x(idx_sonda)));

fprintf('¡Las 7 figuras de análisis se han generado e indexado con éxito!\n');
end
