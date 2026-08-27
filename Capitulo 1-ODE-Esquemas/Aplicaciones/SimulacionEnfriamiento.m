function [T_rk4, T_am] = SimulacionEnfriamiento(N)
% Ecuacion: T' = -k*(T - (Tm + C*sin(w*t)))
% Condicion inicial: T(0) = 100
% Metodos: Runge-Kutta 4 y Adams-Moulton (Predictor-Corrector)

% Parametros del modelo
T0 = 100.0;
Tm = 20.0;
C = 10.0;
w = pi / 12;
k_enf = 0.5;

a = 0.0;
b = 48.0;
h = (b-a)/N;

% Funciones
T_env = @(t) Tm + C .* sin(w .* t);
f = @(t, T) -k_enf .* (T - T_env(t));

% Inicializacion
t = zeros(1, N+1);
T_rk4 = zeros(1, N+1);
T_am = zeros(1, N+1);
tiempos_rk4 = zeros(1, N+1);
tiempos_am = zeros(1, N+1);

t(1) = a;
T_rk4(1) = T0;
T_am(1) = T0;

% RK4
tic
for k = 1:N
    tic;
    k1 = f(t(k), T_rk4(k));
    k2 = f(t(k) + 0.5*h, T_rk4(k) + 0.5*h*k1);
    k3 = f(t(k) + 0.5*h, T_rk4(k) + 0.5*h*k2);
    k4 = f(t(k) + h, T_rk4(k) + h*k3);

    T_rk4(k+1) = T_rk4(k) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    t(k+1) = t(k) + h;

    tiempos_rk4(k+1) = toc * 1e6;
end
toc
% Adams-Moulton
for k = 1:min(4, N+1)
    T_am(k) = T_rk4(k);
end
tic
for k = 4:N
    tic;
    T_pred = T_am(k) + (h/24) * (55*f(t(k), T_am(k)) - 59*f(t(k-1), T_am(k-1)) + ...
                                 37*f(t(k-2), T_am(k-2)) - 9*f(t(k-3), T_am(k-3)));

    T_am(k+1) = T_am(k) + (h/24) * (9*f(t(k+1), T_pred) + 19*f(t(k), T_am(k)) - ...
                                    5*f(t(k-1), T_am(k-1)) + f(t(k-2), T_am(k-2)));
    tiempos_am(k+1) = toc * 1e6;
end
toc
% Solucion analitica
T_exacta = Tm + (k_enf*C/(k_enf^2 + w^2))*(k_enf.*sin(w.*t) - w.*cos(w.*t)) + ...
           (T0 - Tm + (k_enf*C*w)/(k_enf^2 + w^2)) .* exp(-k_enf.*t);

% Graficas
figure('Position', [100, 100, 1000, 700]);

subplot(2,2,1);
plot(t, T_exacta, 'c-', 'LineWidth', 2); hold on;
plot(t, T_rk4, 'b--', 'LineWidth', 1);
title('Solucion Analitica vs RK4');
xlabel('Tiempo (h)');
ylabel('Temperatura (°C)');
legend('Analitica', 'RK4');
grid on; hold off;

subplot(2,2,2);
plot(t, T_exacta, 'c-', 'LineWidth', 2); hold on;
plot(t, T_am, 'r--', 'LineWidth', 1);
title('Solucion Analitica vs Adams-Moulton');
xlabel('Tiempo (h)');
ylabel('Temperatura (°C)');
legend('Analitica', 'Adams-Moulton');
grid on; hold off;

subplot(2,2,[3 4]);
plot(t(2:end), tiempos_rk4(2:end), 'b-', 'LineWidth', 1.5); hold on;
plot(t(5:end), tiempos_am(5:end), 'r-', 'LineWidth', 1.5);
title('Costo Computacional por Iteracion');
xlabel('Tiempo de simulacion (h)');
ylabel('Tiempo (\mus)');
legend('RK4', 'Adams-Moulton');
grid on; hold off;

end
