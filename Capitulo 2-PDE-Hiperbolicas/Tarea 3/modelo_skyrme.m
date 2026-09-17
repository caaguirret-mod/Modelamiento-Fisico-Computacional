function modelo_skyrme ()

    clear; clc; close all;
    
    %=======================================%
    % paso 1:definicion malla, y parametros %
    %=======================================%

    L = pi;
    nx = 400;
    x = linspace(0, L, nx + 1);
    dx = x(2) - x(1);

    r = 0.9; % Coeficiente CFL 

    t = 8;
    dt = dx * r;
    nt = round(t / dt); 
    
    % CORREGIDO: Se agregaron las variables en el orden correcto y sin el '\n' extra
    fprintf('dx=%f, nx=%d, dt=%f, nt=%d, r=%f\n', dx, nx, dt, nt, r);

    %=======================================%
    %      paso 2:condiciones iniciales     %
    %=======================================%

    % Ecuación a modelar para t=0:
    % \theta(x,0) = (pi*x)/L + suma [ A_n * sin((n*pi*x)/L) ]
    
    n = 1:3;
    suma = zeros(size(x));

    for i = n
        suma = suma + (0.2 / (i * i)) * sin((i * pi * x) / L);
    end

    theta_0 = (pi * x) / L + suma;
    %=======================================%
    %      paso 3: vectores de estado       %
    %=======================================%

    theta_pas = theta_0;              % estado en n-1
    theta_presente = zeros(size(x));  % estado en n
    theta_futuro = zeros(size(x));    % estado en n+1

    for i = 2:nx
        theta_presente(i) = theta_0(i) + (r^2 / 2) * (theta_0(i+1) - 2 * theta_0(i) + theta_0(i-1));
    end

    theta_presente(1) = 0;
    theta_presente(end) = pi;

    %=======================================%
    %        paso 4: Bucle leapfrog         %
    %=======================================%

    idx_centro = round(nx / 2) + 1;
    time_vec = zeros(1, nt + 1);
    theta_centro = zeros(1, nt + 1);
    theta_ana_centro = zeros(1, nt + 1);
    error_global = zeros(1, nt + 1);

    theta_centro(1) = theta_0(idx_centro);
    theta_ana_centro(1) = theta_0(idx_centro);
    error_global(1) = 0;

    tic;
    for k = 1:nt-1
        for i = 2:nx
            theta_futuro(i) = 2 * theta_presente(i) - theta_pas(i) + ...
                r^2 * (theta_presente(i+1) - 2 * theta_presente(i) + theta_presente(i-1));
        end
        theta_futuro(1) = 0;
        theta_futuro(end) = pi;

        % rotar para siguiente paso
        theta_pas = theta_presente;
        theta_presente = theta_futuro;
        theta_futuro = zeros(size(x));

        % medidas diagnosticas
        t_actual = k * dt;
        time_vec(k + 1) = t_actual;

        theta_ana_k = pi * x / L ...
            + 0.2 * cos(pi * t_actual / L) .* sin(pi * x / L) ...
            + 0.05 * cos(2 * pi * t_actual / L) .* sin(2 * pi * x / L) ...
            + 0.022 * cos(3 * pi * t_actual / L) .* sin(3 * pi * x / L);

        theta_centro(k + 1) = theta_presente(idx_centro);
        theta_ana_centro(k + 1) = theta_ana_k(idx_centro);
        error_global(k + 1) = max(abs(theta_presente - theta_ana_k));
    end
    tiempo_total = toc;
    fprintf('Tiempo total de integracion = %f segundos\n', tiempo_total);

    % Paso 5: solucion analitica n=3 en t_final para validar
    t_final = nt * dt;
    theta_ana = pi * x / L ...
        + 0.2 * cos(pi * t_final / L) .* sin(pi * x / L) ...
        + 0.05 * cos(2 * pi * t_final / L) .* sin(2 * pi * x / L) ...
        + 0.022 * cos(3 * pi * t_final / L) .* sin(3 * pi * x / L);

    % Graficas
    figure;
    subplot(2,1,1);
    plot(x, theta_0, 'linewidth', 2);
    title('theta en t=0');
    xlabel('x'); ylabel('theta'); grid on;

    subplot(2,1,2);
    plot(x, theta_presente, 'b', 'linewidth', 2); hold on;
    plot(x, theta_ana, 'r--', 'linewidth', 1.5);
    title('numerico vs analitico en t final');
    xlabel('x'); ylabel('theta'); grid on;
    legend('leapfrog','analitica');

    % Error global entre solucion numerica y analitica
    figure;
    plot(time_vec(1:end-1), error_global(1:end-1), 'linewidth', 2);
    title('Error global entre soluciones');
    xlabel('Tiempo t');
    ylabel('|theta_numerica - theta_analitica|_{\infty}');
    grid on;

    % Evolucion temporal en un punto fijo del dominio
    figure;
    plot(time_vec(1:end-1), theta_centro(1:end-1), 'b', 'linewidth', 2); hold on;
    plot(time_vec(1:end-1), theta_ana_centro(1:end-1), 'r--', 'linewidth', 1.5);
    title('Evolucion temporal en x = L/2');
    xlabel('Tiempo t');
    ylabel('\theta(x=L/2,t)');
    legend('numerico', 'analitica');
    grid on;

    % Tiempo de computo vs numero de iteraciones
    figure;
    iter_vals = unique(round(logspace(log10(50), log10(nt), 8)));
    time_iter = zeros(size(iter_vals));

    for j = 1:length(iter_vals)
        n_iter = iter_vals(j);
        theta_pas_j = theta_0;
        theta_pres_j = zeros(size(x));
        theta_fut_j = zeros(size(x));

        for i = 2:nx
            theta_pres_j(i) = theta_0(i) + (r^2 / 2) * (theta_0(i+1) - 2 * theta_0(i) + theta_0(i-1));
        end
        theta_pres_j(1) = 0;
        theta_pres_j(end) = pi;

        tic;
        for k = 1:n_iter-1
            for i = 2:nx
                theta_fut_j(i) = 2 * theta_pres_j(i) - theta_pas_j(i) + ...
                    r^2 * (theta_pres_j(i+1) - 2 * theta_pres_j(i) + theta_pres_j(i-1));
            end
            theta_fut_j(1) = 0;
            theta_fut_j(end) = pi;
            theta_pas_j = theta_pres_j;
            theta_pres_j = theta_fut_j;
            theta_fut_j = zeros(size(x));
        end
        time_iter(j) = toc;
    end

    plot(iter_vals, time_iter, 'o-b', 'linewidth', 2);
    title('Tiempo de computo vs numero de iteraciones');
    xlabel('Numero de iteraciones');
    ylabel('Tiempo de computo (s)');
    grid on;

end

