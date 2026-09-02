%%=========================================================================
%% SCRIPT: Solucion_KdV.m
%% SOLUCIÓN DE LA ECUACIÓN KORTEWEG-DE VRIES (KdV) EN OCTAVE
%% Ecuación: u_t + 6*u*u_x + u_xxx = 0
%%
%% Métodos implementados:
%%   1. Deducibilidad de Esténciles por Suma/Resta de Serie de Taylor:
%%      - Derivada Temporal u_t: (u^{n+1} - u^n) / dt
%%      - Primera Derivada Espacial u_x: (u_{i+1} - u_{i-1}) / (2*dx)
%%      - Tercera Derivada Espacial u_xxx: (u_{i+2} - 2*u_{i+1} + 2*u_{i-1} - u_{i-2}) / (2*dx^3)
%%   2. Método Theta (theta = 0.5 para Crank-Nicolson de 2do orden)
%%   3. Freezing (Congelamiento del coeficiente no lineal u^n)
%%   4. Puntos Fantasma (Ghost Points) para condiciones de frontera periódicas
%%
%% Gráficas incluidas:
%%   - Gráfica 1: Solución Analítica (Solitón Exacto) vs. Solución Numérica
%%   - Gráfica 2: Tiempo de Ejecución (CPU Time vs Número de Nodos N)
%%   - Gráfica 3: Animación en Tiempo Real del Solitón Propagándose
%%=========================================================================

function Solucion_KdV(N_in, theta_in)
    
    clc;
    close all;
    
    %% --------------------------------------------------------------------
    %% PARÁMETROS DE ENTRADA E INTERACTIVIDAD
    %% --------------------------------------------------------------------
    if nargin < 1 || isempty(N_in)
        N = 200; %% Número de puntos espaciales por defecto
    else
        N = N_in;
    end
    
    if nargin < 2 || isempty(theta_in)
        theta = 0.5; %% Método Crank-Nicolson por defecto
    else
        theta = theta_in;
    end
    
    fprintf('=======================================================\n');
    fprintf('   SOLUCION NUMERICA Y ANALITICA DE LA ECUACION KdV\n');
    fprintf('=======================================================\n');
    fprintf(' Parametros:\n');
    fprintf('   Nodos espaciales (N): %d\n', N);
    fprintf('   Parametro Theta     : %.2f\n', theta);
    fprintf('=======================================================\n\n');
    
    %% --------------------------------------------------------------------
    %% DOMINIO FÍSICO Y PARÁMETROS DE LA ONDA SOLITARIA (SOLITÓN)
    %% --------------------------------------------------------------------
    L = 40;                 %% Longitud del dominio espacial [-L/2, L/2]
    x = linspace(-L/2, L/2, N)'; %% Vector de posiciones espaciales (columna)
    dx = x(2) - x(1);       %% Paso espacial h = dx
    
    c = 2.0;                %% Velocidad de propagación del solitón
    x0 = -5.0;              %% Posición inicial del pico del solitón
    
    t_final = 4.0;          %% Tiempo final de simulación
    Nt = 400;               %% Número de pasos temporales
    dt = t_final / Nt;      %% Paso temporal
    
    %% --------------------------------------------------------------------
    %% SOLUCIÓN ANALÍTICA EN t = 0 (CONDICIÓN INICIAL)
    %% u(x,t) = (c/2) * sech^2( (sqrt(c)/2) * (x - c*t - x0) )
    %% --------------------------------------------------------------------
    u_exact_fn = @(x_val, t_val) (c / 2) * (sech( (sqrt(c) / 2) * (x_val - c * t_val - x0) ) .^ 2);
    
    u_num = u_exact_fn(x, 0); %% Condición inicial
    
    %% --------------------------------------------------------------------
    %% CONSTRUCCIÓN DE OPERADORES MATRICIALES (DIFERENCIAS FINITAS Y TAYLOR)
    %% Con Puntos Fantasma (Ghost Points) para Frontera Periódica
    %% --------------------------------------------------------------------
    
    %% Operador D1 (Primera Derivada Espacial u_x por Taylor Centrada)
    %% u_x ≈ (u_{i+1} - u_{i-1}) / (2*dx)
    D1 = zeros(N, N);
    for i = 1:N
        ip1 = mod(i, N) + 1;         %% i + 1 con periodicidad
        im1 = mod(i - 2, N) + 1;     %% i - 1 con periodicidad
        D1(i, ip1) = 1 / (2 * dx);
        D1(i, im1) = -1 / (2 * dx);
    end
    
    %% Operador D3 (Tercera Derivada Espacial u_xxx por Taylor 4 Puntos)
    %% u_xxx ≈ (u_{i+2} - 2*u_{i+1} + 2*u_{i-1} - u_{i-2}) / (2*dx^3)
    D3 = zeros(N, N);
    for i = 1:N
        ip2 = mod(i + 1, N) + 1;     %% i + 2 con punto fantasma periódico
        ip1 = mod(i, N) + 1;         %% i + 1
        im1 = mod(i - 2, N) + 1;     %% i - 1
        im2 = mod(i - 3, N) + 1;     %% i - 2 con punto fantasma periódico
        
        D3(i, ip2) =  1 / (2 * dx^3);
        D3(i, ip1) = -2 / (2 * dx^3);
        D3(i, im1) =  2 / (2 * dx^3);
        D3(i, im2) = -1 / (2 * dx^3);
    end
    
    %% --------------------------------------------------------------------
    %% INTEGRACIÓN TEMPORAL (MÉTODO THETA + FREEZING) CON MEDICIÓN DE TIEMPO
    %% --------------------------------------------------------------------
    fprintf('Ejecutando integracion temporal...\n');
    tic; %% Iniciar cronómetro de CPU
    
    u_history = zeros(N, Nt + 1);
    u_history(:, 1) = u_num;
    
    for n = 1:Nt
        %% Freezing: Se congela u^n como coeficiente en el término no lineal 6*u^n*u_x
        V_n = diag(u_num);
        
        %% Matriz implícita A para u^{n+1}
        %% (I + dt*theta*( 6*diag(u^n)*D1 + D3 )) * u^{n+1} = b
        A = eye(N) + dt * theta * (6 * V_n * D1 + D3);
        
        %% Vector explícito b con el nivel de tiempo actual n
        b = u_num - dt * (1 - theta) * (6 * u_num .* (D1 * u_num) + D3 * u_num);
        
        %% Resolver el sistema lineal A * u^{n+1} = b
        u_num = A \ b;
        
        u_history(:, n + 1) = u_num;
    end
    
    tiempo_ejecucion_principal = toc; %% Guardar tiempo de CPU
    fprintf('Simulacion completada en %.4f segundos.\n\n', tiempo_ejecucion_principal);
    
    %% Solución analítica al tiempo final t_final
    u_exact_final = u_exact_fn(x, t_final);
    
    %% Error relativo norma 2
    err_l2 = norm(u_num - u_exact_final) / norm(u_exact_final);
    fprintf('Error relativo L2 al final (t = %.2f): %.6e\n\n', t_final, err_l2);
    
    %% --------------------------------------------------------------------
    %% ANÁLISIS DE TIEMPO DE EJECUCIÓN (CPU TIME VS N)
    %% --------------------------------------------------------------------
    fprintf('Calculando rendimiento (Tiempo de Ejecucion vs N)...\n');
    N_vector = [50, 100, 150, 200, 250, 300];
    tiempos_cpu = zeros(size(N_vector));
    
    for k = 1:length(N_vector)
        Nk = N_vector(k);
        xk = linspace(-L/2, L/2, Nk)';
        dxk = xk(2) - xk(1);
        uk = u_exact_fn(xk, 0);
        
        D1k = zeros(Nk, Nk);
        D3k = zeros(Nk, Nk);
        for i = 1:Nk
            ip2 = mod(i + 1, Nk) + 1; ip1 = mod(i, Nk) + 1;
            im1 = mod(i - 2, Nk) + 1; im2 = mod(i - 3, Nk) + 1;
            D1k(i, ip1) = 1/(2*dxk);   D1k(i, im1) = -1/(2*dxk);
            D3k(i, ip2) = 1/(2*dxk^3); D3k(i, ip1) = -2/(2*dxk^3);
            D3k(i, im1) = 2/(2*dxk^3); D3k(i, im2) = -1/(2*dxk^3);
        end
        
        tic;
        for n = 1:Nt
            Vk = diag(uk);
            Ak = eye(Nk) + dt * theta * (6 * Vk * D1k + D3k);
            bk = uk - dt * (1 - theta) * (6 * uk .* (D1k * uk) + D3k * uk);
            uk = Ak \ bk;
        end
        tiempos_cpu(k) = toc;
    end
    
    %% ====================================================================
    %% GENERACIÓN DE GRÁFICAS REQUERIDAS
    %% ====================================================================
    
    %% --------------------------------------------------------------------
    %% GRÁFICA 1: SOLUCIÓN ANALÍTICA VS NUMÉRICA (t = t_final)
    %% --------------------------------------------------------------------
    figure(1, 'Name', 'Solucion Analitica vs Numerica KdV', 'NumberTitle', 'off');
    plot(x, u_exact_final, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Solución Analítica (Solitón Exacto)');
    hold on;
    plot(x, u_num, 'b--', 'LineWidth', 1.8, 'DisplayName', sprintf('Solución Numérica (N=%d, \\theta=%.2f)', N, theta));
    grid on;
    xlabel('Posición x', 'FontSize', 12);
    ylabel('Amplitud u(x,t)', 'FontSize', 12);
    title(sprintf('Ecuación KdV: Solución Analítica vs. Numérica a t = %.2f s', t_final), 'FontSize', 13);
    legend('Location', 'northeast', 'FontSize', 11);
    hold off;
    drawnow;
    
    %% --------------------------------------------------------------------
    %% GRÁFICA 2: TIEMPO DE EJECUCIÓN (CPU TIME VS N)
    %% --------------------------------------------------------------------
    figure(2, 'Name', 'Tiempo de Ejecucion vs N', 'NumberTitle', 'off');
    plot(N_vector, tiempos_cpu, 'k-o', 'LineWidth', 2.0, 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    grid on;
    xlabel('Número de Nodos Espaciales (N)', 'FontSize', 12);
    ylabel('Tiempo de Ejecución en CPU (segundos)', 'FontSize', 12);
    title('Desempeño Computacional: Tiempo de Ejecución vs. N', 'FontSize', 13);
    drawnow;
    
    %% --------------------------------------------------------------------
    %% GRÁFICA 3: ANIMACIÓN DE LA PROPAGACIÓN DEL SOLITÓN
    %% --------------------------------------------------------------------
    figure(3, 'Name', 'Animacion del Soliton KdV', 'NumberTitle', 'off');
    
    t_vec = linspace(0, t_final, Nt + 1);
    
    for n = 1:5:Nt+1
        t_curr = t_vec(n);
        u_exact_curr = u_exact_fn(x, t_curr);
        u_num_curr = u_history(:, n);
        
        plot(x, u_exact_curr, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Solución Analítica');
        hold on;
        plot(x, u_num_curr, 'b--', 'LineWidth', 1.8, 'DisplayName', 'Solución Numérica');
        grid on;
        xlim([-L/2, L/2]);
        ylim([-0.2, c/2 + 0.3]);
        xlabel('Posición x', 'FontSize', 12);
        ylabel('Amplitud u(x,t)', 'FontSize', 12);
        title(sprintf('Animación de Onda Solitaria KdV (Solitón) - Tiempo t = %.2f s', t_curr), 'FontSize', 13);
        legend('Location', 'northeast');
        hold off;
        
        drawnow;
        pause(0.02); %% Pausa para fluidez de la animación
    end
    
    %% --------------------------------------------------------------------
    %% GRÁFICA 4: EVOLUCIÓN DEL ERROR EN EL TIEMPO (norma L2 vs t)
    %% --------------------------------------------------------------------
    t_vec = linspace(0, t_final, Nt + 1);
    err_l2_t = zeros(1, Nt + 1);
    for n = 1:Nt+1
        u_exact_n = u_exact_fn(x, t_vec(n));
        err_l2_t(n) = norm(u_history(:, n) - u_exact_n) / norm(u_exact_n);
    end

    figure(4, 'Name', 'Evolucion del Error KdV', 'NumberTitle', 'off');
    plot(t_vec, err_l2_t, 'r-', 'LineWidth', 2.0);
    grid on;
    xlabel('Tiempo t (s)', 'FontSize', 12);
    ylabel('Error Relativo L2', 'FontSize', 12);
    title('Evolución del Error Relativo L2 en el Tiempo', 'FontSize', 13);
    drawnow;

    fprintf('  Error relativo L2 final (N=%d): %.6e\n', N, err_l2);
    fprintf('=======================================================\n');
    fprintf(' Proceso completado exitosamente.\n');
    fprintf('=======================================================\n');
    fprintf('\nPresione cualquier tecla para cerrar...\n');
    pause;
    
end
