% =====================================================================
% OsciladorArmonicoSimple.m
%
% Modelamiento físico y computacional de la ECUACIÓN DE UN OSCILADOR
% ARMÓNICO SIMPLE (basado en el PDF de la primera tarea).
%
% La ecuación diferencial de segundo orden del sistema resorte-masa es:
%
%        m·d²x/dt² + κ·x = 0          (resorte ideal sin fricción)
%
% dividiendo entre la masa m y llamando ω² = κ/m (frecuencia angular
% al cuadrado), la ecuación queda:
%
%        d²x/dt² + ω²·x = 0
%
% Su solución analítica (deducida en el PDF con el ansatz x = e^(r·t)
% y la identidad de Euler) es:
%
%        x(t) = A·cos(ωt) + B·sin(ωt)
%
% En este programa esa ecuación se resuelve de forma NUMÉRICA con dos
% métodos vistos en clase, tomando como referencia los archivos:
%   - PolC.m      -> método de interpolación (dos pasos)
%   - TrunTaylor.m-> método de Taylor truncado
%
% Para poder usar esos métodos (diseñados para EDO de primer orden),
% convertimos la EDO de segundo orden en un SISTEMA de dos EDO de
% primer orden definiendo la velocidad v = dx/dt:
%
%        dx/dt = v
%        dv/dt = -ω²·x
%
% GRÁFICAS QUE SE OBTIENEN:
%   1) Solución analítica vs solución numérica  -> x(t) vs t
%   2) Tiempo total de CPU vs N (rendimiento)   -> loglog
%   3) Tiempo de CPU con N del usuario (barras)
%   4) Error máximo vs N (convergencia)          -> loglog
%
% Uso (todos los argumentos son opcionales, se pueden poner de a uno):
%   OsciladorArmonicoSimple()                        -> valores por defecto
%   OsciladorArmonicoSimple(m, k, A, B, N, orden)
%   OsciladorArmonicoSimple(2.5, 10, 1, 0.7, 500, 4)
% =====================================================================

function OsciladorArmonicoSimple(m, k, A, B, N, orden)

  % Limpiamos la ventana de comandos y cerramos las figuras abiertas.
  % OJO: NO usamos "clear" porque borraría los argumentos m,k,A,B,N,
  % orden que acaban de llegar de la terminal.
  clc; close all;

  % ================= PARÁMETROS FÍSICOS DEL SISTEMA =================
  % El problema original (PDF) define los parámetros:
  %   κ = constante elástica del resorte, m = masa.
  % Si el usuario NO pasa un argumento, usamos el valor por defecto.
  % Nota: la variable local de entrada m se usa también para contar
  % cuántos argumentos llegaron (nargin), por eso la reasignamos.
  if nargin < 1, m = 1.0; end    % masa m por defecto en kg
  if nargin < 2, k = 4.0; end    % constante elástica κ por defecto (N/m)
  if nargin < 3, A = 0.5; end    % amplitud A por defecto (m)
  if nargin < 4, B = 0.3; end    % amplitud B por defecto (m)
  if nargin < 5, N = 400; end    % número de iteraciones por defecto
  if nargin < 6, orden = 4; end  % orden del Taylor por defecto (2..5)

  w = sqrt(k/m);        % frecuencia angular ω = sqrt(κ/m) en rad/s

  % Amplitudes A y B de la solución analítica x(t)=A·cos(ωt)+B·sin(ωt).
  % Son datos del problema y a la vez nos sirven para fijar las
  % condiciones iniciales:
  %   x(0) = A·cos(0) + B·sin(0) = A              -> posición inicial
  %   v(0) = -A·ω·sin(0) + B·ω·cos(0) = B·ω       -> velocidad inicial
  x0 = A;               % posición inicial x(0)
  v0 = B*w;             % velocidad inicial v(0) = dx/dt(0)

  % ================= INTERVALO DE TIEMPO ============================
  % Simulamos un periodo completo del movimiento. Como la frecuencia
  % angular es ω, el periodo vale T = 2π/ω.
  a = 0.0;              % tiempo inicial (s)
  b = 2*pi/w;           % tiempo final = un periodo T (s)

  % =================================================================
  % PARTE 1: SOLUCIÓN NUMÉRICA CON LOS DOS MÉTODOS
  % =================================================================

  % Solución numérica 1: método de interpolación de dos pasos
  % (adaptado de PolC.m). Devuelve el vector de tiempos tI y de
  % posiciones xI.
  [tI, xI] = metodoInterpolacion(N, a, b, w, x0, v0);

  % Solución numérica 2: método de Taylor truncado de orden 4
  % (adaptado de TrunTaylor.m). Devuelve tT y xT.
  [tT, xT] = metodoTaylor(N, a, b, w, x0, v0, orden);

  % Solución analítica x(t)=A·cos(ωt)+B·sin(ωt) evaluada en la misma
  % malla de tiempo (tiene N+1 nodos igual que las numéricas), para
  % poder comparar punto a punto.
  tA = a : (b - a)/N : b;                  % malla de tiempo (s)
  xA = A*cos(w*tA) + B*sin(w*tA);          % solución exacta (m)

  % Error absoluto en cada nodo: |solución exacta - solución numérica|.
  % Esto nos dice qué tan lejos está cada método de la solución real.
  eI = abs(xA - xI);      % error del método de interpolación
  eT = abs(xA - xT);      % error del método de Taylor truncado

  % Mostramos en la consola un resumen con los errores máximos.
  fprintf('---- RESUMEN DEL MODELAMIENTO ----\n');
  fprintf('ω = %.4f rad/s,  periodo T = %.4f s,  pasos N = %d,  h = %.6f s\n', ...
          w, b - a, N, (b - a)/N);
  fprintf('Error máx. interpolación : %.4e\n', max(eI));
  fprintf('Error máx. Truncamiento de Taylor orden %d: %.4e\n', orden, max(eT));

  % =================================================================
  % GRÁFICA 1: SOLUCIÓN ANALÍTICA vs SOLUCIÓN NUMÉRICA
  % =================================================================
  % Dos subplots lado a lado con ajuste manual de posiciones.
  % Para reducir la densidad de marcadores (MarkerIndices no existe en
  % Octave 7.1), se submuestrea el vector numérico cada 'skip' puntos.
  figure(1);
  skip = max(1, round(N / 30));   % ~30 marcadores visibles en total

  % ---- Panel izquierdo: interpolación ----
  subplot(1,2,1, 'Position', [0.08 0.14 0.38 0.74]);
  plot(tA, xA, 'k-', 'LineWidth', 1.8); hold on;
  plot(tI(1:skip:end), xI(1:skip:end), 'b-o', ...
       'LineWidth', 1.0, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
  grid on;
  xlabel('tiempo t (s)', 'FontSize', 9);
  ylabel('x(t) (m)', 'FontSize', 9);
  title('Interpolación', 'FontSize', 11, 'FontWeight', 'bold');
  legend('Analítica', 'Numérica', 'Location', 'northeast', 'FontSize', 8);
  hold off;

  % ---- Panel derecho: Taylor truncado ----
  subplot(1,2,2, 'Position', [0.56 0.14 0.38 0.74]);
  plot(tA, xA, 'k-', 'LineWidth', 1.8); hold on;
  plot(tT(1:skip:end), xT(1:skip:end), 'r-s', ...
       'LineWidth', 1.0, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
  grid on;
  xlabel('tiempo t (s)', 'FontSize', 9);
  ylabel('x(t) (m)', 'FontSize', 9);
  title(sprintf('Trunc. Taylor orden %d', orden), ...
        'FontSize', 11, 'FontWeight', 'bold');
  legend('Analítica', 'Numérica', 'Location', 'northeast', 'FontSize', 8);
  hold off;

  % =================================================================
  % BARRIDO DE N: TIEMPOS Y ERRORES
  % =================================================================
  % Para cada valor de N medimos el tiempo total de ejecución y el error
  % máximo contra la solución analítica. Esto permite construir las
  % gráficas de rendimiento (tiempo vs N) y de convergencia (error vs N).
  Ns = [50 100 200 400 800 1600 3200 6400];   % valores de N a probar
  if ~ismember(N, Ns)
    Ns = [Ns, N];
    Ns = sort(Ns);
  end

  cpuI = zeros(size(Ns));       % tiempo total de ejecución — interpolación
  cpuT = zeros(size(Ns));       % tiempo total de ejecución — Taylor
  errI = zeros(size(Ns));       % error máximo — interpolación
  errT = zeros(size(Ns));       % error máximo — Taylor

  repeticiones = 5;             % repeticiones por N para promediar tiempos

  for i = 1:numel(Ns)
    n = Ns(i);                  % N actual

    % Solución analítica en la malla de N pasos (se re calcula en cada N)
    tA_n = a : (b - a)/n : b;
    xA_n = A*cos(w*tA_n) + B*sin(w*tA_n);

    % ---- Método de interpolación: tiempo y error ----
    t_aux = 0;
    for r = 1:repeticiones
      tic;
      [tI_n, xI_n] = metodoInterpolacion(n, a, b, w, x0, v0);
      t_aux = t_aux + toc;
    end
    cpuI(i) = t_aux / repeticiones;       % tiempo total promedio
    errI(i) = max(abs(xA_n - xI_n));     % error máximo punto a punto

    % ---- Método de Taylor truncado: tiempo y error ----
    t_aux = 0;
    for r = 1:repeticiones
      tic;
      [tT_n, xT_n] = metodoTaylor(n, a, b, w, x0, v0, orden);
      t_aux = t_aux + toc;
    end
    cpuT(i) = t_aux / repeticiones;       % tiempo total promedio
    errT(i) = max(abs(xA_n - xT_n));     % error máximo punto a punto
  end

  % =================================================================
  % GRÁFICA 2: TIEMPO TOTAL DE EJECUCIÓN vs NÚMERO DE ITERACIONES
  % =================================================================
  % Eje X: N (número de iteraciones, escala logarítmica)
  % Eje Y: tiempo total de ejecución en segundos (escala logarítmica)
  % Cada punto corresponde al tiempo promedio de 'repeticiones' corridas
  % del método completo (desde la condición inicial hasta el último paso).
  figure(2);
  loglog(Ns, cpuI, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
  loglog(Ns, cpuT, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 7);
  grid on;
  set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', ':');
  xlabel('Número de iteraciones N');
  ylabel('Tiempo total de ejecución (s)');
  title('Rendimiento: tiempo de CPU vs N');
  legend('Interpolación', sprintf('Truncamiento de Taylor orden %d', orden), ...
         'Location','NorthWest');
  hold off;

  % =================================================================
  % GRÁFICA 3: TIEMPO TOTAL CON EL N DEL USUARIO (barras)
  % =================================================================
  % Diagrama de barras que compara el tiempo de ambos métodos usando el
  % N específico que el usuario eligió. El valor se escribe encima de
  % cada barra para que sea fácil leerlo.
  figure(3);
  bar([cpuI(Ns == N), cpuT(Ns == N)], 0.4);
  grid on;
  set(gca, 'XTickLabel', {'Interpolación', ...
          sprintf('Truncamiento de Taylor orden %d', orden)});
  xlabel('Método numérico');
  ylabel('Tiempo de CPU (s)');
  title(sprintf('Tiempo de CPU con N = %d pasos', N));
  text(1, cpuI(Ns == N), sprintf('%.4f s', cpuI(Ns == N)), ...
       'HorizontalAlignment','center', 'VerticalAlignment','bottom');
  text(2, cpuT(Ns == N), sprintf('%.4f s', cpuT(Ns == N)), ...
       'HorizontalAlignment','center', 'VerticalAlignment','bottom');

  % =================================================================
  % GRÁFICA 4: ERROR GLOBAL vs NÚMERO DE ITERACIONES
  % =================================================================
  % Compara la precisión real de ambos métodos. Eje X: N (escala log),
  % Eje Y: error máximo |x_exacta - x_numérica| (escala log).
  % La pendiente de la recta en log-log indica el orden de convergencia.
  figure(4);
  loglog(Ns, errI, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
  loglog(Ns, errT, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 7);
  grid on;
  set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', ':');
  xlabel('Número de iteraciones N');
  ylabel('Error máximo |x_{exacta} - x_{numérica}|');
  title('Convergencia: error vs N');
  legend('Interpolación', sprintf('Truncamiento de Taylor orden %d', orden), ...
         'Location','NorthEast');
  hold off;

  % Mostramos el resultado también en la consola.
  fprintf('Tiempo de CPU por número de iteraciones:\n');
  fprintf('  %-10s | %-18s | %-18s | %-18s | %-18s\n', 'N', ...
          'T(s) interp', sprintf('T(s) Trunc. Taylor %d', orden), ...
          'Err máx interp', sprintf('Err máx Trunc. Taylor %d', orden));
  fprintf('  %-10s-+-%-18s-+-%-18s-+-%-18s-+-%-18s\n', ...
          '----------', '------------------', '------------------', ...
          '------------------', '------------------');
  for i = 1:numel(Ns)
    fprintf('  %-10d | %14.6f s | %14.6f s | %14.4e | %14.4e\n', ...
            Ns(i), cpuI(i), cpuT(i), errI(i), errT(i));
  end

  % =================================================================
  % FIN: guardamos las figuras como PNG y pausamos.
  % =================================================================
  % Guardamos las figuras en PNG para que se puedan usar en la
  % presentación sin necesidad de sacar captura de pantalla.
  print(1, 'grafica_solucion_vs_analitica.png', '-dpng', '-r150');
  print(2, 'grafica_tiempo_vs_N.png', '-dpng', '-r150');
  print(3, 'grafica_barras_N_usuario.png', '-dpng', '-r150');
  print(4, 'grafica_error_vs_N.png', '-dpng', '-r150');
  fprintf('Figuras guardadas como PNG en la carpeta del proyecto.\n');
  fprintf('Presiona ENTER para cerrar las figuras...\n');
  pause;

end   % ============ fin de la función principal ============

% =====================================================================
% SUBFUNCIÓN: metodoInterpolacion
%
% Resuelve el sistema  dx/dt = v ; dv/dt = -ω²·x  con el método de
% INTERPOLACIÓN de dos pasos (el mismo esquema del archivo PolC.m).
%
% La fórmula del método (para una EDO y'=f(x,y)) es:
%   y(k+1) = y(k) + (h/2)·[ f(x(k-1),y(k-1)) + f(x(k),y(k)) ]
% que es un trapecio explícito de orden 2.
%
% Como el primer paso necesita los valores en dos nodos anteriores,
% arrancamos con un paso de Heun (Runge-Kutta de orden 2), igual que
% hace PolC.m.
% =====================================================================
function [t, x] = metodoInterpolacion(N, a, b, w, x0, v0)

  h = (b - a) / N;                  % paso de tiempo Δt = (b-a)/N

  % Preasignamos los vectores con (N+1) posiciones para que Octave no
  % los vaya agrandando en cada vuelta del bucle (es más rápido).
  t = zeros(1, N + 1);              % vector de tiempos
  x = zeros(1, N + 1);              % vector de posiciones
  v = zeros(1, N + 1);              % vector de velocidades (auxiliar)

  % Condiciones iniciales en el primer nodo.
  t(1) = a;                         % t_0 = a
  x(1) = x0;                        % x(0)
  v(1) = v0;                        % v(0)

  % ---- Arranque: un paso de Heun (RK2) para obtener los nodos 2 ----
  % k1 = f(t1, estado1)  y  k2 = f(t1+h, estado1 + h·k1)
  k1x = v(1);                                   % k1 de dx/dt = v
  k1v = -w^2 * x(1);                            % k1 de dv/dt = -ω²·x
  k2x = v(1) + h * k1v;                         % k2 de la posición
  k2v = -w^2 * (x(1) + h * k1x);                % k2 de la velocidad
  x(2) = x(1) + (h/2) * (k1x + k2x);            % posición en t2
  v(2) = v(1) + (h/2) * (k1v + k2v);            % velocidad en t2
  t(2) = t(1) + h;                              % segundo nodo de tiempo

  % ---- Bucle de interpolación de dos pasos ----
  % Aplicamos la fórmula del trapecio explícito a cada componente del
  % sistema, evaluando f en los dos nodos anteriores (k-1 y k).
  for k = 2:N
    % Para la posición: dx/dt = v, entonces f = v.
    x(k+1) = x(k) + (h/2) * (v(k-1) + v(k));
    % Para la velocidad: dv/dt = -ω²·x, entonces f = -ω²·x.
    v(k+1) = v(k) + (h/2) * (-w^2*x(k-1) - w^2*x(k));
    % Avanzamos el tiempo al siguiente nodo.
    t(k+1) = t(k) + h;
  end

end   % ============ fin de metodoInterpolacion ============

% =====================================================================
% SUBFUNCIÓN: metodoTaylor
%
% Resuelve el sistema  dx/dt = v ; dv/dt = -ω²·x  con el método de
% TAYLOR TRUNCADO de orden 2..5 (el mismo esquema del archivo
% TrunTaylor.m).
%
% La serie de Taylor para la posición x y la velocidad v es:
%
%   x(t+h) = x + h·x' + (h²/2!)·x'' + (h³/3!)·x''' + (h⁴/4!)·x'''' + ...
%   v(t+h) = v + h·v' + (h²/2!)·v'' + (h³/3!)·v''' + (h⁴/4!)·v'''' + ...
%
% Las derivadas se calculan derivando el sistema varias veces:
%   x'   = v         ,  v'   = -ω²·x
%   x''  = -ω²·x     ,  v''  = -ω²·v
%   x''' = -ω²·v     ,  v''' =  ω⁴·x
%   x''''=  ω⁴·x     ,  v''''=  ω⁴·v
% =====================================================================
function [t, x] = metodoTaylor(N, a, b, w, x0, v0, orden)

  h = (b - a) / N;                  % paso de tiempo Δt

  % Preasignamos los vectores para acelerar el bucle.
  t = zeros(1, N + 1);              % vector de tiempos
  x = zeros(1, N + 1);              % vector de posiciones
  v = zeros(1, N + 1);              % vector de velocidades

  % Condiciones iniciales.
  t(1) = a;                         % t_0 = a
  x(1) = x0;                        % x(0)
  v(1) = v0;                        % v(0)

  % ---- Bucle principal del método de Taylor truncado ----
  for k = 1:N
    xi = x(k);                      % posición actual x_k
    vi = v(k);                      % velocidad actual v_k

    % Incremento de la posición: empieza con el término de primer orden.
    incx = h * vi;
    % Incremento de la velocidad: primer orden, dv/dt = -ω²·x.
    incv = h * (-w^2 * xi);

    % Vamos añadiendo los términos de orden superior según el orden
    % pedido (igual que hace TrunTaylor.m con sus if).
    if orden >= 2
      incx = incx + (h^2/2) * (-w^2 * xi);    % x''/2! * h²
      incv = incv + (h^2/2) * (-w^2 * vi);    % v''/2! * h²
    end
    if orden >= 3
      incx = incx + (h^3/6) * (-w^2 * vi);    % x'''/3! * h³
      incv = incv + (h^3/6) * (w^4 * xi);     % v'''/3! * h³
    end
    if orden >= 4
      incx = incx + (h^4/24) * (w^4 * xi);    % x''''/4! * h⁴
      incv = incv + (h^4/24) * (w^4 * vi);    % v''''/4! * h⁴
    end
    if orden >= 5
      incx = incx + (h^5/120) * (w^4 * vi);   % x'''''/5! * h⁵
      incv = incv + (h^5/120) * (-w^6 * xi);  % v'''''/5! * h⁵
    end

    % Actualizamos la posición, la velocidad y el tiempo.
    x(k+1) = xi + incx;             % x_{k+1}
    v(k+1) = vi + incv;             % v_{k+1}
    t(k+1) = t(k) + h;              % t_{k+1} = t_k + h
  end

end   % ============ fin de metodoTaylor ============
