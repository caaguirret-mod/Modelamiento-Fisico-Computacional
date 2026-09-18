import numpy as np
import matplotlib.pyplot as plt

# --- solucionador no lineal ---
def resolver_liouville_no_lineal(Nz, z_min=1.0, z_max=3.0, t_final=1.0, cfl=0.5, eps=0.0, z0=2.0, sigma=0.15, m=1.0, b=1.0):
    coef_nl = 2.0 * m * (b**2)
    z = np.linspace(z_min, z_max, Nz)
    dz = z[1] - z[0]
    
    Nt = int(np.ceil(t_final / (cfl * dz))) + 1
    dt = t_final / (Nt - 1)
    t_vec = np.linspace(0.0, t_final, Nt)
    
    mu2 = (dt / dz)**2
    dt2 = dt**2

    # --- fondo de vacio ---
    rho_bg = -0.5 * np.log(coef_nl) - np.log(z)

    # --- condiciones iniciales ---
    rho = np.zeros((Nz, Nt))
    delta_rho_0 = eps * np.exp(-((z - z0)**2) / (sigma**2))
    rho[:, 0] = rho_bg + delta_rho_0

    # --- paso n=1 (taylor) ---
    for j in range(1, Nz - 1):
        lap_z = (rho[j+1, 0] - 2.0 * rho[j, 0] + rho[j-1, 0]) / (dz**2)
        fuente = coef_nl * np.exp(2.0 * rho[j, 0])
        rho[j, 1] = rho[j, 0] + 0.5 * dt2 * (lap_z - fuente)
        
    # --- condiciones de frontera n=1 ---
    rho[0, 1] = rho_bg[0]
    rho[-1, 1] = rho_bg[-1]

    # --- bucle principal ---
    for n in range(1, Nt - 1):
        for j in range(1, Nz - 1):
            r_curr = rho[j, n]
            r_prev = rho[j, n - 1]
            lap_esp = mu2 * (rho[j+1, n] - 2.0 * r_curr + rho[j-1, n])
            gamma = coef_nl * dt2 * np.exp(2.0 * r_curr)
            
            # --- stencil ---
            # nota mental para steven del futuro: revisar el termino no lineal antes que todo que casi nos explota la simulación
            num = 2.0 * r_curr - r_prev + lap_esp - gamma * (1.0 - 2.0 * r_curr)
            den = 1.0 + 2.0 * gamma
            rho[j, n+1] = num / den

        # --- condiciones de frontera ---
        rho[0, n+1] = rho_bg[0]
        rho[-1, n+1] = rho_bg[-1]

    return z, t_vec, rho, rho_bg, dz, dt

# --- solucionador lineal ---
def resolver_liouville_lineal(Nz, z_min=1.0, z_max=3.0, t_final=0.6, cfl=0.5, eps=0.10, z0=2.0, sigma=0.15):
    z = np.linspace(z_min, z_max, Nz)
    dz = z[1] - z[0]
    Nt = int(np.ceil(t_final / (cfl * dz))) + 1
    dt = t_final / (Nt - 1)
    t_vec = np.linspace(0.0, t_final, Nt)
    mu2 = (dt / dz)**2
    dt2 = dt**2

    # --- condiciones iniciales lineal ---
    delta_rho = np.zeros((Nz, Nt))
    delta_rho[:, 0] = eps * np.exp(-((z - z0)**2) / (sigma**2))
    V_bg = 2.0 / (z**2)

    # --- paso n=1 lineal ---
    for j in range(1, Nz - 1):
        lap_z = (delta_rho[j+1, 0] - 2.0 * delta_rho[j, 0] + delta_rho[j-1, 0]) / (dz**2)
        fuente_lin = V_bg[j] * delta_rho[j, 0]
        delta_rho[j, 1] = delta_rho[j, 0] + 0.5 * dt2 * (lap_z - fuente_lin)

    # --- condiciones de frontera n=1 ---
    delta_rho[0, 1] = 0.0
    delta_rho[-1, 1] = 0.0

    # --- bucle principal lineal ---
    for n in range(1, Nt - 1):
        for j in range(1, Nz - 1):
            curv_esp = mu2 * (delta_rho[j+1, n] - 2.0 * delta_rho[j, n] + delta_rho[j-1, n])
            termino_pot = dt2 * V_bg[j] * delta_rho[j, n]
            delta_rho[j, n+1] = 2.0 * delta_rho[j, n] - delta_rho[j, n-1] + curv_esp - termino_pot

        # --- condiciones de frontera lineal ---
        delta_rho[0, n+1] = 0.0
        delta_rho[-1, n+1] = 0.0

    return z, t_vec, delta_rho

# --- diagnostico del residuo numerico ---
def calcular_residuo_liouville(rho, dz, dt, m=1.0, b=1.0):
    Nz, Nt = rho.shape
    coef_nl = 2.0 * m * (b**2)
    R_num = np.full((Nz, Nt), -coef_nl)

    # --- pasos interiores: diferencia central en tiempo ---
    # nota mental: solo se calcula en pasos interiores para garantizar consistencia del stencil
    for n in range(1, Nt - 1):
        for j in range(1, Nz - 1):
            rho_tt = (rho[j, n+1] - 2.0 * rho[j, n] + rho[j, n-1]) / (dt**2)
            rho_zz = (rho[j+1, n] - 2.0 * rho[j, n] + rho[j-1, n]) / (dz**2)
            R_num[j, n] = np.exp(-2.0 * rho[j, n]) * (rho_tt - rho_zz)

    residuo = R_num - (-coef_nl)
    return R_num, residuo

# --- validacion de vacio (fase 1) ---
z_vac, t_vac, rho_vac, rho_bg_vac, dz_vac, dt_vac = resolver_liouville_no_lineal(Nz=100, t_final=1.0, eps=0.0)

pasos_plot = [1, len(t_vac)//3, 2*len(t_vac)//3, len(t_vac)-1]
fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(12, 4.5))

print("--- metricas de vacio ---")
for idx, n_step in enumerate(pasos_plot):
    t_actual = t_vac[n_step]
    error_Linf = np.max(np.abs(rho_vac[:, n_step] - rho_bg_vac))
    print(f"Paso {n_step:<4} | t = {t_actual:.3f} s | Error L_inf: {error_Linf:.6e}")

ax1a.plot(z_vac, rho_bg_vac, 'k--', lw=1.5, label=r'Analítica $\rho_{\mathrm{bg}}(z)$')
ax1a.plot(z_vac, rho_vac[:, -1], 'b.-', markersize=4, label=r'Numérica ($t = 1.0$ s)')
ax1a.set_title(r'Validación del Vacío $\mathrm{AdS}_2$ ($\varepsilon = 0$)')
ax1a.set_xlabel('Coordenada espacial $z$')
ax1a.set_ylabel(r'Factor conforme $\rho(z)$')
ax1a.grid(True, ls=':', alpha=0.6)
ax1a.legend()

# --- convergencia formal ---
mallas_est = [50, 100, 200, 400]
err_est_Linf, err_est_L2, dz_est_vals = [], [], []

print("\n--- metricas de convergencia estatica ---")
for N in mallas_est:
    z_s, t_s, rho_s, bg_s, dz_val, _ = resolver_liouville_no_lineal(Nz=N, t_final=1.0, eps=0.0)
    err_fin = np.abs(rho_s[:, -1] - bg_s)
    err_est_Linf.append(np.max(err_fin))
    err_est_L2.append(np.sqrt(dz_val * np.sum(err_fin**2)))
    dz_est_vals.append(dz_val)
    print(f"Nz: {N:<4} | Delta z: {dz_val:.5f} | L_inf: {np.max(err_fin):.6e}")

p_est_inf = np.polyfit(np.log(dz_est_vals), np.log(err_est_Linf), 1)[0]
p_est_l2  = np.polyfit(np.log(dz_est_vals), np.log(err_est_L2), 1)[0]
print(f"Pendiente L_inf: {p_est_inf:.3f} | Pendiente L_2: {p_est_l2:.3f}")

ax1b.loglog(dz_est_vals, err_est_Linf, 'ro-', label=f'Norma $L_\\infty$ (Pendiente: {p_est_inf:.2f})')
ax1b.loglog(dz_est_vals, err_est_L2, 'bs--', label=f'Norma $L_2$ (Pendiente: {p_est_l2:.2f})')
ax1b.set_title('Convergencia Espacial Formal')
ax1b.set_xlabel(r'$\Delta z$')
ax1b.set_ylabel('Error numérico residual')
ax1b.grid(True, which="both", ls=":", alpha=0.6)
ax1b.legend()

plt.tight_layout()
fig1.savefig("1_validacion_y_convergencia.png", dpi=300)
plt.close(fig1)

# --- dinamica de la perturbacion (fase 2) ---
Nz_din = 150
z_din, t_din, rho_din, rho_bg_din, dz_din, dt_din = resolver_liouville_no_lineal(Nz=Nz_din, t_final=0.7, eps=0.10)
delta_rho = rho_din - rho_bg_din[:, None]

fig2 = plt.figure(figsize=(13, 5))
ax2a = fig2.add_subplot(121)
pasos_din = [0, len(t_din)//3, 2*len(t_din)//3, len(t_din)-2]
for p in pasos_din:
    ax2a.plot(z_din, delta_rho[:, p], lw=1.6, label=f'$t = {t_din[p]:.3f}$ s')
ax2a.set_title(r'Separación D\'Alembertiana de Modos ($\delta\rho$)')
ax2a.set_xlabel('Coordenada $z$')
ax2a.set_ylabel(r'$\delta\rho(z, t)$')
ax2a.grid(True, ls=':', alpha=0.6)
ax2a.legend()

# --- superficie 3d perturbacion ---
# nota para steven del futuro: aprender a graficar esta vaina que es? media hora acomodando las proyecciones 3D
ax2b = fig2.add_subplot(122, projection='3d')
Z_din_mesh, T_din_mesh = np.meshgrid(z_din, t_din)
surf2 = ax2b.plot_surface(Z_din_mesh, T_din_mesh, delta_rho.T, cmap='magma', edgecolor='none')
ax2b.set_title(r'Superficie Espaciotemporal de $\delta\rho(z, t)$')
ax2b.set_xlabel('$z$')
ax2b.set_ylabel('$t$')
ax2b.set_zlabel(r'$\delta\rho$')
fig2.colorbar(surf2, ax=ax2b, shrink=0.5, aspect=10)

plt.tight_layout()
fig2.savefig("2_dinamica_conforme.png", dpi=300)
plt.close(fig2)

# --- residuo numerico de liouville (fase 3) ---
R_num, residuo = calcular_residuo_liouville(rho_din, dz_din, dt_din)

fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(13, 5))
for p in pasos_din[1:]:
    ax3a.plot(z_din[1:-1], residuo[1:-1, p], lw=1.6, label=f'$t = {t_din[p]:.3f}$ s')
ax3a.set_title(r'Residuo Numérico Local $\mathcal{R}(z, t) = \tilde{R}_{\mathrm{num}} - (-2mb^2)$')
ax3a.set_xlabel('Coordenada $z$')
ax3a.set_ylabel('Error de truncamiento local')
ax3a.grid(True, ls=':', alpha=0.6)
ax3a.legend()

c_map = ax3b.contourf(Z_din_mesh, T_din_mesh, residuo.T, levels=30, cmap='coolwarm')
ax3b.set_title('Mapa Espaciotemporal del Residuo del Esquema')
ax3b.set_xlabel('Espacio $z$')
ax3b.set_ylabel('Tiempo $t$')
fig3.colorbar(c_map, ax=ax3b)

plt.tight_layout()
fig3.savefig("3_residuo_numerico_liouville.png", dpi=300)
plt.close(fig3)

# --- experimento a: lineal vs no lineal ---
amplitudes = [0.01, 0.05, 0.10, 0.20, 0.40]
fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(13, 4.8))

print("\n--- metricas experimento a (lineal vs no lineal) ---")
for eps_val in amplitudes:
    z_a, t_a, rho_nl, bg_a, _, _ = resolver_liouville_no_lineal(Nz=150, t_final=0.6, eps=eps_val)
    d_rho_nl = rho_nl - bg_a[:, None]
    _, _, d_rho_lin = resolver_liouville_lineal(Nz=150, t_final=0.6, eps=eps_val)
    E_nl_t = np.max(np.abs(d_rho_nl - d_rho_lin), axis=0)
    ax4a.plot(t_a, E_nl_t, lw=1.6, label=rf'$\varepsilon = {eps_val}$')
    print(f"Eps: {eps_val:<4.2f} | Error E_NL max: {np.max(E_nl_t):.6e}")

ax4a.set_title(r'Discrepancia $E_{\mathrm{NL}}(t) = \max_z |\delta\rho_{\mathrm{NL}} - \delta\rho_{\mathrm{L}}|$')
ax4a.set_xlabel('Tiempo $t$ [s]')
ax4a.set_ylabel(r'$E_{\mathrm{NL}}$')
ax4a.grid(True, ls=':', alpha=0.6)
ax4a.legend()

# --- experimento c: independencia de fronteras ---
z_max_list = [3.0, 5.0, 8.0]
for zm in z_max_list:
    Nz_c = int(75 * (zm - 1.0))
    z_c, t_c, rho_c, bg_c, _, _ = resolver_liouville_no_lineal(Nz=Nz_c, z_min=1.0, z_max=zm, t_final=0.8, eps=0.10)
    d_rho_c = rho_c - bg_c[:, None]
    idx_centro = np.argmin(np.abs(z_c - 2.0))
    ax4b.plot(t_c, d_rho_c[idx_centro, :], lw=1.6, label=rf'$z_{{\max}} = {zm:.1f}$')

ax4b.set_title(r'Evolución en $z = 2.0$ (Independencia de Fronteras)')
ax4b.set_xlabel('Tiempo $t$ [s]')
ax4b.set_ylabel(r'$\delta\rho(z=2.0, t)$')
ax4b.grid(True, ls=':', alpha=0.6)
ax4b.legend()

plt.tight_layout()
fig4.savefig("4_limites_fisicos_y_frontera.png", dpi=300)
plt.close(fig4)

# --- guardado estructurado de datos ---
# nota para steven del futuro: hacer un hijuemadre repositorio en git antes de perder todo este avance
np.savez("datos_simulacion_liouville.npz", 
         z=z_din, t=t_din, rho=rho_din, rho_bg=rho_bg_din, residuo=residuo)

# nota final: si voy a tener esa falta de criterio para escoger tema mejor me dedico a contar moleculas d eoxigeno
print("\n--- simulacion finalizada y graficas guardadas con exito ---")