import numpy as np
from NEWTON import metodo_newton
from T_EQ.BETA import code_beta, D_linha_w, e_sat, k_linha_a, L_v
from T_EQ.EXP_Y import code_exp_y, m_s, Phi_s, rho_s, sigma_s
from T_EQ import Delta_T_contas, alpha
from TAU_T import partial_rho_vr, rho_vr, contas_tau_T, e_sat_esp, rho_v
from TAU_R.partial_rt_teste import calcular_tau_r
from DGH import Dg_estrela, H_estrela
import matplotlib.pyplot as plt
import math

# ANDREAS
T_a = 18          
T_a_em_k = T_a + 273.15 
T_mar_em_k = 293.15  
T_gota = 19.8     
T_gota_em_k = T_gota + 273.15  
P = 1000          
M_H2O = 18.016e-3  
r_i = 30e-6       
R = 8.31           
M_NaCl = 58.443e-3  
rho_w = 1025       
v_ion = 2         
s = 34/1000       
f = 0.9            
C_ar = 0.0154      
S = 34             
R_atm = 0.082     
g = 9.81          
v_ar = 1.32e-5    
rho_ar = 1.225     
H_s = 6            
T0 = 273.15
P0 = 1013.25
c_ps = 4000       

# Cálculos Iniciais
def func(U_f):
    const = (2 * r_i**2 * g) / (9 * v_ar) * ((rho_w / rho_ar) - 1)
    a = 1 + 0.158 * ((2 * r_i * U_f / v_ar) ** (2 / 3))
    return U_f - const / a

def dfunc(U_f):
    const = (2 * r_i**2 * g) / (9 * v_ar) * ((rho_w / rho_ar) - 1)
    a = 1 + 0.158 * ((2 * r_i * U_f / v_ar) ** (2 / 3))
    d_a = -0.158 * (2 / 3) * (2 * r_i / v_ar) ** (2 / 3) * U_f ** (-1 / 3)
    return 1 - (const * d_a) / (a**2)

U_f0 = 0.1
U_f_sol = metodo_newton.Newton(func, dfunc, U_f0, 1e-6, 50)

tau_f = H_s / (2 * U_f_sol)

e_sat_val = e_sat.calcular_esat(T_a)
L_v_val = L_v.calcular_lv(T_gota)
D_linha_w_val = D_linha_w.calculate_Dw_prime(T_gota_em_k, R, P, r_i, M_H2O, T0, P0, alpha_c=0.036, Delta_w=8e-8)
k_linha_a_val = k_linha_a.calculate_K_linha_a(T_gota, T_gota_em_k, P, r_i, R, T0, P0, M_a=28.9644e-3, alpha_T=0.7, delta_T=2.16e-7, c_pa=1.006e3)

rho_ww = m_s.calcular_rho_ww(T_gota)
m_s_val = m_s.calcular_ms(rho_ww, r_i, s)
v_a = rho_s.calcular_v_a(T_gota, m_s_val, M_NaCl, r_i)
m_w = rho_s.calcular_massa_agua(r_i, rho_ww)
sigma_s_val = sigma_s.calculate_sigma_s(T_gota, m_s_val, m_w)
Phi_s_val = Phi_s.calcular_phi_s(m_s_val, M_H2O, m_w)
rho_s_val = rho_s.calcular_rho_spray(rho_ww, m_s_val, m_w, v_a, M_NaCl)

Delta_T = Delta_T_contas.calculate_delta_T(
    T_a_em_k,
    alpha=alpha.calcular_alpha(T_a_em_k, a=17.502, b=240.97),
    beta=code_beta.calcular_beta(e_sat_val, T_a_em_k, L_v_val, M_H2O, D_linha_w_val, R, k_linha_a_val),
    b=240.97,
    exp_y=code_exp_y.calcular_exp_y(M_H2O, sigma_s_val, T_a_em_k, rho_w, R, v_ion, Phi_s_val, m_s_val, r_i, rho_s_val, M_NaCl),
    f=f)
rho_v_val = rho_v.calcular_rho_v(f, M_H2O, e_sat_val, R, T_a_em_k)

rho_vr_val = rho_vr.calcular_rho_vr(M_H2O,
    e_sat_esp.calcular_esat(T_a),
    code_exp_y.calcular_exp_y(M_H2O, sigma_s_val, T_a_em_k, rho_w, R, v_ion, Phi_s_val, m_s_val, r_i, rho_s_val, M_NaCl),
    R,T_a_em_k)

partial_rho_vr_val = partial_rho_vr.calcular_partial_rho_vr(T_a_em_k, rho_vr_val, a=17.502, b=240.97)

tau_T = contas_tau_T.calculate_tau_T(rho_s_val, r_i, k_linha_a_val, D_linha_w_val, partial_rho_vr_val, L_v_val, c_ps)

T_eq = Delta_T[0] + T_a
T_eq_em_k = T_eq + 273.15

def zeta(r_i):
    term1 = (f - 1)
    term2 = (2 * M_H2O * sigma_s_val) / (R * T_eq_em_k * rho_w * r_i)
    denominator = (4 * np.pi * rho_s_val * (r_i**3 / 3)) - m_s_val
    term3 = (v_ion * Phi_s_val * m_s_val * (M_H2O / M_NaCl)) / denominator
    return term1 - term2 + term3

def dzeta_dr(r_i):
    term1 = (2 * M_H2O * sigma_s_val) / (R * T_eq_em_k * rho_w * (r_i**2))
    denominator = (4 * np.pi * rho_s_val * (r_i**3)) / 3 - m_s_val
    term2 = (v_ion * Phi_s_val * m_s_val * (M_H2O / M_NaCl) * 4 * np.pi * rho_s_val * (r_i**2)) / (denominator**2)
    return term1 - term2

r0 = (1.1) * (3 * m_s_val / (4 * np.pi * rho_s_val)) ** (1 / 3)

r_eq = metodo_newton.Newton(zeta, dzeta_dr, r0, 1e-6, 100)

tau_r = calcular_tau_r(f, M_H2O, sigma_s_val, R, T_a_em_k, rho_w, r_i, rho_s_val, m_s_val, v_ion, Phi_s_val, M_NaCl, D_linha_w_val, e_sat_val, L_v_val, k_linha_a_val, r_eq)[1]

# Simplificações exponenciais para r(t) e T(t) 
def r_simplificado(t):
    return r_eq + (r_i - r_eq) * np.exp(-t / tau_r)

def T_simplificado(t):
    return T_eq_em_k + (T_gota_em_k - T_eq_em_k) * np.exp(-t / tau_T)

# EDO da massa simplificada (apenas a massa varia com o integrador)
def f_m(t, m):
    r = r_simplificado(t)
    T = T_simplificado(t)

    H_t   = H_estrela.calcular_H_estrela(T, S)
    Dg_t  = Dg_estrela.calcular_Dg_estrela(r, T_a_em_k, R_atm)
    vol   = (4 / 3) * np.pi * r**3
    C_gota = m / vol

    dm_dt = 4 * np.pi * r * Dg_t * (C_ar - C_gota / (H_t * R_atm * T))
    return dm_dt

# Método de Euler Explícito Adaptativo (com controle PID) focado na massa
def euler_adaptativo_1d(f_m, m0, t0, t_final, dt_min, dt_max, dt_inicial, tol, K_P, K_I, K_D):
    t_hist = [t0]
    m_hist = [m0]
    dt_hist = [dt_inicial]

    dt = dt_inicial
    dt_prev = dt_min
    n_rej = 0  
    
    e_n_1 = tol
    e_n_2 = tol

    m_n = float(m0)
    t_n = t0

    while t_n < t_final:
        if t_n + dt > t_final:
            dt = t_final - t_n

        # Passo de Euler Explícito 1D
        F1 = f_m(t_n, m_n)
        m_next = m_n + dt * F1

        # Estimativa do erro baseada no incremento relativo do passo
        if m_next != 0:
            e_n = abs(m_next - m_n) / abs(m_next)
        else:
            e_n = 1e-16 # Evita divisão por zero

        # Controle PID 
        if e_n > tol and dt > dt_min: 
            n_rej += 1 
            fator = 1 / e_n 
            if fator > 0.8: 
                fator = 0.8 
            
            dt = max(fator * dt, dt_min) 
            dt_prev = (dt**2) / dt_prev 
            
            continue 
            
        else:
            if e_n == 0:
                e_n = 1e-16

            t_n += dt
            m_n = m_next

            t_hist.append(t_n)
            m_hist.append(m_n)
            dt_hist.append(dt)

            # Eq 9 PID
            fator_P = (e_n_1 / e_n) ** K_P 
            fator_I = (tol / e_n) ** K_I 
            fator_D = ((e_n_1**2) / (e_n * e_n_2)) ** K_D

            dt_next = fator_P * fator_I * fator_D * dt 

            dt = max(dt_next, dt_min) 
            dt = min(dt, dt_max) 

            dt_prev = dt 
            
            e_n_2 = e_n_1 
            e_n_1 = e_n 

    return np.array(t_hist), np.array(m_hist), np.array(dt_hist) 

# Configurações do Euler Adaptativo
tempo_final = tau_f
dt_min      = 1e-6
dt_max      = 1         # Passo máximo permitido
dt_inicial  = 1e-4      # Inicial
tol         = 1e-4      # Tolerância do erro local

# Ganhos do Controlador PID
K_P = 0.075
K_I = 0.175
K_D = 0.01

# Condição inicial da massa 
vol_i = (4 / 3) * np.pi * r_i**3
H_i   = H_estrela.calcular_H_estrela(T_gota_em_k, S)
m_i   = vol_i * C_ar * H_i * R_atm * T_gota_em_k

# Solução com Euler Adaptativo para massa
tempo, massa_euler, dts = euler_adaptativo_1d(
    f_m=f_m, 
    m0=m_i, 
    t0=0.0, 
    t_final=tempo_final, 
    dt_min=dt_min, 
    dt_max=dt_max, 
    dt_inicial=dt_inicial, 
    tol=tol, 
    K_P=K_P, 
    K_I=K_I, 
    K_D=K_D
)

# Reconstrução dos vetores de raio e temperatura
raio_simp = r_simplificado(tempo)
temp_simp = T_simplificado(tempo)

massa_final = massa_euler[-1]
dt_final = dts[-1]

#print(f'raio inicial: {r_i}')
#print(f'massa final: {massa_final}')
#print(f'pontos: {len(tempo)}')


# Gráficos 
plt.rcParams['text.usetex'] = False
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 10), sharex=True)


fig.suptitle(
    f'Euler Explícito | Modelo Simplificado | Adaptativo {len(tempo)} pontos',
    fontsize=13
)

# Raio
ax1.plot(tempo, raio_simp * 1e6, 's-', color="#FF0000", lw=2, ms=4,
         label=f'r_final = {raio_simp[-1] * 1e6:.2f} µm')
ax1.set_ylabel('Raio da Gota (µm)', fontsize=12)
ax1.set_xscale('log')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Temperatura
ax2.plot(tempo, temp_simp - 273.15, 's-', color='#FF0000', lw=2, ms=4,
          label=f'T_final = {temp_simp[-1] - 273.15:.2f} °C')
ax2.set_ylabel('Temperatura da Gota (°C)', fontsize=12)
ax2.set_xscale('log')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

# Massa
ax3.semilogx(tempo, massa_euler, 's-', color='#FF0000', lw=2, ms=4,
              label=f'Massa final: {massa_final:.3e} mol')
ax3.set_ylabel('Massa (mol)', fontsize=12)
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)

# dt
ax4.semilogx(tempo, dts, 's-', color='#FF0000', lw=2, ms=4,
               label=f"dt final: {dt_final:.2e} s")
ax4.set_xlabel('Tempo (s)', fontsize=12)
ax4.set_ylabel('Passo de tempo (s)', fontsize=12)
ax4.grid(True, alpha=0.3, which='both')
ax4.legend(fontsize=10)
ax4.set_yscale('log')

plt.tight_layout()
#plt.show()
plt.savefig("E5_30.png")
