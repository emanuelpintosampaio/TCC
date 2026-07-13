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

#ANDREAS
T_a = 18           # Temperatura do ar em °C ANDREAS(18)
T_a_em_k = T_a + 273.15  # Temperatura do ar em K ANDREAS(291.15)
T_mar_em_k = 293.15  # Temperatura do mar em K ANDREAS(293.15)
T_gota = 19.8      # Temperatura da gota em °C ANDREAS(20.2)
T_gota_em_k = T_gota + 273.15  # Temperatura da gota em K ANDREAS(293.35)
P = 1000           # Pressão em mb
M_H2O = 18.016e-3  # Massa molecular da água em kg/mol
r_i = 100e-6       # Raio em metros
R = 8.31           # Constante universal dos gases
M_NaCl = 58.443e-3  # Massa molecular de sal na água em kg/mol
rho_w = 1025       # Densidade da água do mar kg/m^3
v_ion = 2          # Número de íons por molécula de NaCl dissociada
s = 34/1000        # Salinidade fracionária (34 psu)
f = 0.9            # Umidade relativa fracionária ANDREAS(0.8)
C_ar = 0.0154      # Concentração do gás carbônico no ar
S = 34             # Salinidade do mar
R_atm = 0.082      # Constante universal dos gases
g = 9.81           # Gravidade
v_ar = 1.32e-5     # Viscosidade cinemática do ar
rho_ar = 1.225     # Densidade do ar
H_s = 6            # Altura significativa da onda 
T0 = 273.15
P0 = 1013.25
c_ps = 4000        # Calor específico da spray     

# --- 2. CÁLCULOS PRELIMINARES E TEMPOS DE RELAXAMENTO ---
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
    T_a_em_k, alpha=alpha.calcular_alpha(T_a_em_k, a=17.502, b=240.97),
    beta=code_beta.calcular_beta(e_sat_val, T_a_em_k, L_v_val, M_H2O, D_linha_w_val, R, k_linha_a_val),
    b=240.97, exp_y=code_exp_y.calcular_exp_y(M_H2O, sigma_s_val, T_a_em_k, rho_w, R, v_ion, Phi_s_val, m_s_val, r_i, rho_s_val, M_NaCl),
    f=f)
rho_v_val = rho_v.calcular_rho_v(f, M_H2O, e_sat_val, R, T_a_em_k)
rho_vr_val = rho_vr.calcular_rho_vr(M_H2O, e_sat_esp.calcular_esat(T_a),
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


# --- 3. APROXIMAÇÃO EXPONENCIAL ---
def r_simplificado(t):
    return r_eq + (r_i - r_eq) * np.exp(-t / tau_r)

def T_simplificado(t):
    return T_eq_em_k + (T_gota_em_k - T_eq_em_k) * np.exp(-t / tau_T)


# --- 4. MODELO COMPLETO (SISTEMA DE EDOs) ---
def f_r_T_m(t, vars):
    r, T, m = vars
    rho_ww_novo    = m_s.calcular_rho_ww(T - 273.15)
    m_w_novo       = rho_s.calcular_massa_agua(r, rho_ww_novo)
    sigma_s_novo   = sigma_s.calculate_sigma_s(T - 273.15, m_s_val, m_w_novo)
    Phi_s_novo     = Phi_s.calcular_phi_s(m_s_val, M_H2O, m_w_novo)
    D_linha_w_novo = D_linha_w.calculate_Dw_prime(T, R, P, r, M_H2O, T0, P0, alpha_c=0.036, Delta_w=8e-8)
    L_v_novo       = L_v.calcular_lv(T - 273.15)
    k_linha_a_novo = k_linha_a.calculate_K_linha_a(T - 273.15, T, P, r, R, T0, P0, M_a=28.9644e-3, alpha_T=0.7, delta_T=2.16e-7, c_pa=1.006e3)
    rho_vr_novo    = rho_vr.calcular_rho_vr(
        M_H2O, e_sat_esp.calcular_esat(T - 273.15),
        code_exp_y.calcular_exp_y(M_H2O, sigma_s_novo, T_a_em_k, rho_w, R, v_ion, Phi_s_novo, m_s_val, r, rho_s_val, M_NaCl),
        R, T_a_em_k)

    Y = (2 * M_H2O * sigma_s_novo / (R * T_a_em_k * rho_w * r)) - (v_ion * Phi_s_novo * m_s_val * (M_H2O / M_NaCl) / (m_w_novo - m_s_val))
    den1 = rho_s_val * R * T_a_em_k / (D_linha_w_novo * M_H2O * e_sat_val)
    den2 = rho_s_val * L_v_novo / (k_linha_a_novo * T_a_em_k) * (L_v_novo * M_H2O / (R * T_a_em_k) - 1)
    den  = den1 + den2
    dr_dt = ((f - 1) - Y) / (r * den)

    dT_dt = 3 / (rho_s_val * c_ps * r**2) * (k_linha_a_novo * (T_a_em_k - T) + L_v_novo * D_linha_w_novo * (rho_v_val - rho_vr_novo))

    H_novo  = H_estrela.calcular_H_estrela(T, S)
    Dg_novo = Dg_estrela.calcular_Dg_estrela(r, T_a_em_k, R_atm)
    vol     = (4 / 3) * np.pi * r**3
    C_gota  = m / vol                                       
    dm_dt   = 4 * np.pi * r * Dg_novo * (C_ar - C_gota / (H_novo * R_atm * T))

    return np.array([dr_dt, dT_dt, dm_dt])

def euler_explicito_sistema(f_system, y0, t):
    num_pontos = len(t)
    n_vars     = len(y0)
    y          = np.zeros((num_pontos, n_vars))
    y[0]       = y0
    for n in range(num_pontos - 1):
        dt  = t[n + 1] - t[n]
        F   = f_system(t[n], y[n])
        y[n + 1] = y[n] + dt * F
    return y

# --- 5. DISCRETIZAÇÃO E RESOLUÇÃO ---
dt_fixo    = 1e-5
tempo_final = 1e3
n_steps    = int(tau_f / dt_fixo)
tempo = np.linspace(0, tempo_final, n_steps) # Iniciando em 0.01 para escala log

vol_i = (4 / 3) * np.pi * r_i**3          
H_i   = H_estrela.calcular_H_estrela(T_gota_em_k, S)
m_i = vol_i * C_ar * H_i * R_atm * T_gota_em_k
y0 = np.array([r_i, T_gota_em_k, m_i])

# Resolvendo Modelo Completo (Euler)
sol = euler_explicito_sistema(f_r_T_m, y0, tempo)
raio_euler        = sol[:, 0]
temperatura_euler = sol[:, 1]

# Resolvendo Aproximação Exponencial
raio_simp = r_simplificado(tempo)
temp_simp = T_simplificado(tempo)
print(T_eq)
print(r_eq)
print(tau_T)
print(tau_r)

# --- 6. PLOTAGEM ESTILO ARTIGO ---
import matplotlib.ticker as ticker

plt.rcParams['text.usetex'] = False
fig, ax1 = plt.subplots(figsize=(10, 7))

# Configurações do Eixo X (Tempo: 0.01 a 1000)
ax1.set_xscale('log')
ax1.set_xlabel('Time Since Formation (s)', fontsize=14)
ax1.set_xlim(0.01, 1000)
ax1.set_xticks([0.01, 0.1, 1, 10, 100, 1000])
ax1.get_xaxis().set_major_formatter(ticker.FormatStrFormatter('%g'))

# Cores e estilos de linha padronizados
cor_linha = 'black'
lw = 2

# Plotando Temperatura no eixo esquerdo (ax1: 16 a 21 °C)
ax1.plot(tempo, temperatura_euler - 273.15, '--', color=cor_linha, linewidth=lw, label='Full Model')
ax1.plot(tempo, temp_simp - 273.15, '-', color=cor_linha, linewidth=lw, label='Exponential Approximation')
ax1.set_ylabel('Droplet Temperature (°C)', fontsize=14, color=cor_linha)
ax1.set_ylim(16, 21)
ax1.set_yticks([16, 17, 18, 19, 20, 21])
ax1.tick_params(axis='y', labelcolor=cor_linha, labelsize=12)
ax1.tick_params(axis='x', labelsize=12)

# Criando eixo Y secundário para o Raio (ax2: 50 a 110 µm)
ax2 = ax1.twinx()
ax2.plot(tempo, raio_euler * 1e6, '--', color=cor_linha, linewidth=lw)
ax2.plot(tempo, raio_simp * 1e6, '-', color=cor_linha, linewidth=lw)
ax2.set_ylabel('Droplet Radius (µm)', fontsize=14, color=cor_linha)
ax2.set_ylim(50, 110)
ax2.set_yticks([50, 60, 70, 80, 90, 100, 110])
ax2.tick_params(axis='y', labelcolor=cor_linha, labelsize=12)

# Unificando a legenda na parte inferior
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, loc='lower left', fontsize=17, frameon=False)

# Adicionando a caixa de texto com os parâmetros
textstr = '\n'.join((
    f'$r_0 = {r_i*1e6:.0f}$ µm',
    f'$T_{{mar}} = {T_mar_em_k-273.15:.1f}°C$',
    f'$T_i = {T_gota:.1f}°C$',
    f'$T_a = {T_a:.1f}°C$',
    f'RH = {f*100:.0f}%'
))
props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none')
ax1.text(0.80, 0.95, textstr, transform=ax1.transAxes, fontsize=17,
        verticalalignment='top', bbox=props)

# Adicionando setas indicativas
# Seta esquerda para Temperatura
ax1.annotate('', xy=(0.08, 0.5), xytext=(0.22, 0.5), xycoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=10))

# Seta direita para Raio
ax2.annotate('', xy=(0.92, 0.5), xytext=(0.78, 0.5), xycoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=10))

# Ajustes de borda para simular o estilo do gráfico da imagem
for spine in ax1.spines.values():
    spine.set_linewidth(1.5)

ax1.tick_params(direction='in', length=6, width=1.5, bottom=True, top=True, left=True, right=False)
ax2.tick_params(direction='in', length=6, width=1.5, left=False, right=True)

plt.tight_layout()
plt.show()