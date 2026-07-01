import numpy as np
from NEWTON import metodo_newton
from T_EQ.BETA import code_beta, D_linha_w, e_sat, k_linha_a, L_v
from T_EQ.EXP_Y import code_exp_y, m_s, Phi_s, rho_s, sigma_s
from T_EQ import Delta_T_contas, alpha
from TAU_T import partial_rho_vr, rho_vr, contas_tau_T, e_sat_esp, rho_v
from TAU_R.partial_rt_teste import calcular_tau_r
from DGH import Dg_estrela, H_estrela
import matplotlib.pyplot as plt


T_a         = 18
T_a_em_k    = T_a + 273.15
T_mar_em_k  = 293.15
T_gota      = 19.8
T_gota_em_k = T_gota + 273.15
P           = 1000
M_H2O       = 18.016e-3
r_i         = 500e-6
R           = 8.31
M_NaCl      = 58.443e-3
rho_w       = 1025
v_ion       = 2
s           = 34 / 1000
f           = 0.9
C_ar        = 0.0154
S           = 34
R_atm       = 0.082
g           = 9.81
v_ar        = 1.32e-5
rho_ar      = 1.225
H_s         = 6
T0          = 273.15
P0          = 1013.25
c_ps        = 4000


def func(U_f):
    const = (2 * r_i**2 * g) / (9 * v_ar) * ((rho_w / rho_ar) - 1)
    a = 1 + 0.158 * ((2 * r_i * U_f / v_ar) ** (2 / 3))
    return U_f - const / a

def dfunc(U_f):
    const = (2 * r_i**2 * g) / (9 * v_ar) * ((rho_w / rho_ar) - 1)
    a = 1 + 0.158 * ((2 * r_i * U_f / v_ar) ** (2 / 3))
    d_a = -0.158 * (2 / 3) * (2 * r_i / v_ar) ** (2 / 3) * U_f ** (-1 / 3)
    return 1 - (const * d_a) / (a**2)

U_f_sol = metodo_newton.Newton(func, dfunc, 0.1, 1e-6, 50)
tau_f   = H_s / (2 * U_f_sol)


e_sat_val       = e_sat.calcular_esat(T_a)
L_v_val         = L_v.calcular_lv(T_gota)
D_linha_w_val   = D_linha_w.calculate_Dw_prime(T_gota_em_k, R, P, r_i, M_H2O, T0, P0,
                                                 alpha_c=0.036, Delta_w=8e-8)
k_linha_a_val   = k_linha_a.calculate_K_linha_a(T_gota, T_gota_em_k, P, r_i, R, T0, P0,
                                                  M_a=28.9644e-3, alpha_T=0.7,
                                                  delta_T=2.16e-7, c_pa=1.006e3)
rho_ww          = m_s.calcular_rho_ww(T_gota)
m_s_val         = m_s.calcular_ms(rho_ww, r_i, s)
v_a             = rho_s.calcular_v_a(T_gota, m_s_val, M_NaCl, r_i)
m_w             = rho_s.calcular_massa_agua(r_i, rho_ww)
sigma_s_val     = sigma_s.calculate_sigma_s(T_gota, m_s_val, m_w)
Phi_s_val       = Phi_s.calcular_phi_s(m_s_val, M_H2O, m_w)
rho_s_val       = rho_s.calcular_rho_spray(rho_ww, m_s_val, m_w, v_a, M_NaCl)

Delta_T = Delta_T_contas.calculate_delta_T(
    T_a_em_k,
    alpha   = alpha.calcular_alpha(T_a_em_k, a=17.502, b=240.97),
    beta    = code_beta.calcular_beta(e_sat_val, T_a_em_k, L_v_val, M_H2O,
                                       D_linha_w_val, R, k_linha_a_val),
    b       = 240.97,
    exp_y   = code_exp_y.calcular_exp_y(M_H2O, sigma_s_val, T_a_em_k, rho_w, R,
                                         v_ion, Phi_s_val, m_s_val, r_i, rho_s_val, M_NaCl),
    f       = f)
rho_v_val = rho_v.calcular_rho_v(f, M_H2O, e_sat_val, R, T_a_em_k)
rho_vr_val = rho_vr.calcular_rho_vr(
    M_H2O,
    e_sat_esp.calcular_esat(T_a),
    code_exp_y.calcular_exp_y(M_H2O, sigma_s_val, T_a_em_k, rho_w, R,
                               v_ion, Phi_s_val, m_s_val, r_i, rho_s_val, M_NaCl),
    R, T_a_em_k)
partial_rho_vr_val = partial_rho_vr.calcular_partial_rho_vr(
    T_a_em_k, rho_vr_val, a=17.502, b=240.97)
tau_T = contas_tau_T.calculate_tau_T(
    rho_s_val, r_i, k_linha_a_val, D_linha_w_val, partial_rho_vr_val, L_v_val, c_ps)

T_eq    = Delta_T[0] + T_a
T_eq_em_k = T_eq + 273.15


def zeta(r):
    term1 = (f - 1)
    term2 = (2 * M_H2O * sigma_s_val) / (R * T_eq_em_k * rho_w * r)
    den   = (4 * np.pi * rho_s_val * (r**3 / 3)) - m_s_val
    term3 = (v_ion * Phi_s_val * m_s_val * (M_H2O / M_NaCl)) / den
    return term1 - term2 + term3

def dzeta_dr(r):
    term1 = (2 * M_H2O * sigma_s_val) / (R * T_eq_em_k * rho_w * r**2)
    den   = (4 * np.pi * rho_s_val * r**3) / 3 - m_s_val
    term2 = (v_ion * Phi_s_val * m_s_val * (M_H2O / M_NaCl)
             * 4 * np.pi * rho_s_val * r**2) / (den**2)
    return term1 - term2

r0_newton = 1.1 * (3 * m_s_val / (4 * np.pi * rho_s_val)) ** (1 / 3)
r_eq = metodo_newton.Newton(zeta, dzeta_dr, r0_newton, 1e-6, 100)

tau_r = calcular_tau_r(f, M_H2O, sigma_s_val, R, T_a_em_k, rho_w, r_i, rho_s_val,
                        m_s_val, v_ion, Phi_s_val, M_NaCl,
                        D_linha_w_val, e_sat_val, L_v_val, k_linha_a_val, r_eq)[1]


vol_i = (4 / 3) * np.pi * r_i**3
H_i   = H_estrela.calcular_H_estrela(T_gota_em_k, S)
m_i   = vol_i * C_ar * H_i * R_atm * T_gota_em_k
y0    = np.array([r_i, T_gota_em_k, m_i])


def f_r_T_m(t, vars):
    r, T, m = vars
    rho_ww_n    = m_s.calcular_rho_ww(T - 273.15)
    m_w_n       = rho_s.calcular_massa_agua(r, rho_ww_n)
    sigma_s_n   = sigma_s.calculate_sigma_s(T - 273.15, m_s_val, m_w_n)
    Phi_s_n     = Phi_s.calcular_phi_s(m_s_val, M_H2O, m_w_n)
    D_linha_w_n = D_linha_w.calculate_Dw_prime(T, R, P, r, M_H2O, T0, P0,
                                                 alpha_c=0.036, Delta_w=8e-8)
    L_v_n       = L_v.calcular_lv(T - 273.15)
    k_linha_a_n = k_linha_a.calculate_K_linha_a(T - 273.15, T, P, r, R, T0, P0,
                                                  M_a=28.9644e-3, alpha_T=0.7,
                                                  delta_T=2.16e-7, c_pa=1.006e3)
    rho_vr_n    = rho_vr.calcular_rho_vr(
        M_H2O,
        e_sat_esp.calcular_esat(T - 273.15),
        code_exp_y.calcular_exp_y(M_H2O, sigma_s_n, T_a_em_k, rho_w, R,
                                   v_ion, Phi_s_n, m_s_val, r, rho_s_val, M_NaCl),
        R, T_a_em_k)
    Y    = ((2 * M_H2O * sigma_s_n / (R * T_a_em_k * rho_w * r))
            - (v_ion * Phi_s_n * m_s_val * (M_H2O / M_NaCl) / (m_w_n - m_s_val)))
    den1 = rho_s_val * R * T_a_em_k / (D_linha_w_n * M_H2O * e_sat_val)
    den2 = (rho_s_val * L_v_n / (k_linha_a_n * T_a_em_k)
            * (L_v_n * M_H2O / (R * T_a_em_k) - 1))
    den  = den1 + den2
    dr_dt = ((f - 1) - Y) / (r * den)
    dT_dt = (3 / (rho_s_val * c_ps * r**2)) * (
        k_linha_a_n * (T_a_em_k - T) + L_v_n * D_linha_w_n * (rho_v_val - rho_vr_n))
    H_n   = H_estrela.calcular_H_estrela(T, S)
    Dg_n  = Dg_estrela.calcular_Dg_estrela(r, T_a_em_k, R_atm)
    vol   = (4 / 3) * np.pi * r**3
    dm_dt = 4 * np.pi * r * Dg_n * (C_ar - (m / vol) / (H_n * R_atm * T))
    return np.array([dr_dt, dT_dt, dm_dt])


def r_simplificado(t):
    return r_eq + (r_i - r_eq) * np.exp(-t / tau_r)

def T_simplificado(t):
    return T_eq_em_k + (T_gota_em_k - T_eq_em_k) * np.exp(-t / tau_T)

def f_m(t, m):
    r      = r_simplificado(t)
    T      = T_simplificado(t)
    H_t    = H_estrela.calcular_H_estrela(T, S)
    Dg_t   = Dg_estrela.calcular_Dg_estrela(r, T_a_em_k, R_atm)
    vol    = (4 / 3) * np.pi * r**3
    C_gota = m / vol
    return 4 * np.pi * r * Dg_t * (C_ar - C_gota / (H_t * R_atm * T))


def _grandezas_dinamicas(r, T):
    rho_ww_n    = m_s.calcular_rho_ww(T - 273.15)
    m_w_n       = rho_s.calcular_massa_agua(r, rho_ww_n)
    sigma_s_n   = sigma_s.calculate_sigma_s(T - 273.15, m_s_val, m_w_n)
    Phi_s_n     = Phi_s.calcular_phi_s(m_s_val, M_H2O, m_w_n)
    D_linha_w_n = D_linha_w.calculate_Dw_prime(T, R, P, r, M_H2O, T0, P0,
                                                 alpha_c=0.036, Delta_w=8e-8)
    L_v_n       = L_v.calcular_lv(T - 273.15)
    k_linha_a_n = k_linha_a.calculate_K_linha_a(T - 273.15, T, P, r, R, T0, P0,
                                                  M_a=28.9644e-3, alpha_T=0.7,
                                                  delta_T=2.16e-7, c_pa=1.006e3)
    rho_vr_n    = rho_vr.calcular_rho_vr(
        M_H2O,
        e_sat_esp.calcular_esat(T - 273.15),
        code_exp_y.calcular_exp_y(M_H2O, sigma_s_n, T_a_em_k, rho_w, R,
                                   v_ion, Phi_s_n, m_s_val, r, rho_s_val, M_NaCl),
        R, T_a_em_k)
    Y    = ((2 * M_H2O * sigma_s_n / (R * T_a_em_k * rho_w * r))
            - (v_ion * Phi_s_n * m_s_val * (M_H2O / M_NaCl) / (m_w_n - m_s_val)))
    den1 = rho_s_val * R * T_a_em_k / (D_linha_w_n * M_H2O * e_sat_val)
    den2 = (rho_s_val * L_v_n / (k_linha_a_n * T_a_em_k)
            * (L_v_n * M_H2O / (R * T_a_em_k) - 1))
    den  = den1 + den2
    return k_linha_a_n, D_linha_w_n, L_v_n, rho_vr_n, Y, den


def euler_fixo_sistema(f_system, y0, t):
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for n in range(len(t) - 1):
        dt = t[n + 1] - t[n]
        y[n + 1] = y[n] + dt * f_system(t[n], y[n])
    return y

def rk3_fixo_sistema(f_system, y0, t):
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for n in range(len(t) - 1):
        dt  = t[n + 1] - t[n]
        y_n = y[n];  t_n = t[n]
        F1  = f_system(t_n, y_n)
        y1  = y_n + dt * F1
        F2  = f_system(t_n + dt, y1)
        y2  = (3/4) * y_n + (1/4) * (y1 + dt * F2)
        F3  = f_system(t_n + dt / 2, y2)
        y[n + 1] = (1/3) * y_n + (2/3) * (y2 + dt * F3)
    return y


# subsistema rápido T (Euler com subpassos, r fixo)
def _f_rapido_compl(t, T, r_fixo):
    k_l, D_l, L_v_n, rho_vr_n, Y, den = _grandezas_dinamicas(r_fixo, T)
    return (3 / (rho_s_val * c_ps * r_fixo**2)) * (
        k_l * (T_a_em_k - T) + L_v_n * D_l * (rho_v_val - rho_vr_n))

# subsistema lento [r, m] (Euler macro, T já atualizado)
def _f_lento_compl(t, y_lento, T_final):
    r, m     = y_lento
    k_l, D_l, L_v_n, rho_vr_n, Y, den = _grandezas_dinamicas(r, T_final)
    dr_dt    = ((f - 1) - Y) / (r * den)
    H_n      = H_estrela.calcular_H_estrela(T_final, S)
    Dg_n     = Dg_estrela.calcular_Dg_estrela(r, T_a_em_k, R_atm)
    vol      = (4 / 3) * np.pi * r**3
    dm_dt    = 4 * np.pi * r * Dg_n * (C_ar - (m / vol) / (H_n * R_atm * T_final))
    return np.array([dr_dt, dm_dt])

def _passo_multiescala_e1(r_n, T_n, m_n, t_n, H, M):
    dt = H / M
    T  = T_n
    for j in range(M):
        T = T + dt * _f_rapido_compl(t_n + j * dt, T, r_n)
    T_new   = T
    y_l     = np.array([r_n, m_n])
    y_l_new = y_l + H * _f_lento_compl(t_n, y_l, T_new)
    return y_l_new[0], T_new, y_l_new[1]

def euler_multiescala_completo(r0, T0, m0, t_final, H, M):
    t_l, r_l, T_l, m_l = [0], [r0], [T0], [m0]
    t, r, T, m = 0, r0, T0, m0
    while t < t_final:
        H_step = min(H, t_final - t)
        if H_step < 1e-15: break
        r, T, m = _passo_multiescala_e1(r, T, m, t, H_step, M)
        t += H_step
        t_l.append(t); r_l.append(r); T_l.append(T); m_l.append(m)
    return np.array(t_l), np.array(r_l), np.array(T_l), np.array(m_l)

# ── E2: Euler PID adaptativo — sistema completo 
def euler_pid_completo(f_system, y0, t0, t_final,
                        dt_min, dt_max, dt_ini, tol, K_P, K_I, K_D):
    t_h = [t0];  y_h = [y0];  dt_h = [dt_ini]
    dt = dt_ini;  dt_prev = dt_min
    e_n_1 = tol;  e_n_2 = tol
    y_n = np.array(y0, dtype=float);  t_n = t0
    while t_n < t_final:
        if t_n + dt > t_final: dt = t_final - t_n
        F1     = f_system(t_n, y_n)
        y_next = y_n + dt * F1
        e_r = abs(y_next[0] - y_n[0]) / abs(y_next[0]) if y_next[0] != 0 else 1e-16
        e_T = abs(y_next[1] - y_n[1]) / abs(y_next[1]) if y_next[1] != 0 else 1e-16
        e_m = abs(y_next[2] - y_n[2]) / abs(y_next[2]) if y_next[2] != 0 else 1e-16
        e_n = max(e_r, e_T, e_m)
        if e_n > tol and dt > dt_min:
            fator = min(1 / e_n, 0.8)
            dt    = max(fator * dt, dt_min)
            dt_prev = (dt**2) / dt_prev
            continue
        t_n += dt;  y_n = y_next
        t_h.append(t_n);  y_h.append(y_n);  dt_h.append(dt)
        fator_P = (e_n_1 / e_n)              ** K_P
        fator_I = (tol   / e_n)              ** K_I
        fator_D = (e_n_1**2 / (e_n * e_n_2)) ** K_D
        dt      = np.clip(fator_P * fator_I * fator_D * dt, dt_min, dt_max)
        dt_prev = dt;  e_n_2 = e_n_1;  e_n_1 = e_n
    return np.array(t_h), np.array(y_h), np.array(dt_h)

# ── E3: Euler explícito fixo — sistema simplificado 
def euler_fixo_1d(f_m, m0, t):
    m = np.zeros(len(t))
    m[0] = m0
    for n in range(len(t) - 1):
        dt   = t[n + 1] - t[n]
        m[n + 1] = m[n] + dt * f_m(t[n], m[n])
    return m

# ── E4: Euler multiescala — sistema simplificado 
# r(t) e T(t) analíticos; apenas m integrado via Euler
def _f_lento_simp(t, m, T_new):
    r      = r_simplificado(t)
    H_t    = H_estrela.calcular_H_estrela(T_new, S)
    Dg_t   = Dg_estrela.calcular_Dg_estrela(r, T_a_em_k, R_atm)
    vol    = (4 / 3) * np.pi * r**3
    return 4 * np.pi * r * Dg_t * (C_ar - (m / vol) / (H_t * R_atm * T_new))

def _passo_multiescala_e4(r_n, T_n, m_n, t_n, H, M):
    T_new = T_simplificado(t_n + H)          # subsistema rápido: analítico
    m_new = m_n + H * _f_lento_simp(t_n, m_n, T_new)   # subsistema lento: Euler
    r_new = r_simplificado(t_n + H)
    return r_new, T_new, m_new

def euler_multiescala_simplificado(r0, T0, m0, t_final, H, M):
    t_l, r_l, T_l, m_l = [0], [r0], [T0], [m0]
    t, r, T, m = 0, r0, T0, m0
    while t < t_final:
        H_step = min(H, t_final - t)
        if H_step < 1e-15: break
        r, T, m = _passo_multiescala_e4(r, T, m, t, H_step, M)
        t += H_step
        t_l.append(t); r_l.append(r); T_l.append(T); m_l.append(m)
    return np.array(t_l), np.array(r_l), np.array(T_l), np.array(m_l)

# ── E5: Euler PID adaptativo — sistema simplificado 
def euler_pid_simplificado(f_m, m0, t0, t_final,
                             dt_min, dt_max, dt_ini, tol, K_P, K_I, K_D):
    t_h = [t0];  m_h = [m0];  dt_h = [dt_ini]
    dt = dt_ini;  dt_prev = dt_min
    e_n_1 = tol;  e_n_2 = tol
    m_n = float(m0);  t_n = t0
    while t_n < t_final:
        if t_n + dt > t_final: dt = t_final - t_n
        F1     = f_m(t_n, m_n)
        m_next = m_n + dt * F1
        e_n    = abs(m_next - m_n) / abs(m_next) if m_next != 0 else 1e-16
        if e_n > tol and dt > dt_min:
            fator = min(1 / e_n, 0.8)
            dt    = max(fator * dt, dt_min)
            dt_prev = (dt**2) / dt_prev
            continue
        if e_n == 0: e_n = 1e-16
        t_n += dt;  m_n = m_next
        t_h.append(t_n);  m_h.append(m_n);  dt_h.append(dt)
        fator_P = (e_n_1 / e_n)              ** K_P
        fator_I = (tol   / e_n)              ** K_I
        fator_D = (e_n_1**2 / (e_n * e_n_2)) ** K_D
        dt      = np.clip(fator_P * fator_I * fator_D * dt, dt_min, dt_max)
        dt_prev = dt;  e_n_2 = e_n_1;  e_n_1 = e_n
    return np.array(t_h), np.array(m_h), np.array(dt_h)


dt_fixo    = 1e-5
n_steps    = int(tau_f / dt_fixo)
tempo_fixo = np.linspace(0, tau_f, n_steps)

H_macro = 1e-4
M_sub   = 10

PID = dict(dt_min=1e-6, dt_max=1.0, dt_ini=1e-5, tol=1e-5,
           K_P=0.075, K_I=0.175, K_D=0.01)

# ECONTROL: Euler fixo, completo
sol_ec         = euler_fixo_sistema(f_r_T_m, y0, tempo_fixo)
raio_ec, temp_ec, massa_ec = sol_ec[:, 0], sol_ec[:, 1], sol_ec[:, 2]

# SCONTROL: SSPRK3 fixo, completo
sol_sc         = rk3_fixo_sistema(f_r_T_m, y0, tempo_fixo)
raio_sc, temp_sc, massa_sc = sol_sc[:, 0], sol_sc[:, 1], sol_sc[:, 2]

# E1: Euler subcycling, completo
t_e1, raio_e1, temp_e1, massa_e1 = euler_multiescala_completo(
    r_i, T_gota_em_k, m_i, tau_f, H_macro, M_sub)

# E2: Euler PID, completo
t_e2, sol_e2, dt_e2 = euler_pid_completo(f_r_T_m, y0, 0.0, tau_f, **PID)
raio_e2, temp_e2, massa_e2 = sol_e2[:, 0], sol_e2[:, 1], sol_e2[:, 2]

# E3: Euler fixo, simplificado
massa_e3 = euler_fixo_1d(f_m, m_i, tempo_fixo)
raio_e3  = r_simplificado(tempo_fixo)
temp_e3  = T_simplificado(tempo_fixo)

# E4: Euler subcycling, simplificado
t_e4, raio_e4, temp_e4, massa_e4 = euler_multiescala_simplificado(
    r_i, T_gota_em_k, m_i, tau_f, H_macro, M_sub)

# E5: Euler PID, simplificado
t_e5, massa_e5, dt_e5 = euler_pid_simplificado(f_m, m_i, 0.0, tau_f, **PID)
raio_e5 = r_simplificado(t_e5)
temp_e5 = T_simplificado(t_e5)

print("Pontos:")
print(f"  ECONTROL (Euler fixo, compl.): {len(tempo_fixo)}")
print(f"  SCONTROL (RK3  fixo, compl.): {len(tempo_fixo)}")
print(f"  E1  (Euler sub,  compl.):      {len(t_e1)}")
print(f"  E2  (Euler PID,  compl.):      {len(t_e2)}")
print(f"  E3  (Euler fixo, simp.):       {len(tempo_fixo)}")
print(f"  E4  (Euler sub,  simp.):       {len(t_e4)}")
print(f"  E5  (Euler PID,  simp.):       {len(t_e5)}")


plt.rcParams['text.usetex'] = False
plt.rcParams['font.size']   = 9

COR_COMPLETO = '#0D00FF'   # azul  — sistema completo
COR_SIMPLES  = '#FF0000'   # vermelho — sistema simplificado


fig1, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 7), sharex=True)

# Raio
ax1.plot(tempo_fixo, raio_ec * 1e6, '^-',  color=COR_COMPLETO, lw=2, ms=4,
         label=f'ECONTROL: r final = {raio_ec[-1]*1e6} µm')
ax1.plot(tempo_fixo, raio_sc * 1e6, 'D-',  color=COR_COMPLETO, lw=2, ms=4,
         label=f'SCONTROL: r final = {raio_sc[-1]*1e6} µm')
ax1.plot(t_e1,       raio_e1 * 1e6, 'o-',  color=COR_COMPLETO, lw=2, ms=3,
         label=f'E1: r final = {raio_e1[-1]*1e6} µm')
ax1.plot(tempo_fixo, raio_e3 * 1e6, '^--', color=COR_SIMPLES,  lw=2, ms=4, alpha=0.7,
         label=f'E3: r final = {raio_e3[-1]*1e6} µm')
ax1.plot(t_e4,       raio_e4 * 1e6, 'o--', color=COR_SIMPLES,  lw=2, ms=3, alpha=0.7,
         label=f'E4: r final = {raio_e4[-1]*1e6} µm')
ax1.set_ylabel('Raio (µm)')
ax1.set_xscale('log')
ax1.ticklabel_format(axis='y', useOffset=False)
ax1.legend(fontsize=7.5, loc='best', ncol=2)
ax1.grid(True, alpha=0.3, which='both')

# Temperatura
ax2.plot(tempo_fixo, temp_ec - 273.15, '^-',  color=COR_COMPLETO, lw=2, ms=4,
         label=f'ECONTROL: T final = {temp_ec[-1]} °C')
ax2.plot(tempo_fixo, temp_sc - 273.15, 'D-',  color=COR_COMPLETO, lw=2, ms=4,
         label=f'SCONTROL: T final = {temp_sc[-1]} °C')
ax2.plot(t_e1,       temp_e1 - 273.15, 'o-',  color=COR_COMPLETO, lw=2, ms=3,
         label=f'E1:  T final = {temp_e1[-1]} °C')
ax2.plot(tempo_fixo, temp_e3 - 273.15, '^--', color=COR_SIMPLES,  lw=2, ms=4, alpha=0.7,
         label=f'E3:  T final = {temp_e3[-1]} °C')
ax2.plot(t_e4,       temp_e4 - 273.15, 'o--', color=COR_SIMPLES,  lw=2, ms=3, alpha=0.7,
         label=f'E4:  T final = {temp_e4[-1]} °C')
ax2.set_ylabel('Temperatura (°C)')
ax2.set_xscale('log')
ax2.legend(fontsize=7.5, loc='best', ncol=2)
ax2.grid(True, alpha=0.3, which='both')

# Massa
ax3.plot(tempo_fixo, massa_ec, '^-',  color=COR_COMPLETO, lw=2, ms=4,
         label=f'ECONTROL: m final = {massa_ec[-1]} mol')
ax3.plot(tempo_fixo, massa_sc, 'D-',  color=COR_COMPLETO, lw=2, ms=4,
         label=f'SCONTROL: m final = {massa_sc[-1]} mol')
ax3.plot(t_e1,       massa_e1, 'o-',  color=COR_COMPLETO, lw=2, ms=3,
         label=f'E1: m final = {massa_e1[-1]} mol')
ax3.plot(tempo_fixo, massa_e3, '^--', color=COR_SIMPLES,  lw=2, ms=4, alpha=0.7,
         label=f'E3: m final = {massa_e3[-1]} mol')
ax3.plot(t_e4,       massa_e4, 'o--', color=COR_SIMPLES,  lw=2, ms=3, alpha=0.7,
         label=f'E4: m final = {massa_e4[-1]} mol')
ax3.set_ylabel('Massa (mol)')
ax3.set_xlabel('Tempo (s)')
ax3.set_xscale('log')
ax3.legend(fontsize=7.5, loc='best', ncol=2)
ax3.grid(True, alpha=0.3, which='both')

fig1.tight_layout()
fig1.savefig('all_fixo_sub_30.png', dpi=400)


fig2, (bx1, bx2, bx3, bx4) = plt.subplots(4, 1, figsize=(7, 9), sharex=True)

# Raio
bx1.plot(t_e2, raio_e2 * 1e6, 's-',  color=COR_COMPLETO, lw=2, ms=4,  label=f'E2: r final = {raio_e2[-1]*1e6} µm')
bx1.plot(t_e5, raio_e5 * 1e6, 's--', color=COR_SIMPLES,  lw=2, ms=4,
         alpha=0.7, label=f'E5: r final = {raio_e5[-1]*1e6} µm')
bx1.set_ylabel('Raio (µm)')
bx1.set_xscale('log')
bx1.ticklabel_format(axis='y', useOffset=False)
bx1.legend(fontsize=7.5, loc='best')
bx1.grid(True, alpha=0.3, which='both')

# Temperatura
bx2.plot(t_e2, temp_e2 - 273.15, 's-',  color=COR_COMPLETO, lw=2, ms=4,  label=f'E2: T final = {temp_e2[-1]} °C')
bx2.plot(t_e5, temp_e5 - 273.15, 's--', color=COR_SIMPLES,  lw=2, ms=4,
         alpha=0.7,  label=f'E5: T final = {temp_e5[-1]} °C')
bx2.set_ylabel('Temperatura (°C)')
bx2.set_xscale('log')
bx2.legend(fontsize=7.5, loc='best')
bx2.grid(True, alpha=0.3, which='both')

# Massa
bx3.plot(t_e2, massa_e2, 's-',  color=COR_COMPLETO, lw=2, ms=4,   label=f'E2: m final = {massa_e2[-1]} mol')
bx3.plot(t_e5, massa_e5, 's--', color=COR_SIMPLES,  lw=2, ms=4,
         alpha=0.7,  label=f'E5: m final = {massa_e5[-1]} mol')
bx3.set_ylabel('Massa (mol)')
bx3.set_xscale('log')
bx3.legend(fontsize=7.5, loc='best')
bx3.grid(True, alpha=0.3, which='both')

# Passo de tempo
bx4.plot(t_e2, dt_e2, 's-',  color=COR_COMPLETO, lw=2, ms=4,   label=f'E2: Δt final = {dt_e2[-1]:e} s')
bx4.plot(t_e5, dt_e5, 's--', color=COR_SIMPLES,  lw=2, ms=4,
         alpha=0.7, label=f'E5: Δt final = {dt_e5[-1]:e} s')
bx4.set_ylabel(r'$\Delta t$ (s)')
bx4.set_xlabel('Tempo (s)')
bx4.set_xscale('log')
bx4.set_yscale('log')
bx4.legend(fontsize=7.5, loc='best')
bx4.grid(True, alpha=0.3, which='both')

fig2.tight_layout()
fig2.savefig('all_pid_30.png', dpi=400)
#plt.show()
