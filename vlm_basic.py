#%%
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

from matplotlib import use
use('TkAgg')

#%%
# ---------------------------
# 1) Parâmetros físicos e da asa
# ---------------------------
rho = 1           # densidade do ar (kg/m^3)
U_inf = 1            # velocidade de escoamento livre (m/s)
alpha_deg = 5.0         # ângulo de ataque (graus)
alpha = np.deg2rad(alpha_deg)

# Wing geometry
b = 10.0                # envergadura total (m)
c = 1.0                 # corda (m)
S = b * c               # área da asa (m^2) (retangular)

# Discretização (painéis)
N_span = 20             # número de divisões ao longo da envergadura (no semiplano: considere simetria?)
N_chord = 4             # número de divisões ao longo da corda

# Trailing leg "infinito" approximated length
L_trail = 100 * b     # comprimento do rastro finito (m) - grande para aproximar infinito

#%%
# ---------------------------
# 2) Geração da malha (pontos geométricos)
# ---------------------------
# coordenadas discretas (painéis) no domínio x (chordwise) e y (spanwise)
# x_local: de 0 (bordo de ataque) a c (bordo de fuga)
x_edges = np.linspace(0.0, c, N_chord + 1)       # bordas na direção da corda
y_edges = np.linspace(-b/2.0, b/2.0, N_span + 1) # bordas na direção da envergadura

# Calcule posições do bound vortex (1/4 c) e pontos de controle (3/4 c)
x_bound = 0.25 * c
x_control = 0.75 * c

# Arrays para armazenar posições (centros dos painéis em y e coordenadas do vórtice)
# Vamos indexar j = 0..N_span-1 (spanwise), i = 0..N_chord-1 (chordwise)
panel_centers_y = 0.5 * (y_edges[:-1] + y_edges[1:])   # centro spanwise (por painel)
panel_widths = y_edges[1:] - y_edges[:-1]              # spans locais (uniformes aqui)

# Número total de painéis
N = N_span * N_chord

# Criar arrays de posições para cada painel j (unrolled)
panel_y = np.repeat(panel_centers_y, N_chord)          # repetido por chordwise panels
panel_x_control = np.tile(np.full(N_chord, x_control), N_span)  # x do ponto de controle
panel_x_bound = np.tile(np.full(N_chord, x_bound), N_span)      # x do bound vortex

# coordenadas 3D dos pontos de controle e dos bound vortices
# plano da asa em z=0
ctrl_points = np.vstack([panel_x_control, panel_y, np.zeros_like(panel_y)]).T
bound_points = np.vstack([panel_x_bound, panel_y, np.zeros_like(panel_y)]).T

# Observação: o bound vortex será modelado por um segmento curto (pode-se colocar ele entre
# xi_left e xi_right se desejado). Aqui é modelado como segmento muito pequeno localizado
# em x_bound, com as pernas de rastro saindo do fim desse segmento (aqui tratamos o bound
# como um ponto onde o segmento está centrado em x_bound e as pernas saem do ponto 'end' = bound point).

plt.figure().add_subplot(projection = '3d')
plt.plot(ctrl_points[:, 0], ctrl_points[:, 1], ctrl_points[:, 2], 'ko')
plt.plot(bound_points[:, 0], bound_points[:, 1], bound_points[:, 2], 'ro')
plt.show()


#%%
# ---------------------------
# 3) Biot-Savart para segmento reto (influência de um segmento entre A e B no ponto P)
# ---------------------------
# POTENTIAL FUNCTION: Esta região pode ser extraída como função: induced_velocity_segment(A,B,P)
# Implementamos diretamente (vetorizado por P). Retorna vetor velocidade induzida por unidade de Gamma.
def induced_velocity_segment(A, B, P):
    """
    Induced velocity at point(s) P due to a vortex segment from A to B with circulation Gamma=1.
    A, B: arrays of length 3 (coordinates)
    P: array of shape (M,3) or (3,) (points where velocity is computed)
    returns v: array shape (M,3)
    Formula used (standard finite-segment Biot-Savart):
      v = (Gamma / (4*pi)) * (r1 x r2) / |r1 x r2|^2 * ( r0 . (r1/|r1| - r2/|r2|) )
    where r1 = P - A, r2 = P - B, r0 = B - A
    """
    P = np.atleast_2d(P)
    A = np.asarray(A)
    B = np.asarray(B)
    r1 = P - A
    r2 = P - B
    r0 = B - A
    # cross product r1 x r2 for each P
    cross = np.cross(r1, r2)
    cross_sq = np.sum(cross**2, axis=1)
    # norms of r1 and r2
    r1_norm = np.linalg.norm(r1, axis=1)
    r2_norm = np.linalg.norm(r2, axis=1)
    # To avoid division by zero for points very near vortex, add small eps
    eps = 1e-9
    r1_norm = np.maximum(r1_norm, eps)
    r2_norm = np.maximum(r2_norm, eps)
    denom = cross_sq + eps
    # dot product r0 . (r1/|r1| - r2/|r2|)
    factor = np.dot(r0, (r1.T / r1_norm - r2.T / r2_norm))
    # factor is scalar per point
    # broadcast to vector: v = (1/(4*pi)) * cross / denom * factor[:,None]
    coeff = factor / (4.0 * np.pi * denom)
    v = cross * coeff[:, None]
    return v  # per-unit-Gamma

# ---------------------------
# 4) Montagem da matriz de influência A (N x N) e RHS
# Cada coluna j corresponde a influência do vórtice de ferradura j com Gamma=1 nos N pontos de controle.
# Para cada vórtice j: temos
#   - um segmento bound curto (aqui considerado muito pequeno centralizado em x_bound) -> para modelar
#     rigidez do bound, podemos ignorar um segmento curto e modelar apenas as pernas + um "bound" (melhora se usar pequeno segmento)
#   - duas pernas de rastro (trailing legs): do extremo de bound até x + L_trail (em x) mantendo mesma y
# Implementação simples: representar bound como segmento very small along x (dx_bound), e as pernas: from bound end to (x+L_trail, y, z)
# ---------------------------
dx_bound = 1e-3  # comprimento muito pequeno do segmento bound para evitar singularidade (m)
# Para cada panel j: definimos pontos A (start) and B (end) do bound vortex (orientado em x)
# Vamos definir bound segment from (x_bound - dx/2) to (x_bound + dx/2) at same y
bound_A = np.vstack([np.full(N, x_bound - dx_bound/2.0), panel_y, np.zeros(N)]).T
bound_B = np.vstack([np.full(N, x_bound + dx_bound/2.0), panel_y, np.zeros(N)]).T

# Trailing legs: cada bound vortex tem duas pernas: elas vão para trás (+x) mantendo y (e z=0)
# representaremos cada perna como segmento de bound endpoint to endpoint + [L_trail, 0, 0]
# We will use only one trailing leg per side (but horseshoe has both sides at same y - here trailing goes to +inf)
# For a single horseshoe we have: bound segment (A->B) + trailing from B to B+L_trail*ex
# (this is a simplified horseshoe where trailing legs are parallel and go to +x).
trail_tip = np.array([L_trail, 0.0, 0.0])
trail_end = bound_B + trail_tip  # shape (N,3)

# Pre-allocate influence matrix
A = np.zeros((N, N))

# normal vector of panels: for thin wing it's z-direction
n_hat = np.array([0.0, 0.0, 1.0])

# Construir matriz A: coluna j = influência vertical (n_hat) no ponto i devido à ferradura j (por unidade Gamma)
for j in range(N):
    # vortex j: bound segment from bound_A[j] -> bound_B[j]
    A_seg = bound_A[j]
    B_seg = bound_B[j]
    # trailing segment (apenas uma perna rumando para +inf)
    T_end = trail_end[j]

    # Para a ferradura: velocidade induzida = devido ao bound segment + devido ao trailing segment
    # (no modelo clássico também se inclui a perna que vai ao infinito do outro lado; aqui simplificamos
    #  com um só trailing leg saindo do bound end).
    # Calculamos a influência unitária (Gamma=1) nos N pontos de controle
    v_bound = induced_velocity_segment(A_seg, B_seg, ctrl_points)    # (N,3)
    v_trail = induced_velocity_segment(B_seg, T_end, ctrl_points)   # (N,3)
    v_total = v_bound + v_trail

    # componente normal (z)
    A[:, j] = np.dot(v_total, n_hat)

# ---------------------------
# 5) RHS - condição de não-penetração (U + induced)·n = 0 -> A Γ = -U·n
# n = [0,0,1], U freestream com ângulo alpha: U = [U_inf*cos(alpha), 0, U_inf*sin(alpha)]
Uvec = np.array([U_inf * np.cos(alpha), 0.0, U_inf * np.sin(alpha)])
RHS = -np.dot(Uvec, n_hat) * np.ones(N)   # RHS é mesmo para todos os painéis (asa plana, sem twist)
# Observação: se desejar incluir variação por painel (ex: wing twist), a RHS pode variar.

# ---------------------------
# 6) Resolver sistema linear (A Γ = RHS)
# ---------------------------
# Resolver para vetores Γ (um por painel)
# Pode haver problemas numéricos se A for mal condicionado; aqui usamos solve direto
Gamma = linalg.solve(A, RHS)  # forma (N,)

# ---------------------------
# 7) Cálculo da sustentação e coeficiente de sustentação CL
# Lift L = rho * U∞ * ∫ Γ(y) dy  -> discretizado: sum_j rho * U∞ * Gamma_j * Δy_j
# CL = L / (0.5 * rho * U∞^2 * S)
# ---------------------------
# note: cada painel tem largura panel_widths repeated for chordwise panels
panel_spanwise_widths = np.repeat(panel_widths, N_chord)  # Δy por painel (m)
# total lift
L_total = rho * U_inf * np.sum(Gamma * panel_spanwise_widths)
CL = L_total / (0.5 * rho * U_inf**2 * S)

# ---------------------------
# 8) Cálculo do arrasto induzido via Trefftz-plane discrete:
# D_i = rho * ∫ Γ(y) * w(y) dy  -> discretizado: rho * sum_j Gamma_j * w_j * Δy_j
# onde w_j é downwash (velocidade na direção z) no ponto de controle j devido a todos os vórtices (incluindo o próprio).
# Observação: w_j = sum_k A_jk * Gamma_k  (porque A já é a componente normal por unidade Gamma)
# Assim D = rho * sum_j Gamma_j * w_j * Δy_j
# CDi = D / (0.5 * rho * U∞^2 * S)  -> note o denom U∞^2
# ---------------------------
# Calcule w (induced downwash) em cada controle (per unit sum)
w_induced = A.dot(Gamma)  # (N,) -> isso tem unidades de vel. (m/s)
D_induced = rho * np.sum(Gamma * w_induced * panel_spanwise_widths)
CDi = D_induced / (0.5 * rho * U_inf**2 * S)

# ---------------------------
# 9) Saídas e plots
# ---------------------------
print("------ Resultados ------")
print(f"Ángulo de ataque (deg): {alpha_deg:.3f}")
print(f"Velocidade U_inf (m/s): {U_inf}")
print(f"Envergadura b (m): {b}, corda c (m): {c}, área S (m^2): {S}")
print(f"N_panels (total): {N} (N_chord x N_span = {N_chord} x {N_span})")
print(f"Empuxo total (Lift) L = {L_total:.6f} N")
print(f"Coeficiente de sustentação CL = {CL:.6f}")
print(f"Arrasto induzido D_i = {D_induced:.6f} N")
print(f"Coeficiente de arrasto induzido CDi = {CDi:.6e}")

# Plot da distribuição de Gamma ao longo da envergadura (média chordwise)
# média de Gamma ao longo da corda para cada faixa spanwise
Gamma_mat = Gamma.reshape(N_span, N_chord)  # rows spanwise, cols chordwise
Gamma_span_avg = np.mean(Gamma_mat, axis=1)   # média em chordwise
y_span_centers = panel_centers_y

plt.figure(figsize=(8,5))
plt.plot(y_span_centers, Gamma_span_avg, marker='o')
plt.xlabel('y (m) - spanwise')
plt.ylabel('Circulation Γ (m^2/s)')
plt.title('Distribuição de circulação média por fatia spanwise')
plt.grid(True)
plt.show()

# Plot: downwash w(y)
w_mat = w_induced.reshape(N_span, N_chord)
w_span_avg = np.mean(w_mat, axis=1)
plt.figure(figsize=(8,5))
plt.plot(y_span_centers, w_span_avg, marker='o')
plt.xlabel('y (m) - spanwise')
plt.ylabel('Downwash w (m/s)')
plt.title('Downwash induzida média por fatia spanwise')
plt.grid(True)
plt.show()

# Plot: Gamma chordwise (opcional)
plt.figure(figsize=(8,5))
for idx in range(N_span):
    yv = panel_centers_y[idx]
    plt.plot(np.linspace(0, c, N_chord), Gamma_mat[idx,:], label=f"y={yv:.2f}")
plt.xlabel('x (m) - chordwise index')
plt.ylabel('Γ (m^2/s)')
plt.title('Γ por painel chordwise (cada linha = fatia spanwise)')
plt.grid(True)
plt.tight_layout()
plt.show()

# FIM DO SCRIPT

# ---------------------------
# REGIÕES QUE PODERIAM SER MODULARIZADAS/TRANSFORMADAS EM FUNÇÕES:
# 1) induced_velocity_segment(A,B,P) -> função (já separada aqui)
# 2) montagem da malha (geração de x_edges, y_edges, ctrl_points, bound_points)
# 3) montagem da matriz A (loop sobre vortices e cálculo de influências)
# 4) resolução linear e post-processamento (cálculo de L, D_i, CL, CDi)
# 5) rotinas de plot
# ---------------------------
