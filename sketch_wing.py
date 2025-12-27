import numpy as np
from src.utils import set_aiaa_style, plt, Path
from matplotlib.patches import Arc

set_aiaa_style()
# ======================
# Parâmetros da asa
# ======================
b = 10.0        # envergadura total
cr = 2.0        # corda na raiz
ct = 1.0        # corda na ponta
phi = 25.0      # ângulo de enflechamento [graus] (medido com o eixo y)

phi_rad = np.deg2rad(phi)

# ======================
# Geometria
# ======================
y_root = 0.0
y_tip = b / 2

x_root_le = 0.0
x_tip_le = y_tip * np.tan(phi_rad)

# Bordos
x_le = [x_root_le, x_tip_le]
y_le = [y_root, y_tip]

x_te = [x_root_le + cr, x_tip_le + ct]
y_te = [y_root, y_tip]

# ======================
# Plot
# ======================
# fig, ax = plt.subplots(figsize=(8, 5))
fig, ax = plt.subplots(figsize = (6,4))

# ======================
# Limites manuais (inclui anotações)
# ======================
x_margin = 1.5
y_margin = 1.0

ax.set_xlim(-1.0, x_tip_le + ct + x_margin)
ax.set_ylim(-y_tip - y_margin, y_tip + y_margin)

# Asa
ax.plot(x_le, y_le, 'k')
ax.plot(x_te, y_te, 'k')
ax.plot([x_le[0], x_te[0]], [y_root, y_root], 'k')
ax.plot([x_le[1], x_te[1]], [y_tip, y_tip], 'k')

# Asa espelhada
ax.plot(x_le, [-y for y in y_le], 'k')
ax.plot(x_te, [-y for y in y_te], 'k')
ax.plot([x_le[1], x_te[1]], [-y_tip, -y_tip], 'k')

# ======================
# Eixos
# ======================
ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
# ax.set_xlabel('x')
# ax.set_ylabel('y')

# ======================
# Linha de medida da envergadura
# ======================
ax.annotate(
    '', xy=(x_tip_le + ct + 0.5, y_tip),
    xytext=(x_tip_le + ct + 0.5, -y_tip),
    arrowprops=dict(arrowstyle='<->')
)
ax.text(x_tip_le + ct + 0.7, 0.1, 'b', va='center')

# ======================
# Corda na raiz
# ======================
ax.annotate(
    '', xy=(x_root_le + cr, y_root - 0.4),
    xytext=(x_root_le, y_root - 0.4),
    arrowprops=dict(arrowstyle='<->')
)
ax.text(cr / 2, y_root - 0.6, r'$c_r$', ha='center')

# ======================
# Corda na ponta
# ======================
ax.annotate(
    '', xy=(x_tip_le + ct, y_tip + 0.4),
    xytext=(x_tip_le, y_tip + 0.4),
    arrowprops=dict(arrowstyle='<->')
)
ax.text(x_tip_le + ct / 2, y_tip + 0.6, r'$c_t$', ha='center')

# ======================
# Ângulo de enflechamento (arco)
# ======================
arc = Arc(
    (0, 0),
    width=2.0,
    height=2.0,
    angle=0,
    theta1=90 - phi,
    theta2=90,
    color='k'
)
ax.add_patch(arc)

ax.text(
    0.4 * np.sin(phi_rad),
    1.3 * np.cos(phi_rad),
    r'$\phi$'
)
# ======================
# Referência dos eixos coordenados (x, y)
# ======================
axis_len = 3  # comprimento das setas dos eixos

# eixo x
ax.annotate(
    '', xy=(axis_len, 0),
    xytext=(0, 0),
    arrowprops=dict(arrowstyle='->', linewidth=1.2)
)
ax.text(axis_len + 0.1, 0.1, 'x', va='center')

# eixo y
ax.annotate(
    '', xy=(0, axis_len),
    xytext=(0, 0),
    arrowprops=dict(arrowstyle='->', linewidth=1.2)
)
ax.text(-0.1, axis_len + 0.1, 'y', ha='center')
# ======================
# Ajustes finais
# ======================

# ax.set_aspect('equal', adjustable='box')
# ax.set_title('Vista em planta da asa (plano x–y)')
# ax.grid(True)
# ======================
# Aparência tipo sketch
# ======================
ax.set_axis_off()          # remove eixos, ticks e moldura
ax.set_aspect('equal')     # mantém proporção geométrica
plt.tight_layout()

plt.savefig(Path('images').joinpath('sketch_wing_xy_parameters.pdf'), dpi = 600, format = 'pdf')
plt.show(block = False)

# ======================
# Parâmetros da asa - plano y-z
# ======================
b = 10.0          # envergadura total
theta = 8.0       # ângulo de diedro [graus]

theta_rad = np.deg2rad(theta)

# ======================
# Geometria
# ======================
y_root = 0.0
y_tip = b / 2

z_root = 0.0
z_tip = y_tip * np.tan(theta_rad)

# ======================
# Plot
# ======================
fig, ax = plt.subplots(figsize=(6, 4))

# Asa (meia-asa)
ax.plot([y_root, y_tip], [z_root, z_tip], 'k', linewidth=1.5)

# Asa espelhada
ax.plot([-y_tip, y_root], [z_tip, z_root], 'k', linewidth=1.5)

# Linha de referência (asa sem diedro)
ax.plot([-y_tip, y_tip], [0, 0], '--', color = 'gray', linewidth=0.8)

# ======================
# Ângulo de diedro (arco)
# ======================
arc = Arc(
    (0, 0),
    width=2.0,
    height=2.0,
    angle=0,
    theta1=0,
    theta2=theta,
    color='k'
)
ax.add_patch(arc)

ax.text(
    1.5 * np.cos(theta_rad / 2),
    0.4 * np.sin(theta_rad / 2),
    r'$\theta$'
)

# ======================
# Referência dos eixos coordenados
# ======================
axis_len = 1

# eixo y
ax.annotate(
    '', xy=(axis_len, 0),
    xytext=(0, 0),
    arrowprops=dict(arrowstyle='->', linewidth=1.2)
)
ax.text(axis_len + 0.1, -0.1, 'y', va='center')

# eixo z
ax.annotate(
    '', xy=(0, axis_len),
    xytext=(0, 0),
    arrowprops=dict(arrowstyle='->', linewidth=1.2)
)
ax.text(0, axis_len + 0.1, 'z', ha='center')

# ======================
# Limites manuais
# ======================
margin = 1.0
ax.set_xlim(-y_tip - margin, y_tip + margin)
ax.set_ylim(-0.5, z_tip + margin)

# ======================
# Aparência tipo sketch
# ======================
ax.set_aspect('equal')
ax.set_axis_off()
plt.tight_layout()
plt.savefig(Path('images').joinpath('sketch_wing_yz_parameters.pdf'), dpi = 600, format = 'pdf')

plt.show(block = False)

plt.close('all')
# ======================
# Parâmetros NACA
# ======================
m = 0.04     # arqueamento máximo (ex: NACA 2412 → m = 0.02)
p = 0.4      # posição do arqueamento máximo
n = 400      # resolução

# ======================
# Malha
# ======================
x = np.linspace(0, 1, n)
zc = np.zeros_like(x)

# ======================
# Linha média
# ======================
mask1 = x < p
mask2 = x >= p

zc[mask1] = (m / p**2) * (2*p*x[mask1] - x[mask1]**2)
zc[mask2] = (m / (1 - p)**2) * (1 - 2*p + 2*p*x[mask2] - x[mask2]**2)

# Valor máximo
zmax = m
xmax = p

# ======================
# Plot
# ======================
fig, ax = plt.subplots(figsize=(6, 4))

# Linha média
ax.plot(x, zc, 'k', linewidth=1.8)

# Linha de referência (corda)
ax.plot([0, 1], [0, 0], 'k--', linewidth=0.8)

# ======================
# Evidenciar p
# ======================
ax.plot([p, p], [0, zmax], 'k:', linewidth=1.2)
ax.plot(p, zmax, 'ko', markersize=4)
ax.text(p, -0.005, r'$p$', ha='center', va='top')

# ======================
# Evidenciar m
# ======================
ax.annotate(
    '', xy=(p + 0.06, zmax),
    xytext=(p + 0.06, 0),
    arrowprops=dict(arrowstyle='<->', linewidth=1.2)
)
ax.text(p + 0.08, zmax / 2, r'$m$', va='center')

# ======================
# Referência dos eixos
# ======================
axis_len = 0.12

# eixo x_c
ax.annotate(
    '', xy=(axis_len, 0),
    xytext=(0, 0),
    arrowprops=dict(arrowstyle='->', linewidth=1.2)
)
ax.text(axis_len + 0.01, -0.005, r'$x_c$', va='center')

# eixo z_c
ax.annotate(
    '', xy=(0, axis_len- 0.05),
    xytext=(0, 0),
    arrowprops=dict(arrowstyle='->', linewidth=1.2, color = 'k')
)
ax.text(0, axis_len + 0.005-0.05, r'$z_c$', ha='center')

# ======================
# Limites e aparência
# ======================
ax.set_xlim(-0.08, 1.05)
ax.set_ylim(-0.05, 2 * m)
ax.set_aspect('auto')
ax.set_axis_off()
ax.set_xticks([0,1])

plt.tight_layout()
plt.savefig(Path('images').joinpath('sketch_cambernaca4digits.pdf'), dpi = 600, format = 'pdf')
plt.show()
