import numpy as np
import src.vlm as vlm
import matplotlib.pyplot as plt

"""
Atividades VLM com vórtice de ferradura: 
    - Alihar esteira com velocidade freestream
    - Conferir span multiplicando sustentação
    - Calcular arrasto
    - Condição de simetria
    - Testar diedro
    - Enflechamento
    - Arqueamento
    
Atividade VLM com Kawada-Hardin:
    - Implementar velocidade induzida pela esteira helicoidal
    - Adaptar condições de operação para rotor
"""


#%% Geometria da asa
span = 5
y = np.linspace(-span/2, span/2, 5)
x_le = np.zeros_like(y)
chord = np.ones_like(y)

geometry = np.column_stack([x_le, y, chord])

nspan  = 20
nchord = 5

mesh = vlm.meshPlanar(nspan, nchord, geometry)
vlm.plot_mesh(mesh)
#%% Cl xalpha
alphas_deg = np.arange(0, 6)
alphas = np.deg2rad(alphas_deg)

# alphas = np.deg2rad([3])

Cl = np.zeros_like(alphas)
rho = 1.22
Vtot = 10
q = 0.5*rho*Vtot**2
S = span*np.mean(chord)
c = np.mean(chord)
for i, aoa in enumerate(alphas):
    Vinf = Vtot*np.array([np.cos(aoa), 0, np.sin(aoa)])
    # Resolver sistemas
    A,B,RHS = vlm.influence_coefficients(Vinf, 500, mesh)
    Gamma = np.linalg.solve(A, RHS)
    
    #% Kutta-Joukwski theorem
    
    panels_span =mesh['panels_span'].reshape(Gamma.shape)
    Li = rho*Vtot*Gamma*panels_span
    L = np.sum(Li)
    
    Cl[i]  = L/(q*S*c)
    print(f'alpha = {alphas_deg[i]:.1f}: Cl={Cl[i]:.5f}')

#%% dados xflr5
xf = np.loadtxt('data/retangular-10_0 m_s-case1.txt', skiprows=8, usecols=(0, 2, 3))

plt.figure()
plt.subplot(1,2,1)
plt.plot(xf[:, 0], xf[:, 1], 'bo-')
plt.plot(alphas_deg, Cl, 'ks-')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_L$')
plt.grid()

plt.subplot(1,2,2)
plt.plot(xf[:, 0], xf[:, 2], 'ro-')
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_{Di}$')
plt.grid()

plt.tight_layout()
plt.show()
# %% asa enflechada
# span = 4.0
# c_root = 1
# taper_ratio = 0.1
# N = 20

# y = np.linspace(0, span/2, N)
# chord = c_root * (1 - (1 - taper_ratio) * y/span)
# sweep_deg = 30
# sweep = np.tan(np.radians(sweep_deg))

# x_le = sweep * y


# geometry = np.column_stack([x_le, y, chord])

# X, Y, Z, panels = vlm.meshPlanar(10, 10, geometry)

# vlm.plot_mesh(X, Y, Z)


# %%
