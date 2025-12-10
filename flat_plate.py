import numpy as np
import src.vlm as vlm
#%% Geometria da asa
span = 5
y = np.linspace(-span/2, span/2, 5)
x_le = np.zeros_like(y)
chord = np.ones_like(y)

geometry = np.column_stack([x_le, y, chord])

nspan  = 10
nchord = 5

mesh = vlm.meshPlanar(nspan, nchord, geometry)
vlm.plot_mesh(mesh)
#%% Cl xalpha
alphas = np.deg2rad(np.arange(0, 6))
Cl = np.zeros_like(alphas)
rho = 1.22
Vtot = 10
q = 0.5*rho*Vtot**2
S = span*np.mean(chord)
c = np.mean(chord)
for i, aoa in enumerate(alphas):
    Vinf = Vtot*np.array([np.cos(aoa), np.sin(aoa), 0])
    # Resolver sistemas
    A,B,RHS = vlm.influence_coefficients(Vinf, 500, mesh)
    Gamma = np.linalg.solve(A, RHS)
    
    #% Kutta-Joukwski theorem
    
    panels_span =mesh['panels_span'].reshape(Gamma.shape)
    print(f'{Vtot=}')
    Li = rho*Vtot*Gamma*panels_span
    L = np.sum(Li)
    
    Cl[i]  = L/(q*S*c)

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
