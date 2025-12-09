#%% Imports
import numpy as np
import src.vlm as vlm
#%% Geometria da asa
span = 5
y = np.linspace(-span/2, span/2, 5)
x_le = np.zeros_like(y)
chord = np.ones_like(y)

geometry = np.column_stack([x_le, y, chord])

nspan  = 20
nchord = 10

mesh = vlm.meshPlanar(nspan, nchord, geometry)
vlm.plot_mesh(mesh)
aoa = np.deg2rad(5)
# Vinf = np.array([5, 0, 0])
Vinf = np.array([5*np.cos(aoa), 5*np.sin(aoa), 0])
#%% Resolver sistemas
A,B,RHS = vlm.influence_coefficients(Vinf, 500, mesh)
Gamma = np.linalg.solve(A, RHS)

#%% Kutta-Joukwski theorem
rho = 1.22
Vtot = np.linalg.norm(Vinf)
panels_span =mesh['panels_span'].reshape(Gamma.shape)
print(f'{Vtot=}')
Li = rho*Vtot*Gamma*panels_span
L = np.sum(Li)

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
