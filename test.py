#%% Imports
import numpy as np
import matplotlib.pyplot as plt
import src.vlm as vlm
from matplotlib import use
use('TkAgg')
#%% Geometria da asa
span = 5
y = np.linspace(-span/2, span/2, 5)
x_le = np.zeros_like(y)
chord = np.ones_like(y)

geometry = np.column_stack([x_le, y, chord])

nspan  = 20
nchord = 10

X, Y, Z, panels, panels_id, normals, area = vlm.meshPlanar(nspan, nchord, geometry)

vcp = vlm.vortexAndControlPoint(panels)
vlm.plot_mesh(X, Y, Z, vcp, normals)


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
