#%%
import numpy as np
import src.vlm2d as vl2d
import matplotlib.pyplot as plt

from matplotlib import use

use('TkAgg')

#%% Far awar condition
V = 1
rho = 1
alpha = np.deg2rad(10)
Uinf = np.cos(alpha)*V
Winf = np.sin(alpha)*V
dynrho = 0.5*rho*V**2
 
#%%
cp, vp, n = vl2d.grid(10)

plt.figure()

plt.plot(cp[:,0], cp[:, 1], 'ko', vp[:, 0], vp[:, 1], 'ro')
plt.axis('equal')
plt.grid()
plt.show()

# %%
