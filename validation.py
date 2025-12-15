#%%
import numpy as np
import src.vlm as vlm
from src.utils import set_aiaa_style
import matplotlib.pyplot as plt
from pathlib import Path

from matplotlib import use
# use('TkAgg')
set_aiaa_style()
"""
Atividades VLM com vórtice de ferradura: 
    - Conferir span multiplicando sustentação
    - Calcular arrasto [parcial check]
    - Condição de simetria [check]
    - Testar diedro
    - Enflechamento
    - Arqueamento
    
Atividade VLM com Kawada-Hardin:
    - Implementar velocidade induzida pela esteira helicoidal
    - Adaptar condições de operação para rotor
"""
#%%

#%%
# =============================================================================
# Validation Cases - VLM with horseshoe elements
# =============================================================================
# General parameters
Vinf = 1
rho = 1.225
alpha = np.arange(0, 12, 1)
path = Path('data', 'avlCases')

# Geometric parameters
span = 1
AR = 8 

# Mesh parameters
nspan  = 5
nchord = 3

#%% Rectangular wing
case = 'rectangular_wing'
geo_rec = vlm.geometry(span, AR)
mesh_rec = vlm.meshPlanar(nspan, nchord, geo_rec)
vlm.plot_mesh(mesh_rec, 'Rectangular Wing')
CL, CDi, CDi2 = vlm.run_polar(Vinf=Vinf, rho=rho,alpha=alpha, mesh=mesh_rec)
avlData = np.loadtxt(path/case/'polar.dat')

plt.figure(figsize=(8, 3.5))
plt.subplot(1, 2, 1)
plt.plot(avlData[:,0], avlData[:,1], 'k-o', label = 'AVL')
plt.plot(alpha, CL, 'r--s', label = 'code')
# plt.yticks(avlData[:,1])
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_{L}$')
# plt.grid()
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(avlData[:,0], avlData[:,2], 'k-o', label = 'AVL')
plt.plot(alpha, CDi, 'r--s', label = 'code')
plt.plot(alpha, CDi2, 'm--^', label = 'code 2')
# plt.yticks(avlData[:,1])
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_{D_i}$')
# plt.legend()

plt.tight_layout()
plt.show(block = False)



#%% Wing with taper ratio
geo_taper = vlm.geometry(span, AR, taper=0.6)
mesh_taper = vlm.meshPlanar(nspan, nchord, geo_taper)
vlm.plot_mesh(mesh_taper, r'Wing with taper ratio: $\lambda = 0.6$')

# Wing with sweep
geo_sweep = vlm.geometry(span, AR, sweep=10)
mesh_sweep= vlm.meshPlanar(nspan, nchord, geo_sweep)
vlm.plot_mesh(mesh_sweep, r'Wing with sweep: $\phi = 10^\circ$')

# Wing with twist
geo_twist = vlm.geometry(span, AR, twist=5)
mesh_twist = vlm.meshPlanar(nspan, nchord, geo_twist)
vlm.plot_mesh(mesh_twist, r'Wing with twist: $\theta = 5^\circ$')

plt.close('all')


    