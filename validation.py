#%%
import numpy as np
import src.vlm as vlm
from src.utils import set_aiaa_style, plot_polar
import matplotlib.pyplot as plt
from pathlib import Path

from matplotlib import use
# use('TkAgg')
set_aiaa_style()
"""
Atividades VLM com vórtice de ferradura: 
    - Conferir span multiplicando sustentação [parcial check]
    - Calcular arrasto [parcial check]
    - Condição de simetria [check]
    - Testar diedro [parcial check]
    - Enflechamento [parcial check]
    - Arqueamento
    
Atividade VLM com Kawada-Hardin:
    - Implementar velocidade induzida pela esteira helicoidal
    - Adaptar condições de operação para rotor
"""
#%% functions


#%%
# =============================================================================
# Validation Cases - VLM with horseshoe elements
# =============================================================================
# General parameters
Vinf = 1
rho = 1.225
alpha = np.arange(0, 12, 1)
path = Path('data', 'avlCases')
savefig = False
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
vlmData = vlm.run_polar(Vinf=Vinf, alpha=alpha, mesh=mesh_rec)
plot_polar(vlmData=vlmData, path_case=path.joinpath(case, 'polar.dat'), savefig=savefig)

#%% Wing with taper ratio
plt.close('all')
case = 'taper_wing'
geo_taper = vlm.geometry(span, AR, taper=0.6)
mesh_taper = vlm.meshPlanar(nspan, nchord, geo_taper)
vlm.plot_mesh(mesh_taper, r'Wing with taper ratio: $\lambda = 0.6$')
vlmData = vlm.run_polar(Vinf=Vinf, alpha=alpha, mesh =mesh_taper)
plot_polar(vlmData=vlmData, path_case=path.joinpath(case, 'polar.dat'), savefig=savefig)

#%% Wing with sweep
plt.close('all')
case = 'sweep_wing'
geo_sweep = vlm.geometry(span, AR, sweep=10)
mesh_sweep= vlm.meshPlanar(nspan, nchord, geo_sweep)
vlm.plot_mesh(mesh_sweep, r'Wing with sweep: $\phi = 10^\circ$')
vlmData = vlm.run_polar(Vinf=Vinf, alpha=alpha, mesh=mesh_sweep)
plot_polar(vlmData=vlmData, path_case=path.joinpath(case, 'polar.dat'), savefig=savefig)


#%% Wing with twist
plt.close('all')
case = 'twist_wing'
geo_twist = vlm.geometry(span, AR, twist=5)
mesh_twist = vlm.meshPlanar(nspan, nchord, geo_twist)
vlm.plot_mesh(mesh_twist, r'Wing with twist: $\theta = 5^\circ$')
vlmData = vlm.run_polar(Vinf=Vinf, alpha=alpha, mesh=mesh_twist)
plot_polar(vlmData=vlmData, path_case=path.joinpath(case, 'polar.dat'), savefig=savefig)
pass

    