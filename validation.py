#%%
import numpy as np
import src.vlm as vlm
import matplotlib.pyplot as plt

"""
Atividades VLM com vórtice de ferradura: 
    - Alihar esteira com velocidade freestream
    - Conferir span multiplicando sustentação
    - Calcular arrasto
    - Condição de simetria [check]
    - Testar diedro
    - Enflechamento
    - Arqueamento
    
Atividade VLM com Kawada-Hardin:
    - Implementar velocidade induzida pela esteira helicoidal
    - Adaptar condições de operação para rotor
"""

#%%
# =============================================================================
# Validation Cases - VLM with horseshoe elements
# =============================================================================
# Geometric parameters
span = 5
AR = 8 

# Mesh parameters
nspan  = 5
nchord = 3

# Rectangular wing
geo_rec = vlm.geometry(span, AR)
mesh_rec = vlm.meshPlanar(nspan, nchord, geo_rec)
vlm.plot_mesh(mesh_rec, 'Rectangular Wing')

# Wing with taper ratio
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
#%% run polars

Vinf = 1
rho = 1.225
alpha = np.arange(0, 7, 1)
# Rectangular wing
CL, CDi = vlm.run_polar(Vinf=Vinf, rho=rho,alpha=alpha, mesh=mesh_rec)
# %%
    