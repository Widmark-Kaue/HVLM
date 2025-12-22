from src.vlm2 import np, plt, Wing, VLM
from src.utils import set_aiaa_style, plot_polar, Path

set_aiaa_style()
# General parameters
path = Path('data', 'avlCases')
savefig = True
# Simulation parameters
Vinf = 1
alpha = np.arange(0, 12, 1)
# Geometric parameters
span = 1
AR = 8 
taper = 0.6
sweep = 10
twist = 5
sections = ['Naca 0012', 'Naca 0012']
# Mesh parameters
nspan  = 5
nchord = 3

#%% Rectangular Wing
wing = Wing(span=span, AR=AR, sections=sections,name = 'rectangular_wing')
wing.mesh(nspan=nspan, nchord=nchord)
wing.plot_mesh(True)
wing_sim = VLM(wing)
polar = wing_sim.run_polar(Vinf=Vinf, alpha=alpha, ref='drela')
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)


#%% Wing with taper ratio
plt.close('all')
wing = Wing(span=span, AR = AR, sections=sections, taper=taper, name='taper_wing')
wing.mesh(nspan=nspan, nchord=nchord)
wing.plot_mesh(title=True)
polar = wing_sim.run_polar(Vinf=Vinf, alpha=alpha, ref='drela')
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)



#%% Wing with sweep
plt.close('all')
wing = Wing(span=span, AR = AR, sections=sections, sweep=sweep, name = 'sweep_wing')
wing.mesh(nspan=nspan, nchord=nchord)
wing.plot_mesh(title=True)
polar = wing_sim.run_polar(Vinf=Vinf, alpha=alpha)
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)


#%% Wing with twist
plt.close('all')
wing = Wing(span=span, AR = AR, sections=sections, twist=twist, name = 'twist_wing')
wing.mesh(nspan=nspan, nchord=nchord)
wing.plot_mesh(title=True)
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)

#%% Wing with camber
plt.close('all')
sections = ['Naca 2412', 'Naca 2412']
wing = Wing(span=span, AR=AR, sections=sections,name = 'camber_wing')
wing.mesh(nspan=nspan, nchord=nchord)
wing.plot_mesh(True)
wing_sim = VLM(wing)
polar = wing_sim.run_polar(Vinf=Vinf, alpha=alpha, ref='drela')
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)
pass