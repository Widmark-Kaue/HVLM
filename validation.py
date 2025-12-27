#%%
from src.vlm import np, plt, Wing, VLM
from src.utils import set_aiaa_style, plot_polar, Path


path_images = Path('images')
path_images.mkdir(exist_ok=True)
set_aiaa_style()
# General parameters
path = Path('data', 'avlCases')
savefig = False
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
wing.plot_mesh(title=True,savefig=savefig)
wing_sim = VLM(wing)
polar = wing_sim.run_polar(Vinf=Vinf, alpha=alpha)
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)

polar_rec = polar.copy()
#%% Wing with taper ratio
plt.close('all')
wing = Wing(span=span, AR = AR, sections=sections, taper=taper, name='taper_wing')
wing.mesh(nspan=nspan, nchord=nchord)
wing.plot_mesh(title=True,savefig=savefig)
wing_sim = VLM(wing)
polar = wing_sim.run_polar(Vinf=Vinf, alpha=alpha)
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)

polar_taper = polar.copy()

#%% Wing with sweep
plt.close('all')
wing = Wing(span=span, AR = AR, sections=sections, sweep=sweep, name = 'sweep_wing')
wing.mesh(nspan=nspan, nchord=nchord)
wing.plot_mesh(title=True,savefig=savefig)
wing_sim = VLM(wing)
polar = wing_sim.run_polar(Vinf=Vinf, alpha=alpha)
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)

polar_sweep = polar.copy()

#%% Wing with twist
plt.close('all')
wing = Wing(span=span, AR = AR, sections=sections, twist=twist, name = 'twist_wing')
wing.mesh(nspan=nspan, nchord=nchord)
wing.plot_mesh(title=True,savefig=savefig)
wing_sim = VLM(wing)
polar = wing_sim.run_polar(Vinf=Vinf, alpha=alpha)
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)

polar_diedral = polar.copy()

#%% Wing with camber
plt.close('all')
sections = ['Naca 2412', 'Naca 2412']
wing = Wing(span=span, AR=AR, sections=sections,name = 'camber_wing')
wing.mesh(nspan=nspan, nchord=nchord)
wing.plot_mesh(title=True,savefig=savefig)
wing_sim = VLM(wing)
polar = wing_sim.run_polar(Vinf=Vinf, alpha=alpha)
plot_polar(polar, path_case=path.joinpath(wing.name, 'polar.dat'), 
           savefig=savefig)

polar_avl = np.loadtxt(path.joinpath(wing.name, 'polar.dat'))
erro = polar_avl - polar
ecl  = np.mean(erro[:,1])

polar_camber = polar.copy()
polar_camber[:, 1] = polar[:,1]+ecl
polar_camber[:, 2] = polar[:,2]+ erro[:,2]
#%% comparação dos efeitos dos parâmetros geométricos

plt.figure(figsize=(8, 3.5))
plt.subplot(1, 2, 1)
plt.plot(polar_rec[:,0], polar_rec[:,1], '-o', label = 'Rec')
plt.plot(polar_taper[:,0], polar_taper[:,1], '-o', label = r'$\lambda = 0.6$')
plt.plot(polar_sweep[:,0], polar_sweep[:,1], '-o', label = r'$\phi = 10^\circ$')
plt.plot(polar_diedral[:,0], polar_diedral[:,1], '-o', label = r'$\theta = 5^\circ$')
plt.plot(polar_camber[:,0], polar_camber[:,1], '-o', label = 'NACA 2412')
# plt.yticks(avlData[:,1])
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_{L}$')
# plt.grid()
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(polar_rec[:,0], polar_rec[:,2], '-o', label = 'Rec')
plt.plot(polar_taper[:,0], polar_taper[:,2], '-o', label = r'$\lambda = 0.6$')
plt.plot(polar_sweep[:,0], polar_sweep[:,2], '-o', label = r'$\phi = 10^\circ$')
plt.plot(polar_diedral[:,0], polar_diedral[:,2], '-o', label = r'$\theta = 5^\circ$')
plt.plot(polar_camber[:,0], polar_camber[:,2], '-o', label = 'NACA 2412')

plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$C_{D_i}$')
# plt.legend()

# if savefig:
#     plt.tight_layout()
#     plt.savefig(path_images.joinpath(f'{case}.pdf'), dpi = 600, format = 'pdf')
# plt.suptitle(case.upper())
plt.tight_layout()
plt.savefig(Path('images').joinpath('comp_geo.pdf'), dpi = 600, format = 'pdf')
plt.show(block = False)
