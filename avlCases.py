#%% imports
import numpy as np
import avlwrapper as avl
from pathlib import Path

#%% function
def exportpolar(results:dict, alpha:np.ndarray, path:Path):
    CL = np.zeros(len(alphas))
    CDi = np.zeros(len(alphas))
    for i in range(len(alphas)):
        CL[i] = results[i+1]['SurfaceForces']['wing']['CL']
        CDi[i] = results[i+1]['SurfaceForces']['wing']['CDi']
        
    exportData = np.column_stack([alphas, CL, CDi])
    np.savetxt(path, exportData, header='alpha\tCL\tCDi')

#%% General parameters
# Path to save data
path = Path('data', 'avlCases')
path.mkdir(exist_ok=True)

# Geometric parameters
b = 1           #span
AR = 8 
taper = 0.6
sweep = np.deg2rad(10)
twist = np.deg2rad(5)

S = b**2/AR
cmean = S/b
cr = 2*cmean/(1+taper)
ct = cr*taper

# Mesh parameters
nspan  = 5
nchord = 3

#%% Sweep case
base_case = avl.Case(name='polar',velocity = 1, density = 1.225)
alphas    =np.arange(0, 12)
all_cases = avl.create_sweep_cases(base_case=base_case,
                                   parameters=[{'name': 'alpha',
                                                'values': alphas}])

#%% Rectangular wing
root_section = avl.Section(
    leading_edge_point=avl.Point(0, 0, 0),
    chord = cmean,
    airfoil=avl.NacaAirfoil('0012')
    )
tip_section = avl.Section(
    leading_edge_point=avl.Point(0, b/2, 0),
    chord = cmean,
    airfoil=avl.NacaAirfoil('0012')
)

wing = avl.Surface(name='wing',
                   n_chordwise=nchord,
                   chord_spacing=avl.Spacing.equal,
                   n_spanwise=nspan,
                   span_spacing=avl.Spacing.equal,
                   y_duplicate=0.0,
                   sections=[root_section, tip_section])

WING = avl.Aircraft(name = 'rectangular_wing', 
             reference_area=S,
             reference_chord=cmean,
             reference_span=b,
             reference_point=avl.Point(0,0,0),
             surfaces=[wing])

session = avl.Session(geometry=WING, cases=all_cases)
session.show_geometry()
session.export_run_files(path = path / WING.name)
results = session.run_all_cases()
exportpolar(results, alphas, path / WING.name / 'polar.dat')

# %% Wing with taper ratio
root_section = avl.Section(
    leading_edge_point=avl.Point(0, 0, 0),
    chord = cr,
    airfoil=avl.NacaAirfoil('0012')
    )
tip_section = avl.Section(
    leading_edge_point=avl.Point(0, b/2, 0),
    chord = ct,
    airfoil=avl.NacaAirfoil('0012')
)

wing = avl.Surface(name='wing',
                   n_chordwise=nchord,
                   chord_spacing=avl.Spacing.equal,
                   n_spanwise=nspan,
                   span_spacing=avl.Spacing.equal,
                   y_duplicate=0.0,
                   sections=[root_section, tip_section])

WING = avl.Aircraft(name = 'taper_wing', 
             reference_area=S,
             reference_chord=cmean,
             reference_span=b,
             reference_point=avl.Point(0,0,0),
             surfaces=[wing])

session = avl.Session(geometry=WING, cases=all_cases)
session.show_geometry()
session.export_run_files(path = path / WING.name)
results = session.run_all_cases()
exportpolar(results, alphas, path / WING.name / 'polar.dat')
# %% Wing with sweep
root_section = avl.Section(
    leading_edge_point=avl.Point(0, 0, 0),
    chord = cmean,
    airfoil=avl.NacaAirfoil('0012')
    )
tip_section = avl.Section(
    leading_edge_point=avl.Point(np.tan(sweep)*b/2, b/2, 0),
    chord = cmean,
    airfoil=avl.NacaAirfoil('0012')
)

wing = avl.Surface(name='wing',
                   n_chordwise=nchord,
                   chord_spacing=avl.Spacing.equal,
                   n_spanwise=nspan,
                   span_spacing=avl.Spacing.equal,
                   y_duplicate=0.0,
                   sections=[root_section, tip_section])

WING = avl.Aircraft(name = 'sweep_wing', 
             reference_area=S,
             reference_chord=cmean,
             reference_span=b,
             reference_point=avl.Point(0,0,0),
             surfaces=[wing])

session = avl.Session(geometry=WING, cases=all_cases)
session.show_geometry()
session.export_run_files(path = path / WING.name)
results = session.run_all_cases()
exportpolar(results, alphas, path / WING.name / 'polar.dat')
# %% Wing with twist
root_section = avl.Section(
    leading_edge_point=avl.Point(0, 0, 0),
    chord = cmean,
    airfoil=avl.NacaAirfoil('0012')
    )
tip_section = avl.Section(
    leading_edge_point=avl.Point(0, b/2, np.tan(twist)*b/2),
    chord = cmean,
    airfoil=avl.NacaAirfoil('0012')
)

wing = avl.Surface(name='wing',
                   n_chordwise=nchord,
                   chord_spacing=avl.Spacing.equal,
                   n_spanwise=nspan,
                   span_spacing=avl.Spacing.equal,
                   y_duplicate=0.0,
                   sections=[root_section, tip_section])

WING = avl.Aircraft(name = 'twist_wing', 
             reference_area=S,
             reference_chord=cmean,
             reference_span=b,
             reference_point=avl.Point(0,0,0),
             surfaces=[wing])

session = avl.Session(geometry=WING, cases=all_cases)
session.show_geometry()
session.export_run_files(path = path / WING.name)
results = session.run_all_cases()
exportpolar(results, alphas, path / WING.name / 'polar.dat')