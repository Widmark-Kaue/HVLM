#%% imports
import numpy as np
import avlwrapper as avl

#%% General parameters
# Geometric parameters
b = 5           #span
AR = 8 
taper = 0.6
sweep = 10
twist = 5

S = b**2/AR
cmean = S/b

# Mesh parameters
nspan  = 5
nchord = 3

#%% Create retangular wing
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

session = avl.Session(geometry=wing)

# # check if we have ghostscript
# if 'gs_bin' in session.config.settings:
#     img = session.save_geometry_plot()[0]
#     avl.show_image(img)
# else:
#     session.show_geometry()
