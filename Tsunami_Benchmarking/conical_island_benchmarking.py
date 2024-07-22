# Laboratry experiments exploring the runup of tsunami waves on a conical island have been used to benchmark many tsunami models
# with different formulations. In the following code, the laboratory experiment is replicated, producing elevations at different
# points on the island as a result of an incoming soitary wave.

import scipy.interpolate
import math
from thetis import *
from firedrake import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="firedrake.interpolation")

# Load predefined mesh
mesh2d = Mesh('req_data/conical_mesh_2.msh')
P1_2d = FunctionSpace(mesh2d, 'CG', 1)

# Input bathymetry parameters
height = 0.625
toe_diameter = 7.2
crest_diameter = 2.2
toe_radius = toe_diameter / 2.0
crest_radius = crest_diameter / 2.0
center_x, center_y = 15, 12.5

# Bathymetry assigned according to the parameters detailed above
bathymetry_2d = Function(P1_2d, name='Bathymetry')
x, y = SpatialCoordinate(mesh2d)
r = sqrt((x - center_x)**2 + (y - center_y)**2)
bathymetry_expr = conditional(
    le(r, crest_radius),
    height,
    conditional(
        le(r, toe_radius),
        height - (height / (toe_radius - crest_radius)) * (r - crest_radius),
        0.0
    )
)
bathymetry_2d.interpolate(bathymetry_expr-0.32)


# Set initial elevation to 0
elev_init_expr = 0
elev_init = Function(P1_2d)
elev_init.interpolate(elev_init_expr)

# Set initial velocity to negligibel non-zero
P1v_2d = VectorFunctionSpace(mesh2d, 'CG', 1)
vel_init = Function(P1v_2d)
velocity_expr = as_vector((0.000015,0.0))
vel_init.interpolate(velocity_expr)

# Create Solver - include nonhydrostatic equations
solve_nonhydrostatic_pressure = True
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
options = solver_obj.options
options.simulation_export_time =0.1
options.simulation_end_time = 50
options.swe_timestepper_type = 'CrankNicolson'
options.timestep = 0.01
options.output_directory = 'conical_island_outputs_nh_0.32_f0.0015_v0.001'
options.fields_to_export = ['elev_2d', 'uv_2d']
options.fields_to_export_hdf5 = ['elev_2d', 'uv_2d']
options.use_wetting_and_drying = True
options.wetting_and_drying_alpha = Constant(0.32)
options.horizontal_viscosity = Constant(0.001)
options.linear_drag_coefficient = Constant(0.0015)
if solve_nonhydrostatic_pressure:
    options_nh = options.nh_model_options
    options_nh.solve_nonhydrostatic_pressure = solve_nonhydrostatic_pressure


# Define Solitary wave according to Bousinessq eq.
def sech(x):
    return 1/(np.cosh(x))

# Define wave parameters
h = 0.32
H = 0.045*h
k = np.sqrt(3*H/(4*h**3))
c = np.sqrt(9.81*h)

def solitary_wave(simulation_time):
    """Time-dependent tidal elevation"""
    p = (H / (h * k)) * ((h * np.tanh(k * c * simulation_time)) / (h + H * (1 - (np.tanh(k * c * simulation_time))**2)))+6.25
    eta = H*(sech(k*(c*simulation_time-p)))**2
    return eta

# Assign Boundary conditions
left_bnd_id = 1
right_bnd_id = 2
bot_bnd_id = 3
top_bnd_id = 4
wave_bnd_id = 0
swe_bnd = {}

solitary_wave_const = Constant(solitary_wave(0))
swe_bnd[wave_bnd_id] = {'elev': solitary_wave_const}
solver_obj.bnd_functions['shallow_water'] = swe_bnd

def update_forcings(t_new):
    """Callback function that updates all time dependent forcing fields"""
    uv, elev = solver_obj.fields.solution_2d.subfunctions
    solitary_wave_const.assign(solitary_wave(t_new))

# Assign initial conditions
solver_obj.assign_initial_conditions(elev=elev_init) #, uv=vel_init)

# Create Observation Stations
stations = [('stationA', (15, 8.9)),
            ('stationB', (15, 11.4)),
            ('stationC', (16.1, 12.5)),
            ('stationD', (15, 13.6)),
            ('stationE', (15, 0)),
            ]
for name, (sta_x, sta_y) in stations:
    cb = TimeSeriesCallback2D(solver_obj, ['elev_2d'], sta_x, sta_y, name)
    solver_obj.add_callback(cb)

# Begin Solving
#solver_obj.iterate(update_forcings=update_forcings)