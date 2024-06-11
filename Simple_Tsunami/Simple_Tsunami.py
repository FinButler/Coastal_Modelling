# Import necessary modules
#import os
import scipy.interpolate
import math
from thetis import *
import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="firedrake.interpolation")
warnings.filterwarnings("ignore", category=FutureWarning, module="pyproj.transform")

# Create mesh and define mesh elements
lx = 10e3
ly = 10e3
nx = 20
ny = 20
mesh2d = RectangleMesh(nx, ny, lx, ly)
P1_2d = FunctionSpace(mesh2d, 'CG', 1)

# Define baythymetry
bathymetry_2d = Function(P1_2d, name='Bathymetry')
x = SpatialCoordinate(mesh2d)
bathymetry_expr = 20.0 * (x[0] / lx)
bathymetry_2d.interpolate(bathymetry_expr)

# Define tsunami Wave by inputting a gaussian wave as an initial BC
elev_init = Function(P1_2d)
wave_height = 5.0  # Height of the wave in meters
elev_expr = wave_height * exp(-((x[0] - lx + 50) ** 2) / (2 * 100 ** 2))  # Gaussian wave
elev_init.interpolate(elev_expr)

# Define initial velocity - necessary to include friction
P1v_2d = VectorFunctionSpace(mesh2d, 'CG', 1)
vel_init = Function(P1v_2d)
velocity_expr = as_vector((0.00001,0.0))
vel_init.interpolate(velocity_expr)

# Create solver
solver_obj = solver2d.FlowSolver2d(mesh2d, bathymetry_2d)
options = solver_obj.options
options.simulation_export_time = 360
options.simulation_end_time = 3600
options.swe_timestepper_type = 'CrankNicolson'
options.timestep = 5

# Choose output directory and desired outputs
options.output_directory = 'FB_tsunami_outputs'
options.fields_to_export = ['elev_2d', 'uv_2d']

# Inclusion of wetting and drying as well as friction parameter
options.use_wetting_and_drying = True
options.manning_drag_coefficient = Constant(0.03)

# Displaying Boundary IDs
left_bnd_id = 1
right_bnd_id = 2
bot_bnd_id = 3
top_bnd_id = 4

# Defining boundary conditions and assigning to solver
swe_bnd = {}
swe_bnd[right_bnd_id] = {'elev': Constant(0.0),}
swe_bnd[top_bnd_id] = {'elev': Constant(0.0)}
swe_bnd[bot_bnd_id] = {'elev': Constant(0.0),}
solver_obj.bnd_functions['shallow_water'] = swe_bnd

# Assigning initial conditions to solver and running model
solver_obj.assign_initial_conditions(elev=elev_init, uv=vel_init)
solver_obj.iterate()
