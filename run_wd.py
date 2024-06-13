#  simulation
# ================================================================================
from thetis import *
#Supressing irelevant warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="firedrake.interpolation")
warnings.filterwarnings("ignore", category=FutureWarning, module="firedrake.split()")
# things to alter (relative path to this script)
mesh_file = 'mesh/MonaiValley_B.msh'
bathymetry = 'raw_data/Bathymetry.grd'
# so make sure it matches your boundary forcing zone (i.e. it might be best to be in zone 30
# if modelling the shetlands, even though the centre might be zone 31, dependig on the size of your domain)
# might need to edit this
viscosity = 0.0001
physID = 10
dets = [[4.521, 1.196],
         [4.521, 1.696],
         [4.521, 2.196]]


output_directory = "outputs"
t_end = 50.0
use_wetting_and_drying = True
bottom_friction = 0.03
timestep = 0.05
alpha = 0.01

def boundary_forcing(elev,t):
    import inputwave 
    xvector = mesh2d.coordinates.dat.data
    evector = elev.dat.data
    for i,xy in enumerate(xvector):
        evector[i] = inputwave.val(t)
    
    return evector

def get_bathymetry_from_netcdf(bathymetry_file, mesh2d, minimum_depth=None, abscissa_name='lon', ordinate_name='lat', variable_name='z', negative_netcdf_depth=True):
    from netCDF4 import Dataset as NetCDFFile
    import scipy.interpolate
    #Read data from NetCDF file.
    print_output('Reading bathymetry/topography data from '+bathymetry_file)
    nc = NetCDFFile(bathymetry_file)
    lat = nc.variables[ordinate_name][:]
    lon = nc.variables[abscissa_name][:]
    values = nc.variables[variable_name][:,:]
    #Construct a bathymetry function in the appropriate function space.
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)
    bathymetry_2d = Function(P1_2d, name="bathymetry")
    xvector = mesh2d.coordinates.dat.data
    bvector = bathymetry_2d.dat.data
    assert xvector.shape[0]==bvector.shape[0]
    print_output('Interpolating bathymetry/topography data onto simulation mesh')
    interpolator = scipy.interpolate.RegularGridInterpolator((lat, lon), values)
    if negative_netcdf_depth:
        netcfd_bathy_sign = -1
    else:
        netcfd_bathy_sign = 1
    if minimum_depth != None:
        for i,(y,x) in enumerate(zip(xvector[:][:,1], xvector[:][:,0])):
            bvector[i] = max(netcfd_bathy_sign*interpolator((y, x)), minimum_depth)
    else:
        for i,(y,x) in enumerate(zip(xvector[:][:,1], xvector[:][:,0])):
            bvector[i] = netcfd_bathy_sign*interpolator((y, x))
    print_output('Bathymetry/topography min:'+str(min(bvector)))
    print_output('Bathymetry/topography max:'+str(max(bvector)))
    return bathymetry_2d

   

with timed_stage('reading mesh'):
  mesh2d = Mesh(mesh_file)
output_dir = create_directory(output_directory)

print_output('Loaded mesh '+mesh2d.name)
print_output('Exporting to '+output_dir)
    
t_start = 0
# export interval in seconds
t_export = 0.5

# bathymetry
P1 = FunctionSpace(mesh2d, "CG", 1)
bathymetry2d = get_bathymetry_from_netcdf(bathymetry, mesh2d, abscissa_name='x', ordinate_name='y', variable_name='z')

max_elev = Function(P1, name="max_elev")
max_depth =Function(P1, name="max_depth")
max_elev_file = File(os.path.join(output_directory,'max_elev.pvd'))
max_depth_file = File(os.path.join(output_directory,'max_depth.pvd'))

manning_drag_coefficient = Constant(bottom_friction)

with timed_stage('initialisation'):
  # --- create solver ---
  solverObj = solver2d.FlowSolver2d(mesh2d, bathymetry2d)
  options = solverObj.options
  options.cfl_2d = 1.0
  options.use_nonlinear_equations = True
  options.simulation_export_time = t_export
  options.simulation_end_time = t_end
  options.output_directory = output_dir
  options.check_volume_conservation_2d = True
  options.fields_to_export = ['uv_2d','elev_2d']
  # spatial discretisation
  options.element_family = "dg-dg"
  options.swe_timestepper_type = 'DIRK33'
  #options.timestepper_options.implicitness_theta = 1.0
  #options.timestepper_options.use_semi_implicit_linearization = True
  options.manning_drag_coefficient = manning_drag_coefficient
  options.horizontal_viscosity = Constant(viscosity)
  options.wind_stress = None
  options.timestep = timestep
  options.use_wetting_and_drying = use_wetting_and_drying
  options.verbose = 1
  if use_wetting_and_drying:
    options.wetting_and_drying_alpha = Constant(alpha)
  #options.timestepper_options.solver_parameters = {
  #  'snes_type': 'newtonls',
  #  'snes_rtol': 1e-6,
  #  'snes_linesearch_type': 'bt',
  #  'snes_max_it': 10,
  #  'ksp_type': 'preonly',
  #  'snes_converged_reason': None,
  #  'ksp_converged_reason': None,
  #  'pc_type': 'lu',
  #  'pc_factor_mat_solver_package': 'mumps',
  #}

# boundary conditions
boundary_elev = Function(FunctionSpace(mesh2d, "CG", 1), name='boundary_elev')
solverObj.bnd_functions['shallow_water'] = {
    physID: {'elev': boundary_elev},
    7: {'un': 0.0},
    8: {'un': 0.0},
    9: {'un': 0.0},
}



# Set up usual SWE terms
solverObj.create_equations()


def update_forcings(t):
    t += t_start
    boundary_forcing(boundary_elev, t)
    
    if t>t_start:
        # Update max elev and depth
        eta = Function(P1).project(elev)
        max_elev.dat.data[:] = np.maximum(max_elev.dat.data, eta.dat.data)
        max_depth.dat.data[:] = np.maximum(max_depth.dat.data, eta.dat.data+bathymetry2d.dat.data)
        if (t % t_export == 0):
            max_elev_file.write(max_elev)
            max_depth_file.write(max_depth)
        
update_forcings(0.0)
# set initial condition for elevation, piecewise linear function
solverObj.assign_initial_conditions(uv=Constant((1.0e-6,0.0)), elev=Constant(1.0e-6))

# set up our detectors
cb = DetectorsCallback(solverObj, dets, ['elev_2d'], name='detectors')
solverObj.add_callback(cb,'timestep')


uv, elev = solverObj.timestepper.solution.split()
# Run model

solverObj.iterate(update_forcings=update_forcings)


