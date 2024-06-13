import utm
import uptide
import uptide.tidal_netcdf
import datetime
import scipy.interpolate
import numpy as np
from math import tanh, sin

# edit me!
utm_zone=30
utm_band='U'
fluidity_times = np.loadtxt('../bc_extraction/times.txt')
boundary_coords = np.loadtxt('../bc_extraction/coordinates.txt')
boundary_data = np.loadtxt('../bc_extraction/elev.txt')


def set_tidal_field(elev, t):
    if M2_only:
        constituents = ['M2']
    else:
        constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
    tide = uptide.Tides(constituents)
    tide.set_initial_time(datetime.datetime(2013,11,15,0,0,0))
    tnci = uptide.tidal_netcdf.OTPSncTidalInterpolator(tide, '../../../humber_model/data/gridES2008.nc', '../../../humber_model/data/hf.ES2008.nc')
    tnci.set_time(t)
    mesh2d = elev.function_space().mesh()
    xvector = mesh2d.coordinates.dat.data
    evector = elev.dat.data
    for i,xy in enumerate(xvector):
        lat, lon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band, strict=False)
        try:
            evector[i] = tnci.get_val((lon, lat))
        except uptide.netcdf_reader.CoordinateError:
            # Occurs when given location is on land
            #print(lat, lon, xy)
            #print(xy[0], xy[1])
            evector[i] = 0.

def set_tsunami_field(elev, t):
    #mesh2d = elev.function_space().mesh()
    #xvector = mesh2d.coordinates.dat.data
    #evector = elev.dat.data
    #print(((10000-(t-1000))/5000)*sin((t-1000)/500)) # linearly decaying sine wave
    #for i,xy in enumerate(xvector):
    #    evector[i] = ((10000.-(t-1000.))/5000.)*sin((t-1000.)/500.)    
    #return 
    
    # Set tidal field using nearest neighbour interpolation from boundary points extracted from Fluidity tsunami simulations
    fluidity_dt = fluidity_times[1]-fluidity_times[0]
    if t%fluidity_dt == 0:
        fluidity_index = int(t/fluidity_dt)
        bd = boundary_data[fluidity_index]
        interpolator = scipy.interpolate.NearestNDInterpolator(boundary_coords, bd)
        mesh2d = elev.function_space().mesh()
        xvector = mesh2d.coordinates.dat.data
        evector = elev.dat.data
        for i,xy in enumerate(xvector):
            evector[i] = interpolator((xy[0], xy[1]))
    else:
        fluidity_index = int(t/fluidity_dt)
        boundary_data_1 = boundary_data[fluidity_index]
        boundary_data_2 = boundary_data[fluidity_index+1]
        interpolator_1 = scipy.interpolate.NearestNDInterpolator(boundary_coords, boundary_data_1)
        interpolator_2 = scipy.interpolate.NearestNDInterpolator(boundary_coords, boundary_data_2)
        mesh2d = elev.function_space().mesh()
        xvector = mesh2d.coordinates.dat.data
        evector = elev.dat.data
        for i,xy in enumerate(xvector):
            evector[i] = interpolator_1((xy[0], xy[1])) + (t-fluidity_dt*fluidity_index)*(interpolator_1((xy[0], xy[1]))-interpolator_2((xy[0], xy[1])))/fluidity_dt

