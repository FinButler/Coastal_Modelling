import h5py
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from thetis import *

# Load laboratory data and assign desired outputs to plot
df = pd.read_excel('Conical Island Benchmarking.xlsx', sheet_name='TR A')
options.output_directory = 'conical_island_outputs_nh_0.32_f0.0015_v0.001'


lab_stations = ['g6_a', 'g9_a', 'g16_a', 'g22_a']
stations = ['stationA', 'stationB', 'stationC', 'stationD']

fig, axes = plt.subplots(2, 2, figsize=(12, 6))  # Create a 2x2 grid of subplots

for i, (station, lab_station) in enumerate(zip(stations, lab_stations)):
    station_file = options.output_directory + f'/diagnostic_timeseries_{station}_elev.hdf5'

    values = df[[lab_station]]
    num_groups = len(values) // 13
    reshaped_values = values[:num_groups * 13].values.reshape((num_groups, 13))
    averaged_values = reshaped_values.mean(axis=1)
    averaged_df = pd.DataFrame({'Averaged_Value': averaged_values})

    with h5py.File(station_file, 'r') as h5file:
        time_station = h5file['time'][:].flatten()
        vals_station = h5file['elev'][:].flatten()

    row, col = divmod(i, 2)
    ax = axes[row, col]

    ax.plot(time_station, vals_station, 'k:', zorder=3, label=f'Modelled', lw=1.3)
    ax.plot(np.arange(0, 57.5, 0.5) - 3, (averaged_df['Averaged_Value'] - averaged_values[0]), 'r-', zorder=3,
            label=f'Observed', lw=1.3)

    ax.set_title(f'Elevation at {station}', size='small')
    ax.set_ylabel('Elevation (m)')
    ax.set_xlabel('Time (s)')
    ax.set_ylim(-0.025, 0.025)
    ax.set_xlim(0, 50)
    ax.legend()

fig.tight_layout()
plt.show()