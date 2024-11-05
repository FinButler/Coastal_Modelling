import h5py
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from thetis import *

# Load laboratory data
df = pd.read_excel('req_data/Conical Island Benchmarking.xlsx', sheet_name='TR A')

# Generate output directories based on the pattern
output_directories = [
    f'outputs/con_island_drag_{drag}_visc_{visc}'
    for drag in [0.001, 0.01, 0.1, 1, 10]
    for visc in [0.001, 0.01, 0.1, 1, 10]
]

lab_stations = ['g6_a', 'g9_a', 'g16_a', 'g22_a']
stations = ['stationA', 'stationB', 'stationC', 'stationD']

fig, axes = plt.subplots(2, 2, figsize=(12, 6))

# Define a list of colors and linestyles to cycle through
colors = ['k', 'b', 'g', 'm', 'c', 'y', 'purple', 'orange', 'brown', 'teal']
linestyles = ['-', '--', '-.', ':']
num_styles = len(colors) * len(linestyles)  # Total number of unique color-style combinations

for i, (station, lab_station) in enumerate(zip(stations, lab_stations)):
    values = df[[lab_station]]
    num_groups = len(values) // 13
    reshaped_values = values[:num_groups * 13].values.reshape((num_groups, 13))
    averaged_values = reshaped_values.mean(axis=1)
    averaged_df = pd.DataFrame({'Averaged_Value': averaged_values})

    row, col = divmod(i, 2)
    ax = axes[row, col]

    # Plot the observed values once per station
    ax.plot(np.arange(0, 57.5, 0.5) - 3, (averaged_df['Averaged_Value'] - averaged_values[0]), 'r-', zorder=3,
            label='Observed', lw=1.3)

    # Loop through each output directory and plot the model data
    for j, output_dir in enumerate(output_directories):
        station_file = f'{output_dir}/diagnostic_timeseries_{station}_elev.hdf5'
        color = colors[j % len(colors)]
        linestyle = linestyles[(j // len(colors)) % len(linestyles)]

        try:
            with h5py.File(station_file, 'r') as h5file:
                time_station = h5file['time'][:].flatten()
                vals_station = h5file['elev'][:].flatten()

            ax.plot(time_station, vals_station, color=color, linestyle=linestyle, zorder=3,
                    label=f'Modelled (drag={output_dir.split("_")[-3]}, visc={output_dir.split("_")[-1]})', lw=1.3)
        except FileNotFoundError:
            print(f"File not found: {station_file}")
            continue

    ax.set_title(f'Elevation at {station}', size='small')
    ax.set_ylabel('Elevation (m)')
    ax.set_xlabel('Time (s)')
    ax.set_ylim(-0.025, 0.025)
    ax.set_xlim(0, 50)
    # ax.legend(loc='upper right', fontsize='x-small')

fig.tight_layout()
plt.show()
