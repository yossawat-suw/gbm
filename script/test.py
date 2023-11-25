# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyarrow.feather as feather 
from matplotlib.colors import LinearSegmentedColormap

# Read data
score_meta_4subtypes = feather.read_feather("output/score_4subtype_meta.feather")

# Group and merge data for frequency
freq_table = score_meta_4subtypes.groupby(['AClike', 'NPClike', 'OPClike']).size().reset_index(name='count')
data = pd.merge(score_meta_4subtypes, freq_table, on=['AClike', 'NPClike', 'OPClike'])  

# Normalize counts
data['normalized_count'] = np.interp(data['count'], (data['count'].min(), data['count'].max()), (0.1, 10))

# Create colormap
cmap_radiated = LinearSegmentedColormap.from_list('white_to_blue', ['white', 'blue'])
cmap_control = LinearSegmentedColormap.from_list('white_to_red', ['white', 'red'])

# Separate data based on 'radiation' column
radiated = data[data['radiation'] == 'radiated']
control = data[data['radiation'] == 'control']

# Apply jitter to one of the groups
# Here, we apply a slight offset to the 'control' group along the 'AClike' axis
control_offset = control.copy()
control_offset['AClike'] = control['AClike'] + 1  # Adjust this offset value as needed

# Create 3D scatter plot
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot for 'radiated' and 'control' with offset
scatter1 = ax.scatter3D(radiated['AClike'], radiated['NPClike'], radiated['OPClike'], s=radiated['normalized_count'], c=radiated['MESlike'], cmap=cmap_radiated, marker='o', alpha=0.5)
scatter2 = ax.scatter3D(control_offset['AClike'], control_offset['NPClike'], control_offset['OPClike'], s=control_offset['normalized_count'], c=control_offset['MESlike'], cmap=cmap_control, marker='^', alpha=0.5)

# Labels
ax.set_xlabel('AClike')
ax.set_ylabel('NPClike')
ax.set_zlabel('OPClike')

# Colorbar (for MESlike)
cbar_radiated = fig.colorbar(scatter1, ax=ax, shrink=0.5)
cbar_radiated.set_label('MESlike')
cbar_control = fig.colorbar(scatter2, ax=ax, shrink=0.5)
cbar_control.set_label('MESlike')

# Min and Max Count Annotations
min_count = data['count'].min()
max_count = data['count'].max()
ax.text2D(0.05, 0.95, f"Min Count: {min_count}", transform=ax.transAxes)
ax.text2D(0.05, 0.90, f"Max Count: {max_count}", transform=ax.transAxes)

plt.show()
