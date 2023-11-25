# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pyarrow.feather as feather 

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches



# %%
#Import data
score_meta_4subtypes = feather.read_feather("output/score_4subtype_meta.feather")

# %%
# Calculate the frequency of each unique combination of x, y, and z
freq_table = score_meta_4subtypes.groupby(['AClike', 'NPClike', 'OPClike']).size().reset_index(name='count')

# %%
# Merge the frequencies back into the original DataFrame
data = pd.merge(score_meta_4subtypes, freq_table, on=['AClike', 'NPClike', 'OPClike'])  

# Normalize the counts to be between 1 and 10
data['normalized_count'] = np.interp(data['count'], (data['count'].min(), data['count'].max()), (0.1, 100))


# Adjust the size based on the count
sizes = data['normalized_count']


# %%
# Create a custom colormap from white to purple
cmap = LinearSegmentedColormap.from_list('grey_to_purple', ['grey', 'purple'])

# %%
""" # Create a 3D scatter plot with variable point sizes
fig = plt.figure(figsize=(8, 8)) 
ax = plt.axes(projection ="3d")
# Scatter plot
scatter = ax.scatter3D(data['AClike'], data['NPClike'], data['OPClike'], s=sizes, c=data['MESlike'], cmap=cmap)
# Add labels
ax.set_xlabel('AClike')
ax.set_ylabel('NPClike')
ax.set_zlabel('OPClike')

#plt.tight_layout(pad=100.0)
# Add color bar
#cbar = fig.colorbar(scatter, ax=ax)

plt.show() """

# Separate the data based on the 'radiation' category
radiated = data[data['radiation'] == 'radiated']
control = data[data['radiation'] == 'control']

# Create a 3D scatter plot with variable point sizes
fig = plt.figure(figsize=(10, 10)) 
ax = plt.axes(projection ="3d")

# Scatter plot for 'radiated' and 'control' with different markers/colors
scatter_radiated = ax.scatter3D(radiated['AClike'], radiated['NPClike'], radiated['OPClike'], s=radiated['normalized_count'], c=radiated['MESlike'], cmap=cmap, marker='o')
scatter_control = ax.scatter3D(control['AClike'], control['NPClike'], control['OPClike'], s=control['normalized_count'], c=control['MESlike'], cmap=cmap, marker='^')

# Add labels
ax.set_xlabel('AClike')
ax.set_ylabel('NPClike')
ax.set_zlabel('OPClike')

# Colorbar for MESlike values
cbar = fig.colorbar(scatter_radiated, ax=ax, label='MESlike')

# Create a custom legend for sizes
legend_elements = [mpatches.Circle((0, 0), radius=s, color='gray', label=str(c)) for s, c in zip([10, 50, 100], [20, 100, 500])]
ax.legend(handles=legend_elements, title="Sizes (Counts)")

# Show the plot
plt.show()

