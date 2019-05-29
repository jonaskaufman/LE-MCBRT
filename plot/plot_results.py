import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Plot parameters
colormap = "plasma"
max_density = 1.0

# Nice colorbars (from https://joseph-long.com/writing/colorbars/)
def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

# Specify data paths
densities_path="../source/densities.csv"
doses_path="../source/doses.csv"

# Load data
densities = np.loadtxt(densities_path, delimiter=',')
doses = np.loadtxt(doses_path, delimiter=',')

# Plot data
f, (ax1, ax2) = plt.subplots(1, 2)

ax1.set_title("density")
im1 = ax1.imshow(densities, cmap=colormap, vmin=0.0, vmax=max_density)
colorbar(im1)
ax1.axis('off')

ax2.set_title("dose")
im2 = ax2.imshow(doses, cmap=colormap, vmin=0.0)
colorbar(im2)
ax2.axis('off')

plt.tight_layout()
plt.show()

