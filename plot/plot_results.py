import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

# Plot parameters
colormap = "plasma"
log_scale_dose = False
log_zero_offset = 0.000000001
save = False
DPI = 300;

if 'log' in sys.argv:
    log_scale_dose = True

if 'save' in sys.argv:
    save = True

# Nice colorbars (from https://joseph-long.com/writing/colorbars/)
def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

# Specify data paths
densities_path="densities.csv"
doses_path="doses.csv"

# Load data
densities = np.loadtxt(densities_path, delimiter=',')
doses = np.loadtxt(doses_path, delimiter=',')
if log_scale_dose: doses = np.log(doses + log_zero_offset)

# Plot data
f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(6, 3))

ax1.set_title("density")
im1 = ax1.imshow(densities, cmap=colormap, vmin=0.0, vmax=1.0)
colorbar(im1)
ax1.axis('off')

dose_title = "dose"
if log_scale_dose: dose_title = "ln(dose)"
# TODO: redo log scale using LogNorm colorbar
ax2.set_title(dose_title)
im2 = ax2.imshow(doses, cmap=colormap, vmin=0.0)
colorbar(im2)
ax2.axis('off')

plt.tight_layout()
if save:
    plt.savefig('results.png', dpi=DPI)
else:
    plt.show()

