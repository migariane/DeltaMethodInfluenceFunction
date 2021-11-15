from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

fig = plt.figure()
fig.set_size_inches(8, 6)
ax = fig.gca(projection='3d')

# Make data.
X = np.arange(-1, 1, 0.01)
Y = np.arange(-1, 1, 0.01)
X, Y = np.meshgrid(X, Y)

Z = np.cos(X)**2

# Plot the surface.
ax.plot_wireframe(X, Y, Z, color = "black",
                       linewidth=1, antialiased=False, zorder = 0, alpha = 0.15)
ax.plot([-0.5],[0],zs = [np.cos(0.5)**2 ] , color = "r", zorder = 10, marker='o', markersize=10)
ax.plot([-0.5, -0.5],[0,0],zs = [0, np.cos(0.5)**2] , color = "r", zorder = 1, linewidth = 0.5, linestyle = "dotted")

ax.plot([-0.5, -0.5],[-1,1],zs = [np.cos(0.5)**2, np.cos(0.5)**2] , color = "r", zorder = 1, linewidth = 2, linestyle = "dashed")
ax.plot([0, -1],[0,0],zs = [2*np.cos(0.5)*np.sin(0.5)*(0.5 - 0) + np.cos(0.5)**2, 2*np.cos(0.5)*np.sin(0.5)*(0.5 - 1) + np.cos(0.5)**2] , color = "b", zorder = 1, linewidth = 2, linestyle = "dashed")
ax.quiver(-0.5,0,0,0,-1.2,0, length=0.5, normalize=False, color = "r")
ax.quiver(-0.5,0,0,1.2,0,-0, length=0.5, normalize=False, color = "blue")

# Customize the z axis.
ax.set_zlim(0.0, 1.5)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Get rid of the panes
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))


# Add the labels
ax.set_xlabel(r'$\theta_1$', fontsize = 12)
ax.set_ylabel(r'$\theta_2$', fontsize = 12)
ax.set_zlabel(r'$\phi(\theta_1,\theta_2)$', fontsize = 12)

# make the grid lines transparent
ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)


# Get rid of the ticks
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

ax.text(0, 0, 0.1, r"$v_1$", color = "blue")
ax.text(-0.5, -0.5, 0.1, r"$v_2$", color = "red")

plt.tight_layout()
ax.view_init(elev=5.0, azim=123)
ax.set_box_aspect((1,1,1))

plt.show()
plt.savefig("Figura_2.pdf")
