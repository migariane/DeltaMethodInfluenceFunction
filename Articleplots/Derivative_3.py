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
X = np.arange(-2, 2, 0.01)
Y = np.arange(-2, 2, 0.01)
X, Y = np.meshgrid(X, Y)

phi = lambda X,Y: (X/2)**2  + 0.75

Z = phi(X,Y)

# Plot the surface.
ax.plot_wireframe(X, Y, Z, color = "black",
                       linewidth=1, antialiased=False, zorder = 0, alpha = 0.15)
ax.plot([0],[0],zs = [phi(0,0)] , color = "r", zorder = 10, marker='o', markersize=10)
x1 = -(1.2*np.cos(np.pi*5/6) - 1.2*np.sin(np.pi*5/6))*0.5
y1 = -(1.2*np.cos(np.pi*5/6) + 1.2*np.sin(np.pi*5/6))*0.5
ax.plot([0, 0],[0,0],zs = [0, phi(0,0)] , color = "r", zorder = 1, linewidth = 0.5, linestyle = "dotted")
ax.plot([0,x1], [0,y1],[phi(0,0),phi(x1, y1)], color = "blue", linestyle = "dotted", linewidth = 2)

x1 = -(1.2*np.cos(np.pi*2/3) - 1.2*np.sin(np.pi*2/3))*0.5
y1 = -(1.2*np.cos(np.pi*2/3) + 1.2*np.sin(np.pi*2/3))*0.5
ax.plot([0, 0],[0,0],zs = [0, phi(0,0)] , color = "r", zorder = 1, linewidth = 0.5, linestyle = "dotted")
ax.plot([0,x1], [0,y1],[phi(0,0),phi(x1, y1)], color = "blue", linestyle = "dotted", linewidth = 2)


x1 = -(1.2*np.cos(np.pi) - 1.2*np.sin(np.pi))*0.5
y1 = -(1.2*np.cos(np.pi) + 1.2*np.sin(np.pi))*0.5
ax.plot([0, 0],[0,0],zs = [0, phi(0,0)] , color = "r", zorder = 1, linewidth = 0.5, linestyle = "dotted")
ax.plot([0,x1], [0,y1],[phi(0,0),phi(x1, y1)], color = "blue", linestyle = "dotted", linewidth = 2)


x1 = -(1.2*np.cos(np.pi*2/5) - 1.2*np.sin(np.pi*2/5))*0.5
y1 = -(1.2*np.cos(np.pi*2/5) + 1.2*np.sin(np.pi*2/5))*0.5
ax.plot([0, 0],[0,0],zs = [0, phi(0,0)] , color = "r", zorder = 1, linewidth = 0.5, linestyle = "dotted")
ax.plot([0,x1], [0,y1],[phi(0,0),phi(x1, y1)], color = "blue", linestyle = "dotted", linewidth = 2)


x1 = -(1.2*np.cos(np.pi/2) - 1.2*np.sin(np.pi/2))*0.5
y1 = -(1.2*np.cos(np.pi/2) + 1.2*np.sin(np.pi/2))*0.5
ax.plot([0, 0],[0,0],zs = [0, phi(0,0)] , color = "r", zorder = 1, linewidth = 1, linestyle = "dotted")
ax.plot([0,x1], [0,y1],[phi(0,0),phi(x1, y1)], color = "blue", linestyle = "dotted", linewidth = 2)


x1 = -(1.2*np.cos(np.pi/3) - 1.2*np.sin(np.pi/3))*0.5
y1 = -(1.2*np.cos(np.pi/3) + 1.2*np.sin(np.pi/3))*0.5
ax.plot([0, 0],[0,0],zs = [0, phi(0,0)] , color = "r", zorder = 1, linewidth = 0.5, linestyle = "dotted")
ax.plot([0,x1], [0,y1],[phi(0,0),phi(x1, y1)], color = "blue", linestyle = "dotted", linewidth = 2)


ax.plot([0, 0],[-1,1],zs = [0.75, 0.75] , color = "r", zorder = 1, linewidth = 2, linestyle = "dashed")
ax.quiver(0,0,0,-(1.2*np.cos(np.pi/4) - 1.2*np.sin(np.pi/4)),-(1.2*np.cos(np.pi/4) + 1.2*np.sin(np.pi/4)),-0, length=0.5, normalize=False, color = "r")
ax.quiver(0,0,0,-(1.2*np.cos(np.pi/2) - 1.2*np.sin(np.pi/2)),-(1.2*np.cos(np.pi/2) + 1.2*np.sin(np.pi/2)),-0, length=0.5, normalize=False, color = "blue")
ax.quiver(0,0,0,-(1.2*np.cos(np.pi/3) - 1.2*np.sin(np.pi/3)),-(1.2*np.cos(np.pi/3) + 1.2*np.sin(np.pi/3)),-0, length=0.5, normalize=False, color = "blue")
ax.quiver(0,0,0,-(1.2*np.cos(np.pi*2/5) - 1.2*np.sin(np.pi*2/5)),-(1.2*np.cos(np.pi*2/5) + 1.2*np.sin(np.pi*2/5)),-0, length=0.5, normalize=False, color = "blue")
ax.quiver(0,0,0,-(1.2*np.cos(np.pi*2/3) - 1.2*np.sin(np.pi*2/3)),-(1.2*np.cos(np.pi*2/3) + 1.2*np.sin(np.pi*2/3)),-0, length=0.5, normalize=False, color = "blue")
ax.quiver(0,0,0,-(1.2*np.cos(np.pi) - 1.2*np.sin(np.pi)),-(1.2*np.cos(np.pi) + 1.2*np.sin(np.pi)),-0, length=0.5, normalize=False, color = "blue")
ax.quiver(0,0,0,-(1.2*np.cos(np.pi*5/6) - 1.2*np.sin(np.pi*5/6)),-(1.2*np.cos(np.pi*5/6) + 1.2*np.sin(np.pi*5/6)),-0, length=0.5, normalize=False, color = "blue")

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

ax.text(0, -1, 0, r"$v$", color = "red")

plt.tight_layout()
ax.view_init(elev=35, azim=120)
ax.set_box_aspect((1,1,1))

plt.show()
plt.savefig("Figura_3.pdf")
