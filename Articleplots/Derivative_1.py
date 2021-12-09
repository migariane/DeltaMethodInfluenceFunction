import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


# 100 linearly spaced numbers
x  = np.linspace(-1,2.5,100)
a  = 0.4
b  = 1.5

# the function, which is y = x^2 here

yfun         = lambda x: np.exp(x)
yprime       = lambda x: np.exp(x)
tangent_line = lambda x: yprime(a)*(x - a) + yfun(a)
line         = lambda x: (yfun(b) - yfun(a))/(b - a)*(x - a) + yfun(a)
line2        = lambda x: (yfun(c) - yfun(a))/(c - a)*(x - a) + yfun(a)

# setting the axes at the centre
fig = plt.figure()
fig.set_size_inches(8, 6)
ax = fig.add_subplot(1, 1, 1)
#ax.spines['left'].set_position('zero')
#ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.set_yticks([])
ax.set_xticks([])
plt.ylim(bottom=0, top = yfun(2))
plt.xlim(-1,3)
# plot the function
plt.plot(x,yfun(x), 'k-',zorder=2)
plt.plot(x, line(x), 'r', lw=2, zorder=1, linestyle = "dashed")
#plt.plot(x, line2(x), 'b', lw=2, zorder=1, linestyle = "dashed")
plt.plot(x, tangent_line(x), 'g', lw=2, zorder=1, linestyle = "dotted")
plt.scatter(b, yfun(b), s = 100, c = 'r',zorder=101)
plt.scatter(a, yfun(a), s = 100, c = 'r',zorder=100)
#plt.scatter(c, yfun(c), s = 100, c = 'b',zorder=100)
plt.plot([a, a], [yfun(a), 0], 'r', lw=0.5, zorder=4, linestyle='dashed')
plt.plot([b, b], [yfun(b), 0], 'r', lw=0.5, zorder=4, linestyle='dashed')
plt.plot([-1, a], [yfun(a), yfun(a)], 'r', lw=0.5, zorder=4, linestyle='dashed')
plt.plot([-1, b], [yfun(b), yfun(b)], 'r', lw=0.5, zorder=4, linestyle='dashed')
#plt.plot([c, c], [yfun(c), 0], 'b', lw=0.5, zorder=4, linestyle='dashed')

ax.annotate(r'$\phi(x)$',
            xy=(2, yfun(2)), xycoords='data', fontsize=12, annotation_clip=False,
            horizontalalignment='left', verticalalignment='top')


ax.annotate(r'$\frac{\phi(\theta + h) - \phi(\theta)}{h}\cdot(x - \theta) + \phi(\theta)$', color = "red",
            xy=(2, line(2)), xycoords='data', fontsize=12, annotation_clip=False,
            horizontalalignment='left', verticalalignment='top')


ax.annotate(r"$\mathcal{T}_1(x) =\phi'(\theta)\cdot(x - \theta) + \phi(\theta)$", color = "green",
            xy=(2, tangent_line(2)), xycoords='data', fontsize=12, annotation_clip=False,
            horizontalalignment='left', verticalalignment='top')


ax.annotate(r'$\phi(\theta)$',
            xy=(-1.2, yfun(a)), xycoords='data', fontsize=12, annotation_clip=False,
            horizontalalignment='right', verticalalignment='top')

ax.annotate(r'$\phi(\theta + h)$',
            xy=(-1.2, yfun(b)), xycoords='data', fontsize=12, annotation_clip=False,
            horizontalalignment='right', verticalalignment='top')


ax.annotate(r'$\theta$',
            xy=(a, -0.1), xycoords='data', fontsize=12, annotation_clip=False,
            horizontalalignment='center', verticalalignment='top')

ax.annotate(r'$\theta + h_1$',
            xy=(b, -0.1), xycoords='data', fontsize=12, annotation_clip=False,
            horizontalalignment='center', verticalalignment='top')

#ax.annotate(r'$\theta + h_2$',
#            xy=(c, -0.1), xycoords='data', fontsize=16, annotation_clip=False,
#            horizontalalignment='center', verticalalignment='top')

# show the plot
plt.show()
plt.savefig("Figura_1.pdf")
