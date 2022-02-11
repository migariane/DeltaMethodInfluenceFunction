import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from numpy.random import uniform
from numpy.random import poisson
from numpy import hstack
from statsmodels.distributions.empirical_distribution import ECDF

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams.update({
    'text.latex.preamble': r'\usepackage{amsfonts}'
})

#Uniform variables
#create_xbar = lambda n, m: [np.mean(uniform(low=0, high=1, size=n)) for i in range(m)]
create_xbar = lambda n, m: [np.mean(poisson(1, size=n)) for i in range(m)]

#Get samples
n1 = 5
n2 = 10
n3 = 500
sample1 = create_xbar(n1, 1000)
sample2 = create_xbar(n2, 1000)
sample3 = create_xbar(n3, 1000)

#Transform the variables
sigma_sq  = 1
mu        = 1
normalize = lambda n, xbar: np.sqrt(n/sigma_sq)*(xbar - mu)
sample1   = [normalize(n1, xbar) for xbar in sample1]
sample2   = [normalize(n2, xbar) for xbar in sample2]
sample3   = [normalize(n3, xbar) for xbar in sample3]

#Empirical cdf
x     = np.linspace(-2.5, 2.5, 1000)
ecdf1 = ECDF(sample1)(x)
ecdf2 = ECDF(sample2)(x)
ecdf3 = ECDF(sample3)(x)
cdf   = norm.cdf(x)

# setting the axes at the centre
fig = plt.figure()
fig.set_size_inches(8, 6)
ax = fig.add_subplot(1, 1, 1)
#ax.spines['left'].set_position('zero')
#ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.set_ylabel(r'$F_Z(z) = \mathbb{P}(Z \leq z)$', fontsize = 12)
ax.set_xlabel(r'$z$', fontsize = 12)

# plot the function
plt.plot(x,cdf,'k', zorder=5, alpha = 0.35, lw = 2, label=r"$n = \infty$")
plt.plot(x, ecdf1, 'r', lw=1, zorder=2, linestyle = "dashed", alpha = 0.5, label=r"$n = 10$")
plt.plot(x, ecdf2, 'g', lw=1, zorder=3, linestyle = "dotted", alpha = 0.5, label=r"$n = 20$")
plt.plot(x, ecdf3, 'b', lw=1, zorder=4, linestyle = "dashdot", alpha = 0.5, label=r"$n = 30$")
plt.legend()


# show the plot
plt.show()
plt.savefig("ECDF.pdf")
