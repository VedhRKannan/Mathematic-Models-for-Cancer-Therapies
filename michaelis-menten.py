
############################# Classic Example #####################################################

from scipy.optimize import minimize
from lmfit import Model, Parameters
import pylab
import matplotlib.pyplot as plt
import numpy as np

# Michaelis-Menten constants (M), maximum rate(M.s-1)
Km, Vmax = 0.0033, 20

# Grid of substrate concentration values (M)
E = np.linspace
S = np.linspace(0, 0.0111, 1000)


def rate(S, Km, Vmax):
    return Vmax * S / (Km + S)

# use matplot lib to plot the product on the x axis and the rate of reaction on the y axis

# plt.plot(S, rate(S, Km, Vmax))
# plt.x_label('Substrate')

# pylab.plot(S, rate(S, Km, Vmax))
# pylab.show()


############################## Fitting a Model ###############################
def v(s, v_max, k_m):
    return (v_max * s) / (k_m + s)


data = np.array([[3.6, 1.8, 0.9, 0.45, 0.225, 0.1125, 3.6, 1.8, 0.9, 0.45, 0.225, 0.1125, 3.6, 1.8, 0.9, 0.45, 0.225, 0.1125, 0],
                 [0.004407692, 0.004192308, 0.003553846, 0.002576923, 0.001661538, 0.001064286, 0.004835714, 0.004671429, 0.0039, 0.002857143, 0.00175, 0.001057143, 0.004907143, 0.004521429, 0.00375, 0.002764286, 0.001857143, 0.001121429, 0]]).T

v_real = data[:, 1]
s_real = data[:, 0]


def loss(theta):
    v_max, k_m = theta
    v_pred = v(s_real, v_max, k_m)
    return np.sum((v_real - v_pred)**2)


res = minimize(loss, [1, 1])
print(res.x)

plt.scatter(s_real, v_real)
s_plot = np.linspace(0, 4, 100)
plt.plot(s_plot, v(s_plot, res.x[0], res.x[1]))
plt.xlim([0, 4])
plt.ylim([0, 0.006])
plt.xlabel('[S]')
plt.ylabel('v')
plt.show()
