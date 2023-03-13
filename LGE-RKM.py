import numpy as np
import matplotlib.pyplot as plt

# Define the logistic growth function and EGF-EGFR function


def logistic_growth(N, r, K, EGF_EGFR):
    return r * N * (1 - N/K) * EGF_EGFR


def egf_egfr_conc(egf_egfr, d_func, s_func, t):
    return -d_func(t) * egf_egfr / (s_func(t) + egf_egfr)


# Set the initial conditions
N0 = 1000
EGF_EGFR0 = 0.1

# Set the parameter values
r = 0.2
K = 10000
d = 0.1
s = 0.5

# Set the time step size and simulation time
dt = 0.01
t_max = 50

# Initialize the time, population, and EGF-EGFR concentration arrays
t = np.arange(0, t_max, dt)
N = np.zeros_like(t)
EGF_EGFR = np.zeros_like(t)

# Set the initial values
N[0] = N0
EGF_EGFR[0] = EGF_EGFR0


def s(t):
    return 2 + 0.5*np.sin(2*np.pi*t/10)


def d(t):
    return 0.1 + 0.05*np.sin(2*np.pi*t/5)


# Use the fourth-order Runge-Kutta method to solve the system of differential equations
for i in range(1, len(t)):
    k1_N = logistic_growth(N[i-1], r, K, EGF_EGFR[i-1])
    k1_EGF_EGFR = egf_egfr_conc(EGF_EGFR[i-1], d, s, t[i-1])

    k2_N = logistic_growth(N[i-1] + k1_N * dt/2, r, K,
                           EGF_EGFR[i-1] + k1_EGF_EGFR * dt/2)
    k2_EGF_EGFR = egf_egfr_conc(
        EGF_EGFR[i-1] + k1_EGF_EGFR * dt/2, d, s, t[i-1] + dt/2)

    k3_N = logistic_growth(N[i-1] + k2_N * dt/2, r, K,
                           EGF_EGFR[i-1] + k2_EGF_EGFR * dt/2)
    k3_EGF_EGFR = egf_egfr_conc(
        EGF_EGFR[i-1] + k2_EGF_EGFR * dt/2, d, s, t[i-1] + dt/2)

    k4_N = logistic_growth(N[i-1] + k3_N * dt, r, K,
                           EGF_EGFR[i-1] + k3_EGF_EGFR * dt)

    k4_EGF_EGFR = egf_egfr_conc(
        EGF_EGFR[i-1] + k3_EGF_EGFR * dt, d, s, t[i-1] + dt)

    N[i] = N[i-1] + (k1_N + 2*k2_N + 2*k3_N + k4_N) * dt/6
    EGF_EGFR[i] = EGF_EGFR[i-1] + \
        (k1_EGF_EGFR + 2*k2_EGF_EGFR + 2*k3_EGF_EGFR + k4_EGF_EGFR) * dt/6


plt.figure()
plt.plot(t, N, label='Population size')
plt.plot(t, EGF_EGFR, label='EGF-EGFR concentration')
plt.xlabel('Time')
plt.ylabel('Population size / EGF-EGFR concentration')
plt.legend()
plt.show()
