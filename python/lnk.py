"""
Linear-Nonlinear-Kinetics (LNK) model for moving bar and playback experiments
"""

import numpy as np
import scipy.special

__all__ = ['compute_est', 'simulate_linear scale',
        'simulate_nonlinearity', 'simulate_4state',
        'simulate_nonlinearity_2', 'simulate_4state_2']

# params for filter (15), nonlinearity(5), kinetics(6), scaling(1)
param = [-0.0054, 0.0055, 0.0106, -0.0078, 0.0036,    \
        -0.0118, -0.0045, -0.0038, -0.0014, 9.2057e-4,\
        0.0612, 0.1051, 0.0490, -0.0415, 0.0099,      \
        2.5458, -0.7199, 2.1938, 0.0087, 0.1241,      \
        0, 9.8902, 0, 62.9771, 0, 6.1664,             \
        5.0552]
x_0 = np.array([0, 0, 0, 100]);

def compute_est():
    return 0;

def simulate_linear_scale(lin, p=1):
    LN = lin
    LN_2 = p*lin
    return LN, LN_2

def simulate_nonlinearity(lin, p=param[15:20]):
    # nonlinearity block
    # {lin, p}.shape: (time(T),), (5,)
    LN = np.power(p[0], scipy.special.erf(lin + p[1]) + 1) + p[2]
    LN_2 = p[3]*LN + p[4]
    return LN, LN_2

# final result is from this nonlinearity
def simulate_nonlinearity_2(g, p):
    # nonlinearity block 2
    # {g, p}.shape: (time(T),), (2,)
    u = (scipy.special.erf(g + p[0]) + 1)/2 + p[1]
    return u

def simulate_4state(u, v, p=param[20:26], x_0=x_0):
    # kinetics block
    # {u, v, p, x_0}.shape: (time(T),), (T,), (6,), (4,)
    x_curr = x_0;
    T = u.shape[0]
    X = np.zeros((4, T))
    for t in range(T):
        X[:, t] = x_curr
        x_next = np.zeros((4,))
        x_next[0] = x_curr[0]*(1-0.001*(u[t]+p[4])) + x_curr[2]*0.001*p[3] + x_curr[1]*0.001*p[0]
        x_next[1] = x_curr[1]*(1-0.001*(p[0]+p[1])) + x_curr[0]*0.001*u[t] + x_curr[2]*0.001*p[2]
        x_next[2] = x_curr[2]*(1-0.001*(p[2]+p[3]+p[5])) + x_curr[1]*0.001*p[1] + \
                x_curr[0]*0.001*p[4] + x_curr[3]*0.001*v[t]
        x_next[3] = x_curr[3]*(1-0.001*v[t]) + x_curr[2]*0.001*p[5]
        x_curr = x_next
    return X

# final result is from this kinetics block
def simulate_4state_2(u, p, x_0=x_0):
    # kinetics block
    # {u, p, x_0}.shape: (time(T),), (7,), (4,)
    x_curr = x_0
    T = u.shape[0]
    X = np.zeros((4, T))
    k_a = p[0]
    k_fi = p[1]
    k_fr = p[2]
    k_si = p[3]
    k_sr = p[4]
    b_1 = p[5]
    b_2 = p[6]
    for t in range(T):
        X[:, t] = x_curr
        Q = np.array([
            [-u[t]*k_a-b_1, u[t]*k_a+b_1,              0,              0],
            [            0,        -k_fi,           k_fi,              0],
            [         k_fr,            0,   -(k_fr+k_si),           k_si],
            [            0,            0,  u[t]*k_sr+b_2, -u[t]*k_sr-b_2]])
        x_curr = x_curr + 0.001 * np.matmul(x_curr.reshape(1, 4), Q).reshape(4,)
    return X
