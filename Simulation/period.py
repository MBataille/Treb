# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 12:37:36 2016

@author: martin
"""

import numpy as np
from numpy import sin, cos, pi
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

la = .61
lb = 0.
m = .61
mp = 30.
g = 9.8
h0 = 2.0

def f(x, t):
    alpha, alphap = x

    A = -g*la*sin(alpha)
    B = (I/m)+la*la

    dxdt = [alphap, A/B]
    return dxdt

def period(sol):
    epsilon = 0.001
    alpha_init = sol[0,0]
    for i in range(10,len(sol)):
        if sol[i,0] >= (alpha_init - epsilon) and sol[i,0] <= (alpha_init + epsilon):          
            T = i*0.025
            break
    return T
c_ini = [pi/2, 0.]
t = np.linspace(0., 5., 201)
Im = np.linspace(0.01,0.7,201)
T = []
for I in Im:
    sol = odeint(f, c_ini, t)
    T.append(period(sol))
plt.plot(T, Im)
plt.show()

