# -*- coding: utf-8 -*-
"""
Created on Sat May 28 00:10:15 2016

@author: martin
"""
import numpy as np
from numpy import sin, cos, pi
from scipy.integrate import odeint
import matplotlib.pyplot as plt

la = .8
lb = .5
I = 5.0
m = .61
mp = 30.
g = 9.8
h0 = 3.0

def f(x, t):
    alpha, beta, alphap, betap = x #alpha p = alpha punto

    A = I/m + la*la*(1.+mp/m)
    B = la*lb*cos(alpha-beta)
    C = la*lb*cos(alpha-beta)
    D = lb*lb
    E = -(mp/m)*la*(la*alphap+g)-g*la*sin(alpha)-la*lb*betap*betap*sin(alpha-beta) 
    F = la*lb*alphap*alphap*sin(alpha-beta)-g*lb*sin(beta) 
    det = (I/m + la*la*sin(alpha-beta)*sin(alpha-beta)+la*la*(mp/m))*lb*lb
    det = 1.0/det

    dxdt = [alphap, betap, (D*E-B*F)*det, (-C*E + A*F)*det]
    return dxdt

c_ini = [pi/2, pi/2, 0., 0.]
t = np.linspace(0., 1.15, 101)

sol = odeint(f, c_ini, t)

a = sol[:,0] ; b = sol[:,1]
xa = la*sin(a); ya = -la*cos(a) #-
xb = lb*sin(b) + xa; yb = -lb*cos(b) + ya #-
theta = np.linspace(0., 2*pi, 361)
xc = la*cos(theta) ; yc = la*sin(theta)
xd = [3.]*len(sol); yd = h0+la*(a-c_ini[0])
plt.axes().set_aspect('equal')
plt.xlim(-2,4)
plt.plot(xa, ya, 'go')
plt.plot(xb, yb, 'ro')
plt.plot(xc, yc, 'c-')
plt.plot(xd, yd, 'ko')
for i in range(0,len(xa)):
	plt.plot([xa[i], xb[i]],[ya[i], yb[i]], 'b-')
plt.show()
