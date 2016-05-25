import numpy as np
from numpy import sin, cos, pi
from scipy.integrate import odeint
import matplotlib.pyplot as plt

la = 1.0
lb = 1.0
I = 5.0
m = 1.0
g = 9.8

def f(x,t):
	alpha, beta, alphap, betap = x #alpha p = alpha punto

	A = I/m + la*la
	B = la*lb*cos(alpha-beta)
	C = la*lb*cos(alpha-beta)
	D = lb*lb
	E = g*la*sin(alpha)-la*lb*betap*betap*sin(alpha-beta) 
	F = la*lb*alphap*alphap*sin(alpha-beta)+g*lb*sin(beta) 
	det = (I/m + la*la*sin(alpha-beta)*sin(alpha-beta))*lb*lb

	dxdt = [alphap, betap, (D/det)*E - (B/det)*F, - (C/det)*E+(A/det)*F]
	return dxdt

c_ini = [pi/2, pi/2, 2., 0.]
t = np.linspace(0., 4., 201)

sol = odeint(f, c_ini, t)

a = sol[:,0] ; b = sol[:,1]
print (sol)
xa = la*sin(a); ya = la*cos(a) #-
xb = lb*sin(b) + xa; yb = lb*cos(b) + ya #-
theta = np.linspace(0., 2*pi, 361)
xc = la*cos(theta) ; yc = la*sin(theta)

plt.axes().set_aspect('equal')
plt.plot(xa, ya, 'go')
plt.plot(xb, yb, 'ro')
plt.plot(xc, yc, 'c-')
for i in range(0,len(xa)):
	plt.plot([xa[i], xb[i]],[ya[i], yb[i]], 'b-')
plt.show()