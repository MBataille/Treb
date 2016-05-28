"""
Created on Wed May 25 21:32:37 2016

@author: martin

Simulacion del disco
"""

import numpy as np
from numpy import sin, cos, pi
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

la = 1.0
lb = 1.0
I = 5.0
m = 1.0
g = 9.8

def f(x, t):
    alpha, beta, alphap, betap = x #alpha p = alpha punto

    A = I/m + la*la
    B = la*lb*cos(alpha-beta)
    C = la*lb*cos(alpha-beta)
    D = lb*lb
    E = -g*la*sin(alpha)-la*lb*betap*betap*sin(alpha-beta) 
    F = la*lb*alphap*alphap*sin(alpha-beta)-g*lb*sin(beta) 
    det = (I/m + la*la*sin(alpha-beta)*sin(alpha-beta))*lb*lb
    det = 1.0/det

    dxdt = [alphap, betap, (D*E-B*F)*det, (-C*E + A*F)*det]
    return dxdt

c_ini = [pi, pi/2, 0., 0.]
t = np.linspace(0., 20., 200)

sol = odeint(f, c_ini, t)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1) 
ax1.set_aspect('equal')
ax1.set_xlim(-2,2)
ax1.set_ylim(-2,2)
def animate(i):
    xas = la*sin(sol[i,0])
    yas = -la*cos(sol[i,0])
    xbs = lb*sin(sol[i,1]) + xas
    ybs = -lb*cos(sol[i,1]) + yas
    thisx = [0, xas, xbs]
    thisy = [0, yas, ybs]
    theta = np.linspace(0, 2*pi, 361)
    xcs = la*sin(theta)
    ycs = la*cos(theta)
    ax1.clear()
    ax1.set_xlim(-2,2)
    ax1.set_ylim(-2,2)
    ax1.plot(thisx,thisy, 'bo-')
    ax1.plot(xcs, ycs, 'g-')
#a = sol[:,0] ; b = sol[:,1]
#print (sol)
#xa = la*sin(a); ya = la*cos(a) #-
#xb = lb*sin(b) + xa; yb = lb*cos(b) + ya #-
#theta = np.linspace(0., 2*pi, 361)
#xc = la*cos(theta) ; yc = la*sin(theta)
#
#plt.axes().set_aspect('equal')
#plt.plot(xa, ya, 'go')
#plt.plot(xb, yb, 'ro')
#plt.plot(xc, yc, 'c-')
#for i in range(0,len(xa)):
#	plt.plot([xa[i], xb[i]],[ya[i], yb[i]], 'b-')
ani = animation.FuncAnimation(fig,animate,np.arange(1,len(sol)),interval=50)
plt.show()

resp = input("Guardar [S/N]")
if resp == "S":
    name = input("Nombre del archivo")
    ani.save('{}.mp4'.format(name))
