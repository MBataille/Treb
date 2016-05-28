"""
Created on Wed May 25 21:32:37 2016

@author: martin

Simulacion del disco con peso
"""

import numpy as np
from numpy import sin, cos, pi
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

la = 1
lb = 0.7
I = 10.0
m = .61
mp = 30.
g = 9.8
h0 = 4.5

def abs(ar):
    x = []
    for a in ar:
        if a>0:
            x.append(a)
        else:
            x.append(-a)
    return x

def v(sol, t):
    vel = []
    for i in range(len(t)):
        vel.append( np.sqrt( la*la * sol[i,2]* sol[i,2] +lb*lb*sol[i,3]*sol[i,3] + 2*la*lb*sol[i,2]*sol[i,3]*np.cos( sol[i,0]-sol[i,1] ) ) )
    return vel

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

c_ini = [pi/2, 0, 0., 0.]
t = np.linspace(0., 1.3, 201)

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
    xds = 4
    yds = h0+la*(sol[i,0]-sol[0,0])
    if yds <=0:
        return
    ax1.clear()
    ax1.set_xlim(-4,6)
    ax1.set_ylim(-4,6)
    ax1.plot(thisx,thisy, 'bo-')
    ax1.plot(xcs, ycs, 'g-')
    ax1.plot(xds, yds, 'ko')

ani = animation.FuncAnimation(fig,animate,np.arange(1,len(sol)),interval=25)
plt.show()
fig = plt.figure()
### V_B
ax2 = fig.add_subplot(211)
ax3 = fig.add_subplot(212)
ax2.plot(t, v(sol,t))
lb = la
sol = odeint(f, c_ini, t)
ax3.plot(t, v(sol,t))
plt.show()
resp = input("Guardar [S/N]")
if resp == "S":
    name = input("Nombre del archivo ")
    ani.save('{}.mp4'.format(name))
