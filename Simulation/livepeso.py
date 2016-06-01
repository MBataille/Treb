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

la = .77
lb = .87
I = .50
m = .61
mp = 30.
g = 9.8
h0 = 2.0

def traj(sol, xp, yp, index_exit):
    c_throw = sol[index_exit]
    x0 = lb*sin(c_throw[1]) + la*sin(c_throw[0])
    y0 = -lb*cos(c_throw[1]) - la*cos(c_throw[0])
    t_max = (yp[index_exit]+np.sqrt(yp[index_exit]*yp[index_exit]+2*g*y0))/g
    print(t_max)
    t = np.linspace(0,t_max, (len(sol)-index_exit))
    x = x0+xp[index_exit]*t
    y = y0 + yp[index_exit]*t - 0.5*g*t*t   
    return x,y
    
def traject(sol, xp, yp, time, index_exit):
    c_throw = sol[index_exit]
    x0 = lb*sin(c_throw[1]) + la*sin(c_throw[0])
    y0 = -lb*cos(c_throw[1]) - la*cos(c_throw[0])
    xs = []
    ys = []
    xlast = -1
    for i,t in enumerate(time):
        x = xp[i]*t + x0
        y = -0.5*g*t*t + yp[i]*t + y0   
        if y < 0:
            if xlast == -1:
                xlast = xs[-1]
            xs.append(xlast)
            ys.append(0)
        else:
            xs.append(x)
            ys.append(y)
    return xs, ys
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
        v =  np.sqrt( la*la * sol[i,2]* sol[i,2] +lb*lb*sol[i,3]*sol[i,3] + 2*la*lb*sol[i,2]*sol[i,3]*np.cos( sol[i,0]-sol[i,1] ) )
        vel.append(v )
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

c_ini = [pi/2, pi, 0., 0.]
t = np.linspace(0., 1.5, 201)

sol = odeint(f, c_ini, t)

xp = la*sol[:,2]*cos(sol[:,0]) + lb*sol[:,3]*cos(sol[:,1])
yp = la*sol[:,2]*sin(sol[:,0]) + lb*sol[:,3]*sin(sol[:,1])

ang = np.arctan2(yp, xp)
index_exit = (np.abs(ang-pi/4)).argmin()

xpr , ypr = traj(sol, xp, yp, index_exit)
print(index_exit, np.degrees(ang[index_exit]))    
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1) 
ax1.set_aspect('equal')
ax1.set_xlim(-2,2)
ax1.set_ylim(-2,2)
theta = np.linspace(0,2*pi,361)
def animate(i):

    xas = la*sin(sol[i,0])
    yas = -la*cos(sol[i,0])
    xbs = lb*sin(sol[i,1]) + xas
    ybs = -lb*cos(sol[i,1]) + yas
    thisx = [0, xas, xbs]
    thisy = [0, yas, ybs]
    xcs = la*sin(theta)
    ycs = la*cos(theta)
    xds = 4
    yds = h0+la*(sol[i,0]-sol[0,0]) 
    
    ax1.set_xlim(-2,40)
    ax1.set_ylim(-2,10)
    if i >= index_exit:
        ax1.plot(xpr[:i-index_exit], ypr[:i-index_exit], 'ko-',label='{}'.format(i/100.))
#        ax1.legend()
    else:
        ax1.clear()   
        ax1.set_xlim(-2,40)
        ax1.set_ylim(-2,10)
        
        ax1.plot(xcs, ycs, 'g-')
        ax1.plot(xds, yds, 'ko')
        ax1.plot(thisx,thisy, 'bo-', label='{}'.format(i/100.))
        ax1.legend()
    
ani = animation.FuncAnimation(fig,animate,np.arange(1,len(sol)),interval=50)
plt.show()
fig = plt.figure()
### V_B
ax2 = fig.add_subplot(211)
ax3 = fig.add_subplot(212)
#ax4 = fig.add_subplot(313)
ax2.plot(t, v(sol,t),'b-')
#ax3.plot(t, np.degrees((sol[:,1]%(2*pi))),'r-')
#ax4.plot(t, np.degrees((sol[:,0])%(2*pi)))



ax3.plot(t, ang, 'bo-')
ax3.plot(t, [pi/4]*len(ang), 'r-')

plt.show()
resp = input("Guardar [S/N]")
if resp == "S":
    name = input("Nombre del archivo ")
    ani.save('{}.mp4'.format(name))
