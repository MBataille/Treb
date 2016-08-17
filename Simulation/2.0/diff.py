"""
Created on Wed May 25 21:32:37 2016

@author: martin

Simulacion del disco con peso, Muestra rotacion y trayectoria formato mp4
"""

import numpy as np
from numpy import sin, cos, pi
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def traj(sol, xp, yp, index_exit, params, h = None):
    """Trayectoria del proyectil usando velocidad
    inicial (xp,yp) y posicion inicial (angulos en sol o la
    altura (opcional))"""
    
    la, lb, I, m, mp, g, h0 = params    
    c_throw = sol[index_exit]
    x0 = lb*sin(c_throw[1]) + la*sin(c_throw[0])
    
    if h is None:
        y0 = -lb*cos(c_throw[1]) - la*cos(c_throw[0])
    else:
        y0 = h
    t_max = (yp[index_exit]+np.sqrt(yp[index_exit]*yp[index_exit]+2*g*y0))/g
    t = np.linspace(0,t_max, (len(sol)-index_exit))
    x = x0+xp[index_exit]*t
    y = y0 + yp[index_exit]*t - 0.5*g*t*t  
    return x,y

def abs(ar):
    """Valor absoluto de una lista"""
    x = []
    for a in ar:
        if a>0:
            x.append(a)
        else:
            x.append(-a)
    return x

def v(sol, t, args):
    """Velocidad tangencial en funcion de los angulos y sus derivadas en el tiempo"""
    la, lb, I, m, mp, g, h0 = args
    vel = []
    for i in range(len(t)):
        v =  np.sqrt( la*la * sol[i,2]* sol[i,2] +lb*lb*sol[i,3]*sol[i,3] + 2*la*lb*sol[i,2]*sol[i,3]*np.cos( sol[i,0]-sol[i,1] ) )
        vel.append(v )
    return vel

def f(x, t, args):
    """Resuelve ecuacion diferencial"""
  
    alpha, beta, alphap, betap = x # alpha p = alpha punto
    la, lb, I, m, mp, g, h0 = args
    
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

def f2(x, t, args):
    """Resuelve ecuacion diferencial para el momento de inercia"""
    la, lb, I, m, mp, g, h0 = args
    alpha, alphap = x
    alphapp = (-g*la*sin(alpha))/((I/m)+la*la)
    dxdt = [alphap,alphapp]
    return dxdt

def anim(i, args):
    """Plot dinamico del proyectil en el tiempo"""
    ax1, sol, index_exit, xpr, ypr, la, lb, I, m, mp, g, h0, frames = args
    # Brazo del disco (o corazon)
    xas = la*sin(sol[i,0])
    yas = -la*cos(sol[i,0])
    
    # Posicion del proyectil (cuando rota)
    xbs = lb*sin(sol[i,1]) + xas
    ybs = -lb*cos(sol[i,1]) + yas
    
    # Junta las coordenadas para plotear lineas
    thisx = [0, xas, xbs]
    thisy = [0, yas, ybs]
    
    # Circulo
    theta = np.linspace(0, 360, 361)
    xcs = la*sin(theta)
    ycs = la*cos(theta)
    
    # Contrapeso que cae    
    xds = 4
    yds = h0+la*(sol[i,0]-sol[0,0]) 

    # Time to plot!    
    ax1.clear()  
    ax1.set_xlim(-2,40)
    ax1.set_ylim(-2,10)   
    
    ax1.plot(xcs, ycs, 'g-')
    ax1.plot(xds, yds, 'ko')    
    # Indice en sol en el que el proyectil es disparado
    if i >= index_exit: # Plot de la trayectoria del proyectil
        ax1.plot(xpr[:i-index_exit], ypr[:i-index_exit], 'ko-')

    ax1.plot(thisx,thisy, 'bo-', label='{:04.3f}'.format(float(i)/(frames-1)))
    ax1.legend()

def solution(conf_ini = None, time = None, params = None):
    """ Resuelve ecuacion diferencial con los parametros dados
    (retorna una matriz solucion --sol-- y los parametros)"""
    if conf_ini is None:
        c_ini = [pi/2, pi, 0., 0.] # En esta configuracion, la energia es maxima
    else:
        c_ini = conf_ini

    if time is None:
        frames = 201
        t = np.linspace(0., 1.5, frames) # si el tiempo es +1,5s, 
                                        # hace mas de una rotacion
    else:
        t = time
        frames = len(time)
        
    if params is None:
        args = (.77, .87, .05, .65, 10., 9.8, 2.0)
                #  la,  lb,  I,   m,  mp,   g,   h0
    else:
        args = params + (9.8, 2.0) # params + g,h0
    sol = odeint(f, c_ini, t, (args,)) # Resuelve la ecuacion dif
    return sol, args
    
def inertiaSol(inertia, conf_ini = None, time = None ):
    """Resuelve ec dif para el momento de inercia con parametros dados"""
    if conf_ini is None:
        c_ini = [pi/2, 0.]
    else:
        c_ini = conf_ini

    if time is None:
        frames = 1001
        t = np.linspace(0., 20., frames) 
    else:
        t = time
        frames = len(time)   
    args1 = (.77, .0, inertia)
    args2 = (.65, 10., 9.8, 2.0) 
    sol = odeint(f2, c_ini, t, (args1+args2,)) # Resuelve la ecuacion dif
    return sol
def createAnimation(filename = None):
    """ Produce una animacion en mp4"""
    sol, params = solution() 
    la, lb, I, m, mp, g, h0 = params    

    # Cordenadas del proyectil (cuando esta rotando)
    xp = la*sol[:,2]*cos(sol[:,0]) + lb*sol[:,3]*cos(sol[:,1])
    yp = la*sol[:,2]*sin(sol[:,0]) + lb*sol[:,3]*sin(sol[:,1])
    
    # Determinar el angulo de salida
    ang = np.arctan2(yp, xp)
    index_exit = (np.abs(ang-pi/4)).argmin()
    
    # Trayectoria
    xpr , ypr = traj(sol, xp, yp, index_exit, params)
    
    # Plot!   
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1) 
    ax1.set_aspect('equal')
    params = (ax1, sol, index_exit, xpr, ypr) + params + (201,) # 201 = frames
    # Animate!
    ani = animation.FuncAnimation(fig,anim,np.arange(1,len(sol)),interval=50, fargs = (params,))
    plt.show()
    fig = plt.figure()
    
    if filename is not None:
        ani.save('{}.mp4'.format(filename))
