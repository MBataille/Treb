# -*- coding: utf-8 -*-
"""
Created on Sat May 21 08:59:06 2016

@author: klaus
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
# movimiento de un pendulo amarrado a un disco que gira.

# la = radio del disco que gira
# lb = largo del pendulo amarrado al disco 
# m = masa del objeto del pendulo
# I momento de inercia del disco que gira
# g = gravedad
# y1 = angulo de rotacion del disco
# y2 = angulo de rotacion del pendulo

###############################
# ingreso de valores de los parametros del problema
la = 1.0
lb = 1.0
I = 5.0
m = 1.0
g = 10.0
param = [la, lb, I, m, g]
################################

def b1(y):
    c = -la*lb*np.sin(y[0]-y[1])*y[3]*y[3] + g*la*np.sin(y[0])
    return c

def b2(y):
    c = la*lb*np.sin(y[0]-y[1])*y[2]*y[2] + g*lb*np.sin(y[1])
    return c

def deter(y):
    c = ( I/m + la*la*np.sin(y[0]-y[1])*np.sin(y[0]-y[1]) ) * lb*lb
    return c

def f3(y):
    c = lb*lb * b1(y) - la*lb*np.cos(y[0]-y[1]) * b2(y)
    c = c / deter(y)
    return c

def f4(y):
    c = - la*lb*np.cos(y[0]-y[1]) * b1(y) + (I/m+la*la)*b2(y)
    c = c /deter(y)
    return c

# condicion inicial, supongamos parte en forma horizontal
alfa = np.pi/2 
beta = np.pi/2
gama = 0.0
eta  = 0.0
y0 = np.array( [ alfa , beta , gama , eta ] )
t = np.linspace(0, 4, 201)

def pend(y,t):
    dydt = [ y[2] , y[3] , f3(y), f4(y) ]
    return dydt
    
sol = odeint(pend, y0, t)

xb = la*np.sin(sol[:,0]) + lb * np.sin( sol[:,1])
zb = la * np.cos(sol[:,0]) + lb * np.cos( sol[:,1])

xa = la*np.sin(sol[:,0])
za = la * np.cos(sol[:,0])

fig = plt.figure()
theta = np.linspace(0.0,2*np.pi,361)
xc = la*np.cos(theta) ; zc = la*np.sin(theta)
plt.axes().set_aspect('equal')
plt.plot(xc,zc, 'c-')
plt.plot(xb,zb, 'ro')
plt.plot(xa,za, 'go')
#fig = plt.figure()
for i in range(len(t)):
    plt.plot([xa[i],xb[i]],[za[i],zb[i]], 'b-')
    
plt.show()

 