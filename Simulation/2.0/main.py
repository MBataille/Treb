# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 21:45:17 2016

@author: martin
"""

import diff
import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
import scipy.fftpack

NI = 5
N = 4096
dt = 20. / N

def seno(t,p,b):
    res = b*np.cos(2*np.pi*t/p)
    return res

def period(t,alpha):
    creciente = False
    for i in range(1,len(alpha)-1):
        if not creciente and alpha[i] < alpha[i+1]:
            creciente = True
        elif creciente and alpha[i] > alpha[i+1]:
            break
    index_max = i
    dt = t[-1]/(len(t)-1.)
    return index_max*dt

def solution(I, ang = np.pi/2., t = np.linspace(0,20,N)):
    sol = diff.inertiaSol(I,time = t, conf_ini = [ang, 0.])
    return sol[:,0]

def solEstim():
    g = 9.8
    la = .77
    I = np.linspace(0,1, NI)
    m = 0.04
    k = np.sqrt((g*la)/(I/m + la*la))
    T = 2*np.pi/k
    return I,T
def estim(t,alpha):
    plt.figure()
    T = period(t,alpha)
    
    sp = np.abs( scipy.fftpack.fft(alpha) )[:N/2]
    freq = np.linspace(0.0,1.0/(2*dt),N/2)
        
    indice = np.where(sp==sp.max())
    T_fft = 1./freq[indice]
    plt.plot(t, alpha, 'k-', label = 'alpha(t)')
    plt.plot(t, np.pi/2. * np.cos(2*np.pi*t/T_fft), 'b-', label = 'fft')
    plt.plot(t, np.pi/2. * np.cos(2*np.pi*t/T), 'r-', label = 'algo')
    plt.legend()
    plt.show()

def inertia():
    I = np.linspace(0.1, 1, NI)
    ang = [np.pi/2., np.pi/4., np.pi/8.]
    for c in ang:    
        T = []
        for i in I:
            t = np.linspace(0,20,N)
            sol = diff.inertiaSol(i ,conf_ini = [c, 0.],time = t)
            alpha = sol[:,0]
            T.append(period(t,alpha))
        plt.plot(I, T, label = 'pi/{}'.format(int(np.pi/c)))
    x,y = solEstim()
    plt.plot(x,y, label = 'analytique')

def calcInertia(d, l1, l2, alpha):
    tan = np.tan(alpha)
    cos = np.cos(alpha)
    a = d * (tan + tan*tan*tan/3.) * l2*l2*l2*l2*cos*cos*cos*cos / 4.
    b = d * tan * ((l1*l1*l1*tan*tan*(l1-l2*cos)/3.) - ((l1*l1-l2*l2*cos*cos)*l1*l1*tan*tan/2.)
                    + (l1*(1+tan*tan)*(l1*l1*l1-l2*l2*l2*cos*cos*cos)/3.) 
                    - (tan*tan * (l1*l1*l1*l1-l2*l2*l2*l2*cos*cos*cos*cos))/12.)
    return a + b
    
#la = [0.2159, 0.1905, 0.127, 0.1956, 0.2083, 0.2159]
#ang = [60., 40., 34., 106., 65., 55.]
#d = 600
#print(calcInertia(d,0.2159, 0.2159, np.pi/3))
#print(calcInertia(d, 0.2159, 0.20828, np.pi*55./180.))
#print(calcInertia(d, 0.20828, 0.19558, np.pi*65./180.))
#print(calcInertia())
inertia()
plt.legend()
plt.show()
#estim(np.linspace(0, 20, N), solution(0.5))

