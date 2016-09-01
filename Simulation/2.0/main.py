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

NI = 5000
N = 2048
dt = 20. / N

def period(t,alpha):
    creciente = False
    for i in range(1,len(alpha)-1):
        if not creciente and alpha[i] < alpha[i+1]:
            creciente = True
        elif creciente and alpha[i] > alpha[i+1]:
            break
    index_max = i
    return t[index_max]

def solution(I, ang = np.pi/2., t = np.linspace(0,20,N)):
    sol = diff.inertiaSol(I,time = t, conf_ini = [ang, 0.])
    return sol[:,0]

def solEstim():
    g = 9.8
    la = .77
    I = np.linspace(0,1, NI)
    m = 0.65
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
    I = np.linspace(0, 1, NI)
    ang = [30*np.pi/180., 20*np.pi/180.]
    for c in ang:    
        T = []
        for i in I:
            t = np.linspace(0,5,N)
            sol = diff.inertiaSol(i ,conf_ini = [c, 0.],time = t)
            alpha = sol[:,0]
            T.append(period(t,alpha))
        plt.plot(I, T, label = '{}°'.format(int(round(180.*c/np.pi))))
    x,y = solEstim()
    plt.plot(x,y, label = 'analytique')
    plt.plot([0,1],[1.91,1.91], 'r-', label = 'T')
    plt.xticks(np.arange(I[0], I[-1], 0.1))

def calcInertia(d, l1, l2, alpha):
    tan = np.tan(alpha)
    cos = np.cos(alpha)
    a = d * (tan + tan*tan*tan/3.) * l2*l2*l2*l2*cos*cos*cos*cos / 4.
    b = d * tan * ((l1*l1*l1*tan*tan*(l1-l2*cos)/3.) - ((l1*l1-l2*l2*cos*cos)*l1*l1*tan*tan/2.)
                    + (l1*(1+tan*tan)*(l1*l1*l1-l2*l2*l2*cos*cos*cos)/3.) 
                    - (tan*tan * (l1*l1*l1*l1-l2*l2*l2*l2*cos*cos*cos*cos))/12.)
    return a + b

def recInertia(l):
    d = 540
    a = 0.0254
    b = 2 * a
    s = 0.
    for c in l:
        m = a * b * c * d    
        s += m*(2*a*a + 2*b*b + 2*c*c)/12.
    return s
#t = np.linspace(0, 20, N)
#estim(t, solution(3))
    
#inertia()
#plt.legend()
#plt.show()
l = [0.63, 0.41 ]
print(recInertia(l))
#s = 0
#n = int(input('n '))
#for i in range(n):
#    l1 = float(input('l1 '))
#    l2 = float(input('l2 '))
#    ang = float(input('ang ')) * np.pi/180.
#    s+=calcInertia(600, l1, l2, ang)
#    print(calcInertia(600, l1, l2, ang),)
#print(s)
#estim(np.linspace(0, 20, N), solution(0.5))

