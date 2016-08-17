# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 18:01:25 2016

@author: martin
"""

import numpy as np

l = np.loadtxt('inercia.txt')
print(l.sum()/600.)