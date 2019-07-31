#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import as needed
import mesa_reader as mr
import matplotlib.pylab as plt
import numpy as np
import os
from math import log
from scipy.integrate import cumtrapz
import math

# import function from another file
from ipynb.fs.full.functions import getMaxRadiusProfile
from ipynb.fs.full.functions import plotTInspiral
from ipynb.fs.full.functions import plotTConv

G = 6.67408e-11 # gravitational constant
# change G to cgs units
G = G * 1e3


# In[3]:


directory = '/Users/emilyflynn/Desktop/25M_tests/15M_logs'
p = getMaxRadiusProfile(directory)

m2s = [0.005, 0.1, 0.5, 1, 5, 7.5]

plotTConv(p, True, '15M')

for i in range(len(m2s)):
    s = str(m2s[i]) + 'M'
    plotTInspiral(p, m2s[i], s)

plt.xlim(1e9, 1e14)

plt.xlabel('Radius [cm]')
plt.ylabel('Time [s]')

plt.legend()


# In[ ]:


directory = '/Users/emilyflynn/Desktop/25M_tests/30M_logs'
p = getMaxRadiusProfile(directory)

m2s = [0.05, 0.1, 0.5, 1, 5, 10, 15]

plotTConv(p, True, '30M')

for i in range(len(m2s)):
    s = str(m2s[i]) + 'M'
    plotTInspiral(p, m2s[i], s)

plt.xlim(1e9, 5e14)

plt.xlabel('Radius [cm]')
plt.ylabel('Time [s]')

plt.legend()


# In[ ]:


directory = '/Users/emilyflynn/Desktop/25M_tests/50M_logs'
p = getMaxRadiusProfile(directory)

m2s = [0.1, 0.5, 1, 5, 10, 15, 25]

plotTConv(p, True, '50M')

for i in range(len(m2s)):
    s = str(m2s[i]) + 'M'
    plotTInspiral(p, m2s[i], s)

plt.xlim(1e9, 1e14)

plt.xlabel('Radius [cm]')
plt.ylabel('Time [s]')

plt.legend()


# In[ ]:


directory = '/Users/emilyflynn/Desktop/25M_tests/70M_logs'
p = getMaxRadiusProfile(directory)

m2s = [0.1, 0.5, 1, 5, 10, 15, 25]

plotTConv(p, True, '70M')

for i in range(len(m2s)):
    s = str(m2s[i]) + 'M'
    plotTInspiral(p, m2s[i], s)

plt.xlim(1e9, 1e14)

plt.xlabel('Radius [cm]')
plt.ylabel('Time [s]')

plt.legend()


# In[2]:


directory = '/Users/emilyflynn/Desktop/1.0M_Sun'
p = getMaxRadiusProfile(directory)

m2s = [0.002, 0.005, 0.008, 0.02, 0.05, 0.08, 0.2]

plotTConv(p, True, '1M')

for i in range(len(m2s)):
    s = str(m2s[i]) + 'M'
    plotTInspiral(p, m2s[i], s)

plt.xlim(1e9, 1e14)
plt.ylim(1e3, 1e12)

plt.xlabel('Radius [cm]')
plt.ylabel('Time [s]')

plt.legend()


# In[ ]:




