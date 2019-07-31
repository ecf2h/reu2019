#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# CONTENTS
# getClosestModel
# getMaxRadiusProfile
# plotTInspiral
# plotTConv
# plotEbinds
# plotMaxLum
# plotDragLum
# plotOrbitalEnergyChange
# getTInspiral
# getTConv
# getMaxLum
# getDragLum


# In[ ]:


# find and import the profile closest to your target profile

import os
import mesa_reader as mr

def getClosestModel(modelNumber, workingDirectory):
    """Returns the MESA profile with the closest model number in a specified directory.
    
    Keyword arguments:
    modelNumber -- target/desired model number
    workingDirectory -- target directory
    """
    print('Looking for the profile closest to model #' + str(modelNumber) + ' in ' + str(workingDirectory))
    filenames = []

    for root, dirs, files in os.walk(workingDirectory):
        for file in files:
            if file.endswith("profiles.index"):
                filenames.append(os.path.join(root, file))

    profiles = []
    # import files
    for file in filenames:
        i = mr.MesaProfileIndex(file)
        profiles.append(i)

    # find the closest model number - ugly but functional
    closest = 0
    diff = 1e10
    profilePath = ''
    j = 0
    for index in profiles:
        values = index.model_numbers
        profileNums = index.profile_numbers
        k=0
        for i in values:
            if (abs(modelNumber - i) < diff):
                diff = abs(modelNumber - i)
                closest = i
                og = filenames[j]
                og = og[:-14] # put together the file name given the directory
                profilePath = og + 'profile' + str(profileNums[k]) + '.data'
            k+=1
        j+=1

    print('Actual model number: ' + str(closest))
    print('Difference between target and actual model: ' + str(diff))
    print('File path: ' + str(profilePath))
    print('')

    # import target profile
    p = mr.MesaData(profilePath)
    return p


# In[ ]:


# find the model number with the max radius

def getMaxRadiusProfile(workingDirectory):
    """Returns the profile with the biggest radius in a directory.
    
    Keyword arguments:
    workingDirectory -- target directory
    """
    
    import os
    import mesa_reader as mr
    from heapq import nlargest
    import numpy as np

    filenames = []

    for root, dirs, files in os.walk(workingDirectory):
        for file in files:
            if file.endswith("history.data"):
                filenames.append(os.path.join(root, file))

    # for each file, go in and correct that quotation mark error
    # replace 10.14-2019 with 10.14-2019"
    for file in filenames:
        s = open(file).read()
        s = s.replace('10.14-2019 ', '10.14-2019"')
        f = open(file, 'w')
        f.write(s)
        f.close()

    # for each file, read it into a variable based on which log folder it's in
    hBRExists = False
    hARExists, hRExists = False, False
    htSBExists, htlgTExists = False, False
    hCExists, hFExists, hLExists = False, False, False

    for file in filenames:
        if 'before_remove' in file:
            hBR = mr.MesaData(file)
            hBRExists = True

        elif 'after_remove' in file:
            hAR = mr.MesaData(file)
            hARExists = True

        elif 'remove' in file:
            hR = mr.MesaData(file)
            hRExists = True

        elif 'to_si_burn' in file:
            htSB = mr.MesaData(file)
            htSBExists = True

        elif 'to_lgT_9.9' in file:
            htlgT = mr.MesaData(file)
            hlgTExists = True

        elif 'convert' in file:
            hC = mr.MesaData(file)
            hCExists = True

        elif 'finish' in file:
            hF = mr.MesaData(file)
            hFExists = True

        else:
            hL = mr.MesaData(file)
            hLExists = True

    # frankenstein the data together because MESA is MESA

    hModels = []
    hRadius = []

    if hBRExists:
        hRadius.append(hBR.log_R)
        hModels.append(hBR.model_number)

    if hRExists:
        hRadius.append(hR.log_R)
        hModels.append(hR.model_number)

    if hARExists:
        hRadius.append(hAR.log_R)
        hModels.append(hAR.model_number)

    if htSBExists:
        hRadius.append(htSB.log_R)
        hModels.append(htSB.model_number)

    if htlgTExists:
        hRadius.append(htlgT.log_R)
        hModels.append(htlgT.model_number)

    if hCExists:
        hRadius.append(hC.log_R)
        hModels.append(hC.model_number)

    # these are behaving weirdly
    # if hFExists:
    #     hRadius.append(hF.log_R)
    #     hModels.append(hF.model_number)

    if hLExists:
        hRadius.append(hL.log_R)
        hModels.append(hL.model_number)
    
    
    # find the model number for when the radius is at its maximum
    maxValues = []
    maxModelNumbers = []

    # find the biggest 5 in each history file
    for x in range(len(hRadius)):
        data = hRadius[x]
        maxValues.append(nlargest(1, data))
        index = nlargest(1, range(len(data)), key=lambda idx: data[idx])
        maxModelNumbers.append(hModels[x][index])

    print(maxValues)
    print(maxModelNumbers)
    print('-----')
    print(max(maxValues))
    modelNumber = int(maxModelNumbers[np.argmax(maxValues)])
    print(modelNumber)

    # import desired profile

    # find all profile.index files
    filenames = []

    for root, dirs, files in os.walk(workingDirectory):
        for file in files:
            if file.endswith("profiles.index"):
                filenames.append(os.path.join(root, file))

    profiles = []
    # import files
    for file in filenames:
        i = mr.MesaProfileIndex(file)
        profiles.append(i)

    return getClosestModel(modelNumber, workingDirectory)


# In[1]:


def getR2(m):
    """Return a value for the radius of the secondary given the mass.
    Determined using Wilson & Nordhaus (2019), Eq. 9.
    
    Keyword arguments:
    m -- mass
    """
    import math
    from math import log
    if (m > 0.077):
        r = m**0.92
    elif (m < 0.0026):
        r = 0.10045 # r_jupiter
    else:
        r = 0.117
        r = r - 0.054*(log(m/0.0026)**2)
        r = r + 0.024*(log(m/0.0026)**3)
    return r

def plotTInspiral(p, m2, label):
    """Plot the inspiral time scale.
    
    Keyword arguments:
    p -- MESA profile
    m2 -- mass of the secondary
    label -- text appearing in the legend
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    import os
    from math import log
    from scipy.integrate import cumtrapz
    import math

    # import function from another file
    from ipynb.fs.full.functions import getMaxRadiusProfile

    G = 6.67408e-11 # gravitational constant
    # change G to cgs units
    G = G * 1e3

    coreMass = p.he_core_mass + p.c_core_mass + p.o_core_mass + p.si_core_mass + p.fe_core_mass
    # change from Msuns to grams
    coreMass = coreMass*1.989e33

    radius = p.radius
    radius = radius*69.551e9

    # setting up constants
    r2 = getR2(m2)
    m2 = m2*1.989e33 # units
    r2 = r2*69.551e9 # units
    xi = 4
    rshred = r2 * (2*coreMass/m2)**(1/3)
    k = 4 * xi * math.pi * G * m2

    # density
    rho = p.logRho
    rho = 10**rho

    # masses
    masses = p.mass
    masses = masses*1.989e33

    # keplerian velocity
    vkep_r = np.sqrt(G * masses / radius)

    # dM/dr
    dMdr = np.diff(masses) / np.diff(radius)

    # make all the array sizes the same
    vkep_r = vkep_r[:-1]
    radius = radius[:-1]
    rho = rho[:-1]
    masses = masses[:-1]

    # integrand
    integrand = (dMdr - (masses / radius)) * vkep_r / (k * radius * rho)

    # look through everything in the radius array and compare it to rshred
    # when the value is <= rshred, save that index
    # chop the integrand and radius arrays at that index
    i = 0
    for x in radius:
        if x > rshred:
            i+=1

#     radius = radius[100:i] # cut first
#     integrand = integrand[100:i] # cut first

    # actually integrate
    tInspiral = cumtrapz(y=integrand, x=radius)
    
    radius = radius[100:i] # int first
    tInspiral = tInspiral[100:i] # int first
    
#     radius = radius[:-1] # cut first
    
    radius = np.flip(radius)

    point = plt.plot(radius[0], tInspiral[0], 'o')
    y = plt.getp(point[0], 'color')
    plt.loglog(radius, tInspiral, label=label, color=y, linestyle=':')


# In[ ]:


# if they are from the MESA site it should be (p, False, True)
# if we made them it is (p, True, False)

def plotTConv(p, nearZero, label):
    """Plot the convective time-scale.
    
    Keyword arguments:
    p -- MESA profile
    nearZero -- boolean of whether any values are near zero
    label -- text appearing in legend
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    from math import log
    from scipy.integrate import cumtrapz
    
    r = p.radius
    
    r = r*69.551e9   # change units to cm
    r = r[100:]      # cut off the first 100 points - takes care of surface weirdness

    v = p.log_conv_vel
    
    if nearZero:
        v = [i if i>1e-8 else 1e-8 for i in v] # if having issues with v being too close to 0
    
    v = np.power(10, v) # un-log it
    v = -1/v            # it'll be integrated like this
    v = v[100:]         # cut off first 100 points

    # integrate
    tconv = []
    tconv = cumtrapz(y=v, x=r)
    
    plt.loglog(r[:-1], tconv, label=label)
    return


# In[3]:


def plotEbinds(p, label):
    """Plot the binding energy.
    
    Keyword arguments:
    p -- MESA profile
    label -- text appearing in legend
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    import os
    from math import log
    from scipy.integrate import cumtrapz

    G = 6.67408e-11 # gravitational constant
    # change G to cgs units
    G = G * 1e3
    
    r = p.radius     # bring in the radius

    r = r*69.551e9   # change units to cm
    r = r[100:]      # cut off the first 100 points - takes care of surface weirdness

    m = p.mass
    m = m*1.989e33   # change units to grams
    m = m[100:]      # exclude first 100 points

    # integrate
    integrand = []

    for i in range(len(m)):
        x = G * m[i] / r[i]
        integrand.append(x)

    ebind = cumtrapz(y=integrand, x=m)
    
    plt.loglog(r[:-1], -ebind, label=label)
    return


# In[1]:


def plotMaxLum(p, localSim, label):
    """Plot the maximum luminosity.
    
    Keyword arguments:
    p -- MESA profile
    localSim -- boolean; if the profile was generated locally
    label -- text appearing in legend
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    import os
    from math import log
    import math

    G = 6.67408e-11 # gravitational constant
    # change G to cgs units
    G = G * 1e3
    
    r = p.radius     # bring in the radius

    r = r*69.551e9   # change units to cm
    r = r[100:]      # cut off the first 100 points - takes care of surface weirdness

    # need: beta ~ 5, 4pi, density[r], radius, sound speed c[r]
    # define a constant, k
    k = 5 * 4 * math.pi

    # get density
    rho = p.logRho
    rho = 10**rho
    
    # get sound speed
    if localSim:
        c = p.csound
    
    else:
        # the web files don't have csound
        # csound=sqrt(5/3*pressure[r] / rho[r])
        ks = 5 * p.pressure / 3
        c = np.sqrt(ks / rho)
        
    c = c[100:]
    rho = rho[100:]

    lumMax = k * r**2 * rho * c**3
    plt.loglog(r, lumMax, label=label)
    return


# In[1]:


def getR2(m):
    """Return a value for the radius of the secondary given the mass.
    Determined using Wilson & Nordhaus (2019), Eq. 9.
    
    Keyword arguments:
    m -- mass
    """
    from math import log
    if (m > 0.077):
        r = m**0.92
    elif (m < 0.0026):
        r = 0.10045 # r_jupiter
    else:
        r = 0.117
        r = r - 0.054*(log(m/0.0026)**2)
        r = r + 0.024*(log(m/0.0026)**3)
    return r

def plotDragLum(p, m2, label):
    """Plot the drag luminosity.
    
    Keyword arguments:
    p -- MESA profile
    m2 -- mass of the secondary
    label -- text appearing in legend
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    import os
    from math import log
    from scipy.integrate import cumtrapz
    import math

    # import function from another file
    from ipynb.fs.full.functions import getMaxRadiusProfile

    G = 6.67408e-11 # gravitational constant
    # change G to cgs units
    G = G * 1e3

    coreMass = p.he_core_mass + p.c_core_mass + p.o_core_mass + p.si_core_mass + p.fe_core_mass
    # change from Msuns to grams
    coreMass = coreMass*1.989e33
    
    # get rshred
    r2 = getR2(m2)
    m2 = m2*1.989e33 # units
    r2 = r2*69.551e9 # units
    rshred = r2 * (2*coreMass/m2)**(1/3)

    radius = p.radius
    radius = radius*69.551e9

    masses = p.mass
    masses = masses*1.989e33

    vkep_r = np.sqrt(G * masses / radius)

    rho = p.logRho
    rho = 10**rho

    xi = 4

    rAcc = 2 * G * m2 / vkep_r**2

    dragLum = xi * math.pi * rAcc**2 * rho * vkep_r**3
    
    # stop at rshred
    i = 0
    for x in radius:
        if x > rshred:
            i+=1

    radius = radius[100:i]
    dragLum = dragLum[100:i]

    point = plt.plot(radius[len(radius)-1], dragLum[len(dragLum)-1], 'o')
    y = plt.getp(point[0], 'color')
    plt.loglog(radius, dragLum, label=label, color=y)


# In[ ]:


def plotOrbitalEnergyChange(p, m2, label):
    """Plot the change in orbital energy.
    
    Keyword arguments:
    p -- MESA profile
    m2 -- mass of the secondary
    label -- text appearing in legend
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    import os
    from math import log
    from scipy.integrate import cumtrapz
    import math

    # import function from another file
    from ipynb.fs.full.functions import getMaxRadiusProfile

    G = 6.67408e-11 # gravitational constant
    # change G to cgs units
    G = G * 1e3
    
    coreMass = p.he_core_mass + p.c_core_mass + p.o_core_mass + p.si_core_mass + p.fe_core_mass
    # change from Msuns to grams
    coreMass = coreMass*1.989e33
    
    # get rshred
    r2 = getR2(m2)
    m2 = m2*1.989e33 # units
    r2 = r2*69.551e9 # units
    rshred = r2 * (2*coreMass/m2)**(1/3)

    radius = p.radius
    radius = radius*69.551e9

    masses = p.mass
    masses = masses*1.989e33

    ri = p.radius[0]*69.551e9
    mi = p.initial_mass*1.989e33

    deltaEOrb = (mi / ri) - (masses / radius)
    deltaEOrb = deltaEOrb * G * m2 / 2
    deltaEOrb = abs(deltaEOrb)
    
    # look through everything in the radius array and compare it to rshred
    # when the value is <= rshred, save that index
    # chop the integrand and radius arrays at that index
    i = 0
    for x in radius:
        if x > rshred:
            i+=1

    radius = radius[100:i]
    deltaEOrb = deltaEOrb[100:i]

    point = plt.plot(radius[len(radius)-1], deltaEOrb[len(deltaEOrb)-1], 'o')
    y = plt.getp(point[0], 'color')
    plt.loglog(radius, deltaEOrb, label=label, color=y)


# In[4]:


def getR2(m):
    """Return a value for the radius of the secondary given the mass.
    Determined using Wilson & Nordhaus (2019), Eq. 9.
    
    Keyword arguments:
    m -- mass
    """
    import math
    from math import log
    if (m > 0.077):
        r = m**0.92
    elif (m < 0.0026):
        r = 0.10045 # r_jupiter
    else:
        r = 0.117
        r = r - 0.054*(log(m/0.0026)**2)
        r = r + 0.024*(log(m/0.0026)**3)
    return r

def getTInspiral(p, m2):
    """Return an array containing the inspiral time-scale values.
    The outermost 100 points are removed.
    
    Keyword arguments:
    p -- MESA profile
    m2 -- mass of the secondary
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    import os
    from math import log
    from scipy.integrate import cumtrapz
    import math

    # import function from another file
    from ipynb.fs.full.functions import getMaxRadiusProfile

    G = 6.67408e-11 # gravitational constant
    # change G to cgs units
    G = G * 1e3

    coreMass = p.he_core_mass + p.c_core_mass + p.o_core_mass + p.si_core_mass + p.fe_core_mass
    # change from Msuns to grams
    coreMass = coreMass*1.989e33

    radius = p.radius
    radius = radius*69.551e9

    # setting up constants
    r2 = getR2(m2)
    m2 = m2*1.989e33 # units
    r2 = r2*69.551e9 # units
    xi = 4
    rshred = r2 * (2*coreMass/m2)**(1/3)
    # rshred = r2 * math.pow(2*coreMass/m2, 0.33333333)
    k = 4 * xi * math.pi * G * m2

    # density
    rho = p.logRho
    rho = 10**rho

    # masses
    masses = p.mass
    masses = masses*1.989e33

    # keplerian velocity
    vkep_r = np.sqrt(G * masses / radius)

    # dM/dr
    dMdr = np.diff(masses) / np.diff(radius)

    # make all the array sizes the same
    vkep_r = vkep_r[:-1]
    radius = radius[:-1]
    rho = rho[:-1]
    masses = masses[:-1]

    # integrand
    integrand = (dMdr - (masses / radius)) * vkep_r / (k * radius * rho)

    radius = radius[100:]
    integrand = integrand[100:]
    
    # actually integrate
    tInspiral = cumtrapz(y=integrand, x=radius)
    
    return np.flip(tInspiral)


# In[ ]:


def getTConv(p):
    """Return an array containing the convective time-scale values.
    The outermost 100 points are removed.
    
    Keyword arguments:
    p -- MESA profile
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    from math import log
    from scipy.integrate import cumtrapz
    
    r = p.radius     # bring in the radius
    
    r = r*69.551e9   # change units to cm
    r = r[100:]

    v = p.log_conv_vel
    
    v = [i if i>1e-8 else 1e-8 for i in v] # if having issues with v being too close to 0
    
    v = np.power(10, v) # un-log it
    v = -1/v            # it'll be integrated like this
    v = v[100:]

    # integrate
    tconv = []
    tconv = cumtrapz(y=v, x=r)

    return tconv[:-1]


# In[ ]:


def getMaxLum(p, localSim):
    """Return an array containing the maximum luminosity values.
    The outermost 100 points are removed.
    
    Keyword arguments:
    p -- MESA profile
    localSim -- boolean; if the profile was generated locally
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    import os
    from math import log
    import math

    G = 6.67408e-11 # gravitational constant
    # change G to cgs units
    G = G * 1e3
    
    if localSim:
        r = p.logR       # bring in the radius
        r = 10**r        # un-log it

    else:
        r = p.radius     # bring in the radius

    r = r*69.551e9   # change units to cm

    # need: beta ~ 5, 4pi, density[r], radius, sound speed c[r]
    # define a constant, k
    k = 5 * 4 * math.pi

    # get density
    rho = p.logRho
    rho = 10**rho
    
    # get sound speed
    if localSim:
        c = p.csound
    
    else:
        # the web files don't have csound
        ks = 5 * p.pressure / 3
        c = np.sqrt(ks / rho)

    lumMax = k * r**2 * rho * c**3
    return lumMax[100:-2]


# In[ ]:


def getR2(m):
    """Return a value for the radius of the secondary given the mass.
    Determined using Wilson & Nordhaus (2019), Eq. 9.
    
    Keyword arguments:
    m -- mass
    """
    from math import log
    if (m > 0.077):
        r = m**0.92
    elif (m < 0.0026):
        r = 0.10045 # r_jupiter
    else:
        r = 0.117
        r = r - 0.054*(log(m/0.0026)**2)
        r = r + 0.024*(log(m/0.0026)**3)
    return r

def getDragLum(p, m2):
    """Return an array containing the drag luminosity values.
    The outermost 100 points are removed.
    
    Keyword arguments:
    p -- MESA profile
    m2 -- mass of the secondary
    """
    import mesa_reader as mr
    import matplotlib.pylab as plt
    import numpy as np
    import os
    from math import log
    from scipy.integrate import cumtrapz
    import math

    # import function from another file
    from ipynb.fs.full.functions import getMaxRadiusProfile

    G = 6.67408e-11 # gravitational constant
    # change G to cgs units
    G = G * 1e3

    coreMass = p.he_core_mass + p.c_core_mass + p.o_core_mass + p.si_core_mass + p.fe_core_mass
    # change from Msuns to grams
    coreMass = coreMass*1.989e33
    
    # get rshred
    r2 = getR2(m2)
    m2 = m2*1.989e33 # units
    r2 = r2*69.551e9 # units
    rshred = r2 * (2*coreMass/m2)**(1/3)

    radius = p.radius
    radius = radius*69.551e9

    masses = p.mass
    masses = masses*1.989e33

    vkep_r = np.sqrt(G * masses / radius)

    rho = p.logRho
    rho = 10**rho

    xi = 4

    rAcc = 2 * G * m2 / vkep_r**2

    dragLum = xi * math.pi * rAcc**2 * rho * vkep_r**3

    return dragLum[100:-2]

