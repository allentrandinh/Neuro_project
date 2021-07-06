
# coding: utf-8

# In[15]:

#Without this it complains about a missing display server
import matplotlib
matplotlib.use('Agg')

import sys

from mpmath import *
from matplotlib.pyplot import *
from numpy import double, array, linspace, vectorize, meshgrid


# In[16]:

#If this is commented out, it uses the default precision of 15, see
#https://docs.sympy.org/0.6.7/modules/mpmath/basics.html

#mp.dps = 20#25


# In[3]:

uR = float(sys.argv[1])
# diffusion constant [0, 10]
DR = float(sys.argv[2])

# half-life [0, 24] hours
LifTR = float(sys.argv[5]) * 60 * 60

#translation rate [0,10]
xrange = float(sys.argv[3])


#translation rate [0,10]
tmax = float(sys.argv[4])


kR = log(2.) / LifTR
muR = uR / (2. * DR)
#etaR = kR + uR^2 / (4 * DR);
etaR = -kR - uR**2./(4*DR)
lambdaR = (sqrt(uR**2. + 4. * kR * DR) - uR) / (2.* DR)
dR = kR / lambdaR
#xmax=1/abs(lambdaR)*3
#tmax=LifTR*5

# In[6]:

t = double(linspace(0.1, tmax, 15))*60*60
x = double(linspace(0.1, xrange, 15)) #35
betaR=1;

# In[7]:

def mRNAMathem(y,X):
    res = betaR/(2.*dR)*  exp(muR*X)*(-exp(-dR*y)*(exp(sqrt((-etaR - dR)/DR)*X)* erfc(sqrt(-dR - etaR)*sqrt(y) + X/(2.*sqrt(DR*y)))+exp(-sqrt((-etaR - dR)/DR)*X)*erfc(-sqrt(-dR - etaR)*sqrt(y) + X/(2.*sqrt(DR*y))))     + (exp( sqrt((-etaR)/DR)*X)* erfc(sqrt(-etaR)*sqrt(y) + X/(2.*sqrt(DR*y))) + exp(-sqrt((-etaR)/DR)*X)*   erfc(-sqrt(-etaR)*sqrt(y) + X/(2.*sqrt(DR*y)))))
    return float(res.real)


# In[11]:

vfunc = vectorize(mRNAMathem)


# In[12]:

[X, T] = meshgrid(x, t)


# In[13]:

Z = vfunc(T,X)/mRNAMathem(tmax*3600*24*10000,0)


# In[128]:

#mRNA = array([[float(mRNAMathem(y,X).real) for y in t] for X in x])
#Z = mRNA


# In[21]:

#I changed the width from 7.5 to 9.5


figure(figsize=(9.5, 6))
cp = contourf(X, T/(3600), Z, 30, cmap = cm.Blues) #cm.RdBu_r
cbar = colorbar(cp)
xlabel('distance')
ylabel('time [hours]')
#Save the figure to stdout

savefig(sys.stdout.buffer)
#savefig('ProteinTimeEvolution.png')



# In[ ]:




