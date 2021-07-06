import matplotlib
matplotlib.use('Agg')

import sys

from mpmath import *
from matplotlib.pyplot import *
from numpy import double, array, linspace, vectorize, meshgrid

# In[1]


uR = float(sys.argv[1])

# diffusion constant [0, 10]
DR = float(sys.argv[2])

# half-life [0, 24] hours
LifTR = float(sys.argv[3]) * 60 * 60

SomDentr=float(sys.argv[4])

# protein diffusion coefficient [0,10]
DP = float(sys.argv[5])

#protein velocity [0,5]
uP = float(sys.argv[6])

#protein half life [0,20] days
LifTP = float(sys.argv[7])*24*60*60


#translation rate [0,10]
xrange = float(sys.argv[8])


#translation rate [0,10]
tmax = float(sys.argv[9])


muR = uR / (2. * DR)
muP=uP/(2*DP);

kR=log(2)/(LifTR);
kP=log(2)/(LifTP);

lambdaR=(sqrt(uR**2+4*kR*DR)-uR)/(2*DR);
lambdaP=(sqrt(uP**2+4*kP*DP)-uP)/(2*DP);

dR=kR*2*DR/(sqrt(uR**2+4*kR*DR)-uR);
dRP=kP*2*DP/(sqrt(uP**2+4*kP*DP)-uP);

dP=kP*(DP*lambdaR**2+uP*lambdaR-kP)/(DP*lambdaP*lambdaR+uP*lambdaP-kP);


etaP=kP+uP**2/(4*DP);
etaR = -kR - uR**2./(4*DR);

C1=1/(etaP-DP*(lambdaR+muP)**2);
C4=1/dP-C1;


#xmax=1/abs(lambdaP)*3
#tmax=LifTP*5
betaP=1;


t = double(linspace(0.0000000001, tmax, 20))*3600*24;
x = double(linspace(0, xrange, 15)) #35  ; 

def mRNAMathem(t,x):
    res = betaP/(2.*dRP)*exp(muP*x)*(-exp(-dRP*t)*(exp(sqrt((etaP-dRP)/DP)*x)* erfc(sqrt(-dRP + etaP)*sqrt(t) + x/(2.*sqrt(DP*t)))+exp(-sqrt((etaP - dRP)/DP)*x)*erfc(-sqrt(-dRP+etaP)*sqrt(t) + x/(2.*sqrt(DP*t))))+(exp(sqrt((etaP)/DP)*x)*erfc(sqrt(etaP)*sqrt(t)+x/(2.*sqrt(DP*t))) + exp(-sqrt((etaP)/DP)*x)*erfc(-sqrt(etaP)*sqrt(t) + x/(2.*sqrt(DP*t)))))
    return float(res.real)

def Pdendrite(t,x):
    b=((1/(2*dP)*exp(-lambdaP*x)*erfc(x/(2*sqrt(DP*t))-sqrt(t*etaP))+1/(kP-DP*lambdaR**2-lambdaR*uP)*(exp(-lambdaR*x)-1/2*exp(-lambdaP*x)*erfc(x/(2*sqrt(DP*t))-sqrt(t*etaP)))+C1/2*exp(-1/C1*t+2*muP*x+lambdaR*x)*erfc(x/(2*sqrt(DP*t))+(lambdaR+muP)*sqrt(t*DP))+C4/2*exp((muP+sqrt(etaP/DP))*x)*erfc(x/(2*sqrt(DP*t))+sqrt(t*etaP))+C1*(-exp(-1/C1*t-lambdaR*x)+1/2*exp(-1/C1*t-lambdaR*x)*erfc(x/(2*sqrt(DP*t))-(lambdaR+muP)*sqrt(t*DP)))-1/(2*dP)*exp(-dP*t+muP*x)*(exp(x*sqrt((etaP-dP)/DP))*erfc(x/(2*sqrt(DP*t))+sqrt(t*(etaP-dP)))+exp(-x*sqrt((etaP-dP)/DP))*erfc(x/(2*sqrt(DP*t))-sqrt(t*(etaP-dP))))))
    return float(b.real)

vfunc = vectorize(Pdendrite)
vfunc1 = vectorize(mRNAMathem)


[X, T] = meshgrid(x, t)



# In[13]:

Z = (vfunc(T,X)*SomDentr*1.19+(1-SomDentr)*0.23*vfunc1(T,X))/(Pdendrite(tmax*3600*24*10000,0)*SomDentr*1.19+(1-SomDentr)*0.23*mRNAMathem(tmax*3600*24*10000,0))

figure(figsize=(9.5, 6))
cp = contourf(X, T/(3600*24.), Z, 30, cmap = cm.Oranges)#cm.RdBu_r
cbar = colorbar(cp)
#clabel(cp, inline=True, fontsize=8)
xlabel('distance from soma (μm)')
ylabel('time [days]')
title('Protein')
savefig(sys.stdout.buffer)




# In[ ]:
