import matplotlib
matplotlib.use('Agg')

import sys

from mpmath import *
from matplotlib.pyplot import *
from numpy import double, array, linspace, vectorize, meshgrid

# In[1]

#protein velocity [0,5]
uP = float(sys.argv[2])

# protein diffusion coefficient [0,10]
DP = float(sys.argv[1])

tmax = float(sys.argv[3])*60
xrange = float(sys.argv[4])

Dur = float(sys.argv[5])

x0 = float(sys.argv[6])

LifTP = float(sys.argv[7])*24*60*60

mp.dps = 20
sigma=0.50;
#x0=300;

t0=0.01;
a1=1;
a2=1/(Dur/2*60);
b1=0.1;
b2=1/(Dur/2*60);

def tf(t):
   if t>t0:
       return (t-t0)*(a1*exp(-a2*(t-t0))-b1*exp(-b2*(t-t0)))
   else:
       return 0

vtf=vectorize(tf)
timearray= double(linspace(0, max(1/a2,1/a1)*50, 100));   
Tmaxi=max(vtf(timearray))


#Tmax=9.706
#calculated using equation and given parameters


kP=log(2)/(LifTP);
muP=uP/(2*DP);
etaP=kP+uP**2/(4*DP);
lambdaP=(sqrt(uP**2+4*kP*DP)-uP)/(2*DP);
betaP=0.02301

A1=-(-etaP+a2)/(4*DP);
A2=-(-etaP+b2)/(4*DP);

def Eerf(alpha,beta,gamma):
    return exp(2*alpha*beta)*erf(alpha*gamma+beta/gamma)+exp(-2*alpha*beta)*erf(alpha*gamma-beta/gamma)


def Iminus12(a,b,A,B):
    return 1/2*sqrt(pi)/sqrt(A)*(Eerf(sqrt(A),sqrt(B),sqrt(b))-Eerf(sqrt(A),sqrt(B),sqrt(a))) 


    
def Iplus12(a,b,A,B):
    return-1./A*(sqrt(b)*exp(-A*b-B/b)-sqrt(a)*exp(-A*a-B/a))+1./(4*A)*sqrt(pi)/sqrt(A)*(Eerf(sqrt(A),sqrt(B),sqrt(b))-Eerf(sqrt(A),sqrt(B ),sqrt(a)))+1./(2*A)*sqrt(pi)*sqrt(B)*(Eerf(sqrt(B),sqrt(A),sqrt(a**(-1)))-Eerf(sqrt(B),sqrt(A),sqrt(b**(-1))))
    

    
def P(t,x):
    h=0
    if t>t0:
        h=sqrt(sigma)*exp(sigma*muP**2/2+muP*(x-x0)+sigma*kP/(4*DP))*1./16*DP**(-2)*(a1*exp(-a2*(sigma/(4*DP)+t-t0))*(Iminus12(sigma,sigma+4*DP*(t-t0),A1,(sigma*muP+2*(x-x0))**2/4)*(4*DP*(t-t0)+sigma)-Iplus12(sigma,sigma+4*DP*(t-t0),A1,(sigma*muP+2*(x-x0))**2/4))+b1*exp(-b2*sigma/(4*DP)-b2*(t-t0))*(-Iminus12(sigma,sigma+4*DP*(t-t0),A2,(sigma*muP+2*(x-x0))**2/4)*(4*DP*(t-t0)+sigma)+Iplus12(sigma,sigma+4*DP*(t-t0),A2,(sigma*muP+2*(x-x0))**2/4)));
    return float(h.real)

#t = double(linspace(0, max(1/a2,1/a1)*10, 20));
t = double(linspace(0, tmax, 20));
#x = double(linspace(max(x0-sigma*100,0), x0+sigma*100, 20));
x = double(linspace(max(x0-xrange,0), x0+xrange, 25));



vfunc=vectorize(P)

[X,T]=meshgrid(x,t)   

Z=(2.1/60)/Tmaxi*vfunc(T,X)


figure(figsize=(9.5, 6))
cp = contourf(X, T/60, Z, 10, cmap = cm.RdBu_r)
#the 9000 is to normalize the function. The 55 ensures the time axis is accurate 
cbar = colorbar(cp)
#clabel(cp, inline=True, fontsize=8)
xlabel('distance from soma (μm)')
ylabel('time [minutes]')
savefig(sys.stdout.buffer)


