import math
import matplotlib.pyplot as plt
import numpy as np

G= 6.67428*(10**(-11))
M=1.9885*(10**(30))
r=1496*(10**(8))
vo=1.0
xo=1496*(10**(8))
a= -G*M*(xo**(2))/((r**(3))*(vo**(2)))

fxi=[]
fyi=[]
vxi=[]
vyi=[]
xi=[]
yi=[]

def fx(x,y):
    return a*x
def fy(x,y):
    return a*y


t0=0
tf=365*24*60*60/(xo/vo)
x0=r/xo
xi.append(x0)
y0= 0
yi.append(y0)
vx0=0 
vxi.append(vx0)
vy0= (np.sqrt(G * M / r))/vo   
vyi.append(vy0)
h=(2*24*60*60)/(xo/vo)
t=np.arange(t0,tf,h)

for i in range(len(t)-1):
    vxi.append(vx0+fx(x0, y0)*h)
    vyi.append(vy0+fy(x0, y0)*h)
    xi.append(x0+vxi[i]*h)
    yi.append(y0+vyi[i]*h)
    x0=xi[i+1]
    y0=yi[i+1]
    vx0=vxi[i+1]
    vy0=vyi[i+1]

plt.figure(figsize=(6,6))
plt.plot(xi,yi)
