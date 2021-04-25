import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')
import matplotlib.animation as animation
import seaborn as sns

import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate
def f1(r,v,theta,A_v,t):
    return v
def f2(r,v,theta,A_v,t):
    G=6.67408e-11
    M_s=5*1.989e30                        #YOU CAN CHANGE FACTORE OF MASS TO GET DIFFERENT RESULTS.
    return r*A_v**2-G*M_s/r**2
def f3(r,v,theta,A_v,t):
    return A_v
def f4(r,v,theta,A_v,t):
    return -2*v*A_v/r

def RK4(f1,f2,f3,f4,r,v,theta,A_v,ti,tf,n):
    dt=(tf-ti)/n
    a=np.zeros(shape=(n+1))
    b=np.zeros(shape=(n+1))
    c=np.zeros(shape=(n+1))
    d=np.zeros(shape=(n+1))
    a[0]=r
    b[0]=v
    c[0]=theta
    d[0]=A_v
    t=np.linspace(ti,tf,n+1)
    for i in range(n):
        k11=dt*f1(r,v,theta,A_v,t)
        k21=dt*f2(r,v,theta,A_v,t)
        k31 = dt * f3(r, v, theta, A_v, t)
        k41 = dt * f4(r, v, theta, A_v, t)
        k12 = dt * f1(r+0.5*k11, v+0.5*k21, theta+0.5*k31, A_v+0.5*k41, t+0.5*dt)
        k22 = dt * f2(r+0.5*k11, v+0.5*k21, theta+0.5*k31, A_v+0.5*k41, t+0.5*dt)
        k32 = dt * f3(r + 0.5 * k11, v + 0.5 * k21, theta + 0.5 * k31, A_v + 0.5 * k41, t + 0.5 * dt)
        k42 = dt * f4(r + 0.5 * k11, v + 0.5 * k21, theta + 0.5 * k31, A_v + 0.5 * k41, t + 0.5 * dt)
        k13 = dt * f1(r + 0.5 * k12, v + 0.5 * k22, theta + 0.5 * k32, A_v + 0.5 * k42, t + 0.5 * dt)
        k23 = dt * f2(r + 0.5 * k12, v + 0.5 * k22, theta + 0.5 * k32, A_v + 0.5 * k42, t + 0.5 * dt)
        k33 = dt * f3(r + 0.5 * k12, v + 0.5 * k22, theta + 0.5 * k32, A_v + 0.5 * k42, t + 0.5 * dt)
        k43 = dt * f4(r + 0.5 * k12, v + 0.5 * k22, theta + 0.5 * k32, A_v + 0.5 * k42, t + 0.5 * dt)
        k14 = dt * f1(r + k13, v + k23, theta + k33, A_v + k43, t + dt)
        k24 = dt * f2(r + k13, v + k23, theta + k33, A_v + k43, t + dt)
        k34 = dt * f3(r + k13, v + k23, theta + k33, A_v + k43, t + dt)
        k44 = dt * f4(r + k13, v + k23, theta + k33, A_v + k43, t + dt)
        dr=(k11+2*k12+2*k13+k14)/6
        dv = (k21 + 2 * k22 + 2 * k23 + k24) / 6
        dtheta = (k31 + 2 * k32 + 2 * k33 + k34) / 6
        dA_v = (k41 + 2 * k42 + 2 * k43 + k44) / 6
        r=r+dr
        v=v+dv
        theta=theta+dtheta
        A_v=A_v+dA_v
        a[i+1] = r
        b[i+1] = v
        c[i+1] = theta
        d[i+1] = A_v
    return t,a,b,c,d


r0=1.496e11
v0=0.00
theta0=np.pi/6
A_v0=1.990986e-7

t,a,b,c,d = RK4(f1,f2,f3,f4,r0,v0,theta0,A_v0,0,3.154e7,10000)
n=10000
x=np.zeros(shape=(n+1))
y=np.zeros(shape=(n+1))
for i in range(n+1):
    x[i]=a[i]*np.cos(c[i])
    y[i]=a[i]*np.sin(c[i])

circle=plt.Circle((0,0),696340e4,color='orange',fill=True)
circle1=plt.Circle((1.347e11,6.80e10),4*12742e5,color='b',fill=True)
fig, ax=plt.subplots()
ax.add_artist(circle)
ax.add_artist(circle1)
plt.xlim(-2e11,2e11)
plt.ylim(-2e11,2e11)
ax.plot(x,y,color='b')
plt.show()