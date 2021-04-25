import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

def RK4(f1,ti,tf,y0,n):
    t=np.zeros(n)
    h = (tf - ti) / (n - 1)
    t[0] = ti
    y=np.zeros(shape=(len(y0),n))
    for i in range(len(y0)):
        y[i][0]=y0[i]


    for i in range(n-1):
        a=np.zeros(len(y0))
        b = np.zeros(len(y0))
        c = np.zeros(len(y0))
        e = np.zeros(len(y0))
        for k in range(len(y0)):
            e[k] = y[k][i]
        k1=f1(t[i],e)
        for k in range(len(y0)):
            a[k]=y[k][i]+k1[k]*h/2
        k2=h*f1(t[i] + h/2,a)
        for k in range(len(y0)):
            b[k]=y[k][i]+k2[k]*h / 2
        k3=h*f1(t[i]+h/2,b)
        for k in range(len(y0)):
            c[k]=y[k][i] + k3[k] * h
        k4=h*f1(t[i]+h,c)
        for j in range(len(y0)):
            y[j][i+1]=y[j][i]+((k1[j]+2*(k2[j]+k3[j])+k4[j])*h)/6
        t[i+1]=t[i]+h
    return t,y



def f1(t,y):
    R1=np.array([y[0],y[1],y[2]])
    R2 = np.array([y[3], y[4], y[5]])
    V1 = np.array([y[6], y[7], y[8]])
    V2 = np.array([y[9], y[10], y[11]])
    r=np.linalg.norm(R2-R1)
    A1=G*m2*(R2-R1)/r**3
    A2=G*m1*(R1-R2)/r**3
    dydt=np.array([V1,V2,A1,A2]).reshape(12,)
    return dydt


G=6.67259e-20
m1=1.e26
m2=1.e26
t0=0
tf=480
au=1.496e8
r1=np.array([0.,0.,0.])
r2=np.array([3000,0,0])
v1=np.array([10,20,30]) #230
v2=np.array([0,40,0])


b0=np.array([r1,r2,v1,v2])
y0=b0.reshape(12,)

'''
# For Better Accuracy Use This Block
#t1=np.linspace(t0,tf,2000)
#sol=scipy.integrate.solve_ivp(f1,[t0,tf],y0,method='RK45',t_eval=t1)
#t=sol.t
#y=sol.y
'''


t,y=RK4(f1,t0,tf,y0,1000)
n=len(t)
x=[]
y1=[]
z=[]
for i in range(n):
    x.append((m1*y[0][i]+m2*y[3][i])/(m1+m2))
    y1.append((m1 * y[1][i] + m2 * y[4][i]) / (m1 + m2))
    z.append((m1 * y[2][i] + m2 * y[5][i]) / (m1 + m2))
xg=np.array(x)
yg=np.array(y1)
zg=np.array(z)

fig=plt.figure(figsize=(10,10))

#Here Program Make Three Type Of Graph ,So Remove Comment Block One At Time To Get Plot

ax=fig.add_subplot(1,1,1,projection="3d")
ax.plot(y[0],y[1],y[2],color="darkblue")
ax.plot(y[3],y[4],y[5],color="tab:red")
ax.plot(xg,yg,zg,color="green")
ax.set_title("Motion Relative to Inertial Frame (Same Mass)\n",fontsize=14)
ax.scatter(y[0,-1],y[1,-1],y[2,-1],color="darkblue",marker="o",s=100,label="m1")
ax.scatter(y[3,-1],y[4,-1],y[5,-1],color="tab:red",marker="o",s=100,label="m2")
ax.scatter(xg[-1],yg[-1],zg[-1],color="green",marker="o",s=100,label="Centre Of Mass")


'''
ax=fig.add_subplot(1,1,1,projection="3d")
ax.plot(y[0]-xg,y[1]-yg,y[2]-zg,color="darkblue",label="m1")
ax.plot(y[3]-xg,y[4]-yg,y[5]-zg,color="tab:red",label="m2")
ax.set_title("Motion of m1 and m2 Relative to Centre of Mass\n",fontsize=14)
'''
'''
ax=fig.add_subplot(1,1,1,projection="3d")
ax.plot(y[3]-y[0],y[4]-y[1],y[5]-y[2],color="darkblue",label="m2")
ax.plot(xg-y[0],yg-y[1],zg-y[2],color="tab:red",label="Centre of Mass")
ax.set_title("Motion of m2 and Centre of mass relative to m1\n",fontsize=14)
'''

ax.set_xlabel("x-coordinate",fontsize=14)
ax.set_ylabel("y-coordinate",fontsize=14)
ax.set_zlabel("z-coordinate",fontsize=14)
ax.legend()
plt.show()
