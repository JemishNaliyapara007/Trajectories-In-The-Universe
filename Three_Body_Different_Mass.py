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
    R1=np.array([y[0],y[1]])
    R2 = np.array([y[2], y[3]])
    R3=np.array([y[4],y[5]])
    V1 = np.array([y[6], y[7]])
    V2 = np.array([y[8], y[9]])
    V3 = np.array([y[10], y[11]])
    r12=np.linalg.norm(R2-R1)
    r13 = np.linalg.norm(R3 - R1)
    r23 = np.linalg.norm(R2 - R3)
    A1=G*m2*(R2-R1)/r12**3 + G*m3*(R3-R1)/r13**3
    A2=G*m1*(R1-R2)/r12**3 + G*m3*(R3-R2)/r23**3
    A3 = G * m1 * (R1 - R3) / r13 ** 3 + G * m2 * (R2 - R3) / r23 ** 3
    dydt=np.array([V1,V2,V3,A1,A2,A3]).reshape(12,)
    return dydt


G=6.67259e-20
m1=1.e29
m2=2*1.e29
m3=5*1.e29
t0=0
tf=67000
au=1.4932e8
r1=np.array([0,0])
r2=np.array([300000,0])
r3=np.array([2*300000,0])
v1=np.array([0,0])
v2=np.array([250,250])
v3=np.array([0,0])

m=m1+m2+m3
b0=np.array([r1,r2,r3,v1,v2,v3])
y0=b0.reshape(12,)
#Use Only One Block at Time
#Block:1
# For Reduce Time And Get Better Accuracy Use This Block
t1=np.linspace(t0,tf,1000)
sol=scipy.integrate.solve_ivp(f1,[t0,tf],y0,method='RK45',t_eval=t1)
t=sol.t
y=sol.y
'''
#Block:2
#For RK4 Method, Use This Block 
#It is too much computational complex and take loger time. 
n1=3000000
t,y=RK4(f1,t0,tf,y0,n1)
'''
n=len(t)

x=[]
y1=[]
for i in range(n):
    x.append((m1*y[0][i]+m2*y[2][i]+m3*y[4][i])/m)
    y1.append((m1 * y[1][i] + m2 * y[3][i]+m3*y[5][i])/m)
xg=np.array(x)
yg=np.array(y1)

fig=plt.figure(figsize=(10,10))
#Here Program Make Three Type Of Graph ,So Remove Comment Block One At Time To Get Plot

ax=fig.add_subplot(1,1,1)
ax.plot(y[0],y[1],color="darkblue")
ax.plot(y[2],y[3],color="tab:red")
ax.plot(y[4],y[5],color="green")
ax.plot(xg,yg,color="k")
ax.set_title("Three Body Relative to Inertial Frame\n",fontsize=14)
ax.scatter(y[0,-1],y[1,-1],color="darkblue",marker="o",s=100,label="m1")
ax.scatter(y[2,-1],y[3,-1],color="tab:red",marker="o",s=100,label="m2 = 2*m1 ")
ax.scatter(y[4,-1],y[5,-1],color="green",marker="o",s=100,label="m3 = 5*m1")


'''
ax=fig.add_subplot(1,1,1)
ax.plot(y[0]-xg,y[1]-yg,color="darkblue",label="m1")
ax.plot(y[2]-xg,y[3]-yg,color="tab:red",label="m2 = 2*m1")
ax.plot(y[4]-xg,y[5]-yg,color="green",label="m3 = 5*m1")
ax.set_title("Motion of m1,m2 and m3 Relative to Centre of Mass\n",fontsize=14)
ax.scatter(y[0,0]-xg[0],y[1,0]-yg[0],color="darkblue",marker="o",s=100,label="m1")
ax.scatter(y[2,0]-xg[0],y[3,0]-yg[0],color="tab:red",marker="o",s=100,label="m2")
ax.scatter(y[4,0]-xg[0],y[5,0]-yg[0],color="green",marker="o",s=100,label="m3")
'''

'''
ax=fig.add_subplot(1,1,1)
ax.plot(y[2]-y[0],y[3]-y[1],color="darkblue",label="m2 = 2*m1")
ax.plot(y[4]-y[0],y[5]-y[1],color="green",label="m3 = 5*m1")
ax.plot(xg-y[0],yg-y[1],color="k",label="Centre of Mass")
ax.set_title("Motion of m2 ,m3 and Centre of mass relative to m1\n",fontsize=14)
'''
ax.set_xlabel("x-coordinate",fontsize=14)
ax.set_ylabel("y-coordinate",fontsize=14)
ax.legend()
plt.show()
