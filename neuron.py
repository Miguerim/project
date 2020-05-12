import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation

V_Na=115
V_K=-12
V_L=10.613
g_Na=120
g_K=36
g_L=0.3
C=1
t0=0
t1=1
#print("Dla Ogarniaczy Fizyki przygotowal: Michal Piotrak")
#omega=float(input("Wprowadz czestotliwosc pradu w jednostkach PI (3,14) radianow: "))
I0=5 #float(input("Wprowadz natezenie pradu w miliamperach: "))

def alpha_n(v):
    return 0.01*(-v+10)/(np.exp((-v+10)/10) -1)

def beta_n(v):
    return 0.125*np.exp(-v/80)

def alpha_m(v):
    return 0.1*(-v+25)/(np.exp((-v+25)/10) -1)

def beta_m(v):
    return 4*np.exp(-v/18)
    
def alpha_h(v):
    return 0.07*np.exp(-v/20)

def beta_h(v):
    return 1/(np.exp((-v+30)/10) +1)

def n_inf(v=0.0):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))

def m_inf(v=0.0):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))

def h_inf(v=0.0):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def I(t):
    if t0<t<t1:
        I=I0
    else: I=0.0    
    return I

    
def system(Y,T):
    dY=np.zeros((4,))
    
    v=Y[0]
    n=Y[1]
    m=Y[2]
    h=Y[3]
    

    dY[0]=(1/C)*(I(T)-(g_Na*(m**3)*h*(v-V_Na)+g_K*(n**4)*(v-V_K)+g_L*(v-V_L)))
    dY[1]=alpha_n(v)*(1-n)-beta_n(v)*n
    dY[2]=alpha_m(v)*(1-m)-beta_m(v)*m
    dY[3]=alpha_h(v)*(1-h)-beta_h(v)*h
    
    return dY

Y0=np.array([0.0,n_inf(),m_inf(),h_inf()])

ts=np.linspace(0,100,1000)

Y1=odeint(system,Y0,ts)

Y1[:,0]=Y1[:,0]-65

plt.plot(ts,Y1[:,0])
plt.xlabel('czas [ms]')
plt.ylabel('napiecie [mV]')
plt.show()

#tspike=[]

#for i in range (len(Y1[:,0])):
#    if Y1[i][0]==max(Y1[0:200,0]):
#        tspike.append(ts[i])
#    if Y1[i][0]==max(Y1[200:600,0]):
#        tspike.append(ts[i])
 #   if Y1[i][0]==max(Y1[600:1000,0]):
 #       tspike.append(ts[i])
        
#print(tspike)