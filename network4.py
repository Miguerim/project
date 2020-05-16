import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import animation

N=100

tau_s=1
tau_n=1
R=1

t0=0
t1=10

v_rest=-65

M=[]


for i in range(N):
    M.append([])
    for j in range(N):
        M[i].append(1)
        
M=np.array(M)
        
def I(t):
    I=5   
    return I
        
def system(t,u):
    du=np.zeros((N,))
    
    if u<-30:
        du=-1/tau_s*(u-v_rest)+R*(I(t)-np.dot(u,M))
    else:
        u=v_rest
    return du

u0=v_rest*np.ones(N,)

ts=(0,1000)

sol=solve_ivp(system,ts,u0)
