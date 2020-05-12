import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation

N=100

tau_s=1
tau_n=1
R=1

t0=0
t1=10

v_rest=65

M=[]


for i in range(N):
    M.append([])
    for j in range(N):
        M[i].append(1)
        
M=np.array(M)
        
def I(t):
    if t0<t<t1:
        I=5
    else: I=0.0    
    return I
        
def system(S,T):
    dS=np.zeros((2,N))
    
    v=S[0]
    i=S[1]
    
    dS[0]=-1/tau_s*(v-v_rest)+R*i
    dS[1]= -I(T)+np.dot(M,v)
    
    return dS

S0=np.array([0.0,0])

ts=np.linspace(0,1000,10000)

S1=odeint(system,S0,ts)
