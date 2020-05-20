import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import random as rnd

N=10

C=0.526
g_L=26.3
v_rest=-70.0
v_syn=0.0
theta=-52.0 
t_refrac =2.0
tau_r=0.5
tau_d=2.0 

t0=1
t1=5

M=[]

for i in range(N):
    M.append([])
    for j in range(N):
        M[i].append(rnd.randint(0,1))
        
M=np.array(M)


def t_spike(t,u):
    t_spike=[]
    for i in range(N):
        t_spike.append([])
        if u[i]==v_rest: 
            t_spike[i].append(t)
    return t_spike

def S(t,t_spike):
    S=[]
    for i in range(N):
        S.append([])
        for j in range(len(t_spike[i])):
            s=-1*(np.exp(-(t-t_spike[i][j])/tau_d)-np.exp(-(t-t_spike[i][j])/tau_r))
            S[i].append(s)
    return S
    
def G_out(S):
    G_out=[]
    for i in range(N):
        G_out.append(np.sum(S[i]))
    return G_out
        
def B(t,u):
    B=[]
    for i in range(N):
        B.append(np.multiply(M[i],v_syn-u))
    return B
            
def I_syn(G_out,B):
    I_syn=np.dot(G_out,B)
    return I_syn
    
def I(t):
    I=np.random.uniform(490,510)
    return I


        
def system(t,u):
    du=np.zeros((N,))
    u[u>theta]=v_rest
    du=1/C*(g_L*(v_rest-u)+I_syn(G_out(S(t,t_spike(t,u))),B(t,u))+I(t))
    print(I_syn(G_out(S(t,t_spike(t,u))),B(t,u)))
    return du

u0=v_rest*np.ones(N,)

ts=(0,1)

sol=solve_ivp(system,ts,u0,method='RK23')

plt.plot(sol.t, sol.y[0])

#fig, axs = plt.subplots(N)
#for i in range(N):
#    axs[i].plot(sol.t, sol.y[i])