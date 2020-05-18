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

M=[]

for i in range(N):
    M.append([])
    for j in range(N):
        M[i].append(rnd.randint(0,1))
        
M=np.array(M)

def S(t,t_spike):
    S=[]
    for i in range(N):
        S.append([])
        for j in range(N):
            s=20*(np.exp(-(t-t_spike[i][j])/tau_d)-np.exp(-(t-t_spike[i][j])/tau_r))
            S[i].append(s)
    return S
    
def G_out(S):
    G_out=[]
    for i in range(N):
        g_out=0
        for j in range(N):
            g_out=g_out+S[i][j]
        G_out.append(g_out)
    return G_out
        
def B(t,u):
    B=[]
    for i in range(N):
        B.append([])
        for j in range(N):
            B[i].append(M[i][j]*(v_syn-u[i]))
    return B
            
def I_syn(G_out,B):
    I_syn=np.dot(G_out,B)
    return I_syn
    
def I(t):
    if t>t0 and t<t1:
        I=5
    else:
        I=0
    return I



        
def system(t,u):
    du=np.zeros((N,))
    t_spike=t[u==v_rest]
    du=1/C*(g_L*(v_rest-u)+I_syn(G_out(S(t,t_spike),B(t,u))+I(t))
    u[u>-30]=v_rest
    return du

u0=v_rest*np.ones(N,)

ts=(0,1000)

sol=solve_ivp(system,ts,u0)

        
