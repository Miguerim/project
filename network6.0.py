import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import simps
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

tau_w=93
alpha=0.82
b=1

t0=1
t1=5
P=0.2


#M=[]

#for i in range(N):
#    M.append([])
#    for j in range(N):
#        M[i].append(rnd.randint(0,1))

M=np.zeros((N,N))   
     
for i in range(N):
    for j in range(N):
        if i==0 and j==N-1: M[i][j]=1
        elif i==N-1 and j==0: M[i][j]=1
        elif j==i+1: M[i][j]=1
        elif j==i-1: M[i][j]=1
        else: M[i][j]=0
        
M=np.array(M)

created=0
destroyed=0   
P1=0.2


for i in range(N):
    for j in range(N):
        if M[i][j]==1:
            x=np.random.choice(a=[0,1],p=[1-P,P])
            if x==1:
                list=np.linspace(0,N-1,N)
                list=list.tolist()
                for k in list: 
                    if k==i: list.remove(k)
                M[i][j]=0
                #destroyed+=1
                for l in list:
                    if l==j: list.remove(l)
                y=np.random.choice(list)
                M[i][int(y)]=1
                #created+=1

for i in range(N):
    for j in range(N):
        z=np.random.choice(a=[0,1],p=[1-P1,P1])
        if z==1:
            M[i][j]*=-1
            

t_spike_arr=[]            
for i in range(N):
    t_spike_arr.append([])

def t_spike(t,u):
    for i in range(N):
        if t!=0:
            if u[i]==v_rest: 
                t_spike_arr[i].append(t)
    return t_spike_arr

def S(t,t_spike):
    S=[]
    for i in range(N):
        S.append([])
        for j in range(len(t_spike[i])):
            if t-t_spike[i][j]<tau_d:
                s=1*(np.exp(-(t-t_spike_arr[i][j])/tau_d)-np.exp(-(t-t_spike_arr[i][j])/tau_r))
            else: s=0
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
    I=[]
    for i in range(N):
        I.append(np.random.uniform(490,510))
    return I

def D_w(t):
    Delta_w=[]
    for i in range(N):
        Delta_w.append([])
        for j in range(len(t_spike_arr[i])):
            if t==t_spike_arr[i][j]:
                delta_w=b
            else: delta_w=0
            Delta_w[i].append(delta_w)
    print(np.shape(Delta_w))
    print(Delta_w)
    return Delta_w
        
def system(t,f):
    print(t)
    df=np.zeros((2*N,))
    u=f[:N]
    w=f[N:]
    u[u>theta]=v_rest
    df[N:]=(1/tau_w)*(alpha*(u-v_rest)-w)
    df[:N]=1/C*(g_L*(v_rest-u)+I_syn(G_out(S(t,t_spike(t,u))),B(t,u))+I(t)-w)
    w[u>theta]=w[u>theta]+b
    return df

#u0=v_rest*np.ones(N,)

u0=np.ones(2*N,)
u0[:N]*=v_rest
u0[N:2*N]*=0

ts=(0,0.1)

sol=solve_ivp(system,ts,u0,method='RK23')

plt.plot(sol.t, sol.y[4])
plt.plot(sol.t,sol.y[0])
plt.show()

plt.plot(sol.t,sol.y[14])
plt.plot(sol.t,sol.y[10])
plt.show()


#fig, axs = plt.subplots(N)
#for i in range(N):
#    axs[i].plot(sol.t, sol.y[i])


plt.eventplot(t_spike_arr)
plt.show()


v_hat = (1/N)*np.sum(sol.y,axis=0)
v_hat_2_avg=(1/ts[1])*simps(np.multiply(v_hat,v_hat),sol.t)
v_hat_avg_2=((1/ts[1])*simps(v_hat,sol.t))**2
sigma_hat_2=v_hat_2_avg-v_hat_avg_2

Sigma_2=[]

for i in range(N):
    v_2_avg=(1/ts[1])*simps(np.multiply(sol.y[i],sol.y[i]),sol.t)
    v_avg_2=((1/ts[1])*simps(sol.y[i],sol.t))**2
    sigma_2=v_2_avg-v_avg_2
    Sigma_2.append(sigma_2)
    
mean_sigma_2=(1/N)*np.sum(Sigma_2)

chi=(sigma_hat_2/mean_sigma_2)**(1/2)

print(chi)

#print(created)
#print(destroyed)
