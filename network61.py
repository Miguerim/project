import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import simps
import random as rnd

N=100
q=10

C=0.526
g_L=26.3
v_rest=-70.0
v_syn=0.0
theta=-52.0 
t_refrac =2.0
tau_r=0.5
tau_d=2.0 

tau_w=93e-3
alpha=0.82
b=15
delta=2

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
        elif j>=i+1 and j<=i+q/2: M[i][j]=1
        elif j<=i-1 and j>=i-q/2: M[i][j]=1
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

def Freq(t):
    t_spike_arr=np.array(t_spike_arr)
    freq=[]
    for i in range(N):
        number_spikes=len(t_spike_arr[i]<t)
        freq.append(number_spikes)
    
    Freq=[]
    for i in range(N):
        Freq.append([])
        for j in range(N):
            freq_multiplied=freq[i]*freq[j]
            Freq[i].append(freq_multiplied)
            
    t_spike_arr=t_spike_arr.tolist()
    return Freq
            

def S(t,t_spike):
    S=[]
    for i in range(N):
        S.append([])
        for j in range(len(t_spike[i])):
            if t-t_spike[i][j]<tau_d:
                s=1*(np.exp(-1e-3*(t-t_spike_arr[i][j])/tau_d)-np.exp(-1e-3*(t-t_spike_arr[i][j])/tau_r))
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
    df[:N]=1/C*(g_L*(v_rest-u)+I_syn(G_out(S(t,t_spike(t,u))),B(t,u))+I(t)+g_L*delta*np.exp((u-theta)/delta)-w)
    w[u>theta]=w[u>theta]+b
    return df

#u0=v_rest*np.ones(N,)

u0=np.ones(2*N,)
u0[:N]*=v_rest
u0[N:2*N]*=0

ts=(0,3)

sol=solve_ivp(system,ts,u0,method='RK23')

plt.plot(sol.t, sol.y[4])
plt.plot(sol.t,sol.y[0])
plt.show()

plt.plot(sol.t,sol.y[104])
plt.plot(sol.t,sol.y[100])
plt.show()


#fig, axs = plt.subplots(N)
#for i in range(N):
#    axs[i].plot(sol.t, sol.y[i])


plt.eventplot(t_spike_arr)
plt.show()

Chi=[]

for i in range(1,len(sol.t)):
    v_hat = (1/N)*np.sum(sol.y[:,:i],axis=0)
    v_hat_2_avg=(1/sol.t[i])*simps(np.multiply(v_hat,v_hat),sol.t[:i])
    v_hat_avg_2=((1/sol.t[i])*simps(v_hat,sol.t[:i]))**2
    sigma_hat_2=v_hat_2_avg-v_hat_avg_2
    
    Sigma_2=[]
    
    for j in range(N):
        v_2_avg=(1/sol.t[i])*simps(np.multiply(sol.y[j][:i],sol.y[j][:i]),sol.t[:i])
        v_avg_2=((1/sol.t[i])*simps(sol.y[j][:i],sol.t[:i]))**2
        sigma_2=v_2_avg-v_avg_2
        Sigma_2.append(sigma_2)
        
    mean_sigma_2=(1/N)*np.sum(Sigma_2)
    
    chi=(sigma_hat_2/mean_sigma_2)**(1/2)
    Chi.append(chi)

plt.plot(sol.t[1:],Chi)
plt.show()

#print(created)
#print(destroyed)
A=[]
for t in sol.t:
    n=0
    for i in range(N):
        for j in range(len(t_spike_arr[i])):
            if t==t_spike_arr[i][j]:
                n=n+1
    n=n/N
    A.append(n)
    

plt.plot(sol.t,A)
plt.show()

gradA=np.gradient(A)

lambdaA=(1/N)*np.sum(np.log(abs(gradA)))
Entropy = []
Weights=[]
for t in range(len(sol.t)):
    
    weights=np.zeros(100)

    for i in range(len(A[:t])):
        for j in range(len(weights)):
            if 100*A[i]<=j+1 and 100*A[i]>j:
                weights[j]=weights[j]+1
            
    if np.sum(weights)!=0: weights*=1/np.sum(weights)

    weights=weights.tolist()
    while 0 in weights: weights.remove(0)
    weights=np.array(weights)
    Weights.append(weights)
    entropy=-np.sum(np.multiply(weights,np.log(weights)))
    Entropy.append(entropy)
 
plt.plot(sol.t,Entropy)
plt.show()

plt.hist(Weights[len(sol.t)-1])
plt.show()