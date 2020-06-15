import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation 
from scipy.integrate import solve_ivp
from scipy.integrate import simps
import random as rnd

from matplotlib.animation import FuncAnimation

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

t0=0
t1=10
P=0.2
P1=0.2


M=np.zeros((N,N))   
     
for i in range(N):
    for j in range(N):
        if i==0 and j==N-1: M[i][j]=1
        elif i==N-1 and j==0: M[i][j]=1
        elif j>=i+1 and j<=i+q/2: M[i][j]=1
        elif j<=i-1 and j>=i-q/2: M[i][j]=1
        else: M[i][j]=0
        
M=np.array(M)



connections = 0
for i in range(N):
    for j in range(N):
        if abs(M[i][j]==1):
            connections += 1
print(connections)


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
                
                for l in list:
                    if l==j: list.remove(l)
                    
                y=np.random.choice(list)
                while M[i][int(y)]==1:
                    y=np.random.choice(list)
                M[i][int(y)]=1


for i in range(N):
    for j in range(N):
        z=np.random.choice(a=[0,1],p=[1-P1,P1])
        if z==1:
            M[i][j]*=-1
            

connections = 0
for i in range(N):
    for j in range(N):
        if abs(M[i][j]) == 1:
            connections += 1
print(connections)


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
        if t>t0 and t<t1:
            I.append(np.random.uniform(490,510))
        else: I.append(0)
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

ts=(0,1)


sol=solve_ivp(system,ts,u0,method='RK23')

'''

plt.figure()
plt.plot(sol.t, sol.y[4])
plt.plot(sol.t,sol.y[0])
plt.show()

plt.figure()
plt.plot(sol.t,sol.y[104])
plt.plot(sol.t,sol.y[100])
plt.show()

plt.figure()
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

plt.figure()
plt.plot(sol.t[1:],Chi)
plt.show()



'''
A=[]
for t in sol.t:
    n=0
    for i in range(N):
        for j in range(len(t_spike_arr[i])):
            if t==t_spike_arr[i][j]:
                n=n+1
    n=n/N
    A.append(n)
'''
plt.figure()
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
 
plt.figure() 
plt.plot(sol.t,Entropy)
plt.show()



#U-I
du=np.diff(sol.y[0])
i=C*(du-g_L*(v_rest-sol.y[0][1:])-g_L*delta*np.exp((sol.y[0][1:]-theta)/delta)+sol.y[100][1:])


plt.figure()
plt.plot(i,sol.y[0][1:])
plt.show()


#U_n+1 - U_n
plt.figure()
plt.plot(sol.y[0][1:],sol.y[0][:len(sol.t)-1])
plt.show()



#A_n+1 - A_n
plt.figure()
plt.plot(A[1:],A[:len(A)-1])
plt.show()



#Entropy_n+1 - Entropy_n
plt.figure()
plt.plot(Entropy[1:],Entropy[:len(Entropy)-1])
plt.show()








# Animations:

url = "c://Users/eddyt/Desktop/network/" # where you want to download the animations to

fig = plt.figure() 
ax = plt.axes(xlim=(ts[0], ts[1]), ylim=(min(sol.y[0]), max(sol.y[0])))
line, = ax.plot([], [], lw=2) 

def init(): 
	line.set_data([], []) 
	return line, 

xdata, ydata = [], [] 

def animate(i): 
	x = sol.t[i]
	y = sol.y[0][i]	
	xdata.append(x) 
	ydata.append(y) 
	line.set_data(xdata, ydata) 
	return line, 

plt.title('Neuron Spiking')
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(sol.t), interval=20, blit=True, repeat=False) 
anim.save(url+'neuron_spiking.gif', writer=animation.PillowWriter(fps=30))




fig = plt.figure() 
ax = plt.axes(xlim=(ts[0], ts[1]), ylim=(0, 1))
line, = ax.plot([], [], lw=2) 

def init(): 
	line.set_data([], []) 
	return line, 

xdata, ydata = [], [] 

def animate(i): 
	x = sol.t[i]
	y = Chi[i]
	xdata.append(x) 
	ydata.append(y) 
	line.set_data(xdata, ydata) 
	return line, 

plt.title('Synchronization')
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(len(sol.t)-1), interval=20, blit=True, repeat=False) 
anim.save(url+'synchronization.gif', writer=animation.PillowWriter(fps=30))




fig = plt.figure() 
ax = plt.axes(xlim=(ts[0], ts[1]), ylim=(min(Entropy), max(Entropy)))
line, = ax.plot([], [], lw=2) 

def init(): 
	line.set_data([], []) 
	return line, 

xdata, ydata = [], [] 

def animate(i):
	x = sol.t[i]
	y = Entropy[i]
	xdata.append(x) 
	ydata.append(y) 
	line.set_data(xdata, ydata) 
	return line, 

plt.title('Entropy of the Network')
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(sol.t), interval=20, blit=True, repeat=False) 
anim.save(url+'entropy.gif', writer=animation.PillowWriter(fps=30))




ii = i
fig = plt.figure() 
ax = plt.axes(xlim=(min(ii), max(ii)), ylim=(min(sol.y[0][1:]), max(sol.y[0][1:])))
line, = ax.plot([], [], lw=2) 

def init(): 
	line.set_data([], []) 
	return line, 

xdata, ydata = [], [] 

def animate(i):
	x = ii[i]
	y = sol.y[0][i+1]
	xdata.append(x) 
	ydata.append(y) 
	line.set_data(xdata, ydata) 
	return line, 

plt.title('U-I')
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(len(sol.t)-1), interval=20, blit=True, repeat=False) 
anim.save(url+'U-I.gif', writer=animation.PillowWriter(fps=30))




fig = plt.figure() 
ax = plt.axes(xlim=(min(sol.y[0][1:]), max(sol.y[0][1:])), ylim=(min(sol.y[0][:len(sol.t)-1]), max(sol.y[0][:len(sol.t)-1])))
line, = ax.plot([], [], lw=2) 

def init(): 
	line.set_data([], []) 
	return line, 

xdata, ydata = [], [] 

def animate(i):
	x = sol.y[0][i+1]
	y = sol.y[0][i]
	xdata.append(x) 
	ydata.append(y) 
	line.set_data(xdata, ydata) 
	return line, 

plt.title('U_n+1-U_n')
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(len(sol.t)-1), interval=20, blit=True, repeat=False) 
anim.save(url+'U_n+1-U_n.gif', writer=animation.PillowWriter(fps=30))




fig = plt.figure() 
ax = plt.axes(xlim=(min(A[1:]), max(A[1:])), ylim=(min(A[:len(A)-1]), max(A[:len(A)-1])))
line, = ax.plot([], [], lw=2) 

def init(): 
	line.set_data([], []) 
	return line, 

xdata, ydata = [], [] 

def animate(i):
	x = A[i+1]
	y = A[i]
	xdata.append(x) 
	ydata.append(y) 
	line.set_data(xdata, ydata) 
	return line, 

plt.title('A_n+1 - A_n')
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(len(sol.t)-1), interval=20, blit=True, repeat=False) 
anim.save(url+'A_n+1-A_n.gif', writer=animation.PillowWriter(fps=30))




fig = plt.figure() 
ax = plt.axes(xlim=(min(Entropy[1:]), max(Entropy[1:])), ylim=(min(Entropy[:len(Entropy)-1]), max(Entropy[:len(Entropy)-1])))
line, = ax.plot([], [], lw=2) 

def init(): 
    	line.set_data([], []) 
	return line, 

xdata, ydata = [], [] 

def animate(i):
	x = Entropy[i+1]
	y = Entropy[i]
	xdata.append(x) 
	ydata.append(y) 
	line.set_data(xdata, ydata) 
	return line, 

plt.title('Entropy_n+1 - Entropy_n')
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(len(Entropy)-1), interval=20, blit=True, repeat=False) 
anim.save(url+'Entropy_n+1-Entropy_n.gif', writer=animation.PillowWriter(fps=30))
'''

url = "c://Users/yourName/Desktop/" # where you want to download the animations to

fig = plt.figure() 
ax = plt.axes(xlim=(min(A[1:]), max(A[1:])), ylim=(min(A[:len(A)-1]), max(A[:len(A)-1])))
line, = ax.plot([], [], lw=2) 

def init(): 
	line.set_data([], []) 
	return line, 

xdata, ydata = [], [] 

def animate(i):
    x = A[i+1]
    y = A[i]
    xdata.append(x) 
    ydata.append(y) 
    if i > 25:
        del xdata[0]
        del ydata[0]
    line.set_data(xdata, ydata) 
    return line, 

plt.title('A_n+1 - A_n')
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(len(sol.t)-1), interval=20, blit=True, repeat=False) 
anim.save(url+'A_n+1-A_n7.gif', writer=animation.PillowWriter(fps=30))