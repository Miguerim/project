import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation 
from scipy.integrate import solve_ivp
from scipy.integrate import simps
import random as rnd

from typing import Dict, List
from collections import defaultdict
from matplotlib.figure import Figure
from matplotlib.artist import Artist
from matplotlib.animation import ArtistAnimation

from matplotlib import cbook


N=100
q=10
ts=(0, 0.5)

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

R=2000

M=np.zeros((N,N))   
     
for i in range(N):
    for j in range(N):
        if i==0 and j==N-1: M[i][j]=1
        elif i==N-1 and j==0: M[i][j]=1
        elif j>=i+1 and j<=i+q/2: M[i][j]=1
        elif j<=i-1 and j>=i-q/2: M[i][j]=1
        else: M[i][j]=0
        
M=np.array(M)

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



u0=np.ones(2*N,)
u0[:N]*=v_rest
u0[N:2*N]*=0

sol=solve_ivp(system,ts,u0,method='RK23')


# animation code below
    

spike_times = t_spike_arr
spike_times = np.array(spike_times)
for i in range(N):
    for j in range(len(spike_times[i])):
        spike_times[i][j] = round(1000*spike_times[i][j])



neuron_spikes = np.zeros((N, int(1000 * ts[1])))
        

for i in range(N):
    for j in range(int(1000 * ts[1])):
        if j in spike_times[i]:
            neuron_spikes[i][j]=3
            neuron_spikes[i][j-1]=2
            neuron_spikes[i][j-2]=1
            if j < int(1000* ts[1] - 1):
                neuron_spikes[i][j+1]=2
                if j < int(1000*ts[1] - 2):
                    neuron_spikes[i][j+2]=1
                

'''
neuron_spikes = neuron_spikes.tolist()
for i in range(N):
    neuron_spikes[i] = neuron_spikes[i][45:]    
neuron_spikes = np.asarray(neuron_spikes)
'''
   
def centre(n):
    x0=R*np.cos(2*np.pi*n/N+np.pi/2)
    y0=R*np.sin(2*np.pi*n/N+np.pi/2)
    return [x0,y0]


neurons = []
for i in range(N):
    neurons.append([])
    neurons[i].append(centre(i)[0])
    neurons[i].append(centre(i)[1])
    

url = "c://Users/yourName/Desktop" # where you want to download the animations to


fig = plt.figure(figsize=(10,10))
for i in range(N):
    for j in range(N):
        if M[i][j]>0:
            plt.plot([centre(i)[0],centre(j)[0]],[centre(i)[1],centre(j)[1]],c='blue')
        if M[i][j]<0:
            plt.plot([centre(i)[0],centre(j)[0]],[centre(i)[1],centre(j)[1]],c='red')
for n in range(N):
    plt.scatter(centre(n)[0],centre(n)[1],c=[(0, 0, 0)],linewidths=10, zorder=10)
plt.axis('off')


plt.savefig(url+'animation_background.png', bbox_inches='tight', pad_inches=0)

image_file = cbook.get_sample_data(url+'animation_background.png')    
image = plt.imread(image_file)




fig = plt.figure(figsize=(10,10))


class Camera:
    def __init__(self, figure: Figure) -> None:
        self._figure = figure
        self._offsets: Dict[str, Dict[int, int]] = {
            k: defaultdict(int) for k in [
                'collections', 'patches', 'lines', 'texts', 'artists', 'images'
            ]
        }
        self._photos: List[List[Artist]] = []

    def snap(self) -> List[Artist]:
        frame_artists: List[Artist] = []
        for i, axis in enumerate(self._figure.axes):
            if axis.legend_ is not None:
                axis.add_artist(axis.legend_)
            for name in self._offsets:
                new_artists = getattr(axis, name)[self._offsets[name][i]:]
                frame_artists += new_artists
                self._offsets[name][i] += len(new_artists)
        self._photos.append(frame_artists)
        return frame_artists

    def animate(self, *args, **kwargs) -> ArtistAnimation:
        return ArtistAnimation(self._figure, self._photos, *args, **kwargs)

camera = Camera(fig)
for t in range(int(1000 * ts[1])):
    plt.imshow(image,extent=[-2224, 2224, -2224, 2224], aspect='auto')
    plt.text(-2000, 2000, "t="+str(t)+'ms')
    for i in range(N):
        if neuron_spikes[i][t] == 3:
            plt.scatter(neurons[i][0], neurons[i][1], c=[(1, 1, 0)], linewidths=5)
        elif neuron_spikes[i][t] == 2:
            plt.scatter(neurons[i][0], neurons[i][1], c=[(0.6, 0.6, 0)], linewidths=5)
        elif neuron_spikes[i][t] == 1:
            plt.scatter(neurons[i][0], neurons[i][1], c=[(0.3, 0.3, 0)], linewidths=5)
        elif neuron_spikes[i][t] == 0:
            plt.scatter(neurons[i][0], neurons[i][1], c=[(0, 0, 0)], linewidths=5)
    camera.snap()
    print(str(t+1)+'/'+str(int(1000 * ts[1])))

anim = camera.animate()
anim.save(url+'network_animation.gif', writer=animation.PillowWriter(fps=30))
