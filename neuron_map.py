import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation 
from scipy.integrate import solve_ivp
from scipy.integrate import simps
import random as rnd

N=100
q=10

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

#Y=[]

#for x in np.linspace(-R,R,100000):
#    y1=(R**2-x**2)**0.5
#    y2=-(R**2-x**2)**0.5
#    Y.append([x,y1])
#    Y.append([x,y2])
 
#Y=np.array(Y)
   
def centre(n):
    x0=R*np.cos(2*np.pi*n/N+np.pi/2)
    y0=R*np.sin(2*np.pi*n/N+np.pi/2)
    return [x0,y0]

Centre=[]

for n in range(N):
    x0=R*np.cos(2*np.pi*n/N+np.pi/2)
    y0=R*np.sin(2*np.pi*n/N+np.pi/2)
    Centre.append([x0,y0])

Centre=np.array(Centre)

#def neuron(n):
#    neuron=[]
#    for x in np.linspace(-R-50,R+50,10000):
#        x0=centre(n)[0]
#        y0=centre(n)[1]
#        y1=y0-(25-(x-x0)**2)**0.5
#        y2=y0+(25-(x-x0)**2)**0.5
#        neuron.append([x,y1])
#        neuron.append([x,y2])
#    neuron=np.array(neuron)
#    return neuron


   
plt.figure(figsize=(10,10))

#plt.scatter(Y[:,0],Y[:,1],c='black',linewidths=5)

for i in range(N):
    for j in range(N):
        if M[i][j]>0:
            plt.plot([centre(i)[0],centre(j)[0]],[centre(i)[1],centre(j)[1]],c='blue')
        if M[i][j]<0:
            plt.plot([centre(i)[0],centre(j)[0]],[centre(i)[1],centre(j)[1]],c='red')
            
for n in range(N):
    plt.scatter(centre(n)[0],centre(n)[1],c='green',linewidths=10)
    
plt.show()
    


#line = None

#fig =plt.figure(figsize=(10,10))
#ax=fig.add_subplot(1,1,1)

#def animate(i):
    
#    global line
    
#    if line is not None:
#        line,=None
    
#    line,=ax.scatter(Centre[:,0],Centre[:,1],c='green',linewidths=10)
    
#    return line,

#anime = animation.FuncAnimation(fig,animate)



    




