from numpy import *
import matplotlib.pyplot as plt

L=lambda x,y:sqrt(x**2+(y-1)**2)
 
def simulate(beta=0.9,Theta=30,epsilon=0,num_periods=6,
             time_steps_per_period=60, plot=True):
    
    n=num_periods*time_steps_per_period
    dt=2*pi/time_steps_per_period

    x=zeros(n+1)
    y=zeros(n+1)
    tet=zeros(n+1)
    t=linspace(0,dt*n,n+1)
    
    L1=L(x[0],y[0])
    x[0]=(1+epsilon)*sin(Theta)
    y[0]=1-(1+epsilon)*cos(Theta)
    x[1]=x[0]*(1-(dt**2*0.5*beta/(1-beta))*(1-beta/L1))

    y[1]=y[0]-(y[0]-1)*beta*(dt**2)*(1-beta/L1)/(1-beta)-0.5*beta*dt**2

    for i in range(1,n):
        L1=L(x[i],y[i])

        x[i+1]=2*x[i]-x[i-1]-x[i]*beta*(dt**2)*(1-beta/L1)/(1-beta)
        y[i+1]=2*y[i]-y[i-1]-(y[i]-1)*beta*(dt**2)*(1-beta/L1)/(1-beta)
        -beta*dt**2
    
    tet=arctan(x/(y-1))

    if plot:
        plt.figure()
        #plt.gca().set_aspect('equal')
        plt.plot(x,y)
        plt.show()
        
    return x,y,tet,t



if __name__=="__main__":
    simulate()
    x,y,tet,t=simulate(beta=0.9,Theta=0,epsilon=0,num_periods=6,
             time_steps_per_period=60, plot=True)
    print L(0,0)
