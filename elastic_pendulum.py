from numpy import *

L=lambda x,y:sqrt(x**2+(y-1)**2)
 
def simulate(beta=0.9,Theta=30,epsilon=0,num_periods=6,
             time_steps_per_period=60, plot=True):
    
    n=60*6*2*pi
    dt=2*pi/60

    x=zeros(n+1)
    y=zeros(n+1)
    t=linspace(0,dt*n,n+1)

    x[0]=(1+epsilon)*sin(Theta)
    y[0]=1-(1+epsilon)*cos(Theta)
    x[1]=x[0]*(1-(0.5*beta/(1-beta))*(1-beta/L(x[0],y[0])))
    y[1]=y[0]*(1-(0.5*beta/(1-beta))*(1-beta/L(x[0],y[0])))+
    (0.5*beta/(1-beta))*(1-beta/L(x[0],y[0]))-0.5*beta

    for i in range(1,n):
        x[
    
