from numpy import *

def solver(I,w,V,dt,T,f):
    
    n=int(round(T/dt))
    u = zeros(n+1)
    t=linspace(0,n*dt,n+1)

    u[0]=I
    u[1]=(I*(2-(dt*w)**2)+2*V*dt+f(0)*dt**2)/2
    for i in range(1,n):
        u[i+1]=u[i](2-(dt*w)**2)-u[i-1]+f(i*dt)*dt**2

    return u,t
