from numpy import *
from matplotlib.pyplot import *

def solver(I, a, T, dt, theta):
    """Solve u'=-a*u, u(0)=I, for t in (0,T] with steps of dt."""
    dt = float(dt)            # avoid integer division
    Nt = int(round(T/dt))     # no of time intervals
    T = Nt*dt                 # adjust T to fit time step dt
    u = zeros(Nt+1)           # array of u[n] values
    t = linspace(0, T, Nt+1)  # time mesh

    u[0] = I                  # assign initial condition
    for n in range(0, Nt):    # n=0,1,...,Nt-1
        u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[n]
    return u, t

def decay_plot_error(I, a, T, dt):
    f=lambda x: I*exp(-a*x)
    K=[0,0.5,1]
    for i in K: 
        u,t = solver(I, a, T, dt,i)
        
        plot(t,f(t)-u)

    show()

decay_plot_error(2,1,1,0.1)
