from numpy import *
from matplotlib.pylab import *

def mesh_function(f, t):
    sol = zeros(len(t))

    for i in range(len(t)):
        sol[i]=f(t[i])
    return sol

def F(x):
    if x<=3:
        return exp(-x)
    else:
        return exp(-3*x)

t=linspace(2,4,200)

plot(t,mesh_function(F,t))
show()
