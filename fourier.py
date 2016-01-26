import sympy as sym
from numpy import *
import matplotlib.pyplot as plt
def phi(x):

    a = lambda y: piecewise(y,[y<1./3],[0,1])
    b = lambda y: piecewise(y,[y>2./3],[0,1])

    return a(x)*b(x)


x= linspace(0,1,1000)
y=phi(x)





def fori(N,x):
    
    c,s = sym.mpmath.fourier(phi,[0,1],N)
    print c
    
    S=0
    for i in range(N):
        S=S + c[i]*cos(i*2*pi*x) + s[i]*sin(i*2*pi*x)
        

    return S

plt.plot(x,fori(10,x))
plt.plot(x,y)
plt.show()
