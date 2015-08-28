from numpy import *
from matplotlib.pylab import *

def diferentiate(u, dt):
    n=len(u)

    d=zeros(n)

    d[0]=(u[1]-u[0])/dt
    d[n-1]=(u[n-1]-u[n-2])/dt

    for i in range(1,n-1):
        d[i]=(u[i+1]-u[i-1])/(2*dt)

    return d

def test_diferentiate():
    t = linspace(0,4,40)
    dt=0.1
    f=lambda x: x*x +x+1
    g=lambda x: 2*x +1

    error = sum(abs(diferentiate(f(t),dt)-g(t)))
    plot(t,diferentiate(f(t),dt))
    plot(t,g(t))
    show()

    print error


def diferentiate2(u, dt):
    n=len(u)

    d=zeros(n)

    d[0]=(u[1]-u[0])/dt
    d[n-1]=(u[n-1]-u[n-2])/dt

    d[1:-1] = (u[2:] - u[0:-2])/(2*dt)

    return d

def test_diferentiate2():
    t = linspace(0,1,10000)
    dt=0.0001
    f=lambda x: x*x +x+1
    g=lambda x: 2*x +1

    error = sum(abs(diferentiate2(f(t),dt)[1:-1]-g(t)[1:-1]))
    plot(t,diferentiate2(f(t),dt))
    plot(t,g(t))
    show()

    print error*dt
    
test_diferentiate2()
