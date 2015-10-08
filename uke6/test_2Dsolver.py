from wave2D_solver import *
import numpy as np

def test_constant():

    k=2.718281828459045
    
    #define functions that give costant solution
    I=lambda x,y: k
    q=lambda x,y: 76.3
    
    #define some variables
    Lx = 10
    Ly = 5
    T = 3
    C = 1.
    dt= 0.1
    
    #get constant solution
    Nx = int(round(Lx/dt))
    Ny = int(round(Ly/dt))
    ue =np.zeros((Nx,Ny)) +k

    u,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,mode='scalar')

    assert np.amax(abs(u-ue))<10**(-16)

def test_plug():

    #define function that give
    q=lambda x,y: 1
    
    #define some variables
    Lx = 4
    Ly = 1
    T = 1
    C = 1
    dt= 0.1
    
    #get constant solution
    Nx = int(round(Lx/dt))
    Ny = int(round(Ly/dt))
    ue =np.zeros((Nx,Ny)) 

    #define initial function
    
    def I(x,y):
        if abs(x-Lx/2.)<0.5:
            return 2.1
        else:
            return 0

    
    
    u,x,y,t = solver(I,None,None,q,0,Lx,Ly,dt,T,C,mode='scalar')
    x1 = x+T
    x2 = x-T

    for i in range(Nx):
        for j in range(Ny):
            ue[i,j]=0.5*(I(x1[i],0) + I(x2[i],0))

    assert np.amax(abs(u-ue))<10**(-16)
    

def constructed_bugs():
    k=2.718281828459045

    I=lambda x,y: k
    q=lambda x,y: x+y

    Lx = 10
    Ly = 5
    T = 3
    C = 1.
    dt= 0.5
    
    Nx = int(round(Lx/dt))
    Ny = int(round(Ly/dt))
    ue =np.zeros((Nx,Ny)) +k

    u1,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,mode='scalar',bug='bug1')
    u2,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,mode='scalar',bug='bug2')
    u3,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,mode='scalar',bug='bug3')
    u4,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,mode='scalar',bug='bug4')
    u5,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,mode='scalar',bug='bug5')

    print "bug1 gives error: ",np.amax(abs(u1-ue))
    print "bug2 gives error: ",np.amax(abs(u2-ue))
    print "bug3 gives error: ",np.amax(abs(u3-ue))
    print "bug4 gives error: ",np.amax(abs(u4-ue))
    print "bug5 gives error: ",np.amax(abs(u5-ue))

if __name__ == "__main__":
    #test_constant()
    #constructed_bugs()
    test_plug()
