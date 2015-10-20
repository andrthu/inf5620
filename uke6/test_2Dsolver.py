from wave2D_solver import *
from scitools.std import zeros,cos,amax,pi,sqrt,log,exp,sin

def test_constant():

    k=2.718281828459045
    
    #define functions that give costant solution
    I=lambda x,y: k
    q=lambda x,y: 34 + x
    
    #define some variables
    Lx = 10
    Ly = 4
    T = 3
    C = 1.
    dt= 0.1
    
    #get constant solution
    Nx = int(round(Lx/dt))
    Ny = int(round(Ly/dt))
    ue =zeros((Nx+1,Ny+1)) +k

    u,x,y,t,error=solver(I,None,None,q,0,Lx,Ly,dt,T,C,1,mode='vector')
    print 'Test constant solution. Exercise 3.1'
    print '||abs(u-ue)||: ', amax(abs(u-ue))
    print
    assert amax(abs(u-ue))<10**(-16)

def test_plug():

    #define function that give
    q=lambda x,y: 1
    
    #define some variables
    Lx = 10
    Ly = 3
    T = 1
    C = 1
    dt= 0.1
    
    #get constant solution
    Nx = int(round(Lx/dt))
    Ny = int(round(Ly/dt))
    ue =zeros((Nx+1,Ny+1)) 

    #define initial function
    
    def I(x,y):
        if abs(x-Lx/2.)<0.5:
            return 2.1
        else:
            return 0

    
    #solve numerically
    u,x,y,t,err = solver(I,None,None,q,0,Lx,Ly,dt,T,C,1,mode='scalar')

    #find exact solution at time T
    x1 = x+T
    x2 = x-T
    for i in range(Nx+1):
        for j in range(Ny+1):
            ue[i,j]=0.5*(I(x1[i],0) + I(x2[i],0))
    
    #calculate error at time t=T
    e1 = amax(abs(u-ue))
    

    #Do the exact same for plug wave in y direction:
    def I(x,y):
        if abs(y-Ly/2.)<0.5:
            return 2.1
        else:
            return 0
    #Numerical solution
    u,x,y,t,err = solver(I,None,None,q,0,Lx,Ly,dt,T,C,1,mode='scalar')
    
    #find exact solution at time T
    y1 = x+T
    y2 = x-T    
    for i in range(Nx+1):
        for j in range(Ny+1):
            ue[i,j]=0.5*(I(0,y1[j]) + I(0,y2[j]))


    #calculate error at time t=T
    e2 = amax(abs(u-ue))
    
    print 'Plug error in x ||u(x,y,T)-ue(x,y,T)||: ',e1
    print 'Plug error in y ||u(x,y,T)-ue(x,y,T)||: ',e2
    print
    
    assert e1<10**(-16) and e2<10**(-16)



def test_undampedWaves():
    
    #define constants given in exercise
    A = 1
    mx=7.
    my=2.
    
    #define function that give
    q=lambda x,y: 1
    
    #define some variables
    Lx = 3
    Ly = 1.3
    T = 1
    C = 0.5
    dt= 0.1
    
    #define omega so equation holds
    w=pi*sqrt((mx/Lx)**2 +(my/Ly)**2)
    
    #help varabeles
    kx = pi*mx/Lx
    ky = pi*my/Ly
    
    #Exact solution
    ue = lambda x,y,t: A*cos(x*kx)*cos(y*ky)*cos(t*w)
    
    #initial condition so we get result we want.
    I = lambda x,y: A*cos(x*kx)*cos(y*ky)
    
   
    #factor dt decreeses per step
    step=0.5
    #number of steps I want to do
    val=5
    #array to store errors
    E=zeros(val)
    
    
    
    for i in range(val):
        v='vector'
        #solve eqation
        u,x,y,t,e=solver(I,None,None,q,0,Lx,Ly,dt*step**(i),T,C,1,mode=v,ue=ue)
        
        
        E[i]=e
        
    #find convergence rate between diffrent dt values
    r =zeros(val-1)
    r = log(E[1:]/E[:-1])/log(step)

    
    print "Test convergence for undamped waves:"
    for i in range(val):
        if i==0:
            print "dt: ",dt," Error: ",E[i]
        else:
            print "dt: ",dt*step**(i)," Error: ",E[i], "rate of con.: ", r[i-1]
        
    
    #requiere close to 2 in convergence rate for last r.
    assert abs(r[-1]-2)<0.01 

def test_manufactured():

    A=1
    B=0
    mx=1
    my=1
    
    
    b=1
    c=1.1

    #define some variables
    Lx = 2.
    Ly = 2.
    T = 1
    C = 0.3
    dt= 0.1

    #help varabeles
    kx = pi*mx/Lx
    ky = pi*my/Ly
    w=1

    #Exact solution
    ue = lambda x,y,t: A*cos(x*kx)*cos(y*ky)*cos(t*w)*exp(-c*t)
    I = lambda x,y: A*cos(x*kx)*cos(y*ky)
    V = lambda x,y: -c*A*cos(x*kx)*cos(y*ky)
    q = lambda x,y: x**2+y**2

    f = lambda x,y,t:A*(-b*(c*cos(t*w) + w*sin(t*w))*cos(kx*x)*cos(ky*y) + c**2*cos(kx*x)*cos(ky*y)*cos(t*w) + 2*c*w*sin(t*w)*cos(kx*x)*cos(ky*y) + kx**2*(x**2 + y**2)*cos(kx*x)*cos(ky*y)*cos(t*w) + 2*kx*x*sin(kx*x)*cos(ky*y)*cos(t*w) + ky**2*(x**2 + y**2)*cos(kx*x)*cos(ky*y)*cos(t*w) + 2*ky*y*sin(ky*y)*cos(kx*x)*cos(t*w) - w**2*cos(kx*x)*cos(ky*y)*cos(t*w))*exp(-c*t)
    

    
    
    

    #factor dt decreeses per step
    step=0.5
    #number of steps I want to do
    val=5
    #array to store errors
    E=zeros(val)
    
    
    
    for i in range(val):
        v='vector'
        #solve eqation
        u,x,y,t,e=solver(I,V,f,q,b,Lx,Ly,dt*step**(i),T,C,8,mode=v,ue=ue)
        
        
        E[i]=e
        
    #find convergence rate between diffrent dt values
    r =zeros(val-1)
    r = log(E[1:]/E[:-1])/log(step)

    print
    print "Convergence rates for manufactured solution:" 
    for i in range(val):
        if i==0:
            print "dt: ",dt," Error: ",E[i]
        else:
            print "dt: ",dt*step**(i)," Error: ",E[i], "rate of con.: ", r[i-1]
        
    
    #requiere "close" to 2 in convergence rate for last r.
    assert abs(r[-1]-2)<0.5

    
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
    ue =zeros((Nx+1,Ny+1)) +k

    u1,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,1,mode='scalar',bug='bug1')
    u2,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,1,mode='scalar',bug='bug2')
    u3,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,1,mode='scalar',bug='bug3')
    u4,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,1,mode='scalar',bug='bug4')
    u5,x,y,t=solver(I,None,None,q,0,Lx,Ly,dt,T,C,1,mode='scalar',bug='bug5')
    
    print "Test bugs:" 
    print "bug1 gives error: ",amax(abs(u1-ue))
    print "bug2 gives error: ",amax(abs(u2-ue))
    print "bug3 gives error: ",amax(abs(u3-ue))
    print "bug4 gives error: ",amax(abs(u4-ue))
    print "bug5 gives error: ",amax(abs(u5-ue))
    print

if __name__ == "__main__":
    test_constant()
    constructed_bugs()
    test_plug()
    test_undampedWaves()
    test_manufactured()
