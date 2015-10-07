import numpy as np

def solver(I,V,f,q,b,Lx,Ly,dt,T,C):
    
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1)

    dx = dt/C 
    Nx = int(round(Lx/dx))
    x = np.linspace(0, Nx*dx, Nx+1)
    xv = x[:newaxis]
    
    dy = dt/C 
    Ny = int(round(Ly/dy))
    y = np.linspace(0, Ny*dx, Ny+1)
    yv = y[newaxis:]

    C1 = (dt/dx)**2
    C2 = (dt/dy)**2
    
    F = np.zeros(shape=(Nx,Ny))
    Q = np.zeros(shape=(Nx,Ny))
    V_v = np.zeros(shape=(Nx,Ny))
    I_v = np.zeros(shape=(Nx,Ny))

    u = np.zeros(shape=(Nx,Ny))
    u_1 = np.zeros(shape=(Nx,Ny))
    u_2 = np.zeros(shape=(Nx,Ny))

    #handeling zero sorce term and initial condition
    if f is None or f == 0 :
        f = lambda x,y,t: 0
    if V is None or V == 0:
        V = lambda x,y: 0

    #making variabels to make vectorization easier.
    
    F[:,:] = f(xv[:],yv[:],0)
    Q[:,:] = q(xv[:],yv[:])
    V_v[:,:] = V(xv[:],yv[:])
        
    
    #The first initial condition
    for i in range(1,Nx):
        for j in range(1,Ny):
            u_1[i,j]=I(dx*i,dy*j)

    #The second initial condiyion at inner points
    for i in range(1,Nx-1):
        for j in range(1,Ny-1):
            u[i,j] = step1(u_1,x,y,t,f,V,C1,C2,dt,i,j,b,q)
            

    #Neumann condition at the bouandary
    u[0,0] = u[1,1]
    u[0,Ny-1] = u[1,Ny-2]
    u[Nx-1,0] = u[Nx-2,1]
    u[Nx-1,Ny-1] = u[Nx-2,Ny-2]
    
    #Neumann for the half the bouandary
    for i in range(1,Nx-1):
        u[i,0],u[i,Ny-1] = step1_neuman_y(u_1,x,y,t,f,V,C1,C2,dt,i,Ny,b,q)
           
    #Neumann for the other half of the bouandary 
    for j in range(1,Ny-1):
        u[0,j],u[Nx-1,j] = step1_neuman_x(u_1,x,y,t,f,V,C1,C2,dt,j,Nx,b,q)
    
    #fixing u_1 and u_2 for next timestep.
    for i in range(Nx):
            for j in range(Ny):                
                u_2[i,j]=u_1[i,j]
                u_1[i,j]=u[i,j]
    
    for n in range(1,Nt):
        
        #finding u for t=dt*(n+1) at inner points
        for i in range(1,Nx-1):
            for j in range(1,Ny-1):
                u[i,j] = step(u_1,u_2,x,y,t,f,C1,C2,dt,i,j,n,b,q)
        
        #fixing edgy boundary points with Neumann
        u[0,0] = u[1,1]
        u[0,Ny-1] = u[1,Ny-2]
        u[Nx-1,0] = u[Nx-2,1]
        u[Nx-1,Ny-1] = u[Nx-2,Ny-2]
        
        #half of boundary with Neumann
        for i in range(1,Nx-1):
            u[i,0],u[i,Ny-1]=step_neuman_y(u_1,u_2,x,y,t,f,C1,C2,dt,i,n,Ny,b,q)
        
        #other half of boundary with Neumann
        for j in range(1,Ny-1):
            u[0,j],u[Nx-1,j]=step_neuman_x(u_1,u_2,x,y,t,f,C1,C2,dt,j,n,Nx,b,q)
            

        #fixing u_1 and u_2 for next timestep.
        for i in range(Nx):
            for j in range(Ny):                
                u_2[i,j]=u_1[i,j]
                u_1[i,j]=u[i,j]



    return u




################################################################
################################################################
################################################################
################################################################

def step1(u_1,x,y,t,f,V,C1,C2,dt,i,j,b,q):

    A = (1+0.5*b*dt)*u_1[i,j] + dt*V(x[i],y[j]) + 0.5*dt**2*f(x[i],y[j],t[0])
    
    qx1 = 0.5*(q(x[i],y[j]) + q(x[i+1],y[j]))  
    qx2 = 0.5*(q(x[i],y[j]) + q(x[i-1],y[j]))

    B = 0.5*C1*(qx1*(u_1[i+1,j] - u_1[i,j]) - qx2*(u_1[i,j] - u_1[i-1,j])) 

    qy1 = 0.5*(q(x[i],y[j]) + q(x[i],y[j+1]))  
    qy2 = 0.5*(q(x[i],y[j]) + q(x[i],y[j-1]))
    
    D = 0.5*C2*(qy1*(u_1[i,j+1] - u_1[i,j]) - qy2*(u_1[i,j] - u_1[i,j-1]))

    return (A + B + D)/( 1 + 0.5*b*dt ) 

def step1_neuman_y(u_1,x,y,t,f,V,C1,C2,dt,i,N,b,q):

    A1 = (1+0.5*b*dt)*u_1[i,0]+dt*V(x[i],y[0])+ 0.5*dt**2*f(x[i],y[0],t[0])
    A2 = (1+0.5*b*dt)*u_1[i,N-1]+dt*V(x[i],y[N-1])+0.5*dt**2*f(x[i],y[N-1],t[0])
    
    q11 = 0.5*(q(x[i],y[0]) + q(x[i+1],y[0]))  
    q12 = 0.5*(q(x[i],y[0]) + q(x[i-1],y[0]))

    q21 = 0.5*(q(x[i],y[N-1]) + q(x[i+1],y[N-1]))  
    q22 = 0.5*(q(x[i],y[N-1]) + q(x[i-1],y[N-1]))
    
    B1 = 0.5*C1*(q11*(u_1[i+1,0] - u_1[i,0]) - q12*(u_1[i,0] - u_1[i-1,0]))
    B2 = 0.5*C1*(q21*(u_1[i+1,N-1]-u_1[i,N-1])-q22*(u_1[i,N-1]-u_1[i-1,N-1]))

    D1 = C2*q(x[i],y[0])*(u_1[i,1] - u_1[i,0]) 
    D2 = C2*q(x[i],y[N-1])*(u_1[i,N-2] - u_1[i,N-1]) 

    return (A1+B1+D1)/( 1 + 0.5*b*dt ), (A2+B2+D2)/( 1 + 0.5*b*dt )

def step1_neuman_x(u_1,x,y,t,f,V,C1,C2,dt,j,N,b,q):

    A1 = (1+0.5*b*dt)*u_1[0,j]+dt*V(x[0],y[j])+0.5*dt**2*f(x[0],y[j],t[0])
    A2 = (1+0.5*b*dt)*u_1[N-1,j]+dt*V(x[N-1],y[j])+0.5*dt**2*f(x[N-1],y[j],t[0])
    
    B1 = C1*q(x[0],y[j])*(u_1[1,j] - u_1[0,j]) 
    B2 = C1*q(x[N-1],y[j])*(u_1[N-2,j] - u_1[N-1,j]) 

    q11 = 0.5*(q(x[0],y[j]) + q(x[0],y[j+1]))  
    q12 = 0.5*(q(x[0],y[j]) + q(x[0],y[j-1]))

    q21 = 0.5*(q(x[N-1],y[j]) + q(x[N-1],y[j+1]))  
    q22 = 0.5*(q(x[N-1],y[j]) + q(x[N-1],y[j-1]))

    D1 = 0.5*C2*(q11*(u_1[0,j+1] - u_1[0,j]) - q12*(u_1[0,j] - u_1[0,j-1]))
    D2 = 0.5*C2*(q21*(u_1[N-1,j+1] - u_1[N-1,j])-q22*(u_1[N-1,j]-u_1[N-1,j-1]))

    return (A1+B1+D1)/( 1 + 0.5*b*dt ), (A2+B2+D2)/( 1 + 0.5*b*dt )

def step(u_1,u_2,x,y,t,f,C1,C2,dt,i,j,n,b,q):
    A = (2+ b*dt)*u_1[i,j] - u_2[i,j] + dt**2*f(x[i],y[j],t[n])

    qx1 = 0.5*(q(x[i],y[j]) + q(x[i+1],y[j]))  
    qx2 = 0.5*(q(x[i],y[j]) + q(x[i-1],y[j]))

    B = C1*(qx1*(u_1[i+1,j] - u_1[i,j]) - qx2*(u_1[i,j] - u_1[i-1,j]))
    
    qy1 = 0.5*(q(x[i],y[j]) + q(x[i],y[j+1]))  
    qy2 = 0.5*(q(x[i],y[j]) + q(x[i],y[j-1]))
    
    D = C2*(qy1*(u_1[i,j+1] - u_1[i,j]) -qy2*(u_1[i,j] - u_1[i,j-1]))
    
    return (A + B + D)/(1 + b*dt)

def step_neuman_y(u_1,u_2,x,y,t,f,C1,C2,dt,i,n,N,b,q):
    A1 = (2 + b*dt)*u_1[i,0] - u_2[i,0] + dt**2*f(x[i],y[0],t[n])
    A2 = (2 + b*dt)*u_1[i,N-1] - u_2[i,N-1] + dt**2*f(x[i],y[N-1],t[n])
    
    q11 = 0.5*(q(x[i],y[0]) + q(x[i+1],y[0]))  
    q12 = 0.5*(q(x[i],y[0]) + q(x[i-1],y[0]))

    q21 = 0.5*(q(x[i],y[N-1]) + q(x[i+1],y[N-1]))  
    q22 = 0.5*(q(x[i],y[N-1]) + q(x[i-1],y[N-1]))

    B1 = C1*(q11*(u_1[i+1,0] - u_1[i,0]) - q12*(u_1[i,0] - u_1[i-1,0]))
    B2  = C1*(q21*(u_1[i+1,N-1]-u_1[i,N-1]) - q22*(u_1[i,N-1]-u_1[i-1,N-1]))

    D1 = 2*C2*q(x[i],y[0])*(u_1[i,1] - u_1[i,0])
    D2 = 2*C2*q(x[i],y[N-1])*(u_1[i,N-2] - u_1[i,N-1])

    return (A1+B1+D1)/(1 + b*dt), (A2+B2+D2)/(1 + b*dt)

def step_neuman_x(u_1,u_2,x,y,t,f,C1,C2,dt,j,n,N,b,q):

    A1 = (2 + b*dt)*u_1[0,j] - u_2[0,j]+ dt**2*f(x[0],y[j],t[n])
    A2 = (2 + b*dt)*u_1[N-1,j] - u_2[N-1,j] + dt**2*f(x[N-1],y[j],t[n])
    
    B1 = 2*C1*q(x[0],y[j])*(u_1[1,j] - u_1[0,j])
    B2 = 2*C1*q(x[N-1],y[j])*(u_1[N-2,j] - u_1[N-1,j])

    q11 = 0.5*(q(x[0],y[j]) + q(x[0],y[j+1]))  
    q12 = 0.5*(q(x[0],y[j]) + q(x[0],y[j-1]))

    q21 = 0.5*(q(x[N-1],y[j]) + q(x[N-1],y[j+1]))  
    q22 = 0.5*(q(x[N-1],y[j]) + q(x[N-1],y[j-1]))

    D1 = C2*(q11*(u_1[0,j+1]-u_1[0,j])-q12*(u_1[0,j]-u_1[0,j-1]))
    D2 = C2*(q21*(u_1[N-1,j+1]-u_1[N-1,j])-q22*(u_1[N-1,j]-u_1[N-1,j-1]))

    return (A1+B1+D1)/(1 + b*dt), (A2+B2+D2)/(1 + b*dt)




if __name__ == "__main__":
    
    I = lambda x,y: x+y
    V = lambda x,y: x
    f = lambda x,y,t: 0
    q = lambda x,y: 1
    print solver(I,V,f,q,0,1,1,0.2,1,1)
