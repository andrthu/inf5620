import vectorization as vec 
from scitools.std import linspace, newaxis,zeros,amax

def solver(I,V,f,q,b,Lx,Ly,dt,T,C,mode='scalar',bug=None,ue=None):
    
    Nt = int(round(T/dt))
    t = linspace(0, Nt*dt, Nt+1)
    dt = t[1]-t[0]

    dx = dt/C 
    Nx = int(round(Lx/dx))
    x = linspace(0, Nx*dx, Nx+1)
    dx=x[1]-x[0]
    xv = x[:,newaxis]
    
    dy = dt/C 
    Ny = int(round(Ly/dy))
    y = linspace(0, Ny*dx, Ny+1)
    dy=y[1]-y[0]
    yv = y[newaxis,:]

    C1 = (C)**2
    C2 = (C)**2
    
    

    u = zeros(shape=(Nx+1,Ny+1))
    u_1 = zeros(shape=(Nx+1,Ny+1))
    u_2 = zeros(shape=(Nx+1,Ny+1))

    #handeling zero sorce term and initial condition
    if f is None or f == 0 :
        f = lambda x,y,t: 0
    if V is None or V == 0:
        V = lambda x,y: 0

    
    
    
      
    if mode=='scalar':
        
        
        #The first initial condition
        
        for i in range(Nx+1):
            for j in range(Ny+1):
                u_1[i,j]=I(dx*i,dy*j)
        
        
        ### BUG ### Making a bug for exercise 3,1,4
        if bug=='bug1':
            u_1 = zeros(shape=(Nx+1,Ny+1))
            for i in range(1,Nx+1):
                for j in range(1,Ny+1):
                    u_1[i,j]=I(dx*i,dy*j)
        
        
        #The second initial condiyion at inner points
        for i in range(1,Nx):
            for j in range(1,Ny):
                u[i,j] = step1(u_1,x,y,t,f,V,C1,C2,dt,i,j,b,q)
        
        ### BUG ### Making a bug for exercise 3,1,4
        if bug=='bug2':
            from constructed_bugs import bug_2
            for i in range(1,Nx-1):
                for j in range(1,Ny-1):
                    u[i,j] = bug_2(u_1,x,y,t,f,V,C1,C2,dt,i,j,b,q)

        ### BUG ### Making a bug for exercise 3,1,4
        if bug=='bug3':
            from constructed_bugs import bug_3
            for i in range(1,Nx):
                for j in range(1,Ny):
                    u[i,j] = bug_3(u_1,x,y,t,f,V,C1,C2,dt,i,j,b,q)

        ### BUG ### Making a bug for exercise 3,1,4
        if bug=='bug4':
            from constructed_bugs import bug_4
            for i in range(1,Nx):
                for j in range(1,Ny):
                    u[i,j] = bug_4(u_1,x,y,t,f,V,C1,C2,dt,i,j,b,q)

        #Neumann condition at the bouandary
        u[0,0] = u_1[1,1]
        u[0,Ny] =u_1[1,Ny-1]
        u[Nx,0] = u_1[Nx-1,1]
        u[Nx,Ny] = u_1[Nx-1,Ny-1]
        
        #Neumann for the half the bouandary
        for i in range(1,Nx):
           
            u[i,0],u[i,Ny] = step1_neuman_y(u_1,x,y,t,f,V,C1,C2,dt,i,Ny+1,b,q)
        
        #Neumann for the other half of the bouandary 
        for j in range(1,Ny):
            
            u[0,j],u[Nx,j] = step1_neuman_x(u_1,x,y,t,f,V,C1,C2,dt,j,Nx+1,b,q)
        
        ### BUG ### Making a bug for exercise 3,1,4
        if bug=='bug5':
            for j in range(1,Ny):
                #use current timestep to find the boundary instead
                #of useing previous time step.
                u[0,j],u[Nx,j]=step1_neuman_x(u,x,y,t,f,V,C1,C2,dt,j,Nx+1,b,q)
                
        
        print u  
        print 0
                
        #fixing u_1 and u_2 for next timestep.
        for i in range(Nx+1):
            for j in range(Ny+1):                
                u_2[i,j]=u_1[i,j]
                u_1[i,j]=u[i,j]



    if mode == 'vector':
        u_1[:,:]=I(xv[:],yv[:])
        #print u_1
        error=0
        u=vec.step1(u,u_1,xv,yv,t,f,V,C1,C2,Nx,Ny,dt,b,q)
        if ue!=None:
            err=amax(abs(u-ue(xv[:],yv[:],t[1])))
            if err>error:
                error =err
        
        u_2[:,:]=u_1[:,:]
        u_1[:,:]=u[:,:]
        #print u
        #print 0
        for n in range(1,Nt):
            u=vec.step(u,u_1,u_2,xv,yv,t,f,C1,C2,Nx,Ny,dt,b,q,n)
            u_2[:,:]=u_1[:,:]
            u_1[:,:]=u[:,:]
            #print u
            #print t[n]

            if ue!=None:
                err=amax(abs(u-ue(xv[:],yv[:],t[n+1])))
                if err>error:
                    error =err
        return u,xv,yv,t,error


    for n in range(1,Nt):
        
        #finding u for t=dt*(n+1) at inner points
        for i in range(1,Nx):
            for j in range(1,Ny):
                u[i,j] = step(u_1,u_2,x,y,t,f,C1,C2,dt,i,j,n,b,q)
        
        #fixing edgy boundary points with Neumann
        u[0,0] = u_1[1,1]
        u[0,Ny] = u_1[1,Ny-1]
        u[Nx,0] = u_1[Nx-1,1]
        u[Nx,Ny] = u_1[Nx-1,Ny-1]
        
        #half of boundary with Neumann
        for i in range(1,Nx):           
            u[i,0],u[i,Ny]=step_neuman_y(u_1,u_2,x,y,t,f,C1,C2,dt,i,n,Ny+1,b,q)
            
        #other half of boundary with Neumann
        for j in range(1,Ny):           
            u[0,j],u[Nx,j]=step_neuman_x(u_1,u_2,x,y,t,f,C1,C2,dt,j,n,Nx+1,b,q)
           

        #fixing u_1 and u_2 for next timestep.
        for i in range(Nx+1):
            for j in range(Ny+1):                
                u_2[i,j]=u_1[i,j]
                u_1[i,j]=u[i,j]


        #print u
        #print t[n]
    return u,xv,yv,t




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
    V = lambda x,y: x +0*y
    f = lambda x,y,t: 0
    q = lambda x,y: 1

    b=0
    Lx=1
    Ly=1
    dt=0.2
    T=1
    C=1

    mode='scalar'
    u,x,y,t=solver(I,V,f,q,b,Lx,Ly,dt,T,C,mode)
    print u
    print 


    mode='vector'
    u,x,y,t,error=solver(I,V,f,q,b,Lx,Ly,dt,T,C,mode)
    print u
    print 

"""
solver(I,V,f,q,b,Lx,Ly,dt,T,C,mode='scalar',bug=None)
"""
