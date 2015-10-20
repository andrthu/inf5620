import vectorization as vec 
from numpy import linspace, newaxis,zeros,amax,sqrt,ceil
from scalar import *



"""
2D wave equation solved by finite differences::
  u,xv,yv,t = solver(I,V,f,q,b,Lx,Ly,dt,T,C,q_max,
                       mode='scalar',bug=None,ue=None):
Solve the 2D wave equation u_tt + bu_t = (qu_x)x + (qu_y)y + f(x,y,t) 
on (0,Lx)x(0,Ly) with
du/dn=0 on the boundary and initial condition du/dt=V(x,y).
The mesh points are numbered as (0,0), (1,0), (2,0),
..., (Nx,0), (0,1), (1,1), ..., (Nx, Ny).
dt is the time step. 
T is the stop time for the simulation.
I, V, f are functions: I(x,y), V(x,y), f(x,y,t). V and f
can be specified as None or 0, resulting in V=0 and f=0.

"""
def solver(I,V,f,q,b,Lx,Ly,dt,T,C,q_max,mode='scalar',bug=None,ue=None,
           ani=False):

    
    
    Nt = int(round(T/dt))             #Number of timepoints
    t = linspace(0, Nt*dt, Nt+1)      #t values
    dt = t[1]-t[0]                    #new dt to fit better  
    
    c_max = sqrt(q_max)               #c value for stabilety
    
    
    Nx = int(ceil((Lx*C)/(dt*c_max))) #number of points in x-dir 
    x = linspace(0,Lx,Nx+1)           #x values  
    dx=x[1]-x[0]                      #exact dx
    xv = x[:,newaxis]                 #vectorized axis
    
     
    Ny = int(ceil((Ly*C)/(dt*c_max))) #number of points in y-dir
    y = linspace(0,Ly,Ny+1)           #y values
    dy=y[1]-y[0]                      #exact dy
    yv = y[newaxis,:]                 #vectorized axis
    
    #help variables
    C1 = (dt/dx)**2 
    C2 = (dt/dy)**2
    
    error=0
    
    #Define matrixes to store solution
    u = zeros(shape=(Nx+1,Ny+1))
    u_1 = zeros(shape=(Nx+1,Ny+1))
    u_2 = zeros(shape=(Nx+1,Ny+1))

    #handeling zero sorce term and initial condition
    if f is None or f == 0 :
        f = lambda x,y,t: 0
    if V is None or V == 0:
        V = lambda x,y: 0

    
    
    
    #scalar version
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
        #Handle edges diffrent in vectorized code.
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
                
        
        
                
        #fixing u_1 and u_2 for next timestep.
        for i in range(Nx+1):
            for j in range(Ny+1):                
                u_2[i,j]=u_1[i,j]
                u_1[i,j]=u[i,j]

        #time loop for scalar version
        for n in range(1,Nt):
        
            #finding u for t=dt*(n+1) at inner points
            for i in range(1,Nx):
                for j in range(1,Ny):
                    u[i,j] = step(u_1,u_2,x,y,t,f,C1,C2,dt,i,j,n,b,q)
        
            #fixing edgy boundary points with Neumann
            #Handle edges diffrent in vectorized code.
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



    #vectorized version
    if mode == 'vector':
        #animation variable
        L=[]


        #initial condition
        u_1[:,:]=I(xv[:],yv[:])
        if ani==True:
            L.append(u_1)

        #calculating error variable
        
        #first step
        u=vec.step1(u,u_1,xv,yv,t,f,V,C1,C2,Nx,Ny,dt,b,q,x,y)
        
        if ani==True:
            L.append(u)
            
        #calculating max norm error at timelevel 1
        if ue!=None:
            err=amax(abs(u-ue(xv[:],yv[:],t[1])))
            if err>error:
                error =err
        
        #setting matrixes for next step. 
        u_2[:,:]=u_1[:,:]
        u_1[:,:]=u[:,:]

        
            
            

        #time loop
        for n in range(1,Nt):
            u=vec.step(u,u_1,u_2,xv,yv,t,f,C1,C2,Nx,Ny,dt,b,q,n,x,y)
            u_2[:,:]=u_1[:,:]
            u_1[:,:]=u[:,:]
            
            if ani==True:
                L.append(u)

            if ue!=None:
                err=amax(abs(u-ue(xv[:],yv[:],t[n+1])))
                if err>error:
                    error =err

        
            
        if ani== True:
            return L,x,y,t

        
        
    return u,xv,yv,t,error




