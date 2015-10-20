from numpy import *


def solver(I, V,f, q, L, dt, beta, T, max_q,
           user_action=None,disc='type57'):
    """solving equation u_tt=(qu_x)_x +f with Neumann boundary. """

    
    Nt = int(round(T/dt))
    t = linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = dt*max_q/float(beta)  # Find spacing in x direction so method is stable
    Nx = int(round(L/dx))      # Number of discretization points in x dim.
    x = linspace(0, L, Nx+1)   # Mesh points in space
    C3 =(dt/dx)**2             # Help variable in the scheme
    
    

    #handeling zero sorce term and initial condition
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0


    u   = zeros(Nx+1)   # Solution array at new time level
    u_1 = zeros(Nx+1)   # Solution at 1 time level back
    u_2 = zeros(Nx+1)  # Solution at 2 time levels back
    
   
    

    import time;  t0 = time.clock()  # for measuring CPU time

    
    
    # Load initial condition into u_1
    u_1[:] = I(x[:])
    
    
    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Special formula for first time step
    n = 0
    
    #Vectorized version of expression:
    #u[1,i] = u[0,i] + dt*V(xi) + 0.5*dt**2*f(xi,0) +
    #0.5*(dt/dx)**2(q(xi+dx/2)*(u[0,i+1]-u[0,i])-q(xi-dx/2)*(u[0,i]-u[0,i-1])
    u[1:-1]=u_1[1:-1] + dt*V(x[1:-1]) + 0.5*C3*(q(x[1:-1]+0.5*dx)*(u_1[2:]-u_1[1:-1])-q(x[1:-1]-0.5*dx)*(u_1[1:-1]-u_1[:-2])) + 0.5*dt**2*f(x[1:-1],t[n])
   

    
    #Here are four diffrent types of boundary discretization. 
    #type 54 and 57 are the ones you are asked to implement in exersise
    #13 a). type 3 is the one from c, and type 4 is the one from d. The 
    #diffrence between type54, 57 and 4 is the handeling of the 
    #discretization of (qu_x)_x, therefore for theese types I have made
    #3 if tests that give values to these terms, and then plug them into
    #a general expression that fit all types. For type3 the discretization is
    #handled diffrentley.

    #All handel the u-values u[1,0] and u[1,Nx].
    if disc == 'type54':
        K0=q(0)*(u_1[1]-u_1[0])   
        KN=q(dx*Nx)*(u_1[Nx-1]-u_1[Nx])
        
    if disc == 'type57':
        K0=q(dx*0.5)*(u_1[1]-u_1[0])   
        KN=q(dx*(Nx-0.5))*(u_1[Nx-1]-u_1[Nx])
    
    if disc == 'type4':
        K0=0.5*q(0.5*dx)*(u_1[1]-u_1[0])
        KN=0.5*q(dx*(Nx-0.5))*(u_1[-2]-u_1[-1])
        
        
    if disc !='type3':     
        u[0] = u_1[0] + dt*V(x[0]) + C3*K0 + 0.5*dt**2*f(x[0], t[n])
        u[Nx] = u_1[Nx] + dt*V(x[Nx]) + C3*KN + 0.5*dt**2*f(x[Nx], t[n])
        
    if disc == 'type3':
        u[0] = u[1]
        u[-1]=u[-2]

    
    

    if user_action is not None:
        user_action(u, x, t, 1)

    # Switch variables before next step
    u_2[:] = u_1;  u_1[:] = u
    

    #Time loop starts here. Each step finds u at time n*dt.
    #when loop is finished u holds values for u(x,T).
    for n in range(1, Nt):
        
        #Vectorized version of scheme for inner points. The expression is:
        #u[n+1,i] = 2u[n,i] - u[m-1,i] + 0.5*dt**2*f(xi,n*dt) +
        #(dt/dx)**2(q(xi+dx/2)*(u[n,i+1]-u[n,i])-q(xi-dx/2)*(u[n,i]-u[n,i-1])
        u[1:-1] = - u_2[1:-1] + 2*u_1[1:-1] + dt**2*f(x[1:-1],t[n]) + C3*(q(x[1:-1]+0.5*dx)*(u_1[2:]-u_1[1:-1])-q(x[1:-1]-0.5*dx)*(u_1[1:-1]-u_1[:-2]))
        
          
        
        # Insert boundary conditions
        # Same as above, but diffrent expressions for inner timesteps. 
        if disc == 'type57':
            K0 = 2*q(dx*0.5)*(u_1[1] - u_1[0])  
            KN = 2*q(dx*(Nx - 0.5))*(u_1[Nx-1] - u_1[Nx])
            
        if disc == 'type54':
            K0 = 2*q(0)*(u_1[1] - u_1[0])  
            KN = 2*q(dx*Nx)*(u_1[Nx-1]-u_1[Nx])

        if disc == 'type4':
            K0 = q(0.5*dx)*(u_1[1] - u_1[0])
            KN = q(dx*(Nx-0.5))*(u_1[-2] - u_1[-1])
            
            
        if disc !='type3':
            u[0] = - u_2[0] + 2*u_1[0] + C3*K0 + dt**2*f(x[0], t[n])
            u[Nx] = - u_2[Nx] + 2*u_1[Nx] + C3*KN + dt**2*f(x[Nx], t[n])

        if disc =='type3':
           u[0] = u[1]
           u[-1]=u[-2] 
            
        if user_action is not None:
            if user_action(u, x, t, n+1):
                break
            
       

        # Switch variables before next step
        u_2[:] = u_1;  u_1[:] = u

    cpu_time = t0 - time.clock()

    
    #return u at time T.
    return u,x,t,cpu_time

def test_convergence():
    """Testing solver with q=1+(x-L/2.)**4"""

    
    L=2                        # Size of space domain
    beta=0.9                   # beta constant for secureing stabilety
    T=3                        # Final time
    dt=0.1                     # Initial timestep 
    I = lambda x: cos(pi*x/L)  # Initial condition
    K=pi/L                     # help variable

    
    
    q = lambda x: 1 + (x-L/2.)**4     # Defineing q
    d_q= lambda x: 4*(x-L/2.)**3     # Defineing derivative of q 
    
    ue = lambda x,t:cos(K*x)*cos(t)  # Defineing exact solution

    #Defineing sorceterm so we analytacly get exatct solution with given q.
    f= lambda x,t:-ue(x,t) + q(x)*ue(x,t)*K**2 + d_q(x)*K*sin(K*x)*cos(t)

    #Looking at q we see that it takes maximum at x=0.
    max_q = sqrt(q(0))

    #Testing solver for diffrent dt values in a loop. For each step in the
    #loop I multiply step with previus dt value.
    step=0.9
    
    #number of steps in loop
    val =  10
    
    #Arrays for measureing error for diffrent types of discretizations.
    E1 = zeros(val)
    E2 = zeros(val)
    E3 = zeros(val)
    E4 = zeros(val)

    #solveing the equations for diffrent dt values, and calculate 
    #diffrence between u and exact solution at t=T.
    for i in range(val): 
        
        u1,x1,t,cp=solver(I, None,f, q, L, dt*step**(i), beta, T, max_q,
                                 user_action=None,disc='type54')

        u2,x2,t,c=solver(I, None,f, q, L, dt*step**(i), beta, T, max_q,
                                 user_action=None,disc='type57')

        u3,x3,t,c=solver(I, None,f, q, L, dt*step**(i), beta, T, max_q,
                                 user_action=None,disc='type3')

        u4,x4,t,c=solver(I, None,f, q, L, dt*step**(i), beta, T, max_q,
                                 user_action=None,disc='type4')

        E1[i]=amax(abs(u1-ue(x1,T)))
        E2[i]=amax(abs(u2-ue(x2,T)))
        E3[i]=amax(abs(u3-ue(x3,T)))
        E4[i]=amax(abs(u4-ue(x4,T)))


    #printing out avrage converge rates for the diffrent discretizaton
    #schemes.    
    print 'Average convergent rate q=1 + (x-L/2.)**4:'
    print sum(log(E1[1:]/E1[:-1])/log(step))/(val-1),'type54'
    print sum(log(E2[1:]/E2[:-1])/log(step))/(val-1),'typr57'
    print sum(log(E3[1:]/E3[:-1])/log(step))/(val-1),'type3'
    print sum(log(E4[1:]/E4[:-1])/log(step))/(val-1),'type4 \n'

def test_convergence2():
    """Testing solver with q=1+cos(pi*x/L)"""

    L=2                        # Size of space domain
    beta=0.9                   # beta constant for secureing stabilety
    T=2                        # Final time
    dt=0.1                     # Initial timestep 
    I = lambda x: cos(pi*x/L)  # Initial condition
    K=pi/L                     # help variable

    
    
    q = lambda x: 1 + cos(K*x)     # Defineing q
    d_q= lambda x: -K*sin(K*x)     # Defineing derivative of q 
    
    ue = lambda x,t:cos(K*x)*cos(t)   # Defineing exact solution


    #Defineing sorceterm so we analytacly get exatct solution with given q.
    f= lambda x,t:-ue(x,t) + q(x)*ue(x,t)*K**2 + d_q(x)*K*sin(K*x)*cos(t)

    #Looking at q we see that it takes maximum at x=0.
    max_q = sqrt(q(0))

    #Testing solver for diffrent dt values in a loop. For each step in the
    #loop I multiply step with previus dt value.
    step=0.9

    #number of steps in loop
    val = 10

    #Arrays for measureing error for diffrent types of discretizations.
    E1 = zeros(val)
    E2 = zeros(val)
    E3 = zeros(val)

    #solveing the equations for diffrent dt values, and calculate 
    #diffrence between u and exact solution at t=T.
    for i in range(val): 
        
        u1,x1,t,cp=solver(I, None,f, q, L, dt*step**(i), beta, T, max_q,
                                 user_action=None,disc='type57')

        u2,x2,t,c=solver(I, None,f, q, L, dt*step**(i), beta, T, max_q,
                                 user_action=None,disc='type3')

        u3,x3,t,c=solver(I, None,f, q, L, dt*step**(i), beta, T, max_q,
                                 user_action=None,disc='type4') 

        E1[i]=amax(abs(u1-ue(x1,T)))
        E2[i]=amax(abs(u2-ue(x2,T)))
        E3[i]=amax(abs(u3-ue(x3,T)))
        
        
    
    
    #printing out avrage converge rates for the diffrent discretizaton
    #schemes    
    print 'Average convergent rate, q=1 + cos(pi*x/L):'
    print sum(log(E1[1:]/E1[:-1])/log(step))/(val-1),'type57'
    print sum(log(E2[1:]/E2[:-1])/log(step))/(val-1),'type3'
    print sum(log(E3[1:]/E3[:-1])/log(step))/(val-1),'type4'






def test_constant():
    """Test that linear solution gives exact result. This is not asked for
    in exersise, so I will not comment it."""
    L=2
    beta=0.9
    T=1
    dt=0.01
    I = lambda x: 0*x+1
    K=pi/L 
    disc='type57'

    q = lambda x: 1 +(x-L/2.)**4
    d_q= lambda x: 4*(x-L/2.)**3

    ue = lambda x,t: t + 1 +0*x
    
    max_q = q(0)
    u1,x1,t,cp=solver(I, lambda x: 1,None, q, L, dt, beta, T, max_q,
                                 user_action=None,disc='type54')
    u2,x2,t,cp=solver(I, lambda x: 1,None, q, L, dt, beta, T, max_q,
                                 user_action=None,disc='type57')
    u3,x3,t,cp=solver(I, lambda x: 1,None, q, L, dt, beta, T, max_q,
                                 user_action=None,disc='type3')
    u4,x4,t,cp=solver(I, lambda x: 1,None, q, L, dt, beta, T, max_q,
                                 user_action=None,disc='type4')


    print "Error for linear solution:"
    print amax(abs(ue(x1,T)-u1)),'type54' 
    print amax(abs(ue(x2,T)-u2)),'type57'
    print amax(abs(ue(x3,T)-u3)),'type3'
    print amax(abs(ue(x4,T)-u4)),'type4 \n'

if __name__ == "__main__":
    # Do the tests
    test_constant()   
    test_convergence()
    test_convergence2()
    
#This is the result of running the code as it is with val = 10 in
#test_convergence() and test_convergence2().
"""
terminal> python wave1D_solver.py
Error for linear solution:
8.881784197e-16 type54
8.881784197e-16 type57
8.881784197e-16 type3
8.881784197e-16 type4 

Average convergent rate q=1 + (x-L/2.)**4:
1.90108082555 type54
1.90971602891 typr57
0.872485667178 type3
1.03560912869 type4 

Average convergent rate, q=1 + cos(pi*x/L):
1.34689507223 type57
0.821381710942 type3
1.05822464586 type4
"""

#This is the result of running the code with val = 50 in
#test_convergence() and test_convergence2().
"""
terminal> python wave1D_solver.py
Error for linear solution:
8.881784197e-16 type54
8.881784197e-16 type57
8.881784197e-16 type3
8.881784197e-16 type4 

Average convergent rate q=1 + (x-L/2.)**4:
0.944973346604 type54
1.06980641619 typr57
0.998304775278 type3
0.974988244364 type4 

Average convergent rate, q=1 + cos(pi*x/L):
0.949343099669 type57
0.986696692504 type3
0.990664966219 type4
"""
