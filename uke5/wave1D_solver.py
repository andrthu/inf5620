from numpy import *
from wave_standing import *

def solver(I, V,f, q, L, dt, C, T, max_q,
           user_action=None,disc='type57'):

    Nt = int(round(T/dt))
    t = linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = dt*max_q/float(C)
    Nx = int(round(L/dx))
    x = linspace(0, L, Nx+1)       # Mesh points in space
    C3 =(dt/dx)**2                     # Help variable in the scheme
    

    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0

    u   = zeros(Nx+1)   # Solution array at new time level
    u_1 = zeros(Nx+1)   # Solution at 1 time level back
    u_2 = zeros(Nx+1)   # Solution at 2 time levels back


    import time;  t0 = time.clock()  # for measuring CPU time

    

    # Load initial condition into u_1
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Special formula for first time step
    n = 0
    for i in range(1, Nx):

        K = q(dx*(i+0.5))*(u_1[i+1]-u_1[i]) - q(dx*(i-0.5))*(u_1[i]-u_1[i-1])
        
        u[i] = u_1[i] + dt*V(x[i]) + 0.5*C3*K + 0.5*dt**2*f(x[i], t[n])
    
    if disc == 'type54':
        K0=q(0)*(u_1[1]-u_1[0])   
        KN=q(dx*Nx)*(u_1[Nx-1]-u_1[Nx])
        
    if disc == 'type57':
        K0=q(dx*0.5)*(u_1[1]-u_1[0])   
        KN=q(dx*(Nx-0.5))*(u_1[Nx-1]-u_1[Nx])

    u[0] = u_1[0] + dt*V(x[0]) + C3*K0 + 0.5*dt**2*f(x[0], t[n])
    u[Nx] = u_1[Nx] + dt*V(x[Nx]) + C3*KN + 0.5*dt**2*f(x[Nx], t[n])

    if user_action is not None:
        user_action(u, x, t, 1)

    # Switch variables before next step
    u_2[:] = u_1;  u_1[:] = u

    for n in range(1, Nt):
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            K = q(dx*(i+0.5))*(u_1[i+1]-u_1[i])-q(dx*(i-0.5))*(u_1[i]-u_1[i-1])
            
            u[i] = - u_2[i] + 2*u_1[i] + C3*K + dt**2*f(x[i], t[n])
            
        # Insert boundary conditions

        if disc == 'type57':
            K0=2*q(dx*0.5)*(u_1[1]-u_1[0])  
            KN=2*q(dx*(Nx-0.5))*(u_1[Nx-1]-u_1[Nx])
            
        if disc == 'type54':
            K0=2*q(0)*(u_1[1]-u_1[0])  
            KN=2*q(dx*Nx)*(u_1[Nx-1]-u_1[Nx])

        print dt**2*f(x[0], t[n])
        u[0] = - u_2[0] + 2*u_1[0] + C3*K0 + dt**2*f(x[0], t[n])
        u[Nx] = - u_2[Nx] + 2*u_1[Nx] + C3*KN + dt**2*f(x[Nx], t[n])
        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Switch variables before next step
        u_2[:] = u_1;  u_1[:] = u

    cpu_time = t0 - time.clock()
    return u, x, t, cpu_time

def make_f(q,d_q,dd_q,L):
    
    omega = lambda x: sqrt(q(x))*pi/L
    d_omega = lambda x: pi*d_q(x)/(2.*L*sqrt(q(x)))
    dd_omega = lambda x: dd_q(x)*d_omega(x)-(pi*d_q(x)**2)/(4*L*sqrt(q(x))**3)
    
    ue = lambda x,t: cos(pi*x/L)*cos(t*omega(x))
    help1_f = lambda x,t: sin(pi*x/L)*cos(t*omega(x))
    help2_f = lambda x,t: cos(pi*x/L)*sin(t*omega(x))
    help3_f = lambda x,t: sin(pi*x/L)*sin(t*omega(x))

    u_tt = lambda x,t: -ue(x,t)*omega(x)**2

    u_x = lambda x,t: -pi*help1_f(x,t)/L - t*d_omega(x)*help2_f(x,t)/L
    u_xx1 = lambda x,t: -ue(x,t)*((pi/L)**2 + (t*d_omega(x))**2)
    u_xx2 = lambda x,t: help3_f(x,t)*2*t*pi*d_omega(x)/L
    u_xx3 = lambda x,t: -t*dd_omega(x)*help2_f(x,t)
    u_xx = lambda x,t: u_xx1(x,t)+u_xx2(x,t)+u_xx3(x,t)
    
    return lambda x,t: u_tt(x,t) - d_q(x)*u_x(x,t)-q(x)*u_xx(x,t),ue


if __name__ == "__main__":
    
    L=1
    C=1.1
    T=2
    dt=0.01
    I = lambda x: cos(pi*x/L)
    disc='type57'

    q = lambda x: 1+ (x-L/2.)**4
    d_q= lambda x: 4*(x-L/2.)**3
    dd_q=lambda x: 12*(x-L/2.)**2
    max_q = q(0)

    omega = lambda x: sqrt(q(x))*pi/L
    d_omega = lambda x: pi*d_q(x)/(2.*L*sqrt(q(x)))
    dd_omega = lambda x: dd_q(x)*d_omega(x)-(pi*d_q(x)**2)/(4*L*sqrt(q(x))**3)
    
    ue = lambda x,t: cos(pi*x/L)*cos(t*omega(x))
    help1_f = lambda x,t: sin(pi*x/L)*cos(t*omega(x))
    help2_f = lambda x,t: cos(pi*x/L)*sin(t*omega(x))
    help3_f = lambda x,t: sin(pi*x/L)*sin(t*omega(x))


    u_tt = lambda x,t: -ue(x,t)*omega(x)**2

    u_x = lambda x,t: -pi*help1_f(x,t)/L - t*d_omega(x)*help2_f(x,t)/L
    u_xx1 = lambda x,t: -ue(x,t)*((pi/L)**2 + (t*d_omega(x))**2)
    u_xx2 = lambda x,t: help3_f(x,t)*2*t*pi*d_omega(x)/L
    u_xx3 = lambda x,t: -t*dd_omega(x)*help2_f(x,t)
    u_xx = lambda x,t: u_xx1(x,t)+u_xx2(x,t)+u_xx3(x,t)

    #f = lambda x,t: u_tt(x,t) - d_q(x)*u_x(x,t)-q(x)*u_xx(x,t)
    f,ue=make_f(q,d_q,dd_q,L)
    
    """
    viz(I,None,f, q,L,dt,C,T,ue,max_q,
        -2,2,True,'matplotlib',solver,False,disc)
    """
    q = lambda x: (cos(pi*x/L))**2 +1
    d_q= lambda x: -pi*sin(2*pi*x/L)/L
    dd_q=lambda x: -2*cos(2*pi*x/L)*(pi/L)**2
    max_q = q(0)
    f,ue=make_f(q,d_q,dd_q,L)

    viz(I,None,f, q,L,dt,C,T,ue,max_q,
        -2,2,True,'matplotlib',solver,False,disc)
    
"""
    I, V, f, c, L, dt, C, T, ue, max_q,
    umin, umax,               
    animate=True,             
    tool='matplotlib',        
    solver_function=solver,   
    )
"""
