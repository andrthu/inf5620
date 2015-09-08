from numpy import *

import sympy as sym 

V, t, I, w, dt, b, a = sym.symbols('V t I w dt b a')  # global symbols
f = None  # global variable for the source term in the ODE

def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""

    return sym.diff(u(t), t, t) + w**2*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = DtDt(u,dt)  + w**2*u(t)- f #ode_source_term(u)
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""

    R = 2*(u(dt) - u(0) - dt*V)/(dt**2) + w**2*u(0) - f.subs(t, 0)

    return sym.simplify(R)

def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """

    return (u(t + dt) - 2*u(t) + u(t - dt))/(dt**2)

def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """

    print '=== Testing exact solution: %s ===' % u(t)
    print "Initial conditions u(0)=%s, u'(0)=%s:" % (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sym.simplify(ode_source_term(u))
    

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)

def linear():
    main(lambda t: V*t + I)

def quadratic():
    """
    function that run main with a quadratic polynomial. Notice
    that coeffichent b infront of t**2 is not determined and 
    could be anything
    """

    main(lambda t: b*t**2 + V*t + I)


def cubic():
    """
    function that run main with a cubic polynomial. Notice
    that coeffichent b infront of t**2 and a infront of t**3
    is not determined and could be anything. Reason for also 
    making a cubic function, is to solve exercise 1e).
    """

    main(lambda t: a*t**3 + b*t**2 +V *t + I)

def solver(I,w,V,dt,T,f):
    """
    Solve u'' + w**2*u = f for t in (0,T], u(0)=I and u'(0)=V,
    by a central finite difference method with time step dt.
    """

    n = int(round(T/dt))
    u = zeros(n+1)
    t = linspace(0,n*dt,n+1)
    
    #setting initial conditions u(0)=I and u'(0)=V. 
    #u'(0)=V with central diffrence. 
    u[0] = I
    u[1] = 0.5*(I*(2 - (dt*w)**2) + 2*V*dt + f(0)*dt**2)

    
    #Loop that implements scheme
    for i in range(1,n):
        u[i+1] = u[i]*(2 - (dt*w)**2) - u[i-1] + f(i*dt)*dt**2

    return u,t




if __name__ == '__main__':
    linear()
    quadratic()     
    cubic()
    


"""
terminal> python vib_undamped_verify_mms.py
=== Testing exact solution: I + V*t ===
Initial conditions u(0)=I, u'(0)=V:
residual step1: 0
residual: 0
=== Testing exact solution: I + V*t + b*t**2 ===
Initial conditions u(0)=I, u'(0)=V:
residual step1: 0
residual: 0
=== Testing exact solution: I + V*t + a*t**3 + b*t**2 ===
Initial conditions u(0)=I, u'(0)=V:
residual step1: 2*a*dt
residual: 0
"""
