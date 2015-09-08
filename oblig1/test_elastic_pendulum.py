from numpy import *
import matplotlib.pyplot as plt
from elastic_pendulum import simulate


def test_zero_sulotion():
    """
    Testing that Theta=epsilon=0 gives constant zero
    solutions when inserted in simulate function
    """

    #test for diffrent beta and time steps per period
    N=[50,100,1000]  

    #it only makes phissical sence if beta is in (0,1)
    B=[0.1,0.4,0.7,0.9]

    #loop through the different test values
    for n in N:
        for beta in B:
            
            #use simulate function we have made with arguments
            #given in exercise. Also useing number_of_periods=1
            x,y,tet,t=simulate(beta,0,0,1,n,False)
           
            #find maximal values in x and y solution
            x_max = abs(amax(x))
            y_max = abs(amax(y))
            
            #check that x_max and y_max is 0, with a possabillety
            #for rounding errors.
            assert ((x_max<10**(-16)) and (y_max<10**(-16)))


def test_vertical_motion(dt=2*pi*0.0001,tol=10**(-5)):
    """
    Testing that vertical motion corrosponds with exact solution.
    when Theta=0, epsilon in (0,1), we get pure vertical motion 
    satisfying y'' = -(beta/(1 - beta)y, which has exact solution
    y(t)=Icos(wt), I = -epsilon, and w = sqrt(beta/(1 - beta)). This
    is used to show that simulate is good.

    I set tollerance for my test case to be 10**(-5), this depends
    on size of dt.
    """

    
    
    #test for diffrent epsilons and betas.
    Beta_and_eps=[0.1,0.4,0.7,0.9]
    
    #find timsteps_per_period with dt
    n= int(round(2*pi/dt))

    #loop throgh diffrent testing values.
    for eps in Beta_and_eps:
        for beta in Beta_and_eps:
            
            #solve the equation
            x,y,tet,t = simulate(beta,0,eps,1,n,False)
            
            #create exact solution
            w = sqrt(beta/(1 - beta))
            I = -eps    
            y_exact = I*cos(w*t)

            
            #calculate l2 error of exact solution minus 
            #numerical solution.
            l2_error= sqrt(dt*sum((y-y_exact)**2))
            
            #check if l2 norm is less then tolerence.    
            assert l2_error<tol
