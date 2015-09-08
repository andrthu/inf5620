from numpy import *
from numpy.random import *
from vib_undamped_verify_mms import solver



def test_quadratic(tol = 10**(-10)):
    """
    Testing that constructed solution q(x)=x**2+x+1,
    is computed correctly by solver function.

    Tollerence to fit test.
    """
    #Set up some T,dt and w values for testing. w values
    #choosen so that stabilety requierement dt<2/w holds 
    #for all dt.
    T_list=[1, 4, 10, 20]
    dt_list=[0.5, 0.2, 0.1]
    W_list=[0.2, 0.4, 0.8, 1.2]

    #creating quadratic function
    q=lambda x: x**2+x+1    
       
    #looping over test values
    for T in T_list:
        for dt in dt_list:
            for w in W_list:  

                #creating sorce term
                f=lambda x : 2 + (w**2)*q(x)
            
                #solve equation with solver
                u,t=solver(1,w,1,dt,T,f)

                #check if maximal diffrence between exact and
                #numerical solution less than tlerance.
                assert abs(amax(u-q(t)))<tol

        
