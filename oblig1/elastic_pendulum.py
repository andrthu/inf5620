from numpy import *
import matplotlib.pyplot as plt

from matplotlib import animation

 
def simulate(beta=0.9,Theta=30,epsilon=0,num_periods=6,
             time_steps_per_period=60, plot=True):
    """
    Solve x'' = -C(x,y)*x and y'' = -C(x,y)(y - 1) - beta for t
    in (0,num_periods*P], x(0)=(1+epsilon)*sin(Theta), x'(0)=0,
    y(0)=1 - (1 + epsilon)*cos(Theta) and y'(0)=0 by a central 
    finite difference method with time step dt=P/time_steps_per_period.
    
    Where C(x,y)=beta/(1 - beta)*(1 - beta/L(x,y)), 
    L(x,y) = sqrt(x**2 + (y-1)**2) and P=2*pi.
    """


    n = num_periods * time_steps_per_period  #find number of timesteps
    dt = float(2*pi)/time_steps_per_period   #find timestep
    rad_Theta = pi*Theta/180.                #convert degrees to radians 
    
    #Variable to make calculations simpler
    K = beta/(1 - beta)                      
                                             

    #Create arrays that stores solution
    x=zeros(n+1)
    y=zeros(n+1)
    tet=zeros(n+1)
    t=linspace(0,dt*n,n+1)
    
    
    #Insert first initial condition. Sine and cosine functions
    #take radians as argument.
    x[0] = (1+epsilon)*sin(rad_Theta)
    y[0] = 1 - (1 + epsilon)*cos(rad_Theta)
    
    
    #calculate "euclidian norm" of (x[0],y[0]-1), and make 
    #new temporary help variable K2, to make calculations
    #"cleaner".
    L = sqrt(x[0]**2 + (y[0] - 1)**2)
    K2 = (1 - beta/L)*K
    
    #Find first timestep by useing second initial condition
    #x'(0)=y'(0)=0. centered finite diffrence gives exspression
    #for "x[-1]" and "y[-1]". Use theese in scheme.
    x[1] = x[0]*(1 - 0.5*K2*(dt**2))
    y[1] = y[0] - (y[0] - 1)*K2*0.5*dt**2 - 0.5*beta*dt**2


    #Loop that implements scheme
    for i in range(1,n):
        #calculate "norm" of solution at t=dt*i, and make
        #helpvariable K2, to make calculations "cleaner"
        L = sqrt(x[i]**2+(y[i] - 1)**2)
        K2 = (1-beta/L)*K


        #implementation of scheme
        x[i+1] = x[i]*(2 - K2*dt**2) - x[i-1]
        y[i+1] = 2*y[i] -(y[i]-1)*K2*dt**2 -y[i-1] - beta*dt**2
        
    
    tet=arctan(x/(1-y))

    if plot:
        plt.figure()
        plt.gca().set_aspect('equal')
        plt.plot(x,y)
        plt.title('Motion of pendulum')
        plt.xlabel('x')
        plt.ylabel('y')
        
        plt.figure()
        plt.plot(t, tet*180/pi)
        plt.title('Angle in degrees')
        plt.xlabel('t')
        plt.ylabel('tet')

        plt.show()
        
        if Theta<10:
            non_elastic = Theta*cos(t)
            plt.figure()
            plt.plot(t, tet*180/pi, '-b')
            plt.plot(t, non_elastic, 'r')
            plt.show()
    return x,y,tet,t



def demo(beta,theta):
    x,y,t,tet = simulate(beta, theta, 0, 3, 600, True)


    
    

if __name__=="__main__":
    #x1,y1,t1,tet=simulate()
    #test_simulate()
    #test_vertical_motion(0.9,0.9)
    demo(0.9,30)
