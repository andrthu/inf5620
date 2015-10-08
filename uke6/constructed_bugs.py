

def bug_2(u_1,x,y,t,f,V,C1,C2,dt,i,j,b,q):

    A = (1+0.5*b*dt)*u_1[i,j] + dt*V(x[i],y[j]) + 0.5*dt**2*f(x[i],y[j],t[0])
    
    qx1 = 0.5*(q(x[i],y[j]) + q(x[i+1],y[j]))  
    qx2 = 0.5*(q(x[i],y[j]) + q(x[i-1],y[j]))
    
    B = 0.5*C1*(qx1*(u_1[i+1,j] - u_1[i,j]) - qx2*(u_1[i,j] - u_1[i-1,j])) 

    qy1 = 0.5*(q(x[i],y[j]) + q(x[i],y[j+1]))  
    qy2 = 0.5*(q(x[i],y[j]) + q(x[i],y[j-1]))
    
    ### BUG ### changed a sign         under here
    D = 0.5*C2*(qy1*(u_1[i,j+1] - u_1[i,j]) + qy2*(u_1[i,j] - u_1[i,j-1]))
    
    return (A + B + D)/( 1 + 0.5*b*dt ) 

def bug_3(u_1,x,y,t,f,V,C1,C2,dt,i,j,b,q):

    A = (1+0.5*b*dt)*u_1[i,j] + dt*V(x[i],y[j]) + 0.5*dt**2*f(x[i],y[j],t[0])
    
    qx1 = 0.5*(q(x[i],y[j]) + q(x[i+1],y[j]))  
    qx2 = 0.5*(q(x[i],y[j]) + q(x[i-1],y[j]))
    
    B = 0.5*C1*(qx1*(u_1[i+1,j] - u_1[i,j]) - qx2*(u_1[i,j] - u_1[i-1,j])) 

    qy1 = 0.5*(q(x[i],y[j]) + q(x[i],y[j+1]))  
    qy2 = 0.5*(q(x[i],y[j]) + q(x[i],y[j-1]))

    ### BUG ### changed a sign                         under here
    D = 0.5*C2*(qy1*(u_1[i,j+1] - u_1[i,j]) - qy2*(u_1[i,j] + u_1[i,j-1]))
    
    return (A + B + D)/( 1 + 0.5*b*dt ) 

def bug_4(u_1,x,y,t,f,V,C1,C2,dt,i,j,b,q):

    A = (1+0.5*b*dt)*u_1[i,j] + dt*V(x[i],y[j]) + 0.5*dt**2*f(x[i],y[j],t[0])
    
    ### BUG ### change the evaluated q functions 
    qx1 = 0.5*(q(x[i],y[j]) + q(x[i+1],y[j])) + 28 
    qx2 = 0.5*(q(x[i],y[j]) + q(x[i-1],y[j])) - 73
    
    B = 0.5*C1*(qx1*(u_1[i+1,j] - u_1[i,j]) - qx2*(u_1[i,j] - u_1[i-1,j]))
    
    ### BUG ### change the evaluated q functions 
    qy1 = 0.5*(q(x[i],y[j]) + q(x[i],y[j+1])) -0.5 
    qy2 = 0.5*(q(x[i],y[j]) + q(x[i],y[j-1])) +10.43
    
    D = 0.5*C2*(qy1*(u_1[i,j+1] - u_1[i,j]) - qy2*(u_1[i,j] - u_1[i,j-1]))
    
    return (A + B + D)/( 1 + 0.5*b*dt ) 
