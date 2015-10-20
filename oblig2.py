from numpy import *
from matplotlib.pylab import *
import scipy.linalg as sci
n=4

b1=zeros(n)

A=matrix('1 0 5 0; 0 1 0 0; 1 0 1 0; 0 0 0 1')

#regner ut euclidsk norm til vektoren (x,y)
def norm(x,y):
    return sqrt(sum(x**2+y**2))

#python version av koden fra boka. Kommenterer
#bare rotsolve.
def rothestri(A,b):
    n=len(b)
    
    for k in range(i):
        r = norm(A[k,k],A[k+1,k])
        if r>0:
            c = A[k,k]/r
            s = A[k+1,k]/r
            
            A[[k,k+1],k+1:n+1]=matrix([[c,s],[-s,c]])*A[[k,k+1],k+1:n+1]
        A[k,k]=r
        A[k+1,k]=0

    return sci.solve_triangular(A,b)      

#Loser systemet Ax=b ved hjelp av Givens rotasjoner. For
#generelle nxn matriser A.
def rotsolve(A,b):

    #finner dimensjonen til A
    n=len(b)
    
    #gaar gjennom hver soyle i A, og "nuller ut" alle verdier
    #under diagonalen. Starter med forste soyle i matrise A.
    for i in range(n-1):

        #Gaar gjennom alle verdier i soyle i under diagonalen
        #og "nuller dem ut".
        for j in range(n-1,i,-1):
            
            #Samme framgangsmaate som i boka. Starter
            #med aa finne normen til vektoren (A[i,i],A[j,i]).
            r = float(norm(A[i,i],A[j,i]))

            #Dersom r ikke er 0, lager jeg transformasjonsmatrise
            #[ c, s ]
            #[-s, c ]
            #og ganger denne med rad i og j i A, men bare for verdier
            #i disse radene som ligger til hoyre for soyle i.
            if r>0:
                c = A[i,i]/r
                s = A[j,i]/r

                A[[i,j],i+1:n+1]=matrix([[c,s],[-s,c]])*A[[i,j],i+1:n+1]

            #Gir diagonalen riktig verdi, og "nuller ut" verdi A[j,i].
            A[i,i]=r
            A[j,i]=0
    
    
    return sci.solve_triangular(A,b)


def hilbert(n):
    H = zeros((n,n))

    for i in range(n):
        for j in range(n):
            H[i,j]=1./(i+j+1)

    return H

def metric(x,y):
    return sqrt(sum((x-y)**2))
n=20
error = zeros(n)

for i in range(n):
    x = zeros(i+1)+1
    H=hilbert(i+1)
    b=H*x
    error[i]=metric(x,rotsolve(H,b))


plot(error,'b-*')
show()
