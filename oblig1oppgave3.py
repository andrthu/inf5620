from numpy import *

#dimensjon til matrise
n=8

#lager en matrise med tilfeldige verdier
A=matrix(random.random((n,n)))

#gjor matrisen ovre trianguler ved aa sette verdier 
#under diagonal til 0
A=triu(A)

#lager kopi av A.
U=A.copy()



#Denne nestede lokka inverterer A. Dette gjores ved aa lose Ab(j)=e(j) 
#j=n,...1. Her er b(j) soyle nummer j i den inverse til A(alsaa en ukjent), 
#og e(j) er enhetsvektoren med veri 1 i plass j i vektoren og 0 ellers. 
#Algoritmen starter med aa gjore dette for soyle n i A, og finner den 
#"nederste" verdien i b forst. Naar en verdi i b er regnet ut lagres den
#paa plassen i A, som den utregnede verdien ville hatt i den inverse til A.
#Algoritmen er laget slik at man ikke trenger de verdiene i A som forsvinner
#naar man bytter dem ut med verdiene i den inverse. Dette er mulig siden
#vi vet at den inverse i til ovre triangler matrise er en ovre trianguler
#matrise. 

#k loper over soylene i A. starter med soylen helt til hoyre 
for k in range(n-1,-1,-1):
    #ganger rad k i a med soyle k i A invers. Dette skal bli lik 1.
    #den eneste verdien i produktet som ikke nulles ut er verdiene paa
    #diagonalen.
    U[k,k]=1./U[k,k]
    
    #r looper over verdiene i soyle k. Starter nederst.
    for r in range(k-1,-1,-1):
        #loser ligningen: 
        #A[r,r]*b(k)(r) + A[r,r+1]b(k)(r+1)+...+A[r,k]b(k)(k)=0
        #for b(k)(r) og lagrer verdien i A[r,k] for i>r er
        #b(k)(i) lagret i A[i,k], og vi faar dermed uttrykket under. 
        U[r,k]= -U[r,(r+1):(k+1)]*U[(r+1):(k+1),k]/U[r,r]


#ganger U med A, og ser fra utskrift at matrise produktet
#av A og U blir identiten. Dette betyr at U er den inverse
#til A. 
print U*A
