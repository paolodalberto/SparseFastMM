from multiprocessing import Pool
import scipy
import scipy.io as sio
import math
import numpy
import time

B = sio.mmread("bcsstk32.mtx")
AS = sio.mmread("shufflebcsstk32.mtx")
A = sio.mmread("bcsstk32.mtx")

BC= B.tocsc()
BCC= B.tocsr()
AC= A.tocsr()

b = time.time(); C=A*B; e = time.time(); print(e-b)
b = time.time(); CC=AC*BC; e = time.time(); print(e-b)
b = time.time(); CC=AC*BCC; e = time.time(); print(e-b)


def f(x,N=AC.shape[0]):
    S = N//x[1]
    l = x[0]*S
    #print(x[0],x[1],x[0]==x[1]-1)
    if x[0]==(x[1]-1):
        r = N
    else:
        r = (x[0]+1)*S

    #print(x,l,r)
    b = time.time();
    R = AC[l:r,:]*BCC
    e = time.time(); 
    #print(type(R))
    return [R,e-b] #x[0],e-b,r-l,l,r,S,x]



K =2
while  K<16:
    K = K+1
    print(K)
    X = []
    for k in range(0,K):
        X.append([k,K])
    
    p =  Pool(K)
    b = time.time();
    R = p.map(f, X);
    e = time.time(); print(K,e-b)
    for r in R:
        print(K,r[0].nnz,r[1:])

    

    
