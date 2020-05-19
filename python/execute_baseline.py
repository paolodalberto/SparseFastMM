from multiprocessing import Pool
import scipy
import scipy.io as sio
import math
import numpy
import time


FILES = """
MTX/mult_dcop_03.mtx
MTX/mult_dcop_01.mtx
MTX/mult_dcop_02.mtx
MTX/lp_fit2d.mtx
MTX/bloweya.mtx
MTX/lp_osa_07.mtx
MTX/ex19.mtx
MTX/brainpc2.mtx
MTX/shermanACb.mtx
MTX/cvxqp3.mtx
MTX/case9.mtx
MTX/TSOPF_FS_b9_c6.mtx
MTX/OPF_6000.mtx
MTX/OPF_3754.mtx
MTX/c-47.mtx
MTX/mhd4800a.mtx
MTX/gen4.mtx
MTX/Maragal_6.mtx
MTX/aft01.mtx
MTX/TSOPF_RS_b39_c7.mtx
"""

FILES = """
MTX/mult_dcop_03.mtx
"""

"""
MTX/mult_dcop_01.mtx
MTX/mult_dcop_02.mtx
MTX/lp_fit2d.mtx
MTX/bloweya.mtx
MTX/lp_osa_07.mtx
MTX/ex19.mtx
MTX/brainpc2.mtx
MTX/shermanACb.mtx
MTX/cvxqp3.mtx
MTX/case9.mtx
MTX/TSOPF_FS_b9_c6.mtx
MTX/OPF_6000.mtx
MTX/OPF_3754.mtx
MTX/c-47.mtx
MTX/mhd4800a.mtx
MTX/gen4.mtx
MTX/Maragal_6.mtx
MTX/aft01.mtx
MTX/TSOPF_RS_b39_c7.mtx
"""


T = "/media/ext4/home/paolo/FastMM/Epyc/clSPARSE/bin/Externals/MTX/Bell_Garland/cant/cant.mtx"


import pdb

def timing(A,B, iter=10):
    b = time.time();
    for i in range(0,iter): C=A*B
    e = time.time(); print((e-b)/iter)
    return (e-b)/iter


AC
BCC

def f(x) :
    N=x.shape[1]
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



def explore_execution(AC,BC):
    print("reading")
    B = sio.mmread(T)
    A = sio.mmread(T)

    print("B translating ")
    BC= B.tocsc()
    global BCC= B.tocsr()
    global AC= A.tocsr()

    print("computing")
    if False:
        timing(A,B,3)
        timing(AC,BC,3)
        timing(AC,BCC,3)

    print("Multi")
    K =2
    while  K<16:
        K = K+1
        print(K)
        X = []
        for k in range(0,K):
            X.append([k,K])
        
        p =  Pool(K)
        b = time.time();
        print(X)
        R = p.map(f, X);
        e = time.time(); print(K,e-b)
        for r in R:
            print(K,r[0].nnz,r[1:])
    
        
if __name__ == '__main__':

    for i in FILES.split():
        print(explore_execution(i)) 
    
