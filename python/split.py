

import scipy 
import scipy.io as sio
import ent
import numpy
import matplotlib.pyplot as plt


def visualize(x):
    plt.clf()
    ax1 = plt.subplot(111)
    plt.plot(x)
    plt.show()
def dense(A):
    Hd,Wd = A.shape
    nnz = A.nnz

    DH = numpy.zeros(Hd,dtype='long')
    DW = numpy.zeros(Wd,dtype='long')

    for i in range(0,nnz):
        c = A.col[i]; r = A.row[i]
        DW[c-1] +=1 ; DH[r-1] +=1  

    H = []
    for i in range(Hd):
        H.append([i,DH[i]])
    H = sorted(H, key = lambda x: x[1])
    #import pdb; pdb.set_trace()
    #visualize(sorted(DH.flatten()))


    count_h = len(DH[DH>Wd/10])
    DenseH = numpy.zeros(Wd*count_h).reshape(count_h,Wd)
    
    W = []
    for i in range(Wd):
        W.append([i,DW[i]])
    W = sorted(W, key = lambda x: x[1])
    #import pdb; pdb.set_trace()
    #visualize(sorted(DW))
    
    count_w = len(DW[DW>Hd/10])
    DenseW = numpy.zeros(Hd*count_w).reshape(Hd,count_w)
    
    ## what considered dense ?
    Sparse = A*0.0
    
    

    for i in range(0,nnz):
        c = A.col[i]
        r = A.row[i]
        val = A.data[i]
        c = W[c][0]
        r = H[r][0]
        
        if c<Wd-count_w and  r<Hd-count_h:
            Sparse.col[i] = c
            Sparse.row[i] = r
            Sparse.data[i] = val
        elif c>Wd-count_w and  r<=Hd-count_h:
            DenseW[r,Wd - c-1] = val
        elif r> Hd-count_h:
            DenseH[Hd-r-1,c] = val
    #result matrix plus swaps
    
    

    Sparse.eliminate_zeros()
    return Sparse,DenseW,DenseH,W,H

def split(args):
    A = sio.mmread(args.file)
    TA = ent.measure_a(A,1000,args)
    #print(TA)
    B,DW,DH,W,H = dense(A)
    
    T = ent.measure_a(B,1000,args)
    #print(T)
    if DW.shape[0]* DW.shape[1] >0:
        TD1 = ent.measure_dense_a(DW,1000,args)
        T[1][0] += TD1[1][0]
        T[1][3] += TD1[1][0]
        #print(TD1)
    if DH.shape[0]* DH.shape[1] >0:
        TD2 = ent.measure_dense_a(DH,1000,args)
        #print(TD2)
        T[1][0] += TD2[1][0]
        T[1][3] += TD2[1][0]
        
    #print("Whole  Time COO", TA[1][0],"Time CSR", TA[1][3])
    #print("Sparse Time COO", T[1][0], "Time CSR", T[1][3])

    return 100*((TA[1][0] - T[1][0])/TA[1][0]),100*((TA[1][3] - T[1][3])/TA[1][3])
        
if __name__ == '__main__':


    parser = ent.default_compiler_arg_parser()
    args = parser.parse_args()

    for t in range(args.trials):
        print(split(args))

    
