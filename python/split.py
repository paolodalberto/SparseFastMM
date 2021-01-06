
from multiprocessing import Pool
import scipy
import scipy.io as sio
import ent
import numpy
import matplotlib.pyplot as plt
import time

def visualize(x):
    plt.clf()
    ax1 = plt.subplot(111)
    plt.plot(x)
    plt.show()

######################
## inputs : A COO matrix
## outputs: 
##       S COO matrix
##       W dense
##       H dense
##       sort vectors for permute
######################
def dense(A, x = 2 ):
    Hd,Wd = A.shape #height, width
    nnz   = A.nnz    

    ## histograms by height and width
    DH = numpy.zeros(Hd,dtype='int')
    DW = numpy.zeros(Wd,dtype='int')

    for i in range(0,nnz):
        c = A.col[i]; r = A.row[i]
        DW[c-1] +=1 ; DH[r-1] +=1  

    ## dense part by height
    count_h = len(DH[DH>Wd/10]) if x&1>0 else 0 
    DenseH = numpy.zeros(Wd*count_h,dtype =float).reshape(count_h,Wd)

    count_w = len(DW[DW>Hd/10])
    count_w = count_w if  count_w>= 4 and  x&2>0 else 0 
    DenseW = numpy.zeros(Hd*count_w,dtype=float).reshape(Hd,count_w)

    if count_h ==0 and  count_w == 0:
        return A,DenseW,DenseH,None,None
    H = []
    for i in range(Hd):
        H.append([i,DH[i]])
    H = sorted(H, key = lambda x: x[1])

    
    W = []
    for i in range(Wd):
        W.append([i,DW[i]])
    W = sorted(W, key = lambda x: x[1])

    ## dense part by width

    
    ## what considered dense ?
    Sparse = A*0.0
    
    for i in range(0,nnz):
        c = A.col[i]
        r = A.row[i]
        val = A.data[i]
        c = W[c][0]
        r = H[r][0]
        
        if c<Wd-count_w and  r<Hd-count_h:
            # Sparse
            Sparse.col[i] = c
            Sparse.row[i] = r
            Sparse.data[i] = val
        elif c>Wd-count_w and  r<=Hd-count_h:
            # Dense width 
            DenseW[r,Wd - c-1] = val
        elif r> Hd-count_h:
            # dense height
            DenseH[Hd-r-1,c] = val
    #result matrix plus swaps
    Sparse.eliminate_zeros()
    return Sparse,DenseW,DenseH,W,H

def split_compute(S,W,H,B):
    
    h,w = S.shape 
    Rfw = None
    Rfh = None

    Rfs = S*B
    if W.shape[1]>0:
        Rfw = numpy.dot(W,B[w-W.shape[1]:w])
        Rfs += Rfw
    if H.shape[0]>0: 
        Rfh = numpy.dot(H,B)
        Rfs[h-H.shape[0]:h] += Rfh
        
    
    return Rfs




def measure_local(A, S,DW,DH,W,H, IN=10000, args = None) :
    if args is not None  and args.gpuonly: return "", [0,0]
    
    B = numpy.ones(A.shape[1])
    #S,DW,DH,W,H = dense(A) 
    #print(S.shape, DW.shape, DH.shape)
    b = time.time();
    for ll in range(0,IN):
        Rf = split_compute(S,DW,DH,B)
    e = time.time();
    coo = (e-b)/IN
    #print("COO", IN,e-b, file=sys.stderr)
    nnz = A.nnz
    coogflops = (2*nnz/coo)/1000000000
    AC= S.tocsr()
    b = time.time();
    for ll in range(0,IN):
        Rf = split_compute(AC,DW,DH,B)
    e = time.time(); 
    csr = (e-b)/IN
    #print("CSR", IN,e-b, file=sys.stderr, flush= True)
    csrgflops = (2*nnz/csr)/1000000000
    
    return ("# COO GFLOPS {0:5.2f} CSR GFLOPS {1:5.2f} ".format(coogflops,csrgflops),[coo,coogflops,csr,csrgflops])

def fun(x):
    S,W,H,B = x[1]
    h,w = S.shape
    #print(x[0],S.shape,W.shape, H.shape, B.shape) 
    if x[0]==0:
        return S*B
    elif x[0] ==1 :
        if W.shape[1]>0: return  numpy.dot(W,B[w-W.shape[1]:w])
        else: return None
    elif x[0] ==2 :
        if H.shape[0]>0: return numpy.dot(H,B)
        else: return None
    return None
def measure_local_pool(A, IN=10000, args = None) :
    if args is not None  and args.gpuonly: return "", [0,0]
    P = Pool(16)

    
    B = numpy.ones(A.shape[1])
    S,DW,DH,W,H = dense(A) 
    T = [S,DW,DH,B]
    X = [ [0, T],
          [1, T],
          [2, T]]
    
    
    def split_compute_pool(X,P):
        
        Rx = []
        Rx.extend( P.map(fun, X));
        rs = Rx[0]
        for r in Rx[1:]:
            if r is not None:
                rs += r
        return rs
    

    #print(S.shape, DW.shape, DH.shape)
    b = time.time();
    for ll in range(0,IN):
        Rf = split_compute_pool(X,P)
    e = time.time();
    coo = (e-b)/IN
    #print("COO", IN,e-b, file=sys.stderr)
    nnz = A.nnz
    coogflops = (2*nnz/coo)/1000000000
    AC= S.tocsr()
    b = time.time();
    for ll in range(0,IN):
        Rf = split_compute_pool(X,P)
    e = time.time(); 
    csr = (e-b)/IN
    #print("CSR", IN,e-b, file=sys.stderr, flush= True)
    csrgflops = (2*nnz/csr)/1000000000
    
    return ("COO GFLOPS {0:5.2f} CSR GFLOPS {1:5.2f} ".format(coogflops,csrgflops),[coo,coogflops,csr,csrgflops])

def split_2(A, S,DW,DH,W,H):
    ## read matrix
    TA1 = ent.measure_a(A,100,args)
    
    #print(TA1)
    TA2 = measure_local(A, S,DW,DH,W,H,100,args)
    #print(TA2)
    return TA1, TA2
    
def split(args):
    ## read matrix
    A = sio.mmread(args.file)
    ## estimate of sparse computation time
    TA = ent.measure_a(A,100,args)
    print(TA)
    # creating sparse, dense, dense
    B,DW,DH,W,H = dense(A)

    ## estimate of new sparse 
    T = ent.measure_a(B,100,args)
    #print(T)
    if DW.shape[0]* DW.shape[1] >0:
        TD1 = ent.measure_dense_a(DW,1000,args)
        T[1][0] += TD1[1][0]
        T[1][3] += TD1[1][0]
    if DH.shape[0]* DH.shape[1] >0:
        TD2 = ent.measure_dense_a(DH,1000,args)
        T[1][0] += TD2[1][0]
        T[1][3] += TD2[1][0]
    ## relative measyre
    return 100*((TA[1][0] - T[1][0])/TA[1][0]),100*((TA[1][3] - T[1][3])/TA[1][3])


def fun_2(filename):
    A = sio.mmread(filename)
    S,DW,DH,W,H = dense(A)
    #print(filename,A.nnz,S.shape,DW.shape,DH.shape)
    if DH.shape[0]* DH.shape[1] ==0 and DW.shape[0]* DW.shape[1] ==0:
        return None

    #print("Relative")
    #for t in range(args.trials):
    #    print(split(args))

    COO_1 = []
    CSR_1 = []
    COO_2 = []
    CSR_2 = []
    #print("Simpler")
    #bar = Bar('Processing', max=args.trials)
    for t in range(args.trials):
        #bar.next()
        a,b = split_2(A, S,DW,DH,W,H)
        COO_1.append(a[1][1])
        CSR_1.append(a[1][3])
        COO_2.append(b[1][1])
        CSR_2.append(b[1][3])
    #bar.finish()
    #print(filename,
    #      "COO", max(COO_1),"COO Hybrid",max(COO_2),
    #      "CSR", max(CSR_1),"CSR Hybrid",max(CSR_2))
        
    plt.clf()
    ax1 = plt.subplot(111)
    plt.hist(COO_1,bins=100, alpha=0.5, label="coo-sparse")
    plt.hist(CSR_1,bins=100, alpha=0.5, label="CSR-sparse")
    plt.hist(COO_2,bins=100, alpha=0.5, label="coo-mix")
    plt.hist(CSR_2,bins=100, alpha=0.5, label="CSR-mix")

    plt.xlabel("GFLOPS", size=14)
    plt.ylabel("Count", size=14)
    plt.title("Sparse vs Sparse + Dense ")
    plt.legend(loc='upper right')
    name = args.file[args.file.rfind("/")+1:]
    plt.savefig(name+".png")
    return [filename,
            A.nnz,S.shape,DW.shape,DH.shape,
            "COO", max(COO_1),"COO Hybrid",max(COO_2),
            "CSR", max(CSR_1),"CSR Hybrid",max(CSR_2)]        

if __name__ == '__main__':

    from progress.bar import Bar
    parser = ent.default_compiler_arg_parser()
    args = parser.parse_args()

    files =  args.file.split(",")
    #print(files)
    P = Pool(8)
    Rx = []
    Rx.extend( P.map(fun_2, files));

    Rs = []
    for r in Rx:
        if r is not None:
            print(r)
            Rs.append(r[0])
    print(Rs)
