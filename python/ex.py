import sparsecoo
import pdb
import scipy
import scipy.io as sio
import math


SHUFFLE=False

if SHUFFLE:
    B = sio.mmread("bcsstk32.mtx")
    A = sio.mmread("shufflebcsstk32.mtx")
    C = sio.mmread("shufflebcsstk32.mtx")
else:
    B = sio.mmread("bcsstk32.mtx")
    A = sio.mmread("bcsstk32.mtx")
    C = sio.mmread("bcsstk32.mtx")

    
def orderbyrow(A):
    AL = []
    print("Ordering A Row Major")
    for i in range(0,A.nnz):
        AL.append([A.row[i],A.col[i], i])
    
    AL = sorted(AL, key = lambda x: (x[0],x[1]))
    I = []
    for i in range(0,A.nnz):
        I.append(AL[i][2])
    A.row = A.row[I]
    A.col = A.col[I]
    A.data = A.data[I]
    return A

print("Ordering C by row")
A = orderbyrow(A)

print("Ordering C as A")
C = orderbyrow(C)

print("Ordering B (this will be reorg in the kenel")

B = orderbyrow(B)


DEBUG = False


def sparse_mm(P,C=C,A=A,B=B):


    

    #B.row = B.row[0:IP]
    #B.col = B.col[0:IP]
    #B.data = B.data[0:IP]

    #pdb.set_trace()

    if DEBUG:
        print("Calling ")
        print(type(P),type(A.col),A.col.dtype,type(A.row),A.row.dtype,type(A.data),A.data.dtype,type(A.shape[0]),type(A.shape[1]))
        print(P,C.col.shape,C.row.shape,C.data.shape,C.shape)
        print(C.col[3],C.row[3],C.data[3])
        print(P,A.col.shape,A.row.shape,A.data.shape,A.shape)
        print(type(B.col),type(B.row),type(B.data),type(B.shape[0]),type(B.shape[1]))
        print(B.col.shape,B.row.shape,B.data.shape,B.shape)
        
        print("testing A",
              (A.col>A.shape[1]).any(),(A.col<0).any(),
              (A.row>A.shape[0]).any(),(A.row<0).any()
        )
        
        print("testing B",
              (B.col>B.shape[0]).any(),(B.col<0).any(),
              (B.row>B.shape[1]).any(),(B.row<0).any() 
        )
        
        print("testing C",
              (C.col>C.shape[1]).any(),(C.col<0).any(),
              (C.row>C.shape[0]).any(),(C.row<0).any() 
        )
    
    cr,cc,cv,time,ops =  sparsecoo.pcoomul(
        P,
        C.col,C.row,C.data,C.shape[0],C.shape[1],
        A.col,A.row,A.data,A.shape[0],A.shape[1],
        B.col,B.row,B.data,B.shape[0],B.shape[1]    # transposeing B
    )
    #pdb.set_trace()
    
    #
    if DEBUG: print(A.shape,cr.shape,cc.shape,cv.shape, time)
    
    #pdb.set_trace()
    return scipy.sparse.coo_matrix((cv, (cr, cc))),time,ops

if  True:
    R = []
    for i in range(2,17):
        G,time,ops =sparse_mm(i)
        print(i,G.nnz,time, G.nnz/time[0], ops[0], ops[0]/time[0])
        R.append([i,G.nnz,time, G.nnz/time[0], ops[0], ops[0]/time[0]])

    for r in R:
        print(r)
