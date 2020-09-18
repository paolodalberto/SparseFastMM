import scipy.io as sio
import math
import scipy.spatial.distance as distance
import numpy 
import scipy.stats
import time 
import matplotlib.pyplot as plt
import argparse
import sys 
from collections import namedtuple
from multiprocessing import Pool
from functools import reduce    
import random

from  gpu_script import gpu_execute
from  thread_script import compute_parallel,compute_parallel_2

class Element:
    def __init__(self,row,col,data):
        self.row = row
        self.col = col
        self.data=data
        self.strange = 0.0
        self.prob = 0.0
    def __str__(self):
        return str((self.row, self.col,self.strange, self.prob))

Types = [
    'row', 'col', 'diag', 'rand', #'xXx'
]
TU    = [ 'ro',  'go',  'yo',   'bo',   'r^']
        

def GetNNZ(M, order='row'):
    nnzs = [Element(M.row[i], M.col[i], M.data[i]) for i in range(M.nnz)]
    cnum = M.shape[0]
    if order == 'row':
        return sorted(nnzs, key = lambda tup: (tup.row,tup.col) )
    elif order == 'col':
        return sorted(nnzs, key = lambda tup: (tup.col, tup.row))
        
    elif order == 'diag':
        return sorted(nnzs, key = lambda tup: ((cnum+tup.row-tup.col)%cnum,tup.row))
    elif order == 'rand':
        random.seed(3)
        random.shuffle(nnzs)
        return nnzs
    
    return nnzs


def norm_diff(x : Element,y: Element,
              Mod : int = 8):
    s =0.0
    a = (x.row - y.row)%Mod
    b = (x.col - y.col)%Mod
    
    return a*a+b*b


def distance(w : Element,
             I : list ,
             indexes = list,
             Mod : int = 8):

    x = norm_diff(w,I[indexes[0]],Mod)
    count =1 
    for j in indexes[1:]:
        i = I[j]
        y = norm_diff(w,I[j])

        if y==x:  count +=1
        elif y<x: x= y; count=1
        
    return x/count

def Strange_Set(I : list,
                Mod: int = 8):
    W = [i for i in range(len(I))]
    T = [i for i in range(len(I))]

    hist = {}
    
    for i in W:
        w = I[i]
        T.pop(i)
        D = distance(w,I,T,Mod)
        w.strange = D
        if not D in hist: hist[D] = 0
        hist[D] +=1
        T.insert(i,i)

    N = len(I)
    for i in I:
        i.prob = hist[i.strange]/N
        #if i.prob>1: import pdb;pdb.set_trace()
    

    
    
def H(I : list):
    h = 0
    for i in I:
        if i.prob >0: h += -i.prob*math.log2(i.prob)

        
    return h
    
parameters = [
    ("-f", "--file",       str,None,    'store',"Input matrix"),
    ("-w", "--window",     int,8,    'store',"window 8,16,32"),
    ("-m", "--modulo",     int,8,       'store',"mesh size 8,16,32"),
    ("-png", "--png",      str,"result.png",       'store',"picturemesh size 8,16,32"),
    ("-full", "--full",    bool,False,       'store',"all windows"),
]
    
def default_compiler_arg_parser(params=parameters):
    # NOTE: Not all arguments are passed as separate variables to backend functions.
    # Removing a line here may cause errors until we can completely remove the args parameter from backend functions

    parser = argparse.ArgumentParser()
    for x in params:
        #print("Adding arguments:",x)
        if x[2] is bool:
            if x[0] is not None:
                parser.add_argument(x[0], x[1], default=x[3], action=x[4], help=x[5])
            else:
                parser.add_argument(x[1], default=x[3], action=x[4], help=x[5])
        elif x[0] is not None:
            parser.add_argument(x[0], x[1], type=x[2], default=x[3], action=x[4], help=x[5])
        else:
            parser.add_argument(x[1], type=x[2], default=x[3], action=x[4], help=x[5])
    return parser



def Entropy_Time_Series(t : list ,
                        W : int = 8,
                        M : int = 8,
                        identifier="#"):

    #print(identifier,len(t))
    hs = [] 
    k = 0
    while k<len(t)-W:
        T = t[k:k+W]
        Strange_Set(T,M)
        hoft = H(T)
        #        if hoft < 0:            import pdb;pdb.set_trace()
        #        if k> K and (k//args.window)%K==0:
        #            print(identifier,k,hoft)
        #            for f in T: print(identifier, f)
        #            #pdb.set_trace()
        hs.append(hoft)
        k+= W

    return  hs


if __name__ == '__main__':
    from multiprocessing import Pool
    from multiprocessing import Process 


    parser = default_compiler_arg_parser()
    args = parser.parse_args()



    def f(X):
        #print(t)
        w,m,t = X
        T = GetNNZ(W,t)
        hs = Entropy_Time_Series(T,w,m,t)
        return hs
    
    W = sio.mmread(args.file)
    #print(W.nnz)

    if args.full:
        TYPES = [[8,   8, i] for i in Types] +\
                [[16, 16, i] for i in Types] +\
                [[32, 32, i] for i in Types] 
        
        
    else :
        TYPES = [ [args.window, args.modulo, i] for i in Types] 
    Hs = []
    with Pool(len(TYPES)) as p:
        
        Hs.append(p.map(f,TYPES))

    #import pdb; pdb.set_trace()
    for i in range(len(Hs[0])):
        print(args.file, TYPES[i][0],TYPES[i][1],TYPES[i][2], numpy.mean(Hs[0][i]), numpy.var(Hs[0][i]))
        
    for i in range(1 if not args.full else 3):
        #plt.clf()
        #import pdb; pdb.set_trace()
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
        plt.violinplot(Hs[0][i*len(Types):(i+1)*len(Types)], showmeans=False,showmedians=True)
        
        axes.yaxis.grid(True)
        axes.set_xticks([y + 1 for y in range(len(Types))])
        axes.set_xlabel('Rand')
        axes.set_ylabel('Streaming Entropy')
        plt.setp(axes,xticks=[y + 1 for y in range(len(Types))],xticklabels=Types)
        N = str(TYPES[i*len(Types)]).replace("\s", "").replace(",","x").replace("[","").replace("]","")
        #plt.show() #
        plt.savefig(N+args.png)
                           

