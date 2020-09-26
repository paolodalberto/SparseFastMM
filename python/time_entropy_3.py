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
        
    def __str__(self):
        return str((self.row, self.col))

class Mesh:
    def __init__(self,M : int ):
        self.M = M
        self.x = numpy.zeros(M,dtype=int)
        self.y = numpy.zeros(M,dtype=int)
        
    def add(self,i : Element):
        
        a = (i.row)%self.M
        b = (i.col)%self.M
        
        self.x[a] += 1
        self.y[b] += 1
    def sub(self,i : Element):
        a = (i.row)%self.M
        b = (i.col)%self.M
        self.x[a] -= 1
        self.y[b] -= 1
        

    def x(self):
        return numpy.count_nonzero(self.x)
    def y(self):
        return  numpy.count_nonzero(self.y)
    def entropy(self):
        x = numpy.count_nonzero(self.x)
        y = numpy.count_nonzero(self.y)
        
        return (x+ y)/(2*self.M)


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

        

def Entropy_Time_Series(t : list ,
                        W : int = 8,
                        M : int = 8,
                        identifier="#"):
    mesh = Mesh(M)
    res = []
    
    for k in range(W):
        mesh.add(t[k])
    for k in range(W,len(t)):
        mesh.add(t[k])
        res.append(mesh.entropy())
        mesh.sub(t[k-W])
        
    return  res


parameters = [
    ("-f", "--file",       str,None,    'store',"Input matrix"),
    ("-w", "--window",     int,8,    'store',"window 8,16,32"),
    ("-p", "--pool",       int,7,    'store',"window 8,16,32"),
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



if __name__ == '__main__':
    from multiprocessing import Pool
    from multiprocessing import Process 


    parser = default_compiler_arg_parser()
    args = parser.parse_args()

    MX = {}

    #print("reading")
    W = sio.mmread(args.file)
    
    def g(t):
        #print("sorting",t)
        M = GetNNZ(W,t)
        #print("sorting",t, "done")
        return [t,M] #print(t,MX[t][0])
    Ks = []
    with Pool(len(Types)) as p:
        Ks.append(p.map(g,Types))
    
    for t in Ks[0]:
        
        MX[t[0]] = t[1]
        #print(t[0] ,MX[t[0]][0])
    
    #print("computing")
    del W
    
    def f(X):
        #print(t)
        w,m,t = X
        hs = Entropy_Time_Series(MX[t],w,m,t)
        return hs


    if args.full:
        TYPES = [[8,   8, i] for i in Types] +\
                [[16, 16, i] for i in Types] +\
                [[32, 32, i] for i in Types] 
        
        
    else :
        TYPES = [ [args.window, args.modulo, i] for i in Types] 
    Hs = []
    with Pool(args.pool) as p:
        
        Hs.append(p.map(f,TYPES))

    #import pdb; pdb.set_trace()
    for i in range(len(Hs[0])):
        print(args.file, TYPES[i][0],TYPES[i][1],TYPES[i][2], numpy.mean(Hs[0][i]), numpy.var(Hs[0][i]))
        
    for i in range(1 if not args.full else 3):
        #plt.clf()
        #import pdb; pdb.set_trace()
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
        plt.violinplot(Hs[0][i*len(Types):(i+1)*len(Types)], showmeans=True,showmedians=False)
        
        axes.yaxis.grid(True)
        axes.set_xticks([y + 1 for y in range(len(Types))])
        axes.set_xlabel('Rand')
        axes.set_ylabel('Streaming Entropy')
        plt.setp(axes,xticks=[y + 1 for y in range(len(Types))],xticklabels=Types)
        N = str(TYPES[i*len(Types)]).replace("\s+", "").replace(",","x").replace("[","").replace("]","")
        #plt.show() #
        plt.savefig(N+args.png)




