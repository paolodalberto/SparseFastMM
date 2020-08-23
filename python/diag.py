import scipy.io as sio
import math
import scipy.spatial.distance as distance
import numpy 
import scipy.stats
import time 
import matplotlib.pyplot as plt
import argparse
import sys 
import pdb

parameters = [
    ("-f", "--file",       str,None,    'store',"Input matrix")
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



def diag(Mat: scipy.sparse.coo.coo_matrix) -> scipy.sparse.coo.coo_matrix :
    diag = {}

    M = Mat.shape[0]
    for c in range(M): diag[c] = []
    Result = Mat.copy()


    for i in range(Mat.nnz):
        key = Mat.row[i] - Mat.col[i]
        if key <0 : key = M +key
        diag[key].append((Mat.row[i],Mat.col[i],Mat.data[i]))

    for c in range(M):
        if len(diag[c]) == 0: del diag[c]
        else:                 diag[c] = sorted(diag[c], key = lambda x: x[0])
    
    running = 0
    for k,vs in sorted(diag.items(), key = lambda x: x[0]):
        for v in vs:
            Result.row [running] = v[0]
            Result.col [running] = v[1]
            Result.data[running] = v[2]
            running +=1
    del diag
    return Result


    
if __name__ == '__main__':


    parser = default_compiler_arg_parser()
    args = parser.parse_args()

    W = sio.mmread(args.file)

    R = diag(W)
    pdb.set_trace()
    sio.mmwrite("temp.mtx",R) 
    
    
