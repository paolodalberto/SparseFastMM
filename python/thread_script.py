from multiprocessing import Pool
import scipy
import scipy.io as sio
import math
import numpy
import time
import os
#import pdb
import json
import argparse
import sys

parameters = [
    ("-f", "--file",       str,None,    'store',"Input matrix"),
    ("-d", "--device",     str,'0',    'store',"device"),
    ("-r", "--resfile",    str,"out.json",    'store',"output file"),
    ("-i", "--interval",   int,1000,    'store',"output file"),
    ("-t", "--times",      int,10,      'store',"output file"),
    ("-p", "--poolsize",   int,16,      'store',"output file"),
    
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
MTX/SHUFFLED/0_lp_osa_07.mtx
MTX/SHUFFLED/1_lp_osa_07.mtx
MTX/SHUFFLED/2_lp_osa_07.mtx
MTX/SHUFFLED/3_lp_osa_07.mtx  
MTX/SHUFFLED/4_lp_osa_07.mtx
"""
DEVICE = "CPU"
RESFILE = "result-"+str(DEVICE)+".json"

Template = "Matrix %s K = %2d Time %1.4f time/average %2.4f GFLOPS %2.4f" 

def f(x):
    #print(x)
    AC =x[2] 
    N = AC.shape[0]
    B = x[3]
    IN = x[4]
    
    S = N//x[1]
    l = x[0]*S
    #print(x[0],x[1],x[0]==x[1]-1)
    if x[0]==(x[1]-1):
        r = N
    else:
        r = (x[0]+1)*S

    #print(l,r)
    b = time.time();
    for ll in range(0,IN):
        Rf = AC[l:r,:]*B
    e = time.time(); 
    #print(type(Rf))
    return [(e-b)/IN] #x[0],e-b,r-l,l,r,S,x]
    #return [Rf,e-b] #x[0],e-b,r-l,l,r,S,x]



def compute_parallel(list,
                     DEVICE= 'cpu',
                     R = {'cpu': {}},
                     INTERVAL=5000,
                     TIMES=1,
                     THREADS=16):

    res = R[DEVICE]
    
    for name in list:
        A = sio.mmread(name)
        AC= A.tocsr()
        #pdb.set_trace()
    
        B = numpy.ones(A.shape[1])
        
        
    
        if name not in res:
            res[name] = {}
    
        best = []
        K =1
        while  K<=THREADS:
            
            
            X = []
            for k in range(0,K):
                X.append([k,K,AC,B,INTERVAL])
            with Pool(K) as p:
    
                Rx = []
                b = time.time();
                for k in range(0,TIMES): Rx.extend( p.map(f, X));
                e = time.time(); 
    
                ti = (e-b)/TIMES
                
                
                t = [r[0] for r in Rx] 
                res[name][K] = {
                    'P':K,
                    'nnz':AC.nnz,
                    'shape' : AC.shape,
                    'time':ti,
                    'interval': INTERVAL,
                    'gflops':INTERVAL*2.0*AC.nnz/ti/1000000000,
                    'mean':numpy.mean(t), 'var': numpy.var(t)
                }
                print(Template % (name,  K, res[name][K]['time'],
                                  res[name][K]['time']/res[name][K]['mean']/INTERVAL, res[name][K]['gflops'] ),
                      file=sys.stderr)
                best.append([res[name][K]['gflops'],K,res[name][K]['time']/res[name][K]['mean']/INTERVAL])
                del p
            K = K + (2 if K%2 ==0 else 1)
        #print(name, AC.nnz,sorted(best,key = lambda x: x[0])[-1])
    
        Q = "#P {1:2d} GFLOPS {0:2.3f}  Total/Average {2:2.2f}"
        S = sorted(best,key = lambda x: x[0])[-1]
        
        #print(S)
        result = Q.format(*tuple(S))
        print(result, file= sys.stderr)
        return  (result, S)
        
        
if __name__ == '__main__':


    parser = default_compiler_arg_parser()
    args = parser.parse_args()

    RESFILE = args.resfile
    FILE =  args.file
    DEVICE = args.device
    list = FILES.split()
    if False and os.path.exists(RESFILE):
        print("Exists", RESFILE)
        try :
            with open(RESFILE) as json_file:
                R = json.load(json_file)
                #print(R)
                print("previous results available")
                json_file.close()
                print("read")
        except Exception as e:
            print("Not read",e )
            R= {DEVICE:  {}}
    else:
        R = {DEVICE: {}}

        
    
    #list = os.popen('ls '+ DIRFILE).read().split()
    compute_parallel([FILE], DEVICE,R,args.interval, args.times,args.poolsize)    
