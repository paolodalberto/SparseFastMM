import os
import json
import collections
import argparse
import sys


H = "Device file precision layout M           N           nnz         alpha       beta        Algorithm   GFlops     GBs        msec        iter        verified".split()
H2 = "Device file precision layout m n nnz alpha beta GFlops GBs msec".split()
#print("H", H,file=sys.stderr)
Result = collections.namedtuple(
    "Result",
    H
)

parameters = [
    ("-f", "--file",       str,None,    'store',"Input matrix"),
    ("-d", "--device",     str,'0',    'store',"device"),
    ("-r", "--resfile",    str,"out.json",    'store',"output file"),
    
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


NOTFOUND = """
MTX/loweya.mtx
"""
CANNOTREAD = """
MTX/brainpc2.mtx
MTX/case9.mtx
MTX/TSOPF_FS_b9_c6.mtx
"""



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
MTX/lp_osa_07.mtx
MTX/SHUFFLED/1_lp_osa_07.mtx
MTX/SHUFFLED/2_lp_osa_07.mtx
MTX/SHUFFLED/3_lp_osa_07.mtx
MTX/SHUFFLED/4_lp_osa_07.mtx
"""



DIRFILE = FILES # "MTX/*.mtx"
DEVICE  = "1"
RESFILE = "result-"+str(DEVICE)+".json"
STUBcsr ="""/home/paolo/FastMM/Epyc/rocSPARSE/build/release/clients/staging/rocsparse-bench -d %s -r %s -f csrmv --mtx %s -i 10000  -v 0 | grep -P "(.+\s)+(.+\s)+(.+\s)+(.+\s)+(.+\s)+(.+\s)+(.+\s)+(.+\s)+" | grep -v -i "device" | grep -v with """
STUBcoo ="""/home/paolo/FastMM/Epyc/rocSPARSE/build/release/clients/staging/rocsparse-bench -d %s -r %s -f coomv --mtx %s -i 10000   -v 0 | grep -P "(.+\s)+(.+\s)+(.+\s)+(.+\s)+(.+\s)+(.+\s)+(.+\s)+(.+\s)+" | grep -v -i "device" | grep -v with """


#print("result file", RESFILE)

def gpu_execute(DEVICE : int = 0,
                FILES  : str = FILES,
                R : dict = {0: {'s' : {}, 'd' : {}}}) :

    L = FILES.split()
    res = R[DEVICE]
    #print(L)
             
    for precision in ['s','d']:
        pres = res[precision]
        
        for f in L:
            
            command = STUBcsr %(DEVICE, precision, f)
            try:
                myCmd = os.popen(command).read().split("\n")
                rhead  = ['device','matrix', 'precisio','layout']
                r = [DEVICE,f,precision,'csr']
                
                rhead.extend(myCmd[0].replace("/","").split())
                r.extend(myCmd[1].split())
    
                Result = collections.namedtuple(
                    "Result",
                    rhead
                )
                
                pres['csr'] = Result(*(r))
                print(pres['csr'],file=sys.stderr)
                result = {}
            
                for i in range(0,len(rhead)):
                    try:    result[rhead[i]] = float(r[i])
                    except: result[rhead[i]] = r[i]
                pres['csr'] = result
            except Exception as e:
                print(e, myCmd,file=sys.stderr)
                pres[f+'csr'] = None
    
            try:

                command = STUBcoo %(DEVICE, precision, f)
                myCmd = os.popen(command).read().split("\n")
                rhead  = ['device','matrix', 'precision','layout']
                r = [DEVICE,f,precision,'coo']
                
                rhead.extend(myCmd[0].replace("/","").split())
                r.extend(myCmd[1].split())
                Result = collections.namedtuple(
                    "Result",
                    rhead
                )
                pres['coo'] = Result(*(r))
                print(pres['coo'],file=sys.stderr)
                result = {}
            
                for i in range(0,len(rhead)):
                    try:    result[rhead[i]] = float(r[i])
                    except: result[rhead[i]] = r[i]
                pres['coo'] = result
            except Exception as e:
                print(e, myCmd,file=sys.stderr)
                print(precision,f,"Skip", e,file=sys.stderr)
                pres['coo'] = None
            #print(pres[f])        
    
    return R

if __name__ == '__main__':


    parser = default_compiler_arg_parser()
    args = parser.parse_args()

    RESFILE = args.resfile
    FILES =  args.file
    DEVICE = args.device
    
    #list = os.popen('ls '+ DIRFILE).read().split()
    list = FILES.split()
    if os.path.exists(RESFILE):
        #print("Exists", RESFILE)
        try :
            with open(RESFILE,'r') as json_file:
                #import pdb; pdb.set_trace()
                R = json.load(json_file)
                #print(R)
                #print("previous results available",
                #      len(R[DEVICE]['s']),len(R[DEVICE]['d']))
                json_file.close()
                #print("read")
        except Exception as e:
            print("Not read",e )
            R= {DEVICE: {'s' : {}, 'd' : {}}}
    else:
        R = {DEVICE: {'s' : {}, 'd' : {}}}
    
    count =0
    gpu_execute(DEVICE,FILES,R)
    with open(RESFILE,'w') as outfile:
        json.dump(R, outfile)
        outfile.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
