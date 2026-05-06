'''
Get the simulation data, parse it, and put it in new files
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import os,sys
from types import SimpleNamespace
import re

def isnum(s):
    try:
        a=float(s)
    except ValueError:
        return False
    return True

# function to read in simulation info
def get_params(paramfile):
    f=open(paramfile,'r')
    lines = f.readlines()
    f.close()

    params = [re.split(' is | = ',x.replace('\n','')) for x in lines if ' is ' in x or ' = ' in x]
    params = {x[0]: x[1].split(' ')[0] for x in params}
    to_del = []
    to_add = []
    for k in params.keys():
        if isnum(params[k]):
            params[k] = float(params[k])
        if k.startswith('N'):
            params[k] = int(params[k])
        if type(params[k])==str and ',' in params[k]:
            params[k] = np.array([float(x) for x in params[k].split(',')])
    params = SimpleNamespace(**params)
    return params



##########################################################################################################################################
# READ DATA

folders=['int','nonint']

# collect data
for f in folders:

    # list of different seeds for each system (interacting & noninteracting)
    seeds = [os.path.join(f,x.replace('-param','')) for x in os.listdir(f) if x.endswith('-param')]
    
    # get param info
    p = get_params(seeds[0]+'-param')
    dx = p.Lx/p.Nbinx
    dy = p.Ly/p.Nbiny
    iLx = int(p.Lx+1e-5)

    # displacement data    
    disps_ = [np.mean(np.genfromtxt(x+'-disp',delimiter='\t')[::2,1:],axis=1) for x in seeds]
    disps_ = np.array([np.concatenate([[0],x]) for x in disps_])
    disp_avg = np.mean(disps_,axis=0)
    disp_err = np.std(disps_,axis=0,ddof=1)/np.sqrt(len(seeds))

    # density profile
    profs_ = [np.genfromtxt(x+'-prof',delimiter='\t')[:,1:]/p.NstepProf/p.N/dx for x in seeds]
    prof = sum(profs_)/len(seeds)

    # stress tensor data
    if 'nonint' not in f:

        # thermal/active stress
        sigAs_ = [np.genfromtxt(x+'-sigmaAprof',delimiter='\t')[:,1]/p.NstepProf/p.N/dx for x in seeds]
        sigAs_ = [np.array([x[p.Nbinx*i:p.Nbinx*(i+1)] for i in range(8)]) for x in sigAs_]
        sigA = np.transpose(sum(sigAs_) / len(seeds))
        sigA[:,4:] = sigA[:,4:] / p.dt

        # interaction stress
        sigIKs_ = [np.sum(np.genfromtxt(x+'-sigmaIKprof',delimiter='\t')[:,1:],axis=1)/p.NstepProf/p.dt/p.N for x in seeds]
        sigIKs_ = [np.array([x[iLx*i:iLx*(i+1)] for i in range(3)]) for x in sigIKs_]
        sigIK = np.transpose(sum(sigIKs_) / len(seeds))

        # save stress data
        np.savetxt(os.path.join(f,'sigA_avg'),sigA,delimiter='\t',fmt='%.4e')
        np.savetxt(os.path.join(f,'sigIK_avg'),sigIK,delimiter='\t',fmt='%.4e')


    # save data
    np.savetxt(os.path.join(f,'disp_avg'),disp_avg,delimiter='\t',fmt='%.4e')
    np.savetxt(os.path.join(f,'disp_err'),disp_err,delimiter='\t',fmt='%.4e')
    np.savetxt(os.path.join(f,'prof_avg'),prof,delimiter='\t',fmt='%.4e')


