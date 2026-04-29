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

# list of folders: 4 different types of system, interacting + non-interacting
folders = ['UPBPs','PBPs','ABPs','AOUPs']
folders = list(np.concatenate([[f+'/int',f+'/nonint'] for f in folders]))
subdir1 = 'tf_8000000-dt_0.000250-N_120-Lx_30-Ly_10-Nmode_7'
subdir2 = 'tf_8000000-dt_0.000250-N_120-Lx_30-Ly_10-Nmode_7-k_10'
fs = [os.path.join(f,subdir1) if f.endswith('nonint') else os.path.join(f,subdir2) for f in folders]

# collect data
for f in fs:

    # list of different seeds for each system
    seeds = [os.path.join(f,x.replace('-param','')) for x in os.listdir(f) if x.endswith('-param')]
    
    p = get_params(seeds[0]+'-param')
    dx = p.Lx/p.Nbinx
    dy = p.Ly/p.Nbiny

    # displacement data    
    disps_ = [np.mean(np.genfromtxt(x+'-disp',delimiter='\t')[::2,1:],axis=1) for x in seeds]
    disps_ = np.array([np.concatenate([[0],x]) for x in disps_])
    disp_avg = np.mean(disps_,axis=0)
    disp_err = np.std(disps_,axis=0,ddof=1)/np.sqrt(len(seeds))

    # density profile
    profs_ = [np.genfromtxt(x+'-prof',delimiter='\t')[:,1:]/p.NstepProf*p.Lx*p.Ly/dx/dy/p.N for x in seeds]
    prof = sum(profs_)/len(seeds)

    # current profile
    if 'AOUPs/int' in f: # current profile has extra rows in AOUP data; we exclude here
        Jprofs_ = [np.genfromtxt(x+'-Jprof',delimiter='\t')[:,1:] for x in seeds]
        Jprofs_ = [x[:2*p.Nbinx,:]+x[2*p.Nbinx:,:] for x in Jprofs_]
    else:
        Jprofs_ = [np.genfromtxt(x+'-Jprof',delimiter='\t')[-2*p.Nbinx:,1:] for x in seeds]

    # x and y current
    Jxprofs_ = [x[:p.Nbinx,:] for x in Jprofs_]
    Jyprofs_ = [x[p.Nbinx:,:] for x in Jprofs_]
    Jxprof = sum(Jxprofs_)/len(seeds)
    Jyprof = sum(Jyprofs_)/len(seeds)

    # save data
    np.savetxt(os.path.join(f,'disp_avg'),disp_avg,delimiter='\t',fmt='%.4e')
    np.savetxt(os.path.join(f,'disp_err'),disp_err,delimiter='\t',fmt='%.4e')
    np.savetxt(os.path.join(f,'prof_avg'),prof,delimiter='\t',fmt='%.4e')
    np.savetxt(os.path.join(f,'Jxprof_avg'),Jxprof,delimiter='\t',fmt='%.4e')
    np.savetxt(os.path.join(f,'Jyprof_avg'),Jyprof,delimiter='\t',fmt='%.4e')


