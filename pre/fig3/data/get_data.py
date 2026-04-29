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

folder = os.getcwd()

# get list of discretizations
discs = [os.path.join(folder,x) for x in os.listdir(folder) if os.path.isdir(os.path.join(folder,x)) and not x.startswith('.')]

for disc in discs:

    # get list of seeds for this discretization
    seeds = [os.path.join(disc,x.replace('-param','')) for x in os.listdir(disc) if x.endswith('-param')]

    # get parameter data
    p = get_params(seeds[0]+'-param')

    # displacement data    
    disps_ = [np.mean(np.genfromtxt(x+'-disp',delimiter='\t')[::2,1:],axis=1) for x in seeds]
    disps_ = np.array([np.concatenate([[0],x]) for x in disps_])
    disp_avg = np.mean(disps_,axis=0)
    disp_err = np.std(disps_,axis=0,ddof=1)/np.sqrt(len(seeds))

    # density profile
    profs_ = [np.genfromtxt(x+'-prof',delimiter='\t')[:,1:]/p.NstepProf/p.N/(p.Lx/p.Nbinx)/(p.Ly/p.Nbiny) for x in seeds]
    prof = sum(profs_)/len(seeds)

    # save data
    np.savetxt(os.path.join(disc,'disp_avg'),disp_avg,delimiter='\t',fmt='%.4e')
    np.savetxt(os.path.join(disc,'disp_err'),disp_err,delimiter='\t',fmt='%.4e')
    np.savetxt(os.path.join(disc,'prof_avg'),prof,delimiter='\t',fmt='%.4e')
    
