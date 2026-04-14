import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import os,sys
from types import SimpleNamespace
import re


################################################################################################################
# AUXILIARY FUNCTIONS

def isnum(s):
    try:
        a=float(s)
    except ValueError:
        return False
    return True

# import parameter data from a simulation
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


################################################################################################################
# READ DATA

# get list of seeds
folder = os.getcwd()
seeds = [os.path.join(folder,x.replace('-param','')) for x in os.listdir(folder) if x.endswith('-param') and os.path.getsize(os.path.join(folder,x.replace('param','prof')))>0]

# get parameter info from one seed
p = get_params(seeds[0]+'-param')

# reed profile and displacement data
profs = np.array([np.genfromtxt(seed+'-prof',delimiter='\t')[-p.Nbinx:,1] / (p.NstepProf) for seed in seeds])
FAprofs = np.array([np.genfromtxt(seed+'-FAprof',delimiter='\t')[:,1] / (p.dt*p.NstepProf) for seed in seeds])
Fintprofs = np.array([np.genfromtxt(seed+'-Fintprof',delimiter='\t')[:,1] / (p.dt*p.NstepProf) for seed in seeds])
sigIKprofs = np.array([np.sum(np.genfromtxt(seed+'-sigmaIKprof',delimiter='\t')[:,1:],axis=1) / (p.dt*p.NstepProf) for seed in seeds])
sigAprofs = np.array([np.genfromtxt(seed+'-sigmaAprof',delimiter='\t')[:,1] / (p.dt*p.dt*p.NstepProf) for seed in seeds])
disps = np.array([np.sum(np.genfromtxt(seed+'-disp',delimiter='\t')[::2,1:],axis=1) for seed in seeds])

# average over seeds
prof = np.mean(profs,axis=0)
FAprof = np.mean(FAprofs,axis=0)
Fintprof = np.mean(Fintprofs,axis=0)
sigIKprof = np.mean(sigIKprofs,axis=0)
sigAprof = np.mean(sigAprofs,axis=0)
disp = np.mean(disps,axis=0)

# collate arrays
FAprof = np.transpose(np.array([FAprof[i*p.Nbinx:(i+1)*p.Nbinx] for i in range(2)]),(1,0))
Fintprof = np.transpose(np.array([Fintprof[i*p.Nbinx:(i+1)*p.Nbinx] for i in range(2)]),(1,0))
sigIKprof = np.transpose(np.array([sigIKprof[i*int(p.Lx):(i+1)*int(p.Lx)] for i in range(3)]),(1,0))
sigAprof = np.transpose(np.array([sigAprof[i*p.Nbinx:(i+1)*p.Nbinx] for i in range(4)]),(1,0))

# save averaged data
np.savetxt(os.path.join(folder,'FAprof_avg'),FAprof,delimiter='\t',fmt='%.4e')
np.savetxt(os.path.join(folder,'Fintprof_avg'),Fintprof,delimiter='\t',fmt='%.4e')
np.savetxt(os.path.join(folder,'sigIKprof_avg'),sigIKprof,delimiter='\t',fmt='%.4e')
np.savetxt(os.path.join(folder,'sigAprof_avg'),sigAprof,delimiter='\t',fmt='%.4e')
np.savetxt(os.path.join(folder,'disp_avg'),disp,delimiter='\t',fmt='%.4e')
