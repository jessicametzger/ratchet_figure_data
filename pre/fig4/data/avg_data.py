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

# get seed (assume only 1)
folder = os.getcwd()
seed = [os.path.join(folder,x.replace('-param','')) for x in os.listdir(folder) if x.endswith('-param')]
seed = seed[0]

# get parameter info 
p = get_params(seed+'-param')

# read epr and stress data
epr = np.genfromtxt(seed+'-EPR',delimiter='\t') # epr
sig = np.genfromtxt(seed+'-sigmaIKprof',delimiter='\t') # irving-kirkwood interaction stress
Tst = np.genfromtxt(seed+'-Tprof',delimiter='\t') # thermal stress


################################################################################################################
# CREATE EPR HEATMAP

# average epr over different times
NstoreProf = int(p.tf / p.StoreInterProf + 1e-5)
epr_avg = sum([epr[i*p.Nbinx:(i+1)*p.Nbinx,1:] for i in range(NstoreProf)]) / NstoreProf

# divide by normalization
dx_epr = p.Lx/p.Nbinx
dy_epr = p.Ly/p.Nbiny
T_meas = NstoreProf*p.dt*p.NstepProf # effective amount of time over which EPR is measured
unit_vol = dx_epr*dy_epr
epr_avg = epr_avg / (T_meas * unit_vol * p.N)

# save heatmap
np.savetxt('eprprof_avg',epr_avg,delimiter='\t',fmt='%.4e')


################################################################################################################
# CREATE EPR TIME SERIES

# sum at each time point
epr_series = np.array([np.sum(epr[i*p.Nbinx:(i+1)*p.Nbinx,1:]) for i in range(NstoreProf)])

# divide by normalization
epr_series = epr_series/(T_meas*p.N*unit_vol)

# save time series
np.savetxt('eprseries_avg',epr_series,delimiter='\t',fmt='%.4e')


################################################################################################################
# CREATE STRESS PROFILES

Lx = int(p.Lx)

# sum the stresses over time
sig = sum([sig[i*Lx*3:(i+1)*Lx*3,1:] for i in range(NstoreProf)]) 
Tst = sum([Tst[i*p.Nbinx:(i+1)*p.Nbinx,1:] for i in range(NstoreProf)]) 

# bin array to match same grid as sigmaIK (grid width = 1 rather than dx=Lx/Nbinx)
dx = p.Lx / p.Nbinx
binw = int(1 / dx + 1e-5)
Tst_coarse = sum([Tst[i::binw,:] for i in range(binw)])
Tst_coarse = sum([Tst_coarse[:,i::binw] for i in range(binw)])

# normalize
sig = sig / (NstoreProf*p.NstepProf*p.N*p.dt)
Tst_coarse = Tst_coarse / (NstoreProf*p.NstepProf*p.N)

# average the irving-kirkwood and thermal stresses over y
sig = np.sum(sig,axis=1) / p.Ly
Tst_coarse = np.sum(Tst_coarse,axis=1) / p.Ly

# create single array
stress = np.zeros((Lx, 5))
stress[:,0] = np.linspace(0,Lx-1,Lx) # x points
stress[:,1] = sig[:Lx] # xx component of IK stress
stress[:,2] = sig[Lx:Lx*2] # xy component of IK stress
stress[:,3] = sig[Lx*2:] # yy component of IK stress
stress[:,4] = Tst_coarse # thermal stress

# save stress profiles
np.savetxt('sigprof_avg',stress,delimiter='\t',fmt='%.4e')