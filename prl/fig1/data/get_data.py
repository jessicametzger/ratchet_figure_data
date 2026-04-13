'''
Get the simulation data, parse it, and put it in new files
'''

import numpy as np
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


path = '../C_code'
flags = ['ABP-int','ABP-ni','PBP-int','PBP-ni']
paths = [os.path.join(path,f) for f in flags]

# get parameter object
ps = [get_params(x+'-param') for x in paths]

# get displacement data
disps = [np.genfromtxt(x+'-disp',delimiter='\t') for x in paths]

# get profile data
profs = [np.genfromtxt(x+'-prof',delimiter='\t') for x in paths]

# output destinations
out_paths = ['ABPs/active-Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025',
			 'ABPs/active-ni_Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025',
			 'PBPs/passive-Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025',
			 'PBPs/passive-ni_Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025']

# save output
for k,p in enumerate(ps):
	# save profile data
	np.savetxt(out_paths[k]+'-prof_avg', profs[k], delimiter='\t', fmt='%.4e')

	# save average displacement data
	disp_ = np.zeros((int(p.tf/p.StoreInterDisp+1e-5),2))
	disp_[:,0] = disps[k][::2,0]
	disp_[:,1] = np.mean(disps[k][::2,1:],axis=1) # average over the different particles
	np.savetxt(out_paths[k]+'-disp_avg', disp_, delimiter='\t', fmt='%.4e')

	# copy param info
	f=open(paths[k]+'-param','r')
	s=f.read()
	f.close()

	f=open(out_paths[k]+'-param','w')
	f.write(s)
	f.close()

