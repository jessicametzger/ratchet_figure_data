import numpy as np
import time
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
# READ IN PARAMETER INFO

supf = '/home/jessica/Projects/Ratchet/1D_RTPs/strong_current/new_L10_dt0.0005_v020_a11_a20.8'
N_folders = [os.path.join(supf,x) for x in os.listdir(supf) if x.startswith('N') and os.path.isdir(x)]
N_folders = sorted(N_folders, key=lambda x: int(x.split('/')[-1].replace('N',''))) 

# get parameters of a sample file
sample_seed = [os.path.join(N_folders[0],x) for x in os.listdir(N_folders[0]) if 'param' in x][0]
p = get_params(sample_seed)


################################################################################################################
# AVERAGE DISPLACEMENT OVER TIME

# iterate over different N
for N_folder in N_folders:
    N = int(N_folder.split('/')[-1].replace('N',''))
    print('N =',N)

    files = [x for x in os.listdir(N_folder) if 'disp' in x and '.png' not in x and 'batch' in x]
    amps = list(set(list([x.split('amp-')[1].split('-')[0] for x in files])))
    amps = sorted(amps, key = float)

    # "epsilon" = amp*2/3
    # should give the same list of N*epsilon for each N
    print('\tN epsilon: '+','.join([str(float(x)*N*2/3.) for x in amps]))

    # iterate over interaction strengths
    Js_all = []
    for i,amp in enumerate(amps):
        Nepsilon = N*float(amp)*2/3.

        # get list of particle displacement data files
        files_amp = [os.path.join(N_folder,x) for x in files if 'amp-'+amp in x and 'disp' in x]
        pf = [os.path.join(N_folder,x) for x in os.listdir(N_folder) if 'param' in x][0]
        p = get_params(pf)

        # iterate over different seeds, get list of currents, one for each seed
        Js = []
        for k,f in enumerate(files_amp):
            disp = np.genfromtxt(f, delimiter='\t')
            Js += [np.mean(disp[-N:,2]) / p.L / p.final_time]


        # calculate and save mean and standard deviation of J over the seeds
        Js = np.array(Js)
        mean_J = np.mean(Js)
        std_J = np.sqrt(np.sum((Js-mean_J)**2)/(Js.shape[0]-1))/np.sqrt(Js.shape[0])
        Js_all += [[N, Nepsilon, mean_J, std_J]]

    # save current data for this N
    np.savetxt(os.path.join('N'+str(N),'J_avg'),np.array(Js_all),delimiter='\t')