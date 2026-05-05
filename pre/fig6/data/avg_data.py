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

supf = os.getcwd()
N_folders = [os.path.join(supf,x) for x in os.listdir(supf) if x.startswith('N') and os.path.isdir(x)]
N_folders = sorted(N_folders, key=lambda x: int(x.split('/')[-1].replace('N',''))) 

# get parameters of a sample file
sample_seed = [os.path.join(N_folders[0],x) for x in os.listdir(N_folders[0]) if 'param' in x][0]
p = get_params(sample_seed)


################################################################################################################
# AVERAGE DATA

# iterate over different N
for N_folder in N_folders:
    N = int(N_folder.split('/')[-1].replace('N',''))
    print('N =',N)

    # get list of seeds
    files = [x for x in os.listdir(N_folder) if 'disp' in x and '.png' not in x and 'batch' in x]
    
    # list of interaction strengths
    amps = list(set(list([x.split('amp-')[1].split('-')[0] for x in files])))
    amps = sorted(amps, key = float)

    # "epsilon" = amp*2/3
    # should give the same list of N*epsilon for each N
    print('\tN epsilon: '+','.join([str(float(x)*N*2/3.) for x in amps]))

    # iterate over interaction strengths
    Js_all = []
    for i,amp in enumerate(amps):
        Nepsilon = N*float(amp)*2/3.
        pf = [os.path.join(N_folder,x) for x in os.listdir(N_folder) if 'param' in x][0]
        p = get_params(pf)

        ####################################################################################################
        # AVG CURRENT

        # get list of particle displacement data files
        files_amp = [os.path.join(N_folder,x) for x in files if 'amp-'+amp in x and 'disp' in x and 'avg' not in x]

        # iterate over different seeds, get list of currents, one for each seed
        Js = []
        for k,f in enumerate(files_amp):
            disp = np.genfromtxt(f, delimiter='\t')
            Js += [np.mean(disp[-N:,2]) / p.L / p.final_time]

        # calculate and save mean and standard error of J over the seeds
        Js = np.array(Js)
        mean_J = np.mean(Js)
        err_J = np.sqrt(np.sum((Js-mean_J)**2)/(Js.shape[0]-1))/np.sqrt(Js.shape[0])
        Js_all += [[N, Nepsilon, mean_J, err_J]]

        # NOW, AVG DISPLACEMENT VS TIME
        # iterate over different seeds, get list of currents, one for each seed
        disps = []
        N_disp = int(p.final_time / p.StoreInterDisp+1e-5)
        for k,f in enumerate(files_amp):
            disp = np.genfromtxt(f, delimiter='\t')
            disps += [[np.mean(disp[j*N:(j+1)*N,2])/p.L for j in range(N_disp)]]

        # calculate and save mean and standard deviation of disp over the seeds
        disps = np.array(disps)
        mean_disps = np.mean(disps,axis=1)
        err_disps = np.std(disps,axis=1,ddof=1)/np.sqrt(disps.shape[1])
        disps_ = np.zeros((disps.shape[0],3))
        disps_[:,0] = np.linspace(0,p.final_time,N_disp+1)
        disps_[1:,1] = mean_disps
        disps_[1:,2] = err_disps

        # save displacement vs time
        np.savetxt(os.path.join(N_folder,'amp-'+amp+'-avg_disp'),disps_)


        ####################################################################################################
        # AVG DENSITY PROFILE

        # get files with density profiles for this N, amp
        files_rho = [os.path.join(N_folder,x) for x in os.listdir(N_folder) if 'amp-'+amp in x and 'profile_rho' in x and 'avg' not in x]
        
        # average over all the data
        avg_rho = 0
        counter = 0
        for f in files_rho:
            data = np.genfromtxt(f,delimiter='\t')
            times = np.unique(data[:,0])
            for t in times:
                prof_rho = data[np.where(data[:,0]==t)]
                avg_rho += prof_rho
                counter += 1
        avg_rho = avg_rho/counter
        avg_rho[:,0:2] = avg_rho[:,0:2]
        avg_rho = avg_rho[:,1:]
        
        # save density profile for this N and amp
        np.savetxt(os.path.join(N_folder,'amp-'+amp+'-avg_profile_rho'),avg_rho, delimiter='\t')

        ####################################################################################################
        # AVG MAGNETIZATION PROFILE
        
        # get files with magnetization profiles for this N, amp
        files_m = [os.path.join(N_folder,x) for x in os.listdir(N_folder) if 'amp-'+amp in x and 'profile_m' in x and 'avg' not in x]
        
        # average over all the data
        avg_m = 0
        counter = 0
        for f in files_m:
            data = np.genfromtxt(f,delimiter='\t')
            times = np.unique(data[:,0])
            for t in times:
                prof_m = data[np.where(data[:,0]==t)]
                avg_m += prof_m
                counter += 1
        avg_m = avg_m/counter
        avg_m[:,0:2] = prof_m[:,0:2]
        avg_m = avg_m[:,1:]
    
        # save magnetization profile for this N and amp
        np.savetxt(os.path.join(N_folder,'amp-'+amp+'-avg_profile_m'),avg_m, delimiter='\t')


    # save current data for this N
    np.savetxt(os.path.join('N'+str(N),'J_avg'),np.array(Js_all),delimiter='\t')