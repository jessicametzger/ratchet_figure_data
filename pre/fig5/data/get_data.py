import numpy as np
import os
        
        
def get(paramstr, param, keepstr=False):
    s = paramstr.split(param+' is ')[-1].split('\n')[0]
    if keepstr:
        return s
    return float(s)



supdir = 'N_500-amp_0.2-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12'
folders = [x for x in os.listdir(supdir) if (x.startswith('abp_') or x.startswith('rtp_')) and os.path.isdir(supdir+'/'+x)]

amp = float(supdir.split('amp_')[-1].split('-')[0])

prof_avgs = {'abp': None, 'rtp': None}
disp_avgs = {'abp': None, 'rtp': None}
profs = {'abp': [], 'rtp': []}
disps = {'abp': [], 'rtp': []}
for i,folder in tqdm(enumerate(folders)):
    path = supdir+'/'+folder+'/res-'
    if 'abp' in folder:
        flag = 'abp'
    elif 'rtp' in folder:
        flag = 'rtp'
    
    if not os.path.exists(path+'param'):
        continue
    
    f=open(path+'param')
    param=f.read()
    f.close()

    Lx=get(param,'Lx')
    Ly=get(param,'Ly')
    N=int(get(param,'N'))
    Nbinx=int(get(param,'Nbinx'))
    
    disp = np.genfromtxt(path+'disp',delimiter='\t')
    if disp.shape[0]<10*N:
        continue
    prof = np.genfromtxt(path+'profile_rho',delimiter='\t')
    
    
    Ts = int(prof.shape[0]/Nbinx+1e-5)
    prof_avg = np.zeros((Nbinx,2))
    prof_avg[:,0] = prof[0:Nbinx,1]
    prof_avg[:,1] = sum([prof[Nbinx*i:Nbinx*(i+1),3] for i in range(Ts)])/Ts
    
    profs[flag] += [prof_avg]
    if prof_avgs[flag] is None:
        prof_avgs[flag] = prof_avg
    else:
        prof_avgs[flag][:,1] += prof_avg[:,1]
        
    Ts = int(disp.shape[0]/N+1e-5)
    disp_avg = np.zeros((Ts,2))
    disp_avg[:,0] = disp[::N,0]
    disp_avg[:,1] = sum([disp[i::N,2] for i in range(N)])/(N*Lx)
    
    disps[flag] += [disp_avg.copy()]
    if disp_avgs[flag] is None:
        disp_avgs[flag] = disp_avg.copy()
    else:
        disp_avgs[flag][:,1] += disp_avg.copy()[:,1]
        
for flag in ['abp','rtp']:
    disp_avgs[flag] = np.concatenate((np.zeros((1,2)),disp_avgs[flag]),axis=0)
    for i in range(len(disps[flag])):
        disps[flag][i] = np.concatenate((np.zeros((1,2)),disps[flag][i]),axis=0)
        
    prof_avgs[flag][:,1] = prof_avgs[flag][:,1] / np.trapz(prof_avgs[flag][:,1], prof_avgs[flag][:,0])
    for i in range(len(profs[flag])):
        profs[flag][i][:,1] = profs[flag][i][:,1] / np.trapz(profs[flag][i][:,1], profs[flag][i][:,0])

    np.savetxt(supdir+'/'+flag+'_prof_avg',prof_avgs[flag])


amps = [0, 0.025, 0.05,0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225]

Js_arr = []
for amp in tqdm(amps):
    #supdir = 'N_700-amp_'+str(amp)+'-Lx_60-Ly_15-v0_20-v1_1-landscape_10_20_50'
    supdir = 'N_500-amp_'+str(amp)+'-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12'
    folders = [x for x in os.listdir(supdir) if (x.startswith('abp_') or x.startswith('rtp_')) \
               and os.path.isdir(supdir+'/'+x) and os.path.exists(supdir+'/'+x+'/res-param')]
    
    amp = float(supdir.split('amp_')[-1].split('-')[0])

    Js = {'abp': [], 'rtp': []}
    for i,folder in enumerate(folders):
        path = supdir+'/'+folder+'/res-'
        if 'abp' in folder:
            flag = 'abp'
        elif 'rtp' in folder:
            flag = 'rtp'

        if not os.path.exists(path+'param'):
            continue

        f=open(path+'param')
        param=f.read()
        f.close()

        Lx=get(param,'Lx')
        N=int(get(param,'N'))

        disp = np.genfromtxt(path+'disp',delimiter='\t')
        Js[flag] += [np.mean(disp[-N:,2])/(Lx*disp[-1,0])]
    
    Js_mean = {'abp': None, 'rtp': None}
    Js_std = {'abp': None, 'rtp': None}
    for flag in ['abp','rtp']:
        Js[flag] = np.array(Js[flag])
        Js_mean[flag] = np.mean(Js[flag])
        Js_std[flag] = np.std(Js[flag],ddof=1)/np.sqrt(Js[flag].shape[0])
    Js_arr += [[amp, Js_mean['abp'], Js_mean['rtp'], Js_std['abp'], Js_std['rtp']]]
    
Js_arr = np.array(Js_arr)
print(Js_arr)
np.savetxt('weak_disp_avg',Js_arr)