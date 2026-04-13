import os,sys
import numpy as np
from types import SimpleNamespace
import re

sys.path.append(os.getcwd())
import mydefaults as md


def genfromtxt(fpath,delimiter='\t',dtype=float,hskip=0):
    f=open(fpath,'r')
    data = np.array([x.strip().split(delimiter) for x in f.readlines()[hskip:]]).astype(dtype)
    f.close()
    return data


def getp(s,ps,st=False):
    after = s.split(ps+'_')[1]
    if after[0]=='-':
        ans = '-'+after.split('-')[1]
    else:
        ans = after.split('-')[0]
    if st:
        return ans.split('/')[0]
    return float(ans.split('/')[0])

def get(paramstr, param, keepstr=False):
    s = paramstr.split('\n'+param+' is ')[-1].split('\n')[0]
    if keepstr:
        return s
    return float(s)

def isnum(s):
    try:
        a=float(s)
    except ValueError:
        return False
    return True

def sci(n):
    if n == 0:
        return "0"
    exponent = int(np.floor(np.log10(np.abs(n)))+md.EPS)
    mantissa = n / 10**exponent
    if mantissa == 1:
        return f"$10^{{{exponent}}}$"
    else:
        return f"${mantissa:.0f} \\times 10^{{{exponent}}}$"

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


def avg_disp(files,p, delimiter='\t'):
    '''
    Will return the avg and stderr displacement over time for a list of files

    assume col1 is t, col2 is ids, etc.

    doesn't divide by L
    '''
    if hasattr(p,'tf'):
        tf = p.tf
    elif hasattr(p,'final_time'):
        tf = p.final_time
    else:
        raise RuntimeError("No tf attribute")
    Nstore = int(tf / p.StoreInterDisp + md.EPS)
    nseed = len(files)
    norm = nseed
    data = [np.genfromtxt(x,delimiter=delimiter) for x in files]

    examp = data[0]
    ids = np.unique(examp[:,1])
    N = ids.shape[0]
    ts= np.unique(examp[:,0])

    if ts.shape[0]!=Nstore:
        raise RuntimeError("observing "+str(ts.shape[0])+" time slices but expecting Nstore="+str(Nstore))
    for i,x in enumerate(data):
        if N*Nstore!=x.shape[0]:
            raise RuntimeError("file "+files[i]+" has shape "+str(x.shape[0])+" but expected Nstore="+str(Nstore))

    s=data[0].shape[1]
    if s==3:
        # avg over N
        data_avgd_N = np.concatenate([[[sum(x[i*N:(i+1)*N,2])/float(p.N) for i in range(Nstore)]] for x in data],axis=0)
        
        # sum over seeds
        data_avgd = data_avgd_N/norm
    
        # t, avg, stderr
        ret = np.zeros((Nstore,3))
        ret[:,0] = ts
        ret[:,1] = np.mean(data_avgd,axis=0)
        ret[:,2] = np.std(data_avgd,axis=0,ddof=1)/np.sqrt(nseed)
    else:
        # avg over N
        data_avgd_N = np.concatenate([[[sum(x[i*N:(i+1)*N,2:])/float(p.N) for i in range(Nstore)]] for x in data],axis=0)
        
        # sum over seeds
        data_avgd = data_avgd_N/norm
        
        # t, avg, stderr for each dim
        ret = np.zeros((Nstore+1,1+2*(s-2)))
        ret[1:,0] = ts
        ret[1:,1:(s-1)] = np.mean(data_avgd,axis=0)
        ret[1:,(s-1):] = np.std(data_avgd,axis=0,ddof=1)/np.sqrt(nseed)
        
    return ret
    


def avg_prof_1d(files,p,delimiter='\t'):
    '''
    Will read and average all files in list of filenames ls

    Assumes column 0 is time, columns 1 is x, everything else to be averaged
    '''
    if hasattr(p,'tf'):
        tf = p.tf
    elif hasattr(p,'final_time'):
        tf = p.final_time
    else:
        raise RuntimeError("No tf attribute")
    Nstore = int(tf / p.StoreInterProf + md.EPS)
    Nstep = int(p.NstepProf)
    norm = Nstore*Nstep

    nseed = len(files)
    norm*= nseed
    
    data = [np.genfromtxt(x,delimiter=delimiter) for x in files]
    
    examp = data[0]
    xs = np.unique(examp[:,1])
    size = xs.shape[0]

    # make sure there is the expected amount of data 
    assert np.all([p.Nbin*Nstore==x.shape[0] for x in data])

    # sum all the arrays together
    data_tot = sum([x[:,2:] for x in data])

    # sum it over time
    data_avg = sum([data_tot[size*i:size*(i+1),:] for i in range(Nstore)])
    ret = np.zeros((size,data_avg.shape[1]+1))
    ret[:,0] = xs
    ret[:,1:]= data_avg/norm
    return ret



def avg_prof_2d(files,p,avg_y=False,delimiter='\t', i0=3, dx=None, dy=None, Nstep=None, StoreInter=None): #Nstore,Nstep=1,avg_y=True,Nbinx=None,Nbiny=None
    '''
    Will read and average all files in list of filenames ls

    Assumes column 0 is time, columns 1 is x, column 2 is y, and everything else (after index i0 with default 3) to be averaged

    p is "params" object which includes items Lx, Ly, Nbinx, Nbiny, final_time, StoreInterProfile, and NstepProfile
    '''
    if hasattr(p,'tf'):
        tf = p.tf
    elif hasattr(p,'final_time'):
        tf = p.final_time
    else:
        raise RuntimeError("No tf attribute")
    if StoreInter is None:
        Nstore = int(tf / p.StoreInterProf + md.EPS)
    else:
        Nstore = int(tf / StoreInter + md.EPS)

    if Nstep is None:
        Nstep = int(p.NstepProf)
    else:
        Nstep = int(Nstep)

    if dx is None:
        Nbinx = p.Nbinx
        Nbiny = p.Nbiny
        dx = p.Lx/p.Nbinx
        dy = p.Ly/p.Nbiny
    else:
        Nbinx = int(p.Lx / dx + 1e-10)
        Nbiny = int(p.Ly / dy + 1e-10)
    xs = np.linspace(0,p.Lx-dx, Nbinx)
    ys = np.linspace(0,p.Ly-dy, Nbiny)
    
    # keep a running total of normalization
    norm = Nstore*Nstep

    nseed = len(files)
    norm*= nseed
    
    data = [np.genfromtxt(x,delimiter=delimiter) for x in files]
    
    N = Nbinx
    expected_size = N*Nstore

    # make sure there is the expected amount of data 
    for i in range(len(files)):
        x=data[i]
        if x.shape[0]!=expected_size:
            raise RuntimeError('file '+files[i]+' has length '+str(x.shape[0])+' but expected '+str(expected_size))

    # sum all the array data together
    data_tot = sum([x[:,i0:] for x in data])

    # sum it over time
    data_avg = sum([data_tot[N*i:N*(i+1),:] for i in range(Nstore)])

    # collate in y
    data_avg = np.concatenate([[data_tot[Nbiny*i:Nbiny*(i+1),:]] for i in range(Nbinx)],axis=0)

    # sum it over y
    if avg_y:
        norm *= Nbiny
        data_avg = np.sum(data_avg,axis=1)
        
        ret = np.zeros((xs.shape[0],data_avg.shape[1]+1))
        ret[:,0] = xs
        ret[:,1:]= data_avg/norm
        return ret

    else:
        return data_avg / norm
        
    



# def avg_prof_ddim(files,Nstore,Nstep=1,d=1,avg_dims=None,delimiter='\t'):
#     '''
#     Will read and average all files in list of filenames ls

#     Assumes column 0 is time, columns 1,...,d are coordinates

#     Will average over the avg_dims dimensions (x=1, y=2, ...)

#     Will assume each profile consists of Nstep measurements, so they're divided by Nstep
#     '''
#     # keep a running total of normalization
#     norm = Nstore*Nstep

#     nseed = len(files)
#     norm*= nseed
    
#     data = [np.genfromtxt(x,delimiter=delim) for x in files]
    
#     examp = data[0]
#     coords = [np.unique(examp[:,i+1]) for i in d]
#     sizes = [x.shape[0] for x in coords]

#     # make sure there is the expected amount of data 
#     assert np.all([np.prod(sizes)*Nstore==x.shape[0] for x in data])

#     # sum all the arrays together
#     data_tot = sum(data) 

#     dims_all = range(1,len(sizes)+1)
#     dims_keep = [x for x in dims_all if x not in avg_dims]
#     sizes_keep = list(np.array(sizes)[dims_keep])
#     avgd_arr = np.zeros(tuple(sizes_keep))

#     # NEED TO FINISH THIS LATER