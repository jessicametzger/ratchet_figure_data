import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import os,sys
from types import SimpleNamespace
import re

plt.rcParams.update(
                {'text.usetex' : True,
                'font.size' : 8,
                'font.family' : 'serif',
                'text.latex.preamble' : r"\usepackage{lmodern} \usepackage{amstext}",
                'figure.figsize' : [3.4,3.4*0.7],
                'figure.dpi': 200
                }
)


##########################################################################################################################################
# AUXILIARY FUNCTIONS FOR READING DATA

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

# different types of temperature/potential landscapes
folders = ['data/symm_T-symm_V','data/step_T-symm_V','data/asymm_T-step_V']
ps = []
profs = []
Fprofs = []
Tprofs = []
Js = []
dTdxprofs = []
dps = []
J_stresses = []

for folder in folders:
    seeds = [os.path.join(folder,x.replace('-param','')) for x in os.listdir(folder) if x.endswith('-param') and os.path.getsize(os.path.join(folder,x.replace('-param','-prof')))>0]
    p = get_params(seeds[0]+'-param')
    
    dx = p.Lx/p.Nbinx
    
    # density, force, and temperature profiles
    prof = sum([np.genfromtxt(seed+'-prof',delimiter='\t')[:,1]/p.NstepProf/p.N/dx for seed in seeds])/len(seeds)
    Fprof = sum([np.genfromtxt(seed+'-Fprof',delimiter='\t')[:,1]/p.NstepProf/p.dt/p.N for seed in seeds])/len(seeds)
    Fprof = np.array([Fprof[:p.Nbinx],Fprof[p.Nbinx:]]).T
    Tprof = sum([np.genfromtxt(seed+'-Tprof',delimiter='\t')[:,1]/p.NstepProf/p.N for seed in seeds])/len(seeds)
    
    # displacement data
    disp = np.array([np.mean(np.genfromtxt(seed+'-disp',delimiter='\t')[:,1:],axis=1)/p.Lx for seed in seeds])
    disp_std = np.std(disp,axis=0,ddof=1)/np.sqrt(len(seeds))
    disp = np.mean(disp,axis=0)
    disp = np.array([disp[::2],disp[1::2]]).T
    disp = np.concatenate(([[0,0]],disp),axis=0)
    disp_std = np.array([disp_std[::2],disp_std[1::2]]).T
    disp_std = np.concatenate(([[0,0]],disp_std),axis=0)
    
    xs = np.linspace(0,p.Lx-dx,p.Nbinx)+dx/2.0
    ts = np.arange(0,p.tf+1e-5,p.StoreInterDisp)
    
    # T'(x) data
    dTdxprof = np.copy(Tprof)
    dTdxprof[1:-1] = (Tprof[2:]-Tprof[:-2])/2.0/dx
    dTdxprof[0] = (Tprof[1]-Tprof[-1])/2.0/dx
    dTdxprof[-1] = (Tprof[0]-Tprof[-2])/2.0/dx

    # calculate momentum sources
    if folder==folders[0]:
        x1 = (p.Fxs[0]+p.Fxs[1])/2.0
        x2 = (p.Fxs[1]+p.Fxs[2])/2.0
        x3 = (p.Fxs[2]+p.Fxs[3])/2.0
        dp1 = np.sum(Fprof[(xs>=x1) & (xs<x2)])
        dp2 = np.sum(Fprof[(xs>=x2) & (xs<x3)])
    elif folder==folders[1]:
        dp1 = None
        dp2 = None
    elif folder==folders[2]:
        xmid = (p.Fxs[3] + p.Fxs[4])/2.0
        dp1 = np.sum(Fprof[xs<xmid,0])
        dp2 = np.sum(Fprof[xs>=xmid,0])
    
    # observed current
    Jobs = disp[-1,0]/p.tf
    
    # current inferred from stresses
    binw = 10
    J_stress = (Fprof[:,0]-dTdxprof)/dx
    J_stress = sum([J_stress[i::binw] for i in range(binw)])/binw
    xsb = sum([xs[i::binw] for i in range(binw)])/binw

    ps += [p]
    profs += [prof]
    Fprofs += [Fprof]
    Tprofs += [Tprof]
    Js += [disp[-1,0]/p.tf]
    dTdxprofs += [dTdxprof]
    dps += [(dp1,dp2)]
    J_stresses += [J_stress]
    
    print(np.sum(Fprof[:,0])/p.Lx,disp[-1,0]/p.tf)


##########################################################################################################################################
#  TEMPERATURE/FORCE LANDSCAPE

# linear interpolation over points listed in xs, and y-values listed in fxs
def f_lin(xs,fxs,fs,Lx):
    Nstep = fxs.shape[0]
    fs_ = np.zeros(xs.shape)
    dfdxs = np.zeros(Nstep)
    dfdxs[:-1] = (fs[1:] - fs[:-1]) / (fxs[1:] - fxs[:-1])
    dfdxs[-1] = (fs[0] - fs[-1]) / (fxs[0] + Lx - fxs[-1])
    fs_[xs<fxs[0]] = (xs[xs<fxs[0]] + Lx - fxs[-1]) * dfdxs[-1] + fs[-1]
    for i in range(1,Nstep):
        inds = (xs<fxs[i]) & (xs>=fxs[i-1])
        fs_[inds] = (xs[inds] - fxs[i-1]) * dfdxs[i-1] + fs[i-1]
    inds = (xs>=fxs[-1])
    fs_[inds] = (xs[inds] - fxs[-1]) * dfdxs[-1] + fs[-1]
    return fs_

# cubic interpolation over points listed in xs, and y-values listed in fxs
def f_cub(xs,fxs,fs,Lx):
    cs = np.copy(fs)*0
    ds = np.copy(cs)
    Nstep = fxs.shape[0]
    for i in range(Nstep-1):
        dx = fxs[i+1]-fxs[i]
        if dx==0:
            dx=1
        cs[i] = -3*(fs[i] - fs[i+1])/dx**2
        ds[i] = 2*(fs[i]-fs[i+1])/dx**3
    dx = fxs[0]+Lx-fxs[-1]
    if dx==0:
        dx=1
    cs[-1] = -3*(fs[-1]-fs[0])/dx**2
    ds[-1] = 2*(fs[-1]-fs[0])/dx**3
    
    fs_ = np.zeros(xs.shape)

    inds = xs<fxs[0]
    dxs = xs[inds] + Lx - fxs[-1]
    fs_[inds] = fs[-1] + cs[-1]*dxs**2 + ds[-1]*dxs**3
    for i in range(1,Nstep+1):
        if i==Nstep:
            upper = Lx+fxs[0]
        else:
            upper = fxs[i]
        inds = (xs>=fxs[i-1]) & (xs<upper)
        dxs = xs[inds] - fxs[i-1]
        fs_[inds] = fs[i-1] + cs[i-1]*dxs**2 + ds[i-1]*dxs**3
    return fs_

dx_th = dx*0.02
xs_th = np.linspace(0,ps[0].Lx-dx_th,ps[0].Nbinx*50)
Fs_th = [f_lin(xs_th, p.Fxs, p.Fs, p.Lx) for p in ps]
Ts_th = [f_cub(xs_th, p.Txs, p.Ts, p.Lx) for p in ps]

Us_th = [-np.cumsum(c)*dx_th for c in Fs_th]
for i in range(3):
    Us_th[i] -= np.min(Us_th[i])

p = None


##########################################################################################################################################
# PLOTTING

# create figure/subplots
fig = plt.figure(figsize=(7,3),dpi=200)
gs = fig.add_gridspec(4,4,hspace=0.15,wspace=0.32,height_ratios=[0.5,1,0.5,0.5],width_ratios=[1,1,1,0.15])
axs = [[fig.add_subplot(gs[j,i]) for j in range(4)] for i in [0,1,2]]

# figure colors
cm = plt.get_cmap('viridis')
c_stress = cm(0)
c_ms = cm(0.99)
c_intms = cm(0.5)
c_J = cm(0.84)
c_Jobs = cm(0.35)
c_Fext = cm(0.2)
c_T = '#6aa84fff'

# line size, label spacing, etc.
lw=1.5
lbs = 0.5
htp = 0.6
hl = 1.5
bbx = 0.98
bby = 0.62
edgewidth = 0.5

# plot temperature & potential landscapes
for i in range(3):
    axs[i][0].plot(xs_th, Us_th[i], color=c_Jobs, label=r'$V(x)$')
    axs[i][0].plot(xs_th, Ts_th[i], color=c_T, label=r'$T(x)$')
axs[-1][0].legend(frameon=False, bbox_to_anchor=(bbx, bby), loc='center left', handletextpad=htp,labelspacing=lbs,handlelength=hl)

# plot thermal stress and external force profiles
for i in range(3):
    axs[i][1].fill_between(xs, Fprofs[i][:,0], lw=0, color=c_intms, alpha=0.8, label=None, zorder=3)
    axs[i][1].plot(xs, Fprofs[i][:,0], color=c_ms, label=r'$-\rho \partial_x V$',zorder=3)
    axs[i][1].plot(xs, Tprofs[i], color=c_stress, label=r'$\rho T$', zorder=4)
    axs[i][1].set_ylim(min(np.min(Fprofs[i][:,0]),np.min(Tprofs[i]))*1.5, max(np.max(Fprofs[i][:,0]),np.max(Tprofs[i]))*1.3)

# label momentum sources
i=0
axs[i][1].text(ps[i].Fxs[1]+0.5,np.min(Fprofs[i][:,0])*1.25, r'$'+'{:.2f}'.format(dps[i][0])+r'$', color=c_intms, zorder=4,va='top',ha='center',
              bbox=dict(boxstyle='square,pad='+str(.1),edgecolor='gray',facecolor='white',lw=edgewidth,alpha=0.5))
axs[i][1].text(ps[i].Fxs[2], np.max(Fprofs[i][:,0])*0.2, r'$'+'{:.2f}'.format(dps[i][1])+r'$', color=c_intms, zorder=4,va='bottom',ha='center',
              bbox=dict(boxstyle='square,pad='+str(.1),edgecolor='gray',facecolor='white',lw=edgewidth,alpha=0.5))
axs[i][1].set_ylim(ymin=np.min(Fprofs[i][:,0])*2.5)
i=-1
axs[i][1].text((ps[i].Fxs[0]+ps[i].Fxs[3])/2.0,-np.max(Fprofs[i][:,0])*0.1, r'$'+'{:.2f}'.format(dps[i][0])+r'$', color=c_intms, zorder=4,va='top',ha='center')
axs[i][1].text((ps[i].Fxs[4]+ps[i].Fxs[7])/2.0-3, np.max(Fprofs[i][:,0])*0.1, r'$'+'{:.2f}'.format(dps[i][1])+r'$', color=c_intms, zorder=4,va='bottom',ha='center')
axs[-1][1].legend(frameon=False, bbox_to_anchor=(bbx, bby), loc='center left', handletextpad=htp,labelspacing=lbs,handlelength=hl)

# plot forces
for i in range(3):
    axs[i][2].plot([-1,xs[-1]+1],[0,0],color='k',lw=0.5)
    axs[i][2].plot(xs, Fprofs[i][:,0], color=c_ms, label=r'$-\rho \partial_x V$')
    axs[i][2].plot(xs, -dTdxprofs[i], color=c_stress, ls='--', label=r'$-\partial_x (\rho T)$')
axs[-1][2].legend(frameon=False, bbox_to_anchor=(bbx, bby), loc='center left', handletextpad=htp,labelspacing=lbs,handlelength=hl)

# plot current, from observing particle trajectories & from measuring forces
for i in range(3):
    axs[i][3].plot(xsb,J_stresses[i],color=c_J, lw=lw, label=r'$J_{\rm stress}$',zorder=1)
    axs[i][3].plot(xs,[Js[i]]*xs.shape[0],color=c_Jobs,label=r'$J_{\rm obs}$',ls='--')
    axs[i][3].set_xlabel(r'$x$')
    axs[i][3].xaxis.set_label_coords(0.94,-0.075)
    axs[i][3].tick_params(axis='x',which='both',direction='in',pad=1,length=3)
    axs[i][3].set_ylim(0,np.abs(Js[i])*2)
    axs[i][3].plot([-10,1e3],[0,0],lw=lw*0.5,color='black')
axs[-1][3].legend(frameon=False, bbox_to_anchor=(bbx, bby), loc='center left',
              handletextpad=htp,labelspacing=lbs,handlelength=hl)

# axis label/limit/etc.
for i in range(3):
    for j in range(4):
        axs[i][j].set_xlim(0,ps[i].Lx)
        axs[i][j].tick_params(axis='y', which='both', length=2, direction='in',pad=2)
    for j in range(3):
        axs[i][j].set_xticks([])
    axs[i][0].set_ylim(0,max(np.max(Ts_th[i]),np.max(Us_th[i]))*1.5)
    axs[i][-1].set_xticks([0,10,20,30])
    

# add (a), (b), etc.
pad = 1.25
fontsize=8.0
labels=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)']
count=0
for i in range(3):
    for j in range(4):
        bbox = axs[i][j].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        pad_x = pad / 72 / width
        pad_y = pad / 72 / height
        inside_pad = pad/fontsize
        text = axs[i][j].text(pad_x,1-pad_y,labels[count],ha='left',va='top',transform=axs[i][j].transAxes,
                              bbox=dict(boxstyle='square,pad='+str(inside_pad),edgecolor='black',facecolor='white',lw=edgewidth))
        count += 1
        text.set_clip_on(True)

plt.subplots_adjust(top=0.93,bottom=0.04,right=0.97,left=0.06)
plt.savefig('blowtorch_stress_appendix_fig-tf_{:d}-dt_{:.6f}.pdf'.format(int(ps[i].tf),ps[i].dt),dpi=1000)
