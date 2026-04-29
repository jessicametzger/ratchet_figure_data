'''
Create the figure, using the averaged data from data/get_data.py
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os,sys
from types import SimpleNamespace
import re
from scipy.ndimage import convolve

plt.rcParams.update(
                {'text.usetex' : True,
                'font.size' : 8,
                'font.family' : 'serif',
                'text.latex.preamble' : r"\usepackage{lmodern} \usepackage{amstext}",
                'figure.figsize' : [3.4,3.4*0.7],
                'figure.dpi': 200
                }
)

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


# bin 2d array
def bin_arr(arr,binw):
    arr_=np.copy(arr)
    arr_ = sum([arr_[i::binw,:] for i in range(binw)])/binw
    arr_ = sum([arr_[:,i::binw] for i in range(binw)])/binw
    return arr_


##########################################################################################################################################
# READ DATA

folder = 'data'

# get list of discretizations
discs = [os.path.join(folder,x) for x in os.listdir(folder) if os.path.isdir(os.path.join(folder,x)) and not x.startswith('.')]

# get data for each discretization; put into dict
disp_avgs = {x.split('/')[-1]: None for x in discs}
disp_errs = {x.split('/')[-1]: None for x in discs}
prof_avgs = {x.split('/')[-1]: None for x in discs}
ps = {x.split('/')[-1]: None for x in discs}
for disc in discs:
    p = get_params([os.path.join(disc,x) for x in os.listdir(disc) if x.endswith('-param')][0])
    ps[disc.split('/')[-1]] = p
    disp_avgs[disc.split('/')[-1]] = np.genfromtxt(os.path.join(disc,'disp_avg'),delimiter='\t')
    disp_errs[disc.split('/')[-1]] = np.genfromtxt(os.path.join(disc,'disp_err'),delimiter='\t')
    prof_avgs[disc.split('/')[-1]] = np.genfromtxt(os.path.join(disc,'prof_avg'),delimiter='\t')

ts = np.linspace(0,p.tf,11)


##########################################################################################################################################
# TEMPERATURE LANDSCAPE

# cubic activity landscape function    
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

dx = p.Lx/p.Nbinx
xs_data = np.arange(0,p.Lx-1e-5,dx)
Ts_data = f_cub(xs_data, p.Txs, p.Ts, p.Lx)



##########################################################################################################################################
# J THEORY

# spatial discretization
dx = 0.00025
xs_th = np.arange(0,p.Lx-1e-5,dx)+dx/2.0

# temperature landscape with this spatial discretization
Ts = f_cub(xs_th, p.Txs, p.Ts, p.Lx)

# harmonic interaction force kernel
xs_ker = np.arange(-1, 1, dx)+dx/2.
Fs = xs_ker.copy()
Fs[xs_ker>0] =  1-xs_ker[xs_ker>0]
Fs[xs_ker<0] = -1-xs_ker[xs_ker<0]

# loop over different discretization alpha; calculate linear coefficient of J
alphas_plot = np.linspace(0,1,20)
J1s = []
for alpha in alphas_plot:
    Ts_alpha = Ts**alpha
    Ts_1minusalpha = Ts**(1-alpha)
    
    c = 1/np.sum(dx/Ts_alpha)
    kappa = 1 / np.sum(dx/Ts_1minusalpha)
    rho = kappa/Ts_1minusalpha
    rho_conv = convolve(rho,-Fs)*dx

    J1s += [-c * dx * np.sum(kappa*rho_conv/Ts)]

J1s = np.array(J1s)


##########################################################################################################################################
# PLOTTING

# colors
cm = plt.get_cmap('plasma')
def color(disc):
    disc_ = float(disc)
    return cm(disc_*0.92)
    
# create figure and subplots
fig = plt.figure(figsize=(7, 7*.22),dpi=600)
gs = gridspec.GridSpec(2, 3, height_ratios=[.5, 1], hspace=0, wspace=0.13)
ax1 = fig.add_subplot(gs[1,0])
ax0 = fig.add_subplot(gs[0,0],sharex=ax1)
ax2 = fig.add_subplot(gs[:,1])
ax3 = fig.add_subplot(gs[:,2])
axs = [ax0,ax1,ax2,ax3]

# plot temperature landscape
ax0.plot(xs_th, Ts, color='green')
ax0.set_ylabel(r'$T(x)$',rotation='horizontal')
ax0.yaxis.set_label_coords(-0.08,0.7)
ax0.set_ylim(ymin=0,ymax=35)
ax0.set_yticks([0,20])
ax0.tick_params(axis='x', which='both', length=0, labelbottom=False)

# plot density profile data; only include legend for alternating alpha
discs = [x.split('/')[-1] for x in discs]
for disc in sorted(discs):
    label=None
    if disc=='0.0':
        label=r'$\alpha$=0.0'
    elif disc=='0.25':
        label=r'$\alpha$=0.25'
    elif disc=='0.5':
        label=r'$\alpha$=0.5'
    elif disc=='0.625':
        ax1.plot([-1,-1],[0,1],color='white',label=r'$\,$',alpha=0)
    elif disc=='0.75':
        label=r'$\alpha$=0.75'
    elif disc=='1.0':
        label=r'$\alpha$=1'
    ax1.plot(xs_data, prof_avgs[disc], color=color(disc),label=label)
ax1.set_xlim(xmin=0,xmax=p.Lx)
ax1.set_ylim(ymin=0,ymax=np.max(prof_avgs['0.0'])*1.3)
ax1.legend(frameon=False,ncol=2,loc='lower left',
          handletextpad=0.4,labelspacing=0.25,handlelength=1.5, columnspacing=1.5)
ax1.set_ylabel(r'$\rho_{\rm s}(x)$',rotation='horizontal')
ax1.yaxis.set_label_coords(-0.08,0.75)
ax1.set_yticks([0,.025,.05])
ax1.set_xticks([0,5,10,15])
ax1.set_xlabel(r'$x$')
ax1.xaxis.set_label_coords(0.93,-0.02)

# plot displacement vs. time data (only for alternating alpha)
ax2.plot([-1,p.tf+1],[0,0],color='black',lw=0.5)
for disc in ['0.0','0.25','0.5',None,'0.75','1.0']:
    if disc is None:
        ax2.plot([-1,-1],[0,1],color='white',label=' ',alpha=0)
        continue
    ax2.fill_between(ts,
                     (disp_avgs[disc]-3*disp_errs[disc])/p.Lx, 
                     (disp_avgs[disc]+3*disp_errs[disc])/p.Lx,
                     color=color(disc),zorder=1, alpha=0.25,linewidth=0)
    ax2.plot(ts, disp_avgs[disc]/p.Lx, label=r'$\alpha$='+disc, color=color(disc),zorder=2)
    ax2.errorbar(ts, disp_avgs[disc]/p.Lx,
                 color=color(disc),zorder=2,marker='o',markersize=3) #label=r'$\alpha$='+disc, 
ax2.legend(frameon=False,ncol=2,loc='lower center',bbox_to_anchor=(0.4,-0.025),
          handletextpad=0.4,labelspacing=0.25,handlelength=1.5, columnspacing=1.5)
ax2.set_xlim(xmin=0,xmax=p.tf)
ax2.set_ylim(ymax=4,ymin=-14)
ax2.set_yticks([0,-5,-10],labels=['0',r'$-5$',r'$-10$'])
ax2.set_xticks([0,50000,100000,150000])
ax2.set_xlabel(r'$t$')
ax2.set_ylabel(r'$\frac{\Delta x}{L_x}$',rotation='horizontal')
ax2.yaxis.set_label_coords(-0.06,0.9)
ax2.xaxis.set_label_coords(0.95,-0.02)

# plot currrent vs. alpha
for i,disc in enumerate(discs):
    ax3.errorbar([float(disc)], [disp_avgs[disc][-1]/p.Lx/p.tf], 
                yerr=[3*disp_errs[disc][-1]/p.Lx/p.tf],
                marker='o',ls='',color=color(disc))#,markeredgecolor='black',
                # markeredgewidth=0.75)
ax3.plot([-0.1,1.1],[0,0],color='black',lw=0.5)
ax3.plot(alphas_plot, [J1 * p.k / p.Ly for J1 in J1s], color='black')
ax3.set_xlim(xmin=-0.05,xmax=1.08)
ax3.set_ylim(ymax=1.75e-5)
ax3.set_xticks([0,0.25,0.5,0.75,1],labels=['0','0.25','0.5','0.75','1'])
ax3.set_xlabel(r'$\alpha$')
ax3.xaxis.set_label_coords(.99,-0.02)
ax3.set_ylabel(r'$J$',rotation='horizontal')
ax3.set_yticks([-6e-5,-4e-5,-2e-5,0],labels=[r'$-6$',r'$-4$',r'$-2$',r'$0$'])
ax3.text(-0.035,-6.8e-5,r'$\times 10^{-5}$')
ax3.yaxis.set_label_coords(-0.05,0.93)

# adjust tick params
for ax in axs:
    ax.tick_params(which='both',direction='in',pad=1.5,length=2.5)
ax0.tick_params(axis='x',which='both',length=0)

# add (a), (b), etc.
edgewidth = 0.5
pad = 1.25
fontsize=8.0
labels=['(a)','(b)','(c)','(d)']
for i,ax in enumerate([ax0,ax1,ax2,ax3]):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    pad_x = pad / 72 / width
    pad_y = pad / 72 / height
    inside_pad = pad/fontsize
    print(pad_x,pad_y)
    text = ax.text(pad_x,1-pad_y,labels[i],ha='left',va='top',transform=ax.transAxes,
                       bbox=dict(boxstyle='square,pad='+str(inside_pad),edgecolor='black',facecolor='white',lw=edgewidth))
    text.set_clip_on(True)

plt.subplots_adjust(top=0.975,bottom=0.08,right=0.985,left=0.05)

# save figure
plt.savefig('disc_current_3panel.pdf',dpi=1000)