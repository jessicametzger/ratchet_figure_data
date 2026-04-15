import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
from types import SimpleNamespace
import os,sys
import re

plt_rcparams = {'text.usetex' : True,
                'font.size' : 8,
                'font.family' : 'serif',
                'text.latex.preamble' : r"\usepackage{lmodern} \usepackage{amstext}",
                'figure.figsize' : [3.4,3.4*0.7],
                'figure.dpi': 200
                }
plt.rcParams.update(plt_rcparams)

def isnum(s):
    try:
        a=float(s)
    except ValueError:
        return False
    return True

# function to read the parameter info of a simulation
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

# function to reconstruct cubic temperature landscape
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
    

#########################################################################################################################################
# IMPORT DATA

epr_series = np.genfromtxt('data/eprseries_avg',delimiter='\t')
epr = np.genfromtxt('data/eprprof_avg',delimiter='\t')
sig = np.genfromtxt('data/sigprof_avg',delimiter='\t')

p = get_params('data/1337-param')
dx_prof = p.Lx/p.Nbinx 
dy_prof = p.Ly/p.Nbiny
xs_prof = np.linspace(0,p.Lx - dx_prof, p.Nbinx)
ys_prof = np.linspace(0,p.Lx - dy_prof, p.Nbiny)
NstoreProf = int(p.tf / p.StoreInterProf + 1e-5)

# temperature landscape
Ts = f_cub(xs_prof, p.Txs, p.Ts, p.Lx)


#########################################################################################################################################
# PLOTTING

# color scheme
cm = plt.get_cmap('viridis')
c_stot = cm(0)
c_sid = cm(0.35)
c_sIK = cm(0.75)

fig = plt.figure(figsize=(3.4, 3.4*.5))
gs = gridspec.GridSpec(3, 2, height_ratios=[.45, .85, 1], width_ratios = [1,0.04], hspace=0, wspace=0.1)
ax2 = fig.add_subplot(gs[2, 0])
ax1 = fig.add_subplot(gs[1, 0],sharex=ax2)
ax0 = fig.add_subplot(gs[0, 0],sharex=ax2)
axs = [ax0,ax1,ax2]

# colormap for epr density
rng = np.max(np.abs(epr))
colors1 = plt.cm.ocean(np.linspace(0., 1, 128))
colors2 = plt.cm.hot(np.flip(np.linspace(0, 1, 128)))
colors_ = np.vstack((colors1, colors2))
cmap = colors.LinearSegmentedColormap.from_list('my_colormap', colors_)
epr_plot = axs[1].imshow(epr.T, extent=(0,p.Lx,0,p.Ly), origin='lower', cmap=cmap, 
                         vmin=-rng, vmax=rng,aspect='auto') 
cbar_ax = fig.add_subplot(gs[1, 1])
cbar_ax_inset = inset_axes(cbar_ax, width="100%", height="90%", 
                           bbox_to_anchor=(-0.45, -0.1, 1, 1),  # (x0, y0, width, height) in axes coordinates
                           bbox_transform=cbar_ax.transAxes,    # Transform the bbox to the axes coordinate system
                           borderpad=0)
cbar = fig.colorbar(epr_plot, cax=cbar_ax_inset, pad=0.1)
cbar.set_label(r'$\hat{s}(\mathbf{r})$',rotation='horizontal')
cbar.ax.yaxis.set_label_coords(1,1.25)
cbar.ax.set_yticks([-0.005,0,0.005],labels=['-0.005',r'$0$',r'$0.005$'])
cbar.ax.tick_params(axis='y', which='both', length=2, pad=1)
cbar.ax.set_ylim(ymin=-0.0052)

for spine in cbar_ax.spines.values():
    spine.set_visible(False)
for spine in cbar_ax_inset.spines.values():
    spine.set_visible(False)

cbar_ax.set_xticks([])
cbar_ax.set_xticklabels([])
cbar_ax.set_yticks([])
cbar_ax.set_yticklabels([])
    
# Ensure colorbar outline is visible
cbar_ax.spines['top'].set_visible(False)
cbar_ax.spines['right'].set_visible(False)
cbar_ax.spines['bottom'].set_visible(False)
cbar_ax.spines['left'].set_visible(False)
cbar_ax_inset.spines['top'].set_visible(True)
cbar_ax_inset.spines['right'].set_visible(True)
cbar_ax_inset.spines['bottom'].set_visible(True)
cbar_ax_inset.spines['left'].set_visible(True)

# stress plot
axs[2].plot(sig[:,0]+1/2., sig[:,1]+sig[:,4], color=c_stot,ls='--',label=r'$-\mathbf{\sigma}^{xx}_{\textnormal{\tiny tot}}$')
axs[2].plot(sig[:,0]+1/2., sig[:,4], color=c_sid,label=r'$\rho T$')
axs[2].plot(sig[:,0]+1/2., sig[:,1], color=c_sIK, label=r'$-\mathbf{\sigma}^{xx}_{\textnormal{\tiny IK}}$')

# entropy production time series inset
ax_inset = inset_axes(axs[1], width="18%", height="35%",
                      bbox_to_anchor=(0.65, -0.6, 1.44, 1.6),
                      bbox_transform=axs[1].transAxes, loc='upper left')
ax_inset.plot(np.linspace(0,p.tf,NstoreProf+1), [0]+list(epr_series), color='black')
ax_inset.text(.01*p.tf,np.max(epr_series*.7),r'$S(t)$')
ax_inset.text(0.9*p.tf,np.max(epr_series*0.05),r'$t$')
ax_inset.tick_params(axis='x', which='both', pad=1,length=2)
ax_inset.tick_params(axis='y', which='both', pad=1,length=2)
ax_inset.set_xticks([0,1e5],labels=['0',r'$10^5$'])
ax_inset.set_xlim(xmin=0,xmax=p.tf)
ax_inset.set_ylim(ymin=0)
ax_inset.xaxis.set_label_coords(0.5,-0.08)
ax_inset.set_yticks([0,2e6],labels=['0',r'$2\hspace{-0.2em}\times\hspace{-0.2em}10^6$'])

fig.add_subplot(gs[0, 1]).axis('off')

# T(x) plot
axs[0].plot(xs_prof+dx_prof/2.0, Ts, color='#6aa84fff', lw=2)

# axis limits
axs[2].set_xlim(xmin=0,xmax=p.Lx)
axs[0].set_ylim(ymax=17,ymin=0)
axs[2].set_ylim(ymin=0)
ylim = axs[2].get_ylim()
axs[2].set_ylim(ymax=ylim[1]*1.15)

# legend for stress plot
axs[2].legend(loc='lower right', bbox_to_anchor=(1.26,-0.1),handletextpad=0.5,
              frameon=False,handlelength=1.5)

# axis tick params
axs[0].set_yticks([0,10],labels=['0','10'])
axs[0].set_xticks([])
axs[0].set_xticklabels([])
axs[0].tick_params(axis='x', which='both', length=0,direction='in')
axs[0].tick_params(axis='y', length=2, pad=2,direction='in')
axs[1].tick_params(axis='x', which='both', length=0,direction='in')
axs[2].set_xticks([0,10,20,30,40],labels=[0,10,20,30,40])
axs[2].tick_params(axis='x', which='both', length=3, pad=2,direction='in')
axs[2].tick_params(axis='y', which='both', length=3, pad=2,direction='in')
axs[1].tick_params(axis='y', which='both', length=3, pad=2,direction='in')
axs[1].set_yticks([0,5,10],labels=[0,5,10])
axs[2].set_yticks([0,0.006])

# axis labels
axs[2].set_xlabel(r'$x$')
axs[1].set_ylabel(r'$y$', rotation='horizontal')
axs[0].set_ylabel(r'$T(x)$',rotation='horizontal')

# position axis labels
axs[2].xaxis.set_label_coords(0.97,-0.03)
axs[0].yaxis.set_label_coords(-0.07,0.75)
axs[1].yaxis.set_label_coords(-0.1,0.45)

plt.setp(axs[0].get_xticklabels(), visible=False)
plt.setp(axs[1].get_xticklabels(), visible=False)

# plot geometry
gs.update(hspace=0)
plt.subplots_adjust(left=0.091,bottom=0.07,right=0.9,top=0.98)

plt.savefig('EPR_tot.pdf',dpi=1000)
plt.savefig('EPR_tot.png',dpi=1000)


