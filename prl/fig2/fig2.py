import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
import os,sys

sys.path.append(os.getcwd())
import read_data as rd
import mydefaults as md

plt.rcParams.update(md.plt_rcparams)


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
    

epr_file = 'test-EPR'
sig_file = 'test-sigmaIKprof'
T_file = 'test-Tprof'
prof_file = 'test-prof'
param_file = 'test-param'

p = rd.get_params(param_file)
dx_prof = p.Lx/p.Nbinx 
dy_prof = p.Ly/p.Nbiny
xs_prof = np.linspace(0,p.Lx - dx_prof, p.Nbinx)
ys_prof = np.linspace(0,p.Lx - dy_prof, p.Nbiny)
NstoreProf = md.fint(p.tf / p.StoreInterProf)

# temperature landscape
Ts = f_cub(xs_prof, p.Txs, p.Ts, p.Lx)

#########################################################################################################################################
# EPR DENSITY

epr = np.genfromtxt(epr_file, delimiter='\t')

# avg over time (not x or y)
# turn into epr PER UNIT VOL, per particle
T_meas = NstoreProf*p.dt*p.NstepProf # effective amount of time over which EPR is measured
unit_vol = dx_prof*dy_prof
epr_tot_img = sum([epr[i*p.Nbinx:(i+1)*p.Nbinx,1:].T for i in range(NstoreProf)])/(T_meas*p.N*unit_vol)

plt.hist(epr_tot_img.flatten())
plt.show()

# sum over x and y, then sum cumulatively over time
series_EPR_tot = np.array([np.sum(epr[i*p.Nbinx:(i+1)*p.Nbinx,1:]) for i in range(NstoreProf)])
series_EPR_tot = np.cumsum(series_EPR_tot)

#########################################################################################################################################



#########################################################################################################################################
# IRVING-KIRKWOOD STRESS

sig = np.genfromtxt(sig_file, delimiter='\t')[:,3:] / (p.NstepProf*p.N*p.dt)
xs_sig = np.linspace(0,p.Lx-1,md.fint(p.Lx))
ys_sig = np.linspace(0,p.Ly-1,md.fint(p.Ly))
sig_tot = sum([sig[i*md.fint(p.Lx*3):(i+1)*md.fint(p.Lx*3),:] for i in range(NstoreProf)]) / NstoreProf
sig_tot_vs_x = np.zeros((md.fint(p.Lx),3))
sig_tot_vs_x[:,0] = np.sum(sig_tot[:md.fint(p.Lx)],axis=1) / p.Ly
sig_tot_vs_x[:,1] = np.sum(sig_tot[md.fint(p.Lx):md.fint(p.Lx*2)],axis=1) / p.Ly
sig_tot_vs_x[:,2] = np.sum(sig_tot[md.fint(p.Lx*2):md.fint(p.Lx*3)],axis=1) / p.Ly

#########################################################################################################################################



#########################################################################################################################################
# THERMAL STRESS

Tprof = np.genfromtxt(T_file, delimiter='\t')[:,3:]
Tprof = Tprof / (p.NstepProf*p.N*dx_prof)
Tprof = sum([Tprof[i*p.Nbinx:(i+1)*p.Nbinx] for i in range(NstoreProf)]) / NstoreProf # avg over time
Tprof = np.sum(Tprof,axis=1)/p.Ly # sum over y

# make coarse version to sum with irving-kirkwood stress
binw = md.fint(1 / dx_prof)
Tprof_coarse = sum([Tprof[i::binw] for i in range(binw)]) / binw

#########################################################################################################################################



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

rng = np.max(np.abs(epr_tot_img))
colors1 = plt.cm.ocean(np.linspace(0., 1, 128))
colors2 = plt.cm.hot(np.flip(np.linspace(0, 1, 128)))
colors_ = np.vstack((colors1, colors2))
cmap = colors.LinearSegmentedColormap.from_list('my_colormap', colors_)
epr_plot = axs[1].imshow(epr_tot_img, extent=(0,p.Lx,0,p.Ly), origin='lower', cmap=cmap, 
                         vmin=-rng, vmax=rng,aspect='auto') 
cbar_ax = fig.add_subplot(gs[1, 1])
cbar_ax_inset = inset_axes(cbar_ax, width="100%", height="90%", 
                           bbox_to_anchor=(-0.45, -0.1, 1, 1),  # (x0, y0, width, height) in axes coordinates
                           bbox_transform=cbar_ax.transAxes,    # Transform the bbox to the axes coordinate system
                           borderpad=0)
cbar = fig.colorbar(epr_plot, cax=cbar_ax_inset, pad=0.1)
cbar.set_label(r'$\hat{s}(\mathbf{r})$',rotation='horizontal')
cbar.ax.yaxis.set_label_coords(1,1.25)
cbar.ax.set_yticks([-1,0,1],labels=['-1',r'$0$',r'$1$'])
cbar.ax.tick_params(axis='y', which='both', length=2, pad=1)
cbar.ax.set_ylim(ymin=-1,ymax=1)

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


axs[2].plot(xs_sig+1/2., sig_tot_vs_x[:,0]+Tprof_coarse, color=c_stot,ls='--',label=r'$-\mathbf{\sigma}^{xx}_{\textnormal{\tiny tot}}$')
axs[2].plot(xs_prof+dx_prof/2., Tprof, color=c_sid,label=r'$\rho T$')
axs[2].plot(xs_sig+1/2., sig_tot_vs_x[:,0], color=c_sIK, label=r'$-\mathbf{\sigma}^{xx}_{\textnormal{\tiny IK}}$')

ax_inset = inset_axes(axs[1], width="18%", height="35%",
                      bbox_to_anchor=(0.65, -0.6, 1.44, 1.6),
                      bbox_transform=axs[1].transAxes, loc='upper left')
ax_inset.plot(np.linspace(0,p.tf,NstoreProf+1), [0]+list(series_EPR_tot), color='black')
ax_inset.text(.01*p.tf,np.max(series_EPR_tot*.7),r'$S(t)$')
ax_inset.text(0.9*p.tf,np.max(series_EPR_tot*0.05),r'$t$')
ax_inset.tick_params(axis='x', which='both', pad=1,length=2)
ax_inset.tick_params(axis='y', which='both', pad=1,length=2)
ax_inset.set_xticks([0,1e5],labels=['0',r'$10^5$'])
ax_inset.set_xlim(xmin=0,xmax=p.tf)
ax_inset.set_ylim(ymin=0)
ax_inset.xaxis.set_label_coords(0.5,-0.08)
ax_inset.set_yticks([0,2e5],labels=['0',r'$2\hspace{-0.2em}\times\hspace{-0.2em}10^5$'])

fig.add_subplot(gs[0, 1]).axis('off')

xs = np.linspace(0,p.Lx,400)
ys = np.linspace(0,p.Ly,400)
dx_th = p.Lx/400.
axs[0].plot(xs_prof+dx_prof/2.0, Ts, color='#6aa84fff', lw=2)

axs[2].set_xlim(xmin=0,xmax=p.Lx)
axs[0].set_ylim(ymax=17,ymin=0)
axs[2].set_ylim(ymin=0)
ylim = axs[2].get_ylim()
axs[2].set_ylim(ymax=ylim[1]*1.15)

axs[2].legend(loc='lower right', bbox_to_anchor=(1.26,-0.1),handletextpad=0.5,
              frameon=False,handlelength=1.5)

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

axs[2].set_xlabel(r'$x$')
axs[1].set_ylabel(r'$y$', rotation='horizontal')
axs[0].set_ylabel(r'$T(x)$',rotation='horizontal')

axs[2].xaxis.set_label_coords(0.97,-0.03)
axs[0].yaxis.set_label_coords(-0.07,0.75)
axs[1].yaxis.set_label_coords(-0.1,0.45)

plt.setp(axs[0].get_xticklabels(), visible=False)
plt.setp(axs[1].get_xticklabels(), visible=False)

gs.update(hspace=0)  # Increase this value to add more space
plt.subplots_adjust(left=0.091,bottom=0.07,right=0.9,top=0.98)

plt.savefig(param_file.split('-')[0]+'_EPR_tot.pdf',dpi=1000)
plt.savefig(param_file.split('-')[0]+'_EPR_tot.png',dpi=1000)


