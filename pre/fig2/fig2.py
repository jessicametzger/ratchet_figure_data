'''
Create the figure, using the averaged data from data/get_data.py
'''

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

folders=['data/int','data/nonint']

# interacting data
disp_avg = np.genfromtxt(os.path.join(folders[0],'disp_avg'),delimiter='\t')
disp_err = np.genfromtxt(os.path.join(folders[0],'disp_err'),delimiter='\t')
prof_avg = np.genfromtxt(os.path.join(folders[0],'prof_avg'),delimiter='\t')
sigm_avg = np.genfromtxt(os.path.join(folders[0],'sigm_avg'),delimiter='\t')
sigIK_avg = np.genfromtxt(os.path.join(folders[0],'sigIK_avg'),delimiter='\t')

# non-interacting data
disp_NI_avg = np.genfromtxt(os.path.join(folders[1],'disp_avg'),delimiter='\t')
disp_NI_err = np.genfromtxt(os.path.join(folders[1],'disp_err'),delimiter='\t')
prof_NI_avg = np.genfromtxt(os.path.join(folders[1],'prof_avg'),delimiter='\t')

# get param info
p = get_params([os.path.join(folders[0],x) for x in os.listdir(folders[0]) if x.endswith('-param')][0])
dx = p.Lx / p.Nbinx
dy = p.Ly / p.Nbiny

# get x and y axis data
dx = p.Lx/p.Nbinx
xs_data = np.arange(0,p.Lx-1e-5, dx)
xs_sig = np.linspace(0,p.Lx-1,int(p.Lx))
dx_sig = 1.0

# bin thermal stress array to match irving-kirkwood stress (coarser grid)
binw = int(dx_sig / dx + 1e-5)
sigm_avg = sum([sigm_avg[i::binw,:] for i in range(binw)])/binw

# t axis data
ts = np.arange(0,p.tf+p.StoreInterDisp/2.0,p.StoreInterDisp)


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

xs = np.linspace(0,p.Lx-dx,p.Nbinx) + dx/2.0
fofx = f_cub(xs, p.Txs, p.Ts, p.Lx)



##########################################################################################################################################
# PLOTTING

# colors
cm1=plt.get_cmap('plasma')
cm = plt.get_cmap('viridis')
c_stot = cm(0.1)
c_stoty= cm1(0.5)
c_stotxy= cm1(0.62)
c_stotyx= cm1(0.85)
c_sA1 = cm(0.35)
c_sIK = cm(0.75)
c_sA2 = cm(0.99)

# create figure and subplots
fig=plt.figure(figsize=(7,3.4*0.55),dpi=200)
gs = fig.add_gridspec(100, 3, hspace=0, wspace=0.08, height_ratios=[1]*100,width_ratios=[.9,.6,1])
ax4 = fig.add_subplot(gs[71:100,2])
ax3 = fig.add_subplot(gs[0:65,2],sharex=ax4)
ax2 = fig.add_subplot(gs[:,1])
ax1 = fig.add_subplot(gs[33:,0])
ax0 = fig.add_subplot(gs[0:27,0],sharex=ax1)
axs = [ax0,ax1,ax2,ax3,ax4]

# plot temperature landscape
ax0.plot(xs_data,fofx,color='green',lw=2.)

# plot density landscape & overdamped theory
ax1.plot(xs_data+dx/2.0, (1/fofx)/np.trapezoid(1/fofx,xs_data), color='tan',lw=1.,label=r'$\propto 1/T(x)$')
ax1.plot(xs_data+dx/2.0, prof_NI_avg, color='darkgoldenrod',lw=2., label=r'$\rho_{\rm s}(x),\,U_{\rm int}=0$')
ax1.plot(xs_data+dx/2.0, prof_avg, color='brown', lw=2., label=r'$\rho_{\rm s}(x),\,U_{\rm int}\neq 0$')

# plot xx stress components
ax3.plot(xs_sig+dx_sig/2.0, sigIK_avg[:,0]+sigm_avg[:,0],lw=2., color=c_stot,
         label=r'$-\sigma_{\rm tot}^{xx}$')
ax3.plot([-1,-1],[0,0],color='white',alpha=0,label=' ')
ax3.plot(xs_sig+dx_sig/2.0, sigIK_avg[:,0], color=c_sIK, lw=2., label=r'$-\sigma_{\rm IK}^{xx}$')
ax3.plot([-1,-1],[0,0],color='white',alpha=0,label=' ')
ax3.plot(xs_sig+dx_sig/2.0, sigm_avg[:,0], color=c_sA1, lw=2., label=r'$\langle (v^x)^2\rangle$')
ax3.plot([-1,p.Lx+1],[0,0],lw=0.3,color='black')

# plot total xx and yy stresses
ax4.plot(xs_sig+dx_sig/2.0, sigIK_avg[:,0]+sigm_avg[:,0],lw=2., color=c_stot,
         label=r'$-\sigma_{\rm tot}^{xx}$')
ax4.plot(xs_sig+dx_sig/2.0, sigIK_avg[:,2]+sigm_avg[:,2],lw=2., color=c_stoty, 
         label=r'$-\sigma_{\rm tot}^{yy}$')

# plot particle displacement over time (ends up being zero on average)
ax2.plot([0,p.tf],[0,0],color='black',lw=.75)
ax2.fill_between(ts,(disp_NI_avg-3*disp_NI_err)/p.Lx, (disp_NI_avg+3*disp_NI_err)/p.Lx,
                color='darkgoldenrod',alpha=0.2,linewidth=0.5,zorder=0)
ax2.plot(ts,disp_NI_avg/p.Lx,color='darkgoldenrod',lw=1,marker='o',
         markersize=4,zorder=3,label=r'$U_{\rm int}=0$',markerfacecolor='white')
ax2.fill_between(ts,(disp_avg-3*disp_err)/p.Lx, (disp_avg+3*disp_err)/p.Lx,
                color='darkred',alpha=0.2,linewidth=0.5,zorder=0)
ax2.plot(ts,disp_avg/p.Lx,color='darkred',lw=1,marker='o',markersize=4,zorder=3,label=r'$U_{\rm int}\neq 0$')

for ax in axs[:2]:
    ax.set_ylim(ymin=0)

# adjust ticks
for ax in axs:
    ax.tick_params(axis='both',which='both',direction='in',length=2,pad=2)
for ax in [ax0,ax3,ax4]:
    ax.tick_params(axis='x', which='both', labelbottom=False)
ax2.tick_params(axis='y',which='both',pad=-22)

# positioning y ticks
for ax in axs[3:]:
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")

# x axis limits
for ax in axs:
    ax.set_xlim(xmin=0, xmax=p.Lx)
ax2.set_xlim(xmin=0,xmax=p.tf)
    
# y axis limits
ax0.set_ylim(ymax=35)
ax1.set_ylim(ymax=0.06)
ax2.set_ylim(ymax=1.45,ymin=-1.1)
ax3.set_ylim(ymax=0.24,ymin=0)
ax4.set_ylim(0.16,0.24)

# axis ticks & labels
ax0.set_yticks([0,10,20])
ax1.set_yticks([0,0.02,0.04])
ax1.set_xticks([0,10,20,30,40,50],labels=['0','','20','','40',''])
ax2.set_xticks([0,5e3,1e4,1.5e4],labels=['0','',r'$10^4$',''])
ax2.set_yticks([-1,-.5,0,.5,1])

# legends
ax1.legend(frameon=False,loc='upper right',bbox_to_anchor=(1,1),handletextpad=0.5, handlelength=1.75)
ax2.legend(frameon=False,loc='upper right',bbox_to_anchor=(1.04,.93),handletextpad=0.5, handlelength=1.75)
ax3.legend(frameon=False,ncol=3,loc='upper right',bbox_to_anchor=(1.02,1.06), columnspacing=1.5,handletextpad=0.5, handlelength=1.75)
ax4.legend(frameon=False,ncol=2,loc='upper right',bbox_to_anchor=(1,1.1), columnspacing=1.5,handletextpad=0.5, handlelength=1.75)

# axis labels
ax0.set_ylabel(r'$T(x)$',rotation='horizontal')
ax0.yaxis.set_label_coords(0.9,0.2)
ax1.set_xlabel(r'$x$')
ax1.xaxis.set_label_coords(0.97,-.03)
ax2.set_ylabel(r'$\Delta x(t)/L_x$',rotation='horizontal')
ax2.yaxis.set_label_coords(0.6,0.9)
ax2.set_xlabel(r'$t$')
ax2.xaxis.set_label_coords(0.9,-.02)

# add (a), (b), etc.
edgewidth = 0.5
pad = 1.25
fontsize=8.0
labels=['(a)','(b)','(c)','(d)','(e)','(f)']
for i,ax in enumerate([ax0,ax1,ax2,ax3,ax4]):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    pad_x = pad / 72 / width
    pad_y = pad / 72 / height
    inside_pad = pad/fontsize
    text = ax.text(pad_x,1-pad_y,labels[i],ha='left',va='top',transform=ax.transAxes,
                       bbox=dict(boxstyle='square,pad='+str(inside_pad),edgecolor='black',facecolor='white',lw=edgewidth))
    text.set_clip_on(True)

# figure geometry
plt.subplots_adjust(left=0.035,bottom=0.07,top=0.99,right=0.955)

# save figure
plt.savefig('ud_pbp_eos_large.pdf',dpi=1000)