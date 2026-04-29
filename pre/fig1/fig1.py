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

# list of folders: 4 different types of system, interacting + non-interacting
folders = ['UPBPs','PBPs','ABPs','AOUPs']
folders = list(np.concatenate([[f+'/int',f+'/nonint'] for f in folders]))
subdir1 = 'tf_8000000-dt_0.000250-N_120-Lx_30-Ly_10-Nmode_7'
subdir2 = 'tf_8000000-dt_0.000250-N_120-Lx_30-Ly_10-Nmode_7-k_10'
fs = [os.path.join('data/'+f,subdir1) if f.endswith('nonint') else os.path.join('data/'+f,subdir2) for f in folders]

# collect data
disp_avgs = []
disp_errs = []
profs = []
Jxprofs = []
Jyprofs = []
for f in fs:
	disp_avgs += [np.genfromtxt(os.path.join(f,'disp_avg'),delimiter='\t')]
	disp_errs += [np.genfromtxt(os.path.join(f,'disp_err'),delimiter='\t')]
	profs += [np.genfromtxt(os.path.join(f,'prof_avg'),delimiter='\t')]
	Jxprofs += [np.genfromtxt(os.path.join(f,'Jxprof_avg'),delimiter='\t')]
	Jyprofs += [np.genfromtxt(os.path.join(f,'Jyprof_avg'),delimiter='\t')]

# get param data for one
p = get_params([os.path.join(f,x) for x in os.listdir(f) if x.endswith('-param')][0])

# get bounds for plots
vmax = max([np.max(np.abs(x-1)) for x in profs])
disp_max = max([np.max(disp_avgs[i]+3*disp_errs[i]) for i in range(len(disp_avgs))])
disp_min = min([np.min(disp_avgs[i]-3*disp_errs[i]) for i in range(len(disp_avgs))])

# get x and y axis data
dx = p.Lx/p.Nbinx
dy = p.Ly/p.Nbiny
xgs = np.linspace(0,p.Lx-dx,p.Nbinx)
ygs = np.linspace(0,p.Ly-dy,p.Nbiny)
xgs, ygs = np.meshgrid(xgs,ygs)

# t axis data
ts = np.arange(0,p.tf+p.StoreInterDisp/2.0,p.StoreInterDisp)


##########################################################################################################################################
# ACTIVITY LANDSCAPE

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
fofx = f_cub(xs, p.Dxs, p.Ds, p.Lx)


##########################################################################################################################################
# PLOTTING

# plot geometry
gapL = 5  # gap to left of plots
gapR = 3  # gap to right of plots
gapT = 3  # gap above plots
gapB = 4  # gap below plots
h0 = 8    # height of v(x), T(x), etc. plots
h1 = 10   # height of heatmap plots
h2 = 10   # height of displacement plots
wsp = 3   # horizontal spacing width
hsp0 = 2  # vertical spacing width below v(x), T(x), etc. plots
hsp1 = 2  # vertical spacing between heatmap plots
hsp2 = 5  # vertical spacing between heatmap plots and displacement plots
caxw = 2  # width of colorbar axis
caxh = 14 # height of colorbar axis

# auxiliary geometric parameters
w1 = h1*3
w2 = h1*2-2
caxtop = int((h0+h1*2+hsp1*2-caxh)/2.0+1e-5)
htot = gapT+gapB+h0+h1*2+h2+hsp1*2+hsp2
wtot = gapL+gapR+w1*4+wsp*4+caxw
aspect = htot/wtot # overall aspect ratio

tickpad = 1.2
ticklen = 2

# colors
cmrho = plt.get_cmap('Spectral_r')
cdisps = ['slateblue','palevioletred']
cJ = 'black'

# create figure and subplots
fig = plt.figure(figsize=(7,7*aspect))
gs = fig.add_gridspec(h0+2*h1+h2+hsp1*2+hsp2,4*w1+wsp*4+caxw,wspace=0,hspace=0)
axs = [[fig.add_subplot(gs[:h0,:w1]), fig.add_subplot(gs[:h0,w1+wsp:2*w1+wsp]), fig.add_subplot(gs[:h0,w1*2+wsp*2:w1*3+wsp*2]), fig.add_subplot(gs[:h0,w1*3+wsp*3:w1*4+wsp*3])],
       [fig.add_subplot(gs[h0+hsp1:h0+h1+hsp1,:w1]), fig.add_subplot(gs[h0+hsp1:h0+h1+hsp1,w1+wsp:2*w1+wsp]), fig.add_subplot(gs[h0+hsp1:h0+h1+hsp1,w1*2+wsp*2:w1*3+wsp*2]), fig.add_subplot(gs[h0+hsp1:h0+h1+hsp1,w1*3+wsp*3:w1*4+wsp*3])],
       [fig.add_subplot(gs[h0+h1+hsp1*2:h0+h1*2+hsp1*2,:w1]), fig.add_subplot(gs[h0+h1+hsp1*2:h0+h1*2+hsp1*2,w1+wsp:2*w1+wsp]), fig.add_subplot(gs[h0+h1+hsp1*2:h0+h1*2+hsp1*2,w1*2+wsp*2:w1*3+wsp*2]), fig.add_subplot(gs[h0+h1+hsp1*2:h0+h1*2+hsp1*2,w1*3+wsp*3:w1*4+wsp*3])],
       [fig.add_subplot(gs[h0+h1*2+hsp1*2+hsp2:,2:w1-2]), fig.add_subplot(gs[h0+h1*2+hsp1*2+hsp2:,w1+wsp+2:2*w1+wsp-2]), fig.add_subplot(gs[h0+h1*2+hsp1*2+hsp2:,w1*2+wsp*2+2:w1*3+wsp*2-2]), fig.add_subplot(gs[h0+h1*2+hsp1*2+hsp2:,w1*3+wsp*3+2:w1*4+wsp*3-2])]]

# first, do axis titles, labels, limits, etc.
titles=['UBPs','OBPs','ABPs','AOUPs']
for j in range(4):

	# system label, fluctuation landscape plot limits
    axs[0][j].set_title(titles[j],pad=3)
    axs[0][j].set_ylim(0,5)

    # heatmap x axis
    axs[2][j].set_xlabel(r'$x$')
    axs[2][j].xaxis.set_label_coords(0.5,-0.15)

    # displacement plot
    axs[3][j].set_xlim(0,p.tf)
    axs[3][j].set_ylim(disp_min*1.5,disp_max*1.2)
    axs[3][j].set_xticks([0,3e6,6e6],labels=['0',r'$3\times 10^6$',r'$6\times 10^6$'])
    axs[3][j].set_xlabel(r'$t$')
    axs[3][j].xaxis.set_label_coords(0.5,-0.2)
    axs[3][j].set_yticks([0,2500])

    # x-axis of heatmap/fluctuation landscape plots
    for i in range(3):
        axs[i][j].set_xlim(0,p.Lx)
        axs[i][j].set_xticks([0,p.Lx/3,2*p.Lx/3,p.Lx],labels=None if i==2 else ['']*4)

    # y-axis of heatmap plots
    for i in range(1,3):
        axs[i][j].set_ylim(0,p.Ly)
        axs[i][j].set_yticks([0,p.Ly/2,p.Ly],labels=None if j==0 else ['']*3)

    # axis ticks
    for i in range(4):
        axs[i][j].tick_params(which='both',direction='in',length=ticklen,pad=tickpad)

# y axis label for leftmost plots
for i in range(1,3):
    axs[i][0].set_ylabel(r'$y$',rotation='horizontal')
    axs[i][0].yaxis.set_label_coords(-0.135,0.5)
axs[3][0].set_ylabel(r'$\Delta x(t)$',rotation='horizontal')
axs[3][0].yaxis.set_label_coords(-0.15,0.7)

# adjust figure geometry
plt.subplots_adjust(left=gapL/wtot,right=1-gapR/wtot,bottom=gapB/htot,top=1-gapT/htot)

# figure out how to scale arrows --- OD PBPs should be much smaller; choose empirical scale to match UD
scale=0.12
binw = 5
scale1 = scale * np.sqrt(np.mean(bin_arr(Jxprofs[2],binw).T**2+bin_arr(Jyprofs[2],binw).T**2 + bin_arr(Jxprofs[3],binw).T**2 + bin_arr(Jyprofs[3],binw).T**2)) / np.sqrt(np.mean(bin_arr(Jxprofs[0],binw).T**2+bin_arr(Jyprofs[0],binw).T**2 + bin_arr(Jxprofs[1],binw).T**2 + bin_arr(Jyprofs[1],binw).T**2))

# for (a), (b), etc. tags
ew = 0.5
pad = 1.25
extra=2
fontsize = 8.0
tags=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)']

# plot data
labels=['Interact.','Non-interact.']
fofx_labels = [r'$T(x)$', r'$T(x)$', r'$v(x)$', r'$D(x)$']
count = 0
for j in range(4):
	# fluctuation landscape display
    axs[0][j].plot(xs, fofx, color='#6aa84fff')
    
    # heatmaps (interact & non-interact)
    axs[1][j].imshow(profs[j*2].T,origin='lower',extent=(0,p.Lx,0,p.Ly),cmap=cmrho,vmax=1.5,vmin=0.5)
    axs[2][j].imshow(profs[j*2+1].T,origin='lower',extent=(0,p.Lx,0,p.Ly),cmap=cmrho,vmax=1.5,vmin=0.5)

    # current arrows
    binw = 10
    axs[1][j].quiver(bin_arr(xgs,binw), bin_arr(ygs,binw), bin_arr(Jxprofs[j*2],binw).T,bin_arr(Jyprofs[j*2],binw).T,color=cJ,scale=scale if j!=1 else scale1)
    axs[2][j].quiver(bin_arr(xgs,binw), bin_arr(ygs,binw), bin_arr(Jxprofs[j*2+1],binw).T,bin_arr(Jyprofs[j*2+1],binw).T,color=cJ,scale=scale if j!=1 else scale1)

    # displacement data (with error bars)
    axs[3][j].fill_between(ts, disp_avgs[j*2]-disp_errs[j*2]*3, disp_avgs[j*2]+disp_errs[j*2]*3, color=cdisps[0], lw=0,alpha=0.3)
    axs[3][j].fill_between(ts, disp_avgs[j*2+1]-disp_errs[j*2+1]*3, disp_avgs[j*2+1]+disp_errs[j*2+1]*3, color=cdisps[1], lw=0,alpha=0.3)
    axs[3][j].plot(ts, disp_avgs[j*2], color=cdisps[0], label='Interact.')
    axs[3][j].plot(ts, disp_avgs[j*2+1], color=cdisps[1], label='Non-interact.')

    if j==0:
        axs[3][j].legend(frameon=False,handlelength=1,handletextpad=0.5,loc='upper right')

    # fluctuation landscape labels
    bbox = axs[0][j].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width,height = bbox.width,bbox.height
    pad_x = (pad+extra) / 72 / width
    pad_y = (pad+extra) / 72 / height
    inside_pad = pad/fontsize
    text = axs[0][j].text(1-pad_x,1-pad_y,fofx_labels[j],ha='right',va='top',transform=axs[0][j].transAxes,color='#6aa84fff')

    # system labels (interact or non-interact)
    for i in range(2):
        bbox = axs[i][j].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        pad_x = (pad+extra) / 72 / width
        pad_y = (pad+extra) / 72 / height
        inside_pad = pad/fontsize
        text = axs[i+1][j].text(pad_x,pad_y,labels[i],ha='left',va='bottom',transform=axs[i+1][j].transAxes,
                              bbox=dict(boxstyle='square,pad='+str(inside_pad),edgecolor=cdisps[i],facecolor=cdisps[i],lw=ew,alpha=0.7))

    # (a), (b), etc. labels
    for i in range(4):
        bbox = axs[i][j].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        pad_x = pad / 72 / width
        pad_y = pad / 72 / height
        inside_pad = pad/fontsize
        text = axs[i][j].text(pad_x,1-pad_y,tags[count],ha='left',va='top',transform=axs[i][j].transAxes,
                              bbox=dict(boxstyle='square,pad='+str(inside_pad),edgecolor='black',facecolor='white',lw=ew))
        count += 1

# colormap 
cax = fig.add_subplot(gs[caxtop:caxtop+caxh,w1*4+wsp*4:w1*4+wsp*4+caxw])
cax.imshow(np.array([np.linspace(1-vmax,1+vmax,200)]).T,origin='lower',cmap=cmrho,extent=(0,1,0,2))
cax.yaxis.tick_right()
cax.set_aspect('auto')
cax.set_xticks([])
cax.set_title(r'$\rho/\overline{\rho}$')
cax.tick_params(which='both',axis='both',direction='in',pad=tickpad,length=ticklen)

# save figure
plt.savefig('ratchet_exceptions_summary_fofx.pdf',dpi=1000)
