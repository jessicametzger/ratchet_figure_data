import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import os,sys

sys.path.append(os.getcwd())
import read_data as rd
import mydefaults as md

plt.rcParams.update(md.plt_rcparams)


################################################################################################################
# AUXILIARY FUNCTIONS

# Bin 2d array 
def bin_arr(arr,binw):
    arr_=np.copy(arr)
    arr_ = sum([arr_[i::binw,:] for i in range(binw)])/binw
    arr_ = sum([arr_[:,i::binw] for i in range(binw)])/binw
    return arr_

# activity landscape function
def v(xs,p):
    cs = np.copy(p.vs)*0
    ds = np.copy(cs)
    for i in range(p.Nstep_v-1):
        dx = p.vxs[i+1]-p.vxs[i]
        if dx==0:
            dx=1
        cs[i] = -3*(p.vs[i] - p.vs[i+1])/dx**2
        ds[i] = 2*(p.vs[i]-p.vs[i+1])/dx**3
    dx = p.vxs[0]+p.Lx-p.vxs[-1]
    if dx==0:
        dx=1
    cs[-1] = -3*(p.vs[-1]-p.vs[0])/dx**2
    ds[-1] = 2*(p.vs[-1]-p.vs[0])/dx**3
    
    ans = np.copy(xs)*0
    inds = xs<p.vxs[0]
    dxs = xs[inds]+p.Lx-p.vxs[-1]
    ans[inds] = p.vs[-1] + cs[-1]*dxs**2 + ds[-1]*dxs**3
    for i in range(1,p.Nstep_v+1):
        if i==p.Nstep_v:
            upper = p.Lx+p.vxs[0]
        else:
            upper = p.vxs[i]
        inds = (xs>=p.vxs[i-1]) & (xs<upper)
        dxs = xs[inds] - p.vxs[i-1]
        ans[inds] = p.vs[i-1] + cs[i-1]*dxs**2 + ds[i-1]*dxs**3
    return ans


################################################################################################################
# READ DATA

p = rd.get_params('param')
dx = p.Lx/p.Nbinx
pd = p.__dict__
vxs = pd['vx1,vx2,...']
vs = pd['v1,v2,...']
p.vs = vs
p.vxs = vxs

Jprof = np.genfromtxt('Jprof_avg',delimiter='\t')
FAprof = np.genfromtxt('FAprof_avg',delimiter='\t')
Fintprof = np.genfromtxt('Fintprof_avg',delimiter='\t')
sigIKprof = np.genfromtxt('sigIKprof_avg',delimiter='\t')
sigAprof = np.genfromtxt('sigAprof_avg',delimiter='\t')
disp = np.genfromtxt('disp_avg',delimiter='\t')


################################################################################################################
# SHIFT, BIN, ETC. DATA

shiftx = 5
ishiftx = md.fint(shiftx/dx)
ishiftx1 = md.fint(shiftx)
Jprof = np.concatenate((Jprof[-ishiftx:,:],Jprof[:-ishiftx,:]),axis=0)
FAprof = np.concatenate((FAprof[-ishiftx:,:],FAprof[:-ishiftx,:]),axis=0)
Fintprof = np.concatenate((Fintprof[-ishiftx:,:],Fintprof[:-ishiftx,:]),axis=0)
sigIKprof = np.concatenate((sigIKprof[-ishiftx1:,:],sigIKprof[:-ishiftx1,:]),axis=0)
sigAprof = np.concatenate((sigAprof[-ishiftx:,:],sigAprof[:-ishiftx,:]),axis=0)

dsigAdxprof = np.copy(sigAprof)
dsigAdxprof[1:-1,:] = (sigAprof[2:,:]-sigAprof[:-2,:])/2.0/dx
dsigAdxprof[0] = (sigAprof[1,:]-sigAprof[-1,:])/2.0/dx
dsigAdxprof[-1] = (sigAprof[0,:]-sigAprof[-2,:])/2.0/dx

dsigIKdxprof = np.copy(sigIKprof)
dsigIKdxprof[1:-1,:] = (sigIKprof[2:,:]-sigIKprof[:-2,:])/2.0
dsigIKdxprof[0] = (sigIKprof[1,:]-sigIKprof[-1,:])/2.0
dsigIKdxprof[-1] = (sigIKprof[0,:]-sigIKprof[-2,:])/2.0

xs_ = np.linspace(0,p.Lx-dx,p.Nbinx)
xs1_ = np.linspace(0,p.Lx-1,md.fint(p.Lx))
xs = xs_+dx/2.0
xs1 = xs1_+0.5
dx1 = 1
binw = md.fint(1/dx)
print('bin size = ',binw)

vs = v(xs,p)
vs = np.concatenate((vs[-ishiftx:],vs[:-ishiftx]),axis=0)

# make copies on coarser grid (to match sigIK)
FAprof1 = sum([FAprof.copy()[i::binw,:] for i in range(binw)])
Fintprof1 = sum([Fintprof.copy()[i::binw,:] for i in range(binw)])
sigAprof1 = sum([sigAprof.copy()[i::binw,:] for i in range(binw)])
dsigAdxprof1 = sum([dsigAdxprof.copy()[i::binw,:] for i in range(binw)])

# calculate size of momentum sources
dp1 = np.sum(FAprof[(xs>=p.vxs[0]+shiftx) & (xs<=p.vxs[1]+shiftx),0])*dx
dp2 = np.sum(FAprof[(xs>=p.vxs[2]+shiftx) & (xs<=p.vxs[3]+shiftx),0])*dx

# calculate J due to stress
Js = FAprof[:,0]+Fintprof[:,0]-dsigAdxprof[:,0]

# stress J -> coarser grid
binw = 10
xsbin = sum([xs_[i::binw] for i in range(binw)])/binw + dx*binw/2.0
Jsbin = sum([Js[i::binw] for i in range(binw)])/binw

# observed J
Jobs = disp[-1]/p.tf/p.N


################################################################################################################
# PLOT

fig = plt.figure(figsize=(3.4,3.4*.7),dpi=200)
gs = fig.add_gridspec(4,2,hspace=0.15,height_ratios=[0.5,1,0.5,0.5],width_ratios=[1,0.4])
axs = [fig.add_subplot(gs[0,0]), fig.add_subplot(gs[1,0]), fig.add_subplot(gs[2,0]),fig.add_subplot(gs[3,0])]

# colors
cm = plt.get_cmap('viridis')
c_stot = cm(0)
c_sA = cm(0.35)
c_sIK = cm(0.75)
c_J = cm(0.84)
c_FA = cm(0.99)
c_intFA = cm(0.5)
c_Jobs = c_sA
c_v = '#6aa84fff'
lw=1.5
lbs = 0.5
htp = 0.7
hl = 1.5

# v(x) plot
axs[0].fill_between([vxs[0]+shiftx,vxs[1]+shiftx],[-1,-1],[500,500],color=c_intFA,alpha=0.3,linewidth=0)
axs[0].fill_between([vxs[2]+shiftx,vxs[3]+shiftx],[-1,-1],[500,500],color=c_intFA,alpha=0.3,linewidth=0)
axs[0].plot(xs,vs,color=c_v,label=r'$v(x)$')
axs[0].set_ylim(0,np.max(p.vs)*1.5)
axs[0].text((vxs[0]+vxs[1])/2.0+shiftx,np.max(p.vs)*1.4,r'$\mathcal{R}_1$',ha='center',va='top')
axs[0].text((vxs[2]+vxs[3])/2.0+shiftx,np.max(p.vs)*1.4,r'$\mathcal{R}_2$',ha='center',va='top')
axs[0].legend(frameon=False, bbox_to_anchor=(0.99, 0.62), loc='center left', handletextpad=htp,labelspacing=lbs,handlelength=hl)

# stress plot
axs[1].plot(xs1, sigIKprof[:,0]+sigAprof1[:,0],lw=lw, color=c_stot,  label=r'$-\sigma_{\textnormal{\tiny tot}}^{xx}$', zorder=4)
axs[1].plot(xs , sigAprof[:,0]/dx,             lw=lw, color=c_sA,    label=r'$-\sigma_{\textnormal{\tiny A}}^{xx}$')
axs[1].plot(xs1, sigIKprof[:,0],               lw=lw, color=c_sIK,   label=r'$-\sigma_{\textnormal{\tiny IK}}^{xx}$')
axs[1].fill_between(xs,0*xs,FAprof[:,0]/dx,    lw=0 , color=c_intFA, label=None, alpha=0.8,zorder=3)
axs[1].plot(xs , FAprof[:,0]/dx,               lw=lw, color=c_FA,    label=r'$\delta F_{\textnormal{\tiny A}}^x$',zorder=3)
axs[1].text((vxs[0]+vxs[1])/2.0+shiftx+1,-np.max(FAprof[:,0]/dx)*0.1, r'$\Delta p_1='+'{:.1f}'.format(dp1)+r'$', color=c_intFA, zorder=4,va='top',ha='center')
axs[1].text((vxs[2]+vxs[3])/2.0+shiftx-3, np.max(FAprof[:,0]/dx)*0.1, r'$\Delta p_2='+'{:.1f}'.format(dp2)+r'$', color=c_intFA, zorder=4,va='bottom',ha='center')
axs[1].legend(frameon=False, bbox_to_anchor=(0.99, 0.62), loc='center left', handletextpad=htp,labelspacing=lbs,handlelength=hl)
axs[1].set_ylim(-np.max(FAprof[:,0]/dx)*0.8)

# force plot
axs[2].plot(xs,FAprof[:,0], color=c_FA, lw=lw,label=r'$\delta F_{\textnormal{\tiny A}}^x$')
axs[2].plot(xs,Fintprof[:,0]-dsigAdxprof[:,0],color=c_stot,ls='--',lw=lw,zorder=2,label=r'$\partial_x \sigma_{\textnormal{\tiny tot}}^{xx}$')
axs[2].legend(frameon=False, bbox_to_anchor=(0.99, 0.5), loc='center left',
              handletextpad=htp,labelspacing=lbs,handlelength=hl)

# current plot
axs[3].plot(xsbin,Jsbin/p.Lx,color=c_J, lw=lw, label=r'$J_{\rm stress}$',zorder=1)
axs[3].plot(xs,[Jobs/p.Lx]*xs.shape[0],color=c_Jobs,label=r'$J_{\rm obs}$',ls='--')
axs[3].legend(frameon=False, bbox_to_anchor=(0.99, 0.5), loc='center left',
              handletextpad=htp,labelspacing=lbs,handlelength=hl)
axs[3].set_ylim(Jobs*2/p.Lx,0)
axs[3].set_xlabel(r'$x$')
axs[3].xaxis.set_label_coords(0.89,-0.075)
axs[3].tick_params(axis='x',which='both',direction='in',pad=1,length=3)

# adjust axis limits
for ax in axs:
    ax.set_xlim(0,p.Lx)
for ax in axs[:3]:
    ax.set_xticks([])
    
# adjust axis tick params
for ax in axs:
    ax.tick_params(axis='y', which='both', length=2, direction='in',pad=2)

# add (a), (b), etc.
edgewidth = 0.5
pad = 1.25
fontsize=8.0
labels=['(a)','(b)','(c)','(d)']
for i,ax in enumerate(axs):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    pad_x = pad / 72 / width
    pad_y = pad / 72 / height
    inside_pad = pad/fontsize
    text = ax.text(pad_x,1-pad_y,labels[i],ha='left',va='top',transform=ax.transAxes,
                       bbox=dict(boxstyle='square,pad='+str(inside_pad),edgecolor='black',facecolor='white',lw=edgewidth))
    text.set_clip_on(True)

plt.subplots_adjust(top=0.995,bottom=0.05,right=1.185,left=0.115)

plt.savefig('Lx_40-Ly_15-N_400-vs_1.0_30.0-xs_0.0_4.0_14.0_30.0-k_50.0-tf_700000-dt_0.000500-stress_fig-binw{:d}.pdf'.format(binw),dpi=1000)
plt.savefig('Lx_40-Ly_15-N_400-vs_1.0_30.0-xs_0.0_4.0_14.0_30.0-k_50.0-tf_700000-dt_0.000500-stress_fig-binw{:d}.png'.format(binw),dpi=1000)
