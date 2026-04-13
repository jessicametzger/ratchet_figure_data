import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os,sys
from pdf2image import convert_from_path
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

# check if something is a number
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


# read active data
ac_path = 'data/ABPs'
ac_int_path= ac_path+'/active-Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025'    # interacting
ac_ni_path = ac_path+'/active-ni_Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025' # non-interacting

a_prof = np.genfromtxt(ac_int_path+'-prof_avg',delimiter='\t')   # density profile (active, interacting)
a_ni_prof = np.genfromtxt(ac_ni_path+'-prof_avg',delimiter='\t') # density profile (active, non-interacting)
a_disp = np.genfromtxt(ac_int_path+'-disp_avg',delimiter='\t')   # net displacement (active, interacting)
a_ni_disp = np.genfromtxt(ac_ni_path+'-disp_avg',delimiter='\t') # net displacement (active, non-interacting)

# read passive data
pa_path = 'data/PBPs'
pa_int_path = pa_path+'/passive-Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025'   # interacting
pa_ni_path = pa_path+'/passive-ni_Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025' # non-interacting

p_prof = np.genfromtxt(pa_int_path+'-prof_avg',delimiter='\t')   # density profile (passive, interacting)
p_ni_prof = np.genfromtxt(pa_ni_path+'-prof_avg',delimiter='\t') # density profile (passive, non-interacting)
p_disp = np.genfromtxt(pa_int_path+'-disp_avg',delimiter='\t')   # net displacement (passive, interacting)
p_ni_disp = np.genfromtxt(pa_ni_path+'-disp_avg',delimiter='\t') # net displacement (passive, non-interacting)

# import params (same for active and passive, except v <-> T)
p = get_params(ac_int_path+'-param')

pd = p.__dict__
vxs = pd['vx1,vx2,...']
vs = pd['v1,v2,...']
p.vs = vs
p.vxs = vxs


##########################################################################################################################################
# ACTIVITY LANDSCAPE

# activity landscape function
def v(xs,p):
    bs = np.copy(p.vs)*0
    for i in range(p.Nstep_v-1):
        dx = p.vxs[i+1]-p.vxs[i]
        if dx==0:
            dx=1
        bs[i] = (p.vs[i+1] - p.vs[i])/dx
    dx = p.vxs[0]+p.Lx-p.vxs[-1]
    if dx==0:
        dx=1
    bs[-1] = (p.vs[0]-p.vs[-1])/dx
    
    ans = np.ones(xs.shape)
    inds = xs<p.vxs[0]
    dxs = xs[inds]+p.Lx-p.vxs[-1]
    ans[inds] = p.vs[-1] + bs[-1]*dxs
    for i in range(1,p.Nstep_v+1):
        if i==p.Nstep_v:
            upper = p.Lx+p.vxs[0]
        else:
            upper = p.vxs[i]
        inds = (xs>=p.vxs[i-1]) & (xs<upper)
        dxs = xs[inds] - p.vxs[i-1]
        ans[inds] = p.vs[i-1] + bs[i-1]*dxs
    return ans

# v(x) profile and 1/v(x)
xs_theory = np.linspace(0,p.Lx,1000)
vs = v(xs_theory,p)
vs_norm = np.trapezoid(vs, xs_theory)
inv_vs = 1/vs
inv_vs_norm = np.trapezoid(inv_vs, xs_theory)
dx_theory = xs_theory[1]-xs_theory[0]


##########################################################################################################################################
# DENSITY PROFILES

dx = p.Lx/float(p.Nbinx)
xs = np.linspace(0,p.Lx-dx,p.Nbinx)
norm = p.N*p.NstepProf*dx

# normalization for profiles from requiring integral of rho is 1
a_norm = np.sum(a_prof[:,1])*dx
a_ni_norm = np.sum(a_ni_prof[:,1])*dx
p_norm = np.sum(p_prof[:,1])*dx
p_ni_norm = np.sum(p_ni_prof[:,1])*dx


##########################################################################################################################################
# PLOTTING

# colors
active_color='#ff9900ff'
passive_color='cornflowerblue'

fig=plt.figure(figsize=(3.4,3.4*0.95))

# Create a GridSpec with 4 rows and 1 column
gs = gridspec.GridSpec(5, 2, height_ratios=[2.2,0.5,1, 0.25, 1], width_ratios = [1,0.33], hspace=0)

# Create the three subplots that share the x-axis
ax2 = fig.add_subplot(gs[1, :])
ax0 = fig.add_subplot(gs[2, :],sharex=ax2)
ax1 = fig.add_subplot(gs[4, 0])
axs = [ax0,ax1,ax2]

# add schematic panel
pdf_path = 'vofx-bacteria-schematic-1.pdf'
pages = convert_from_path(pdf_path,dpi=1000)
ax = fig.add_subplot(gs[0, :])
ax.imshow(pages[0])
ax.axis('off')

# plot v(x)
axs[2].plot(xs_theory, vs, color='#6aa84f',lw=1.5)

# plot cumulative density profile
axs[0].plot(xs+dx/2., a_prof[:,1]/a_norm, color=active_color, label=r'$\rho_{\textnormal{\tiny ABP}}$',lw=1.5)
axs[0].plot(xs+dx/2., p_prof[:,1]/p_norm, color=passive_color, label=r'$\rho_{\textnormal{\tiny OBP}}$',lw=1.5)

# plot theory 1/v -- indistinguishable from measured non-interacting density profile
axs[0].plot(xs_theory, inv_vs/inv_vs_norm, color='black',zorder=-1,lw=1.5,ls='--', 
            label=r'$\rho_{\textnormal{\tiny NI}}\propto \frac{1}{v}=\frac{1}{T}$')

# plot displacements
axs[1].plot(a_disp[:,0], 100*a_disp[:,1]/80., color=active_color, lw=1.5,label=r'$\Delta x_{\textnormal{\tiny ABP}}(t)$')
axs[1].plot(a_ni_disp[:,0], 100*a_ni_disp[:,1]/80., color=active_color, lw=1.5, ls='--',label=r'$\Delta x_{\textnormal{\tiny ABP,NI}}(t)$')
axs[1].plot(p_disp[:,0], 100*p_disp[:,1]/80., color=passive_color, lw=1.5,label=r'$\Delta x_{\textnormal{\tiny OBP}}(t)$')
axs[1].plot(p_ni_disp[:,0], 100*p_ni_disp[:,1]/80., color=passive_color, lw=1.5, ls='--',label=r'$\Delta x_{\textnormal{\tiny OBP,NI}}(t)$')

# create legend in each figure
legends  = [axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.45),frameon=False,
              handletextpad=0.4,labelspacing=0.5,handlelength=1.75)]
legends += [axs[0].legend(loc='center left', frameon=False, bbox_to_anchor=(-0.01,0.4),
              handletextpad=0.4,labelspacing=0.5,handlelength=1.75,ncols=2,columnspacing=0.9)]

# shorten/adjust axis ticks
axs[2].tick_params(axis='x',which='both',length=0)
plt.setp(axs[2].get_xticklabels(), visible=False)
axs[0].set_xticks([])
for ax in axs:
    ax.tick_params(axis='both',which='both',direction='in',pad=2)

# manually set yticks
axs[0].set_xticks([0,20,40,60])
axs[0].set_yticks([0,0.02])
axs[1].set_yticks([0,10],[r'0',r'$10$'])
axs[1].set_xticks([0,4e5,8e5],labels=['0',r'$4\times 10^5$',r'$8\times 10^5$'])
axs[2].set_yticks([0,10,20])

# set axis labels 
axs[0].set_xlabel(r'$x$')
axs[0].set_ylabel(r'$\rho(x)$',rotation='horizontal')
axs[1].set_xlabel(r'$t$')
axs[1].set_ylabel(r'$\frac{\Delta x(t)}{L_x}$',rotation='horizontal')
axs[2].set_ylabel(r'$v(x)$ or $T(x)$',rotation='horizontal',color='#4b7a36')

# position axis labels
axs[0].xaxis.set_label_coords(0.98,-0.05)
axs[0].yaxis.set_label_coords(-0.05,0.65)
axs[1].xaxis.set_label_coords(0.98,-0.05)
axs[1].yaxis.set_label_coords(-0.07,0.7)
axs[2].yaxis.set_label_coords(0.18,0.22)

# set axis limits
axs[0].set_xlim(xmin=0,xmax=p.Lx)
axs[0].set_ylim(ymin=0,ymax=np.max(a_prof[:,1]/a_norm)*1.15)
axs[1].set_xlim(xmin=0,xmax=np.amax(a_disp[:,0]))
axs[2].set_ylim(ymin=0,ymax=28)

# Adjust subplot parameters to remove vertical spacing between top three plots
plt.subplots_adjust(left=0.09,bottom=0.05,top=0.99,right=0.99)

plt.savefig('ratchet-exception-fig.pdf',dpi=1000)
plt.savefig('ratchet-exception-fig.png',dpi=1000)


##########################################################################################################################################