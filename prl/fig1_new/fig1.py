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
ac_int_path= 'active-Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025'    # interacting
ac_ni_path = 'active-ni_Lx80_Ly50_v11_vcl2.5_vcr22.5_vr35_P1_N800_dt0.00025' # non-interacting

a_prof = np.genfromtxt(ac_int_path+'-prof_avg',delimiter='\t')   # density profile (active, interacting)
a_ni_prof = np.genfromtxt(ac_ni_path+'-prof_avg',delimiter='\t') # density profile (active, non-interacting)
a_disp = np.genfromtxt(ac_int_path+'-disp_avg',delimiter='\t')   # net displacement (active, interacting)
a_ni_disp = np.genfromtxt(ac_ni_path+'-disp_avg',delimiter='\t') # net displacement (active, non-interacting)

# read passive data
pa_int_path = 'passive-Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025'   # interacting
pa_ni_path = 'passive-ni_Lx80_Ly50_T11_Tcl2.5_Tcr22.5_Tr35_P1_N800_dt0.00025' # non-interacting

p_prof = np.genfromtxt(pa_int_path+'-prof_avg',delimiter='\t')   # density profile (passive, interacting)
p_ni_prof = np.genfromtxt(pa_ni_path+'-prof_avg',delimiter='\t') # density profile (passive, non-interacting)
p_disp = np.genfromtxt(pa_int_path+'-disp_avg',delimiter='\t')   # net displacement (passive, interacting)
p_ni_disp = np.genfromtxt(pa_ni_path+'-disp_avg',delimiter='\t') # net displacement (passive, non-interacting)

# import params (same for active and passive, except v <-> T)
p = get_params(ac_int_path+'-param')



##########################################################################################################################################
# ACTIVITY LANDSCAPE

# step activity pattern
def v(x,v0,v1,v_center1, v_center2,v_right,Lx,P):
    motif_L = Lx/P
    n_motif = int(x/motif_L)
    x_in_motif = x - motif_L*n_motif
    
    if x_in_motif>v_right:
        return v0
    elif x_in_motif<v_center1:
        return v0 + (v1-v0)*x_in_motif/v_center1
    elif x_in_motif<v_center2:
        return v1
    else:
        return v1 - (v1-v0)*(x_in_motif-v_center2)/(v_right-v_center2)

# parameters of activity landscape
dxR = p.v_right - p.v_center2
plateau = p.v_center2 - p.v_center1
vpL = (p.v1-p.v0)/p.v_center1
vpR = (p.v0-p.v1)/dxR
v_params = [p.v0,p.v1,p.v_center1,p.v_center2,p.v_right,p.Lx,p.P]

# the original x=0 point will be shifted to Lx*shift_factor
shift_factor = 0.5
shift = shift_factor*p.Lx
rshift = p.Lx-shift

# v(x) profile and 1/v(x)
xs_theory = np.linspace(0,p.Lx,1000)
vs = np.array([v(x,*v_params) for x in xs_theory])
vs_norm = np.trapezoid(vs, xs_theory)
inv_vs = 1/vs
inv_vs_norm = np.trapezoid(inv_vs, xs_theory)
dx_theory = xs_theory[1]-xs_theory[0]

# shifted v(x) profile and 1/v(x)
shifted_vs = np.concatenate((vs[xs_theory>=rshift], vs[xs_theory<rshift]),axis=0)
shifted_inv_vs = np.concatenate((inv_vs[xs_theory>=rshift], inv_vs[xs_theory<rshift]),axis=0)


##########################################################################################################################################
# DENSITY PROFILES

dx = p.Lx/float(p.Nbinx)
norm = p.N*p.NstepProfile*dx

# normalization for profiles from requiring integral of rho is 1
a_norm = np.sum(a_prof[:,1])*dx
a_ni_norm = np.sum(a_ni_prof[:,1])*dx
p_norm = np.sum(p_prof[:,1])*dx
p_ni_norm = np.sum(p_ni_prof[:,1])*dx

# shift the x axis
original_xs = a_prof[:,0]
shifted_xs = np.concatenate((original_xs[np.where(original_xs>=rshift)],original_xs[np.where(original_xs<rshift)]), axis=0)

# shift the density profiles along x
a_prof = np.concatenate((a_prof[np.where(a_prof[:,0]>=rshift)],a_prof[np.where(a_prof[:,0]<rshift)]), axis=0)
p_prof = np.concatenate((p_prof[np.where(p_prof[:,0]>=rshift)],p_prof[np.where(p_prof[:,0]<rshift)]), axis=0)
a_ni_prof = np.concatenate((a_ni_prof[np.where(a_ni_prof[:,0]>=rshift)],a_ni_prof[np.where(a_ni_prof[:,0]<rshift)]), axis=0)
p_ni_prof = np.concatenate((p_ni_prof[np.where(p_ni_prof[:,0]>=rshift)],p_ni_prof[np.where(p_ni_prof[:,0]<rshift)]), axis=0)

# boundaries of regions with nonzero v'(x)
left_inds = np.where(shifted_xs<=p.v_center1)
right_inds = np.where((shifted_xs>=p.v_center2) & (shifted_xs<=p.v_right)) 
outside_inds = np.where((shifted_xs>=p.v_center1) & (shifted_xs<p.v_center2))
left_inds = (np.concatenate((np.array([left_inds[0][0]-1]), left_inds[0])),)
right_inds = (np.concatenate((np.array([right_inds[0][0]-1]), right_inds[0])),)

# activity landscape in data coordinates
vs_data = np.array([v(x+dx/2.,*v_params) for x in original_xs])
vps_data = 0*vs_data
vps_data[left_inds] = vpL
vps_data[right_inds] = vpR
inv_vs_data = 1/vs_data
vs_data_norm = np.sum(vs_data)*dx
inv_vs_data_norm = np.sum(inv_vs_data)*dx
shifted_vs_data = np.concatenate((vs_data[np.where(original_xs>=rshift)],vs_data[np.where(original_xs<rshift)]), axis=0)
shifted_inv_vs_data = np.concatenate((inv_vs_data[np.where(original_xs>=rshift)], inv_vs_data[np.where(original_xs<rshift)]), axis=0)



##########################################################################################################################################
# PLOTTING

# colors
active_color='#ff9900ff'
passive_color='cornflowerblue'

fig=plt.figure(figsize=(3.4,3.4*0.84))

# Create a GridSpec with 4 rows and 1 column
gs = gridspec.GridSpec(5, 2, height_ratios=[1.6,0.5,1, 0.25, 1], width_ratios = [1,0.33], hspace=0)

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
axs[2].plot(xs_theory, shifted_vs, color='#6aa84f',lw=1.5)

# plot cumulative density profile
axs[0].plot(original_xs+dx/2., a_prof[:,1]/a_norm, color=active_color, label=r'$\rho_{\textnormal{\tiny ABP}}$',lw=1.5)
axs[0].plot(original_xs+dx/2., p_prof[:,1]/p_norm, color=passive_color, label=r'$\rho_{\textnormal{\tiny OBP}}$',lw=1.5)

# plot theory 1/v -- indistinguishable from measured non-interacting density profile
axs[0].plot(xs_theory, shifted_inv_vs/inv_vs_norm, color='black',zorder=-1,lw=1.5,ls='--', 
            label=r'$\rho_{\textnormal{\tiny NI}}\propto \frac{1}{v}=\frac{1}{T}$')

# plot displacements
axs[1].plot(a_disp[:,0], a_disp[:,1]/80., color=active_color, lw=1.5,label=r'$\Delta x_{\textnormal{\tiny ABP}}(t)$')
axs[1].plot(a_ni_disp[:,0], a_ni_disp[:,1]/80., color=active_color, lw=1.5, ls='--',label=r'$\Delta x_{\textnormal{\tiny ABP,NI}}(t)$')
axs[1].plot(p_disp[:,0], p_disp[:,1]/80., color=passive_color, lw=1.5,label=r'$\Delta x_{\textnormal{\tiny OBP}}(t)$')
axs[1].plot(p_ni_disp[:,0], p_ni_disp[:,1]/80., color=passive_color, lw=1.5, ls='--',label=r'$\Delta x_{\textnormal{\tiny OBP,NI}}(t)$')

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