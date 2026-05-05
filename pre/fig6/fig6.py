import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os,sys
from pdf2image import convert_from_path
from types import SimpleNamespace
import re
from scipy.ndimage import convolve

plt_rcparams = {'text.usetex' : True,
                'font.size' : 8,
                'font.family' : 'serif',
                'text.latex.preamble' : r"\usepackage{lmodern} \usepackage{amstext}",
                'figure.figsize' : [3.4,3.4*0.7],
                'figure.dpi': 200
                }
plt.rcParams.update(plt_rcparams)


################################################################################################################
# AUXILIARY FUNCTIONS

def isnum(s):
    try:
        a=float(s)
    except ValueError:
        return False
    return True

# import parameter data from a simulation
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



################################################################################################################
# IMPORT DATA

# we will create a figure for this N*epsilon
N_epsilon = 15
N_epsilon_i = int(N_epsilon / 3.0 - 1 + 1e-5) # which index in list 3,6,9,...

N_folders = [os.path.join('data',x) for x in os.listdir('data') if x.startswith('N') and os.path.isdir(os.path.join('data',x))]
N_folders = sorted(N_folders, key=lambda x: int(x.split('/')[-1].replace('N',''))) 
Ns = [int(x.split('/')[-1].replace('N','')) for x in N_folders]

# get current, density, and magnetization profile data
Js_all_N = []
disps_all_N = []
rhos_all_N = []
ms_all_N = []
for N_folder in N_folders:
    N_int = int(N_folder.split('/')[-1].replace('N',''))

    # current data
    Js_all = np.genfromtxt(os.path.join(N_folder,'J_avg'),delimiter='\t')
    Js_all_N += [Js_all[N_epsilon_i,1:]]

    # displacement data
    file = [os.path.join(N_folder,x) for x in os.listdir(N_folder) if x.endswith('-avg_disp') and int(N_int*float(x.split('amp-')[1].split('-')[0])*2/3.+1e-5)==N_epsilon][0]
    disps_all_N += [np.genfromtxt(file)]

    # density profile
    file = [os.path.join(N_folder,x) for x in os.listdir(N_folder) if x.endswith('-avg_profile_rho') and int(N_int*float(x.split('amp-')[1].split('-')[0])*2/3.+1e-5)==N_epsilon][0]
    rho = np.genfromtxt(file)
    norm = np.trapezoid(rho[:,1],rho[:,0])
    rho[:,1] /= norm
    rhos_all_N += [rho]

    # magnetization profile
    file = [os.path.join(N_folder,x) for x in os.listdir(N_folder) if x.endswith('-avg_profile_m') and int(N_int*float(x.split('amp-')[1].split('-')[0])*2/3.+1e-5)==N_epsilon][0]
    m = np.genfromtxt(file)
    m[:,1] /= norm
    ms_all_N += [m]

dx = rhos_all_N[0][1,0] - rhos_all_N[0][0,0]

# get parameters of a sample file
sample_seed = [os.path.join(N_folders[0],x) for x in os.listdir(N_folders[0]) if 'param' in x][0]
p = get_params(sample_seed)

################################################################################################################
# CALCULATE THEORETICAL CURRENT/DENSITY/MAGNETIZATION, PERTURBATIVELY

# convolution routine
# computes integral of arr(x+y)*ker(y), with PBC
# Must flip ker because of how we write force calculation in the theory
def conv(xs, ker, arr):
    return convolve(arr, ker, mode='wrap')*(xs[1]-xs[0])
    
# to append to beginning of indefinite integral arrays
stub=np.array([0])

# activity landscape
def v(x,v0,a1,a2,P,L):
    return v0*(1 + 0.5*a1*np.sin(2*np.pi*x*P/L) + .5*a2*np.sin(4*np.pi*x*P/L))

# give parameters, return density profile and current expansion coefficients
def integrate(v0,a1,a2,L,P,n_orders,sigma,res):
    
    res = int(res)
    
    # number of points in the domain and their spacing
    Nx = res * P
    dx = L/float(Nx)
    
    # evenly spaced Nx points from 0 to L, excluding L; shifted to the middle of the interval
    xs = np.arange(0, L, L/float(Nx))+dx/2.
    
    # sum arrays numpy-style
    def npsum(arrs):
        if len(arrs)==0:
            return xs*0
        return sum(arrs)
    
    vs = v(xs, v0, a1, a2, P, L)
    
    C = 1/np.trapezoid(1/vs, xs)
    rho0_ = C/vs
    
    # create array for the interaction kernel
    xs_kernel = np.arange(-sigma, sigma-1e-10, dx)+dx/2.
    # make sure the interaction kernel has odd shape for symmetry of convolution
    if xs_kernel.shape[0] % 2 == 0:
        xs_kernel = np.arange(-sigma-dx/2., sigma+dx/2., dx)+dx/2.
    
    # interaction potential that integrates to zero 
    U_hat = (3/(4*sigma))*(1 - xs_kernel**2/sigma**2)
    f_hat = (3/2)*xs_kernel/sigma**3
    
    # make sure interaction kernel is zero outside range
    U_hat[np.where(np.abs(xs_kernel)>=sigma)] = 0
    f_hat[np.where(np.abs(xs_kernel)>=sigma)] = 0
    
    # convolution version
    f_conv_rhon = [conv(xs, f_hat, rho0_)]
    Jn = [0]
    omegan = [C]
    rhon = [rho0_]
    mn = [0*xs]
    
    # integrate the EOM to calculate the density, current, etc. order by order
    for n in range(1,n_orders+1):

        # current
        Jn_integrand = rhon[0]*npsum([rhon[k]*f_conv_rhon[n-k-1] for k in range(n)])
        Jn += [np.trapezoid(Jn_integrand, xs)]

        # magnetization
        mn += [(Jn[n] - npsum([rhon[k]*f_conv_rhon[n-k-1] for k in range(n)]))/vs]

        # condition that rhon integrates to zero
        mnint = np.concatenate((stub,np.cumsum(dx*(mn[n][1:]+mn[n][:-1])/2.)))
        omegan_integrand = (npsum([mn[k]*f_conv_rhon[n-k-1] for k in range(n)]) + mnint)/vs
        omegan += [np.trapezoid(omegan_integrand, xs)*omegan[0]]

        # density
        rhon += [-omegan_integrand + omegan[n]/vs]
        

        # save density convolved w interaction kernel for future steps
        f_conv_rhon += [conv(xs, f_hat, rhon[n])]
    
    return [xs, Jn, rhon, mn]


# calculate current
n_orders=6
res = p.L*1000 # number of points in domain of numerical integrator
xs_th,Jns,rhons,mns = integrate(p.v0,p.a1,p.a2,p.L,p.P,n_orders,p.sigma,res)
dx_th = xs_th[1]-xs_th[0]

# should return the following values:
# Jns = [0.0, 0.00016024329392465422, -8.134845248146115e-07, -1.4033438024323694e-07, 4.648862326420472e-10, -4.5923807968686305e-11]

# sum to get theory
J_th = sum([Jns[i]*N_epsilon**i for i in range(n_orders)])
rho_th = sum([rhons[i]*N_epsilon**i for i in range(n_orders)])
m_th = sum([mns[i]*N_epsilon**i for i in range(n_orders)])

################################################################################################################
# PLOTTING

cmap = plt.get_cmap('plasma')
cmap_N = plt.get_cmap('plasma')

# create figure and subplots
fig=plt.figure(figsize=(3.4,3.5))
gs = gridspec.GridSpec(7, 1, height_ratios=[.4,.1,1,.1, 1,.2,.6], hspace=0)
v_ax = fig.add_subplot(gs[0])
ax0 = fig.add_subplot(gs[2])
ax1 = fig.add_subplot(gs[4])
ax2 = fig.add_subplot(gs[6])
axs = [ax0,ax1,ax2]

# plot the density profiles for panel (b)
for k in range(len(rhos_all_N)):
    axs[0].plot(rhos_all_N[k][:,0]+dx/2., rhos_all_N[k][:,1], color=cmap_N(k/float((len(N_folders)-1))),lw=2,
             label=r'$\overline{\rho}$='+'{:.1f}'.format(float(Ns[k])/10.0))
             
    # plot the theory in the middle so the legend appears in the bottom left
    if k==int(len(rhos_all_N)/2.-1+1e-5):
        rho_th_plot, = axs[0].plot(xs_th+dx_th/2., rho_th, color='black',lw=2,label=r'theory',zorder=4,ls='--')

# create legend for panel (b)
top_legend = axs[0].legend(loc='upper left',frameon=False, ncol=2, bbox_to_anchor=(0,0.9))
ylim = axs[0].get_ylim()
axs[0].set_ylim(ymax=ylim[1]*1.1)
axs[0].set_ylabel(r'$\rho_{\rm s}(x)$',rotation='horizontal')
axs[0].yaxis.set_label_coords(-0.05,0.8)
axs[0].set_yticks([0,0.1,0.2],labels=['0','0.1','0.2'])
axs[0].set_xlim(0,10)

# create v(x) for panel (a) and theoretical rho_0(x) for panel (c) - first, with the x grid used in the data
xs_data = rhos_all_N[0][:,0] + dx/2.
v_inv_data = 1/(p.v0*(1 + 0.5*p.a1*np.sin(2*np.pi*xs_data/p.L) + 0.5*p.a2*np.sin(4*np.pi*xs_data/p.L)))
v_inv_data_norm = np.trapezoid(v_inv_data, xs_data)
v_inv_data = v_inv_data/v_inv_data_norm

# v(x) for panel (a) and theoretical rho_0(x) in the x grid used for theory
vs_theory = (p.v0*(1 + 0.5*p.a1*np.sin(2*np.pi*xs_th/p.L) + 0.5*p.a2*np.sin(4*np.pi*xs_th/p.L)))
v_inv_theory = 1/vs_theory
v_inv_theory_norm = np.trapezoid(v_inv_theory, xs_th)
v_inv_theory = v_inv_theory/v_inv_theory_norm

# plot v(x) in panel (a)
v_ax.plot(xs_th+dx_th/2., vs_theory, color='green',lw=2)
v_ax.set_xlim(xmin=0,xmax=p.L)
v_ax.set_ylim(ymin=0,ymax=v_ax.get_ylim()[1]*1.4)
v_ax.set_ylabel(r'$v(x)$',rotation='horizontal')
v_ax.yaxis.set_label_coords(-0.05,0.7)
v_ax.set_yticks([0,20],labels=['0','20'])

# plot rho(x)-rho_0(x) and m(x)-m_0(x) for panel c
for k in range(len(rhos_all_N)):
    axs[1].plot(rhos_all_N[k][:,0]+dx/2., rhos_all_N[k][:,1]-v_inv_data, color=cmap_N(k/float((len(N_folders)-1))),lw=2,
             label=r'$\overline{\rho}$='+'{:.1f}'.format(float(Ns[k])/10.0))

    axs[1].plot(ms_all_N[k][:,0]+dx/2., ms_all_N[k][:,1], color=cmap_N(k/float((len(N_folders)-1))),lw=2,ls='--')

# panel c: theory version
rho_th_plot, = axs[1].plot(xs_th+dx_th/2., rho_th-v_inv_theory, color='black',lw=2,label=r'theory, $\rho_{\rm s}-\rho_0$')
m_th_plot, = axs[1].plot(xs_th+dx_th/2., m_th, color='black',lw=2,ls='--',label=r'theory, $m_{\rm s}-m_0$')

# plot x displacement vs time
k=4
tf = disps_all_N[k][-1,0]
axs[2].plot([0]+list(disps_all_N[k][39::40,0]), [0]+list(disps_all_N[k][39::40,1]), marker='o', markersize=5, ls='-', lw=0, 
            color=cmap_N(k/float((len(N_folders)-1))), label=r'$\overline{\rho}$='+'{:.1f}'.format(float(Ns[k])/10.0))

axs[2].plot([0,tf],[0,J_th*tf], color='black', label='theory', lw=2,zorder=0)

# panel (a) axis ticks and labels    
v_ax.set_xticks([])
v_ax.set_xticklabels([])
v_ax.tick_params(axis='x', which='both', length=0, direction='in')
v_ax.tick_params(axis='y', which='both', length=3, direction='in')

# panel (b) axis ticks, limits, and labels
axs[0].set_xticks([])
axs[0].set_xticklabels([])
axs[0].set_ylim(ymin=0)
axs[0].tick_params(axis='x', which='both', length=0)
axs[0].tick_params(axis='y',direction='in',pad=2,length=3)

# panel (c) axis ticks, limits, and labels
axs[1].set_xlim(xmin=0,xmax=p.L)
axs[1].set_xticks([0,2.5,5,7.5], ['0','2.5','5','7.5'])
axs[1].set_xlabel(r'$x$')
axs[1].set_yticks([-0.1,0],labels=['-0.1','0'])
axs[1].xaxis.set_label_coords(0.95,-0.04)
axs[1].text(-1.54,0.065,r'$\rho_{\rm s}-\rho_0,$')
axs[1].text(-1.63,0.04,r'$m_{\rm s}-m_0$')
axs[1].tick_params(axis='y',direction='in',pad=2,length=3)
axs[1].tick_params(axis='x',pad=2,direction='in',length=3)
axs[1].legend(handles=[rho_th_plot, m_th_plot], loc='lower left',frameon=False)
fig.subplots_adjust(hspace=0, wspace=0)

# panel (d) axis ticks, limits, and labels
axs[2].set_xlim(xmin=0,xmax=tf)
axs[2].set_ylim(ymin=0)
axs[2].legend(frameon=False,loc='upper left',bbox_to_anchor=(0.05,1))
axs[2].tick_params(axis='y',direction='in',pad=2,length=3)
axs[2].tick_params(axis='x',pad=2,direction='in',length=3)
axs[2].set_xlabel(r'$t$')
axs[2].set_ylabel(r'$\frac{\langle \Delta x(t)\rangle}{L_x}$',rotation='horizontal')
axs[2].xaxis.set_label_coords(0.92,-0.06)
axs[2].yaxis.set_label_coords(-0.08,0.75)
axs[2].set_xticks([0,1e6,2e6,3e6,4e6],labels=[0,1,2,3,4])
axs[2].text(tf*0.88,200,r'$\times 10^6$')

# add (a), (b), etc.
edgewidth = 0.5
pad = 1.25
fontsize=8.0
labels=['(a)','(b)','(c)','(d)']
for i,ax in enumerate([v_ax]+axs):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    pad_x = pad / 72 / width
    pad_y = pad / 72 / height
    inside_pad = pad/fontsize
    text = ax.text(pad_x,1-pad_y,labels[i],ha='left',va='top',transform=ax.transAxes,
                   bbox=dict(boxstyle='square,pad='+str(inside_pad),edgecolor='black',facecolor='white',lw=edgewidth))
text.set_clip_on(True)

plt.setp(axs[0].get_xticklabels(), visible=False)
plt.setp(v_ax.get_xticklabels(), visible=False)

plt.subplots_adjust(left=0.15,right=0.98,bottom=0.04,top=0.98)

plt.savefig('Nepsilon'+str(N_epsilon)+'_prof_disp.pdf',dpi=400)
