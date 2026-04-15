import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os,sys
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


N_folders = [os.path.join('data',x) for x in os.listdir('data') if x.startswith('N') and os.path.isdir(x)]
N_folders = sorted(N_folders, key=lambda x: int(x.split('/')[-1].replace('N',''))) 

Js_all_N = []
for N_folder in N_folders:
    N = int(N_folder.split('/')[-1].replace('N',''))
    Js_all = np.genfromtxt(os.path.join(N_folder,'J_avg'),delimiter='\t')
    Js_all_N += [Js_all]


# get parameters of a sample file
sample_seed = [os.path.join(N_folders[0],x) for x in os.listdir(N_folders[0]) if 'param' in x][0]
p = get_params(sample_seed)


################################################################################################################
# CALCULATE THEORETICAL CURRENT, PERTURBATIVELY

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
xs,Jns,rhons,mns = integrate(p.v0,p.a1,p.a2,p.L,p.P,n_orders,p.sigma,res)

# should return the following values:
# Jns = [0.0, 0.00016024329392465422, -8.134845248146115e-07, -1.4033438024323694e-07, 4.648862326420472e-10, -4.5923807968686305e-11]


################################################################################################################
# PLOTTING

# colors
cmap = plt.get_cmap('viridis')
cmap_N = plt.get_cmap('plasma')

# 2 plots side-by-side
fig = plt.figure(figsize=(8.6/2.54,8.6*0.4/2.54))
gs = gridspec.GridSpec(1, 2, width_ratios=[0.85, 1], wspace=0, hspace=0)
axs = [fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1])]

# add pre-made schematic to left plot
pdf_path = "mean-field-schematic-thin.pdf"
pages = convert_from_path(pdf_path, dpi=1000)  # Increase dpi for higher resolution
image = pages[0]
axs[0].imshow(image,extent=(0,1,0,1),aspect='equal')
axs[0].set_xlim(0,1)
axs[0].set_ylim(0.07,1)
axs[0].axis('off')

# plot simulation current
Ns = [x[0,0] for x in Js_all_N]
lines = []
for i in range(len(Js_all_N)):
    lines += [plt.errorbar(Js_all_N[i][:,1], Js_all_N[i][:,2],
                         color=cmap_N(i/float(len(Js_all_N)-1)),label=str(int(Js_all_N[i][1,0])),fmt='o',ls='',markersize=5,
                         mec='black',mew=0.5)]

# plot theoretical current
Nepsilons_grid = np.linspace(0, 21.5, 100)

# 1st order theoretical current
axs[1].plot(Nepsilons_grid, Jns[1]*Nepsilons_grid, color='gray', lw=2)

# 6th order theoretical current
axs[1].plot(Nepsilons_grid, sum([Jns[n]*Nepsilons_grid**n for n in range(len(Jns))]), color='black', lw=2)#, label='theory')

# axis labels
axs[1].set_xlabel(r'$N\varepsilon$')
axs[1].set_ylabel(r'$J$',rotation='horizontal')
axs[1].yaxis.set_label_coords(0.05,.5)
axs[1].xaxis.set_label_coords(0.94,-0.02)

# legend (which N)
legend = axs[1].legend(ncol=2,frameon=False,loc='lower right',title=r'$N$=',
                    bbox_to_anchor=(1.035,-0.035),handletextpad=0,columnspacing=0.75)
legend.get_title().set_position((-20, 0))

# axis limits
axs[1].set_ylim(ymin=0,ymax=0.0026)
axs[1].set_xlim(xmin=0,xmax=22)

# axis ticks
axs[1].set_yticks([0.001,0.002],labels=['0.001','0.002'])
axs[1].set_xticks([0,9,18],labels=['0','9','18'])
axs[1].tick_params(axis='y',direction='in',pad=-24)
axs[1].tick_params(axis='x',pad=2,direction='in')

plt.subplots_adjust(left=0,right=0.99,bottom=0.1,top=0.99,wspace=0)

plt.savefig('J_all_no_title.pdf',dpi=1000)