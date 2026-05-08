import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import rgb_to_hsv
from matplotlib import patches
import matplotlib.gridspec as gridspec
import os,sys
import copy
import matplotlib
from matplotlib.patches import FancyBboxPatch
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

cm = plt.get_cmap('plasma')


##########################################################################################################################################
# AUXILIARY FUNCTIONS FOR READING DATA

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

def avg_prof_1d(files,p,delimiter='\t'):
    '''
    Will read and average all files in list of filenames ls

    Assumes column 0 is time, columns 1 is x, everything else to be averaged
    '''
    if hasattr(p,'tf'):
        tf = p.tf
    elif hasattr(p,'final_time'):
        tf = p.final_time
    else:
        raise RuntimeError("No tf attribute")
    Nstore = int(tf / p.StoreInterProfile + 1e-5)
    Nstep = int(p.NstepProfile)
    norm = Nstore*Nstep

    nseed = len(files)
    norm*= nseed
    
    data = [np.genfromtxt(x,delimiter=delimiter) for x in files]
    
    examp = data[0]
    xs = np.unique(examp[:,1])
    size = xs.shape[0]

    # make sure there is the expected amount of data 
    assert np.all([p.Nbin*Nstore==x.shape[0] for x in data])

    # sum all the arrays together
    data_tot = sum([x[:,2:] for x in data])

    # sum it over time
    data_avg = sum([data_tot[size*i:size*(i+1),:] for i in range(Nstore)])
    ret = np.zeros((size,data_avg.shape[1]+1))
    ret[:,0] = xs
    ret[:,1:]= data_avg/norm
    return ret

path = 'data'
seeds= [os.path.join(path,x.replace('-param','')) for x in os.listdir(path) if x.endswith('-param')]
p = get_params(seeds[0]+'-param')

profs = [avg_prof_1d([x+'-profile_rho'],p) for x in seeds]
prof = sum(profs) / len(profs)
ms = [avg_prof_1d([x+'-profile_m'],p) for x in seeds]
m = sum(ms) / len(ms)

dx_data = prof[1,0]-prof[0,0]

disps = []
for x in seeds:
    disps += [np.genfromtxt(x+'-disp')]

# further parse displacement data
disps1 = [np.stack([x[i::p.N] for i in range(p.N)]) for x in disps]
disps1 = np.concatenate(disps1,axis=0)
disps2 = np.zeros((disps1.shape[1],3))
disps2[:,0] = disps1[0,:,0]
disps2[:,1] = np.mean(disps1[:,:,2],axis=0)
disps2[:,2] = np.std(disps1[:,:,2],ddof=1,axis=0)/np.sqrt(disps1.shape[0])
disps2[:,1:] = disps2[:,1:] / p.L


##########################################################################################################################################
#  ACTIVITY/POTENTIAL LANDSCAPE

def vofx(xs, p):
    vs = np.ones(xs.shape)
    for i in range(p.As.shape[0]):
        vs += p.As[i]*np.sin(2*np.pi*(i+1)*xs*p.P/p.L + 2*np.pi*p.dAs[i]*p.P)
    return vs*p.v0
def vpofx(xs, p):
    vps = np.zeros(xs.shape)
    for i in range(p.As.shape[0]):
        vps += p.As[i]*(2*np.pi*(i+1)*p.P/p.L)*np.cos(2*np.pi*(i+1)*xs*p.P/p.L + 2*np.pi*p.dAs[i]*p.P)
    return vps*p.v0

def Uofx(xs, p):
    Us = np.zeros(xs.shape)
    for i in range(p.Bs.shape[0]):
        Us += p.Bs[i]*np.sin(2*np.pi*(i+1)*xs*p.P/p.L + 2*np.pi*p.dBs[i]*p.P)
    return Us
def Upofx(xs, p):
    Ups = np.zeros(xs.shape)
    for i in range(p.Bs.shape[0]):
        Ups += p.Bs[i]*(2*np.pi*(i+1)*p.P/p.L)*np.cos(2*np.pi*(i+1)*xs*p.P/p.L + 2*np.pi*p.dBs[i]*p.P)
    return Ups
def Uppofx(xs, p):
    Upps = np.zeros(xs.shape)
    for i in range(p.Bs.shape[0]):
        Upps -= p.Bs[i] * (2*np.pi*(i+1)*p.P/p.L)**2 * np.sin(2*np.pi*(i+1)*xs*p.P/p.L + 2*np.pi*p.dBs[i]*p.P)
    return Upps

xs = np.linspace(0,p.L,400)
dx = xs[1]-xs[0]
vs = vofx(xs,p)
Us = Uofx(xs,p)
vps = vpofx(xs,p)
Ups = Upofx(xs,p)
Upps= Uppofx(xs,p)

##########################################################################################################################################
# CALCULATE THEORETICAL DENSITY, MAGNETIZATION, CURRENT, ETC.

# to append to beginning of indefinite integral arrays
stub=np.array([0])
    
def Teff(xs, Ups, vs, mu, alpha):
    return vs ** 2 / (mu * alpha) - mu * Ups ** 2 / alpha
    
def mueff(xs, Ups, Upps, vs, vps, mu, alpha):
    return mu / (1 - mu * Upps / alpha + mu * Ups * vps / (vs * alpha))
    
def Veff(xs, Us, Ups, vs, vps, mu, alpha):
    Veffs = Us - vs ** 2 / (2 * mu * alpha)
    Veffs_nonlocal_integrand = (Ups ** 2) * vps / vs
    Veffs_nonlocal = np.concatenate(
        (stub, np.cumsum((xs[1:] - xs[:-1]) * (Veffs_nonlocal_integrand[1:] + Veffs_nonlocal_integrand[:-1]) / 2.)),
        axis=0)
    return Veffs + (mu / alpha) * Veffs_nonlocal
    

def Phi(xs, Us, Ups, Upps, vs, vps, mu, alpha):
    stub = np.array([0])

    Teffs = vs ** 2 / (mu * alpha) - mu * Ups ** 2 / alpha
    mueffs = mu / (1 + mu * Upps / alpha - mu * Ups * vps / (vs * alpha))
    Veffs = Us - vs ** 2 / (2 * mu * alpha)
    Veffs_nonlocal_integrand = (Ups ** 2) * vps / vs
    Veffs_nonlocal = np.concatenate(
        (stub, np.cumsum((xs[1:] - xs[:-1]) * (Veffs_nonlocal_integrand[1:] + Veffs_nonlocal_integrand[:-1]) / 2.)),
        axis=0)
    Veffs = Veffs + (mu / alpha) * Veffs_nonlocal
    Veffps = Ups - vs * vps / (mu * alpha) + (mu / alpha) * Veffs_nonlocal_integrand

    # calculate van Kampen effective nonlocal potential Phi
    Phis_integrand = Veffps / Teffs
    Phis = np.concatenate((stub, np.cumsum((xs[1:] - xs[:-1]) * (Phis_integrand[1:] + Phis_integrand[:-1]) / 2.)),
                          axis=0)

    return Phis

def J_exact(xs, Us, Ups, Upps, vs, vps, mu, alpha):
    '''
    Calculate the exact current in a noninteracting system with
    potential U(x) and activity landscape v(x), using the formula
    from van Kampen 1988.
    '''
    stub = np.array([0])

    # calculate effective fields given U, x, mu, alpha
    # these fields are s.t. active system maps onto inhomogenous passive system
    Teffs = vs ** 2 / (mu * alpha) - mu * Ups ** 2 / alpha
    mueffs = mu / (1 - mu * Upps / alpha + mu * Ups * vps / (vs * alpha))
    Veffs = Us - vs ** 2 / (2 * mu * alpha)
    Veffs_nonlocal_integrand = (Ups ** 2) * vps / vs
    Veffs_nonlocal = np.concatenate(
        (stub, np.cumsum((xs[1:] - xs[:-1]) * (Veffs_nonlocal_integrand[1:] + Veffs_nonlocal_integrand[:-1]) / 2.)),
        axis=0)
    Veffs = Veffs + (mu / alpha) * Veffs_nonlocal
    Veffps = Ups - vs * vps / (mu * alpha) + (mu / alpha) * Veffs_nonlocal_integrand

    # calculate van Kampen effective nonlocal potential Phi
    Phis_integrand = Veffps / Teffs
    Phis = np.concatenate((stub, np.cumsum((xs[1:] - xs[:-1]) * (Phis_integrand[1:] + Phis_integrand[:-1]) / 2.)),
                          axis=0)

    # calculate exact expression
    inner_integrand = np.exp(Phis) / mueffs
    inner_integral = np.concatenate(
        (stub, np.cumsum((xs[1:] - xs[:-1]) * (inner_integrand[1:] + inner_integrand[:-1]) / 2.)), axis=0)

    denom_integrand = (np.exp(-Phis) / Teffs) * (
                (inner_integral[-1] - inner_integral) + inner_integral * np.exp(Phis[-1]))
    denom = np.trapezoid(denom_integrand, xs)

    J_exact = (1 - np.exp(Phis[-1])) / denom
    return J_exact

J = J_exact(xs, Us, Ups, Upps, vs, vps, 1, p.alpha)

Teffs = Teff(xs,Ups,vs,1,p.alpha)
mueffs= mueff(xs,Ups,Upps,vs,vps,1,p.alpha)
Veffs = Veff(xs,Us,Ups,vs,vps,1,p.alpha)
Phis  = Phi(xs,Us,Ups,Upps,vs,vps,1,p.alpha)

rho0_T0 = -J*np.sum(np.exp(Phis)/mueffs)*dx / (np.exp(Phis[-1]) - 1)
rho_exact = (np.exp(-Phis)/Teffs)*(rho0_T0 - J*np.cumsum(np.exp(Phis)/mueffs)*dx)
m_exact = (J + rho_exact*Ups)/vs



##########################################################################################################################################
# PLOTTING

# create figure/subplots
fig=plt.figure(figsize=(3.34,3.34*1.3),dpi=200)
gs = fig.add_gridspec(7, 1, hspace=0, wspace=0, height_ratios = [1,.15,1,.15,1.5,0.3,1.2])
ax0 = fig.add_subplot(gs[0,:])
ax1 = fig.add_subplot(gs[2,:],sharex=ax0)
ax2 = fig.add_subplot(gs[4,:],sharex=ax1)
ax3 = fig.add_subplot(gs[6,:])
axs = [ax0,ax1,ax2,ax3]

# plot activity/potential landscape
axs[0].plot(xs+dx/2., vs, color='green',label=r'$v(x)$')
axs[0].plot(xs+dx/2., Us+1, color='gray',label=r'$U(x)$')
axs[0].text(1,1.25,r'$v(x)$',color='green')
axs[0].text(1,0.7,r'$U(x)$',color='gray')
axs[0].set_ylim(ymax=np.max(vs)*1.4)

# plot effective temperature and pseudopotential
axs[1].plot(xs+dx/2., Teffs, color=cm(0.), label=r'$T_{\rm eff}(x)$')
axs[1].plot(xs+dx/2., Phis, color=cm(0.4), label=r'$\Phi(x)$')
axs[1].set_xlabel(r'$x$')
axs[1].text(1,0.65,r'$T_{\rm eff}(x)$',color=cm(0))
axs[1].text(1,-.1,r'$\Phi(x)$',color=cm(0.4))
axs[1].set_ylim(ymax=np.max(Teffs)*1.4)

# plot density and magnetization, simulation & theory
axs[2].scatter(prof[:,0]+dx_data/2.0,prof[:,1]/200,color=cm(0.7),lw=0,marker='o',s=20,label=r'$\rho_{\rm s}(x)$, data')
axs[2].plot(xs+dx/2.0, rho_exact,color='black',ls='--')
axs[2].scatter(m[:,0]+dx_data/2.0,m[:,1]/200,color=cm(0.9),lw=0,marker='o',s=20,label=r'$m_{\rm s}(x)$, data')
axs[2].plot(xs+dx/2.0, m_exact,color='black',ls='--',label='theory')
axs[2].set_xlabel(r'$x$')
axs[2].xaxis.set_label_coords(0.9,-0.05)
axs[2].legend(frameon=False,loc='upper left',bbox_to_anchor=(0.06,1))

for ax in axs[0:3]:
    ax.set_xlim(0,p.L)

# plot displacement vs. time, simulation & theory
axs[3].errorbar(disps2[::10,0], disps2[::10,1], yerr=disps2[::10,2]*3, color=cm(0.4), lw=0, marker='o', label=r'data')
axs[3].plot(disps2[:,0], disps2[:,0]*J, color='black', ls='--', label=r'theory')
axs[3].set_ylabel(r'$\frac{\langle \Delta x(t)\rangle}{L}$',rotation='horizontal')
axs[3].set_xlabel(r'$t$')
axs[3].xaxis.set_label_coords(0.9,-0.05)
axs[3].yaxis.set_label_coords(-0.07,0.85)
axs[3].set_xlim(0,p.final_time)
axs[3].set_ylim(-2,0.75)
axs[3].legend(frameon=False)

for ax in axs:
    ax.tick_params(axis='both',which='both',direction='in',length=2,pad=2)
for ax in axs[:2]:
    ax.tick_params(axis='x', which='both', labelbottom=False)
axs[2].set_xticks([0,5,10,15,20])


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
    print(pad_x,pad_y)
    text = ax.text(pad_x,1-pad_y,labels[i],ha='left',va='top',transform=ax.transAxes,
                       bbox=dict(boxstyle='square,pad='+str(inside_pad),edgecolor='black',facecolor='white',lw=edgewidth))
    text.set_clip_on(True)

plt.subplots_adjust(bottom=0.04,top=0.99,right=0.96)

plt.savefig('disorder-U-v.pdf',dpi=1000)