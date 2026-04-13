import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os,sys

sys.path.append(os.getcwd())
import read_data as rd
import mydefaults as md

plt.rcParams.update(md.plt_rcparams)


################################################################################################################
# READ DATA

alphas = ['1.00','2.00','4.00','8.00']
currents_all = []
currents_all_th = []
for alpha in alphas:
    currents_all += [[float(alpha),np.genfromtxt('current_'+alpha)]]
    currents_all_th += [[float(alpha),np.genfromtxt('current_th_'+alpha)]]
        

################################################################################################################
# ACTIVITY AND POTENTIAL LANDSCAPES

# get parameters of one sample simulation
p = rd.get_params('param')

# trigonometric v and U fields
def v(x,v0,a1,a2,d, P,L):
    return v0 * (1 + .5 * a1 * np.sin(2 * np.pi * x*P / L + 2*np.pi*d) + .5 * a2 * np.sin(4 * np.pi * x*P / L + 4*np.pi*d))
def vp(x,v0,a1,a2,d, P,L):
    return v0 * (.5*a1*2*P*np.pi*(1/L)*np.cos(2*np.pi*x*P/L + d*2*np.pi) + .5*a2*4*P*np.pi*(1/L)*np.cos(4*np.pi*x*P/L + d*4*np.pi))    
def U(x,epsilon,a1,a2,b1,b2,P,L):
    return epsilon*( b1 * np.sin(2 * np.pi * x*P / L) +  b2 * np.sin(4 * np.pi * x*P / L))
def Up(x,epsilon,a1,a2,b1,b2,P,L):
    return epsilon*( b1 * (2*P*np.pi/L)*np.cos(2 * np.pi * x*P / L) + b2 * (4*P*np.pi/L)*np.cos(4 * np.pi * x*P / L))
def Upp(x,epsilon,a1,a2,b1,b2,P,L):
    return epsilon*(- b1 * (2*P*np.pi/L)*(2*P*np.pi/L)*np.sin(2 * np.pi * x*P / L) - b2 * (4*P*np.pi/L)*(4*P*np.pi/L)*np.sin(4 * np.pi * x*P / L))


################################################################################################################
# EFFECTIVE TEMPERATURE, PSEUDOPOTENTIAL, ETC. FIELDS


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

# construct effective mu, T, V, etc. profiles
dx = p.L/(5000.*p.P)
xs = np.arange(0, p.L+dx, dx)
vs = v(xs, p.v0, p.a1, p.a2, p.d, p.P, p.L)
vps = vp(xs, p.v0, p.a1, p.a2, p.d, p.P, p.L)
Us = U(xs,p.epsilon,p.a1,p.a2,p.b1,p.b2,p.P,p.L)
Ups = Up(xs,p.epsilon,p.a1,p.a2,p.b1,p.b2,p.P,p.L)
Upps = Upp(xs,p.epsilon,p.a1,p.a2,p.b1,p.b2,p.P,p.L)
J_exact = J_exact(xs, Us, Ups, Upps, vs, vps, 1, p.alpha)
Phi0s = Phi(xs, Us*0, Ups*0, Upps*0, vs, vps, 1, p.alpha)
Phis = Phi(xs, Us, Ups, Upps, vs, vps, 1, p.alpha)
mueffs = mueff(xs, Ups, Upps, vs, vps, 1, p.alpha)
Teffs = Teff(xs, Ups, vs, 1, p.alpha)


################################################################################################################
# PLOT

cm=plt.get_cmap('viridis')
colors=[cm(i/float(len(alphas)-1)) for i in range(len(alphas))]

fig = plt.figure(figsize=(8.6/2.54, 8.6*.4/2.54))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], hspace=0, wspace=0.05)

# Create the subplots on the left side with shared x-axis
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

# Create the subplot on the right side
ax4 = fig.add_subplot(gs[:, 1])

# plot currents (simulation)
for i,currents in enumerate(currents_all):
    alpha,currents = currents
    ax4.errorbar(currents[:,0],currents[:,1],color=colors[i], mec='black',mew=.75,
    			label=r'$\tau^{-1}$='+str(int(alpha)),fmt='o',markersize=4)
ylim=plt.ylim()

# plot currents (theory)
for i,currents in enumerate(currents_all_th):
    alpha,currents = currents
    ax4.plot(currents[:,0],currents[:,1],color=colors[i],lw=2)
ax4.plot([-100,-101],[0,0],ls='-',color='black',lw=1.85)#,label='theory')
ax4.set_ylim(*ylim)
ax4.set_ylim(ymax=ylim[1]*1.05)
ax4.set_xlim(xmin=0,xmax=.9)
ax4.set_xlabel(r'$\varepsilon$')
ax4.set_ylabel(r'$J$',rotation='horizontal')
ax4.set_ylim(ymin=0)
ax4.legend(frameon=False,loc=(-0.04,.29),handletextpad=-0.1)

# Move y-axis to the right, set ticks/labels
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.set_xticks([0,0.25,0.5,0.75],['0','0.25','0.5','0.75'])
ax4.set_yticks([0,0.003,0.006],['0','0.003','0.006'])
ax4.yaxis.set_label_coords(1.05,.95)
ax4.xaxis.set_label_coords(0.97,-0.01)

# plot v(x) and V(x)
ax1.plot(xs+dx/2., vs, color='#0d800f',lw=1.25,label=r'$v(x)$')
ax1.plot(xs+dx/2., Us*0, color='gray',lw=1,ls='--')
ax1.plot(xs+dx/2., Us*1.5, color='purple',lw=1.25,label=r'$U(x)$')
adjust=0.1 # the arrow isn't long enough for some reason
ax1.annotate('', (6.25,-adjust), xytext=(6.25,np.max(Us)*1.5+adjust), arrowprops=dict(arrowstyle='<->', lw=0.75, color='black', shrinkA=0, shrinkB=0, mutation_scale=3,connectionstyle='arc3'))
ax1.text(7.5,np.max(Us)*.5,r'$\varepsilon$')
ax1.set_yticks([0])
ax1.text(p.L*0.025,vs[0]*0.35,r'$V(x)$',color='black')
ax1.text(p.L*0.15,vs[0]*1.15,r'$v(x)$',color='green')
ax1.yaxis.set_label_coords(-0.1,0.25)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

# plot pseudo-potential
ax2.plot(xs+dx/2., Phi0s, color='gray', lw=1, ls='--', label=r'$\Phi_{U=0}$')
ax2.plot(xs+dx/2., Phis, color='purple', lw=1.25, label=r'$\Phi$')
ax2.arrow(10,np.max(Phi0s)*0.77,11,0,
          color='purple',head_width=0.1,head_length=0.75)
ax2.text(13,np.max(Phi0s)*0.35,r'$J$',color='purple')
ax2.set_yticks([0])
ax2.text(p.L*0.025,np.max(Phi0s)*0.5,r'$\Phi(x)$',color='black')
ax2.yaxis.set_label_coords(-0.1,0.16)
ax2.set_xlabel(r'$x$')
ax2.xaxis.set_label_coords(0.95,-0.04)
ax2.set_xlim(xmin=0,xmax=p.L)
ax2.set_xticks([0,10,20],labels=[r'$0$',r'$L_x$',r'$2L_x$'])

# adjust axis ticks
for ax in [ax1,ax2,ax4]:
    ax.tick_params(axis='both',which='both',direction='in',pad=2)

plt.subplots_adjust(bottom=0.11,left=0.04,right=.91,top=.97)


plt.savefig('noninteract_simple.pdf', dpi=1000)
