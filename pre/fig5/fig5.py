import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import copy

# matplotlib options
plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
params = {'text.usetex' : True,
          'font.size' : 8,
          'font.family' : 'lmodern'
          }
plt.rcParams.update(params)


#########################################################################################################################################
# IMPORT DATA - FUNCTIONS
    
# for reading param info
def get(paramstr, param, keepstr=False):
    s = paramstr.split(param+' is ')[-1].split('\n')[0]
    if keepstr:
        return s
    return float(s)

# read in simulation data (abp and rtp, avg over multiple seeds)
def get_prof_disp(supdir):
    folders = [x for x in os.listdir(supdir) if (x.startswith('abp_') or x.startswith('rtp_')) and os.path.isdir(supdir+'/'+x)]

    # amp = interaction strength
    amp = float(supdir.split('amp_')[-1].split('-')[0])

    # where to leave data
    prof_avgs = {'abp': None, 'rtp': None}
    disp_avgs = {'abp': None, 'rtp': None}
    profs = {'abp': [], 'rtp': []}
    disps = {'abp': [], 'rtp': []}

    # iterate through seeds (abp and rtp)
    for i,folder in enumerate(folders):
        path = supdir+'/'+folder+'/res-'

        # which type (abp or rtp)
        if 'abp' in folder:
            flag = 'abp'
        elif 'rtp' in folder:
            flag = 'rtp'

        # get param info
        f=open(path+'param')
        param=f.read()
        f.close()
        Lx=get(param,'Lx')
        Ly=get(param,'Ly')
        N=int(get(param,'N'))
        Nbinx=int(get(param,'Nbinx'))

        # import data
        disp = np.genfromtxt(path+'disp',delimiter='\t')
        prof = np.genfromtxt(path+'profile_rho',delimiter='\t')

        # how many time points in profile data
        Ts = int(prof.shape[0]/Nbinx+1e-5)

        # averaged density profile
        prof_avg = np.zeros((Nbinx,2))
        prof_avg[:,0] = prof[0:Nbinx,1]
        prof_avg[:,1] = sum([prof[Nbinx*i:Nbinx*(i+1),3] for i in range(Ts)])/Ts

        # save avg profile to aggregated data
        profs[flag] += [prof_avg]
        if prof_avgs[flag] is None:
            prof_avgs[flag] = prof_avg
        else:
            prof_avgs[flag][:,1] += prof_avg[:,1]

        # how many time points in displacement data
        Ts = int(disp.shape[0]/N+1e-5)

        # averaged displacement
        disp_avg = np.zeros((Ts,3)) # columns: time, avg displacement, stderr displacement
        disp_avg[:,0] = disp[::N,0]
        disp_avg[:,1] = sum([disp[i::N,2] for i in range(N)])/(N*Lx)

        # save avg displacement to aggregated data
        disps[flag] += [disp_avg.copy()]
        if disp_avgs[flag] is None:
            disp_avgs[flag] = disp_avg.copy()
        else:
            disp_avgs[flag][:,1] += disp_avg.copy()[:,1]
            

    # average aggregated data over different seeds
    for flag in ['abp','rtp']:

        # divide average by number of seeds
        disp_avgs[flag][:,1] = disp_avgs[flag][:,1]/len(disps[flag])

        # concatenate zero to beginning of displacement arrays
        disp_avgs[flag] = np.concatenate((np.zeros((1,3)),disp_avgs[flag]),axis=0)
        for i in range(len(disps[flag])):
            disps[flag][i] = np.concatenate((np.zeros((1,3)),disps[flag][i]),axis=0)

        # calculate standard error in last column
        disp_avgs[flag][:,2] = np.array([np.std([disps[flag][i][j,1] for i in range(len(disps[flag]))],ddof=1) / np.sqrt(len(disps[flag])) for j in range(disps[flag][0].shape[0])])

        # divide by integral to normalize
        prof_avgs[flag][:,1] = prof_avgs[flag][:,1] / np.trapezoid(prof_avgs[flag][:,1], prof_avgs[flag][:,0])
        for i in range(len(profs[flag])):
            profs[flag][i][:,1] = profs[flag][i][:,1] / np.trapezoid(profs[flag][i][:,1], profs[flag][i][:,0])

    # together
    disps_tog = {flag: np.array(disps[flag])[:,:,1] for flag in ['abp','rtp']}
            
    return profs,prof_avgs,disps,disp_avgs,disps_tog




#########################################################################################################################################
# IMPORT DATA - EXAMPLE SNAPSHOT/DENSITIES/DISPLACEMENTS (PANELS B-D)

# which directory to use as example snapshot/density profile/displacement (for panels b-d)
supdir = 'data/N_500-amp_0.075-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12'
profs,prof_avgs,disps,disp_avgs,disps_tog = get_prof_disp(supdir)
disp_stderr = {flag: np.std(disps_tog[flag],axis=0,ddof=1)/np.sqrt(len(disps[flag])) for flag in ['abp','rtp']}

# get parameter info
sds = [x for x in os.listdir(supdir) if x.startswith('abp_') and os.path.isdir(os.path.join(supdir,x))]
f=open(os.path.join(supdir,sds[0]+'/res-param'))
param=f.read()
f.close()
Lx=get(param,'Lx')
N=int(get(param,'N'))
Ly=get(param,'Ly')
final_time=get(param,'final_time')

# get particle location snapshot for that run
pos = np.genfromtxt(os.path.join(supdir,sds[0]+'/res-pos'),delimiter='\t')
pos = pos[-N:,:]

# which directory to use for non-interacting example density profile/displacement (panels c-d)
supdir = 'data/N_500-amp_0-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12'
profs_NI,prof_avgs_NI,disps_NI,disp_avgs_NI,disps_tog_NI=get_prof_disp(supdir)
disp_stderr_NI = {flag: np.std(disps_tog_NI[flag],axis=0,ddof=1)/np.sqrt(len(disps_NI[flag])) for flag in ['abp','rtp']}

# get array of x points
xs_data = prof_avgs['abp'][:,0]

# shift x axis for display
x0 = 18.5
inds1 = np.where(xs_data>x0)
inds2 = np.where(xs_data<=x0)

# shift density profiles
prof_avgs_ = prof_avgs.copy()
prof_avgs_NI_ = prof_avgs_NI.copy()
for flag in ['abp','rtp']:
    prof_avgs_[flag][:,1] = np.concatenate([prof_avgs[flag].copy()[inds1[0],1], prof_avgs[flag].copy()[inds2[0],1]])
    prof_avgs_NI_[flag][:,1] = np.concatenate([prof_avgs_NI[flag].copy()[inds1[0],1], prof_avgs_NI[flag].copy()[inds2[0],1]])
prof_avgs = prof_avgs_ 
prof_avgs_NI = prof_avgs_NI_

# shift particle positions
pos[:,2] = np.remainder(pos[:,2]-x0, Lx)



#########################################################################################################################################
# IMPORT DATA - AVG CURRENT FOR EACH INTERACTION STRENGTH (PANEL E)

# list of interaction strengths
amps = [0, 0.025, 0.05,0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225]

Js_arr = []
for amp in amps:

    # directory name for that interaction strength
    supdir = 'data/N_500-amp_'+str(amp)+'-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12'

    # list of seeds within that directory (both abp and rtp)
    folders = [x for x in os.listdir(supdir) if (x.startswith('abp_') or x.startswith('rtp_')) \
               and os.path.isdir(supdir+'/'+x) and os.path.exists(supdir+'/'+x+'/res-param')]

    # dict to store current data
    Js = {'abp': [], 'rtp': []}

    # iterate through seeds
    n_seeds_abp = 0
    n_seeds_rtp = 0
    for i,folder in enumerate(folders):
        path = supdir+'/'+folder+'/res-'

        # which system (abp or rtp)
        if 'abp' in folder:
            flag = 'abp'

            # the example data for amp=0,0.075 (panels b-d) has extra seeds. 
            # for panel (e), use same number of seeds (100) as other interaction strengths.
            n_seeds_abp += 1
            if n_seeds_abp > 100:
                continue

        elif 'rtp' in folder:
            flag = 'rtp'

            # the example data for amp=0,0.075 (panels b-d) has extra seeds. 
            # for panel (e), use same number of seeds (100) as other interaction strengths.
            n_seeds_rtp += 1
            if n_seeds_rtp > 100:
                continue

        # get param info
        f=open(path+'param')
        param=f.read()
        f.close()
        Lx=get(param,'Lx')
        N=int(get(param,'N'))

        # read in displacement data
        disp = np.genfromtxt(path+'disp',delimiter='\t')

        # average displacement for that seed
        Js[flag] += [np.mean(disp[-N:,2])/(Lx*disp[-1,0])]
    
    # aggregate data (mean and std)
    Js_mean = {'abp': None, 'rtp': None}
    Js_std = {'abp': None, 'rtp': None}
    for flag in ['abp','rtp']:
        Js[flag] = np.array(Js[flag])
        Js_mean[flag] = np.mean(Js[flag])
        Js_std[flag] = np.std(Js[flag],ddof=1)/np.sqrt(Js[flag].shape[0])
    Js_arr += [[amp, Js_mean['abp'], Js_mean['rtp'], Js_std['abp'], Js_std['rtp']]]
    
disp = np.array(Js_arr)

# theoretical current expansion coefficients (Eq. 98 of paper)
Jns=[0.0, -3.0876260512230874e-05]


#########################################################################################################################################
# ACTIVITY LANDSCAPE

def v_arr(xs,p):
    '''
    Piecewise cubic activity landscape
    '''
    
    v0, v1, vc1, vc2, vr, Lx = p
    dv = v1-v0
    dxR = vr-vc2
    
    c1    = 3*dv / pow(vc1,2)
    d1    = -2*dv / pow(vc1,3)
  
    a2    = -dv * pow(vc2,2) * (3*dxR + 2*vc2) / pow(dxR,3) + v1
    b2    = 6*dv * vc2 * (dxR + vc2) / pow(dxR,3)
    c2    = -3*dv * (dxR + 2*vc2) / pow(dxR,3)
    d2    = 2*dv / pow(dxR,3)
    
    vs = np.zeros(xs.shape)
    vs[xs>=vr] = v0
    vs[xs<vc1] = v0 + c1*pow(xs[xs<vc1],2) + d1*pow(xs[xs<vc1],3)
    vs[(xs<=vc2) & (xs>=vc1)] = v1
    vs[(xs>vc2) & (xs<vr)] = a2 + b2*xs[(xs>vc2) & (xs<vr)] + c2*pow(xs[(xs>vc2) & (xs<vr)],2) + d2*pow(xs[(xs>vc2) & (xs<vr)],3)
    return vs

# get activity landscape
vc1 = get(param,'v_center1')
vc2 = get(param,'v_center2')
vr  = get(param,'v_right')
v0  = get(param,'v0')
v1  = get(param,'v1')
p = v0, v1, vc1, vc2, vr, Lx
vs_data = v_arr(xs_data, p)

# shift v(x)
vs_data = np.concatenate([vs_data[inds1],vs_data[inds2]])


#########################################################################################################################################
# PLOTTING

# create figure and sub-panels
fig = plt.figure(figsize=(7,7*0.3),dpi=600)
gs = fig.add_gridspec(32, 5, hspace=0, wspace=0,width_ratios=[1,0.1,1,0.1,1])
ax2 = fig.add_subplot(gs[19:,0])
ax0 = fig.add_subplot(gs[0:8,0])
ax1 = fig.add_subplot(gs[10:17,0])
ax3 = fig.add_subplot(gs[:,2])
ax4 = fig.add_subplot(gs[:,4])
axs = [ax0,ax1,ax2,ax3,ax4]

####################################################
# PANEL A - activity landscape

ax0.plot(xs_data,vs_data,color='green')

ax0.set_yticks([0,10,20])
ax0.set_ylabel(r'$v(x)$',rotation='horizontal')
ax0.yaxis.set_label_coords(-0.062,0.7)
ax0.set_ylim(ymax=38)

####################################################
# PANEL B - simulation snapshot

# get window extent bbox and convert from inches to display coords
xlim = ax1.get_xlim()
ylim = ax1.get_ylim()
display_coords = ax1.transData.transform([(xlim[0],ylim[0]), (xlim[1], ylim[1])])
w,h = display_coords[1]-display_coords[0] # width, height of axis in PIXELS (will not change)
aspect=w/h
ax1.set_xlim(xmin=0,xmax=Lx)
ax1.set_ylim(ymin=0,ymax=Lx*h/w)

# display simulation shapshot of particles (make sure size is correct)
ppp = plt.gcf().dpi / 72 # px per inch * inch per pt = px per pt
s = (np.pi/4.)*((w/Lx)/ppp)**2
ax1.scatter(pos[:,2],pos[:,3],color='orange',s=s,linewidths=0,alpha=0.5) #edgecolors='#ed6d05', linewidths=0

# plot arrows for self-propulsion
vs_pos = v_arr(np.remainder(pos[:,2]+x0, Lx),p)
f=0.05
ax1.quiver(pos[:,2],pos[:,3],vs_pos*np.cos(pos[:,4])*f,vs_pos*np.sin(pos[:,4])*f,
           color='black',headwidth=4.5)

ax1.set_ylabel(r'$y$',rotation='horizontal')
ax1.yaxis.set_label_coords(-0.03,0.55)
ax1.set_yticks([])

####################################################
# PANEL C - density profile

# plot density profile
ax2.plot(prof_avgs_NI['abp'][:,0], prof_avgs_NI['abp'][:,1], color='black', ls='--', label='NI')
ax2.plot(prof_avgs['abp'][:,0], prof_avgs['abp'][:,1], color='orange', label='ABP')
ax2.plot(prof_avgs['rtp'][:,0], prof_avgs['rtp'][:,1], color='purple', label='RTP')

ax2.set_xlabel(r'$x$')
ax2.set_xticks([0,5,10,15])
ax2.set_yticks([0,0.1])
ax2.xaxis.set_label_coords(0.92,-0.03)
ax2.set_ylabel(r'$\rho_{\rm s}(x)$',rotation='horizontal')
ax2.yaxis.set_label_coords(-0.063,0.75)
ax2.legend(frameon=False,loc='upper right')

####################################################
# PANEL D - particle displacement vs. time

# delta x = 0 line
ax3.plot([-100000,-100000],[-1,0],markersize=4,marker='o',color='black',ls='--',markerfacecolor='white',label='NI')

# non-interacting simulation displacement vs. time
ax3.plot(disp_avgs_NI['abp'][:,0], disp_avgs_NI['abp'][:,1],ls='--',color='orange',zorder=2, 
         marker='o',markersize=4,markerfacecolor='white')
ax3.plot(disp_avgs_NI['rtp'][:,0], disp_avgs_NI['rtp'][:,1],ls='--',color='purple',zorder=2, 
         marker='o',markersize=4,markerfacecolor='white')
ax3.fill_between(disp_avgs_NI['abp'][:,0],
                 disp_avgs_NI['abp'][:,1]+disp_stderr_NI['abp']*3,
                 disp_avgs_NI['abp'][:,1]-disp_stderr_NI['abp']*3,
                 zorder=0, alpha=0.25, color='orange',ls='--')
ax3.fill_between(disp_avgs_NI['rtp'][:,0],
                 disp_avgs_NI['rtp'][:,1]+disp_stderr_NI['rtp']*3,
                 disp_avgs_NI['rtp'][:,1]-disp_stderr_NI['rtp']*3,
                 zorder=0, alpha=0.25, color='purple',ls='--')

# interacting simulation displacement vs. time
ax3.plot(disp_avgs['abp'][:,0], disp_avgs['abp'][:,1],label='ABP',color='orange',zorder=2, marker='o',markersize=4)
ax3.plot(disp_avgs['rtp'][:,0], disp_avgs['rtp'][:,1],label='RTP',color='purple',zorder=2, marker='o',markersize=4)
ax3.fill_between(disp_avgs['abp'][:,0],
                 disp_avgs['abp'][:,1]+disp_stderr['abp']*3,
                 disp_avgs['abp'][:,1]-disp_stderr['abp']*3,
                 color='orange', zorder=1, alpha=0.25)
ax3.fill_between(disp_avgs['rtp'][:,0],
                 disp_avgs['rtp'][:,1]+disp_stderr['rtp']*3,
                 disp_avgs['rtp'][:,1]-disp_stderr['rtp']*3,
                 color='purple', zorder=1, alpha=0.25)

ax3.set_xticks([0,5000,10000,15000])
ax3.set_xlabel(r'$t$')
ax3.xaxis.set_label_coords(0.95,-0.03)
ax3.set_ylabel(r'$\Delta x(t)/L_x$',rotation='horizontal')
ax3.yaxis.set_label_coords(.15,0.35)
ax3.set_xlim(xmin=0,xmax=final_time)
ax3.set_ylim(ymax=1.2)
ax3.tick_params(axis='y', which='both', pad=1)
ax3.legend(frameon=False,loc='lower center',bbox_to_anchor=(0.2,0))

####################################################
# PANEL E - current vs. interaction strength

# simulation data
ax4.errorbar(disp[:,0]*N/Ly, disp[:,1], yerr=disp[:,3]*3, color='orange', label='ABP', marker='o',markersize=5)
ax4.errorbar(disp[:,0]*N/Ly, disp[:,2], yerr=disp[:,4]*3, color='purple', label='RTP', marker='o',markersize=5)

# simulation data - non-interacting point
ax4.errorbar(disp[:1,0]*N/Ly,disp[:1,1],yerr=disp[:1,3]*3,color='orange', marker='o', markersize=5, markerfacecolor='white')
ax4.errorbar(disp[:1,0]*N/Ly,disp[:1,2],yerr=disp[:1,4]*3,color='purple', marker='o', markersize=5, markerfacecolor='white')

# plot theoretical current (just linear part)
xlim = ax4.get_xlim()
ylim = ax4.get_ylim()
epsilons = np.linspace(0,disp[-1,0]*N/Ly,100)
plt.plot(epsilons, sum([epsilons**n * Jns[n] for n in range(len(Jns))]),
         color='black',label=r'$J_{1} N \varepsilon/L_y$')
ax4.plot([xlim[0]-10,xlim[1]+10],[0,0],color='black',lw=0.5)
ax4.set_ylim(*ylim)
ax4.set_xlim(*xlim)

ax4.set_xticks([0,5,10,15,20])
ymin = np.min(disp[:,1])*2
ymax = np.max(disp[:,2])*1.1
ax4.set_yticks([-0.0004,0,0.0004,0.0008],labels=['$-4$', '0', '4','8'])
ax4.text(20.7,ymin*0.8,r'$\times 10^{-4}$')
ax4.set_xlabel(r'$N \varepsilon / L_y$')
ax4.xaxis.set_label_coords(0.96,0.-0.02)
ax4.set_ylabel(r'$J$',rotation='horizontal')
ax4.yaxis.set_label_coords(1.05,0.8)
ax4.set_ylim(ymin=ymin,ymax = ymax)
ax4.set_xlim(xmax=25)
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.legend(frameon=False,loc='center left',bbox_to_anchor=(0,0.7))



####################################################
# OTHER PLOTTING PARAMETERS

# axis ticks
for ax in axs:
    ax.tick_params(axis='both',which='both',direction='in',pad=2)
for ax in axs[:2]:
    ax.set_xticks([0,5,10,15],labels=['']*4)

# axis limits
for ax in axs[:3]:
    ax.set_xlim(xmin=0,xmax=Lx)
    ax.set_ylim(ymin=0)


# add (a), (b), etc.
edgewidth = 0.5
pad = 1.25
fontsize=8.0
labels=['(a)','(b)','(c)','(d)','(e)']
for i,ax in enumerate(axs):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    pad_x = pad / 72 / width
    pad_y = pad / 72 / height
    inside_pad = pad/fontsize
    text = ax.text(pad_x,1-pad_y,labels[i],ha='left',va='top',transform=ax.transAxes,
                       bbox=dict(boxstyle='square,pad='+str(inside_pad),edgecolor='black',facecolor='white',lw=edgewidth))
    text.set_clip_on(True)

# figure geometry
plt.subplots_adjust(hspace=0,left=0.039,bottom=0.088,right=0.975,top=0.99)
    
plt.savefig('weak_rtp_abp_2d_current.pdf',dpi=400)
