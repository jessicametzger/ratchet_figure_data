import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from tqdm import tqdm
import matplotlib.gridspec as gridspec
from IPython.display import Image,display,HTML
import copy

#Direct input 
plt.rcParams['text.latex.preamble']=r"\usepackage{lmodern}"
#Options
params = {'text.usetex' : True,
          'font.size' : 8,
          'font.family' : 'lmodern'
          }
plt.rcParams.update(params)

def v_arr(xs,p):
    '''
    Piecewise cubic
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
        
    
def get(paramstr, param, keepstr=False):
    s = paramstr.split(param+' is ')[-1].split('\n')[0]
    if keepstr:
        return s
    return float(s)


def data_to_marker_size(data_size, ax):

    display_coords = ax.transData.transform([(0, 0), (1, 1)])
    data_coords = ax.transAxes.transform([(0, 0), (1, 1)])

    # Calculate the scale between data and display coordinates (numpy division -> (dx1/dx2,dy1/dy2))
    scale = (display_coords[1] - display_coords[0]) / (data_coords[1] - data_coords[0])

    # Apply the scale to the marker size in the data coordinate
    # assume that aspect ratio of ax has been made equal
    marker_size_display = data_size / (scale[1] * scale[0])

    # Convert from pixels to points (1 point = 1/72 inch)
    # And then to points^2
    ppp = plt.gcf().dpi / 72 # px per inch * inch per pt = px per pt
    marker_size_points = marker_size_display / ppp**2 # px squared
    
    return marker_size_points



def get_prof_disp(supdir):
    folders = [x for x in os.listdir(supdir) if (x.startswith('abp_') or x.startswith('rtp_')) and os.path.isdir(supdir+'/'+x)]

    amp = float(supdir.split('amp_')[-1].split('-')[0])

    prof_avgs = {'abp': None, 'rtp': None}
    disp_avgs = {'abp': None, 'rtp': None}
    profs = {'abp': [], 'rtp': []}
    disps = {'abp': [], 'rtp': []}
    for i,folder in tqdm(enumerate(folders)):
        path = supdir+'/'+folder+'/res-'
        if 'abp' in folder:
            flag = 'abp'
        elif 'rtp' in folder:
            flag = 'rtp'

        if not os.path.exists(path+'param'):
            continue

        f=open(path+'param')
        param=f.read()
        f.close()

        Lx=get(param,'Lx')
        Ly=get(param,'Ly')
        N=int(get(param,'N'))
        Nbinx=int(get(param,'Nbinx'))

        disp = np.genfromtxt(path+'disp',delimiter='\t')
        if disp.shape[0]<10*N:
            continue
        prof = np.genfromtxt(path+'profile_rho',delimiter='\t')

        Ts = int(prof.shape[0]/Nbinx+1e-5)
        prof_avg = np.zeros((Nbinx,2))
        prof_avg[:,0] = prof[0:Nbinx,1]
        prof_avg[:,1] = sum([prof[Nbinx*i:Nbinx*(i+1),3] for i in range(Ts)])/Ts

        profs[flag] += [prof_avg]
        if prof_avgs[flag] is None:
            prof_avgs[flag] = prof_avg
        else:
            prof_avgs[flag][:,1] += prof_avg[:,1]

        Ts = int(disp.shape[0]/N+1e-5)
        disp_avg = np.zeros((Ts,3))
        disp_avg[:,0] = disp[::N,0]
        disp_avg[:,1] = sum([disp[i::N,2] for i in range(N)])/(N*Lx)

        disps[flag] += [disp_avg.copy()]
        if disp_avgs[flag] is None:
            disp_avgs[flag] = disp_avg.copy()
        else:
            disp_avgs[flag][:,1] += disp_avg.copy()[:,1]
            

    for flag in ['abp','rtp']:
        disp_avgs[flag] = np.concatenate((np.zeros((1,3)),disp_avgs[flag]),axis=0)
        for i in range(len(disps[flag])):
            disps[flag][i] = np.concatenate((np.zeros((1,3)),disps[flag][i]),axis=0)            
        disp_avgs[flag][:,1] = disp_avgs[flag][:,1]/len(disps[flag])
        disp_avgs[flag][:,2] = np.array([np.std([disps[flag][i][j,1] for i in range(len(disps[flag]))],ddof=1) / np.sqrt(len(disps[flag])) for j in range(disps[flag][0].shape[0])])

        prof_avgs[flag][:,1] = prof_avgs[flag][:,1] / np.trapezoid(prof_avgs[flag][:,1], prof_avgs[flag][:,0])
        for i in range(len(profs[flag])):
            profs[flag][i][:,1] = profs[flag][i][:,1] / np.trapezoid(profs[flag][i][:,1], profs[flag][i][:,0])

    disps_tog = {flag: np.array(disps[flag])[:,:,1] for flag in ['abp','rtp']}
            
    return profs,prof_avgs,disps,disp_avgs,disps_tog


supdir = 'N_500-amp_0.075-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12'
profs,prof_avgs,disps,disp_avgs,disps_tog = get_prof_disp(supdir)
disp_stderr = {flag: np.std(disps_tog[flag],axis=0,ddof=1)/np.sqrt(len(disps[flag])) for flag in ['abp','rtp']}

np.savetxt(os.path.join(supdir,supdir+'_prof_avg_abp'),prof_avgs['abp'],delimiter='\t')
np.savetxt(os.path.join(supdir,supdir+'_prof_avg_rtp'),prof_avgs['rtp'],delimiter='\t')
np.savetxt(os.path.join(supdir,supdir+'_disp_avg_abp'),disp_avgs['abp'],delimiter='\t')
np.savetxt(os.path.join(supdir,supdir+'_disp_avg_rtp'),disp_avgs['rtp'],delimiter='\t')

sds = [x for x in os.listdir(supdir) if x.startswith('abp_') and os.path.isdir(os.path.join(supdir,x))]
f=open(os.path.join(supdir,sds[0]+'/res-param'))
param=f.read()
f.close()
Lx=get(param,'Lx')
N=int(get(param,'N'))
Ly=get(param,'Ly')
final_time=get(param,'final_time')

pos = np.genfromtxt(os.path.join(supdir,sds[0]+'/res-pos'),delimiter='\t')
pos = pos[-N:,:]

supdir = 'N_500-amp_0-Lx_20-Ly_5-v0_20-v1_1-landscape_2.5_5_12'
profs_NI,prof_avgs_NI,disps_NI,disp_avgs_NI,disps_tog_NI=get_prof_disp(supdir)
disp_stderr_NI = {flag: np.std(disps_tog_NI[flag],axis=0,ddof=1)/np.sqrt(len(disps_NI[flag])) for flag in ['abp','rtp']}

np.savetxt(os.path.join(supdir,supdir+'_prof_avg_abp'),prof_avgs_NI['abp'],delimiter='\t')
np.savetxt(os.path.join(supdir,supdir+'_prof_avg_rtp'),prof_avgs_NI['rtp'],delimiter='\t')
np.savetxt(os.path.join(supdir,supdir+'_disp_avg_abp'),disp_avgs_NI['abp'],delimiter='\t')
np.savetxt(os.path.join(supdir,supdir+'_disp_avg_rtp'),disp_avgs_NI['rtp'],delimiter='\t')


prof_avgs_ = prof_avgs.copy()
prof_avgs_NI_ = prof_avgs_NI.copy()
pos_ = pos.copy()



xs_data = prof_avgs['abp'][:,0]

vc1 = get(param,'v_center1')
vc2 = get(param,'v_center2')
vr  = get(param,'v_right')
v0  = get(param,'v0')
v1  = get(param,'v1')
p = v0, v1, vc1, vc2, vr, Lx

vs_data = v_arr(xs_data, p)

x0 = 18.5
inds1 = np.where(xs_data>x0)
inds2 = np.where(xs_data<=x0)

vs_data = np.concatenate([vs_data[inds1],vs_data[inds2]])

prof_avgs = {}
prof_avgs_NI = {}
for flag in ['abp','rtp']:
    prof_avgs[flag] = prof_avgs_[flag].copy()
    prof_avgs_NI[flag] = prof_avgs_NI_[flag].copy()
    prof_avgs[flag][:,1] = np.concatenate([prof_avgs_[flag].copy()[inds1[0],1], prof_avgs_[flag].copy()[inds2[0],1]])
    prof_avgs_NI[flag][:,1] = np.concatenate([prof_avgs_NI_[flag].copy()[inds1[0],1], prof_avgs_NI_[flag].copy()[inds2[0],1]])

pos = pos_.copy()
pos[:,2] = np.remainder(pos[:,2]-x0, Lx)



disp = np.genfromtxt('weak_disp_avg')

Jns=[0.0, -3.0876260512230874e-05]#, 2.9904413400135847e-08, 2.8860723957634118e-09, -1.011495607160336e-10, 2.8840815109664217e-12, -7.958150939476408e-14, 2.20549317977401e-15, -6.192579157336801e-17, 1.7621364434826948e-18, -5.071566087590012e-20, 1.4729224486014256e-21, -4.307947796110543e-23, 1.266741154139194e-24]



fig = plt.figure(figsize=(7,7*0.3),dpi=600)
gs = fig.add_gridspec(32, 5, hspace=0, wspace=0,width_ratios=[1,0.1,1,0.1,1])
ax2 = fig.add_subplot(gs[19:,0])
ax0 = fig.add_subplot(gs[0:8,0])
ax1 = fig.add_subplot(gs[10:17,0])
ax3 = fig.add_subplot(gs[:,2])
ax4 = fig.add_subplot(gs[:,4])
axs = [ax0,ax1,ax2,ax3,ax4]

ax0.plot(xs_data,vs_data,color='green')

# get window extent bbox and convert from inches to display coords
xlim = ax1.get_xlim()
ylim = ax1.get_ylim()
display_coords = ax1.transData.transform([(xlim[0],ylim[0]), (xlim[1], ylim[1])])
w,h = display_coords[1]-display_coords[0] # width, height of axis in PIXELS (will not change)
aspect=w/h
ax1.set_xlim(xmin=0,xmax=Lx)
ax1.set_ylim(ymin=0,ymax=Lx*h/w)

ppp = plt.gcf().dpi / 72 # px per inch * inch per pt = px per pt
s = (np.pi/4.)*((w/Lx)/ppp)**2
ax1.scatter(pos[:,2],pos[:,3],color='orange',s=s,linewidths=0,alpha=0.5) #edgecolors='#ed6d05', linewidths=0
vs_pos = v_arr(np.remainder(pos[:,2]+x0, Lx),p)
f=0.05
ax1.quiver(pos[:,2],pos[:,3],vs_pos*np.cos(pos[:,4])*f,vs_pos*np.sin(pos[:,4])*f,
           color='black',headwidth=4.5)


ax2.plot(prof_avgs_NI['abp'][:,0], prof_avgs_NI['abp'][:,1], color='black', ls='--', label='NI')
# ax2.plot(prof_MF[:,0], prof_MF[:,1], color='black', label='MF')
ax2.plot(prof_avgs['abp'][:,0], prof_avgs['abp'][:,1], color='orange', label='ABP')
ax2.plot(prof_avgs['rtp'][:,0], prof_avgs['rtp'][:,1], color='purple', label='RTP')


ax3.plot([-100000,-100000],[-1,0],markersize=4,marker='o',color='black',ls='--',markerfacecolor='white',label='NI')
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
# ylim = ax3.get_ylim()
# for i in range(len(disps['abp'])):
#     ax3.plot(disps['abp'][i][:,0], disps['abp'][i][:,1],color='orange',alpha=0.1,zorder=1,lw=1)
#     ax3.plot(disps['rtp'][i][:,0], disps['rtp'][i][:,1],color='purple',alpha=0.1,zorder=1,lw=1)
# ax3.set_ylim(*ylim)

ax4.errorbar(disp[:,0]*N/Ly, disp[:,1], yerr=disp[:,3]*3, color='orange', label='ABP', marker='o',markersize=5)
ax4.errorbar(disp[:,0]*N/Ly, disp[:,2], yerr=disp[:,4]*3, color='purple', label='RTP', marker='o',markersize=5)
ax4.errorbar(disp[:1,0]*N/Ly,disp[:1,1],yerr=disp[:1,3]*3,color='orange', marker='o', markersize=5, markerfacecolor='white')
ax4.errorbar(disp[:1,0]*N/Ly,disp[:1,2],yerr=disp[:1,4]*3,color='purple', marker='o', markersize=5, markerfacecolor='white')


xlim = ax4.get_xlim()
ylim = ax4.get_ylim()
epsilons = np.linspace(0,disp[-1,0]*N/Ly,100)
plt.plot(epsilons, sum([epsilons**n * Jns[n] for n in range(len(Jns))]),
         color='black',label=r'$J_{1} N \varepsilon/L_y$')
ax4.plot([xlim[0]-10,xlim[1]+10],[0,0],color='black',lw=0.5)
ax4.set_ylim(*ylim)
ax4.set_xlim(*xlim)

for ax in axs:
    ax.tick_params(axis='both',which='both',direction='in',pad=2)
    
for ax in axs[:2]:
    ax.set_xticks([0,5,10,15],labels=['']*4)
    # ax.tick_params(axis='x', which='both', length=0, labelbottom=False)

    
ax0.set_yticks([0,10,20])
ax0.set_ylabel(r'$v(x)$',rotation='horizontal')
ax0.yaxis.set_label_coords(-0.062,0.7)
ax0.set_ylim(ymax=38)

ax1.set_ylabel(r'$y$',rotation='horizontal')
ax1.yaxis.set_label_coords(-0.03,0.55)
ax1.set_yticks([])

ax2.set_xlabel(r'$x$')
ax2.set_xticks([0,5,10,15])
ax2.set_yticks([0,0.1])
ax2.xaxis.set_label_coords(0.92,-0.03)
ax2.set_ylabel(r'$\rho_{\rm s}(x)$',rotation='horizontal')
ax2.yaxis.set_label_coords(-0.063,0.75)
ax2.legend(frameon=False,loc='upper right')

ax3.set_xticks([0,5000,10000,15000])
ax3.set_xlabel(r'$t$')
ax3.xaxis.set_label_coords(0.95,-0.03)
ax3.set_ylabel(r'$\Delta x(t)/L_x$',rotation='horizontal')
ax3.yaxis.set_label_coords(.15,0.35)
ax3.set_xlim(xmin=0,xmax=final_time)
ax3.set_ylim(ymax=1.2)
ax3.tick_params(axis='y', which='both', pad=1)
ax3.legend(frameon=False,loc='lower center',bbox_to_anchor=(0.2,0))

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

for ax in axs[:3]:
    ax.set_xlim(xmin=0,xmax=Lx)
    ax.set_ylim(ymin=0)


# add (a), (b), etc.
edgewidth = 0.5
pad = 1.25
fontsize=8.0
labels=['(a)','(b)','(c)']
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


plt.subplots_adjust(hspace=0,left=0.039,bottom=0.088,right=0.975,top=0.99)
    
plt.savefig('weak_rtp_abp_2d_current.pdf',dpi=400)
plt.close()