import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os,sys

sys.path.append(os.getcwd())
import read_data as rd
import mydefaults as md

plt.rcParams.update(md.plt_rcparams)


# L=10, v0=20, a1=1, a2=0.8
Jns = [0.0, 0.00016024329392465422, -8.134845248146115e-07, -1.4033438024323694e-07, 4.648862326420472e-10, -4.5923807968686305e-11]
radius = 24.01 # approximate radius of convergence

L=10
v0=20
a1=1
a2=0.8
dt=0.0005

cmap = plt.get_cmap('viridis')
cmap_N = plt.get_cmap('plasma')

N_folders = [x for x in os.listdir(os.getcwd()) if x.startswith('N') and os.path.isdir(x)]
N_folders = sorted(N_folders, key=lambda x: int(x.split('/')[-1].replace('N',''))) 

Js_all_N = []
for N_folder in N_folders:
    N = int(N_folder.split('/')[-1].replace('N',''))
    Js_all = np.genfromtxt(os.path.join(N_folder,'J_avg'),delimiter='\t')
    Js_all_N += [Js_all]

fig = plt.figure(figsize=(8.6/2.54,8.6*0.4/2.54))
gs = gridspec.GridSpec(1, 2, width_ratios=[0.85, 1], wspace=0, hspace=0)
axs = [fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1])]

pdf_path = "mean-field-schematic-thin.pdf"
pages = convert_from_path(pdf_path, dpi=1000)  # Increase dpi for higher resolution
image = pages[0]
axs[0].imshow(image,extent=(0,1,0,1),aspect='equal')
axs[0].set_xlim(0,1)
axs[0].set_ylim(0.07,1)
axs[0].axis('off')

Ns = [x[0,0] for x in Js_all_N]
lines = []
for i in range(len(Js_all_N)):
    lines += [plt.errorbar(Js_all_N[i][:,1], Js_all_N[i][:,2],
                         color=cmap_N(i/float(len(Js_all_N)-1)),label=str(int(Js_all_N[i][1,0])),fmt='o',ls='',markersize=5,
                         mec='black',mew=0.5)]
Nepsilons = [3,6,9,12,15,18,21,24]
Nepsilons_grid = np.linspace(0, 21.5, 100)
axs[1].plot(Nepsilons_grid, Jns[1]*Nepsilons_grid, color='gray', lw=2)
axs[1].plot(Nepsilons_grid, sum([Jns[n]*Nepsilons_grid**n for n in range(len(Jns))]), color='black', lw=2)#, label='theory')
axs[1].set_xlabel(r'$N\varepsilon$')
axs[1].set_ylabel(r'$J$',rotation='horizontal')

legend = axs[1].legend(ncol=2,frameon=False,loc='lower right',title=r'$N$=',
                    bbox_to_anchor=(1.035,-0.035),handletextpad=0,columnspacing=0.75)
legend.get_title().set_position((-20, 0))
axs[1].set_ylim(ymin=0,ymax=0.0026)

axs[1].set_xlim(xmin=0,xmax=22)
axs[1].set_yticks([0.001,0.002],labels=['0.001','0.002'])
axs[1].set_xticks([0,9,18],labels=['0','9','18'])
axs[1].tick_params(axis='y',direction='in',pad=-24)
axs[1].tick_params(axis='x',pad=2,direction='in')
axs[1].yaxis.set_label_coords(0.05,.5)
axs[1].xaxis.set_label_coords(0.94,-0.02)

plt.subplots_adjust(left=0,right=0.99,bottom=0.1,top=0.99,wspace=0)

plt.savefig('J_all_no_title.pdf',dpi=1000)

