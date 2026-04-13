import numpy as np
import matplotlib.pyplot as plt
import time
import os,sys
from tqdm import tqdm
from matplotlib import gridspec

sys.path.append(os.getcwd())
import read_data as rd
import mydefaults as md

plt.rcParams.update(md.plt_rcparams)


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

alphas = ['1.00','2.00','4.00','8.00']

currents_all = []
currents_all_th = []
for alpha in tqdm(alphas):
    alpha_str = alpha

    # SIMULATION CURRENT

    epsilons = [x.split('epsilon')[-1] for x in os.listdir('alpha'+str(alpha)) if 'epsilon' in x]
    epsilons = sorted(epsilons, key=float)
    
    currents = []
    for epsilon in epsilons:
        
        run_name=os.path.join(os.getcwd(),'alpha'+alpha_str+'/epsilon'+str(epsilon)+'/res')
        p = rd.get_params(run_name+'-param')

        disp = np.loadtxt(dispf, delimiter='\t')
        start_disp=np.zeros((int(p.N),3))
        start_disp[:,1] = np.arange(0,p.N,1).astype(int)
        disp = np.concatenate((start_disp,disp),axis=0)
        disp[:,2] = disp[:,2]/p.L
        times = np.unique(disp[:,0])
        ids = np.unique(disp[:,1])
        displacement_avg=[]
        for t in times[-1:]:
            now_data = disp[np.where(disp[:,0]==t)[0],:]
            displacement_avg += [np.mean(now_data[:,2])]
        displacement_avg = np.array(displacement_avg)
        J_sim = displacement_avg[-1]/times[-1]
        currents += [[epsilon,J_sim]]

    currents = sorted(currents, key = lambda x: x[0])
    currents = np.array(currents)
    currents_all += [[alpha,currents]]

    np.savetxt('current_'+alpha_str,currents)


    # THEORY CURRENT
        
    epsilons_fine = np.arange(0,.88,0.01)
    if int(alpha)==4:
        epsilons_fine = np.arange(0,.851,.01)
    currents_th = []
    for epsilon in epsilons_fine:

        # system domain
        xs = np.arange(0,L+L/50000.,L/50000.)
    
        # activity landscape
        dx = L/(5000.*P)
        xs = np.arange(0, L+dx, dx)
        vs = v(xs, p.v0, p.a1, p.a2, p.d, p.P,p.L)
        vps = vp(xs, p.v0, p.a1, p.a2, p.d, p.P,p.L)
        Us = U(xs,p.epsilon,p.a1,p.a2,p.b1,p.b2,p.P,p.L)
        Ups = Up(xs,p.epsilon,p.a1,p.a2,p.b1,p.b2,p.P,p.L)
        Upps = Upp(xs,p.epsilon,p.a1,p.a2,p.b1,p.b2,p.P,p.L)
        J_exact = functions.J_exact(xs, Us, Ups, Upps, vs, vps, 1, p.alpha)
        currents_th += [[p.epsilon, J_exact]]
        
    currents_th = sorted(currents_th, key=lambda x: x[0])
    currents_th = np.array(currents_th)
    currents_all_th += [[alpha,currents_th]]
    
    np.savetxt('current_th_'+alpha_str,currents_th)