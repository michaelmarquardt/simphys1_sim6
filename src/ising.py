from numpy import *
from matplotlib.pyplot import *
import matplotlib
from cising import *
import argparse
import pickle 
import sys, os

# Command line arguments
parser = argparse.ArgumentParser()

parser.add_argument( '--L', type=int, nargs='+', default=[4, 64], help='--L <L1> <L2> ...\tsystem sizes L to simulate.' )

args = parser.parse_args()

# Beginning of the program
L = 10
Ts = arange(1.0, 5.1, 0.1)
num_sweeps = 1000

def compute_energy(sigma):
    L = sigma.shape[0]
    shifted_ixs = range(1,L) + [0]
    E = -(sigma*sigma[shifted_ixs,:]).sum()
    E -= (sigma*sigma[:,shifted_ixs]).sum()
    return E

def compute_magnetization(sigma):
    return sigma.sum()

def compute_act_error(x):
    x = asarray(x)
    N = len(x)
    xmean = x.mean()
    xvar = x.var()
    acfs = []
    tau_ints = []
    k = 0
    tau_int = 0.5
    while k < 6*tau_int:
        acf = ((x[:N-k]*x[k:]).mean() - xmean*xmean) / xvar
        tau_int += acf
        N_eff = N/(2*tau_int)
        err_tau = tau_int*sqrt(12./N_eff)

        acfs.append(acf)
        tau_ints.append(tau_int)
        k += 1

    err_x = sqrt(xvar/N*2.0*tau_int)
    return xmean, xvar, err_x, tau_int, err_tau, N_eff, \
        array(acfs), array(tau_ints)
        
def binder_parameter(ms):
    '''
    computes the binder parameter U=1-1/3*<mu^4>/<mu^2>^2
    '''
    return 1-1./3.*(ms**4).mean()*(ms**2).mean()**-2
    
    
##################################################
## EXACT SUMMATION
##################################################
def exact_sum(L, Ts):
    # we compute the mean energy and magnetization for all
    # temperatures at once!
    ws = zeros_like(Ts)
    Es = zeros_like(Ts)
    ms = zeros_like(Ts)
    # beta is a NumPy array with len(Ts) elements
    beta = 1./Ts

    V = float(L*L)

    sigma = ones((L, L), dtype=int)
    # the bit pattern of the integer "state" is used to generate all
    # possible states of sigma
    for state in range(2**(L*L)):
        # read out the bitpattern
        for i in range(L):
            for j in range(L):
                k = i*L + j
                if state & 2**k > 0:
                    sigma[i,j] = 1
                else:
                    sigma[i,j] = -1

        if state%10000==0: print state

        # compute energy and magnetization of this state
        E = compute_energy(sigma)
        mu = compute_magnetization(sigma)

        # this is a vector operation, as beta is a vector
        w = exp(-beta*E)
        ws += w
        Es += E/V*w
        ms += abs(mu)/V*w

    Emeans = Es/ws
    mmeans = ms/ws
    return Emeans, mmeans

# Main program
'''
print "exact summation"
Ts = arange(1.0, 5.1, 0.1)
Emeans, mmeans = exact_sum(4, Ts)
for i in range(len(Ts)):
    print "\tT = {} E = {} m = {}".format(Ts[i], Emeans[i], mmeans[i])

figure(0)
subplot(211, title='Energy vs. Temperature')
plot(Ts, Emeans, 'o-', label='exact')
legend()

subplot(212, title='Magnetization vs. Temperature')
plot(Ts, mmeans, 'o-', label='exact')
legend()
'''


##################################################
## MONTE CARLO
##################################################
def monte_carlo_ising(L, T, num_sweeps):
    V = L*L
    beta = 1.0/T

    # generate random configuration
    sigma = random.randint(0, 2, (L, L), dtype=int32)
    sigma *= 2
    sigma -= 1

    E = compute_energy(sigma)
    mu = sigma.sum()

    Es = []
    ms = []
    
    # Simulation
    for sweep in range(num_sweeps):
        E, mu = spin_flip(beta, E, mu, sigma)

        Es.append(E/float(V))
        ms.append(abs(mu)/float(V))

        print "\tT = {} {:5}/{:5}\r".format(T, sweep, num_sweeps),
        sys.stdout.flush()

    Emean, _, Eerr, tauE, _, _, _, _ = compute_act_error(array(Es))
    mmean, _, merr, tauM, _, _, _, _ = compute_act_error(array(ms))
    U = binder_parameter(array(ms))
    print "\rT = {:.1f} tau_E = {:.2f} tau_M = {:.2f} E = {:.3f}+/-{:.3f} m = {:.3f}+/-{:.3f}"\
        .format(T, tauE, tauM, Emean, Eerr, mmean, merr)

    return Emean, Eerr, mmean, merr, sigma, U

# Main program
plotname = 'L'
for L in args.L:
    print "MC (L={})".format(L)
    
    # Generating names
    plotname += '_{}'.format(L)
    datafilename = '../dat/data_L{}'.format(L)
    
    Emeans = []
    Eerrs = []
    mmeans = []
    merrs = []
    sigmas = []
    Us = []

    for T in Ts:
        Emean, Eerr, mmean, merr, sigma, U = monte_carlo_ising(L, T, num_sweeps)
        Emeans.append(Emean)
        Eerrs.append(Eerr)
        mmeans.append(mmean)
        merrs.append(merr)
        sigmas.append(sigma)
        Us.append(U)

    # Write data to file
    datafile = open(datafilename, 'w')
    pickle.dump([Emeans, Eerrs, mmeans, merrs, sigma, Us], datafile)
    datafile.close()
    
    # Plot E and mu
    fig1 = figure(0)    
    subplot(311, title='Energy vs. Temperature')
    errorbar(Ts, Emeans, yerr=Eerrs, fmt='o-', label='MC L={}'.format(L))
    legend()

    subplot(312, title='Magnetization vs. Temperature')
    errorbar(Ts, mmeans, yerr=merrs, fmt='o-', label='MC L={}'.format(L))
    legend()
    
    subplot(313, title='Binder Parameter vs. Temperature')
    plot(Ts, Us, 'o-', label='MC L={}'.format(L))
    legend()
    
    # Plot the final states
    fig2 = figure('Final states')
    numplots = len(sigmas)
    cols = int(ceil(sqrt(numplots)))
    rows = int(ceil(numplots/float(cols)))
    for i in range(numplots):
        subplot(rows, cols, i+1, title='T={}'.format(Ts[i]))
        axis('off')
        imshow(sigmas[i], interpolation='nearest')
    matplotlib.rcParams.update({'font.size': 6})
    tight_layout()
    fig2.savefig('../dat/states_L{}.pdf'.format(L))
    close(fig2)

fig1.savefig('../dat/E_mu_'+plotname+'.pdf')
close(fig1)
