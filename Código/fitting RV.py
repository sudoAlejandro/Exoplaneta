# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
import phoebe

# Subscribe to messages printed to the screen
logger = phoebe.get_basic_logger(clevel='INFO')

# Step 1: Create synthetic radial velocity curves and add some noise. We'll use
# those as simulated data to fit

# Define a system
mybundle = phoebe.Bundle('binary')
mybundle['ecc'] = 0.34
mybundle['per0'] = 53.0
mybundle['vgamma'] = -10.0
mybundle['sma'] = 10.0
mybundle['incl'] = 67.0
mybundle['q'] = 0.66
mybundle['t0@orbit'] = 1.0

# Generate a time base to compute the RVs on, and set the uncertainties on the
# radial velocity curve
time = np.sort(np.random.uniform(low=0, high=10*mybundle['period'], size=100))
sigma = np.ones(len(time))

# --- Primary component
pos1, velo1, btime1, ptime1 = mybundle.get_orbit('primary', time=time)
noise1 = np.random.normal(scale=sigma)

# --- Secondary component
pos2, velo2, btime2, ptime2 = mybundle.get_orbit('secondary', time=time)
noise2 = np.random.normal(scale=sigma)

# Save the RVs to a file
np.savetxt('system.rvs', np.column_stack([time, velo1[2]+noise1, sigma,
                                                    velo2[2]+noise2, sigma]))

# Step 2: Create a new system, and load the previously generated radial velocity
# curves
mybundle = phoebe.Bundle('binary')

# Add the RV curve as data to the primary and secondary
mybundle.rv_fromfile('system.rvs', columns=['time', 'rv', 'sigma'],
                     objref='primary', method='dynamical')
mybundle.rv_fromfile('system.rvs', columns=['time', '', '', 'rv', 'sigma'],
                     objref='secondary', method='dynamical')

# We don't want to fit the semi-major axis because it correlates heavily with
# inclination angle. A better parameter is the projected semi-major axis. This
# is not a standard parameter so we'll have to add it. Additionally, we'll
# remove the semi-major axis as a free parameter but let be derived from asini
# instead
mybundle.add_parameter('asini@orbit', replaces='sma')

# Define priors
mybundle.set_prior('ecc', distribution='uniform', lower=0, upper=1)
mybundle.set_prior('per0', distribution='uniform', lower=0, upper=360)
mybundle.set_prior('vgamma', distribution='uniform', lower=-30, upper=10)
mybundle.set_prior('incl', distribution='uniform', lower=0, upper=90)
mybundle.set_prior('q', distribution='uniform', lower=0.5, upper=1)
mybundle.set_prior('sma', distribution='uniform', lower=6, upper=16)
mybundle.set_prior('t0@orbit', distribution='uniform', lower=0, upper=2)
mybundle.set_prior('asini', distribution='uniform', lower=0, upper=15)

# Mark the parameters that we want to include in the fit
mybundle.set_adjust('ecc')
mybundle.set_adjust('per0')
mybundle.set_adjust('vgamma')
mybundle.set_adjust('asini')
mybundle.set_adjust('incl')
mybundle.set_adjust('q')
mybundle.set_adjust('t0@orbit')

# Add fitting options: we want to use the emcee package, and we'll use 500
# iterations on 50 walkers.
mybundle.add_fitting(context='fitting:emcee', label='mcmc', iters=500, walkers=50)

# To speed up the computations, we'll use the MPI framework with 6 threads
mpi = phoebe.ParameterSet('mpi', np=6)

# Run the fitting algorithm
mybundle.run_fitting(fittinglabel='mcmc', mpi=mpi)

# Plot the history of the probabilities of chain
plt.figure()
mybundle['mcmc@feedback'].plot_logp()

# Only select the last 200 iterations, cut off in lnprob and set the resulting
# chain as the posteriors on the parameters
mybundle['mcmc@feedback'].modify_chain(lnproblim=-40, burnin=300)
mybundle.accept_feedback('mcmc')

# Plot the posteriors and the correlations
mybundle['mcmc@feedback'].plot_summary()

# Print out the report of the feedback
print(mybundle['mcmc@feedback'])

# Plot the best model    
plt.figure()
mybundle.plot_obs('rv01', fmt='ko', phased=True)
mybundle.plot_obs('rv02', fmt='ks', phased=True)
mybundle.plot_syn('rv01', 'r-', lw=2, phased=True)
mybundle.plot_syn('rv02', 'b-', lw=2, phased=True)

plt.show()