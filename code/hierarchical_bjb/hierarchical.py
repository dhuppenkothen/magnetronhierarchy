import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
import copy

rng.seed(0)

def get_samples():
  """
  Load samples of within-burst hyperparameters
  Actually, for now, just generate some fake ones
  Returns a list of 2D numpy arrays. Each element
  of the list is the posterior samples for one burst.
  """
  num_bursts = 20
  num_params = 2   # Number of within-burst hyperparameters
  samples = []

  for i in range(0, num_bursts):
    samples.append(rng.randn() + 0.1*rng.randn(10000, num_params))

  return samples

# A 'global' dataset (really posterior samples from
# the non-hierarchical model)
samples = get_samples()

# How many hyperparameters are there?
num_params = 4

# Some idea of how big the Metropolis proposals should be
jump_sizes = np.hstack([20.*np.ones(num_params//2),\
                         np.log(1E6)*np.ones(num_params//2)])

def from_prior():
  """
  A function to generate parameter values from the prior.
  Returns a numpy array of parameter values.
  """
  return np.hstack([-10. + 20.*rng.rand(num_params//2),\
                      np.log(1E-3) + np.log(1E6)*rng.rand(num_params//2)])

def log_prior(params):
  """
  Evaluate the (log of the) prior distribution
  """
  mu = params[0:(num_params//2)]
  log_sig = params[(num_params//2):]

  # Minus infinity, if out of bounds
  if np.logical_or(np.any(mu < -10.), np.any(mu > 10.)):
    return -np.Inf
  if np.logical_or(np.any(log_sig < np.log(1E-3)), np.any(log_sig > np.log(1E3))):
    return -np.Inf

  return 0.

def log_likelihood(params):
  """
  Log of the marginal likelihood for the hyperparameters
  The marginal likelihood is:
  p(data | alpha) = \int p(theta | alpha) p(data | theta) dtheta
                  = Z \int p(theta | alpha)/pi(theta) p(data | theta) pi(theta)/Z dtheta
  """
  mu = params[0:(num_params//2)]
  log_sig = params[(num_params//2):]
  sig = np.exp(log_sig)

  logL = 0.
  for i in range(len(samples)):
    interim_prior = 1./200
    for j in range(len(mu)):
      integrand = np.exp(-0.5*(samples[i][j] - mu[j])**2/sig[j]**2)/sig[j]\
                         /interim_prior
      logL += np.log(np.mean(integrand))

  if np.isinf(logL) or np.isnan(logL):
    logL = -1E300

  return logL


def proposal(params):
  """
  Generate new values for the parameters, for the Metropolis algorithm.
  """
  # Copy the parameters
  new = copy.deepcopy(params)

  # Which one should we change?
  which = rng.randint(num_params)
  new[which] += jump_sizes[which]*10.**(1.5 - 6.*rng.rand())*rng.randn()
  return new

#if __name__ == '__main__':
#  # Generate or load samples
#  samples = get_samples()

#  # Display the samples
#  for i in range(0, len(samples)):
#    plt.hist(samples[i][:,0], 100, alpha=0.2)
#  plt.show()

#  # Calculate log likelihood as a function of mu
#  mu = np.linspace(-10., 10., 201)
#  logL = np.zeros(mu.shape)
#  for i in range(0, len(logL)):
#    logL[i] = log_likelihood(mu[i], samples)
#  plt.plot(mu, np.exp(logL - logL.max()))
#  plt.show()

