import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

rng.seed(0)

def get_samples():
  """
  Load samples of within-burst hyperparameters
  Actually, for now, just generate some fake ones
  Returns a list of 2D numpy arrays. Each element
  of the list is the posterior samples for one burst.
  """
  num_bursts = 3
  num_params = 1   # Number of within-burst hyperparameters
  samples = []

  for i in range(0, num_bursts):
    samples.append(rng.randn() + rng.randn(10000, num_params))

  return samples

def log_likelihood(mu, samples):
  """
  Log of the marginal likelihood for the hyperparameters
  The marginal likelihood is:
  p(data | alpha) = \int p(theta | alpha) p(data | theta) dtheta
                  = Z \int p(theta | alpha)/pi(theta) p(data | theta) pi(theta)/Z dtheta
  """
  logL = 0.
  for i in range(len(samples)):
    pi = 1./200
    integrand = np.exp(-0.5*(samples[i][0] - mu)**2)/pi
    logL += np.log(np.mean(integrand))

  return logL




if __name__ == '__main__':
  # Generate or load samples
  samples = get_samples()

  for i in range(0, len(samples)):
    plt.hist(samples[i][:,0], 100, alpha=0.2)
  plt.show()

  print(log_likelihood(0., samples))

