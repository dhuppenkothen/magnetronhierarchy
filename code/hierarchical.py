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
  num_bursts = 20
  num_params = 1   # Number of within-burst hyperparameters
  samples = []

  for i in range(0, num_bursts):
    samples.append(rng.randn() + rng.randn(10000, num_params))

  return samples


if __name__ == '__main__':
  samples = get_samples()
  plt.hist(samples[0][:,0], 100)
  plt.show()

