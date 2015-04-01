import h5py

from combine_gauss import gaussian_combination

with h5py.File('sgr1550_sample_hypers.dat', 'r') as f:
    means = f['means'][:] # MxD
    covars = f['covars'][:] # MxDxD

mu, Sigma = gaussian_combination(means, covars)

# If had raw samples in an MxSxD array xx, could do:
# gauss_from_ensemble_samples(xx)
