import numpy as np

# TODO test this a little bit!

def gaussian_combination(mus, Sigmas):
    """
    Report the Gaussian that proportional to the product of M Gaussians

    All inputs and outputs are numpy arrays with the given sizes.

    Inputs:
        mus M,D,
     Sigmas M,D,D,

    Outputs:
         mu D,
      Sigma D,D,

    Sigmas can be Inf or a diagonal of Inf's to have an improper broad Gaussian
    with no effect on the product.

    Sigmas can be negative to divide by a Gaussian instead of multiplying.
    Negative Sigmas can lead to a negative final covariance however.
    """
    # Yes, I just naively take matrix inverses. These matrices better be well
    # conditioned, so numerics shouldn't be a problem. And currently I'm dealing
    # with small matrices, so computational cost is irrelevant.
    M, D = mus.shape
    prec = np.zeros_like(Sigmas[0])
    mu_tmp = np.zeros_like(mus[0])
    for mm in range(M):
        prec_m = np.linalg.inv(Sigmas[mm])
        prec += prec_m
        mu_tmp += np.dot(prec_m, mus[mm])
    Sigma = np.linalg.inv(prec)
    mu = np.dot(Sigma, mu_tmp)

    return (mu, Sigma)

def gauss_from_ensemble_samples(xx):
    """
    Estimate a Gaussian that's a product of distributions, from samples

    Inputs:
        xx M,S,D Ensemble of size M. Each has S, D-dimensional samples
           or a list of M arrays of size S_m,D

    We estimate M Gaussians, then find the Gaussian proportional to the product
    of the Gaussians.
    """
    M = len(xx)
    S0, D = xx[0].shape
    Sigmas = np.zeros((M,D,D))
    mus = np.zeros((M,D))
    for mm in range(M):
        # TODO consider smoothing diagonal
        mus[mm] = xx[mm].mean(0)
        Sigmas[mm] = np.cov(xx[mm], rowvar=False)

    return gaussian_combination(mus, Sigmas)

# TODO if the agents in the ensemble used Gaussian priors we want to divide out,
# we'll need to get those from somewhere, and add them to the list, with
# negative Sigmas. If there's a global non-flat prior we want to include, we'll
# have to include that in the product too.

