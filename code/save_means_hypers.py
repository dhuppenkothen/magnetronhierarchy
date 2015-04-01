

import glob
import h5py
import numpy as np


def save_means_covars(datadir="./", outfile="sgr1550_sample_hypers.dat"):
    ### find all posterior files
    pfiles = glob.glob(datadir+"*posterior.txt")
    mean_all, cov_all = [], []
    for p in pfiles:
        data = np.loadtxt(p)
        d_hyper = data[:,3:3+9]
        d_hyper[:,-3:] = np.log(d_hyper[:,-3:])
        m =  np.mean(d_hyper, axis=0)
        c =  np.cov(d_hyper, rowvar=False)
        mean_all.append(m)
        cov_all.append(c)

    mean_all = np.array(mean_all)
    cov_all = np.array(cov_all)
    f = h5py.File(datadir+outfile, "w")
    f.create_dataset("means", data=mean_all)
    f.create_dataset("covars", data=cov_all)
    f.close()
    return



