import glob
import numpy as np
import combine_gauss

all_data = []
# Use all bursts:
filenames = glob.glob('sgr1550_hypers/*.txt')
# Or just the top few longest ones:
# filenames = ['090122120+016.086_hypers.txt', '090122187+292.952_hypers.txt',
#         '090122173+241.347_hypers.txt', '090122283+230.654_hypers.txt',
#         '090122187+262.576_hypers.txt', '090122194+052.410_hypers.txt',
#         '090122173+037.024_hypers.txt', '090122187+156.963_hypers.txt',
#         '090122187+229.756_hypers.txt', '090122187-000.051_hypers.txt',
#         '090122218+157.280_hypers.txt', '090122283+131.840_hypers.txt',
#         '090122194+058.836_hypers.txt', '090122187+248.037_hypers.txt',
#         '090122104+258.828_hypers.txt', '090122218+183.986_hypers.txt',
#         '090122120+106.245_hypers.txt', '090122218+121.883_hypers.txt',
#         '090122120+203.953_hypers.txt', '090122194-006.278_hypers.txt']
# filenames = map(lambda x: 'sgr1550_hypers/' + x, filenames)
for fn in filenames:
    all_data.append(np.loadtxt(fn))
# Each element of all_data is Sx9, where S~2000 samples
# all_data is now M bursts long

mu, Sigma = combine_gauss.gauss_from_ensemble_samples(all_data)

