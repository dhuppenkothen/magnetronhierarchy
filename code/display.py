from pylab import *
import numpy as np

data = np.loadtxt('090122283_+096.17100_data.txt')
posterior_sample = atleast_2d(loadtxt('posterior_sample.txt'))

ion()
for i in xrange(0, posterior_sample.shape[0]):
  hold(False)
  plot(data[:,0], data[:,1], 'bo')
  hold(True)
  plot(data[:,0], posterior_sample[i, -data.shape[0]:], 'r-')
  ylim([0, 1.2*data[:,1].max()])
  draw()
ioff()
show()

hist(posterior_sample[:,12], 200)
xlabel('Number of Bursts')
show()

pos = posterior_sample[:, 13:113]
pos = pos[pos != 0.]
hist(pos, 1000)
xlabel('Time')
title('Positions of Bursts')
show()

