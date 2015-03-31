from pylab import *
import numpy as np

data = np.loadtxt('090122283_+096.17100_eventfile.dat')
posterior_sample = atleast_2d(loadtxt('posterior_sample.txt'))

ion()
hold(False)
for i in xrange(0, posterior_sample.shape[0]):
  counts, bin_edges = np.histogram(data[:,0], 100)
  dt = bin_edges[1]-bin_edges[0]
  countrate = counts/dt
  plot(bin_edges[:-1]+dt/2., countrate, color="black")
  hold(True)
  plot(data[:,0], posterior_sample[i, -data.shape[0]:], 'r-')
  ylim([0, 1.1*countrate.max()])
  draw()
ioff()
show()

hist(posterior_sample[:,12], 20)
xlabel('Number of Bursts')
show()

pos = posterior_sample[:, 13:113]
pos = pos[pos != 0.]
hist(pos, 1000)
xlabel('Time')
title('Positions of Bursts')
show()

