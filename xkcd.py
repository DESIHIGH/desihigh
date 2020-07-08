import pylab      as     pl
import pandas     as     pd
import numpy      as     np

from   matplotlib import pyplot as plt


def xkcd_plot():
  plt.xkcd()

  fig = plt.figure()
  ax  = fig.add_subplot(1, 1, 1)
  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')

  # plt.xticks([])
  # plt.yticks([])

  return fig, ax


if __name__ == '__main__':
  fig, ax = xkcd_plot()

  dat = pd.read_csv('hubble.dat', sep='\s+')

  ax.plot(dat['r'], dat['v'])
  
  ax.set_xlabel('Distance from us [Mpc]')
  ax.set_ylabel('Recession velocity [km/s]')

  plt.tight_layout()
  
  plt.show()
