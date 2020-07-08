import pylab      as     pl
import pandas     as     pd
import numpy      as     np

from   matplotlib import pyplot as plt
from   scipy      import stats


def plot_xkcd():
  plt.xkcd()

  fig = plt.figure()
  ax  = fig.add_subplot(1, 1, 1)

  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')

  return  fig, ax


