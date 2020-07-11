import pylab      as     pl
import pandas     as     pd
import numpy      as     np
import matplotlib.pyplot as plt

from   matplotlib import pyplot as plt
from   scipy      import stats


def plot_xkcd():
  plt.xkcd()

  plt.rcParams["font.family"] = "xkcd-font"
  
  fig = plt.figure(figsize=(10, 7.5))
  ax  = fig.add_subplot(1, 1, 1)

  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')

  return  fig, ax


