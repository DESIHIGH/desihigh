import requests
import matplotlib.pyplot       as     plt


def cutout(ax, ra, dec):  
  url       = 'http://legacysurvey.org/viewer/jpeg-cutout?ra={:.4f}&dec={:.4f}&layer=decals-dr7&pixscale=0.27&bands=grz'.format(ra, dec)

  image     = plt.imread("http://matplotlib.sourceforge.net/_static/logo2.png")
  
  ax.imshow(image)
  
