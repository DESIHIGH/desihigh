import requests
import matplotlib.pyplot       as     plt

from PIL import Image


def cutout(ax, ra, dec):  
  url       = 'http://legacysurvey.org/viewer/jpeg-cutout?ra={:.4f}&dec={:.4f}&layer=decals-dr7&pixscale=0.27&bands=grz'.format(ra, dec)

  image     = Image.open(url)
  
  ax.imshow(image)
  
