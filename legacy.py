import requests

from   IPython.display         import Image
from   io                      import BytesIO


def cutout(ax, ra, dec):  
  url       = 'http://legacysurvey.org/viewer/jpeg-cutout?ra={:.4f}&dec={:.4f}&layer=decals-dr7&pixscale=0.27&bands=grz'.format(ra, dec)
  response  = requests.get(url)
  img       = Image.open(BytesIO(response.content))
  
  ax.imshow(img)
  
