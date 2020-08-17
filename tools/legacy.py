import requests
import matplotlib.pyplot       as     plt

from   PIL                     import  Image
from   io                      import  BytesIO

def cutout(ax, ra, dec):  
  url       = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra={:.4f}&dec={:.4f}&layer=dr8&pixscale=1.00'.format(ra, dec)

  response  = requests.get(url)
  img       = Image.open(BytesIO(response.content))
  
  ax.imshow(img)
  
