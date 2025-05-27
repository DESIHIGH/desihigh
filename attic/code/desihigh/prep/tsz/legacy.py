import requests
import matplotlib.pyplot       as      plt

from   PIL                     import  Image
from   io                      import  BytesIO

def cutout(ax, ra, dec):  
  url       = 'https://www.legacysurvey.org/viewer/cutout.jpg?ra={:.4f}&dec={:.4f}&layer=dr8&pixscale=1.00'.format(ra, dec)

  try:
    # Timeout in seconds.
    response  = requests.get(url, timeout=30.)
    img       = Image.open(BytesIO(response.content))
  
    ax.imshow(img)

  except:
    print('Failed to retrieve {}'.format(url))

  return  url
  

if __name__ == '__main__':
  fig, ax = plt.subplots(1,1, figsize=(5,5,))
  
  url     = cutout(ax, 223.6459, 23.6432)

  print(url)
