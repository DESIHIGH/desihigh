import numpy as np
import healpy as hp

from   desitarget.geomask import brick_names_touch_hp
from   desiutil import depend, brick

nside = 32

# items = brick_names_touch_hp(nside=nside, numproc=16)

# list of names of bricks that touch HEALPixel 0.
# print(items[0])

# nside = 32
items = ['0436p000', '0438p000', '0441p000', '0443p000', '0446p000', '0448p000', '0451p000', '0453p000', '0456p000', '0458p000', '0461p000', '0463p000', '0436p002', '0438p002', '0441p002', '0443p002', '0446p002', '0448p002', '0451p002', '0453p002', '0456p002', '0458p002', '0461p002', '0463p002', '0436p005', '0438p005', '0441p005', '0443p005', '0446p005', '0448p005', '0451p005', '0453p005', '0456p005', '0458p005', '0461p005', '0463p005', '0436p007', '0438p007', '0441p007', '0443p007', '0446p007', '0448p007', '0451p007', '0453p007', '0456p007', '0458p007', '0461p007', '0463p007', '0436p010', '0438p010', '0441p010', '0443p010', '0446p010', '0448p010', '0451p010', '0453p010', '0456p010', '0458p010', '0461p010', '0463p010', '0436p012', '0438p012', '0441p012', '0443p012', '0446p012', '0448p012', '0451p012', '0453p012', '0456p012', '0458p012', '0461p012', '0463p012', '0436p015', '0438p015', '0441p015', '0443p015', '0446p015', '0448p015', '0451p015', '0453p015', '0456p015', '0458p015', '0461p015', '0463p015', '0436p017', '0438p017', '0441p017', '0443p017', '0446p017', '0448p017', '0451p017', '0453p017', '0456p017', '0458p017', '0461p017', '0463p017', '0436p020', '0438p020', '0441p020', '0443p020', '0446p020', '0448p020', '0451p020', '0453p020', '0456p020', '0458p020', '0461p020', '0463p020', '0436p022', '0438p022', '0441p022', '0443p022', '0446p022', '0448p022', '0451p022', '0453p022', '0456p022', '0458p022', '0461p022', '0463p022', '0436p025', '0438p025', '0441p025', '0443p025', '0446p025', '0448p025', '0451p025', '0453p025', '0456p025', '0458p025', '0461p025', '0463p025']

# nside = 64.
# items = ['0443p000', '0446p000', '0448p000', '0451p000', '0453p000', '0456p000', '0443p002', '0446p002', '0448p002', '0451p002', '0453p002', '0456p002', '0443p005', '0446p005', '0448p005', '0451p005', '0453p005', '0456p005', '0443p007', '0446p007', '0448p007', '0451p007', '0453p007', '0456p007', '0443p010', '0446p010', '0448p010', '0451p010', '0453p010', '0456p010', '0443p012', '0446p012', '0448p012', '0451p012', '0453p012', '0456p012']

area  = len(items) * 0.25 * 0.25

print(area, hp.pixelfunc.nside2pixarea(nside, degrees=True))



bricktable        = brick.Bricks(bricksize=0.25).to_table()

bricktable['A']   = ['{:04d}'.format(x) for x in np.array(bricktable['RA1']  * 10.).astype(np.int)]
bricktable['B']   = ['{:04d}'.format(x) for x in np.array(bricktable['DEC1'] * 10.).astype(np.int)]

bricktable['C']   = ['{:04d}'.format(x) for x in np.array(bricktable['RA2']  * 10.).astype(np.int)]
bricktable['D']   = ['{:04d}'.format(x) for x in np.array(bricktable['DEC2'] * 10.).astype(np.int)]

bricktable['B']   = [x.replace('-', 'm') if '-' in x else 'p' + x[1:] for x in bricktable['B']] 
bricktable['D']   = [x.replace('-', 'm') if '-' in x else 'p' + x[1:] for x in bricktable['B']]

bricktable['EXT'] = ['                                    '] * len(bricktable)

for i, _ in enumerate(bricktable['EXT']):
  bricktable['EXT'][i] = bricktable['A'][i] + bricktable['B'][i] + '-' + bricktable['C'][i] + bricktable['D'][i]

print(bricktable)
