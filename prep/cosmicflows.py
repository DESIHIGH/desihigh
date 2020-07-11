import numpy             as np
import pandas            as pd
import pylab             as pl
import matplotlib.pyplot as plt


cflows      = pd.read_csv('../CosmicFlows3.txt', comment='#')
cflows      = cflows[['Name', 'Dist', 'RAJ', 'DeJ']]

names       = set(cflows['Name'])
names       = list(names)

ngc         = [x for x in names if type(x) == type('string')]
ngc         = [x for x in ngc   if x[:3] == 'NGC']

isin        = np.isin(cflows['Name'], ngc)

ngc         = cflows[isin]

lens        = np.array([len(x) for x in np.array(ngc['DeJ']).tolist()])

ngc['DEC']  = np.array([x[-2:]   for x in np.array(ngc['DeJ']).astype(np.str).tolist()]).astype(np.float) / 3600.
ngc['DEC'] += np.array([x[-4:-2] for x in np.array(ngc['DeJ']).astype(np.str).tolist()]).astype(np.float) / 60.
ngc['DEC']  = np.array([x[-6:-4] for x in np.array(ngc['DeJ']).astype(np.str).tolist()])

print(lens)

# pl.plot(cflows['RAJ'], cflows['DeJ'], c='k', marker='.', lw=0.0)

# pl.xlabel('RAJ')
# pl.ylabel('DeJ')

# plt.tight_layout()
# pl.show()
