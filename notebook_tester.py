import os
import papermill as pm

pm.execute_notebook(
   './Intro.ipynb',
   './pmout/pmout.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './DESI.ipynb',
   './pmout/pmout.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './DesigningDESI.ipynb',
   './pmout/pmout.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)
