import os
import papermill as pm

pm.execute_notebook(
   './Intro.ipynb',
   './pmout/Intro.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './DESI.ipynb',
   './pmout/DESI.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './DesigningDESI.ipynb',
   './pmout/DesigningDESI.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './SnowWhiteDwarf.ipynb',
   './pmout/SnowWhiteDwarf.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Clusters.ipynb',
   './pmout/Clusters.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)
