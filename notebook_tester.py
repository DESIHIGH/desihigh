import os
import papermill as pm

pm.execute_notebook(
   './Intro.ipynb',
   './pmout/Intro.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Espanol/Intro_es.ipynb',
   './pmout/Intro_es.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './DESI.ipynb',
   './pmout/DESI.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Espanol/DESI_es.ipynb',
   './pmout/DESI_es.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './DesigningDESI.ipynb',
   './pmout/DesigningDESI.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Espanol/DesigningDESI_es.ipynb',
   './pmout/DesigningDESI_es.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './SnowWhiteDwarf.ipynb',
   './pmout/SnowWhiteDwarf.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Espanol/SnowWhiteDwarf_es.ipynb',
   './pmout/SnowWhiteDwarf_es.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Clusters.ipynb',
   './pmout/Clusters.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Espanol/Clusters_es.ipynb',
   './pmout/Clusters_es.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './FromMayaToDESI.ipynb',
   './pmout/FromMayaToDESI.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Espanol/FromMayaToDESI_es.ipynb',
   './pmout/FromMayaToDESI_es.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)
