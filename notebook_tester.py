import os
import papermill as pm

pm.execute_notebook(
   './Intro.ipynb',
   'desihigh/pmout/Intro.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Espanol/Intro_es.ipynb',
   'desihigh//pmout/Intro_es.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './DESI.ipynb',
   'desihigh//pmout/DESI.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

#pm.execute_notebook(
#   './Espanol/DESI_es.ipynb',
#   'desihigh/pmout/DESI_es.ipynb',
#   parameters=dict(alpha=0.6, ratio=0.1)
#)

pm.execute_notebook(
   './DesigningDESI.ipynb',
   'desihigh/pmout/DesigningDESI.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './Espanol/DesigningDESI_es.ipynb',
   'desihigh/pmout/DesigningDESI_es.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

pm.execute_notebook(
   './SnowWhiteDwarf.ipynb',
   'desihigh/pmout/SnowWhiteDwarf.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

#pm.execute_notebook(
#   './Espanol/SnowWhiteDwarf_es.ipynb',
#   'desihigh/pmout/SnowWhiteDwarf_es.ipynb',
#   parameters=dict(alpha=0.6, ratio=0.1)
#)

pm.execute_notebook(
   './Clusters.ipynb',
   'desihigh/pmout/Clusters.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

#pm.execute_notebook(
#   './Espanol/Clusters_es.ipynb',
#   'desihigh/pmout/Clusters_es.ipynb',
#   parameters=dict(alpha=0.6, ratio=0.1)
#)

pm.execute_notebook(
   './FromMayaToDESI.ipynb',
   'desihigh/pmout/FromMayaToDESI.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)

#pm.execute_notebook(
#   './Espanol/FromMayaToDESI_es.ipynb',
#   'desihigh/pmout/FromMayaToDESI_es.ipynb',
#   parameters=dict(alpha=0.6, ratio=0.1)
#)

pm.execute_notebook(
   './SupernovaeBrain.ipynb.ipynb',
   'desihigh/pmout/SupernovaeBrain.ipynb',
   parameters=dict(alpha=0.6, ratio=0.1)
)
