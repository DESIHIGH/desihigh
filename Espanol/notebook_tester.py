import os
import papermill as pm

notebooks =['Intro',\
            'DESI',\
            'DesigningDESI',\
            'SnowWhiteDwarf',\
            'Clusters',\
            'FromMayaToDESI',\
            'SupernovaeBrain',\
            'nbody']
    
class TestSuite(object):
    def test_all(self, tmpdir):
        for nb in notebooks:
            pm.execute_notebook(
                './{}_es.ipynb'.format(nb),
                '{}/{}_es.ipynb'.format(tmpdir, nb),
                parameters=dict(alpha=0.6, ratio=0.1)
            )

            return
