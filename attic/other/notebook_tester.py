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
                './{}.ipynb'.format(nb),
                '{}/{}.ipynb'.format(tmpdir, nb),
                parameters=dict(alpha=0.6, ratio=0.1)
            )
