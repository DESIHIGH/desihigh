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
        for lang in ['', '_es']:
            for nb in notebooks:
                pm.execute_notebook(
                    './{}{}.ipynb'.format(nb, lang),
                    '{}/{}{}.ipynb'.format(tmpdir, nb, lang),
                    parameters=dict(alpha=0.6, ratio=0.1)
                )

            os.chdir('./Espanol/')
