
from distutils.core import setup, Extension

ORBSUBS = ['orbfit1.c','orbfit2.c', 'transforms.c','dms.c','ephem_earth.c']
NR      = ['nrutil.c', 'ludcmp.c', 'lubksb.c', 'gaussj.c', 'mrqmin_orbit.c',
           'mrqcof_orbit.c', 'covsrt.c', 'gasdev.c', 'ran1.c']

module1= Extension('bk_orbfit',
                   sources=['orbfitmodule.c']+ORBSUBS+NR,
                   libraries=['m'])

setup(name='bk_orbfit',
      version='0.1',
      description='Python binding on B&K Orbfit.',
      author='JJ Kavelaars',
      author_email='jj.kavelaars@nrc-cnrc.gc.ca',
      url='http://salish.dao.nrc.ca/~kavelaar/',
      ext_modules=[module1])

                   
                   
                   
