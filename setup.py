
from distutils.core import setup, Extension

ORBSUBS = ['orbfit1.c','transforms.c','dms.c','ephem_earth.c']
NR      = ['nrutil.c', 'ludcmp.c', 'lubksb.c', 'gaussj.c', 'mrqmin_orbit.c',
           'mrqcof_orbit.c', 'covsrt.c', 'gasdev.c', 'ran1.c']

module1= Extension('orbfit',
                   sources=['orbfitmodule.c']+ORBSUBS+NR,
                   libraries=['m'])

setup(name='orbfit',
      version='1.0',
      description='KBO orbfit return from G. Berstein',
      author='JJ Kavelaars',
      author_email='jjk@hia.nrc.ca',
      url='http://salish.dao.nrc.ca/~kavelaar/',
      ext_modules=[module1])

                   
                   
                   
