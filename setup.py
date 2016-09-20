import os
import sys
from setuptools import setup, find_packages, Extension

version = open('mp_ephem/__version__.py').read().split("=")[1].strip()

dependencies = ['astropy >= 1.0',
                'numpy >= 1.6.1',
                'six']

if sys.version_info[0] > 2:
    print('mp_ephem package is only compatible with Python version 2.7+, not yet with 3.x')
    sys.exit(-1)

setup(name='mp_ephem',
      version=version,
      url='http://github.com/OSSOS/liborbfit',
      author='''JJ Kavelaars (jjk@uvic.ca),
              Michele Bannister (micheleb@uvic.ca)''',
      maintainer='M Bannister and JJ Kavelaars',
      maintainer_email='jjk@uvic.ca',
      description="BK Orbfit",
      long_description='See http://www.ossos-survey.org/ for science details.',
      classifiers=['Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Development Status :: 3 - Alpha',
                   'Programming Language :: Python :: 2 :: Only',
                   'Operating System :: MacOS :: MacOS X',
                   'Environment :: X11 Applications',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   ],
      package_data={'mp_ephem': ['data/*']},
      install_requires=dependencies,
      packages=find_packages(exclude=['test', ]),
      ext_modules=[Extension('mp_ephem.orbit', [
          'src/aeiderivs.c',
          'src/covsrt.c',
          'src/dms.c',
          'src/ephem_earth.c',
          'src/gasdev.c',
          'src/gaussj.c',
          'src/lubksb.c',
          'src/ludcmp.c',
          'src/mrqcof_orbit.c',
          'src/mrqmin_orbit.c',
          'src/nrutil.c',
          'src/orbfit1.c',
          'src/orbfit2.c',
          'src/orbfitmodule.c',
          'src/ran1.c',
          'src/orbfit.h',
          'src/nrutil.h',
          'src/ephem_types.h',
          'src/ephem_read.h',
          'src/transforms.c'], extra_compile_args=['-Wno-unused-variable'])],
      )
