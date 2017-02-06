import glob
import sys
from setuptools import setup, find_packages, Extension
from mp_ephem import __version__

version = __version__.version

dependencies = ['astropy >= 1.0',
                'numpy >= 1.6.1',
                'six']

if sys.version_info[0] > 2:
    print('mp_ephem package is only compatible with Python version 2.7+, not yet with 3.x')
    sys.exit(-1)

sources = [x for x in glob.glob('src/*.c')]

setup(name='mp_ephem',
      version=version,
      url='http://github.com/OSSOS/liborbfit',
      author='''JJ Kavelaars (jjk@uvic.ca), Michele Bannister (micheleb@uvic.ca)''',
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
          ext_modules=[Extension('mp_ephem.orbit', sources,
                             extra_compile_args=['-Wno-unused-variable'],
                             )],
      )
