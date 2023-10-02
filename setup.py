import glob
import sys
from setuptools import setup, find_packages, Extension

version = "0.13.0"

dependencies = ['astropy', 'pyerfa', 'numpy', 'requests']

sources = [x for x in glob.glob('orbfit/*.c')]
headers = [x for x in glob.glob('orbfit/*.h')]

setup(name='mp_ephem',
      version=version,
      url='http://github.com/OSSOS/liborbfit',
      author='''JJ Kavelaars (jjk@uvic.ca), Michele Bannister (micheleb@uvic.ca)''',
      maintainer='M Bannister and JJ Kavelaars',
      maintainer_email='jjk@uvic.ca',
      description="BK Orbfit",
      long_description='See http://www.ossos-survey.org/ for science details.',
      long_description_content_type='text/x-rst',
      classifiers=['Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Development Status :: 3 - Alpha',
                   'Programming Language :: Python :: 3 :: Only',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   ],
      scripts = ['scripts/minor_planet_ephemeris.py'],
      package_data={'mp_ephem': ['data/*']},
      install_requires=dependencies,
      packages=find_packages(where='src'),
      package_dir = {"": "src"},
      headers = headers,
      ext_modules=[Extension('mp_ephem.orbfit', sources,
                             extra_compile_args=['-Wno-unused-variable'],
                             )],
      )
