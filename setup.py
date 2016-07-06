from distutils.core import setup, Extension

setup(packages=['mp_ephem'],
      package_dir={'mp_ephem': 'mp_ephem'},
      package_data={'mp_ephem': ['data/binEphem.405_32',
                                 'data/binEphem.405_64',
                                 'data/observatories.dat']},
      version='1.0',
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
          'src/transforms.c'], extra_compile_args=['-Wno-unused-variable'])],
      requires=['astropy', 'numpy', 'six']
      )
