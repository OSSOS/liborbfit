
from distutils.core import setup

setup(name='bk_orbfit',
      version='0.1',
      description='Python binding on B&K Orbfit.',
      data_files=[('','data/binEphem.405'),
                  ('','data/observatories.dat')],
      author='JJ Kavelaars',
      package_dir=['bk_orbfit'],
      author_email='jj.kavelaars@nrc-cnrc.gc.ca',
      url='http://salish.dao.nrc.ca/~kavelaar/',
      ext_modules=[module1])

                   
                   
