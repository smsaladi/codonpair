
import os
import os.path
from setuptools import setup

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

version = '0.1.1'
# For travis/pypi-test
if os.environ.get('TRAVIS') == 'true' and os.environ.get('TRAVIS_TAG') == '':
      N, M = os.environ['TRAVIS_JOB_NUMBER'].split('.')
      version = "{v}-a{N}.dev{M}".format(v=version, N=N, M=M)

setup(name='codonpair',
      version=version,
      description='A Python implementation to calculate codon pair score',
      url='http://github.com/smsaladi/codonpair',
      author='Shyam Saladi',
      author_email='saladi@caltech.edu',
      license='MIT',
      packages=['codonpair'],
      install_requires=['numpy', 'pandas', 'biopython'],
      package_data={'codonpair': ['data/*']},
      entry_points={'console_scripts': ['cps=codonpair.codonpair:main']},
      long_description=long_description,
      long_description_content_type='text/markdown',
      tests_require = [
            'pytest',
      ],
      test_suite="pytest",
      zip_safe=True,
      include_package_data=True
)
