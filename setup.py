from setuptools import setup

setup(name='codonpair',
      version='1.0a',
      description='A Python implementation to calculate codon pair score',
      url='http://github.com/smsaladi/codonpair',
      author='Shyam Saladi',
      author_email='saladi@caltech.edu',
      license='MIT',
      packages=['cps'],
      install_requires=['numpy', 'pandas', 'biopython'],
      package_data={'cps': ['data/*.tbd']},
      entry_points={'console_scripts': ['cps=cps.cps:main']},
      zip_safe=True,
      include_package_data=True)
