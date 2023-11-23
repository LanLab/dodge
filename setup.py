from setuptools import setup,find_packages
from dodge import __version__

def readme():
    with open('README.md') as f:
        return f.read()


setup(name='dodge',
      version=__version__,
      description='Dynamic Outbreak Detection for Genomic Epidemiology',
      long_description=readme(),
      classifiers=[
          'License :: OSI Approved :: GPLv3',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Intended Audience :: Science/Research',
      ],
      keywords='genomic epidemiology, outbreak detection, foodborne pathogen, genomic surveillance',
      url='https://github.com/LanLab/dodge',
      author='Michael Payne',
      author_email='michael.payne@unsw.edu.au',
      license='GPLv3',
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': ['dodge=dodge.dodge:main','dodgedists=dodge.dodgedists:main'],
      },
      zip_safe=False)
