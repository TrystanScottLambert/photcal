from distutils.core import setup

setup(
  name = 'photcal',
  packages = ['photcal'],
  version = '0.0.2',
  license='MIT',
  description = 'Package for photometrically calibrating images with survey catalogs.',
  author = 'Trystan Scott Lambert and Dejene Zwedie',
  author_email = 'trystanscottlambert@gmail.com',
  url = 'https://github.com/TrystanScottLambert',
  download_url = 'https://github.com/TrystanScottLambert/photcal/archive/refs/tags/v0.0.2.tar.gz',
  keywords = ['astronomy', 'photometry', 'calibration'],
  install_requires=[
    'numpy',
    'astropy',
    'matplotlib',
   ],

  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
  ],
)
