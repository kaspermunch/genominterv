
import setuptools

setuptools.setup(name='genominterv',
      version='1.0.1',
      author='Kasper Munch',
      description='Utilities for working with intervals in separate chromosomes.',
      # long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/kaspermunch/genominterv',
      packages=setuptools.find_packages(),
      python_requires='>=3.6',
      install_requires=[
      'pandas>=1.0',
      'numpy>=1.1',
      'statsmodels>=0.8',
      ])
