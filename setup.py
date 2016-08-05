from setuptools import setup
from setuptools.extension import Extension



with open('gpbo/VERSION') as version_file:
    version = version_file.read().strip()

def readme():
    with open('gpbo/README.rst') as f:
        return f.read()

setup(name='gpbo',
      version=version,
      description='a package',
      long_description=readme(),
      url='https://github.com/markm541374/gpbo',
      author='markm541374',
      license='MIT',
      packages=['gpbo','gpbo.core','gpbo.examples','gpbo.test','gpbo.config'],
      package_dir={'gpbo':'gpbo'},
      package_data={'gpbo':['cproj/*','VERSION','README.rst']},
      zip_safe=False)
