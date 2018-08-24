from setuptools import setup
from setuptools.extension import Extension



with open('gpbo/VERSION') as version_file:
    version = version_file.read().strip()

def readme():
    with open('README.md') as f:
        return f.read()

compile_flags = ['-O3']
from numpy import get_include
extensions = [
     Extension(name ="gpbo/core/ESutils",
              sources = ["gpbo/core/ESutils.c"],
              include_dirs = ['.','core',get_include()],
              extra_compile_args=compile_flags
    ),
    Extension(name ="gpbo/core/GPdc",
              sources = ["gpbo/core/GPdc.c"],
              include_dirs = ['.','core',get_include()],
              extra_compile_args=compile_flags
    ),
    Extension(name ="gpbo/core/PES",
              sources = ["gpbo/core/PES.c"],
              include_dirs = ['.','core',get_include()],
              extra_compile_args=compile_flags
    ),
    Extension(name ="gpbo/core/eprop",
              sources = ["gpbo/core/eprop.c"],
              include_dirs = ['.','core',get_include()],
              extra_compile_args=compile_flags
    ),
    Extension(name ="gpbo/core/slice",
              sources = ["gpbo/core/slice.c"],
              include_dirs = ['.','core',get_include()],
              extra_compile_args=compile_flags
    ),
    Extension(name ="gpbo/core/acquisitions",
              sources = ["gpbo/core/acquisitions.c"],
              include_dirs = ['.','core',get_include()],
              extra_compile_args=compile_flags
    ),
    Extension(name ="gpbo/core/optutils",
              sources = ["gpbo/core/optutils.c"],
              include_dirs = ['.','core',get_include()],
              extra_compile_args=compile_flags
    )
]

setup(name='gpbo',
      version=version,
      description='a package',
      long_description=readme(),
      url='https://github.com/markm541374/gpbo',
      author='markm541374',
      license='MIT',
      packages=['gpbo','gpbo.core','gpbo.examples','gpbo.exps'],
      package_dir={'gpbo':'gpbo'},
      package_data={'gpbo':['cproj/*','VERSION','README.rst']},
      install_requires=['numpy','scipy','tqdm','direct','matplotlib','pandas','emcee','cvxopt','cma','sklearn'],
      ext_modules= extensions,
      zip_safe=False)
