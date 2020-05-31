import numpy
import os, sys 
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate=True
from setuptools import setup, Extension


setup(
    name='pygl',
    version='1.0.0',
    url='https://github.com/rajeshrinet/pygl',
    author='The PyGL team',
    license='MIT',
    description='python library for numerical simulation of field theories',
    long_description='PyGL is a library for numerical simulation of field theories\
    in non-equilibrium statistical mechanics',
    platforms='works on all platforms (such as LINUX, macOS, and Microsoft Windows)',
    ext_modules=cythonize([ Extension("pygl/*", ["pygl/*.pyx"],
        include_dirs=[numpy.get_include()],
        )],
        compiler_directives={"language_level": sys.version_info[0]},
        ),
    libraries=[],
    packages=['pygl'],
    package_data={'pygl': ['*.pxd']},
)
