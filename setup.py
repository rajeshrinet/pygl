import numpy
import os, sys 
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate=True
from setuptools import setup, Extension


setup(
    name='pymaft',
    version='1.0.0',
    url='https://github.com/rajeshrinet/pymaft',
    author='The PyMAFT team',
    license='MIT',
    description='python library for numerical simulation of field theories',
    long_description='PyGL is a library for numerical simulation of field theories\
    in non-equilibrium statistical mechanics',
    platforms='works on all platforms (such as LINUX, macOS, and Microsoft Windows)',
    ext_modules=cythonize([ Extension("pymaft/*", ["pymaft/*.pyx"],
        include_dirs=[numpy.get_include()],
        )],
        compiler_directives={"language_level": sys.version_info[0]},
        ),
    libraries=[],
    packages=['pymaft'],
    install_requires=['cython','numpy','scipy'],
    extras_require={
            'plotting': ['matplotlib'],
            'notebook': ['jupyter', 'nbconvert']},
    package_data={'pymaft': ['*.pxd']},
)
