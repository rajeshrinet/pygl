import numpy
import os, sys 
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate=True
from setuptools import setup, Extension


setup(
    name='pymaft',
    version='1.0.1',
    url='https://github.com/rajeshrinet/pymaft',
    author='The PyMAFT team',
    license='MIT',
    description='PyMAFT is a numerical library for simulations of Models of Active Field Theories in Python.',
    long_description='PyMAFT is a numerical library for simulations of Models of Active Field Theories in Python.\
                      It constructs differentiation matrices using finite-difference and spectral methods. \
                      It also allows to solve Stokes equation using a spectral method, which satisfies compressibility exactly. \
                       The library currently offers support for doing field theoretical simulation and a direct numerical simulation of the Stokes equation \
                       in both two and three space ',
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
