import numpy
import os, sys 
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate=True
from setuptools import setup, Extension


setup(
    name='pygl',
    version='2.0.2',
    url='https://github.com/rajeshrinet/pygl',
    author='The PyGL team',
    license='MIT',
    description='PyGL is a numerical library for simulations of field theories in Python.',
    long_description='PyGL is a numerical library for simulations of field theories in Python. \
                      It constructs differentiation matrices using finite-difference and spectral methods. \
                      It also allows to solve Stokes equation using a spectral method, which satisfies compressibility exactly. \
                       The library currently offers support for doing field theoretical simulation and a direct numerical simulation of the Stokes equation \
                       in both two and three space ',
    platforms='tested on LINUX and macOS',
    ext_modules=cythonize([ Extension("pygl/*", ["pygl/*.pyx"],
        include_dirs=[numpy.get_include()],
        )],
        compiler_directives={"language_level": sys.version_info[0]},
        ),
    libraries=[],
    packages=['pygl'],
    install_requires=['cython','numpy','scipy'],
    extras_require={
            'plotting': ['matplotlib'],
            'notebook': ['jupyter', 'nbconvert']},
    package_data={'pygl': ['*.pxd']},
)
