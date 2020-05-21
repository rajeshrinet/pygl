import numpy
import os, sys, os.path, tempfile, subprocess, shutil
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True


setup(
    name='pylandau',
    version='1.0.0',
    url='https://gitlab.com/rajeshrinet/',
    author='Rajesh Singh',
    author_email='rajeshrinet@gmail.com',
    license='MIT',
    description='python library for numerical simulation of fields',
    long_description='pylandau is a library for numerical simulation of fields',
    platforms='tested on LINUX',
    ext_modules=cythonize([ Extension("pylandau/*", ["pylandau/*.pyx"],
        include_dirs=[numpy.get_include()],
        )]),
    libraries=[],
    packages=['pylandau'],
    package_data={'pylandau': ['*.pxd']}
)


