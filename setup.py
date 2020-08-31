import numpy
import os, sys, re
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate=True
from setuptools import setup, Extension


with open("README.md", "r") as fh:
    long_description = fh.read()


cwd = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(cwd, 'pygl', '__init__.py')) as fp:
    for line in fp:
        m = re.search(r'^\s*__version__\s*=\s*([\'"])([^\'"]+)\1\s*$', line)
        if m:
            version = m.group(2)
            break
    else:
        raise RuntimeError('Unable to find own __version__ string')


setup(
    name='pygl',
    version=version,
    url='https://github.com/rajeshrinet/pygl',
    author='The PyGL team',
    license='MIT',
    description='PyGL is a numerical library for statistical field theory in Python',
    long_description='PyGL is a numerical library for statistical field theory in Python',
    #long_description=long_description,
    #long_description_content_type='text/markdown',
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
