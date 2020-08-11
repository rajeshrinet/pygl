![Self-propulsion of active droplets](https://raw.githubusercontent.com/rajeshrinet/pygl/master/examples/banner.jpeg)
## PyGL: field theories without time-reversal symmetry in Python [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rajeshrinet/pygl/master?filepath=examples) ![CI](https://github.com/rajeshrinet/pygl/workflows/CI/badge.svg) ![Notebooks](https://github.com/rajeshrinet/pygl/workflows/Notebooks/badge.svg) [![Documentation Status](https://readthedocs.org/projects/pygl/badge/?version=latest)](https://pygl.readthedocs.io/en/latest/?badge=latest)[![PyPI version](https://badge.fury.io/py/pygl.svg)](https://badge.fury.io/py/pygl) [![Downloads](https://pepy.tech/badge/pygl)](https://pepy.tech/project/pygl)

[About](#about) |  [Documentation](https://pygl.readthedocs.io/en/latest/) | [News](#news) | [Installation](#installation) | [Examples](#examples) | [Publications ](#publications)| [Support](#support) | [License](#license)


## About
[PyGL](https://github.com/rajeshrinet/pygl) is a numerical library for simulations of field theories in Python. GL corresponds to the [Ginzburg–Landau theory](https://en.wikipedia.org/wiki/Ginzburg%E2%80%93Landau_theory). The library constructs differentiation matrices using finite-difference and spectral methods. The library also allows to compute fluid flow on a 2D or 3D grid from the solution of the Stokes equation using a spectral solver. 
 

## Installation
Clone (or download) the repository and use a terminal to install using

```
>> git clone https://github.com/rajeshrinet/pygl.git
>> cd pygl
>> python setup.py install
``` 

PyGL requires the following software 


- Python 2.6+ or Python 3.4+
- [Cython 0.25.x+](http://docs.cython.org/en/latest/index.html) |  [Matplotlib 2.0.x+](https://matplotlib.org) | [NumPy 1.x+](http://www.numpy.org) | [SciPy 1.1.x+](https://www.scipy.org/) 

## Pip
Alternatively install latest PyPI version

```
>> pip install pygl 
```


## Examples

See the [examples folder](https://github.com/rajeshrinet/pygl/tree/master/examples) for a list of examples. 

## Publications
* [Hydrodynamically interrupted droplet growth in scalar active matter](https://doi.org/10.1103/PhysRevLett.123.148005). Rajesh Singh and Michael E. Cates. Phys. Rev. Lett. 123, 148005 (2019).

* [Self-propulsion of active droplets without liquid-crystalline order](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.032024). Rajesh Singh, Elsen Tjhung, and Michael E. Cates. Phys. Rev. Research 2, 032024(R) (2020).


## News
* Our paper has been highlighted in the Journal Club for Condensed Matter Physics with a [commentary](https://doi.org/10.36471/JCCM_March_2020_01).


## Support
Please use the [issue tracker](https://github.com/rajeshrinet/pygl/issues) on GitHub.

## License
We believe that openness and sharing improves the practice of science and increases the reach of its benefits. This code is released under the [MIT license](http://opensource.org/licenses/MIT). Our choice is guided by the excellent article on [Licensing for the scientist-programmer](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002598). 


