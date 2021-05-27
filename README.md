![Self-propulsion of active droplets](https://raw.githubusercontent.com/rajeshrinet/pygl/master/examples/banner.jpeg)
## PyGL: statistical field theory in Python [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rajeshrinet/pygl/master?filepath=examples) ![CI](https://github.com/rajeshrinet/pygl/workflows/CI/badge.svg) ![Notebooks](https://github.com/rajeshrinet/pygl/workflows/Notebooks/badge.svg) [![Documentation Status](https://readthedocs.org/projects/pygl/badge/?version=latest)](https://pygl.readthedocs.io/en/latest/?badge=latest)[![PyPI version](https://badge.fury.io/py/pygl.svg)](https://badge.fury.io/py/pygl) [![Downloads](https://pepy.tech/badge/pygl)](https://pepy.tech/project/pygl)

[About](#about) |  [Documentation](https://pygl.readthedocs.io/en/latest/) | [News](#news) | [Installation](#installation) | [Examples](#examples) | [Publications ](#publications)| [Support](#support) | [License](#license)


## About
[PyGL](https://github.com/rajeshrinet/pygl) is a numerical library for statistical field theory in Python. The name GL corresponds to the [Ginzburg–Landau theory](https://en.wikipedia.org/wiki/Ginzburg%E2%80%93Landau_theory). The library has been specifically designed to study field theories without time-reversal symmetry. The library can be used to study models of statistical physics of various symmetries and conservation laws. In particular, we allow models with mass and momentum conservations. The library constructs differentiation matrices using finite-difference and spectral methods. To study the role of momentum conservation, the library also allows computing fluid flow from the solution of the Stokes equation. 

![Self-propulsion of active droplets](https://raw.githubusercontent.com/rajeshrinet/pystokes-misc/master/gallery/droplets/ssi.gif)

The above simulation is done using PyGL. It shows steady-state microphase separation (phase separation arrested to a length-scale) in an active scalar field theory. A self-shearing instability interrupts the growth of droplets by splitting them. Read more: https://arxiv.org/abs/1907.04819

 

## Installation

### From a checkout of this repository

Install PyGL and required dependencies using

```
>> git clone https://github.com/rajeshrinet/pygl.git
>> cd pygl
>> pip install -r requirements.txt
>> python setup.py install
``` 

Install PyGL and its dependencies in a `pygl` environment:

```
>> git clone https://github.com/rajeshrinet/pygl.git
>> cd pygl
>> make env
>> conda activate pygl
>> make
```

### Pip
Alternatively, install the latest PyPI version

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



