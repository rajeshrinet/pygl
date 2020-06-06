## PyMAFT: Field theoretic simulations in python [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/rajeshrinet/pymaft/master?filepath=binder)
[About](#about) | [News](#news) | [Installation](#installation) | [Examples](#examples) | [Publications ](#publications)| [Support](#support) | [License](#license)

## About
[PyMAFT](https://github.com/rajeshrinet/pymaft) is a numerical library for simulations of **M**odels of **A**ctive **F**ield **T**heories in Python.


![Self-shearing instability in active scalar field theory](examples/ssi.gif)

The above simulation is done using PyMAFT. It shows that the droplet growth is interrupted via a self-shearing instability for contractile stress in active model H. Read more: https://arxiv.org/abs/1907.04819

## News
* Our paper has been highlighted in the Journal Club for Condensed Matter Physics with a [commentary](https://doi.org/10.36471/JCCM_March_2020_01).


## Installation
Clone (or download) the repository and use a terminal to install using

```
git clone https://github.com/rajeshrinet/pymaft.git
cd pymaft
python setup.py install
``` 

PyMAFT requires the following software 

## Pip
```
pip install pymaft 
```

- Python 2.6+ or Python 3.4+
- [Cython 0.25.x+](http://docs.cython.org/en/latest/index.html) |  [Matplotlib 2.0.x+](https://matplotlib.org) | [NumPy 1.x+](http://www.numpy.org) | [SciPy 1.1.x+](https://www.scipy.org/) 


## Examples

See the [examples folder](https://github.com/rajeshrinet/pymaft/tree/master/examples) for a list of examples. 

## Publications
* [Hydrodynamically interrupted droplet growth in scalar active matter](https://doi.org/10.1103/PhysRevLett.123.148005). Rajesh Singh and Michael E. Cates. Phys. Rev. Lett. 123, 148005 (2019).

* [Self-propulsion of active droplets without liquid-crystalline order](https://arxiv.org/abs/2004.06064). Rajesh Singh, Elsen Tjhung, and Michael E. Cates. arXiv:2004.06064.  


## Support
Please use the [issue tracker](https://github.com/rajeshrinet/pymaft/issues) on GitHub.

## License
We believe that openness and sharing improves the practice of science and increases the reach of its benefits. This code is released under the [MIT license](http://opensource.org/licenses/MIT). Our choice is guided by the excellent article on [Licensing for the scientist-programmer](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002598). 


