# README

# What is this AlphaSpace2?
AlphaSpace2 is a surface topographical mapping tool.
Based on the algorithm of original AlphaSpace published
[here](http://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00103),
the alphaspace2 is a rewritten implementation with multiple new features
are added for a more friendly user interface and performance boost.

The documentation can be found at

Key features:
1. Clear and pythonic object based API for manipulation of pockets, alpha atoms, and communities. 
2. Integrated trajectory loading and interpretation from **mdtraj**.
3. Create Chimera script for publication quality visualization
4. Provide new beta score for beta atoms

# Getting started

## Requirements
```
python
cython
numpy
scipy
mdtraj
```

## Installation
You can install the package with pip command:

```
pip install alphaspace2
```

To install the latest version

```
git clone https://github.com/lenhsherr/AlphaSpace2.git\
cd AlphaSpace2
pip install -e ./
```

