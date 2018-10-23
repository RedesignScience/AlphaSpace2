# README

# What is this AlphaSpace?
AlphaSpace is a surface topographical mapping tool.
Currently the use of this software is limited to group member of Yingkai Zhang's lab in NYU

Based on the algorithm of original AlphaSpace published [here](http://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00103), the current 2.0+ version is rewritten, multiple new features are added for a more friendly user interface and performance boost. 

The documentation can be found at [https://lenhsherr.github.io/AlphaSpace/](https://lenhsherr.github.io/AlphaSpace/)

Key features:
1. Clear and pythonic object based API for manipulation of pockets, alpha atoms, and communities. 
2. Integrated trajectory loading and interpretation from **mdtraj**.
3. Native support for multi processing for multi snapshot trajectories. 
4. Integration of Jupyter notebook visualization of protein structure from nglview, and additional visualization of pockets and alpha atoms. 

# Getting started

## Requirements
```python 3.6
numpy
scipy
jupyter notebook
nglview
mdtraj
ipywidgets
```

## Installation
Download the zip package and unzip it in you preferred installation location.
Navigate to the AlphaSpace folder, where you should see this README.md file

You can install the AlphaSpace by entering:
```
pip install .
```
The dependencies should be automatically installed through pip.


### Using jupyter notebook

For interactive session in jupyter notebook, type in
```
 jupyter notebook 
```
To start, you can open the [FCTM_tutorial](examples/1_FCTM.ipynb) in the [examples](examples) folder.

You might also want to checkout [Pocket Comparison](examples/2_Compare.ipynb)