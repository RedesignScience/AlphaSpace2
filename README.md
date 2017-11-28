# README

# What is this AlphaSpace?
AlphaSpace is a surface topographical mapping tool.
Currently the use of this software is limited to group member of Yingkai Zhang's lab in NYU

Based on the algorithm of original AlphaSpace published [here](http://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00103), the current 2.0+ version is rewritten, multiple new features are added for a more friendly user interface and performance boost. 

Key features:
1. Clear and pythonic object based API for manipulation of pockets, alpha atoms, and communities. 
2. Integrated trajectory loading and interpretation from **mdtraj**.
3. Native support for multi processing for multi snapshot trajectories. 
4. Integration of Jupyter notebook visualization of protein structure from nglview, and additional visualization of pockets and alpha atoms. 

# Getting started

## Dependencies
```python 3.6+

jupyter notebook

nglview>=1.0

mdtraj 

ipywidgets 7.0
```

## Installation
Download the zip package and unzip it in you preferred installation location.
Navigate to the AlphaSpace folder, where you should see this README.md file

Install by entering:
```
pip install -e .
```
The dependencies should be automatically installed through pip.



## Activate ipython widget environment
To enable jupyter notebook visualization, you may need to enable ipython widget extension. 

```
jupyter-nbextension enable --py --sys-prefix widgetsnbextension
jupyter-nbextension enable --py --sys-prefix nglview
python -m ipykernel install --sys-prefix
```

## Running first AlphaSpace session

### Using jupyter notebook

For interactive session in jupyter notebook, type in
```
 jupyter notebook 
```
To start, you can open the [FCTM_tutorial](examples/FCTM_tutorial.ipynb) in the [examples](examples) folder.

### Using commandline __experimental__
AlphaSpace also support command line usage, if you run AlphaSpace on a pdb file the mapping will be outputted in the working directory as PDB file and chimera .py file. You will be able to open them in chimera later. Alternatively, a binary pickled file can be saved, this allows for direct loading of results in python.   
To run it in command line mode, do
```
 python run_alphaspace.py -i [input_file] -o [output_file] -c [OPTIONAL:config_file_path] —-chimera/pickle
```


## Changes from V1.0
~~1. Pocket community definition has been adjusted to be user defined, or based on percentage of the whole surface.~~ 
2. Pocket score are now calculated based on polar and non-polar lining atom SASA ratio, instead of direct number of atoms ratio.  
## Bugs and issue
Known bugs:
SASA calculation may occasionally fail due to internal bug in Shrake algorithms. 

For any feature request or reporting any bugs, please create an issue.