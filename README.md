# README

# What is this AlphaSpace?
AlphaSpace is a surface topographical mapping tool.
Currently the use of this software is limited to group member of Yingkai Zhang's lab in NYU

Based on the algorithm of original AlphaSpace published here[old paper title ](none), the current 2.0+ version is rewritten, multiple new features are added for a more friendly user interface and performance boost. 

Key features:
1. Clear and pythonic object based API for manipulation of pockets, alpha atoms, and communities. 
2. Integrated trajectory loading and interpretation from *mdtraj*.
3. Native support for multi processing for multi snapshot trajectories. 
4. Integration of Jupyter notebook visualization of protein structure from nglview, and additional visualization of pockets and alpha atoms. 

# Getting started

## Dependencies
python 3.6+
jupyter notebook
nglview >=1.0b2
mdtraj 
ipywidgets 7.0

## Installation
Download the zip package and unzip it in you preferred installation location.
Navigate to the AlphaSpace folder, where you should see setup.py. 

Install by entering:
```
pip install -e .
```

The dependencies should be automatically installed through pip. If not, you can install them separately like this. 

```
pip install nglview==1.0.b5
pip install ipywidgets==7
conda install -c conda-forge mdtraj
```

## Activate ipython widget environment
To enable jupyter notebook visualization, you need to enable ipython widget extension. 

```
jupyter-nbextension enable --py --sys-prefix widgetsnbextension
jupyter-nbextension enable --py --sys-prefix nglview
python -m ipykernel install --sys-prefix
```

## Running first AlphaSpace session

For interactive session in jupyter notebook, type in
```
 jupyter notebook 
```
in your terminal to start the notebook.
You can import AlphaSpace using
```
import AlphaSpace
```

You can find a demo notebook in the example folder, which contains the general practice for mapping a protein. 

AlphaSpace also support command line usage, if you run AlphaSpace on a pdb file the mapping will be outputted in the working directory as PDB file and chimera .py file. You will be able to open them in chimera later. Alternatively, a binary pickled file can be saved, this allows for direct loading of results in python.   
To run it in command line mode, do
```
 python alphaspace.py -i [input_file] -o [OPTIONAL:output_directory] -c [OPTIONAL:config_file_path] —-chimera/pickle
```

## Using custom option and parameters
You can use your own config.ini file, by using -c [config file path]
The default config file location is 
`$ALPHASPACE/alphaspace/config.ini`
This file follows the naming convention of V1.0 AS_Param.py and AS_Option.py, they are combined to a .ini file by adding [options] and [parameters]. 



## Changes from V1.0
1. Pocket community definition has been adjusted to be user defined, or based on percentage of the whole surface. 
2. Pocket score are now calculated based on polar and non-polar lining atom SASA ratio, instead of direct number of atoms ratio.  
## Bugs and issue
Known bugs:
SASA calculation may occasionally fail due to internal bug in Shrake algorithms. 

For any feature request or reporting any bugs, please create an issue.