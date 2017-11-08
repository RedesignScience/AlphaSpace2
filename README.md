# AlphaSpace readme
# What is this AlphaSpace?
AlphaSpace is a surface topographical mapping tool.
Currently the use of this software is limited to group member of Yingkai Zhang's lab in NYU

# Getting started
##Python version
3.6+

## Dependencies
jupyter notebook
nglview 1.0
colour
mdtraj
ipywidgets 7.0

## Installation

We suggest using conda and pip wto install packages.

```
pip install nglview==1.0.b5
pip install colour
pip install ipywidgets==7
conda upgrade --all
conda install -c conda-forge mdtraj
```

## Activate ipython widget environment

```
jupyter-nbextension enable --py --sys-prefix widgetsnbextension
jupyter-nbextension enable --py --sys-prefix nglview
python -m ipykernel install --sys-prefix
```

## Running first AlphaSpace session

Type in `jupyter notebook` in your terminal to start the notebook.
You can import AlphaSpace using
 `from alphaspace import *`
Make sure the AS folder is in your script directory

## Bugs and issue

For any feature request or reporting any bugs, please create an issue.