# README

# What is this AlphaSpace?
AlphaSpace is a surface topographical mapping tool.
Currently the use of this software is limited to group member of Yingkai Zhang's lab in NYU

# Getting started
## Python version
3.6+

## Dependencies
jupyter notebook
nglview >=1.0
mdtraj 
ipywidgets 7.0

## Installation
Download the zip package and unzip it in you preferred installation location.
Navigate to the AlphaSpace folder, where you should see setup.py. 

Install by entering:
```
pip install -e .
```

The dependencies should be automatically installed through pip. 

```
pip install nglview==1.0.b5
pip install colour
pip install ipywidgets==7
conda upgrade --all
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
 `jupyter notebook` 
in your terminal to start the notebook.
You can import AlphaSpace using
 `import AlphaSpace `

To run it in command line mode, do
` python alphaspace.py -i [input_file] -o [output_directory] -c [config_file_path]`

## Using custom option and parameters
You can use your own config.ini file, by using -c [config file path]
The default config file location is 
`$ALPHASPACE/alphaspace/config.ini`

## Bugs and issue

For any feature request or reporting any bugs, please create an issue.