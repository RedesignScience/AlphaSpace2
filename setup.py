from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='alphaspace',
      version='2.0',
      description='Protein topographical mapping tool',
      url='http://github.com/lenhsherr/alphaspace',
      author='Haotian Li',
      author_email='hl2368@nyu.edu',
      license='MIT',
      packages=['alphaspace'],

      package_data={
      	'': ['*.txt', '*.md'],
        # And include any *.msg files found in the 'hello' package, too:
        'alphaspace.tests.bcl2.lig': ['*.pdb'],
        'alphaspace.tests.bcl2.prot': ['*.pdb'],
          'Examples.mdm2_p53': ['*.pdb'],

    },

      install_requires={
          'numpy',
          'scipy',
          'cython',
          'jupyter',
          'nglview',
          'mdtraj',
          'ipywidgets',
          'networkx',

      },
      include_package_data=True,
      zip_safe=False,

      scripts=['scripts/run_alphaspace.py'],



      )