from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='alphaspace2',
      version='0.1.2',
      description='Protein topographical mapping tool',
      long_description="""
      AlphaSpace2 is a surface topographical mapping tool.
        Based on the algorithm of original AlphaSpace published
        [here](http://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00103),
        the alphaspace2 is a rewritten implementation with multiple new features
        are added for a more friendly user interface and performance boost.
      """,
      url='http://github.com/lenhsherr/alphaspace2',
      author='Haotian Li',
      author_email='hl2368@nyu.edu',
      license='gpl-3.0',
      packages=['alphaspace2'],
      scripts=['bin/alphaspace2'],
      include_package_data=True,
      zip_safe=False,

      install_requires=['mdtraj',
                        'cython',
                        'configargparser',
                        'scipy',
                        ],
      classifiers=[
          "Programming Language :: Python :: 3",
          "Operating System :: OS Independent",
      ]

      )
