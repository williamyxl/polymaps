from setuptools import setup

setup(
    name='polymaps',
    version='0.0.22',    
    description='A simple Python package to preprocess and postprocess LAMMPS simulations for polymers and small molecules',
    url='https://gitlab.com/doublylinkedlist/polymaps/',
    author='Xiaoli Yan',
    author_email='xyan11@uic.edu',
    license='MIT',
    packages=['polymaps'],
    install_requires=['jupyter==1.0.0',
                      'jupyterlab==4.0.4',
                      'numpy==1.25.0',
                      'pandas==2.0.2',
                      'matplotlib==3.7.1',
                      'plotly==5.15.0'
                      ],

    
)

