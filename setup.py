'''
Sets up exoVista_wrapper, a package to make it easier to interface with
exoVista data using python
'''

from setuptools import setup, find_packages
setup(
    name='exoVista_wrapper',
    version='0.00001',
    packages=find_packages(include=['exoVista_wrapper', 'exoVista_wrapper.*']),
    install_requires=[
        'astropy',
        'matplotlib',
        'numpy',
        'pandas',
        'tqdm',
    ]
)
