#!/usr/bin/env python

from setuptools import setup, find_packages

version = '0.1dev'

print '''------------------------------
Installing nfRnaReport version {}
------------------------------
'''.format(version)


setup(
    name='nfRnaReport',
    version=version,
    author='lx Gui',
    author_email='guilixuan@gmail.com',
    keywords=['bioinformatics', 'NGS', 'QC'],
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    scripts=['scripts/nfRnaReport', ],
    install_requires=[
        'Jinja2',
        'configparser',
    ]

)


print '''------------------------------
nfRnaReport installation complete!
------------------------------
'''
