"""
Project: 
File: setup

Created by Scrat on 11.05.2017
"""

from setuptools import setup


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='pySTARMA',
    version='0.1.1',
    description='',
    long_description=readme,
    author='Andreas Wolf',
    author_email='andreas.wolf.ke@gmail.com',
    url='https://github.com/scrat-online/pySTARMA.git',
    license=license,
    packages=['pySTARMA'],
    install_requires=['numpy', 'pandas', 'prettytable'],
)