# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
setup(
    name='bookofnumbers',
    description='Various simple number modules',
    author='Jack Briody'
    author_email='jackbriody@gmail.com',
    url='https://jmbriody.github.io'
    version="0.1.0",
    packages=find_packages(),
    include_package_data=False,
    install_requires=["pytest"],
)
