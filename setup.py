# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
setup(
    name='bookofnumbers',
    description='Various simple number modules',
    author='Jack Briody',
    author_email='jackbriody@gmail.com',
    license='MIT',
    url='https://github.com/jmbriody/bookofnumbers',
    version="0.1.0",
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    keywords='algebra logic cdnf Quine McCluskey',
    include_package_data=False,
    # install_requires=["pytest"],
    extras_require={
        'test': ['pytest']
    },
)
