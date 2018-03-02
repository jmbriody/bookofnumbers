# -*- coding: utf-8 -*-
import os
from setuptools import setup, find_packages, Command

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info ./.tox')

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
    cmdclass={
        'clean': CleanCommand,
    },
)
