#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import sys
sys.path.append('biofrost')
import codecs

from pathlib import Path
from setuptools import setup, find_packages

from biofrost.version import __version__


def read(infile):
    return codecs.open(Path(__file__).parent / infile).read()


setup(
    name='biofrost',
    python_requires='>3.7.0',
    version=__version__,
    url='https://github.com/Kevinzjy/biofrost',
    description='The magic bifrost bridge of bioinformatic',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    author='Jinyang Zhang',
    author_email='zhangjinyang@biols.ac.cn',
    maintainer='Jinyang Zhang',
    maintainer_email='zhangjinyang@biols.ac.cn',
    license='BSD-3',
    keywords='bioinformatics',
    packages=find_packages(exclude=['docs', 'tests']),
    entry_points={
        'console_scripts': [
            'biofrost=biofrost.main:main',
        ]
    },
    include_package_data=True,
    zip_safe=True,
    install_requires=[
        'argparse>=1.2.1',
        'matplotlib>=3.5.1',
        'pandas>=1.1.5',
        'scikit-learn>=1.0.2',
        'pysam>=0.19.1',
    ],
    extras_require={
        'bert': ['bert-serving-server>=1.8.6', 'bert-serving-client>=1.8.6', 'pytorch-transformer', 'flair'],
        'vision': ['opencv-python>=4.0.0', 'imagehash>=4.0', 'image', 'peakutils'],
    },
    test_suite="nose.collector",
    tests_require=['nose==1.3.7'],
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
    ],
)