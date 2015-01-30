#!/usr/bin/env python

from setuptools import setup
setup(name='graph_frags',
      version='0.1',
      packages=['graph_frags'],
      tests_require=[
        'nose'],
      test_suite='nose.collector'
      )
