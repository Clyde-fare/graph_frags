#!/usr/bin/env python

from setuptools import setup

setup(name='graph_frags',
      version='0.1',
      packages=['graph_frags', 'graph_frags.tests'],
      data_files=[('graph_frags/tests', ['graph_frags/tests/test1.xyz'])],
      )
