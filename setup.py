#! /usr/bin/env python

from setuptools import setup

setup(name='Halo tools',
      version='0.1',
      description='Tools to compute properties of Dark Matter halos',
      author='Nicolas Garavito',
      author_email='jngaravitoc@email.arizona.edu',
      install_requieres=['numpy', 'scipy', 'astropy', 'pygadgetreade'],
      packages=['halo'],
     )
