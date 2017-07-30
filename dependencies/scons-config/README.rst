============
scons-config
============

Description
===========

While SCons is an excellent build system, I find that its configuration scheme
is a little lacking. For packages that have complex setups, with variations
across platforms and between versions, configuration can be a very difficult
process. scons-config is intended to reduce the pain of locating and validating
installed packages on a machine.

Dependencies
============

There are two dependencies.

  Python 2.6 or 2.7
    Good old Python, version 2.6 or 2.7.

  SCons
    A ``make`` alternative. Download at "http://www.scons.org".

Installation
============

Run the ``setup.py`` script as per distutils installation::

  sudo python setup.py install

Example SConstruct
==================

There is an example ``SConstruct`` file in ``sconsconfig``, it is intended to be
used as a template in your projects. Hopefully it's not too obfuscated and can
be easily understood.

Available Packages
==================

The packages available for configuration can be found under
``sconsconfig/packages``. They can be added to the SConstruct file in the same
way as demonstrated in the example.

Download Packages
=================

Packages can be setup to be downloaded from a source and automatically built
according to the needs of your project. TODO: Needs documentation.
