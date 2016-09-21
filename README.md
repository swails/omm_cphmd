# Build Status

[![Linux Build Status](https://travis-ci.org/swails/omm_cphmd.svg?branch=master)](https://travis-ci.org/swails/omm_cphmd.svg)

# What is it?

This project provides a basic C++ library for working with native Amber files.

In particular, it can:

* Create an OpenMM ``System`` object from Amber topology files (note,
  CHAMBER-style topology files are *not* supported)
* Read and write NetCDF trajectory and restart files

# License

This code is distributed under the MIT license. See ``LICENSE`` for more
details.

# How to use it

## Prerequisites

1. You need OpenMM installed in order to build this library (and present in the
   default library search paths).
2. You also need a NetCDF library installed in order to build NetCDF file
   support.

## Building

1. Run ``./configure [options] <compiler>``, where the supported compilers are
   either ``gnu``, ``clang``, or ``intel``. Use the ``--help`` flag for a full
   list of options.
2. Install via the command ``make install``.

This will create a dynamic libamber.so and static libamber.a library that you
can use to link to your own programs. It provides a header file (``Amber.h``)
that needs to be included to provide access to the API.

The entire API lives inside the ``Amber`` namespace.

## Testing

There are a set of unit tests covering most of the code and all of its core
functionality. ``make test`` will run these tests (but you need to run ``make
install`` first).
