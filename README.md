# pturbdb

Parallel Turbulence Database Framework

Copyright (c) 2018-2019 Jason Graham <jgraham@compukix.net>

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.

## Overview

Contained here is a working prototype for the Parallel Turbulence Database
(pturbdb) library. The goal of pturbdb is to provide a high-level framework for
constructing a database from a sequence of temporally ordered fields (the target
fields are turbulent fields, but can really be anything) from which complex
analyses may be easily perform on a desktop machine or a large-scale distruted
system.

At the core of pturbdb is the parallel fields class (PField) which provides
support for operating on distributed arrays. The idea is to have constructs
which may perform operations from basic arithmetic to spatial derivatives (and
others) while taking care of the parallel complexity and removing it from the
user application. This class may also be used as the base construct other
applications and can serve as the core for a full PDE solver.

pturbdb currently uses HDF5 as the data reader/writer and performs all I/O
operations in parallel using MPI-IO. By adding other data reader/writer the
library may also be used for other database types. Currently data caching is
provided by the FieldCache template library to minimize the I/O time needed for
the PCHIP temporal interpolation.

A couple of notes and things to consider:

 - Since it is a working prototype things may not be fully optimized
   or may not necessarily have the best design in some places (have
   some room to grow.)

 - Since the target usage was for constructing and analyzing database
   fields, which is an I/O bound problem, most CPU optimizations have
   not been implemented.

 - The PField class should probably be reimplented as a class template
   which will give much more flexibility. This would allow not only
   its use with various standard C++ types, but also things like
   Blitz++, Eigen, etc. vectors or arrays may also be used. For
   example, the pfield_math module (see pfield_math.hpp) uses a vector
   of fields. I believe it would be better to have a field of
   vectors. Complex data is also handled similarly.

## Build and Install

The full installation instructions can be found in INSTALL. Autotools
is used as the build system, so generally

    % ./configure
    % make
    % make install

will be all that is needed. See *./configure --help* for additional
information.
