Pulsar-Simint
=============

Pulsar-Simint is a library providing Pulsar bindings for the Simint integral
library.

Contents
========
- [Obtaining Pulsar-Simint](#obtaining-pulsar-simint)
- [Building Pulsar-Simint](#building-pulsar-simint)


Obtaining Pulsar-Simint
-----------------------

Pulsar-Simint can be obtained from the official GitHub repository, located
:link:[here](https://github.com/pulsar-chem/Pulsar-Simint).  To obtain the
source use the usual git clone command:

```.sh
git clone https://github.com/pulsar-chem/Pulsar-Simint
```
At the time of writing, Pulsar-Simint will not build Simint for you.  Thus
you will also have to obtain and build Simint.  Simint can be downloaded from
:link:[here](http://www.bennyp.org/research/simint/).  Documentation for how to
build Simint is also located there.

Building Pulsar-Simint
----------------------

Building Pulsar-Simint is similar to building most other Pulsar supermodules.
Configuring is done via the CMake command:

```.sh
cmake -Bbuild -H. -DOPTION1=Value1 -DOPTION2=...
```

Pulsar-Simint tries to honor most CMake variables.  Some influential variables
are:
| Option Name            | Default                 | Description |
|:----------------------:|:-----------------------:|-------------|
| CMAKE_C_COMPILER       | N/A | This is the C compiler to use.  By default CMake will try to find a C compiler on your system. Typically this means it will find  `/bin/gcc`.  |
| CMAKE_CXX_COMPILER     | N/A | Similar to above, except for the C++ compiler. |
| CMAKE_C_FLAGS          | N/A | Any additional flags you want to pass to the C compiler. |
| CMAKE_CXX_FLAGS | N/A | Any additional flags you want to pass to the C++ compiler. |
| CMAKE_BUILD_TYPE | Release | Can be used to enable debug builds.  Valid choices are: `Release`, `Debug`, or `RelWithDebInfo`. |
| CMAKE_PREFIX_PATH | N/A | Used to tell CMake additional places to look for dependencies.  CMake will look for executables in `CMAKE_PREFIX_PATH/bin/`, libraries in `CMAKE_PREFIX_PATH/lib/`, *etc*. |
| CMAKE_INSTALL_PREFIX | `/usr` | The root directory where the final project will be installed following usual GNU conventions.  *i.e.* libraries will be installed to `CMAKE_INSTALL_PREFIX/lib`, header files to `CMAKE_INSTALL_PREFIX/include`, *etc.* |

Additionally, you will need to tell Pulsar-Simint where Simint is located.  To
this end Pulsar-Simint defines the option `SIMINT_PATH` which should point to
the root directory of the Simint installation tree.

After a successful configuration the remainder of the build should be possible
via:
```.sh
cd build && make
make install
```
