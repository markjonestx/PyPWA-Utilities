PWA Utilities
=============

This package is a collection of various utilities for use alongside PyPWA.
The utilities were pulled from pwa2000 and ppgen along with just the source
and header files needed to compile these resources.

This package includes:
* GAMP
* HGAMP
* VAMP
* PPGEN


Compiling and Installing
------------------------

Create a build directory either in the repository or outside of the repository:
```
$ mkdir build
$ cd build
```

Run cmake pointing to the original repository:
```
$ cmake ../path/to/source/repo
```

Or optionally with an install prefix:
```
$ cmake -DCMAKE_INSTALL_PREFIX=/home/user/.local/ ../path/to/source/repo
```

Then run make and install:
```
$ make
$ make install
```
