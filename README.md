# MG-evolution

Copyright (c) 2020-2024 Farbod Hassani (MG-evolution) and Julian Adamek (gevolution)
(University of Oslo -- University of Zurich)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
  
**The Software is provided "as is", without warranty of any kind, expressed or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the Software or the use or other dealings in the
Software.**

## Compilation and usage

Before compilation, make sure that all required external libraries are
installed:

* LATfield2 [version 1.1](https://github.com/daverio/LATfield2.git)
* FFTW version 3
* GNU Scientific Library (GSL) including CBLAS
* HDF5

Make sure that the include paths are set properly, or add them to the
makefile. Also check the compiler settings in the makefile. The code is
compiled by typing:

    make

A typical command to run a simulation looks like this:

    mpirun -np 16 ./mgevolution -n 4 -m 4 -s settings.ini

For further information, please refer to the User Manual (manual.pdf)

## Credits

If you use MG-evolution for academic purposes, we kindly ask you to cite

*Farbod Hassani, Lucas Lombriser, Mon.Not.Roy.Astron.Soc. 497 (2020) 2, 1885-1894*
in your articles.

For bug reports and other important feedback you can contact the authors,
farbod.hassani@astro.uio.no (for questions related to MG-evolution)
julian.adamek@uzh.ch (for queries related to gevolution)
developers@latfield.org (for queries related to the LATfield2 library)
