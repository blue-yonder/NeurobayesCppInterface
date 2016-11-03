# NeuroBayes C++ Interface

Interface to NeuroBayes(R) for academic purpose

The NeuroBayes package is available under the NeuroBayes license from Blue Yonder GmbH.
Please contact neurobayes@blue-yonder.com for general information and for information about the license.

This repository provides the C++ Interface of the NeuroBayes package. You can only build this project if you
obtained the NeuroBayes licence and NeuroBayes Core library from Blue Yonder GmbH.

The NeuroBayes C++ Interface is usually already shipped with the NeuroBayes package.
The NeuroBayes Core library does only depenend on a very small set of system libraries like libc and libgfortran,
on the other hand the NeuroBayes C++ Interface depends on libstdc++.
It may be necessary to recompile the C++ Interface in the future to support new platforms,
hence Blue Yonder GmbH decided to release the source code of the C++ Interface under the MIT licence for the convinience
of the NeuroBayes user community (in particular for high energy physics experiments).


# Build process and installation
  * Extract the NeuroBayes Core shared library named libNeuroBayesCore\_shared.so from the package you obtained from Blue Yonder GmbH and copy it into the previously empty directory core.
  * cd build/
  * cmake ../
  * make


## Building rpm, deb or tgz packages
Optionally you can build rpm, deb and tgz packages using
  * make package
Additional third-party libraries are required to do so.


## Maxnode
The NeuroBayes Core library is optimized for a maximum number of nodes.
The maximum number of nodes is printed by the NeuroBayes Teacher
at the start of the program e.g. execute the minimaltest executable and search for NB\_MAXNODE.
The default number is 100. If your library is optimized for a different number of nodes,
you have to provide the correct number using an environment variable named MAXNODE\_PREPRO.
during the build process e.g. MAXNODE\_PREPRO=200 cmake ../


## Tests
To verify the correctness of your executable, you can execute the following test suits:
  * cd build
  * ./minimaltest
  * ./interfacetest

The minimaltest provides also a minimal example to get you started.
