This directory contains C++ library files related to calculation of the Gumbel parameters for pairwise sequence alignment.

Usage with "make".

One way to use this library is with the "make" command. The following assumes you have "make" and a C++ compiler suitably installed. If you use the command line to enter the "cpp" directory and type "make", it should create a library file called "libalp.a". How to use the library is shown in the example directory. If you enter this directory and type "make", it should compile the test program: this will work only if it can find the header and library files. In "example/Makefile", the -I flag to the C preprocessor adds a directory to search for headers ("sls_alignment_evaluer.hpp"), the -L flag to the linker adds a directory to search for libraries, and -lalp specifies linking to "libalp".


Please see the URL
http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/index/software.html#6
for further information.
