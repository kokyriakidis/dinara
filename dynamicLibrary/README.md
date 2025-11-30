# Directory dinara /dynamicLibrary

This directory builds the Dinara dynamic library `dinara.so`.
This library is used in three ways:

* It is linked in by the dinara dynamic executable `dinaraDynamic`.
* It can be imported by a python script via `import dinara` to provide Dinara Python bindings.
* It can be statically linked in by other C++ code outside Dinara that uses Dinara as a library.
