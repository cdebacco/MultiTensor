# MultiTensor
Multilayer network tensor factorization, for community detection, link prediction and measure layer interdependence.

Copyright (c) 2016 Caterina De Bacco

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Dependencies:
Boost Graph Library: http://www.boost.org/doc/libs/1_62_0/libs/graph/doc/index.html

## How to compile the code:
Multitensor should first be complied modifing appropriately the Makefile included.
You need to specify the file to be compiled by assigning BIN. Example: if you want to compile em.cpp then assign BIN=em .
LIBS and LIBS should point to the directories where you installed boost. You can also arbitrarily modify the compiler's flags under CXXFLAGS.

Then just type from the command line:

make

## How to run run the code:
Type in the command line the name of the binary file (the same you used assigned to BIN) + all the command line options you want to use. Example:

./em -k 2 -l 3 -a "adjacency.dat" -E "_endfile.dat" 


