# MultiTensor
Multilayer network tensor factorization, for community detection, link prediction and measure layer interdependence.

Copyright (c) 2016 Caterina De Bacco

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Dependencies:
Boost Graph Library: http://www.boost.org/doc/libs/1_62_0/libs/graph/doc/index.html


## Requirements:
Need to make a directory called `data` inside the folde where the `.cpp`  and `Makefile` are stored. Just ype from the command line, inside this folder: 
`mkdir data`

## How to compile the code:
Multitensor should first be complied modifing appropriately the Makefile included.
You need to specify the file to be compiled by assigning BIN. Example: if you want to compile em.cpp then assign BIN=em .
LIBS and LIBS should point to the directories where you installed boost. You can also arbitrarily modify the compiler's flags under CXXFLAGS.

Then just type from the command line:

`make`

## How to run run the code:
Type in the command line the name of the binary file (the same you used assigned to BIN) + all the command line options you want to use. Example:

`./em -k 2 -l 3 -a "adjacency.dat" -E "_endfile.dat" `

# Required arguments

`-a` : Adjacency matrix file
`-f` : Folder where the adjacency input and output are/will be stored (inside `data` folder).

# Optional arguments

`-E` : Output end of file where the paramters' files will be stored. Example: `-E="_abc.dat" ` output files will be `u_K4_abc.dat`,`v_K4_abc.dat`,`w_K4_abc.dat` (assuming that k=4). Default value is `-E=".dat"`.
`-i` : Initialization flag: if `i=0` than parametrs are randomly initialized; if `i=1` the membership vectors u and v and w are initialized form file; if `i=2` only w is initialized from file; if `i=3` only u and v are initialized from file, w instead is random.
`-w` : End of the file where the parameters can be initialized from, in case initialization variable is greater than 0.
`-l` : Number of layers, default is 1.
`-k` : Number of communities, default is 4.
`-r` : Number of different realizations, the final parameters will be the one correspondinf to the realization leading to the max likelihood. Default is 1.
`-t` : Max iteration time. Default is 500.
`-e` : Convergence tolerance. Default is 0.1 .
`-y` : Decision variable for convergence. Default is 10.
`-z` : Seed for random real numbers.
`-s` : Seed for random integer numbers.
