# MultiTensor: Python code
Python version.
Copyright (c) 2016 Caterina De Bacco

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Dependencies:
Needs three main Python modules to be downloaded:

* `numpy` : https://docs.scipy.org/doc/numpy-1.10.1/user/install.html
* `networkx` : https://networkx.github.io/documentation/development/install.html
* `argparse` : https://pypi.python.org/pypi/argparse

## What's included:
- `main.py` : General version of the algorithm. Considers both directed and undirected weigthed multilayer networks with any community structures (non-diagonal or restricted diagonal affinity matrices W).
- `MultiTensor.py` : Contains the class definition of a Multilayer network with all the member functions required.
- `tools.py` : Contains non-class functions.

Use the version that most resembles your network, i.e. if you have an undirected network set the flag '-u=1'. If you also know that the partition is assortative then use the flag '-A=1'.

## Requirements:
Need to make a directory called `data` outside the folder where the `python` folder is stored. Just ype from the command line, inside that folder: 
* `mkdir data`

## How to compile the code:
The Python code does not need to be compiled (but is slower than the compiled c++ version).
Just run it!

## How to run run the code:
Type in the command line the name of the binary file (the same you used assigned to BIN) + all the command line options you want to use. Example:

`python main.py -k=2 -l=4 -a="adjacency.dat" -E="_endfile.dat" `

As an example type on the command line:

`python main.py`

It will use the sample adjacency file contained in `../data`, which is an undirected, unweighted network with `L=4` and `K=5`. 

### Required arguments

- `-a` : Adjacency matrix file
- `-f` : Folder where the adjacency input and output are/will be stored (inside `data` folder).

### Optional arguments

- `-E` : Output end of file where the paramters' files will be stored. Example: `-E="_abc.dat" ` output files will be `u_K4_abc.dat`,`v_K4_abc.dat`,`w_K4_abc.dat` (assuming that k=4). Default value is `-E=".dat"`.
- `-i` : Initialization flag: if `i=0` than parametrs are randomly initialized; if `i=1` the membership vectors u and v and w are initialized form file; if `i=2` only w is initialized from file; if `i=3` only u and v are initialized from file, w instead is random. Default is `i=0`.

* `-w` : End of the file where the parameters can be initialized from, in case initialization variable is greater than 0.

* `-l` : Number of layers, default is 4.
* `-k` : Number of communities, default is 5.
* `-r` : Number of different realizations, the final parameters will be the one corresponding to the realization leading to the max likelihood. Default is 1.
* `-t` : Max iteration time. Default is 500.
* `-e` : Convergence tolerance. Default is 0.1 .
* `-g` : Error added when intializing parameters from file. Default is 0.1 .
* `-o` : Flag to output adjacency matrix. Default is 0 (False).
* `-y` : Decision variable for convergence. Default is 2.
* `-z` : Seed for random real numbers.
* `-A` : Flag to call the (faster) restricted assortative version (purely diagonal affinity matrix).
* `-u` : Flag to call the undirected network, default is 0 (False).

## Input format.
The multilayer adjacency matrix should be formatted as an edge list with L+3 columns:

`E node1 node2 3 0 0 1`

The first columns tells the algorithm that the row denotes an edge; the second and third are the source and target nodes of that edge, respectively; l+3 column tells if there is that edge in the l-th layer and the weigth (must be integer). In this example the edge node1 --> node2 exists in layer 1 with weight 3 and in layer 4 with weight 1, but not in layer 2 and 3.

## Output.
Three files will be generated: the two NxK membership matrices `U` and `V`, and the KxK layer affinity matrix `W`. Supposing that K=4 and `E=".dat"` the output files will be inside `data` folder with names:
- `u_K4.dat`
- `v_K4.dat`
- `w_K4.dat`

The first line outputs the Max Likelihood among the realizations.
For the membership files, the subsequent lines contain L+1 columns: the first one is the node label, the follwing ones are the (not normalized) membership vectors' entries.
For the affinity matrix file, the subsequent lines start with the number of the layer and then the matrix for that layer.
For the restricted assortative version only the diagonal entries of the affinity matrix are printed. The first entry of each row is the layer index.



