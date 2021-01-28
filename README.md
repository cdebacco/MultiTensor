# MultiTensor
Multilayer network tensor factorization, for community detection, link prediction and measure layer interdependence. 

**New version**: a new updated and efficent `cpp` and `python`, version can be found at the [github/MPI-IS/multitensor](https://github.com/MPI-IS/multitensor). At this link you can find the documentation and usage example. 

The codes in this repository are therefore no longer maintained, all the future new updates will be uploaded to the new repository [github/MPI-IS/multitensor](https://github.com/MPI-IS/multitensor).

Implements the algorithm described in:

[1] De Bacco, C., Power, E. A., Larremore, D. B., & Moore, C. (2017). *Community detection, link prediction, and layer interdependence in multilayer networks.* Physical Review E, 95(4), 042317.

If you use this code please cite [[1]](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.042317).  
A _preprint_ version can be found [here](http://cdebacco.com/files/multitensor.pdf) or [here](https://arxiv.org/abs/1701.01369).

If you are further interested in MultiTensor **extensions**:
  - incorporating **node attributes**, we have implemented [`MTCOV`](https://github.com/mcontisc/MTCOV).  
  [[2]](https://www.nature.com/articles/s41598-020-72626-y) Contisciani M., Power E. & De Bacco C. (2020). _Community detection with node attributes in multilayer networks_, Scientific Reports 10, 15736 (2020).
  - incorporating **reciprocity**,  we have implemented [`CRep`](https://github.com/mcontisc/CRep).   
  [[3]](https://arxiv.org/abs/2012.08215) Safdari H., Contisciani M. & De Bacco C. (2020). _A generative model for reciprocity and community detection in networks_, arXiv:2012.08215.

Copyright (c) 2016 [Caterina De Bacco](http://cdebacco.com/).

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
