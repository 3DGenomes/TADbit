

+-----------------------+-+
|                       | |
| Current version: v0.1 | |
|                       | |
+-----------------------+-+


TADdyn is a Python library that allows to model and explore single or time-series 3C-based data. 
These datasets are constituted by interaction matrices that describe distinct stages of naturally 
occurring or induced cellular process such as the cell trans-differentiation, or reprogramming. 
With TADdyn the user can load at once the raw and normalised interaction binned matrices (Hi-C like matrices) 
at each of the experimental stages, build 4D models, and finally, extract structural properties from the models. 
The 4D models reproduce the expected interaction patterns at the experimental time-points, 
and also describe the structural modifications at intermediate moments (between stages) under the hypothesis 
that the changes occurring between consecutive experimental time-points are smooth. To do this, 
TADdyn is designed as a combination of restraint-based modelling, and steered Langevin dynamics of Physics-based 
chromatin models. 

Documentation
*************

**Install LAMMPS as a shared library**
   | 1 - Download lammps
   | git clone -b stable https://github.com/lammps/lammps.git mylammps
   
   | 2 - Download the colvar modified version
   | git clone https://github.com/david-castillo/colvars

   | 3 - Update the user-defined colvars library
   | ./colvars/update-colvars-code.sh ./mylammps/

   | 4 - Compile colvars library
   | cd ./mylammps/lib/colvars
   | make -f Makefile.g++

   | 5 - Install lammps as a shared library
   | cd ../../src/
   | make yes-user-colvars
   | make yes-molecule
   | make serial mode=shlib
   | make yes-python
   
   | export LD_LIBRARY_PATH="/complete-path-to-mylammps/mylammps/src/"
   
   | cd ../../

**Install packages**
   | conda update python
   | conda install -c r r-base
   | conda install -y scipy           # scientific computing in python
   | conda install -y numpy           # scientific computing in python
   | conda install -y matplotlib      # to produce plots
   | conda install -y jupyter         # this notebook :)
   | conda install -y -c https://conda.anaconda.org/bcbio pysam # to deal with SAM/BAM files
   | conda install -y -c https://conda.anaconda.org/salilab imp # for 3D modeling
   | conda install -y pip             # yet another python package manager
   | conda install -y -c bioconda mcl # for clustering
   | conda uninstall Pebble
   | pip install Pebble==4.3.10

**Install TADdyn**
   | 1 - Download TADdyn from the Github repository
   | git clone https://github.com/david-castillo/TADbit.git -b TADdyn TADdyn

   | 2 - Install TADdyn
   | cd TADdyn
   | python setup.py install
   | cd ..

**Try TADdyn**
   | cd test/Sox2
   | python test_TADdyn_on_Sox2.py
   
Citation
********
Please, cite this article if you use TADdyn.

Marco Di Stefano, Ralph Stadhouders, Irene Farabella, David Castillo, François Serra, Thomas Graf, Marc A. Marti-Renom.
**Dynamic simulations of transcriptional control during cell reprogramming reveal spatial chromatin caging.**
*bioRxiv* 642009; `doi: https://doi.org/10.1101/642009`_

Methods implemented in TADdyn
-----------------------------
In the actual implementation, TADdyn relies on TADbit for the preparation of the 3C-based datasets from mapping to normalization, 
and on LAMMPS [Plimpton]_ for the implementation of the simulations.

Bibliography
************

.. [Plimpton] Plimpton, S. Fast Parallel Algorithms for Short-Range Molecular Dynamics. J Comp Phys 117, 1-19 (1995) and Fiorin, G., Klein, M.L. & Hénin, J. Using collective variables to drive molecular dynamics simulations. Molecular Physics 111, 3345-3362 (2013).