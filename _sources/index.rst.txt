.. Tadbit documentation master file, created by
   sphinx-quickstart on Tue Jan 15 18:23:49 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.




Welcome to TADbit’s documentation!
==================================

TADbit is a computational package that deals with 3C-based data (more
specifically Hi-C data). It handles all the steps from the alignment of
paired-end reads to the detection of Topologically Associating Domain
(TAD) borders, compartments and three-dimensional modeling of chromatin
based on interaction matrices. The different steps:

-  preprocessing of paired-end reads from Hi-C experiment
-  alignment of the reads using
   `GEM <http://algorithms.cnag.cat/wiki/The_GEM_library>`__
-  fitering of the mapped reads
-  construction of interaction matrices and normalization
-  the identification of TAD borders and compartments
-  the 1D and 2D comparison of TADs
-  the three-dimensional modeling of the chromatin

The library has been designed to be used by researchers with minimal
level of expertise in computer science. A set of the all-in-one scripts
(called **TADbit tools**) allow to run a full analysis in one single
command line.

On the other hand, advanced users may produce their own scripts using
the TADbit API.

Documentation index:
~~~~~~~~~~~~~~~~~~~~

.. epigraph::

  .. toctree::
     :maxdepth: 3
  
     install
  
     tutorial
     
     tools
  
     reference/index
  
     biblio

