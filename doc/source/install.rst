Install Tadbit on GNU/Linux
***************************

.. contents::

*Note: for now this section is only relative to installation under Debian-linux*


Tadbit requires python>=2.7 as well as several dependencies that are listed below.

Dependencies
============

Python libraries
----------------

**Required:**
::

  apt-get install python-scipy
  apt-get install python-numpy

Accessory (but **highly** recommended):

::

  apt-get install python-matplotlib

IMP - 3D modelling
------------------

Determination of the three dimensional structure of a given locus is done with the IMP package.

In particular, Tadbit works with version 2.0.1 (but may also work with other versions).

Installation instruction for IMP can be found at their website (http://salilab.org/imp/nightly/doc/html/ ). 

**Also we do not provide support for IMP installation**, we make available here the minimal steps to install IMP on Ubuntu  (12.04 and 13.04):

::

    sudo apt-get install cmake
    sudo apt-get install libboost1.49-all-dev
    sudo apt-get install libhdf5-dev
    sudo apt-get install swig
    sudo apt-get install libcgal-dev
    sudo apt-get install python-dev


Download tarball from http://salilab.org/imp/ and uncompress:

::

   wget http://salilab.org/imp/get.php?pkg=2.0.1/download/imp-2.0.1.tar.gz -O imp-2.0.1.tar.gz
   tar xzvf imp-2.0.1.tar.gz
   cd imp-2.0.1

compile (*Note:* the `-j` option stands for the number of CPUs you want to assign to the task, the higher the faster).

::

   cmake . -DCMAKE_BUILD_TYPE=Release -DIMP_MAX_CHECKS=NONE -DIMP_MAX_LOG=SILENT
   make -j4 

Once all is done, open the file setup_environment.sh in your imp directory and copy the first lines to your `~/.bashrc` file. These lines should look like:

::

  LD_LIBRARY_PATH="/SOMETHING/imp-2.0.1/lib:/SOMETHING/imp-2.0.1/lib:/SOMETHING/imp-2.0.1/src/dependency/RMF/:$LD_LIBRARY_PATH"
  export LD_LIBRARY_PATH

  PYTHONPATH="/SOMETHING/imp-2.0.1/lib:/SOMETHING/imp-2.0.1/lib:/SOMETHING/imp-2.0.1/src/dependency/RMF/:$PYTHONPATH"
  export PYTHONPATH


**Important note**

  Do not copy the lines above, copy them from setup_environment.sh, where ``SOMETHING`` is replaced by your real path to IMP.


MCL - clustering
----------------

MCL program is used for clustering 3-dimensional models can be downloaded from http://micans.org/mcl/ .

*Note: in case MCL is not found by tadbit, an alternative clustering method will be used. Nevertheless we strongly recommend to use MCL.*

On Debian/Ubuntu machines you can install MCL from the repositories:

::

  sudo apt-get install mcl


Chimera - visualization
-----------------------

Chimera is a program for visualization and analysis of molecular structures. It is used in TADbit to visualize chromatin strands.

Software is available at: http://www.cgl.ucsf.edu/chimera/

*This software is only needed for the visualization of thee dimensional models from inside Tadbit.*


Tadbit
======

Once done Tadbit can be downloaded from here:

https://github.com/3DGenomes/tadbit.git

Download the latest version as tar.gz and install it as follows:

::

  wget https://github.com/3DGenomes/tadbit/archive/master.zip -O tadbit.zip
  unzip tadbit.zip
  cd tadbit-master/src
  sudo PYTHONPATH=$PYTHONPATH python setup.py install
  

Finally it is a good thing to test if every thing is working fine.
Go to the test directory and run:

::

  cd ../test
  python test_all.py

