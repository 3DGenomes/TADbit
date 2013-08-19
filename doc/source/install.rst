Download and Install Tadbit for Python
**************************************

.. contents::

*Note: for now this section is only relative to installation under Debian-linux*

GNU/Linux
=========

Tadbit requires python>=2.7 as well as several dependencies that are listed below.

Only needed for drawing graph:

* **python-matplotlib**
* **python-numpy**
* **python-scipy**

Install Python libraries
------------------------

**Required:**
::

  apt-get install python-scipy
  apt-get install python-numpy

Accessory (but **highly** recommended):

::

  apt-get install python-matplotlib

IMP
---

Determination of the three dimensional structure of a given locus is done with the IMP package.

In particular, Tadbit works with version 2.0.1 (but may also work with other versions).

Installation instruction for IMP can be found at their website (http://salilab.org/imp/nightly/doc/html/ ). 

**Also we do not provide support for IMP installation**, we make available here the minimal steps to install IMP on Ubuntu machines (tested on 12.04 and 13.04):*

::

    sudo apt-get install cmake
    sudo apt-get install libboost1.49-all-dev
    sudo apt-get install libhdf5-dev
    sudo apt-get install swig
    sudo apt-get install libcgal-dev
    sudo apt-get install python-dev


Download tarball from http://salilab.org/imp/get.php?pkg=2.0.1/download/imp-2.0.1.tar.gz. And uncompress:


::

   tar xzvf imp-2.0.1.tar.gz
   cd imp-2.0.1

compile (*Note:* the `-j` option stands for the number of CPUs you want to assign to the task, the higher the faster).

::

   cmake . -DCMAKE_BUILD_TYPE=Release -DIMP_MAX_CHECKS=NONE -DIMP_MAX_LOG=SILENT
   make -j4 

Once all is done open the file setup_environment.sh in your imp directory and copy the first lines in your `~/.bashrc` file. These lines should look like:

::

  LD_LIBRARY_PATH="/something/imp-2.0.1/lib:/something/imp-2.0.1/lib:/something/imp-2.0.1/src/dependency/RMF/:$LD_LIBRARY_PATH"

  export LD_LIBRARY_PATH

  PYTHONPATH="/something/imp-2.0.1/lib:/something/imp-2.0.1/lib:/something/imp-2.0.1/src/dependency/RMF/:$PYTHONPATH"

  export PYTHONPATH


MCL
---

MCL program is used for clustering 3-dimensional models can be downloaded from http://micans.org/mcl/ .

*Note: in case MCL is not found by tadbit, an alternative clustering method will be used. Nevertheless we strongly recommend to use MCL.*

On Debian/Ubuntu machines you can install MCL from the repositories:

::

  sudo apt-get install mcl



Installing Tadbit
==================

Once done EcoloPy can be downloaded from here:

https://github.com/gui11aume/tadbit.git

Download the package as tar.gz and uncompress it:

::

  tar xzvf tadbit-tadbit-master.tar.gz

once done, go in the tadbit/src directory and install it:

::

  sudo python setup.py install

finally it is a good thing to test if every thing is working fine.

Go to the test directory and run:

::

  python test_all.py

