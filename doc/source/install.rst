Download and Install Tadbit for Python
**************************************

.. contents::

*Note: for now this section is only relative to installation under Debian-linux*

GNU/Linux
=========

Tadbit requires python>=2.7 as well as several dependencies:

* **blablaabl**

---------------------------------------------------------

Only needed for drawing graph:

* **python-matplotlib**
* **python-numpy**

Install Python libraries:
-------------------------

**Required:**
::

  apt-get install python-scipy

Accessory:

::

  apt-get install python-matplotlib python-scipy

IMP
---

Tadbit needs a specific version of IMP, to download it proceed as this:

::

   git clone https://github.com/salilab/imp
   cd imp
   git reset --hard 847e65d44da7d06718bcad366b09264c818752d5


*Note: if you do not have git installed you can install it from your linux repositories (git-core)*

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

