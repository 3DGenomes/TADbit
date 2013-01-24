Download and Install Tadbit for R
*********************************

.. contents::

*Note: for now this section is only relative to installation under Debian-linux*

GNU/Linux
=========

Tadbit requires R>=??? as well as several dependencies:

* **blablaabl**


Installing Tadbit
==================

Once done EcoloPy can be downloaded from here:

https://github.com/gui11aume/tadbit.git

Download the package as tar.gz and uncompress it:

::

  tar xzvf tadbit-tadbit-master.tar.gz

once done, go in the tadbit/src directory and install it:

::

  sudo R CMD INSTALL /full_path_to_tadbit

finally it is a good thing to test if every thing is working fine.

Go to the test directory and run:

::

  R --vanilla < test.R 

