Download and Install
********************

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

