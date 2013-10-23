TADBit
******


Install python wrapper:
=======================

* Download tadbit
* Go to your tadbit/src directory
::

   cd src
   python setup.py build
   sudo python setup.py install

If something goes wrong try:
::

   python setup.py clean --all
   python setup.py build
   sudo python setup.py install


Test python wrapper:
--------------------

To test the installation:
::

   cd test
   python test_all.py


Links
=====
`Travis CI <https://travis-ci.org/#!/tkf/emacs-jedi>`_ |build-status|

.. |build-status|
   image:: https://secure.travis-ci.org/fransua/tadbit.png
           ?branch=master
   :target: http://travis-ci.org/fransua/tadbit
   :alt: Build Status

