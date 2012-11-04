
PYTHON INSTALL:
---------------

 * Download tadbit
 * Go to your tadbit directory

   cd src
   python setup.py build
   sudo python setup.py install

If something goes wrong try:

   python setup.py clean --all
   python setup.py build
   sudo python setup.py install


TEST PYTHON INSTALL:
--------------------

To test the installation:

    cd test
    python test_all.py

