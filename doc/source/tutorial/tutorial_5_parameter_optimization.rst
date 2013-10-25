
Parameter optimization for IMP.
===============================


*Recover data from previous section by loading the previously saved chromosome and its Hi-C data.*

Import the necessary libraries.

.. code:: python

    # Libraries
    from pytadbit import load_chromosome # to load chromosomes
    from pytadbit.imp.impoptimizer import IMPoptimizer

First, load the chromosome from previous tutorial (:ref:`run_tadbit`).

.. code:: python

    # Load the chromosome
    my_chrom = load_chromosome('some_path.tdb')

Next, load Hi-C data for each experiment (Hi-C data is not saved inside chromosome objects because of their size):

.. code:: python

    # Loop over experiments in chromosome and load Hi-C data.
    res = 100000
    
    for exp in my_chrom.experiments:
        try:
            exp.load_hic_data('../../scripts/sample_data/HIC_{0}_{1}_{1}_{2}_obs.txt'.format(
                              exp.name, my_chrom.name, res))
        except AttributeError:
            print 'file not found for experiment: ' + exp.name
            continue
        print exp


.. parsed-literal::

    Experiment k562 (resolution: 100Kb, TADs: 43, Hi-C rows: 639, normalized: None)
    Experiment gm06690 (resolution: 100Kb, TADs: 31, Hi-C rows: 639, normalized: None)
    file not found for experiment: k562+gm06690
    file not found for experiment: batch_gm06690_k562


The log indicates that experiment "k562+gm06690" had no file. Such experiment was built ad-hoc in our previous tutorial and needs to be created again by summing the Hi-C matrices from the individual experiments.

.. code:: python

    # Load Hi-C of the individual experiments and put it into the sum experiment BR+TR1+TR2
    my_chrom.experiments['k562+gm06690'].load_hic_data(
                  (my_chrom.experiments['k562'] + my_chrom.experiments['gm06690']).hic_data, 
                  'k562+gm06690')
    exp = my_chrom.experiments['k562+gm06690']
    
    print my_chrom.experiments


.. parsed-literal::

    [Experiment k562 (resolution: 100Kb, TADs: 43, Hi-C rows: 639, normalized: None), Experiment gm06690 (resolution: 100Kb, TADs: 31, Hi-C rows: 639, normalized: None), Experiment k562+gm06690 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None), Experiment batch_gm06690_k562 (resolution: 100Kb, TADs: 37, Hi-C rows: 639, normalized: None)]


Optimization of IMP 3D modeling parameters
------------------------------------------


In the previous tutorial we found a specific TAD (region 406 to 448) that seemed quite conserved accross different cell types.

Next, we will optimize the three IMP parameters for this TAD. The IMP parameters to optimize are maximal distance between two non-interacting particles (maxdist), Upper-bound Z-score (upfreq) and Lower-bound Z-score (lowfreq). For details see Bau & Marti-Renom. METHODS, 2012, 58:300-306.

.. code:: python

    optimizer = IMPoptimizer(exp, 405, 449, n_models=100, n_keep=50)

.. parsed-literal::

    Experiment k562+gm06690 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None)
    405 449


.. parsed-literal::

    /usr/local/lib/python2.7/dist-packages/TADBit-0.1-py2.7-linux-x86_64.egg/pytadbit/experiment.py:529: UserWarning: WARNING: normalizing according to visibility method
      warn('WARNING: normalizing according to visibility method')


.. code:: python

    # Optimize parameters. Be aware that this step is CPU intensive. If you want to se the progress, set verbose=True.
    optimizer.run_grid_search(n_cpus=8, lowfreq_range=(-1, -0.2, 0.1), upfreq_range=(0.0, 1.0, 0.1), 
                              scale_range=(0.005, 0.009, 0.001), maxdist_range=(350, 500, 50), verbose=False)


.. parsed-literal::

    /usr/lib/python2.7/dist-packages/numpy/lib/function_base.py:2124: RuntimeWarning: invalid value encountered in divide
      return c/sqrt(multiply.outer(d,d))


.. note::
   The above warning is given when a small matrix is loaded. TADBit has a filtering function that is applied to all Hi-C matrices with the aim of removing entire rows with very low counts. Those rows/colums are treated then for modeling as "missing-data" points. This flitering function can only be applied for relatively large matrices.


Visualize the results
---------------------


.. code:: python

    optimizer.write_result('results.log')
.. code:: python

    # Visualize the results of the optimization.
    optimizer.plot_2d()

We can also ask to mark on the plot the best N combination of parameters with the "show_best" parameter.

.. code:: python

    # Visualize the results of the optimization and mark the best 10 parameter sets
    optimizer.plot_2d(show_best=20)

One can also visualize the parameter optimization according to ne of the three optimization parameters.

.. code:: python

    # Visualize the results of the optimization based on the lowfreq parameter.
    optimizer.plot_2d(axes=('upfreq', 'lowfreq', 'maxdist', 'scale'),show_best=10)

.. code:: python

    optimizer.plot_2d(skip={"scale":0.008}, show_best=10)

TADBit also provides the possibility to view it all together in a 3D plot (note that, while here its a static image, inside matplotlib GUI you would be able to turn around and zoom):

.. code:: python

    # Visualize the results of the optimization using a 3D representation with the three optimization parameters in the axis.
    optimizer.plot_3d(axes=('maxdist', 'upfreq', 'lowfreq', 'scale'))

.. code:: python

    optimizer.run_grid_search(n_cpus=8, lowfreq_range=(-1, -0.4, 0.1), upfreq_range=(0.0, 0.1, 0.05), 
                              scale_range=(0.006, 0.008, 0.001), maxdist_range=(375, 475, 25), verbose=False)

.. code:: python

    optimizer.plot_2d(show_best=100)
.. code:: python

    optimizer.write_result('results.log')
.. code:: python

    optimizer2 = IMPoptimizer(exp, 405, 449, n_models=100, n_keep=50)
.. code:: python

    optimizer2.load_from_file('results.log')
.. code:: python

    optimizer2.results.keys()[125]
.. code:: python

    optimizer2.plot_2d(show_best=20)