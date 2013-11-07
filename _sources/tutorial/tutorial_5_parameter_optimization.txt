
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
    exp = my_chrom.experiments['gm06690']
    
    print my_chrom.experiments


.. parsed-literal::

    [Experiment k562 (resolution: 100Kb, TADs: 43, Hi-C rows: 639, normalized: None), Experiment gm06690 (resolution: 100Kb, TADs: 31, Hi-C rows: 639, normalized: None), Experiment k562+gm06690 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None), Experiment batch_gm06690_k562 (resolution: 100Kb, TADs: 37, Hi-C rows: 639, normalized: None)]


Optimization of IMP 3D modeling parameters
------------------------------------------


In the previous tutorial we found a specific TAD (region 406 to 448) that seemed quite conserved accross different cell types.

Next, we will optimize the three IMP parameters for this TAD. The IMP parameters to optimize are maximal distance between two non-interacting particles (maxdist), Upper-bound Z-score (upfreq) and Lower-bound Z-score (lowfreq). For details see Bau & Marti-Renom. METHODS, 2012, 58:300-306.

.. code:: python

    optimizer = IMPoptimizer(exp, 100, 200, n_models=50, n_keep=25, cutoff=1000)

.. parsed-literal::

    Experiment gm06690 (resolution: 100Kb, TADs: 31, Hi-C rows: 639, normalized: None)
    100 200


.. parsed-literal::

    /usr/local/lib/python2.7/dist-packages/pytadbit/experiment.py:521: UserWarning: WARNING: normalizing according to visibility method
      warn('WARNING: normalizing according to visibility method')


.. code:: python

    # Optimize parameters. Be aware that this step is CPU intensive. If you want to se the progress, set verbose=True.
    optimizer.run_grid_search(n_cpus=8, lowfreq_range=(-1, 0, 0.2), upfreq_range=(0.2, 0.8, 0.2), 
                              scale_range=[0.005], maxdist_range=(300, 700, 200), verbose=True)


.. parsed-literal::

        1   0.2 -1 300 0.005 0.742589665946
        2   0.2 -0.8 300 0.005 0.743731102622
        3   0.2 -0.6 300 0.005 0.743129144311
        4   0.2 -0.4 300 0.005 0.745745969261
        5   0.2 -0.2 300 0.005 0.747609817096
        6   0.2 0 300 0.005 0.745233469167
        7   0.4 -1 300 0.005 0.770054963037
        8   0.4 -0.8 300 0.005 0.76942242474
        9   0.4 -0.6 300 0.005 0.770251379695
       10   0.4 -0.4 300 0.005 0.771380386735
       11   0.4 -0.2 300 0.005 0.770437793549
       12   0.4 0 300 0.005 0.770763164744
       13   0.6 -1 300 0.005 0.76005248793
       14   0.6 -0.8 300 0.005 0.760415040959
       15   0.6 -0.6 300 0.005 0.761461637946
       16   0.6 -0.4 300 0.005 0.762443522726
       17   0.6 -0.2 300 0.005 0.762557747165
       18   0.6 0 300 0.005 0.764553556895
       19   0.8 -1 300 0.005 0.724943456545
       20   0.8 -0.8 300 0.005 0.729637474149
       21   0.8 -0.6 300 0.005 0.729382490436
       22   0.8 -0.4 300 0.005 0.727957455579
       23   0.8 -0.2 300 0.005 0.727951762476
       24   0.8 0 300 0.005 0.722713060451
       25   0.2 -1 500 0.005 0.745980082303
       26   0.2 -0.8 500 0.005 0.747126241596
       27   0.2 -0.6 500 0.005 0.746372473251
       28   0.2 -0.4 500 0.005 0.745176160556
       29   0.2 -0.2 500 0.005 0.743964799592
       30   0.2 0 500 0.005 0.743067944566
       31   0.4 -1 500 0.005 0.765954987658
       32   0.4 -0.8 500 0.005 0.765244183252
       33   0.4 -0.6 500 0.005 0.764158662204
       34   0.4 -0.4 500 0.005 0.766008963176
       35   0.4 -0.2 500 0.005 0.766867625799
       36   0.4 0 500 0.005 0.763996590231
       37   0.6 -1 500 0.005 0.758411236105
       38   0.6 -0.8 500 0.005 0.757062187495
       39   0.6 -0.6 500 0.005 0.755222961808
       40   0.6 -0.4 500 0.005 0.756663183488
       41   0.6 -0.2 500 0.005 0.754922984664
       42   0.6 0 500 0.005 0.757588493623
       43   0.8 -1 500 0.005 0.727671000083
       44   0.8 -0.8 500 0.005 0.729532867757
       45   0.8 -0.6 500 0.005 0.727060609187
       46   0.8 -0.4 500 0.005 0.724870624667
       47   0.8 -0.2 500 0.005 0.725864379595
       48   0.8 0 500 0.005 0.725582367188
       49   0.2 -1 700 0.005 0.734145482851
       50   0.2 -0.8 700 0.005 0.740273068246
       51   0.2 -0.6 700 0.005 0.738971131719
       52   0.2 -0.4 700 0.005 0.737529831999
       53   0.2 -0.2 700 0.005 0.740222403839
       54   0.2 0 700 0.005 0.734277725033
       55   0.4 -1 700 0.005 0.760666138241
       56   0.4 -0.8 700 0.005 0.767745409564
       57   0.4 -0.6 700 0.005 0.76672114195
       58   0.4 -0.4 700 0.005 0.766339698097
       59   0.4 -0.2 700 0.005 0.766485574888
       60   0.4 0 700 0.005 0.765149451534
       61   0.6 -1 700 0.005 0.748673149088
       62   0.6 -0.8 700 0.005 0.751237098688
       63   0.6 -0.6 700 0.005 0.752105091377
       64   0.6 -0.4 700 0.005 0.755911312093
       65   0.6 -0.2 700 0.005 0.75376621202
       66   0.6 0 700 0.005 0.752738662488
       67   0.8 -1 700 0.005 0.719467207186
       68   0.8 -0.8 700 0.005 0.718419703596
       69   0.8 -0.6 700 0.005 0.719658609213
       70   0.8 -0.4 700 0.005 0.72112707188
       71   0.8 -0.2 700 0.005 0.717889000058
       72   0.8 0 700 0.005 0.722267543013


.. note::
   The above warning is given when a small matrix is loaded. TADBit has a filtering function that is applied to all Hi-C matrices with the aim of removing entire rows with very low counts. Those rows/colums are treated then for modeling as "missing-data" points. This flitering function can only be applied for relatively large matrices.


Optimizing from Experiment
~~~~~~~~~~~~~~~~~~~~~~~~~~


The exact same as above can be done from Experiment objects directly:

.. code:: python

    optimizer = exp.optimal_imp_parameters(100, 200, n_cpus=8, n_models=50, n_keep=25, cutoff=1000,
                                           lowfreq_range=(-1, 0, 0.2), upfreq_range=(0.2, 0.8, 0.2), 
                                           scale_range=[0.005], maxdist_range=(300, 700, 200), verbose=False)


.. parsed-literal::

    Experiment gm06690 (resolution: 100Kb, TADs: 31, Hi-C rows: 639, normalized: visibility)
    100 200


Visualize the results
---------------------


.. code:: python

    optimizer.write_result('results.log')
.. code:: python

    # Visualize the results of the optimization.
    optimizer.plot_2d()



.. image:: ../nbpictures/tutorial_5_parameter_optimization_20_0.png


We can also ask to mark on the plot the best N combination of parameters with the "show_best" parameter.

.. code:: python

    # Visualize the results of the optimization and mark the best 10 parameter sets
    optimizer.plot_2d(show_best=20)



.. image:: ../nbpictures/tutorial_5_parameter_optimization_22_0.png


.. code:: python

    axes_range = [[float(i) for i in optimizer.scale_range],
                        [float(i) for i in optimizer.maxdist_range],
                        [float(i) for i in optimizer.upfreq_range],
                        [float(i) for i in optimizer.lowfreq_range]]
    
    print axes_range
    [round(i, 3) for i in axes_range[3]]
    result = optimizer._result_to_array()
    wax = [round(i, 3) for i in axes_range[0]]
    zax = [round(i, 3) for i in axes_range[1]]
    xax = [round(i, 3) for i in axes_range[3]]
    yax = [round(i, 3) for i in axes_range[2]]
    sort_result = sorted([(result[i, j, k, l], wax[i], zax[j], xax[l], yax[k])
                                  for i in range(len(wax))
                                  for j in range(len(zax))
                                  for k in range(len(yax))
                                  for l in range(len(xax))
                                  if not np.isnan(result[i, j, k, l])
                                  ], key=lambda x: x[0],
                                 reverse=True)[0]
    print sort_result

.. parsed-literal::

    [[0.005], [300.0, 500.0, 700.0], [0.2, 0.4, 0.6, 0.8], [-1.0, -0.8, -0.6, -0.4, -0.2, 0.0]]
    (0.77138038673540765, 0.005, 300.0, -0.4, 0.4)


One can also visualize the parameter optimization according to ne of the three optimization parameters.

.. code:: python

    # Visualize the results of the optimization based on the lowfreq parameter.
    optimizer.plot_2d(axes=('upfreq', 'lowfreq', 'maxdist', 'scale'),show_best=10)



.. image:: ../nbpictures/tutorial_5_parameter_optimization_25_0.png


.. code:: python

    optimizer.plot_2d(skip={"scale":0.005}, show_best=10)



.. image:: ../nbpictures/tutorial_5_parameter_optimization_26_0.png


TADBit also provides the possibility to view it all together in a 3D plot (note that, while here its a static image, inside matplotlib GUI you would be able to turn around and zoom):

.. code:: python

    # Visualize the results of the optimization using a 3D representation with the three optimization parameters in the axis.
    optimizer.plot_3d(axes=('scale', 'maxdist', 'upfreq', 'lowfreq'))



.. image:: ../nbpictures/tutorial_5_parameter_optimization_28_0.png


.. code:: python

    optimizer.run_grid_search(n_cpus=8, lowfreq_range=(-1., -0.0, 0.1), upfreq_range=(0.3, 0.6, 0.05), 
                              scale_range=[0.005], maxdist_range=[200,250,300,350], verbose=False)

.. code:: python

    optimizer.scale_range
    optimizer.maxdist_range




.. parsed-literal::

    ['200', '250', '300', '350', '500', '700']



.. code:: python

    optimizer.plot_2d(show_best=100)


.. image:: ../nbpictures/tutorial_5_parameter_optimization_31_0.png


.. code:: python

    optimizer.write_result('results.log')
.. code:: python

    optimizer2 = IMPoptimizer(exp, 100, 200, n_models=50, n_keep=25, cutoff=1000)

.. parsed-literal::

    Experiment gm06690 (resolution: 100Kb, TADs: 31, Hi-C rows: 639, normalized: visibility)
    100 200


.. code:: python

    optimizer2.load_from_file('results.log')
.. code:: python

    optimizer2.results.keys()[105]



.. parsed-literal::

    ('0.005', '300', '0.2', '-0.4')



.. code:: python

    optimizer2.plot_2d(show_best=20)


.. image:: ../nbpictures/tutorial_5_parameter_optimization_36_0.png


Retrieve best parameters
------------------------


Once done, best results can be returned as a dictionary to be used for modeling (see next section of the tutorial)

.. code:: python

    config = optimizer.get_best_parameters_dict(reference='gm cell from Job Dekker 2009')
    
    print config


.. parsed-literal::

    {'maxdist': 250.0, 'upfreq': 0.4, 'kforce': 5, 'reference': 'gm cell from Job Dekker 2009', 'lowfreq': -0.1, 'scale': 0.005}

