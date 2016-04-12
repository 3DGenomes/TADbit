
.. _getting_start:

Getting started
===============

Hi-C data format
----------------

Hi-C data are usually represented as symmetric matrices in a tab
separated file, as in the example below:

::

    chrT_001      chrT_002        chrT_003        chrT_004        chrT_005        chrT_006
    chrT_001      629     164     88      105     10      35
    chrT_002      86      612     175     110     40      29
    chrT_003      159     216     437     105     43      73
    chrT_004      100     111     146     278     70      42
    chrT_005      16      36      26      79      243     39
    chrT_006      19      50      53      42      37      224

TADbit allows to load most of this kind of matrices. A Hi-C matrix is
loaded as this *(note: this example loads human chromosome 19 from*
[Lieberman-Aiden2009]\_ *results.)*:

.. code:: python

    from pytadbit import Chromosome
      
    # initiate a chromosome object that will store all Hi-C data and analysis
    my_chrom = Chromosome(name='chr19', centromere_search=True, 
                          species='Homo sapiens', assembly='NCBI36')
    
    # load Hi-C data
    my_chrom.add_experiment('k562', cell_type='wild type', exp_type='Hi-C', identifier='k562',
                            project='TADbit tutorial',
                            hic_data="../../scripts/sample_data/HIC_k562_chr19_chr19_100000_obs.txt", resolution=100000)
    my_chrom.add_experiment('gm06690',  cell_type='cancer', exp_type='Hi-C', identifier='gm06690',
                            project='TADbit tutorial',
                            hic_data="../../scripts/sample_data/HIC_gm06690_chr19_chr19_100000_obs.txt", resolution=100000)



.. ansi-block::

    /usr/local/lib/python2.7/dist-packages/pytadbit/parsers/hic_parser.py:93: UserWarning: WARNING: non integer values
      warn('WARNING: non integer values')


.. warning::
   TADbit assumes that Hi-C data matrices starts at chromosome position 1. If your matrix do not represents the full chromosome length, fill the missing columns and rows with zeros.

.. note::
   by default TADbit will not search for centromeres, however, as, in this example, we know that we are working with chromosome 19, we can ask TADbit search for it as explained in :ref:`forbidden_regions`

.. note::
   TADbit parser is also able to read compressed data directly (gzip files).

.. note::
   Many optional descriptive arguments can be passed to the chromosome and to the experiments. This arguments may be used by TADbit when generating the output.

Unrecognized data format
~~~~~~~~~~~~~~~~~~~~~~~~

In some cases users may found data that are not readable by TADbit. In
such case, a parser function can be written and passed as an argument to
the ``add_experiment`` function:

.. code:: python

    def read_strange_file(f_handle):
        """
        reads from file handler (already openned)
        """
        nums = []
        for line in f_handle:
            if not line:
                continue
            # modify the following line to fit your parsing needings
            ##
            values = line.split()
            ##
            # feed the "num" list with the list of values you parsed
            nums.append([int(v) for v in values])
        return tuple([nums[j][i] for i in xrange(len(nums)) for j in xrange(len(nums))]), len(nums)


And call it as follow:

.. code:: python

    other_chrom = Chromosome(name='An other chromosome')
    other_chrom.add_experiment('First Hi-C experiment', hic_data='../../src/test/data/hESC_chr19-rep1.txt',
                               parser=read_strange_file, resolution=20000)


Experiment objects
------------------

Experiments, when loaded, are stored in a special kind of list attached
to chromosome objects:

.. code:: python

    my_chrom.experiments




.. ansi-block::

    [Experiment k562 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None),
     Experiment gm06690 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None)]



A specific Experiment can be accessed either by its name or by its
position in :class:``pytadbit.chromosome.ExperimentList`` :

.. code:: python

    my_chrom.experiments[0] == my_chrom.experiments["k562"]




.. ansi-block::

    True



Each Experiment is an independent object with a list of associated functions 
(see :class:`pytadbit.experiment.Experiment`).

.. _exp_operations:

Basic manipulation of Experiments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two Hi-C experiments can be summed up easily, resulting in a new Hi-c
experiment contatining the sum of the interaction counts of the summed
experiments:

.. code:: python

    exp = my_chrom.experiments["k562"] + my_chrom.experiments["gm06690"]
    print exp



.. ansi-block::

    Experiment k562+gm06690:
       resolution        : 100Kb
       TADs              : None
       Hi-C rows         : 639
       normalized        : None
       identifier        : k562+gm06690
       cell type         : wild type+cancer
       restriction enzyme: UNKNOWN
       project           : TADbit tutorial
    


.. ansi-block::

    WARNING: experiments should be normalized before being summed


*Note the last warning asking to normalize experiments before summing
them. Indeed normalizing the sum of different experiment should result
in bayesing the final result towards the experiment with more raw Hi-C
interaction counts.*

The resulting experiment (which default name is the concatenation of the
summed experiments) can be added to the experiments of a given
chromosome.

.. code:: python

    my_chrom.add_experiment(exp)
    print my_chrom.experiments



.. ansi-block::

    [Experiment k562 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None), Experiment gm06690 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None), Experiment k562+gm06690 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None)]


Hi-C matric visualization
~~~~~~~~~~~~~~~~~~~~~~~~~

To quickly view how does the interaction matrix look like, experiment objects have the :func:`pytadbit.experiment.Experiment.view` function

.. code:: python

    exp.view()


.. ansi-block::

    /usr/local/lib/python2.7/dist-packages/pytadbit/experiment.py:1090: RuntimeWarning: divide by zero encountered in log2
      img = axe.imshow(fun(matrix), origin='lower', vmin=vmin, vmax=vmax,



.. image:: ../nbpictures//tutorial_1_general_27_1.png




.. ansi-block::

    <matplotlib.image.AxesImage at 0x7f6fa312e7d0>



This plot shows the log2 interaction counts, resulting from the given
Hi-C experiment, and it can be drawn from individual experiments or from
the chromosome object.

.. code:: python

    my_chrom.visualize([('k562', 'gm06690'), 'k562+gm06690'])



.. image:: ../nbpictures//tutorial_1_general_29_0.png


Note how we pass the list of experiments to show to the :func:`pytadbit.chromosome.Chromosome.visualize` we ask to view 3 experiments, the first two being grouped in a single matrix. This view is useful to compare experiments, and do not suppose a lose of information as Hi-C matrices are symetric.

.. _run_tadbit:

Find Topologically Associating Domains
--------------------------------------

Once an experiment has been loaded, the location of Topologically Associating Domains (TADs) can be estimated as:

.. code:: python

    my_chrom.find_tad('k562', n_cpus=8)
    my_chrom.find_tad('gm06690', n_cpus=8)


:func:`pytadbit.chromosome.Chromosome.find_tad` is called from our Chromosome object but is applied to a 
specific experiment. Therefore, TADs found by TADBbit will be associated to this specific experiment. 
They can be accessed as following:

.. code:: python

    exp = my_chrom.experiments["k562"]
    exp.tads




.. ansi-block::

    {1: {'brk': 5.0, 'end': 5.0, 'score': 2.0, 'start': 0.0},
     2: {'brk': 5.0, 'end': 5.0, 'score': 2.0, 'start': 0.0},
     3: {'brk': 11.0, 'end': 11.0, 'score': 2.0, 'start': 6.0},
     4: {'brk': 31.0, 'end': 31.0, 'score': 7.0, 'start': 12.0},
     5: {'brk': 45.0, 'end': 45.0, 'score': 5.0, 'start': 32.0},
     6: {'brk': 56.0, 'end': 56.0, 'score': 5.0, 'start': 46.0},
     7: {'brk': 69.0, 'end': 69.0, 'score': 7.0, 'start': 57.0},
     8: {'brk': 82.0, 'end': 82.0, 'score': 3.0, 'start': 70.0},
     9: {'brk': 103.0, 'end': 103.0, 'score': 4.0, 'start': 83.0},
     10: {'brk': 108.0, 'end': 108.0, 'score': 7.0, 'start': 104.0},
     11: {'brk': 128.0, 'end': 128.0, 'score': 6.0, 'start': 109.0},
     12: {'brk': 183.0, 'end': 183.0, 'score': 6.0, 'start': 129.0},
     13: {'brk': 194.0, 'end': 194.0, 'score': 4.0, 'start': 184.0},
     14: {'brk': 236.0, 'end': 236.0, 'score': 4.0, 'start': 195.0},
     15: {'brk': 244.0, 'end': 244.0, 'score': 5.0, 'start': 237.0},
     16: {'brk': 329.0, 'end': 329.0, 'score': 4.0, 'start': 324.0},
     17: {'brk': 347.0, 'end': 347.0, 'score': 4.0, 'start': 330.0},
     18: {'brk': 352.0, 'end': 352.0, 'score': 2.0, 'start': 348.0},
     19: {'brk': 377.0, 'end': 377.0, 'score': 8.0, 'start': 353.0},
     20: {'brk': 383.0, 'end': 383.0, 'score': 3.0, 'start': 378.0},
     21: {'brk': 399.0, 'end': 399.0, 'score': 6.0, 'start': 384.0},
     22: {'brk': 412.0, 'end': 412.0, 'score': 9.0, 'start': 400.0},
     23: {'brk': 432.0, 'end': 432.0, 'score': 3.0, 'start': 413.0},
     24: {'brk': 445.0, 'end': 445.0, 'score': 6.0, 'start': 433.0},
     25: {'brk': 471.0, 'end': 471.0, 'score': 4.0, 'start': 446.0},
     26: {'brk': 477.0, 'end': 477.0, 'score': 5.0, 'start': 472.0},
     27: {'brk': 485.0, 'end': 485.0, 'score': 7.0, 'start': 478.0},
     28: {'brk': 500.0, 'end': 500.0, 'score': 3.0, 'start': 486.0},
     29: {'brk': 505.0, 'end': 505.0, 'score': 7.0, 'start': 501.0},
     30: {'brk': 521.0, 'end': 521.0, 'score': 2.0, 'start': 506.0},
     31: {'brk': 530.0, 'end': 530.0, 'score': 5.0, 'start': 522.0},
     32: {'brk': 553.0, 'end': 553.0, 'score': 4.0, 'start': 531.0},
     33: {'brk': 562.0, 'end': 562.0, 'score': 6.0, 'start': 554.0},
     34: {'brk': 569.0, 'end': 569.0, 'score': 5.0, 'start': 563.0},
     35: {'brk': 593.0, 'end': 593.0, 'score': 6.0, 'start': 570.0},
     36: {'brk': 608.0, 'end': 608.0, 'score': 7.0, 'start': 594.0},
     37: {'brk': 638.0, 'end': 638.0, 'score': 10.0, 'start': 609.0}}



The "tads" variable returned in this example is a dictionary of TADs,
each of each is in turn a new dictionary containing information about
the start and end positions of a TAD.

"start" and "end" values correspond respectively to the start and end
positions of the given TAD in the chromosome (note that this numbers
have to be multiplied by the resolution of the experiment,
"exp.resolution"). The "brk" key corresponds to the value of "end", all
"brk" together corresponds to all TAD's boundaries.

TADs to text
~~~~~~~~~~~~

TADs can also be seen or saved to a file using this write function:

.. code:: python

    exp.write_tad_borders()


.. ansi-block::

    #      start   end score
    2        1.0   6.0   2.0
    3        1.0   6.0   2.0
    4        7.0  12.0   2.0
    5       13.0  32.0   7.0
    6       33.0  46.0   5.0
    7       47.0  57.0   5.0
    8       58.0  70.0   7.0
    9       71.0  83.0   3.0
    10      84.0 104.0   4.0
    11     105.0 109.0   7.0
    12     110.0 129.0   6.0
    13     130.0 184.0   6.0
    14     185.0 195.0   4.0
    15     196.0 237.0   4.0
    16     238.0 245.0   5.0
    17     325.0 330.0   4.0
    18     331.0 348.0   4.0
    19     349.0 353.0   2.0
    20     354.0 378.0   8.0
    21     379.0 384.0   3.0
    22     385.0 400.0   6.0
    23     401.0 413.0   9.0
    24     414.0 433.0   3.0
    25     434.0 446.0   6.0
    26     447.0 472.0   4.0
    27     473.0 478.0   5.0
    28     479.0 486.0   7.0
    29     487.0 501.0   3.0
    30     502.0 506.0   7.0
    31     507.0 522.0   2.0
    32     523.0 531.0   5.0
    33     532.0 554.0   4.0
    34     555.0 563.0   6.0
    35     564.0 570.0   5.0
    36     571.0 594.0   6.0
    37     595.0 609.0   7.0
    38     610.0 639.0  10.0
    


Another way to view TADs is using the matrix visualization:

TADs in interaction matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    my_chrom.visualize(exp.name, paint_tads=True)



.. image:: ../nbpictures//tutorial_1_general_43_0.png


In this case we can also put them side by side, view a given region, and
use nrmalized data instead of the log2 of raw data:

.. code:: python

    my_chrom.visualize([('k562', 'gm06690')], paint_tads=True, focus=(490,620), normalized=True)


.. ansi-block::

    /usr/local/lib/python2.7/dist-packages/pytadbit/experiment.py:1090: RuntimeWarning: invalid value encountered in log2
      img = axe.imshow(fun(matrix), origin='lower', vmin=vmin, vmax=vmax,



.. image:: ../nbpictures//tutorial_1_general_45_1.png


*Note that the width of the line is proportional to the score of the TAD
border.*


.. _density_plot:


TADs in density plots
~~~~~~~~~~~~~~~~~~~~~

Finally TAD bourders can be seen using the density plot summary:

.. code:: python

    my_chrom.tad_density_plot('k562')



.. image:: ../nbpictures//tutorial_1_general_50_0.png


In this plot, each grey-filled-arc represents a TAD. The eight of the
TAD is proportional to the relative amount of interactions in this TAD.
If this relative amount of interactions is higher than 1 (dark grey
TADs) the number of interactions inside the TAD is higher than expected
accordin to its size (a black horizontal line, highlights this
theoretical threshold).

TAD borders in this plot are marked by colored triangles. The color
code, from blue to red, represents how confident is TADbit about the
identification of the border.

Finding TADs in related Hi-C experiments
----------------------------------------

TADbit also allows to search for TADs in a chromosome using the information of several Hi-C experiments. To do this experiments do not need to be summed up (like in :ref:`exp_operations`), as the find_tad function has a batch_mode:

.. code:: python

    my_chrom.find_tad(['k562', 'gm06690'], batch_mode=True, n_cpus=8)
    print my_chrom.experiments



.. ansi-block::

    [Experiment k562 (resolution: 100Kb, TADs: 37, Hi-C rows: 639, normalized: visibility), Experiment gm06690 (resolution: 100Kb, TADs: 34, Hi-C rows: 639, normalized: visibility), Experiment k562+gm06690 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None), Experiment batch_gm06690_k562 (resolution: 100Kb, TADs: 35, Hi-C rows: 639, normalized: visibility)]


.. note::
   In this case a new experiment is also created (e.g. "batch_gm06690_k562"), it contains the TADs detected when taking into account the two experiments "gm06690" and "k562".


.. _forbidden_regions:

Forbidden regions and centromeres
---------------------------------

By default TADbit does not put limitation in sizes of TADs, owever it may make sense to set a maximum TAD size of 5 Mb in humans [Dekker2013]_. In this case we can set this cutoff when first defining the chromosome (see the option in :class:`pytadbit.Chromosome`), or just use the function :class:`pytadbit.Chromosome.set_max_tad_size`:

.. code:: python

    my_chrom.set_max_tad_size(3000000)

Once TADs are detected by the core :func:`pytadbit.tadbit.tadbit` function, TADbit checks that they are not larger than a given value (5 Mb here). If a TAD is larger than this value, it will be marked with a 
**negative score**, and will be automatically excluded from the main TADbit functions.

.. code:: python

    my_chrom.visualize('k562', paint_tads=True)



.. image:: ../nbpictures//tutorial_1_general_61_0.png


Another optional inspection performed by TADbit is the presence of
centromeric regions. TADbit assumes that the larger gap found in a Hi-C
matrix corresponds to the centromere. In the case of this example we set
the argument "search\_centromere=True", but the default is False, as it
may give unexpected results in telocentric chromosomes.

The search for centromere is updated and refined each time a new
experiment is linked to a given Chromosome. Typically, TADs calculated
by the :func:``pytadbit.tadbit.tadbit`` function include centromeric
regions; if a centromere is found, TADbit will split the TAD containing
it into two TADs, one ending before the centromere and one starting
after.

Saving and restoring data
-------------------------

In order to avoid having to calculate TAD positions each time, TADbit allows to save and load Chromosome 
objects, with all the associated experiments. To save a Chromosome object:

.. code:: python

    my_chrom.save_chromosome("some_path.tdb", force=True)


And to load it:

.. code:: python

    from pytadbit import load_chromosome
    
    my_chrom = load_chromosome("some_path.tdb")
    
    print my_chrom.experiments



.. ansi-block::

    [Experiment k562 (resolution: 100Kb, TADs: 37, Hi-C rows: 639, normalized: None), Experiment gm06690 (resolution: 100Kb, TADs: 34, Hi-C rows: 639, normalized: None), Experiment k562+gm06690 (resolution: 100Kb, TADs: None, Hi-C rows: 639, normalized: None), Experiment batch_gm06690_k562 (resolution: 100Kb, TADs: 35, Hi-C rows: 639, normalized: None)]


*Note: while information about TADs can be saved, in order to save disk space, raw Hi-C data are not stored in this way but can be loaded again for each experiment:*

.. code:: python

    expr = my_chrom.experiments["k562"]
    
    expr.load_hic_data("../../scripts/sample_data/HIC_k562_chr19_chr19_100000_obs.txt")
