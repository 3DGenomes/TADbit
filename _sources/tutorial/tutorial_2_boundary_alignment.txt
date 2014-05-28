
Alignment of TAD boundaries
===========================


.. contents::
   :depth: 3


TADbit allows to use the information from different Hi-C experiments and to put it together in order to 
decide whether some TAD boundaries are conserved or not.

Following with the example in the previous section (:ref:`getting_start`), we will load one extra experiment 
(from the same works of [Lieberman-Aiden2009]_):

.. code:: python

    from pytadbit import Chromosome
    
    # initiate a chromosome object that will store all Hi-C data and analysis
    my_chrom = Chromosome(name='My fisrt chromosome', centromere_search=True)
    
    # load Hi-C data
    my_chrom.add_experiment('First Hi-C experiment', hic_data="../../scripts/sample_data/HIC_k562_chr19_chr19_100000_obs.txt", resolution=100000)
    my_chrom.add_experiment('Second Hi-C experiment', hic_data="../../scripts/sample_data/HIC_gm06690_chr19_chr19_100000_obs.txt", resolution=100000)
    
    # run core tadbit function to find TADs, on each experiment
    my_chrom.find_tad('First Hi-C experiment')
    my_chrom.find_tad('Second Hi-C experiment')
       
    print my_chrom.experiments


.. ansi-block::

    [Experiment First Hi-C experiment (resolution: 100Kb, TADs: 37, Hi-C rows: 639, normalized: visibility), Experiment Second Hi-C experiment (resolution: 100Kb, TADs: 34, Hi-C rows: 639, normalized: visibility)]


.. ansi-block::

    /usr/local/lib/python2.7/dist-packages/pytadbit/parsers/hic_parser.py:93: UserWarning: WARNING: non integer values
      warn('WARNING: non integer values')


We now have loaded two Hi-C experiments, both at 100 Kb resolution, and have predicted the location of TADs in each of them (42 TADs detected in the first experiment and 31 in the second).

Aligning boundaries
-------------------


To align TAD boundaries several algorithms have been implemented 
(see :func:`pytadbit.chromosome.Chromosome.align_experiments`); our recommendation, however, is to use 
the default "reciprocal" method (:func:`pytadbit.boundary_aligner.reciprocally.reciprocal`). 

*Note: If the align_experiments function is run with no argument, by default all the loaded experiments will be aligned.*

Continuing with the example, the two loaded experiments are aligned as follow:


.. code:: python

    my_chrom.align_experiments(names=["First Hi-C experiment", "Second Hi-C experiment"])
    
    print my_chrom.alignment

.. ansi-block::

    {('First Hi-C experiment', 'Second Hi-C experiment'): Alignment of boundaries (length: 56, number of experiments: 2)}


All the alignments done between the experiments belonging to the same chromosome are stored under the 
alignment dictionary attached to the Chromosome object. Each alignment is an object itself 
(see :class:`pytadbit.alignment.Alignment`)


Check alignment consistency through randomization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


In order to check that the alignment makes sense and that it does not correspond to a random association of boundaries, the “randomize” parameter can be set to True when aligning:

.. code:: python

    score, pval = my_chrom.align_experiments(randomize=True, rnd_method="interpolate",
                                             rnd_num=1000)
    
    print 'score:', score
    print 'p-value:', pval

.. ansi-block::

    score: 0.223214285714
    p-value: 0.0


Alignment objects
-----------------


Visualization
~~~~~~~~~~~~~


The first function to call to check the quality of the generated alignments is the 
:func:`pytadbit.alignment.Alignment.write_alignment`:

.. code:: python

    ali = my_chrom.alignment[('First Hi-C experiment', 'Second Hi-C experiment')]
    
    print ali

.. ansi-block::

    Alignment shown in 100 Kb (2 experiments) (scores: [34m0[0m [34m1[0m [34m2[0m [36m3[0m [0m4[0m [1m5[0m [33m6[0m [33m7[0m [35m8[0m [35m9[0m [31m10[0m)
     First Hi-C experiment:|     [34m6[0m|     [34m6[0m| ---- |    [34m12[0m| ---- | ---- |    [33m32[0m| ---- |    [1m46[0m|    [1m57[0m|    [33m70[0m| ---- |    [36m83[0m| ---- | ---- |   [0m104[0m|   [33m109[0m| ---- | ---- |   [33m129[0m| ---- | ---- | ---- |   [33m184[0m|   [0m195[0m|   [0m237[0m|   [1m245[0m|   [0m330[0m|   [0m348[0m|   [34m353[0m| ---- |   [35m378[0m|   [36m384[0m| ---- |   [33m400[0m| ---- |   [35m413[0m| ---- |   [36m433[0m|   [33m446[0m|   [0m472[0m|   [1m478[0m|   [33m486[0m| ---- |   [36m501[0m|   [33m506[0m| ---- |   [34m522[0m|   [1m531[0m|   [0m554[0m|   [33m563[0m|   [1m570[0m|   [33m594[0m|   [33m609[0m| ---- |   [31m639[0m
    Second Hi-C experiment:|     [1m5[0m| ---- |    [1m12[0m| ---- |    [35m18[0m|    [36m28[0m| ---- |    [33m42[0m|    [36m47[0m|    [33m57[0m| ---- |    [31m79[0m| ---- |    [1m86[0m|    [31m98[0m| ---- | ---- |   [35m115[0m|   [0m126[0m| ---- |   [1m131[0m|   [33m145[0m|   [0m164[0m| ---- |   [35m195[0m| ---- |   [33m246[0m| ---- | ---- | ---- |   [33m375[0m| ---- | ---- |   [33m397[0m| ---- |   [33m402[0m|   [1m413[0m|   [33m431[0m| ---- | ---- | ---- |   [35m478[0m|   [35m486[0m|   [33m498[0m| ---- | ---- |   [33m510[0m| ---- |   [0m531[0m|   [31m554[0m|   [1m563[0m|   [33m569[0m|   [33m593[0m|   [31m609[0m|   [0m624[0m|   [31m639[0m
    


The different colors, corresponding to the TADbit confidence in detecting the boundaries, show how conserved the boundaries are between (in this case) cell types.

Alignment can also be viewed using matplotlib (already mention in :ref:`density_plot`):

.. code:: python

    ali.draw()


.. image:: ../nbpictures/tutorial_2_boundary_alignment_17_0.png


*Note that this function can also be zoomed in.*

The get\_column function
~~~~~~~~~~~~~~~~~~~~~~~~


The :func:`pytadbit.alignment.Alignment.get_column` function allows to select specific columns of an alignment. 

To select, for example, the third column of an alignment:


.. code:: python

    ali.get_column(3)



.. ansi-block::

    [(2, [>-<, >1100<])]



The first element of the tuple is the column index, while the two values of the second element of the tuple 
are the TADs associated to the aligned boundaries in that column. Note that TAD objects are represented 
between the '>' and '<' symbols (see: :class:`pytadbit.alignment.TAD`).

The :func:`pytadbit.alignment.Alignment.get_column` function can also take as an argument a function, in 
order to select a column (or several) depending on a specific condition. For example, to select all the 
boundaries with a score higher than 7:



.. code:: python

    cond1 = lambda x: x['score'] > 7

and to the get the selected columns:

.. code:: python

    ali.get_column(cond1=cond1)




.. ansi-block::

    [(55, [>63800<, >63800<])]



resulting, in the selection of these 3 columns.

To add a second condition, e.g. to select only the columns after the 50th column of the alignment:

.. code:: python

    cond2 = lambda x: x['pos'] > 50
    ali.get_column(cond1=cond1, cond2=cond2)



.. ansi-block::

    [(55, [>63800<, >63800<])]



Finally, to be more flexible, this conditions can be applied to only a given number of experiments (in this example of a pairwise alignment, it does not make a lot of sense):

.. code:: python

    ali.get_column(cond1=cond1, cond2=cond2, min_num=1)



.. ansi-block::

    [(53, [>60800<, >60800<]), (55, [>63800<, >63800<])]


