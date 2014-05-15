
TAD clustering
==============




.. contents::
   :depth: 3

TADbit allows to compare the topology of TADs using directly a Hi-C matrix. This comparison is done using 
the same methodology used for protein structure comparison with contact map overlap (CMO) [DiLena2010]_.

The function that allows such comparison is :func:`pytadbit.tad_clustering.tad_cmo.optimal_cmo`. The 
comparisons are pairwise; the main idea of using it is to pull together several of these pairwise comparisons 
in order to find group of TADs with specific characteristics.



Compare two TADs
----------------


To compare TADs, a Chromosome with defined TADs is needed. Thus, following with the example chromosome used 
in :ref:`getting_start` (following the example up to :ref:`run_tadbit`), these re the steps to follow:


.. code:: python

    from pytadbit import Chromosome
    my_chrom = Chromosome(name='My first chromosome')
    my_chrom.add_experiment('First Hi-C experiment', hic_data="../../scripts/sample_data/HIC_k562_chr19_chr19_100000_obs.txt", resolution=100000)
    my_chrom.find_tad('First Hi-C experiment')
    


.. ansi-block::

    /usr/local/lib/python2.7/dist-packages/pytadbit/parsers/hic_parser.py:89: UserWarning: WARNING: non integer values
      warn('WARNING: non integer values')
    
    WARNING: removing columns having less than 24.165 counts: (detected threshold)
        8    9   10   12  245  246  247  248  249  250  251  252  253  254  255
      256  257  258  259  260  261  262  263  264  265  266  267  268  269  270
      271  272  273  274  275  276  277  278  279  280  281  282  283  284  285
      286  287  288  289  290  291  292  293  294  295  296  297  298  299  300
      301  302  303  304  305  306  307  308  309  310  311  312  313  314  315
      316  317  318  319  320  321  322  323  324  639


Once done, all the TADs for this chromosome will be defined. To select two TADs:


.. code:: python

    tad1 = list(my_chrom.iter_tads('First Hi-C experiment'))[31]
    tad2 = list(my_chrom.iter_tads('First Hi-C experiment'))[35]


And to align them:



.. code:: python

    
    from pytadbit.tad_clustering.tad_cmo import optimal_cmo
    align1, align2, score = optimal_cmo(tad1[1], tad2[1], max_num_v=8, long_nw=True, long_dist=True, method='frobenius')


The optimal_cmo function returns two alignments corresponding to the sequence of gaps needed for each TAD to 
be aligned with the other. The score element contains three values, an alignment score that depends on the 
method used to align, and the p-value and rho value of the Spearman correlation between the two Hi-C matrices.


.. code:: python

    from pytadbit.tad_clustering.tad_cmo import optimal_cmo

Following is the output of the example script *'clustering.py'*:

.. figure::  ../pictures/clustering.png
   :align:   center
   :width:   900

