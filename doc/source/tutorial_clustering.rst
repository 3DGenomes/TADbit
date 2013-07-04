TAD clustering
**************

.. contents::
   :depth: 3

Tadbit allows to compare the topology of TADs using directly the Hi-C matrix. This comparison is done using the same methodology as for protein structure comparison using contact map overlap (CMO) [DiLena2010]_.

The function that allows such comparison is :func:`pytadbit.tad_clustering.tad_cmo.optimal_cmo`. The comparisons are pairwise however the main idea of using it is to pull together several of these pairwise comparisons in order to find group of TADs with specific characteristics.

Compare two TADs
----------------

To compare TADs we need a Chromosome with defined TADs, we can thus follow with our example chromosome used in :ref:`getting_start` (following the example up to :ref:`run_tadbit`), basically we need to archieve these steps:

::

   from pytadbit import Chromosome
   my_chrom = Chromosome(name='My first chromosome')
   my_chrom.add_experiment('First Hi-C experiment', xp_handler="sample_data/HIC_k562_chr19_chr19_100000_obs.txt", resolution=100000)
   my_chrom.find_tads('First Hi-C experiment')

Once done we have the TADs defined for this chromosome. Lets select two TADs:

::

   tad1 = list(my_chrom.iter_tads('First Hi-C experiment'))[41]
   tad2 = list(my_chrom.iter_tads('First Hi-C experiment'))[39]

And we align them:

::

   from pytadbit.tad_clustering.tad_cmo import optimal_cmo

   align1, align2, score = optimal_cmo(tad1[1], tad2[1], max_num_v=8, long_nw=True, long_dist=True, method='frobenius')
   
optimal_cmo function returns two alignments corresponding to the sequence of gaps needed for each TAD to be aligned with the other. The score element contains three values, an alignment score that depends on the method used to align, and the p-value and rho value of a Spearman correlation between the two Hi-C matrices.

::

   from pytadbit.tad_clustering.tad_cmo import optimal_cmo



Here the output of the example script *'clustering.py'*:

.. figure::  pictures/clustering.png
   :align:   center
   :width:   900
