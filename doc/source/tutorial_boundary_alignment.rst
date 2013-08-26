=============================
 Alignment of TAD boundaries
=============================

.. contents::
   :depth: 3


Tadbit allows to use information about different Hi-C experiments and put them together in order to decide whether some TAD boundaries are conserved or not.

Following with the example in previous section (:ref:`getting_start`), we will load one extra experiment (from the same works of [Lieberman-Aiden2009]_)

::

   from pytadbit import Chromosome
  
   # initiate a chromosome object that will store all Hi-C data and analysis
   my_chrom = Chromosome(name='My fisrt chromosome')

   # load Hi-C data
   my_chrom.add_experiment('First Hi-C experiment', xp_handler="sample_data/HIC_k562_chr19_chr19_100000_obs.txt", resolution=100000)
   my_chrom.add_experiment('Second Hi-C experiment', xp_handler="sample_data/HIC_gm06690_chr19_chr19_100000_obs.txt", resolution=100000)

   # run core tadbit function to find TADs, on each experiment
   my_chrom.find_tad('First Hi-C experiment')
   my_chrom.find_tad('Second Hi-C experiment')
   
   print my_chrom.experiments


Here the output of the last print statement:

::

   [Experiment First Hi-C experiment (resolution: 100Kb, TADs: 42, Hi-C rows: 639),
   Experiment Second Hi-C experiment (resolution: 100Kb, TADs: 31, Hi-C rows: 639)]   

We thus have loaded two Hi-C experiments, both at 100 Kb resolution, and have predicted the location of TADs in each (42 in the case of the first experiment and 31 in the former). 

Aligning boundaries
===================

In order to align TAD boundaries several algorithms have been implemented (see :func:`pytadbit.chromosome.Chromosome.align_experiments`), however our recommendation is to use the default "reciprocal" method (:func:`pytadbit.boundary_aligner.reciprocally.reciprocal`). 

Following with the example, we align these two experiments:

::

   my_chrom.align_experiments(names=["First Hi-C experiment", "Second Hi-C experiment"])

   print my_chrom.alignment

In this case we may run align_experiments function with no argument, by default all loaded experiments will be aligned.

The result of the print statement is:

:: 

   {('First Hi-C experiment', 'Second Hi-C experiment'): Alignment of boundaries (length: 60, number of experiments: 2)}

Under the alignment dictionary attached to Chromosomes, are stored all the alignments done between the experiments belonging to this same chromosome. Each alignment is an object (see :class:`pytadbit.alignment.Alignment`)


Check alignment consistency through randomization
-------------------------------------------------

In order to check that the alignment makes sense and that it do not corresponds to a random association of boundaries, we can set the "randomize" parameter to True when aligning.:

::

   score, pval = my_chrom.align_experiments(randomize=True, rnd_method="interpolate", rnd_num=1000)
   
For randomization two methodologies are available, "shuffle" or "interpolate" (the default). The first one, is based on shuffling the observed TADs while the former one draws the distribution of TAD lengths and generates random set of TADs according to this distribution (see:* :func:`pytadbit.alignment.randomization_test` *for more details).*

In the case of our example, the score of the alignment is 0.27 and its p-value < 0.001.


Alignment objects
=================

Visualization
-------------

The first function we may want to call in order to see how well does our Hi-C experiment align, is the :func:`pytadbit.alignment.Alignment.write_alignment`:

::

   ali = my_chrom.alignment[('First Hi-C experiment', 'Second Hi-C experiment')]
   
.. raw:: html

   <?xml version="1.0" encoding="UTF-8" ?>
            <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
            <!-- This file was created with the aha Ansi HTML Adapter. http://ziz.delphigl.com/tool_aha.php -->
            <html xmlns="http://www.w3.org/1999/xhtml">
            <head>
            <meta http-equiv="Content-Type" content="application/xml+xhtml; charset=UTF-8" />
            <title>stdin</title>
            </head>
            <h1></h1>
            <body>
            <pre>Alignment shown in Kb (2 experiments) (scores: <span>0</span> <span style="color:blue;">1</span> <span style="color:blue;">2</span> <span style="color:purple;">3</span> <span style="color:purple;">4</span> <span style="color:teal;">5</span> <span style="color:teal;">6</span> <span style="color:olive;">7</span> <span style="color:olive;">8</span> <span style="color:red;">9</span> <span style="color:red;">10</span>)
     First Hi-C experiment :|   <span style="color:blue;">500</span>|  <span style="color:teal;">1200</span>| ---- | ---- |  <span style="color:olive;">3100</span>| ---- |  <span style="color:teal;">4500</span>| ---- |  <span style="color:purple;">5800</span>|  <span style="color:teal;">6900</span>|  <span style="color:blue;">7700</span>| ---- | ---- | <span style="color:olive;">10300</span>| <span style="color:purple;">10800</span>| <span style="color:purple;">11400</span>| <span style="color:blue;">12400</span>| ---- | <span style="color:blue;">13100</span>| <span style="color:purple;">13600</span>| <span style="color:olive;">14400</span>| <span style="color:teal;">16300</span>| <span style="color:teal;">18300</span>| <span style="color:blue;">18800</span>| <span style="color:olive;">19400</span>| <span style="color:red;">24400</span>| <span style="color:red;">32900</span>| <span style="color:purple;">34700</span>| <span style="color:teal;">35500</span>| <span style="color:olive;">37700</span>| <span style="color:purple;">38300</span>| ---- | <span style="color:purple;">39900</span>| ---- | <span style="color:red;">41200</span>| ---- | <span style="color:purple;">43400</span>| <span style="color:teal;">44600</span>| <span style="color:purple;">45200</span>| <span style="color:purple;">45700</span>| <span style="color:purple;">47100</span>| <span style="color:purple;">47700</span>| <span style="color:olive;">48500</span>| <span style="color:purple;">49700</span>| <span style="color:teal;">50500</span>| ---- | <span style="color:purple;">52300</span>| <span style="color:olive;">53000</span>| <span style="color:teal;">55300</span>| <span style="color:teal;">56200</span>| ---- | <span style="color:teal;">59300</span>| <span style="color:olive;">60800</span>| ---- | <span style="color:red;">63800</span>
     Second Hi-C experiment:|   <span style="color:purple;">400</span>|  <span style="color:teal;">1100</span>|  <span style="color:olive;">1700</span>|  <span style="color:blue;">2600</span>| ---- |  <span style="color:teal;">4100</span>|  <span style="color:blue;">4600</span>|  <span style="color:olive;">5600</span>| ---- | ---- |  <span style="color:red;">7800</span>|  <span style="color:teal;">8500</span>|  <span style="color:red;">9700</span>| ---- | ---- | <span style="color:red;">11400</span>| ---- | <span style="color:teal;">12600</span>| ---- | ---- | ---- | ---- | ---- | ---- | <span style="color:red;">19400</span>| <span style="color:red;">24500</span>| ---- | ---- | ---- | <span style="color:red;">37700</span>| ---- | <span style="color:teal;">39600</span>| ---- | <span style="color:teal;">40100</span>| <span style="color:teal;">41200</span>| <span style="color:teal;">42900</span>| ---- | ---- | ---- | ---- | ---- | <span style="color:red;">47700</span>| <span style="color:olive;">48500</span>| <span style="color:teal;">49700</span>| ---- | <span style="color:olive;">50900</span>| ---- | <span style="color:purple;">53000</span>| <span style="color:olive;">55300</span>| <span style="color:teal;">56200</span>| <span style="color:olive;">56800</span>| <span style="color:teal;">59200</span>| <span style="color:red;">60800</span>| <span style="color:purple;">62300</span>| <span style="color:red;">63800</span>
   </pre></body></html>

Here, in different colors (corresponding to tadbit confidence in the detection of the boundaries), we can see how conserved are the boundaries, in this case, between cell types.


The get_column function
-----------------------

:func:`pytadbit.alignment.Alignment.get_column` allows to select specific columns of an alignment. 

Let say that, for example, that we want the third column of an alignment:

::

   ali.get_column(3)

This will return:

::

   [(3, [>-<, >2600<])]

The first element of the tuple corresponds to the column index, and the two values of the second element are the TADs associated to the aligned boundaries in that column. Note that TAD objects are represented between '>' and '<' (see: :class:`pytadbit.alignment.TAD`).

Also :func:`pytadbit.alignment.Alignment.get_column` may take as argument specific function in order to test and select columns.

For example, if we now want all boundaries with scores higher than 7, we might right this kind of function:

::

   cond1 = lambda x: x['score'] > 7
   
and apply it to the get_column call:

::

   ali.get_column(cond1=cond1)

resulting in this 3 columns selected:

::

   [(24, [>19400<, >19400<]), (25, [>24400<, >24500<]), (54, [>63800<, >63800<])]

We may also want to add a second condition, let say that we want only columns after 50th column of the alignment:

::

   cond2 = lambda x: x['pos'] > 50
   ali.get_column(cond1=cond1, cond2=cond2)

Resulting in:

::

   [(54, [>63800<, >63800<])]

Finally, in order to be more flexible, we may want to apply this conditions to only a given number of experiments (of course, in this case, of a pairwise alignment, it does not make a lot of sense):

::

   ali.get_column(cond1=cond1, cond2=cond2, min_num=1)

Will result in:

::

   [(52, [>60800<, >60800<]), (54, [>63800<, >63800<])]


