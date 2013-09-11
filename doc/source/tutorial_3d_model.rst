Three-dimensional modeling
**************************

.. contents::
   :depth: 3


Data filtering
==============

The first step in the structure determination method by IMP, as in any integrative approach, is to check 
for possible biases in the experimental data the could give rise to artifacts in the final 3D models. In 
TADBit, this is done by studying the distribution of the sum of interactions per raw/columns in the Hi-C 
matrix. According to this distribution, some columns may be removed if the sum of their interactions 
significantly lower than the median of median of all the interactions.
This step is performed automatically within TADBbit each time new experimental data are loaded. To ensure 
that only the otliers columns are removed, TADBit checks if the selected cutoff corresponds is a concave 
down in the region between zero and the median of the overall distribution. The rejected columns are stored 
in the variable Experiment._zeros, which represents the columns to be skipped in the consecutive steps.
The following example shows the cutoff (dashed red line, autmatically determined with the above described 
method) for the Human chromosome 19 data (data preloaded in the 'exp' variable). According to the fit showed 
in the plot, all the columns corresponding the Hi-C raw data with total count of interaction below the dashed 
red line (~46 in this case) will be discarded. Their corresponding particles will be treated as 
'missing particles' (particles with no experimental data).

::

  from pytadbit.utils import filter_with_polynomial_fit

  filter_with_polynomial_fit(exp.get_hic_matrix(), draw_hist=True)

.. figure::  pictures/example_filtering.png
   :align:   center


Data normalization
==================
To account for possible biases in the experimental data, the Z-score of the 5C data is computed. However, 
since the Z-score requires the data to be normally distributed, 5C data are log10 transformed prior the 
computation of the Z-score. The input data for IMP consist then in the resulting Z-score matrix computed 
on the experimental data.

The filtered Hi-C data (stored in :class:`pytadbit.experiment.Experiment`) need then to be normalized before entering the IMP pipeline. The normalization is carried out in two steps: (i) generation of the weight for each pair of interactions, based on the interaction count in the corresponding row/column; (ii) Z-score ( `z-score <http://en.wikipedia.org/wiki/Standard_score#Calculation_from_raw_score>`_) compitation for each interaction pair.

Calculation of weights
----------------------

Weights can be calculated according to two formulas (see 
:class:`pytadbit.experiment.Experiment.normalize_hic`); however, in the context of three-dimensional modeling, 
the "over_tot" method is recommended, as the distribution of values generated is closer to normal.

Hi-C interaction counts are thus normalized according to the following formula:

.. math::

  weight(I, J) = \frac{\sum^N_{i=0}{(matrix(i, J))} \times \sum^N_{j=0}{(matrix(I, j))}}{\sum^N_{i=0}{\sum^N_{j=0}{(matrix(i, j))}}}


where 'matrix' is the raw data (count of interactions), and N is the number of rows/columns.

The result is stored in a new matrix, called 'weight'. The values that will be used in the next step are the 
multiplication of this weights per the raw data.


Calculation of the z-score
--------------------------

Z-scores are computed according to classical formula (:math:`\frac{x-\mu}{\sigma}`) on the log10 transformed 
values of the normalized data (see above):

.. math::

  zscore(I, J) = \frac{log_{10}(weight(I, J) \times matrix(I, J)) - mean(log_{10}(weight \times matrix))}{stddev(log_{10}(weight \times matrix))}

**Important: values on the diagonal are not taken into account in this calculation.**

Dealing with zeros
^^^^^^^^^^^^^^^^^^

In an Hi-C interaction matrix, a value of zero means that the corresponding fragments of chromatin were 
not crosslinked in the experiment. Since this could be caused either by an experimental artifact or by 
the fact that the two corresponding fragments were never found close enough to be crosslinked, we assume 
that :math:`log_{10}(0) = 0`, in the calculation of the mean and stddev, and equal to -1 in the calculation
 of the z-score.



TADBit tutorial #1. From Hi-C data to TAD 3D models
===================================================



This tutorial will walk users in the necessary steps for modeling a 3D genomic domain from a Hi-C fly dataset (REF). First, we will load a "Chromosome" object from which we take one "Experiment" object ('exp'). From this Experiment object we can model a given region using IMP.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



1. Loading data
~~~~~~~~~~~~~~~


Import the necessary libraries and set some variables:

*****
    WE WILL HAVE TO TELL IN THE MANUAL HOW TO USE THIS TUTORIAL IN THEIR OWN COMPUTERS, WITH GRAPHS, ETC...
    REMEMBER TO CHANGE THE NAMES OF THE LIBRARIES THAT YOU MAY HAVE CHANGED!
*****

::

    # Libraries
    from pytadbit import Chromosome
    from pytadbit.tad_clustering.tad_cmo import optimal_cmo
    from pytadbit.imp.structuralmodels import load_structuralmodels
    from pytadbit.imp.CONFIG import CONFIG

Definition of a Chromosome object:

::

    crm = '2R' # Chromosome name
    my_chrom = Chromosome(crm) # Create the chromosome object

Load all experiments done on Drosophila's chromosome 2R (Hi-C matrices), which include Corces' technical 
and biolobical replicates. Then visualize the data: 


::

    # Load three different experimental data named TR1, TR2 and BR
    datasets = ['TR1', 'TR2', 'BR']
    for xnam in datasets:
        my_chrom.add_experiment(xnam, 
                                resolution=10000, 
                                xp_handler='/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
    my_chrom.experiments
    
    # Visualize the experimental data for each of the three experiments
    for xnam in datasets:
        my_chrom.visualize(xnam)


.. parsed-literal::

    /usr/local/lib/python2.7/dist-packages/pytadbit/chromosome.py:468: RuntimeWarning: divide by zero encountered in log2
      img = axe.imshow(fun(matrix), origin='lower', vmin=vmin, vmax=vmax,

.. image:: pictures/Tadbit_for_IMP_fransua_8_1.png

.. image:: pictures/Tadbit_for_IMP_fransua_8_2.png

.. image:: pictures/Tadbit_for_IMP_fransua_8_3.png


2. Calculating TAD borders
~~~~~~~~~~~~~~~~~~~~~~~~~~


Calculate TAD borders for the three expriments


::

    # Identify TADs of the experimental data saved in my_chrom (about 20' execution time)
    for xnam in datasets:
        my_chrom.find_tad(xnam) # Finding TADs can take significant CPU time (by default it uses all available CPUs)
        my_chrom.visualize(xnam, paint_tads=True)

.. image:: pictures/Tadbit_for_IMP_fransua_11_0.png

.. image:: pictures/Tadbit_for_IMP_fransua_11_1.png

.. image:: pictures/Tadbit_for_IMP_fransua_11_2.png

Calculating TADs takes significant CPU time. Thus, it is a good idea to save the results and reload them when
 needed:

::

    # Save the TAD definition for all experiments. To re-load a chromosome use the load_chromosome() function.
    my_chrom.save_chromosome('my_chrom.tdb')

Align the three experiments (TAD borders) and visualize the alignment:

::

    # Align all experiments
    my_chrom.align_experiments()
    my_chrom.alignment


.. parsed-literal::

    {('BR',
      'TR1',
      'TR2'): Alignment of boundaries (length: 210, number of experiments: 3)}



::

    # Visualize the result of the alignment
    ali = my_chrom.alignment[('BR','TR1','TR2')]
    ali.write_alignment()


.. raw:: html

   <pre>Alignment shown in 10 Kb (3 experiments) (scores: <span class="ansiblue">0</span> <span class="ansiblue">1</span> <span class="ansiblue">2</span> <span class="ansicyan">3</span> 4 5</span> <span class="ansiyellow">6</span> <span class="ansiyellow">7</span> <span class="ansipurple">8</span> <span class="ansipurple">9</span> <span class="ansired">10</span>)
   TR1:| ---- | ---- | ---- | ---- | ---- |   <span class="ansired">155</span>|   <span class="ansipurple">160</span>|   <span class="ansired">165</span>|   <span class="ansired">184</span>|   <span class="ansired">191</span>|   <span class="ansired">196</span>|   <span class="ansired">203</span>|   210|   <span class="ansired">215</span>|   <span class="ansired">226</span>|   <span class="ansired">237</span>|   <span class="ansired">251</span>|   <span class="ansipurple">263</span>|   <span class="ansipurple">269</span>|   <span class="ansiyellow">283</span>|   <span class="ansiyellow">289</span>| ---- | ---- |   <span class="ansired">334</span>|   <span class="ansired">341</span>|   <span class="ansiyellow">355</span>|   <span class="ansipurple">361</span>| ---- |   <span class="ansired">371</span>|   <span class="ansipurple">379</span>|   <span class="ansired">387</span>|   <span class="ansired">395</span>|   <span class="ansipurple">405</span>|   <span class="ansired">448</span>|   <span class="ansipurple">454</span>|   <span class="ansiyellow">461</span>| ---- |   <span class="ansired">477</span>|   <span class="ansired">483</span>|   <span class="ansiyellow">497</span>|   <span class="ansiyellow">502</span>|   <span class="ansired">508</span>|   <span class="ansired">517</span>|   <span class="ansiyellow">527</span>|   <span class="ansiyellow">532</span>|   <span class="ansired">542</span>|   <span class="ansired">547</span>|   <span class="ansired">553</span>|   <span class="ansiyellow">558</span>|   <span class="ansired">569</span>|   <span class="ansiyellow">587</span>|   <span class="ansiyellow">596</span>|   <span class="ansiyellow">601</span>|   <span class="ansired">629</span>|   <span class="ansired">644</span>|   <span class="ansired">652</span>|   <span class="ansired">669</span>|   <span class="ansipurple">678</span>|   700|   705|   <span class="ansiyellow">711</span>|   <span class="ansiyellow">716</span>|   <span class="ansiyellow">721</span>|   <span class="ansiyellow">728</span>|   <span class="ansiyellow">733</span>|   <span class="ansired">746</span>|   <span class="ansiyellow">752</span>|   <span class="ansiyellow">757</span>|   <span class="ansired">772</span>|   <span class="ansired">777</span>| ---- |   789|   <span class="ansired">794</span>|   <span class="ansired">802</span>|   <span class="ansired">808</span>|   <span class="ansired">817</span>|   <span class="ansired">822</span>|   <span class="ansired">830</span>|   <span class="ansipurple">835</span>|   <span class="ansipurple">849</span>|   <span class="ansired">856</span>|   <span class="ansired">864</span>|   <span class="ansired">872</span>|   878</span>|   <span class="ansired">883</span>|   <span class="ansired">892</span>|   <span class="ansipurple">901</span>|   <span class="ansipurple">906</span>|   <span class="ansipurple">912</span>|   <span class="ansiyellow">931</span>|   <span class="ansipurple">936</span>|   <span class="ansired">945</span>|   <span class="ansired">974</span>|   <span class="ansired">981</span>|   <span class="ansired">987</span>|  <span class="ansired">1003</span>|  1009</span>|  <span class="ansipurple">1014</span>|  <span class="ansired">1036</span>|  <span class="ansired">1050</span>|  <span class="ansired">1063</span>|  <span class="ansipurple">1069</span>|  <span class="ansipurple">1075</span>|  <span class="ansipurple">1084</span>|  1090|  1096</span>|  <span class="ansired">1102</span>|  <span class="ansired">1108</span>|  <span class="ansiyellow">1113</span>|  <span class="ansipurple">1120</span>|  <span class="ansired">1126</span>|  1139</span>|  <span class="ansiyellow">1145</span>|  1173</span>|  <span class="ansiyellow">1181</span>|  <span class="ansiyellow">1189</span>|  <span class="ansipurple">1197</span>|  <span class="ansipurple">1202</span>| ---- |  <span class="ansired">1210</span>|  <span class="ansired">1222</span>|  <span class="ansired">1246</span>|  <span class="ansired">1265</span>|  <span class="ansipurple">1271</span>|  <span class="ansiyellow">1277</span>|  <span class="ansired">1289</span>|  <span class="ansired">1298</span>|  <span class="ansiyellow">1326</span>|  <span class="ansiyellow">1331</span>|  1336|  <span class="ansiyellow">1342</span>|  <span class="ansipurple">1349</span>|  <span class="ansired">1354</span>|  <span class="ansired">1361</span>|  <span class="ansipurple">1367</span>|  <span class="ansired">1375</span>|  <span class="ansired">1401</span>|  <span class="ansired">1406</span>|  <span class="ansipurple">1427</span>|  <span class="ansipurple">1432</span>|  <span class="ansired">1449</span>|  <span class="ansiyellow">1454</span>|  <span class="ansiyellow">1459</span>|  <span class="ansired">1467</span>|  <span class="ansired">1474</span>|  <span class="ansired">1500</span>|  1510</span>|  <span class="ansipurple">1516</span>|  <span class="ansipurple">1525</span>|  <span class="ansipurple">1533</span>|  <span class="ansipurple">1538</span>|  1558|  1565|  <span class="ansired">1612</span>|  <span class="ansired">1621</span>|  <span class="ansired">1646</span>|  <span class="ansipurple">1652</span>|  <span class="ansiyellow">1680</span>|  <span class="ansiyellow">1685</span>|  <span class="ansired">1690</span>|  <span class="ansired">1696</span>|  <span class="ansipurple">1706</span>| ---- |  <span class="ansipurple">1717</span>|  <span class="ansired">1726</span>|  <span class="ansired">1736</span>|  <span class="ansired">1741</span>|  <span class="ansired">1748</span>|  <span class="ansiyellow">1754</span>|  <span class="ansired">1760</span>|  <span class="ansired">1795</span>|  <span class="ansiyellow">1800</span>|  <span class="ansired">1809</span>|  <span class="ansipurple">1838</span>|  <span class="ansired">1844</span>|  <span class="ansired">1849</span>|  <span class="ansired">1858</span>|  <span class="ansired">1868</span>|  <span class="ansiyellow">1875</span>|  <span class="ansipurple">1880</span>|  1885</span>|  1890</span>|  <span class="ansipurple">1895</span>|  <span class="ansired">1924</span>|  <span class="ansired">1929</span>|  <span class="ansired">1940</span>|  1946|  1953</span>|  <span class="ansiyellow">1958</span>|  <span class="ansiyellow">1963</span>|  <span class="ansired">1971</span>|  <span class="ansired">1978</span>|  <span class="ansired">1985</span>|  <span class="ansired">1994</span>|  <span class="ansired">1999</span>|  <span class="ansired">2007</span>|  <span class="ansired">2018</span>|  <span class="ansiyellow">2024</span>|  <span class="ansired">2040</span>|  2045</span>|  2052|  2057</span>|  <span class="ansipurple">2062</span>|  <span class="ansipurple">2067</span>|  <span class="ansiyellow">2075</span>|  <span class="ansired">2081</span>|  <span class="ansired">2089</span>|  <span class="ansired">2096</span>| ---- |  <span class="ansired">2114</span>
   TR2:| ---- | ---- | ---- | ---- | ---- |   <span class="ansired">155</span>|   <span class="ansipurple">160</span>|   <span class="ansired">165</span>|   <span class="ansired">184</span>|   <span class="ansired">191</span>|   <span class="ansired">196</span>|   <span class="ansired">203</span>|   <span class="ansiyellow">210</span>|   <span class="ansired">215</span>| ---- |   <span class="ansipurple">237</span>|   <span class="ansired">251</span>|   <span class="ansipurple">263</span>|   <span class="ansipurple">269</span>|   <span class="ansired">283</span>|   <span class="ansired">289</span>|   <span class="ansipurple">294</span>|   <span class="ansipurple">299</span>|   <span class="ansired">334</span>|   <span class="ansired">341</span>|   355</span>|   <span class="ansired">361</span>|   <span class="ansiyellow">367</span>|   <span class="ansiyellow">372</span>|   <span class="ansired">379</span>|   <span class="ansired">387</span>|   <span class="ansipurple">395</span>|   <span class="ansipurple">405</span>|   <span class="ansired">448</span>|   <span class="ansired">454</span>|   <span class="ansipurple">461</span>| ---- |   <span class="ansired">477</span>|   <span class="ansired">483</span>|   497</span>|   502</span>|   <span class="ansired">508</span>|   <span class="ansired">517</span>|   527</span>|   532</span>|   <span class="ansired">542</span>|   <span class="ansired">547</span>|   <span class="ansired">553</span>|   <span class="ansipurple">558</span>|   <span class="ansired">569</span>|   <span class="ansipurple">587</span>|   <span class="ansiyellow">596</span>|   <span class="ansiyellow">601</span>|   <span class="ansired">629</span>|   <span class="ansired">644</span>|   <span class="ansired">652</span>|   <span class="ansired">669</span>|   <span class="ansipurple">678</span>|   703</span>| ---- |   <span class="ansipurple">711</span>|   <span class="ansipurple">716</span>|   721</span>|   <span class="ansipurple">728</span>|   <span class="ansipurple">733</span>|   <span class="ansired">746</span>|   <span class="ansired">752</span>|   <span class="ansired">757</span>|   <span class="ansired">772</span>|   <span class="ansired">777</span>| ---- |   789|   <span class="ansired">794</span>|   <span class="ansired">802</span>|   <span class="ansired">808</span>|   <span class="ansipurple">817</span>|   <span class="ansired">822</span>|   <span class="ansired">830</span>|   <span class="ansipurple">835</span>|   <span class="ansired">849</span>|   <span class="ansipurple">856</span>|   <span class="ansired">864</span>|   <span class="ansired">872</span>|   <span class="ansiyellow">878</span>|   <span class="ansired">883</span>|   <span class="ansired">892</span>|   <span class="ansiyellow">901</span>|   <span class="ansiyellow">906</span>|   <span class="ansiyellow">912</span>|   <span class="ansired">931</span>|   <span class="ansired">936</span>|   <span class="ansired">945</span>|   <span class="ansiyellow">973</span>|   <span class="ansired">981</span>|   <span class="ansired">987</span>|  <span class="ansired">1003</span>|  <span class="ansipurple">1009</span>|  <span class="ansired">1014</span>|  <span class="ansired">1036</span>|  <span class="ansired">1050</span>|  <span class="ansiyellow">1063</span>|  <span class="ansiyellow">1069</span>|  <span class="ansiyellow">1075</span>|  <span class="ansired">1084</span>| ---- |  1094|  <span class="ansired">1102</span>|  <span class="ansired">1108</span>|  <span class="ansipurple">1113</span>|  <span class="ansipurple">1120</span>|  <span class="ansired">1126</span>|  <span class="ansiyellow">1139</span>|  <span class="ansired">1145</span>|  <span class="ansiyellow">1174</span>|  1180</span>|  <span class="ansipurple">1189</span>|  <span class="ansipurple">1197</span>|  <span class="ansipurple">1202</span>| ---- |  <span class="ansired">1210</span>|  <span class="ansired">1222</span>|  <span class="ansired">1246</span>|  <span class="ansired">1265</span>|  <span class="ansiyellow">1271</span>|  1277</span>|  <span class="ansired">1289</span>|  <span class="ansired">1298</span>|  <span class="ansiyellow">1326</span>|  <span class="ansiyellow">1331</span>|  <span class="ansiyellow">1336</span>|  <span class="ansiyellow">1342</span>|  <span class="ansipurple">1349</span>|  <span class="ansired">1354</span>|  <span class="ansired">1361</span>|  <span class="ansiyellow">1367</span>|  <span class="ansired">1375</span>|  <span class="ansired">1401</span>|  <span class="ansired">1406</span>|  <span class="ansired">1427</span>|  <span class="ansired">1432</span>|  <span class="ansired">1449</span>|  <span class="ansipurple">1454</span>|  <span class="ansipurple">1459</span>|  <span class="ansired">1467</span>|  <span class="ansired">1474</span>|  <span class="ansired">1500</span>| ---- |  <span class="ansipurple">1516</span>|  <span class="ansired">1525</span>|  <span class="ansiyellow">1533</span>|  <span class="ansiyellow">1538</span>|  <span class="ansired">1559</span>|  <span class="ansired">1565</span>|  <span class="ansired">1612</span>|  <span class="ansired">1621</span>|  <span class="ansired">1646</span>|  <span class="ansired">1652</span>|  <span class="ansiyellow">1681</span>| ---- |  <span class="ansired">1690</span>|  <span class="ansired">1696</span>|  <span class="ansipurple">1706</span>| ---- |  <span class="ansipurple">1717</span>|  <span class="ansired">1726</span>|  <span class="ansired">1736</span>|  <span class="ansired">1741</span>|  <span class="ansired">1748</span>|  <span class="ansipurple">1754</span>|  <span class="ansired">1760</span>|  <span class="ansired">1795</span>|  <span class="ansiyellow">1801</span>|  <span class="ansiyellow">1806</span>|  <span class="ansipurple">1838</span>|  <span class="ansired">1844</span>|  <span class="ansired">1849</span>|  <span class="ansired">1858</span>|  <span class="ansired">1868</span>|  <span class="ansipurple">1875</span>|  <span class="ansired">1880</span>|  1885</span>|  1890</span>|  <span class="ansiyellow">1895</span>|  <span class="ansired">1924</span>|  <span class="ansired">1929</span>|  <span class="ansired">1940</span>|  1946</span>|  1953</span>|  1958</span>|  1963</span>|  <span class="ansired">1971</span>|  <span class="ansired">1978</span>|  <span class="ansired">1985</span>|  <span class="ansired">1994</span>|  <span class="ansired">1999</span>|  <span class="ansired">2007</span>|  <span class="ansired">2018</span>|  <span class="ansipurple">2024</span>|  <span class="ansired">2040</span>|  <span class="ansiyellow">2045</span>|  <span class="ansiyellow">2052</span>|  <span class="ansiyellow">2057</span>|  <span class="ansired">2062</span>|  <span class="ansired">2067</span>|  <span class="ansired">2075</span>|  <span class="ansired">2081</span>|  <span class="ansired">2089</span>|  <span class="ansired">2096</span>| ---- |  <span class="ansired">2114</span>
   BR :|     <span class="ansired">7</span>|    <span class="ansiyellow">15</span>|    <span class="ansiyellow">23</span>|    <span class="ansired">28</span>|    <span class="ansired">33</span>|   <span class="ansired">155</span>|   <span class="ansired">160</span>|   <span class="ansired">165</span>|   <span class="ansired">184</span>|   <span class="ansired">191</span>|   <span class="ansired">196</span>|   <span class="ansired">203</span>| ---- |   <span class="ansired">215</span>|   <span class="ansired">229</span>|   <span class="ansired">237</span>|   <span class="ansipurple">251</span>|   <span class="ansired">263</span>|   <span class="ansired">269</span>|   <span class="ansired">283</span>|   <span class="ansired">289</span>|   <span class="ansired">294</span>| ---- |   <span class="ansired">334</span>|   <span class="ansired">341</span>| ---- |   <span class="ansired">361</span>| ---- |   <span class="ansiyellow">371</span>|   <span class="ansired">379</span>|   <span class="ansired">387</span>|   <span class="ansired">395</span>|   <span class="ansired">405</span>|   <span class="ansired">448</span>|   <span class="ansiyellow">454</span>|   <span class="ansicyan">460</span>|   <span class="ansiyellow">466</span>|   <span class="ansired">477</span>|   <span class="ansired">483</span>|   <span class="ansired">498</span>| ---- |   <span class="ansired">508</span>|   <span class="ansired">517</span>|   <span class="ansired">527</span>|   <span class="ansired">533</span>|   <span class="ansired">542</span>|   <span class="ansired">547</span>|   <span class="ansired">553</span>|   <span class="ansipurple">558</span>|   <span class="ansired">569</span>|   <span class="ansired">587</span>|   <span class="ansired">597</span>|   <span class="ansired">602</span>|   <span class="ansired">629</span>|   <span class="ansired">644</span>|   <span class="ansired">652</span>|   <span class="ansired">669</span>|   <span class="ansired">678</span>|   <span class="ansipurple">701</span>|   <span class="ansiyellow">706</span>|   <span class="ansipurple">711</span>|   <span class="ansired">716</span>|   <span class="ansiyellow">721</span>|   <span class="ansired">728</span>|   <span class="ansired">733</span>|   <span class="ansired">746</span>|   <span class="ansipurple">752</span>|   <span class="ansipurple">757</span>|   <span class="ansipurple">773</span>|   <span class="ansipurple">778</span>|   <span class="ansipurple">784</span>|   <span class="ansiyellow">789</span>|   <span class="ansiyellow">794</span>|   <span class="ansired">802</span>|   <span class="ansired">808</span>|   <span class="ansired">817</span>|   <span class="ansired">822</span>|   <span class="ansipurple">830</span>|   <span class="ansipurple">835</span>|   <span class="ansired">849</span>|   856</span>|   <span class="ansired">864</span>|   <span class="ansired">872</span>|   878</span>|   <span class="ansired">883</span>|   <span class="ansired">892</span>|   <span class="ansired">901</span>|   <span class="ansipurple">906</span>|   <span class="ansired">912</span>|   <span class="ansired">931</span>|   <span class="ansired">940</span>|   <span class="ansired">945</span>|   <span class="ansiyellow">974</span>|   <span class="ansired">981</span>|   <span class="ansired">987</span>|  <span class="ansired">1003</span>| ---- |  1014</span>|  <span class="ansired">1036</span>|  <span class="ansired">1050</span>|  <span class="ansired">1063</span>|  <span class="ansired">1069</span>|  <span class="ansired">1075</span>|  <span class="ansired">1084</span>|  1089|  <span class="ansiyellow">1094</span>|  <span class="ansired">1102</span>|  <span class="ansired">1108</span>|  <span class="ansired">1113</span>|  <span class="ansired">1120</span>|  <span class="ansired">1126</span>|  <span class="ansipurple">1139</span>|  <span class="ansired">1145</span>|  <span class="ansired">1173</span>|  <span class="ansired">1181</span>|  <span class="ansired">1189</span>| ---- |  <span class="ansipurple">1200</span>|  <span class="ansipurple">1205</span>|  <span class="ansired">1210</span>|  <span class="ansired">1222</span>|  <span class="ansired">1246</span>|  <span class="ansired">1265</span>|  <span class="ansired">1270</span>|  <span class="ansired">1275</span>|  <span class="ansiyellow">1289</span>|  <span class="ansired">1298</span>|  <span class="ansired">1325</span>|  <span class="ansired">1330</span>|  <span class="ansired">1335</span>|  <span class="ansired">1342</span>|  <span class="ansired">1349</span>|  <span class="ansired">1354</span>|  <span class="ansired">1361</span>|  1366</span>|  <span class="ansired">1375</span>|  <span class="ansiyellow">1401</span>|  <span class="ansired">1406</span>|  <span class="ansired">1427</span>|  <span class="ansired">1432</span>|  <span class="ansired">1449</span>|  <span class="ansired">1454</span>|  <span class="ansired">1459</span>|  <span class="ansired">1467</span>|  <span class="ansired">1474</span>|  <span class="ansired">1500</span>| ---- |  <span class="ansired">1516</span>|  <span class="ansiyellow">1525</span>|  <span class="ansired">1533</span>|  <span class="ansired">1538</span>|  1558</span>|  1565</span>|  <span class="ansired">1612</span>|  <span class="ansired">1621</span>|  <span class="ansipurple">1646</span>|  <span class="ansired">1652</span>|  <span class="ansiyellow">1681</span>| ---- |  <span class="ansired">1690</span>|  <span class="ansired">1696</span>|  <span class="ansiyellow">1704</span>|  <span class="ansiyellow">1709</span>|  <span class="ansired">1717</span>|  <span class="ansired">1726</span>|  <span class="ansired">1736</span>|  <span class="ansired">1741</span>|  <span class="ansired">1748</span>|  <span class="ansiyellow">1754</span>|  <span class="ansired">1760</span>|  <span class="ansired">1795</span>|  <span class="ansired">1800</span>| ---- | ---- |  <span class="ansired">1844</span>|  <span class="ansired">1849</span>|  <span class="ansired">1858</span>|  <span class="ansired">1868</span>|  <span class="ansiyellow">1875</span>|  <span class="ansiyellow">1880</span>|  1886| ---- |  <span class="ansired">1895</span>|  <span class="ansired">1924</span>|  <span class="ansired">1929</span>|  <span class="ansired">1940</span>|  <span class="ansiyellow">1946</span>|  <span class="ansiyellow">1953</span>|  <span class="ansiyellow">1958</span>|  <span class="ansiyellow">1963</span>|  <span class="ansired">1971</span>|  <span class="ansired">1978</span>|  <span class="ansired">1985</span>|  <span class="ansired">1994</span>|  <span class="ansired">1999</span>|  <span class="ansired">2007</span>|  <span class="ansiyellow">2018</span>|  <span class="ansiyellow">2024</span>|  <span class="ansiyellow">2040</span>|  2045</span>| ---- |  2057</span>|  <span class="ansipurple">2062</span>|  <span class="ansiyellow">2067</span>|  <span class="ansired">2074</span>|  <span class="ansired">2081</span>|  <span class="ansired">2089</span>|  <span class="ansired">2096</span>|  <span class="ansired">2106</span>|  <span class="ansired">2114</span>
   
   </pre>
    

The region between boundaries 405 and 448 (430Kb) has well conserved boundaries (score of 9 and 10) and will 
therefore be selected for the modeling:

::

    # Get the TAD number for each of the selected TADs in each experiment using the brk (upper boundary)
    brkp = 448
    ntad = {}
    ltad = {}
    for xnam in datasets:
        ltad[xnam] = my_chrom.experiments[xnam].tads
        for tad_k, tad_v in ltad[xnam].iteritems():
            for k, v in tad_v.iteritems():
                if k == 'brk' and v == brkp:
                    ntad[xnam] = tad_k
                    print "Experiment "+ xnam + " has selected TAD " + str(ntad[xnam])


.. parsed-literal::

    Experiment TR1 has selected TAD 25
    Experiment TR2 has selected TAD 27
    Experiment BR has selected TAD 29


::

    # alternative way to do the same
    # I know you should have passed a hard time to write your loop, but I think this way is more elegant
    # and it works even if one of the TADs does not end exactly at 448 (it could be that the column is 448 448 449).
    brkp = 448
    ntad = {}
    columns = ali.get_column(lambda x: x['brk'] == brkp, 
                             min_num=1) # selection of columns having TADs ending at genomic bin '448'
                                        # Can also be done by counting columns, this one is the 32: 
                                        # >>> ali.get_column(32)
    column= columns[0]                  # only one column
    for tad in column[1]:               # second element of a column is the list of TADs (first is the column number)
        ntad[tad['exp'].name] = tad['index']
        print "Experiment "+ tad['exp'].name + " has selected TAD " + str(tad['index'])


.. parsed-literal::

    Experiment TR1 has selected TAD 25
    Experiment TR2 has selected TAD 27
    Experiment BR has selected TAD 29

Compare all-against-all the selected tads:

::

    # We can now compare the selected TADs
    for i in range(0, len(datasets)):
        for j in range(i+1, len(datasets)):
            xnam1 = datasets[i]        
            xnam2 = datasets[j]
            hic1 = my_chrom.get_tad_hic(my_chrom.experiments[xnam1].tads[ntad[xnam1]],xnam1)
            hic2 = my_chrom.get_tad_hic(my_chrom.experiments[xnam2].tads[ntad[xnam2]],xnam2)
            align1, align2, score = optimal_cmo(hic1, hic2, max_num_v=8, long_nw=True, long_dist=True, method='frobenius')
            print "{:>3}[{}] .vs. {:>3}[{}] CC: {:.3}".format(xnam1, ntad[xnam1], xnam2, ntad[xnam2], score['rho'])


.. parsed-literal::

    TR1[25] .vs. TR2[27] CC: 0.949
    TR1[25] .vs.  BR[29] CC: 0.93
    TR2[27] .vs.  BR[29] CC: 0.945

All the comparisons have a high correlation coefficient (>0.92); thus, the three experiments for the selected TADs can be merged into a single experiment:

::

    # Sum all datasets into a single experimetal set
    exp = my_chrom.experiments['BR'] + my_chrom.experiments['TR1'] + my_chrom.experiments['TR2']
    my_chrom.experiments.append(exp)
    exp


.. parsed-literal::

    Experiment BR+TR1+TR2 (resolution: 10Kb, TADs: None, Hi-C rows: 2115)


Re-calculate all the TADs for the sum of the experiments:

::

    # Re-calculate the TAD borders for the new dataset with the sum of the individual matrices
    xnam = 'BR+TR1+TR2'
    my_chrom.find_tad(xnam) # Finding TADs can take significant CPU time (by default it uses all available CPUs)
    my_chrom.visualize(xnam, paint_tads=True)

.. image:: pictures/Tadbit_for_IMP_fransua_25_0.png


::

    # Is the selected TAD also found in the new TADs from matrix joining the three experiments?
    strp = 406 # starting point
    endp = 448 # end point
    ltad[xnam] = my_chrom.experiments[xnam].tads
    for tad_k, tad_v in ltad[xnam].iteritems():
        if tad_v['start'] == strp and tad_v['end'] == endp:
            ntad[xnam] = tad_k
            print "Experiment "+ xnam + " has selected TAD " + str(ntad[xnam])


.. parsed-literal::

    Experiment BR+TR1+TR2 has selected TAD 24

The selected TAD (TAD 24) is also found in the new TADs from the matrix joining the three experiments. 
Visualize the selected TAD:

***** 
    I ALSO WORRY FOR THE BLUE LINES THAT GO TO THE DIAGONAL. ARE THOSE CORRECT? MISSING VALUES? => yes, data is a bit wired...
*****

::

    # Visualize a specific TAD
    my_tad = my_chrom.experiments[xnam].tads[ntad[xnam]]
    my_chrom.visualize(xnam, tad=my_tad)

.. image:: pictures/Tadbit_for_IMP_fransua_28_0.png


.. parsed-literal::

    <matplotlib.image.AxesImage at 0x3f5f290>



3. Model the 3D structure of a selected TAD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Run the 3D modeling on the selected TAD. The generated models are stored in a dictionary whose keys are consecutive integers, ordered by the IMP Objective Function value (smaller to larger). All the models are based on specific sets of experimental data for which TADBit modeling parameters need to be optimized (see Tutorial #2). Optimizing the parameters takes significant CPU time and thus we have pre-optimized them for the datasets in this tutorial. The specific parameters are stored in a python dictionary called CONFIG. The one used in this tutorial are CONFIG['dmel_01'].

******
    IF I USE ONLY 406 to 448 I GET AN ERROR ON nmin, SEEMS THAT CANNNOT MODEL SMALL REGIONS? TO CORRECT! => indeed it is the step of "polynomial-fit" filtering that was craching, this steps needs a minimum of columns to be able to fit something... This step is now skipped when number of columns is too low.
****** 

::

    # 3D modeling configuration.
    CONFIG


.. parsed-literal::

    {'dmel_01': {'kforce': 5,
      'lowfreq': -0.7,
      'lowrdist': 100,
      'maxdist': 600,
      'reference': 'victor corces dataset 2013',
      'upfreq': 0.3}}



::

    # Build 3D models based on the Hi-C data. This is done by IMP. 
    models = exp.model_region(406, 448, n_models=500, n_keep=100, n_cpus=8, keep_all=True, config=CONFIG['dmel_01'])
    print models


.. parsed-literal::

    /usr/local/lib/python2.7/dist-packages/pytadbit/utils/hic_filtering.py:145: UserWarning: WARNING: Too few data to filter columns. SKIPPING...
      warn('WARNING: Too few data to filter columns. SKIPPING...')


.. parsed-literal::

    StructuralModels with 100 models (objective function range: 922965 - 932558)
       (corresponding to the best models out of 500 models).
      IMP modelling used this parameters:
       - maxdist     : 600
       - upfreq      : 0.3
       - reference   : victor corces dataset 2013
       - kforce      : 5
       - lowfreq     : -0.7
       - lowrdist    : 100
      Models where clustered into 0 clusters

Select a subset of all the generated models:

::

    # Select top 10 models
    models.define_best_models(10)
    print "Lowest 10 IMP OF models:" 
    print models
    
    # Select top 50 models
    models.define_best_models(100)
    print "Lowest 50 IMP OF models:" 
    print models


.. parsed-literal::

    Lowest 10 IMP OF models:
    StructuralModels with 10 models (objective function range: 922965 - 923673)
       (corresponding to the best models out of 500 models).
      IMP modelling used this parameters:
       - maxdist     : 600
       - upfreq      : 0.3
       - reference   : victor corces dataset 2013
       - kforce      : 5
       - lowfreq     : -0.7
       - lowrdist    : 100
      Models where clustered into 0 clusters
    Lowest 50 IMP OF models:
    StructuralModels with 100 models (objective function range: 922965 - 932558)
       (corresponding to the best models out of 500 models).
      IMP modelling used this parameters:
       - maxdist     : 600
       - upfreq      : 0.3
       - reference   : victor corces dataset 2013
       - kforce      : 5
       - lowfreq     : -0.7
       - lowrdist    : 100
      Models where clustered into 0 clusters

For each model, the IMP objective function value, the random number generator feed used in IMP, the xyz coordinates and the log file of the search for the best conformation (i.e. the 3D conformation that best satisfies the input restraints - corresponding to the lowest IMP objective function) are stored. Data from a model can be retreived as followed:


::

    # Get the data for the lowest IMP OF model (number 0) in the set of models
    model = models[0]
    print model


.. parsed-literal::

    IMP model of 43 particles with: 
     - Final objective function value: 922965.537593
     - random initial value: 479
     - first coordinates:
            X      Y      Z
          119   -500    282
           37   -476    212
          136   -459    187
    

The IMP objective function can be plotted as function of each optimization step during the modeling:


::

    # Get the IMP OF of the stored model in "model"
    model.objective_function(log=True, smooth=True)

.. image:: pictures/Tadbit_for_IMP_fransua_38_0.png

It is important to verify whether the set of selected models (usually the top 20%, based on their IMP objective function value) has a good correlation with the input Hi-C data. High correlation values mean that the select set of 3D models represent well the input experimental data. The correlation is calculated with the following function:

::

    # Calculate the correlation coefficient between the selected models and the original Hi-C matrix
    models.correlate_with_real_data(models=range(50), plot=True, cutoff=300)

.. image:: pictures/Tadbit_for_IMP_fransua_40_0.png


4. Model analysis
~~~~~~~~~~~~~~~~~



4.1 Model clustering
^^^^^^^^^^^^^^^^^^^^

The selected models are clustered based on their structural similarity. Clusters are numbered from larger (more models) to smaller (less models):

::

    # Cluster models based on structural similarity
    models.cluster_models(fact=0.90, dcutoff=200)
    print models.clusters


.. parsed-literal::

    {0: [57, 39, 70, 43, 63, 44, 42, 12, 9, 0, 15, 30, 62, 17, 8, 28, 75, 35, 38, 56, 34, 54, 36, 40, 73, 1, 45, 59, 14, 64, 58, 55, 53, 51, 27, 2], 1: [50, 41, 88, 84, 77, 89, 61, 67, 94, 81, 79, 46, 82, 47, 71, 68, 48, 60, 66, 95], 2: [3, 76, 22, 5, 4, 65, 13, 23, 7, 83, 78, 21, 24], 3: [29, 33, 32, 37, 25, 20], 4: [80, 86, 87, 74, 69, 72], 5: [98, 91, 85, 93], 6: [10, 11, 16, 6], 7: [97, 92, 99], 8: [90, 96], 9: [52, 49], 10: [19, 18], 11: [31, 26]}

Generated clusters can be plotted for ease of visualization. The "y" axis of the plot shows the IMP objective function. The largest cluster (the one numbered with a "0") is expected to have the lowest IMP objective function:

******
     NEED TO SAVE THE CLUSTERING AND DENDOGRAM DATA INTO A FILE. => format? always or on demand?
******

::

    # Plot the resulting clusers
    cl = models.cluster_analysis_dendrogram(color=True)

.. image:: pictures/Tadbit_for_IMP_fransua_46_0.png

Similarity betwen clusters can be shown for a limited number of them (5 in this example):

::

    # Show de dendogram for only the 5 top clusters and no colors
    cl = models.cluster_analysis_dendrogram(n_best_clusters=5)

.. image:: pictures/Tadbit_for_IMP_fransua_48_0.png


4.2 Models consistency
^^^^^^^^^^^^^^^^^^^^^^


To assess how "deterministic" a cluster is, for each particle the percentage of models (in the cluster) that superimpose the particle within a given cut-off (pre-set cut-offs of 50, 100, 150 and 200 nm) can be calculated. The lower the consistency value (in %) the less deterministic (i.e. more variable) the corresponfing region of models within the selected cluster. This measure can be taken as a proxy of variability across the model:

***** 
    NEED TO SAVE THE CONSISTENCY DATA INTO A FILE. => Done, on demand only
*****

::

    # Calculate the consistency plot for all models in the first cluster (cluster 0)
    models.model_consistency(cluster=0)

.. image:: pictures/Tadbit_for_IMP_fransua_51_0.png


4.3 DNA density plots
^^^^^^^^^^^^^^^^^^^^^


From the 3D models, the DNA density (or local compactness, in bp/nm) can be calculated as the ratio of the bin size (in base pairs) and the distances (nm) between consequtive particles in the models. The higher the density the more compact DNA for the region. Values are calculated using running averages:

******
     NEED TO SAVE THE DENSITY DATA INTO A FILE. => done on demand only
******

::

    # Calculate a DNA density plot
    models.density_plot()

.. image:: pictures/Tadbit_for_IMP_fransua_54_0.png

***** 
    I NEED TO UNDESTAND HOW YOU GET THIS STANDARD DEVIATION (MAYBE USE STANDAR ERROR?). => It's simply two times the standard deviation, but it moves a lot from one model to the other.
    ALSO, I DID TRIED TO USE ONLY ONE STEPS (10) AND DID NOT WORK. WHY? => FIXED. before only tuples were accepted.
*****

::

    # Get a similar plot for only the top cluster and show the standar deviation for a specific running window (steps)
    models.density_plot(error=True, steps=(2, 5))

.. image:: pictures/Tadbit_for_IMP_fransua_56_0.png


4.4 Models contact map
^^^^^^^^^^^^^^^^^^^^^^


Given a set of selected models (either from a cluster or a list), the percentage of pairs of particles within a for a distance cut-off can be calculated. The result can then be represented as a heat-map similar to the ones used to represent Hi-C data:

******
    ALSO REMOVE THE Z LABEL => all?
    NEED TO SAVE THE CONTACT MAP DATA INTO A FILE. => done.
******

::

    # Get a contact map for the top 50 models at a distance cut-off of 300nm
    models.contact_map(models=range(50), cutoff=300)

.. image:: pictures/Tadbit_for_IMP_fransua_59_0.png

The goal of IMP is to find a 3D structure (or an ensemble of structures) that best satisfies the original Hi-C matrix. To verify the how well the generated models represent the input experimental data, the contact map produced above can be compared to the original Hi-C input matrix:

::

    # Correlate the contact map with the original input Hi-C matrix
    models.correlate_with_real_data(models=range(50), plot=True, cutoff=300)

.. image:: pictures/Tadbit_for_IMP_fransua_61_0.png


4.5 Calculating distances between particles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Sometimes is useful to get the distribution of distances between pairs of particles in the models (or in a sub-set of the models): 

******
    WE ALSO NEED TO CALCULATE ANGLES BETWEEN THREE PARTICLES.
    WOULD BE NICE TO ALSO HAVE THE MEDIAN FOR A DISTANCE (NOT JUST THE AVERAGE). => hmmm the name of the function is bad, what is returned is the median... you can also get the full list of value if you want to compute the mean for example. I haven't done the mean. But if you wish I can.
    ALL OUTPUT NEEDS TO HAVE A MAXIMUM OF TWO DECIMALS. => At the time to write them imto a file OK. But for the return of a function, I do not think this is a good idea, if you want it you just do "round(result, 3)"
    THERE IS A NEED ALSO TO HAVE AN OPTION FOR SAVING THE ACTUAL DISTANCES IN A FILE. => format? usually the format is a three columns file: P1 P2 Distance
******

::

    # Get the average distance between particles 13 and 23 in all kept models
    models.average_3d_dist(13, 23, plot=False)


.. parsed-literal::

    353.18443338679435


Plot of the distribution used to get the mean value:

::

    # Plot the distance distributions between particles 13 and 23 in all kept models
    models.average_3d_dist(13, 23, plot=True)

.. image:: pictures/Tadbit_for_IMP_fransua_66_0.png

To use only the 10 first models (lowest energy), or the models belonging to a cluster (example cluster 0):

::

    # Plot the distance distributions between particles 13 and 23 in the top 10 models
    models.average_3d_dist(13, 23, models=range(10))

.. image:: pictures/Tadbit_for_IMP_fransua_68_0.png


::

    # Plot the distance distributions between particles 13 and 23 in the models from cluster 0
    models.average_3d_dist(13, 23, plot=True, cluster=0)

.. image:: pictures/Tadbit_for_IMP_fransua_69_0.png


4.6 Visualizing the 3D models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Specific 3D models can be saved in two formats:
    - CMM format, which can be directly loaded into Chimera for visualization.
    - XYZ format, which is a simple format that can be useful for further analysis that require coordinates.

******
    CAN GET GET THE OPTION TO WRITE CMMs AND XYZ FILES FOR A RANGE OF MODELS OR ALL MODELS IN A CLUSTER?
    WE NEED TO PREPARE A PAIR OF FUNCTIONS THAT WRITE MODELS AND SCRIPTS FOR CHIMERA
        1. SINGLE MODEL TO VISUALIZE AS A "tube"
        2. ENTIRE CLUSTER TO VISUALIZE AS "surface"
    DAVIDE CAN HELP IN PREPARING THE CHIMERA PYTHON SCRIPTS => YES!!! :)
******

::

    # Write a CMM file for the top model
    models.write_cmm(directory="./", model_num=0, models=None, cluster=None)
    # Write a XYZ file for the top model
    models.write_xyz(directory="./", model_num=0, models=None, cluster=None)


Save and load your analysis
---------------------------


By saving the analysis, some of the most expensive calculations will not need to be repeated:

::

    # Save your entire analysis and models
    models.save_models('dmel_01.models')

To load them:

::

    # Load the models
    models = load_structuralmodels('dmel_01.models')
    print models


.. parsed-literal::

    StructuralModels with 100 models (objective function range: 922965 - 932558)
       (corresponding to the best models out of 500 models).
      IMP modelling used this parameters:
       - maxdist     : 600
       - upfreq      : 0.3
       - reference   : victor corces dataset 2013
       - kforce      : 5
       - lowfreq     : -0.7
       - lowrdist    : 100
      Models where clustered into 12 clusters

