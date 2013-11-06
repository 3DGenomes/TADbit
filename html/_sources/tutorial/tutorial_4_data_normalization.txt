
Hi-C data processing
====================


Data filtering
--------------


To model chromatin structure, we need to ensure that our data is clean enough. The first step is thus to draw the distribution of the sum of interactions per raw/columns in the Hi-C matrix. According to this distribution, we may remove some columns if they present a suspiciously low count of interaction.

Here an example, where "exp" is an preloaded Experiment corresponding to human's 19th chromosome:

.. code:: python

    from pytadbit import Chromosome
    from pytadbit.utils.hic_filtering import filter_by_mean
    
    my_chrom = Chromosome('19')
    my_chrom.add_experiment('gm', resolution=10000, 
                            hic_data='../../scripts/sample_data/HIC_gm06690_chr19_chr19_100000_obs.txt')
    
    exp = my_chrom.experiments[0]
    
    zeroes = filter_by_mean(exp.get_hic_matrix(), draw_hist=True)



.. image:: ../nbpictures/tutorial_4_data_normalization_3_0.png


Than, according to the fit represented above, we would discard all columns in the Hi-C raw data having cumulative count of interaction below the dashed red line in the graph above (~46). This columns will be removed from the modeling, and their associated particles will have no experimental data.

*This step is done automatically within tadbit each time an experiment is loaded. In order to ensure that we do remove outlier columns, tadbit checks if this root corresponds to a* **concave down** *region and if it stands* **between zero and the median** *of the overall distribution. The result of these "bad" columns is stored in the variable Experiment._zeros, that represents the columns to be skipped in the consecutive steps.*

Data normalization
------------------


Hi-C data stored in :class:`pytadbit.experiment.Experiment` might be normalized in order to be used by IMP.
This normalization is achieve in two steps, first we generate weight for each pair of interactions, depending on the interaction count in the corresponding row/column, second we calculate the `z-score <http://en.wikipedia.org/wiki/Standard_score#Calculation_from_raw_score>`_ of each of these interaction pairs.

Calculation of weights
~~~~~~~~~~~~~~~~~~~~~~


Weights can be calculated according to two formulas (see :class:`pytadbit.experiment.Experiment.normalize_hic`), however, in the context of three-dimensional modeling, the "over_tot" method is recommended, as the distribution of values generated is closer to normal.

Hi-C interaction count are thus normalized according to this formula:

.. math::

  weight(I, J) = \frac{\sum^N_{i=0}{(matrix(i, J))} \times \sum^N_{j=0}{(matrix(I, j))}}{\sum^N_{i=0}{\sum^N_{j=0}{(matrix(i, j))}}}


"matrix", being our row data (count of interactions), N the number of rows/columns.

The result is stored in a new matrix, called weight. The values that will be used in the next step are the multiplication of this weights per the raw data.

Calculation of the z-score
~~~~~~~~~~~~~~~~~~~~~~~~~~


Z-scores are computed according to classical formula (:math:`\frac{x-\mu}{\sigma}`), over the decimal logarithm values of the normalized data (see above). Ending in this formula:

.. math::

  zscore(I, J) = \frac{log_{10}(weight(I, J) \times matrix(I, J)) - mean(log_{10}(weight \times matrix))}{stddev(log_{10}(weight \times matrix))}

**Important: values on the diagonal are not taken into account for this calculus.**

Dealing with zeros
^^^^^^^^^^^^^^^^^^


A zero in an Hi-C interaction matrix, means that the given two fragments of DNA were never found close enough to be crosslinked together. However such values are also highly suspicious to be artifacts. 

Right now we assume that :math:`log_{10}(0) = 0`, in the calculation of the mean and stddev, and equal to -1 in the calculation of the z-score itself.
