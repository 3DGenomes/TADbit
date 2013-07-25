Three-dimensional modeling
**************************

.. contents::
   :depth: 3


Data filtering
==============

To model chromatin structure, we need to ensure that our data is clean enough. The first step is thus to draw the distribution of the sum of interactions per raw/columns in the Hi-C matrix. According to this distribution, we may remove some columns if they present a suspiciously low count of interaction.

Here an example, where "exp" is an preloaded Experiment corresponding to human's 19th chromosome:

::

  from pytadbit.utils import filter_with_polynomial_fit

  filter_with_polynomial_fit(exp.get_hic_matrix(), draw_hist=True)

.. figure::  pictures/example_filtering.png
   :align:   center

Than, according to the fit represented above, we would discard all columns in the Hi-C raw data having cumulative count of interaction below the dashed red line in the graph above (~46). This columns will be removed from the modeling, and their associated particles will have no experimental data.

*This step is done automatically within tadbit each time an experiment is loaded. In order to ensure that we do remove outlier columns, tadbit checks if this root corresponds to a* **concave down** *region and if it stands* **between zero and the median** *of the overall distribution. The result of these "bad" columns is stored in the variable Experiment._zeros, that represents the columns to be skipped in the consecutive steps.*


Data normalization
==================

Hi-C data stored in :class:`pytadbit.experiment.Experiment` might be normalized in order to be used by IMP.
This normalization is achieve in two steps, first we generate weight for each pair of interactions, depending on the interaction count in the corresponding row/column, second we calculate the `z-score <http://en.wikipedia.org/wiki/Standard_score#Calculation_from_raw_score>`_ of each of these interaction pairs.

Calculation of weights
----------------------

Weights can be calculated according to two formulas (see :class:`pytadbit.experiment.Experiment.normalize_hic`), however, in the context of three-dimensional modeling, the "over_tot" method is recommended, as the distribution of values generated is closer to normal.

Hi-C interaction count are thus normalized according to this formula:

.. math::

  weight(I, J) = \frac{\sum^N_{i=0}{(matrix(i, J))} \times \sum^N_{j=0}{(matrix(I, j))}}{\sum^N_{i=0}{\sum^N_{j=0}{(matrix(i, j))}}}


"matrix", being our row data (count of interactions), N the number of rows/columns.

The result is stored in a new matrix, called weight. The values that will be used in the next step are the multiplication of this weights per the raw data.


Calculation of the z-score
--------------------------

Z-scores are computed according to classical formula (:math:`\frac{x-\mu}{\sigma}`), over the decimal logarithm values of the normalized data (see above). Ending in this formula:

.. math::

  zscore(I, J) = \frac{log_{10}(weight(I, J) \times matrix(I, J)) - mean(log_{10}(weight \times matrix))}{stddev(log_{10}(weight \times matrix))}

**Important: values on the diagonal are not taken into account for this calculus.**

Dealing with zeros
^^^^^^^^^^^^^^^^^^

A zero in an Hi-C interaction matrix, means that the given two fragments of DNA were never found close enough to be crosslinked together. However such values are also highly suspicious to be artifacts. 

Right now we assume that :math:`log_{10}(0) = 0`, in the calculation of the mean and stddev, and equal to -1 in the calculation of the z-score itself.



How to get ThreeDeeModels?
==========================


Here we load a Chromosome object, from which we take one Experiment object ('exp'). 

From this Experiment object we can model a given region using IMP.

The result of these few line is a ThreeDeeModels object, which will have all the function you asked me to implement (... yes will... like, in the future :S )
In[1]:

.. code:: python

    from pytadbit import Chromosome

I define my chromosome
In[2]:

.. code:: python

    crm = '2R'
    crmbit = Chromosome('2R')

I load all experiments done on Drosophila's chromosome 2R (Hi-C matrices), and sum the Hi-C matrices (Corces' technical and biolobical replicates) into a single experiment
In[3]:

.. code:: python

    for xnam in ['TR2', 'TR1', 'BR']:
        crmbit.add_experiment(xnam, resolution=10000, 
                              xp_handler='/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
    
    exp = crmbit.experiments['TR1'] + crmbit.experiments['TR2'] + crmbit.experiments['BR']

Finally run the IMP modelling on a given region (this region crresponds to the one Davide shows at meeting with Guillaume)
In[4]:

.. code:: python

    models = exp.model_region(190, 295, n_models=500, n_keep=200, n_cpus=8)


.. parsed-literal::

    processing model #100
    processing model #200
    processing model #300
    processing model #500
    processing model #400


Playing around with models
--------------------------


models are stored in a dictionary which keys are number (the lowest the less energy).
Thus to have a look to the best model we just type:
In[5]:

.. code:: python

    print models


.. parsed-literal::

    ThreeDeeModels with 200 models (energy range: 4135345.18085-4203587.19953)
       (corresponding to the best models out of 200 models).
      Models where clustered into 0 clusters

Thus for each model is stored, the final energy, the random initial number used with IMP, the coordinates xyz and the log of the search for the best conformation lowering the energy.

Each can be reached like this:

In[6]:

.. code:: python

    model = models[0]
    print model



.. parsed-literal::

    IMP model of 105 particles with: 
     - Final energy: 4135345.18085
     - random initial value: 430
     - first coordinates:
            X      Y      Z
          438    380    331
          426    357    357
          411    337    375
    


Objective function
~~~~~~~~~~~~~~~~~~


We want to plot the objective function for this best model:
In[7]:

.. code:: python

    models.objective_function_model(0, log=False, smooth=False)

.. image:: pictures/Tadbit_for_IMP_notebook_16_0.png

... perhaps nicer with log (note that it can be done using the IMPmodel object directely):
In[8]:

.. code:: python

    model = models[0]
    model.objective_function(log=True, smooth=True)

.. image:: pictures/Tadbit_for_IMP_notebook_18_0.png


Clustering models
-----------------


First we run the clustering. The result of this will be stored inside the ThreeDeeModels object.
In[9]:

.. code:: python

    models.cluster_models(fact=0.75, dcutoff=200)
    print models.clusters


.. parsed-literal::

    {0: [0, 1, 10, 104, 11, 113, 117, 12, 121, 123, 127, 131, 132, 14, 142, 145, 150, 155, 157, 160, 163, 167, 17, 170, 171, 172, 177, 182, 187, 190, 191, 197, 2, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 32, 33, 34, 38, 4, 40, 41, 43, 44, 46, 48, 5, 52, 6, 60, 61, 62, 68, 7, 71, 74, 77, 8, 85, 86, 88, 89, 9, 91, 92, 93, 95, 97, 99], 1: [101, 107, 108, 109, 110, 112, 114, 115, 116, 118, 119, 120, 122, 124, 125, 126, 129, 130, 133, 134, 135, 136, 137, 138, 139, 140, 141, 161, 179, 185, 49, 51, 59, 63, 66, 69, 75, 76, 79, 80, 84, 87], 2: [13, 15, 16, 169, 173, 18, 19, 192, 196, 20, 30, 31, 36, 37, 42, 45, 47, 55, 56, 57, 67, 72, 78, 82, 90, 94, 96], 3: [143, 147, 148, 151, 154, 158, 159, 162, 164, 166, 168, 175, 176, 180, 181], 4: [146, 184, 193, 50, 53, 58, 64, 70, 73, 81, 83], 5: [100, 102, 103, 105, 106, 111, 98], 6: [128, 149, 152, 153, 178, 189, 65], 7: [183, 186, 188, 195, 198, 199], 8: [144, 156, 165, 174, 194], 9: [35, 39, 54]}


Plot clusters
~~~~~~~~~~~~~


We can plot everything (The 10 clusters found):
In[10]:

.. code:: python

    cl = models.cluster_analysis_dendrogram(color=True)

.. image:: pictures/Tadbit_for_IMP_notebook_24_0.png

Or just 6 of them (without this colors that no one understands...)
In[11]:

.. code:: python

    cl = models.cluster_analysis_dendrogram(n_best_clusters=6)

.. image:: pictures/Tadbit_for_IMP_notebook_26_0.png


Distance between particles
--------------------------



Between 2 particles
~~~~~~~~~~~~~~~~~~~


We can just quickly get a value of the distance between particle 13 and 23
In[12]:

.. code:: python

    models.average_3d_dist(13, 23, plot=False)


.. parsed-literal::

    318.63149011844644


This by default, is calculated over the ensemble of models we have. Lets plot the distribution used to get this mean value:
In[13]:

.. code:: python

    models.average_3d_dist(13, 23, plot=True)

.. image:: pictures/Tadbit_for_IMP_notebook_32_0.png

We may also want to use only the 10 first models, or the models belonging to cluster number 0:
In[14]:

.. code:: python

    models.average_3d_dist(13, 23, models=range(10))

.. image:: pictures/Tadbit_for_IMP_notebook_34_0.png

In[15]:

.. code:: python

    models.average_3d_dist(13, 23, plot=True, cluster=0)

.. image:: pictures/Tadbit_for_IMP_notebook_35_0.png


Density plot
~~~~~~~~~~~~


Using distances between particle, we can plot now the density (bp per nm) of our chromosomic region.
In[16]:

.. code:: python

    models.density_plot(models=None)

.. image:: pictures/Tadbit_for_IMP_notebook_38_0.png

In[17]:

.. code:: python

    models.density_plot(cluster=0, error=True, steps=(5,10))

.. image:: pictures/Tadbit_for_IMP_notebook_39_0.png


Contact Map
~~~~~~~~~~~


In[18]:

.. code:: python

    models.contact_map_consistency(models=None, cluster=None, cutoff=150)

.. image:: pictures/Tadbit_for_IMP_notebook_41_0.png


Consistency Plot
----------------


