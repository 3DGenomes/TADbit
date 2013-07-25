
How to get ThreeDeeModels?
==========================


Here we load a Chromosome object, from which we take one Experiment object ('exp'). 

From this Experiment object we can model a given region using IMP.

The result of these few line is a ThreeDeeModels object, which will have all the function you asked me to implement (... yes will... like, in the future :S )


.. code:: python

    from pytadbit import Chromosome

I define my chromosome


.. code:: python

    crm = '2R'
    crmbit = Chromosome('2R')

I load all experiments done on Drosophila's chromosome 2R (Hi-C matrices), and sum the Hi-C matrices (Corces' technical and biolobical replicates) into a single experiment


.. code:: python

    for xnam in ['TR2', 'TR1', 'BR']:
        crmbit.add_experiment(xnam, resolution=10000, 
                              xp_handler='/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
    
    exp = crmbit.experiments['TR1'] + crmbit.experiments['TR2'] + crmbit.experiments['BR']

Finally run the IMP modelling on a given region (this region crresponds to the one Davide shows at meeting with Guillaume)


.. code:: python

    models = exp.model_region(190, 295, n_models=500, n_keep=250, n_cpus=8)


.. parsed-literal::

    processing model #100
    processing model #200
    processing model #300
    processing model #400
    processing model #500


Playing around with models
--------------------------


models are stored in a dictionary which keys are number (the lowest the less energy).
Thus to have a look to the best model we just type:


.. code:: python

    print models


.. parsed-literal::

    ThreeDeeModels with 250 models (energy range: 4135345-4212180)
       (corresponding to the best models out of 250 models).
      Models where clustered into 0 clusters

Thus for each model is stored, the final energy, the random initial number used with IMP, the coordinates xyz and the log of the search for the best conformation lowering the energy.

Each can be reached like this:



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
------------------


We want to plot the objective function for this best model:


.. code:: python

    models.objective_function_model(0, log=False, smooth=False)

.. image:: pictures/Tadbit_for_IMP_notebook_16_0.png

... perhaps nicer with log (note that it can be done using the IMPmodel object directely):


.. code:: python

    model = models[0]
    model.objective_function(log=True, smooth=True)

.. image:: pictures/Tadbit_for_IMP_notebook_18_0.png


Clustering models
-----------------


First we run the clustering. The result of this will be stored inside the ThreeDeeModels object.


.. code:: python

    models.cluster_models(fact=0.75, dcutoff=200)
    print models.clusters


.. parsed-literal::

    {0: [0, 1, 10, 11, 113, 117, 12, 121, 123, 127, 131, 132, 14, 142, 145, 150, 155, 157, 160, 163, 167, 17, 170, 171, 172, 177, 182, 187, 19, 190, 191, 197, 2, 21, 212, 214, 219, 22, 226, 228, 23, 24, 246, 25, 26, 27, 28, 29, 3, 32, 33, 34, 36, 38, 4, 40, 41, 42, 43, 44, 45, 46, 48, 5, 52, 56, 6, 60, 61, 62, 67, 68, 7, 71, 72, 74, 77, 8, 85, 86, 88, 89, 9, 91, 92, 93, 94, 95, 97, 99], 1: [101, 107, 108, 109, 110, 112, 114, 115, 116, 118, 119, 120, 122, 124, 125, 126, 128, 129, 130, 133, 134, 135, 136, 137, 139, 140, 141, 161, 179, 185, 189, 49, 51, 59, 63, 66, 69, 75, 76, 79, 80, 84, 87, 90, 96], 2: [144, 146, 169, 173, 174, 184, 192, 193, 194, 200, 206, 208, 209, 210, 215, 220, 222, 225, 227, 230, 231, 233, 237, 239, 240, 241, 244, 37, 50, 53, 58, 64, 65, 70, 73, 78, 81, 83], 3: [104, 143, 147, 148, 151, 154, 158, 159, 162, 164, 166, 168, 175, 176, 180, 181, 201, 211, 216, 218, 221, 229, 234, 242, 243, 245, 247, 249], 4: [138, 178, 183, 186, 188, 195, 198, 199, 202, 203, 207, 213, 217, 224, 54], 5: [13, 15, 16, 18, 196, 20, 30, 31, 47, 55, 57], 6: [100, 102, 103, 105, 106, 111, 82, 98], 7: [223, 232, 235, 236, 238, 248], 8: [149, 152, 153], 9: [156, 165], 10: [204, 205], 11: [35, 39]}


Plot clusters
-------------


We can plot everything (The 12 clusters found):


.. code:: python

    cl = models.cluster_analysis_dendrogram(color=True)

.. image:: pictures/Tadbit_for_IMP_notebook_24_0.png

Or just 6 of them (without this colors that no one understands...)


.. code:: python

    cl = models.cluster_analysis_dendrogram(n_best_clusters=6)

.. image:: pictures/Tadbit_for_IMP_notebook_26_0.png


Distance between 2 particles
----------------------------


We can just quickly get a value of the distance between particle 13 and 23


.. code:: python

    models.average_3d_dist(13, 23, plot=False)


.. parsed-literal::

    315.29332979218623


This by default, is calculated over the ensemble of models we have. Lets plot the distribution used to get this mean value:


.. code:: python

    models.average_3d_dist(13, 23, plot=True)

.. image:: pictures/Tadbit_for_IMP_notebook_31_0.png

We may also want to use only the 10 first models, or the models belonging to cluster number 0:


.. code:: python

    models.average_3d_dist(13, 23, models=range(10))

.. image:: pictures/Tadbit_for_IMP_notebook_33_0.png



.. code:: python

    models.average_3d_dist(13, 23, plot=True, cluster=0)

.. image:: pictures/Tadbit_for_IMP_notebook_34_0.png


Density plot
------------


Using distances between particle, we can plot now the density (bp per nm) of our chromosomic region.


.. code:: python

    models.density_plot(models=None)

.. image:: pictures/Tadbit_for_IMP_notebook_37_0.png



.. code:: python

    models.density_plot(cluster=0, error=True, steps=(5,20))

.. image:: pictures/Tadbit_for_IMP_notebook_38_0.png


Contact Map
-----------




.. code:: python

    models.contact_map_consistency(models=None, cluster=None, cutoff=150)

.. image:: pictures/Tadbit_for_IMP_notebook_40_0.png


Consistency Plot
----------------




.. code:: python

    models.model_consistency(cluster=0, cutoffs=(50, 100, 150, 200))

.. image:: pictures/Tadbit_for_IMP_notebook_42_0.png



.. code:: python

    

