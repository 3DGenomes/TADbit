=======================================
Summary of TADbit classes and functions
=======================================


Alignment module
----------------

   - `generate_shuffle_tads <http://3dgenomes.github.io/tadbit/reference/reference_boundary_alignment.html#pytadbit.alignment.generate_shuffle_tads>`_: Returns a shuffle version of a given list of TADs

   - `randomization_test <http://3dgenomes.github.io/tadbit/reference/reference_boundary_alignment.html#pytadbit.alignment.randomization_test>`_: Return the probability that original alignment is better than an                                             alignment of randomized boundaries.

   - `generate_rnd_tads <http://3dgenomes.github.io/tadbit/reference/reference_boundary_alignment.html#pytadbit.alignment.generate_rnd_tads>`_: Generates random TADs over a chromosome of a given size according to a given                                             distribution of lengths of TADs.

TAD class
+++++++++
                      Specific class of TADs, used only within Alignment objects.
                      It is directly inheriting from python dict.
                      a TAD these keys:
                      
                      - 'start': position of the TAD
                      - 'end': position of the TAD
                      - 'score': of the prediction of boundary
                      - 'brk': same as 'end'
                      - 'pos': in the alignment (column number)
                      - 'exp': Experiment this TAD belongs to
                      - 'index': of this TAD within all TADs in the Experiment

Alignment class
+++++++++++++++
    Alignment of TAD borders

      - `draw <http://3dgenomes.github.io/tadbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.draw>`_ [#first]_: Draw alignments as a plot.

      - `get_column <http://3dgenomes.github.io/tadbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.get_column>`_: Get a list of column responding to a given characteristic.

      - `itercolumns <http://3dgenomes.github.io/tadbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.itercolumns>`_: Iterate over columns in the alignment

      - `iteritems <http://3dgenomes.github.io/tadbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.iteritems>`_: Iterate over experiment names and aligned boundaries

      - `itervalues <http://3dgenomes.github.io/tadbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.itervalues>`_: Iterate over experiment names and aligned boundaries

      - `write_alignment <http://3dgenomes.github.io/tadbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.write_alignment>`_: Print alignment of TAD boundaries between different experiments.                                             Alignment are displayed with colors according to the tadbit                                             confidence score for each boundary.

Boundary_aligner aligner module
-------------------------------

   - consensusize:                           Given two alignments returns a consensus alignment. Used for the generation                                             of multiple alignments

   - `align <http://3dgenomes.github.io/tadbit/reference/reference_aligner.html#pytadbit.boundary_aligner.aligner.align>`_: Align Topologically Associating Domain borders. Supports multiple alignment                                             by building a consensus TAD sequence and aligning each experiment to it.

Boundary_aligner globally module
--------------------------------

   - `needleman_wunsch <http://3dgenomes.github.io/tadbit/reference/reference_aligner.html#pytadbit.boundary_aligner.globally.needleman_wunsch>`_: Align two lists of TAD boundaries using a Needleman-Wunsh implementation

Boundary_aligner reciprocally module
------------------------------------

   - find_closest_reciprocal:                Function to check the needleman_wunsch algorithm.

   - `reciprocal <http://3dgenomes.github.io/tadbit/reference/reference_aligner.html#pytadbit.boundary_aligner.reciprocally.reciprocal>`_: Method based on reciprocal closest boundaries (bd). bd1 will be aligned                                             with bd2 (closest boundary from bd1) if and only if bd1 is the closest                                             boundary of bd2 too (and of course if the distance between bd1 and bd2 is                                             lower than max_dist).

Chromosome module
-----------------

   - `load_chromosome <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.load_chromosome>`_: Load a Chromosome object from a file. A Chromosome object can be saved with                                             the save_chromosome function.

ChromosomeSize class
++++++++++++++++++++
                      Chromosome size in base pairs

ExperimentList class
++++++++++++++++++++
                      Inherited from python built in list, modified for tadbit
                      Experiment.
                      
                      Mainly, `getitem`, `setitem`, and `append` were modified in order to
                      be able to search for experiments by index or by name, and to add
                      experiments simply using Chromosome.experiments.append(Experiment).
                      
                      The whole ExperimentList object is linked to a Chromosome instance
                      (Chromosome).

AlignmentDict class
+++++++++++++++++++
                      dict of Alignment
                      
                      Modified getitem, setitem, and append in order to be able to search
                      alignments by index or by name.
                      
                      linked to a Chromosome

RelativeChromosomeSize class
++++++++++++++++++++++++++++
                      Relative Chromosome size in base pairs.
                      
                      Only used for TAD alignment randomization.

Chromosome class
++++++++++++++++
    A Chromosome object designed to deal with Topologically Associating Domains
    predictions from different experiments, in different cell types for a given
    chromosome of DNA, and to compare them.

      - `add_experiment <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.add_experiment>`_: Add a Hi-C experiment to Chromosome

      - `align_experiments <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.align_experiments>`_: Align the predicted boundaries of two different experiments. The                                             resulting alignment will be stored in the self.experiment list.

      - `find_tad <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.find_tad>`_: Call the tadbit function to calculate the                                             position of Topologically Associated Domain boundaries

      - `get_experiment <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.get_experiment>`_: Fetch an Experiment according to its name.                                             This can also be done directly with Chromosome.experiments[name].

      - `get_tad_hic <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.get_tad_hic>`_: Retrieve the Hi-C data matrix corresponding to a given TAD.

      - `iter_tads <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.iter_tads>`_: Iterate over the TADs corresponding to a given experiment.

      - `save_chromosome <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.save_chromosome>`_: Save a Chromosome object to a file (it uses load from                                             the cPickle). Once saved, the object can be loaded with                                             load_chromosome.

      - `set_max_tad_size <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.set_max_tad_size>`_: Change the maximum size allowed for TADs. It also applies to the                                             computed experiments.

      - `tad_density_plot <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.tad_density_plot>`_ [#first]_: Draw an summary of the TAD found in a given experiment and their density                                             in terms of relative Hi-C interaction count.

      - `visualize <http://3dgenomes.github.io/tadbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.visualize>`_ [#first]_: Visualize the matrix of Hi-C interactions of a given experiment

Experiment module
-----------------

Experiment class
++++++++++++++++
    Hi-C experiment.

      - `get_hic_matrix <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.get_hic_matrix>`_: Return the Hi-C matrix.

      - `get_hic_zscores <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.get_hic_zscores>`_: Normalize the Hi-C raw data. The result will be stored into                                             the private Experiment._zscore list.

      - `load_hic_data <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.load_hic_data>`_: Add a Hi-C experiment to the Chromosome object.

      - `load_tad_def <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.load_tad_def>`_: Add the Topologically Associated Domains definition detection to Slice

      - `model_region <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.model_region>`_ [#second]_: Generates of three-dimentional models using IMP, for a given segment of                                             chromosome.

      - `normalize_hic <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.normalize_hic>`_: Normalize the Hi-C data. This normalization step does the same of                                             the tadbit function (default parameters),                                                                                          It fills the Experiment.norm variable with the Hi-C values divided by                                             the calculated weight.                                                                                          The weight of a given cell in column i and row j corresponds to the                                             square root of the product of the sum of column i by the sum of row                                             j.                                                                                          normalization is done according to this formula:

      - `optimal_imp_parameters <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.optimal_imp_parameters>`_ [#second]_: Find the optimal set of parameters to be used for the 3D modeling in                                             IMP.

      - `print_hic_matrix <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.print_hic_matrix>`_: Return the Hi-C matrix as string

      - `set_resolution <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.set_resolution>`_: Set a new value for the resolution. Copy the original data into                                             Experiment._ori_hic and replace the Experiment.hic_data                                             with the data corresponding to new data                                             (compare_condition).

      - `view <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.view>`_ [#first]_: Visualize the matrix of Hi-C interactions

      - `write_interaction_pairs <http://3dgenomes.github.io/tadbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.write_interaction_pairs>`_: Creates a tab separated file with all the pairwise interactions.

Imp imp_modelling module
------------------------

   - `generate_3d_models <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.imp_modelling.generate_3d_models>`_ [#second]_: This function generates three-dimensional models starting from Hi-C data.                                             The final analysis will be performed on the n_keep top models.

Imp impmodel module
-------------------

   - `load_impmodel_from_xyz <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.load_impmodel_from_xyz>`_: Loads an IMPmodel object using an xyz file of the form:

   - `load_impmodel_from_cmm <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.load_impmodel_from_cmm>`_: Loads an IMPmodel object using an cmm file of the form:

IMPmodel class
++++++++++++++
    A container for the IMP modeling results.

      - `accessible_surface <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.accessible_surface>`_ [#first]_: Calculates a mesh surface around the model (distance equal to input                                             **radius**) and checks if each point of this mesh could be replaced by                                             an object (i.e. a protein) of a given **radius**                                                                                          Outer part of the model can be excluded from the estimation of                                             accessible surface, as the occupancy outside the model is unkown (see                                             superradius option).

      - `center_of_mass <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.center_of_mass>`_: Gives the center of mass of a model

      - `contour <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.contour>`_: Total length of the model

      - `cube_side <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.cube_side>`_: Calculates the side of a cube containing the model.

      - `cube_volume <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.cube_volume>`_: Calculates the volume of a cube containing the model.

      - `distance <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.distance>`_: Calculates the distance between one point of the model and an external                                             coordinate

      - `inaccessible_particles <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.inaccessible_particles>`_: Gives the number of loci/particles that are accessible to an object                                             (i.e. a protein) of a given size.

      - `longest_axe <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.longest_axe>`_: Gives the distance between most distant particles of the model

      - `min_max_by_axis <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.min_max_by_axis>`_: Calculates the minimum and maximum coordinates of the model

      - `objective_function <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.objective_function>`_ [#first]_: This function plots the objective function value per each Monte-Carlo                                             step.

      - `persistence_length <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.persistence_length>`_: Calculates the persistence length (Lp) of given section of the model.                                             Persistence length is calculated according to [Bystricky2004] :

      - `radius_of_gyration <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.radius_of_gyration>`_: Calculates the radius of gyration or gyradius of the model                                                                                          Defined as:

      - `shortest_axe <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.shortest_axe>`_: Minimum distance between two particles in the model

      - `view_model <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.view_model>`_ [#first]_: Visualize a selected model in the three dimensions. (either with Chimera                                             or through matplotlib).

      - `write_cmm <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.write_cmm>`_ [#second]_: Save a model in the cmm format, read by Chimera                                             (http://www.cgl.ucsf.edu/chimera).                                                                                          **Note:** If none of model_num, models or cluster parameter are set,                                             ALL the models will be written.

      - `write_xyz <http://3dgenomes.github.io/tadbit/reference/reference_imp_model.html#pytadbit.imp.impmodel.IMPmodel.write_xyz>`_ [#second]_: Writes a xyz file containing the 3D coordinates of each particle in the                                             model.                                                                                          **Note:** If none of model_num, models or cluster parameter are set,                                             ALL the models will be written.

Imp impoptimizer module
-----------------------

IMPoptimizer class
++++++++++++++++++
    This class optimizes a set of paramaters (scale, maxdist, lowfreq and
    upfreq) in order to maximize the correlation between the models generated
    by IMP and the input data.

      - `get_best_parameters_dict <http://3dgenomes.github.io/tadbit/reference/reference_imp_optimizer.html#pytadbit.imp.impoptimizer.IMPoptimizer.get_best_parameters_dict>`_: 

      - `load_from_file <http://3dgenomes.github.io/tadbit/reference/reference_imp_optimizer.html#pytadbit.imp.impoptimizer.IMPoptimizer.load_from_file>`_: Loads the optimized parameters from a file generated with the function:                                             pytadbit.imp.impoptimizer.IMPoptimizer.write_result.                                             This function does not overwrite the parameters that were already                                             loaded or calculated.

      - `load_grid_search <http://3dgenomes.github.io/tadbit/reference/reference_imp_optimizer.html#pytadbit.imp.impoptimizer.IMPoptimizer.load_grid_search>`_: Loads one file or a list of files containing pre-calculated Structural                                             Models (keep_models parameter used). And correlate each set of models                                             with real data. Usefull to run different correlation on the same data                                             avoiding to re-calculate each time the models.

      - `plot_2d <http://3dgenomes.github.io/tadbit/reference/reference_imp_optimizer.html#pytadbit.imp.impoptimizer.IMPoptimizer.plot_2d>`_ [#first]_: A grid of heatmaps representing the result of the optimization.

      - `plot_3d <http://3dgenomes.github.io/tadbit/reference/reference_imp_optimizer.html#pytadbit.imp.impoptimizer.IMPoptimizer.plot_3d>`_: A grid of heatmaps representing the result of the optimization.

      - `run_grid_search <http://3dgenomes.github.io/tadbit/reference/reference_imp_optimizer.html#pytadbit.imp.impoptimizer.IMPoptimizer.run_grid_search>`_ [#second]_: This function calculates the correlation between the models generated                                             by IMP and the input data for the four main IMP parameters (scale,                                             maxdist, lowfreq and upfreq) in the given ranges of values.

      - `write_result <http://3dgenomes.github.io/tadbit/reference/reference_imp_optimizer.html#pytadbit.imp.impoptimizer.IMPoptimizer.write_result>`_: This function writes a log file of all the values tested for each                                             parameter, and the resulting correlation value.                                                                                          This file can be used to load or merge data a posteriori using                                             the function pytadbit.imp.impoptimizer.IMPoptimizer.load_from_file

Imp structuralmodels module
---------------------------

   - `load_structuralmodels <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.load_structuralmodels>`_: Loads StructuralModels from a file                                             (generated with                                             save_models).

StructuralModels class
++++++++++++++++++++++
    This class contains three-dimensional models generated from a single Hi-C
    data. They can be reached either by their index (integer representing their
    rank according to objective function value), or by their IMP random intial
    number (as string).

      - `align_models <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.align_models>`_: Three-dimensional aligner for structural models.

      - `angle_between_3_particles <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.angle_between_3_particles>`_: Calculates the angle between 3 particles.                                                                                                                                       Given three particles A, B and C, the angle g (angle ACB, shown below):

      - `average_model <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.average_model>`_: Builds and returns an average model representing a given group of models

      - `centroid_model <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.centroid_model>`_: Estimates and returns the centroid model of a given group of models.

      - `cluster_analysis_dendrogram <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.cluster_analysis_dendrogram>`_ [#first]_: Representation of the clustering results. The length of the leaves if                                             proportional to the final objective function value of each model. The                                             branch widths are proportional to the number of models in a given                                             cluster (or group of clusters, if it is an internal branch).

      - `cluster_models <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.cluster_models>`_: This function performs a clustering analysis of the generated models                                             based on structural comparison. The result will be stored in                                             StructuralModels.clusters                                                                                          Clustering is done according to a score of pairwise comparison                                             calculated as:

      - `contact_map <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.contact_map>`_ [#first]_ [#second]_: Plots a contact map representing the frequency of interaction (defined                                             by a distance cutoff) between two particles.

      - `correlate_with_real_data <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.correlate_with_real_data>`_ [#first]_: Plots the result of a correlation between a given group of models and                                             original Hi-C data.

      - `deconvolve <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.deconvolve>`_ [#first]_: This function performs a deconvolution analysis of a given froup of models.                                             It first clusters models based on structural comparison (dRMSD), and                                             then, performs a differential contact map between each possible pair                                             of cluster.

      - `define_best_models <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.define_best_models>`_: Defines the number of top models (based on the objective function) to                                             keep. If keep_all is set to True in                                             generate_3d_models or in                                             model_region, then the full set                                             of models (n_models parameter) will be used, otherwise only the n_keep                                             models will be available.

      - `density_plot <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.density_plot>`_ [#first]_ [#second]_: Plots the number of nucleotides per nm of chromatin vs the modeled                                             region bins.

      - `dihedral_angle <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.dihedral_angle>`_: Calculates the dihedral angle between 2 planes formed by 4 particles.

      - `fetch_model_by_rand_init <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.fetch_model_by_rand_init>`_: Models are stored according to their objective function value (first                                             best), but in order to reproduce a model, we need its initial random                                             number. This method helps to fetch the model corresponding to a given                                             initial random number stored under                                             StructuralModels.models[N]['rand_init'].

      - `get_contact_matrix <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.get_contact_matrix>`_: Returns a matrix with the number of interactions observed below a given                                             cutoff distance.

      - `interactions <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.interactions>`_ [#first]_ [#second]_: Plots, for each particle, the number of interactions (particles closer                                             than the guiven cut-off). The value given is the average for all models.

      - `median_3d_dist <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.median_3d_dist>`_ [#first]_: Computes the median distance between two particles over a set of models

      - `model_consistency <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.model_consistency>`_ [#first]_ [#second]_: Plots the particle consistency, over a given set of models, vs the                                             modeled region bins. The consistency is a measure of the variability                                             (or stability) of the modeled region (the higher the consistency value,                                             the higher stability).

      - `objective_function_model <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.objective_function_model>`_ [#first]_: This function plots the objective function value per each Monte-Carlo                                             step

      - `particle_coordinates <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.particle_coordinates>`_: Returns the mean coordinate of a given particle in a group of models.

      - `save_models <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.save_models>`_ [#second]_: Saves all the models in pickle format (python object written to disk).

      - `view_centroid <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.view_centroid>`_: shortcut for                                             view_models(tool='plot', show='highlighted', highlight='centroid')

      - `view_models <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.view_models>`_ [#first]_: Visualize a selected model in the three dimensions (either with Chimera                                             or through matplotlib).

      - `walking_angle <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.walking_angle>`_ [#first]_ [#second]_: Plots the angle between successive loci in a given model or set of                                             models. In order to limit the noise of the measure angle is calculated                                             between 3 loci, between each are two other loci. E.g. in the scheme                                             bellow, angle are calculated between loci A, D and G.

      - `walking_dihedral <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.walking_dihedral>`_ [#first]_: Plots the dihedral angle between successive planes. A plane is formed by                                             3 successive loci.

      - `zscore_plot <http://3dgenomes.github.io/tadbit/reference/reference_imp_structuralmodels.html#pytadbit.imp.structuralmodels.StructuralModels.zscore_plot>`_ [#first]_: Generate 3 plots. Two heatmaps of the Z-scores used for modeling, one                                             of which is binary showing in red Z-scores higher than upper cut-off;                                             and in blue Z-scores lower than lower cut-off. Last plot is an histogram                                             of the distribution of Z-scores, showing selected regions. Histogram                                             also shows the fit to normal distribution and the result of a test of                                             normality (D'Agostino and Pearson's test, that combines skew and                                             kurtosis).

Parsers hic_parser module
-------------------------

   - autoreader:                             Auto-detect matrix format of HiC data file.

   - `read_matrix <http://3dgenomes.github.io/tadbit/reference/reference_parser.html#pytadbit.parsers.hic_parser.read_matrix>`_: Read and checks a matrix from a file (using                                             autoreader) or a list.

Parsers tad_parser module
-------------------------

   - `parse_tads <http://3dgenomes.github.io/tadbit/reference/reference_parser.html#pytadbit.parsers.tad_parser.parse_tads>`_: Parse a tab separated value file that contains the list of TADs of a given                                             experiment. This file might have been generated whith the                                             print_result_R or with the R binding for tadbit

Tad_clustering tad_cmo module
-----------------------------

   - core_nw:                                Core of the fast Needleman-Wunsch algorithm that aligns matrices

   - virgin_score:                           Fill a matrix with zeros, except first row and first column filled with     multiple values of penalty.

   - core_nw_long:                           Core of the long Needleman-Wunsch algorithm that aligns matrices

   - `optimal_cmo <http://3dgenomes.github.io/tadbit/reference/reference_clustering.html#pytadbit.tad_clustering.tad_cmo.optimal_cmo>`_: Calculates the optimal contact map overlap between 2 matrices

Tadbit module
-------------

   - `tadbit <http://3dgenomes.github.io/tadbit/reference/reference_tadbit.html#pytadbit.tadbit.tadbit>`_: The TADbit algorithm works on raw chromosome interaction count data.                                             The normalization is neither necessary nor recommended,                                             since the data is assumed to be discrete counts.                                                                                          TADbit is a breakpoint detection algorithm that returns the optimal                                             segmentation of the chromosome under BIC-penalized likelihood. The                                             model assumes that counts have a Poisson distribution and that the                                             expected value of the counts decreases like a power-law with the                                             linear distance on the chromosome. This expected value of the counts                                             at position (i,j) is corrected by the counts at diagonal positions                                             (i,i) and (j,j). This normalizes for different restriction enzynme                                             site densities and 'mappability' of the reads in case a bin contains                                             repeated regions.

   - `batch_tadbit <http://3dgenomes.github.io/tadbit/reference/reference_tadbit.html#pytadbit.tadbit.batch_tadbit>`_ [#second]_: Use tadbit on directories of data files.                                             All files in the specified directory will be considered data file. The                                             presence of non data files will cause the function to either crash or                                             produce aberrant results.                                                                                          Each file has to contain the data for a single unit/chromosome. The                                             files can be separated in sub-directories corresponding to single                                             experiments or any other organization. Data files that should be                                             considered replicates have to start with the same characters, until                                             the character sep. For instance, all replicates of the unit                                             'chr1' should start with 'chr1\_', using the default value of sep.                                                                                          The data files are read through read.delim. You can pass options                                             to read.delim through the list read_options. For instance                                             if the files have no header, use read_options=list(header=FALSE) and if                                             they also have row names, read_options=list(header=FALSE, row.names=1).                                                                                          Other arguments such as max_size, n_CPU and verbose are passed to                                             tadbit.

Utils extraviews module
-----------------------

   - `compare_models <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.extraviews.compare_models>`_: Plots the difference of contact maps of two group of structural models.

   - `plot_3d_model <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.extraviews.plot_3d_model>`_ [#first]_: Given a 3 lists of coordinates (x, y, z) plots a three-dimentional model                                             using matplotlib

   - `color_residues <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.extraviews.color_residues>`_: Function to color residues from blue to red.

   - `plot_2d_optimization_result <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.extraviews.plot_2d_optimization_result>`_ [#first]_: A grid of heatmaps representing the result of the optimization.

   - colorize:                               Colorize with ANSII colors a string for printing in shell. this acording to                                             a given number between 0 and 10

   - `tad_border_coloring <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.extraviews.tad_border_coloring>`_: Colors TAD borders from blue to red (bad to good score). TAD are displayed                                             in scale of grey, from light to dark grey (first to last particle in the                                             TAD)

   - `tad_coloring <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.extraviews.tad_coloring>`_: Colors TADs from blue to red (first to last TAD). TAD borders are displayed                                             in scale of grey, from light to dark grey (again first to last border)

   - augmented_dendrogram [#first]_:         

   - `chimera_view <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.extraviews.chimera_view>`_ [#first]_: Open a list of .cmm files with Chimera (http://www.cgl.ucsf.edu/chimera)                                             to view models.

   - `plot_3d_optimization_result <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.extraviews.plot_3d_optimization_result>`_: Displays a three dimensional scatter plot representing the result of the                                             optimization.

   - nicer:                                  writes resolution number for human beings.

Utils hic_filtering module
--------------------------

   - filter_by_mean:                         fits the distribution of Hi-C interaction count by column in the matrix to                                             a polynomial. Then searches for the first possible

   - `hic_filtering_for_modelling <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.hic_filtering.hic_filtering_for_modelling>`_: Main filtering function, to remove artefactual columns in a given Hi-C                                             matrix

   - filter_by_zero_count:                   fits the distribution of Hi-C interaction count by column in the matrix to                                             a polynomial. Then searches for the first possible

Utils tadmaths module
---------------------

   - zscore:                                 Calculates the log10, Z-score of a given list of values.

   - `calinski_harabasz <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.tadmaths.calinski_harabasz>`_: Implementation of the CH score [CalinskiHarabasz1974], that has shown to be                                             one the most accurate way to compare clustering methods                                             [MilliganCooper1985] [Tibshirani2001].                                                                                          The CH score is:

   - newton_raphson:                         Newton-Raphson method as defined in:                                             http://www.maths.tcd.ie/~ryan/TeachingArchive/161/teaching/newton-raphson.c.html                                             used to search for the persistence length of a given model.

Interpolate class
+++++++++++++++++
                      simple linear interpolation

Utils three_dim_stats module
----------------------------

   - `square_distance <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.square_distance>`_: Calculates the square distance between two particles.

   - `dihedral <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.dihedral>`_: Calculates dihedral angle between 4 points in 3D (array with x,y,z)

   - `generate_circle_points <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.generate_circle_points>`_: Returns list of 3d coordinates of points on a circle using the                                             Rodrigues rotation formula.                                                                                          see *Murray, G. (2013). Rotation About an Arbitrary Axis in 3 Dimensions*                                             for details

   - mass_center:                            Transforms coordinates according to the center of mass

   - rotate_among_y_axis:                    Rotate and object with a list of x, y, z coordinates among its center of                                             mass

   - `calc_eqv_rmsd <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.calc_eqv_rmsd>`_: Calculates the RMSD, dRMSD, the number of equivalent positions and a score                                             combining these three measures. The measure are done between a group of                                             models in a one against all manner.

   - get_center_of_mass:                     get the center of mass of a given object with list of x, y, z coordinates

   - find_angle_rotation_improve_x:          Finds the rotation angle needed to face the longest edge of the molecule

   - fast_square_distance:                   Calculates the square distance between two coordinates.

   - `angle_between_3_points <http://3dgenomes.github.io/tadbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.angle_between_3_points>`_: Calculates the angle between 3 particles                                                                                          Given three particles A, B and C, the angle g (angle ACB, shown below):

   - generate_sphere_points:                 Returns list of 3d coordinates of points on a sphere using the                                             Golden Section Spiral algorithm.

   - build_mesh:                             Main function for the calculation of the accessibility of a model.


.. [#first] functions generating plots

.. [#second] functions writing text files

