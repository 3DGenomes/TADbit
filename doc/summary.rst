=======================================
Summary of TADbit classes and functions
=======================================

**Note**:
  - *1: for functions generating plots*
  - *2: for functions generating out files*

Imp imp_modelling module
------------------------

   - **generate_3d_models** (2):             This function generates three-dimensional models starting from Hi-C data.                                             The final analysis will be performed on the n_keep top models.

Parsers hic_parser module
-------------------------

   - **autoreader**:                         Auto-detect matrix format of HiC data file.

   - **read_matrix**:                        Read and checks a matrix from a file (using                                             autoreader) or a list.

Imp structuralmodels module
---------------------------

   - **load_structuralmodels**:              Loads StructuralModels from a file                                             (generated with                                             save_models).

StructuralModels class
++++++++++++++++++++++
    This class contains three-dimensional models generated from a single Hi-C
    data. They can be reached either by their index (integer representing their
    rank according to objective function value), or by their IMP random intial
    number (as string).

      - **angle_between_3_particles**:       Calculates the angle between 3 particles.                                                                                                                                       Given three particles A, B and C, the angle g (angle ACB, shown below):

      - **average_model**:                   Builds and returns an average model representing a given group of models

      - **centroid_model**:                  Estimates and returns the centroid model of a given group of models.

      - **cluster_analysis_dendrogram** (1): Representation of the clustering results. The length of the leaves if                                             proportional to the final objective function value of each model. The                                             branch widths are proportional to the number of models in a given                                             cluster (or group of clusters, if it is an internal branch).

      - **cluster_models**:                  This function performs a clustering analysis of the generated models                                             based on structural comparison. The result will be stored in                                             StructuralModels.clusters                                                                                          Clustering is done according to a score of pairwise comparison                                             calculated as:

      - **contact_map** (1,2):               Plots a contact map representing the frequency of interaction (defined                                             by a distance cutoff) between two particles.

      - **correlate_with_real_data** (1):    Plots the result of a correlation between a given group of models and                                             original Hi-C data.

      - **deconvolve** (1):                  This function performs a clustering analysis of the generated models                                             based on structural comparison (dRMSD).                                             Then, performs a differential contact map between each possible pair                                             of cluster.

      - **define_best_models**:              Defines the number of top models (based on the objective function) to                                             keep. If keep_all is set to True in                                             generate_3d_models or in                                             model_region, then the full set                                             of models (n_models parameter) will be used, otherwise only the n_keep                                             models will be available.

      - **density_plot** (1,2):              Plots the number of nucleotides per nm of chromatin vs the modeled                                             region bins.

      - **dihedral_angle**:                  Calculates the dihedral angle between 2 planes formed by 4 particles.

      - **fetch_model_by_rand_init**:        Models are stored according to their objective function value (first                                             best), but in order to reproduce a model, we need its initial random                                             number. This method helps to fetch the model corresponding to a given                                             initial random number stored under                                             StructuralModels.models[N]['rand_init'].

      - **get_contact_matrix**:              Returns a matrix with the number of interactions observed below a given                                             cutoff distance.

      - **median_3d_dist** (1):              Computes the median distance between two particles over a set of models

      - **model_consistency** (1,2):         Plots the particle consistency, over a given set of models, vs the                                             modeled region bins. The consistency is a measure of the variability                                             (or stability) of the modeled region (the higher the consistency value,                                             the higher stability).

      - **objective_function_model** (1):    This function plots the objective function value per each Monte-Carlo                                             step

      - **particle_coordinates**:            Returns the mean coordinate of a given particle in a group of models.

      - **save_models** (2):                 Saves all the models in pickle format (python object written to disk).

      - **view_models** (1):                 Visualize a selected model in the three dimensions (either with Chimera                                             or through matplotlib).

      - **walking_angle** (1,2):             Plots the angle between successive loci in a given model or set of                                             models. In order to limit the noise of the measure angle is calculated                                             between 3 loci, between each are two other loci. E.g. in the scheme                                             bellow, angle are calculated between loci A, D and G.

      - **walking_dihedral** (1):            Plots the dihedral angle between successive planes. A plane is formed by                                             3 successive loci.

      - **zscore_plot** (1):                 Generate 3 plots. Two heatmaps of the Z-scores used for modeling, one                                             of which is binary showing in red Z-scores higher than upper cut-off;                                             and in blue Z-scores lower than lower cut-off. Last plot is an histogram                                             of the distribution of Z-scores, showing selected regions.

Tad_clustering tad_cmo module
-----------------------------

   - **core_nw**:                            Core of the fast Needleman-Wunsch algorithm that aligns matrices

   - **core_nw_long**:                       Core of the long Needleman-Wunsch algorithm that aligns matrices

   - **optimal_cmo**:                        Calculates the optimal contact map overlap between 2 matrices

Utils three_dim_stats module
----------------------------

   - **square_distance**:                    Calculates the square distance between two particles.

   - **dihedral**:                           Calculates dihedral angle between 4 points in 3D (array with x,y,z)

   - **generate_circle_points**:             Returns list of 3d coordinates of points on a circle using the                                             Rodrigues rotation formula.                                                                                          see *Murray, G. (2013). Rotation About an Arbitrary Axis in 3 Dimensions*                                             for details

   - **calc_eqv_rmsd**:                      Calculates the RMSD, dRMSD, the number of equivalent positions and a score                                             combining these three measures. The measure are done between a group of                                             models in a one against all manner.

   - **fast_square_distance**:               Calculates the square distance between two coordinates.

   - **angle_between_3_points**:             Calculates the angle between 3 particles                                                                                          Given three particles A, B and C, the angle g (angle ACB, shown below):

   - **generate_sphere_points**:             Returns list of 3d coordinates of points on a sphere using the                                             Golden Section Spiral algorithm.

   - **build_mesh**:                         Main function for the calculation of the accessibility of a model.

Utils extraviews module
-----------------------

   - **compare_models**:                     Plots the difference of contact maps of two group of structural models.

   - **plot_3d_model** (1):                  Given a 3 lists of coordinates (x, y, z) plots a three-dimentional model                                             using matplotlib

   - **color_residues**:                     Function to color residues from blue to red.

   - **plot_2d_optimization_result**:        A grid of heatmaps representing the result of the optimization.

   - **colorize**:                           Colorize with ANSII colors a string for printing in shell. this acording to                                             a given number between 0 and 10

   - **tad_border_coloring**:                Colors TAD borders from blue to red (bad to good score). TAD are displayed                                             in scale of grey, from light to dark grey (first to last particle in the                                             TAD)

   - **tad_coloring**:                       Colors TADs from blue to red (first to last TAD). TAD borders are displayed                                             in scale of grey, from light to dark grey (again first to last border)

   - **augmented_dendrogram** (1):           

   - **chimera_view** (1):                   

   - **plot_3d_optimization_result**:        Displays a three dimensional scatter plot representing the result of the                                             optimization.

   - **nicer**:                              writes resolution number for human beings.

Boundary_aligner reciprocally module
------------------------------------

   - **reciprocal**:                         Method based on reciprocal closest boundaries (bd). bd1 will be aligned                                             with bd2 (closest boundary from bd1) if and only if bd1 is the closest                                             boundary of bd2 too (and of course if the distance between bd1 and bd2 is                                             lower than max_dist).

   - **find_closest_reciprocal**:            Function to check the needleman_wunsch algorithm.

Utils tadmaths module
---------------------

   - **zscore**:                             Calculates the log10, Z-score of a given list of values.

   - **calinski_harabasz**:                  Implementation of the CH score [CalinskiHarabasz1974], that has shown to be                                             one the most accurate way to compare clustering methods                                             [MilliganCooper1985] [Tibshirani2001].                                                                                          The CH score is:

Interpolate class
+++++++++++++++++
                      simple linear interpolation

Parsers tad_parser module
-------------------------

   - **parse_tads**:                         Parse a tab separated value file that contains the list of TADs of a given                                             experiment. This file might have been generated whith the                                             print_result_R or with the R binding for tadbit

Imp impoptimizer module
-----------------------

IMPoptimizer class
++++++++++++++++++
    This class optimizes a set of paramaters (scale, maxdist, lowfreq and
    upfreq) in order to maximize the correlation between the models generated
    by IMP and the input data.

      - **get_best_parameters_dict**:        

      - **load_from_file**:                  Loads the optimized parameters from a file generated with the function:                                             pytadbit.imp.impoptimizer.IMPoptimizer.write_result.                                             This function does not overwrite the parameters that were already                                             loaded or calculated.

      - **plot_2d**:                         A grid of heatmaps representing the result of the optimization.

      - **plot_3d**:                         A grid of heatmaps representing the result of the optimization.

      - **run_grid_search**:                 This function calculates the correlation between the models generated                                             by IMP and the input data for the four main IMP parameters (scale,                                             maxdist, lowfreq and upfreq) in the given ranges of values.

      - **write_result**:                    This function writes a log file of all the values tested for each                                             parameter, and the resulting correlation value.                                                                                          This file can be used to load or merge data a posteriori using                                             the function pytadbit.imp.impoptimizer.IMPoptimizer.load_from_file

Tadbit module
-------------

   - **tadbit**:                             The TADbit algorithm works on raw chromosome interaction count data.                                             The normalization is neither necessary nor recommended,                                             since the data is assumed to be discrete counts.                                                                                          TADbit is a breakpoint detection algorithm that returns the optimal                                             segmentation of the chromosome under BIC-penalized likelihood. The                                             model assumes that counts have a Poisson distribution and that the                                             expected value of the counts decreases like a power-law with the                                             linear distance on the chromosome. This expected value of the counts                                             at position (i,j) is corrected by the counts at diagonal positions                                             (i,i) and (j,j). This normalizes for different restriction enzynme                                             site densities and 'mappability' of the reads in case a bin contains                                             repeated regions.

   - **batch_tadbit** (2):                   Use tadbit on directories of data files.                                             All files in the specified directory will be considered data file. The                                             presence of non data files will cause the function to either crash or                                             produce aberrant results.                                                                                          Each file has to contain the data for a single unit/chromosome. The                                             files can be separated in sub-directories corresponding to single                                             experiments or any other organization. Data files that should be                                             considered replicates have to start with the same characters, until                                             the character sep. For instance, all replicates of the unit                                             'chr1' should start with 'chr1\_', using the default value of sep.                                                                                          The data files are read through read.delim. You can pass options                                             to read.delim through the list read_options. For instance                                             if the files have no header, use read_options=list(header=FALSE) and if                                             they also have row names, read_options=list(header=FALSE, row.names=1).                                                                                          Other arguments such as max_size, n_CPU and verbose are passed to                                             tadbit.

Imp impmodel module
-------------------

   - **load_impmodel_from_xyz**:             Loads an IMPmodel object using an xyz file of the form:

   - **load_impmodel_from_cmm**:             Loads an IMPmodel object using an cmm file of the form:

IMPmodel class
++++++++++++++
    A container for the IMP modeling results.

      - **accessible_surface** (1):          Calculates a mesh surface around the model (distance equal to input                                             **radius**) and checks if each point of this mesh could be replaced by                                             an object (i.e. a protein) of a given **radius**                                                                                          Outer part of the model can be excluded from the estimation of                                             accessible surface, as the occupancy outside the model is unkown (see                                             superradius option).

      - **center_of_mass**:                  Gives the center of mass of a model

      - **contour**:                         Total length of the model

      - **cube_side**:                       Calculates the side of a cube containing the model.

      - **cube_volume**:                     Calculates the volume of a cube containing the model.

      - **distance**:                        Calculates the distance between one point of the model and an external                                             coordinate

      - **inaccessible_particles**:          Gives the number of loci/particles that are accessible to an object                                             (i.e. a protein) of a given size.

      - **longest_axe**:                     

      - **min_max_by_axis**:                 Calculates the minimum and maximum coordinates of the model

      - **objective_function** (1):          This function plots the objective function value per each Monte-Carlo                                             step.

      - **radius_of_gyration**:              Calculates the radius of gyration or gyradius of the model                                                                                          Defined as:

      - **shortest_axe**:                    Minimum distance between two particles in the model

      - **view_model** (1):                  Visualize a selected model in the three dimensions. (either with Chimera                                             or through matplotlib).

      - **write_cmm** (2):                   Save a model in the cmm format, read by Chimera                                             (http://www.cgl.ucsf.edu/chimera).                                                                                          **Note:** If none of model_num, models or cluster parameter are set,                                             ALL the models will be written.

      - **write_xyz** (2):                   Writes a xyz file containing the 3D coordinates of each particle in the                                             model.                                                                                          **Note:** If none of model_num, models or cluster parameter are set,                                             ALL the models will be written.

Chromosome module
-----------------

   - **load_chromosome**:                    Load a Chromosome object from a file. A Chromosome object can be saved with                                             the save_chromosome function.

ChromosomeSize class
++++++++++++++++++++
                      This is an integer.
                      
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
                      This is an integer.
                      
                      Relative Chromosome size in base pairs.

Chromosome class
++++++++++++++++
    A Chromosome object designed to deal with Topologically Associating Domains
    predictions from different experiments, in different cell types for a given
    chromosome of DNA, and to compare them.

      - **add_experiment**:                  Add a Hi-C experiment to Chromosome

      - **align_experiments**:               Align the predicted boundaries of two different experiments. The                                             resulting alignment will be stored in the self.experiment list.

      - **find_tad**:                        Call the tadbit function to calculate the                                             position of Topologically Associated Domain boundaries

      - **get_experiment**:                  Fetch an Experiment according to its name.                                             This can also be done directly with Chromosome.experiments[name].

      - **get_tad_hic**:                     Retrieve the Hi-C data matrix corresponding to a given TAD.

      - **iter_tads**:                       Iterate over the TADs corresponding to a given experiment.

      - **save_chromosome**:                 Save a Chromosome object to a file (it uses load from                                             the cPickle). Once saved, the object can be loaded with                                             load_chromosome.

      - **set_max_tad_size**:                Change the maximum size allowed for TADs. It also applies to the                                             computed experiments.

      - **visualize** (1):                   Visualize the matrix of Hi-C interactions of a given experiment

Experiment module
-----------------

Experiment class
++++++++++++++++
    Hi-C experiment.

      - **get_hic_matrix**:                  Return the Hi-C matrix.

      - **get_hic_zscores**:                 Normalize the Hi-C raw data. The result will be stored into                                             the private Experiment._zscore list.

      - **load_hic_data**:                   Add a Hi-C experiment to the Chromosome object.

      - **load_tad_def**:                    Add the Topologically Associated Domains definition detection to Slice

      - **model_region** (2):                Generates of three-dimentional models using IMP, for a given segment of                                             chromosome.

      - **normalize_hic**:                   Normalize the Hi-C data. This normalization step does the same of                                             the tadbit function (default parameters),                                                                                          It fills the Experiment.norm variable with the Hi-C values divided by                                             the calculated weight.                                                                                          The weight of a given cell in column i and row j corresponds to the                                             square root of the product of the sum of column i by the sum of row                                             j.                                                                                          normalization is done according to this formula:

      - **optimal_imp_parameters** (2):      Find the optimal set of parameters to be used for the 3D modeling in                                             IMP.

      - **print_hic_matrix**:                Return the Hi-C matrix as string

      - **set_resolution**:                  Set a new value for the resolution. Copy the original data into                                             Experiment._ori_hic and replace the Experiment.hic_data                                             with the data corresponding to new data                                             (compare_condition).

      - **write_interaction_pairs**:         Creates a tab separated file with all the pairwise interactions.

Boundary_aligner globally module
--------------------------------

   - **equal**:                              

   - **needleman_wunsch**:                   Align two lists of TAD boundaries.

   - **virgin_score**:                       creates empty matrix

Utils hic_filtering module
--------------------------

   - **filter_by_mean**:                     fits the distribution of Hi-C interaction count by column in the matrix to                                             a polynomial. Then searches for the first possible

   - **hic_filtering_for_modelling**:        Main filtering function, to remove artefactual columns in a given Hi-C                                             matrix

   - **filter_by_zero_count**:               fits the distribution of Hi-C interaction count by column in the matrix to                                             a polynomial. Then searches for the first possible

Alignment module
----------------

   - **generate_shuffle_tads**:              Returns a shuffle version of a given list of TADs

   - **randomization_test**:                 Return the probability that original alignment is better than an                                             alignment of randomized boundaries.

   - **generate_rnd_tads**:                  Generates random TADs over a chromosome of a given size according to a given                                             distribution of lengths of TADs.

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

      - **draw** (1):                        Draw alignments as a plot.

      - **get_column**:                      Get a list of column responding to a given characteristic.

      - **itercolumns**:                     Iterate over columns in the alignment

      - **iteritems**:                       Iterate over experiment names and aligned boundaries

      - **itervalues**:                      Iterate over experiment names and aligned boundaries

      - **write_alignment**:                 Print alignment of TAD boundaries between different experiments.                                             Alignment are displayed with colors according to the tadbit                                             confidence score for each boundary.

Boundary_aligner aligner module
-------------------------------

   - **consensusize**:                       Given two alignments returns a consensus alignment. Used for the generation                                             of multiple alignments

   - **align**:                              Align Topologically Associating Domain borders. Supports multiple alignment                                             by building a consensus TAD and aligning each TAD to it.

