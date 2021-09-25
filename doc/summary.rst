=======================================
Summary of TADbit classes and functions
=======================================


Root module
-----------

   - `get_dependencies_version <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.get_dependencies_version>`_: Check versions of TADbit and all dependencies, as well and retrieves system                                             info. May be used to ensure reproducibility.

Alignment module
----------------

   - `generate_rnd_tads <http://3dgenomes.github.io/TADbit/reference/reference_boundary_alignment.html#pytadbit.alignment.generate_rnd_tads>`_: Generates random TADs over a chromosome of a given size according to a given                                             distribution of lengths of TADs.

   - `generate_shuffle_tads <http://3dgenomes.github.io/TADbit/reference/reference_boundary_alignment.html#pytadbit.alignment.generate_shuffle_tads>`_: Returns a shuffle version of a given list of TADs

   - `randomization_test <http://3dgenomes.github.io/TADbit/reference/reference_boundary_alignment.html#pytadbit.alignment.randomization_test>`_: Return the probability that original alignment is better than an                                             alignment of randomized boundaries.

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

      - `draw <http://3dgenomes.github.io/TADbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.draw>`_ [#first]_: Draw alignments as a plot.

      - `get_column <http://3dgenomes.github.io/TADbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.get_column>`_: Get a list of column responding to a given characteristic.

      - `itercolumns <http://3dgenomes.github.io/TADbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.itercolumns>`_: Iterate over columns in the alignment

      - `iteritems <http://3dgenomes.github.io/TADbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.iteritems>`_: Iterate over experiment names and aligned boundaries

      - `itervalues <http://3dgenomes.github.io/TADbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.itervalues>`_: Iterate over experiment names and aligned boundaries

      - `write_alignment <http://3dgenomes.github.io/TADbit/reference/reference_boundary_alignment.html#pytadbit.alignment.Alignment.write_alignment>`_: Print alignment of TAD boundaries between different experiments.                                             Alignments are displayed with colors according to the TADbit                                             confidence score for each boundary.

Boundary_aligner aligner module
-------------------------------

   - `align <http://3dgenomes.github.io/TADbit/reference/reference_aligner.html#pytadbit.boundary_aligner.aligner.align>`_: Align Topologically Associating Domain borders. Supports multiple-alignment                                             by building a consensus TAD sequence and aligning each experiment to it.

   - consensusize:                           Given two alignments returns a consensus alignment. Used for the generation                                             of multiple alignments

Boundary_aligner globally module
--------------------------------

   - `needleman_wunsch <http://3dgenomes.github.io/TADbit/reference/reference_aligner.html#pytadbit.boundary_aligner.globally.needleman_wunsch>`_: Align two lists of TAD boundaries using a Needleman-Wunsh implementation

Boundary_aligner reciprocally module
------------------------------------

   - `reciprocal <http://3dgenomes.github.io/TADbit/reference/reference_aligner.html#pytadbit.boundary_aligner.reciprocally.reciprocal>`_: Method based on reciprocal closest boundaries (bd). bd1 will be aligned                                             with bd2 (closest boundary from bd1) if and only if bd1 is the closest                                             boundary of bd2 too (and of course if the distance between bd1 and bd2 is                                             lower than max_dist).

   - find_closest_reciprocal:                Function to check the needleman_wunsch algorithm.

Chromosome module
-----------------

   - `load_chromosome <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.load_chromosome>`_: Load a Chromosome object from a file. A Chromosome object can be saved with                                             the save_chromosome function.

AlignmentDict class
+++++++++++++++++++
                      dict of Alignment
                      
                      Modified getitem, setitem, and append in order to be able to search
                      alignments by index or by name.
                      
                      linked to a Chromosome

ExperimentList class
++++++++++++++++++++
                      Inherited from python built in list, modified for TADbit
                      Experiment.
                      
                      Mainly, `getitem`, `setitem`, and `append` were modified in order to
                      be able to search for experiments by index or by name, and to add
                      experiments simply using Chromosome.experiments.append(Experiment).
                      
                      The whole ExperimentList object is linked to a Chromosome instance
                      (Chromosome).

Chromosome class
++++++++++++++++
    A Chromosome object designed to deal with Topologically Associating Domains
    predictions from different experiments, in different cell types for a given
    chromosome of DNA, and to compare them.

      - `add_experiment <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.add_experiment>`_: Add a Hi-C experiment to Chromosome

      - `align_experiments <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.align_experiments>`_: Align the predicted boundaries of two different experiments. The                                             resulting alignment will be stored in the self.experiment list.

      - `find_tad <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.find_tad>`_: Call the tadbit function to calculate the                                             position of Topologically Associated Domain boundaries

      - `get_experiment <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.get_experiment>`_: Fetch an Experiment according to its name.                                             This can also be done directly with Chromosome.experiments[name].

      - `get_tad_hic <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.get_tad_hic>`_: Retrieve the Hi-C data matrix corresponding to a given TAD.

      - `iter_tads <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.iter_tads>`_: Iterate over the TADs corresponding to a given experiment.

      - `save_chromosome <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.save_chromosome>`_: Save a Chromosome object to a file (it uses load from                                             the pickle). Once saved, the object can be loaded with                                             load_chromosome.

      - `set_max_tad_size <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.set_max_tad_size>`_: Change the maximum size allowed for TADs. It also applies to the                                             computed experiments.

      - `tad_density_plot <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.tad_density_plot>`_ [#first]_: Draw an summary of the TAD found in a given experiment and their density                                             in terms of relative Hi-C interaction count.

      - `visualize <http://3dgenomes.github.io/TADbit/reference/reference_chromosome.html#pytadbit.chromosome.Chromosome.visualize>`_ [#first]_: Visualize the matrix of Hi-C interactions of a given experiment

Experiment module
-----------------

   - `load_experiment_from_reads <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.load_experiment_from_reads>`_: Loads an experiment object from TADbit-generated read files, that are lists                                             of pairs of reads mapped to a reference genome.

Experiment class
++++++++++++++++
    Hi-C experiment.

      - `filter_columns <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.filter_columns>`_ [#first]_: Call filtering function, to remove artifactual columns in a given Hi-C                                             matrix. This function will detect columns with very low interaction                                             counts. Filtered out columns will be stored in the dictionary Experiment._zeros.

      - `get_hic_matrix <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.get_hic_matrix>`_: Return the Hi-C matrix.

      - `get_hic_zscores <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.get_hic_zscores>`_: Normalize the Hi-C raw data. The result will be stored into                                             the private Experiment._zscore list.

      - `load_hic_data <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.load_hic_data>`_: Add a Hi-C experiment to the Chromosome object.

      - `load_norm_data <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.load_norm_data>`_: Add a normalized Hi-C experiment to the Chromosome object.

      - `load_tad_def <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.load_tad_def>`_: Add the Topologically Associated Domains definition detection to Slice

      - `model_region <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.model_region>`_ [#second]_: Generates of three-dimensional models using IMP, for a given segment of                                             chromosome.

      - `normalize_hic <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.normalize_hic>`_: Normalize the Hi-C data. This normalization step does the same of                                             the tadbit function (default parameters),                                                                                          It fills the Experiment.norm variable with the Hi-C values divided by                                             the calculated weight.                                                                                          The weight of a given cell in column i and row j corresponds to the                                             square root of the product of the sum of column i by the sum of row                                             j.                                                                                          normalization is done according to this formula:

      - `optimal_imp_parameters <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.optimal_imp_parameters>`_ [#second]_: Find the optimal set of parameters to be used for the 3D modeling in                                             IMP.

      - `print_hic_matrix <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.print_hic_matrix>`_: Return the Hi-C matrix as string

      - `set_resolution <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.set_resolution>`_: Set a new value for the resolution. Copy the original data into                                             Experiment._ori_hic and replace the Experiment.hic_data                                             with the data corresponding to new data                                             (compare_condition).

      - `view <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.view>`_ [#first]_: Visualize the matrix of Hi-C interactions

      - `write_interaction_pairs <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.write_interaction_pairs>`_: Creates a tab separated file with all the pairwise interactions.

      - `write_json <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.write_json>`_: Save hic matrix in the json format, read by TADkit.

      - `write_tad_borders <http://3dgenomes.github.io/TADbit/reference/reference_experiment.html#pytadbit.experiment.Experiment.write_tad_borders>`_ [#second]_: Print a table summarizing the TADs found by tadbit. This function outputs                                             something similar to the R function.

Hic_data module
---------------

   - isclose:                                https://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python/33024979#33024979

HiC_data class
++++++++++++++
    This may also hold the print/write-to-file matrix functions

      - `add_sections <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.add_sections>`_: Add genomic coordinate to HiC_data object by getting them from a FASTA                                             file containing chromosome sequences. Orders matters.

      - `add_sections_from_fasta <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.add_sections_from_fasta>`_: Add genomic coordinate to HiC_data object by getting them from a FASTA                                             file containing chromosome sequences

      - `cis_trans_ratio <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.cis_trans_ratio>`_: Counts the number of interactions occurring within chromosomes (cis) with                                             respect to the total number of interactions

      - `find_compartments <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.find_compartments>`_ [#first]_ [#second]_: Search for A/B compartments in each chromosome of the Hi-C matrix.                                             Hi-C matrix is normalized by the number interaction expected at a given                                             distance, and by visibility (one iteration of ICE). A correlation matrix                                             is then calculated from this normalized matrix, and its first                                             eigenvector is used to identify compartments. Changes in sign marking                                             boundaries between compartments.                                             Result is stored as a dictionary of compartment boundaries, keys being                                             chromosome names.

      - `find_compartments_beta <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.find_compartments_beta>`_ [#first]_ [#second]_: Search for A/B compartments in each chromosome of the Hi-C matrix.                                             Hi-C matrix is normalized by the number interaction expected at a given                                             distance, and by visibility (one iteration of ICE). A correlation matrix                                             is then calculated from this normalized matrix, and its first                                             eigenvector is used to identify compartments. Changes in sign marking                                             boundaries between compartments.                                             Result is stored as a dictionary of compartment boundaries, keys being                                             chromosome names.

      - `get_hic_data_as_csr <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.get_hic_data_as_csr>`_: Returns a scipy sparse matrix in Compressed Sparse Row format of the Hi-C data in the dictionary

      - `get_matrix <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.get_matrix>`_: returns a matrix.

      - `load_biases <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.load_biases>`_: Load biases, decay and bad columns from pickle file

      - `save_biases <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.save_biases>`_: Save biases, decay and bad columns in pickle format (to be loaded by                                             the function load_hic_data_from_bam)

      - `sum <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.sum>`_: Sum Hi-C data matrix                                             WARNING: parameters are not meant to be used by external users

      - `write_compartments <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.write_compartments>`_ [#second]_: Write compartments to a file.

      - `write_cooler <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.write_cooler>`_: writes the hic_data to a cooler file.

      - `write_coord_table <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.write_coord_table>`_: writes a coordinate table to a file.

      - `write_matrix <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.write_matrix>`_: writes the matrix to a file.

      - `yield_matrix <http://3dgenomes.github.io/TADbit/reference/reference_hic_data.html#pytadbit.hic_data.HiC_data.yield_matrix>`_: Yields a matrix line by line.                                             Bad row/columns are returned as null row/columns.

Mapping module
--------------

   - eq_reads:                               Compare reads accounting for multicontacts

   - `get_intersection <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.get_intersection>`_: Merges the two files corresponding to each reads sides. Reads found in both                                             files are merged and written in an output file.                                                                                          Dealing with multiple contacts:                                             - a pairwise contact is created for each possible combnation of the                                             multicontacts.                                             - if no other fragment from this read are mapped than, both are kept                                             - otherwise, they are merged into one longer (as if they were mapped                                             in the positive strand)

   - gt_reads:                               Compare reads accounting for multicontacts

   - `merge_2d_beds <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.merge_2d_beds>`_: Merge two result files (file resulting from get_intersection or from                                             the filtering) into one.

   - merge_bams:                             Merge two bam files with samtools into one.

Mapping analyze module
----------------------

   - `correlate_matrices <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.analyze.correlate_matrices>`_ [#first]_ [#second]_: Compare the interactions of two Hi-C matrices at a given distance,                                             with Spearman rank correlation.                                                                                          Also computes the SCC reproducibility score as in HiCrep (see                                             https://doi.org/10.1101/gr.220640.117).

   - `eig_correlate_matrices <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.analyze.eig_correlate_matrices>`_ [#first]_ [#second]_: Compare the interactions of two Hi-C matrices using their 6 first                                             eigenvectors, with Pearson correlation

   - fragment_size [#first]_:                Plots the distribution of dangling-ends lengths

   - get_reproducibility:                    Compute reproducibility score similarly to HiC-spector                                             (https://doi.org/10.1093/bioinformatics/btx152)

   - `hic_map <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.analyze.hic_map>`_ [#first]_ [#second]_: function to retrieve data from HiC-data object. Data can be stored as                                             a square matrix, or drawn using matplotlib

   - `insert_sizes <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.analyze.insert_sizes>`_ [#first]_: Deprecated function, use fragment_size

   - `plot_distance_vs_interactions <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.analyze.plot_distance_vs_interactions>`_ [#first]_: Plot the number of interactions observed versus the genomic distance between                                             the mapped ends of the read. The slope is expected to be around -1, in                                             logarithmic scale and between 700 kb and 10 Mb (according to the prediction                                             of the fractal globule model).

   - `plot_genomic_distribution <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.analyze.plot_genomic_distribution>`_ [#first]_ [#second]_: Plot the number of reads in bins along the genome (or along a given                                             chromosome).

   - `plot_iterative_mapping <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.analyze.plot_iterative_mapping>`_ [#first]_: Plots the number of reads mapped at each step of the mapping process (in the                                             case of the iterative mapping, each step is mapping process with a given                                             size of fragments).

   - `plot_strand_bias_by_distance <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.analyze.plot_strand_bias_by_distance>`_ [#first]_: Classify reads into four categories depending on the strand on which each                                             of its end is mapped, and plots the proportion of each of these categories                                             in function of the genomic distance between them.                                                                                          Only full mapped reads mapped on two diferent restriction fragments (still                                             same chromosome) are considered.                                                                                          The four categories are:                                             - Both read-ends mapped on the same strand (forward)                                             - Both read-ends mapped on the same strand (reverse)                                             - Both read-ends mapped on the different strand (facing), like extra-dangling-ends                                             - Both read-ends mapped on the different strand (opposed), like extra-self-circles

Mapping filter module
---------------------

   - `apply_filter <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.filter.apply_filter>`_ [#second]_: Create a new file with reads filtered

   - `filter_reads <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.filter.filter_reads>`_ [#second]_: Filter mapped pair of reads in order to remove experimental artifacts (e.g.                                             dangling-ends, self-circle, PCR artifacts

Mapping full_mapper module
--------------------------

   - fast_fragment_mapping:                  Maps FASTQ reads to an indexed reference genome with the knowledge of                                             the restriction enzyme used (fragment-based mapping).

   - `full_mapping <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.full_mapper.full_mapping>`_: Maps FASTQ reads to an indexed reference genome. Mapping can be done either                                             without knowledge of the restriction enzyme used, or for experiments                                             performed without one, like Micro-C (iterative mapping), or using the                                             ligation sites created from the digested ends (fragment-based mapping).

   - transform_fastq:                        Given a FASTQ file it can split it into chunks of a given number of reads,                                             trim each read according to a start/end positions or split them into                                             restriction enzyme fragments

Mapping restriction_enzymes module
----------------------------------

   - `map_re_sites <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.restriction_enzymes.map_re_sites>`_: map all restriction enzyme (RE) sites of a given enzyme in a genome.                                             Position of a RE site is defined as the genomic coordinate of the first                                             nucleotide after the first cut (genomic coordinate starts at 1).                                                                                                                                       In the case of HindIII the genomic coordinate is this one:                                                                                          123456 789

   - iupac2regex:                            Convert target sites with IUPAC nomenclature to regex pattern

   - religateds:                             returns the resulting list of all possible sequences after religation of two                                             digested and repaired ends.

   - `repaired <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.mapping.restriction_enzymes.repaired>`_: returns the resulting sequence after reparation of two digested and repaired                                             ends, marking dangling ends.

   - identify_re:                            Search most probable restriction enzyme used in the Hi-C experiment.                                             Uses binomial test and some heuristics.

   - map_re_sites_nochunk:                   map all restriction enzyme (RE) sites of a given enzyme in a genome.                                             Position of a RE site is defined as the genomic coordinate of the first                                             nucleotide after the first cut (genomic coordinate starts at 1).                                                                                                                                       In the case of HindIII the genomic coordinate is this one:                                                                                          123456 789

Modelling impmodel module
-------------------------

   - `load_impmodel_from_cmm <http://3dgenomes.github.io/TADbit/reference/reference_modelling_impmodel.html#pytadbit.modelling.impmodel.load_impmodel_from_cmm>`_: Loads an IMPmodel object using an cmm file of the form:

   - `load_impmodel_from_xyz <http://3dgenomes.github.io/TADbit/reference/reference_modelling_impmodel.html#pytadbit.modelling.impmodel.load_impmodel_from_xyz>`_: Loads an IMPmodel object using an xyz file of the form:

IMPmodel class
++++++++++++++
    A container for the IMP modeling results.

      - `objective_function <http://3dgenomes.github.io/TADbit/reference/reference_modelling_impmodel.html#pytadbit.modelling.impmodel.IMPmodel.objective_function>`_ [#first]_: This function plots the objective function value per each Monte-Carlo                                             step.

Modelling structuralmodel module
--------------------------------

IMPmodel class
++++++++++++++
    A container for the IMP modeling results.

      - accessible_surface [#first]_:        Calculates a mesh surface around the model (distance equal to input                                             **radius**) and checks if each point of this mesh could be replaced by                                             an object (i.e. a protein) of a given **radius**                                                                                          Outer part of the model can be excluded from the estimation of                                             accessible surface, as the occupancy outside the model is unknown (see                                             superradius option).

      - center_of_mass:                      Gives the center of mass of a model

      - contour:                             Total length of the model

      - cube_side:                           Calculates the side of a cube containing the model.

      - cube_volume:                         Calculates the volume of a cube containing the model.

      - `distance <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.distance>`_: Calculates the distance between one point of the model and an external                                             coordinate

      - inaccessible_particles:              Gives the number of loci/particles that are accessible to an object                                             (i.e. a protein) of a given size.

      - longest_axe:                         Gives the distance between most distant particles of the model

      - min_max_by_axis:                     Calculates the minimum and maximum coordinates of the model

      - persistence_length:                  Calculates the persistence length (Lp) of given section of the model.                                             Persistence length is calculated according to [Bystricky2004] :

      - radius_of_gyration:                  Calculates the radius of gyration or gyradius of the model                                                                                          Defined as:

      - shortest_axe:                        Minimum distance between two particles in the model

      - view_model [#first]_:                Visualize a selected model in the three dimensions. (either with Chimera                                             or through matplotlib).

      - `write_cmm <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.write_cmm>`_ [#second]_: Save a model in the cmm format, read by Chimera                                             (http://www.cgl.ucsf.edu/chimera).                                                                                          **Note:** If none of model_num, models or cluster parameter are set,                                             ALL the models will be written.

      - `write_xyz <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.write_xyz>`_ [#second]_: Writes a xyz file containing the 3D coordinates of each particle in the                                             model.                                             Outfile is tab separated column with the bead number being the                                             first column, then the genomic coordinate and finally the 3                                             coordinates X, Y and Z                                                                                          **Note:** If none of model_num, models or cluster parameter are set,                                             ALL the models will be written.

      - `write_xyz_babel <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.write_xyz_babel>`_ [#second]_: Writes a xyz file containing the 3D coordinates of each particle in the                                             model using a file format compatible with babel                                             (http://openbabel.org/wiki/XYZ_%28format%29).                                             Outfile is tab separated column with the bead number being the                                             first column, then the genomic coordinate and finally the 3                                             coordinates X, Y and Z                                             **Note:** If none of model_num, models or cluster parameter are set,                                             ALL the models will be written.

Modelling structuralmodels module
---------------------------------

   - `load_structuralmodels <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.load_structuralmodels>`_: Loads StructuralModels from a file                                             (generated with                                             save_models).

StructuralModels class
++++++++++++++++++++++
    This class contains three-dimensional models generated from a single Hi-C
    data. They can be reached either by their index (integer representing their
    rank according to objective function value), or by their IMP random intial
    number (as string).

      - `accessibility <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.accessibility>`_ [#first]_ [#second]_: Calculates a mesh surface around the model (distance equal to input                                             **radius**) and checks if each point of this mesh could be replaced by                                             an object (i.e. a protein) of a given **radius**                                                                                          Outer part of the model can be excluded from the estimation of                                             accessible surface, as the occupancy outside the model is unkown (see                                             superradius option).

      - `align_models <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.align_models>`_: Three-dimensional aligner for structural models.

      - `angle_between_3_particles <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.angle_between_3_particles>`_: Calculates the angle between 3 particles.                                                                                                                                       Given three particles A, B and C, the angle g (angle ACB, shown below):

      - `average_model <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.average_model>`_: Builds and returns an average model representing a given group of models

      - `centroid_model <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.centroid_model>`_: Estimates and returns the centroid model of a given group of models.

      - `cluster_analysis_dendrogram <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.cluster_analysis_dendrogram>`_ [#first]_: Representation of the clustering results. The length of the leaves if                                             proportional to the final objective function value of each model. The                                             branch widths are proportional to the number of models in a given                                             cluster (or group of clusters, if it is an internal branch).

      - `cluster_models <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.cluster_models>`_: This function performs a clustering analysis of the generated models                                             based on structural comparison. The result will be stored in                                             StructuralModels.clusters                                                                                          Clustering is done according to a score of pairwise comparison                                             calculated as:

      - `contact_map <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.contact_map>`_ [#first]_ [#second]_: Plots a contact map representing the frequency of interaction (defined                                             by a distance cutoff) between two particles.

      - `correlate_with_real_data <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.correlate_with_real_data>`_ [#first]_: Plots the result of a correlation between a given group of models and                                             original Hi-C data.

      - `deconvolve <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.deconvolve>`_ [#first]_: This function performs a deconvolution analysis of a given froup of models.                                             It first clusters models based on structural comparison (dRMSD), and                                             then, performs a differential contact map between each possible pair                                             of cluster.

      - `define_best_models <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.define_best_models>`_: Defines the number of top models (based on the objective function) to                                             keep. If keep_all is set to True in                                             generate_3d_models or in                                             model_region, then the full set                                             of models (n_models parameter) will be used, otherwise only the n_keep                                             models will be available.

      - `density_plot <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.density_plot>`_ [#first]_ [#second]_: Plots the number of nucleotides per nm of chromatin vs the modeled                                             region bins.

      - `dihedral_angle <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.dihedral_angle>`_: Calculates the dihedral angle between 2 planes formed by 5 particles                                             (one common to both planes).

      - `fetch_model_by_rand_init <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.fetch_model_by_rand_init>`_: Models are stored according to their objective function value (first                                             best), but in order to reproduce a model, we need its initial random                                             number. This method helps to fetch the model corresponding to a given                                             initial random number stored under                                             StructuralModels.models[N]['rand_init'].

      - `get_contact_matrix <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.get_contact_matrix>`_: Returns a matrix with the number of interactions observed below a given                                             cutoff distance.

      - `get_persistence_length <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.get_persistence_length>`_ [#first]_ [#second]_: Calculates the persistence length (Lp) of given section of the model.                                             Persistence length is calculated according to [Bystricky2004] :

      - `infer_unrestrained_particle_coords <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.infer_unrestrained_particle_coords>`_: if a given particle (and direct neighbors) have no restraints. Infer                                             the coordinates by linear interpolation using closest particles with                                             restraints.

      - `interactions <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.interactions>`_ [#first]_ [#second]_: Plots, for each particle, the number of interactions (particles closer                                             than the given cut-off). The value given is the average for all models.

      - `median_3d_dist <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.median_3d_dist>`_ [#first]_: Computes the median distance between two particles over a set of models

      - `model_consistency <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.model_consistency>`_ [#first]_ [#second]_: Plots the particle consistency, over a given set of models, vs the                                             modeled region bins. The consistency is a measure of the variability                                             (or stability) of the modeled region (the higher the consistency value,                                             the higher stability).

      - `objective_function_model <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.objective_function_model>`_ [#first]_: This function plots the objective function value per each Monte-Carlo                                             step

      - `particle_coordinates <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.particle_coordinates>`_: Returns the mean coordinate of a given particle in a group of models.

      - `save_models <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.save_models>`_ [#second]_: Saves all the models in pickle format (python object written to disk).

      - `view_centroid <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.view_centroid>`_: shortcut for                                             view_models(tool='plot', show='highlighted', highlight='centroid')

      - `view_models <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.view_models>`_ [#first]_: Visualize a selected model in the three dimensions (either with Chimera                                             or through matplotlib).

      - `walking_angle <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.walking_angle>`_ [#first]_ [#second]_: Plots the angle between successive loci in a given model or set of                                             models. In order to limit the noise of the measure angle is calculated                                             between 3 loci, between each are two other loci. E.g. in the scheme                                             bellow, angle are calculated between loci A, D and G.

      - `walking_dihedral <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.walking_dihedral>`_ [#first]_ [#second]_: Plots the dihedral angle between successive planes. A plane is formed by                                             3 successive loci.

      - `zscore_plot <http://3dgenomes.github.io/TADbit/reference/reference_modelling_structuralmodels.html#pytadbit.modelling.structuralmodels.StructuralModels.zscore_plot>`_ [#first]_: Generate 3 plots. Two heatmaps of the Z-scores used for modeling, one                                             of which is binary showing in red Z-scores higher than upper cut-off;                                             and in blue Z-scores lower than lower cut-off. Last plot is an histogram                                             of the distribution of Z-scores, showing selected regions. Histogram                                             also shows the fit to normal distribution.

Parsers bed_parser module
-------------------------

   - parse_bed:                              simple BED and BEDgraph parser that only checks for the fields 1, 2, 3 and 5                                             (or 1, 2 and 3 if 5 not availbale).

   - parse_mappability_bedGraph:             parse BEDgraph containing mappability.                                             GEM mappability file obtained with:                                                                                          gem-indexer -i hg38.fa -o hg38                                             gem-mappability -I hg38.gem -l 50 -o hg38.50mer -T 8                                             gem-2-wig -I hg38.gem -i hg38.50mer.mappability -o hg38.50mer                                             wigToBigWig hg38.50mer.wig hg38.50mer.sizes hg38.50mer.bw                                             bigWigToBedGraph hg38.50mer.bw  hg38.50mer.bedGraph

Parsers cooler_parser module
----------------------------

   - cooler_file:                            Cooler file wrapper.

   - close:                                  Copy remaining buffer to file, index the pixelsand complete information

   - create_bins:                            Write bins to cooler file.

   - prepare_matrix:                         Prepare matrix datasets to be written as chunks.

   - write_bins:                             Write the bins table.

   - write_indexes:                          Write the indexes from existing bins and pixels.

   - write_info:                             Write the file description and metadata attributes.

   - write_iter:                             Write bin1, bin2, value to buffer. When the chunk number changes the buffer                                             is written to the h5py file.

   - write_regions:                          Write the regions table.

   - write_weights:                          Write the weights in the bins table.

   - is_cooler:                              Check if file is a cooler and contains the wanted resolution

   - parse_cooler:                           Read matrix stored in cooler

   - parse_header:                           Read matrix header stored in cooler

   - rlencode:                               Run length encoding.                                             Based on http://stackoverflow.com/a/32681075, which is based on the rle                                             function from R.                                                                                          Parameters                                             ----------                                             x : 1D array_like                                             Input array to encode                                             dropna: bool, optional                                             Drop all runs of NaNs.                                                                                          Returns                                             -------                                             start positions, run lengths, run values

Parsers genome_parser module
----------------------------

   - `parse_fasta <http://3dgenomes.github.io/TADbit/reference/reference_parser.html#pytadbit.parsers.genome_parser.parse_fasta>`_: Parse a list of fasta files, or just one fasta.                                                                                          WARNING: The order is important

   - get_gc_content:                         Get GC content by bins of a given size. Ns are nottaken into account in the                                             calculation, only the number of Gs and Cs over As, Ts, Gs and Cs

Parsers hic_bam_parser module
-----------------------------

   - bed2D_to_BAMhic:                        function adapted from Enrique Vidal <enrique.vidal@crg.eu> scipt to convert                                             2D beds into compressed BAM format.                                                                                          Gets the *_both_filled_map.tsv contacts from TADbit (and the corresponding                                             filter files) and outputs a modified indexed BAM with the following fields:                                                                                          - read ID                                             - filtering flag (see codes in header)                                             - chromosome ID of the first pair of the contact                                             - genomic position of the first pair of the contact                                             - MAPQ set to 0                                             - pseudo CIGAR with sequence length and info about current copy (P: first copy, S: second copy)                                             - chromosome ID of the second pair of the contact                                             - genomic position of the second pair of the contact                                             - mapped length of the second pair of the contact                                             - sequence is missing (*)                                             - quality is missing (*)                                             - TC tag indicating single (1) or multi contact (3 6

   - get_biases_region:                      Retrieve biases, decay, and bad bins from a dictionary, and re-index it                                             according to a region of interest.

   - get_filters:                            get all filters

Parsers hic_parser module
-------------------------

   - `read_matrix <http://3dgenomes.github.io/TADbit/reference/reference_parser.html#pytadbit.parsers.hic_parser.read_matrix>`_: Read and checks a matrix from a file (using                                             autoreader) or a list.

   - load_hic_data_from_bam:                 

   - `load_hic_data_from_reads <http://3dgenomes.github.io/TADbit/reference/reference_parser.html#pytadbit.parsers.hic_parser.load_hic_data_from_reads>`_: 

   - abc_reader:                             Read matrix stored in 3 column format (bin1, bin2, value)

   - autoreader:                             Auto-detect matrix format of HiC data file.

   - is_asymmetric:                          Helper functions for the autoreader.

   - is_asymmetric_dico:                     Helper functions for the optimal_reader

   - optimal_reader:                         Reads a matrix generated by TADbit.                                             Can be slower than autoreader, but uses almost a third of the memory

   - symmetrize:                             Make a matrix symmetric by summing two halves of the matrix

   - symmetrize_dico:                        Make an HiC_data object symmetric by summing two halves of the matrix

AutoReadFail class
++++++++++++++++++
                      Exception to handle failed autoreader.

Parsers map_parser module
-------------------------

   - parse_map:                              Parse map files                                                                                          Keep a summary of the results into 2 tab-separated files that will contain 6                                             columns: read ID, Chromosome, position, strand (either 0 or 1), mapped                                             sequence lebgth, position of the closest upstream RE site, position of                                             the closest downstream RE site.                                                                                          The position of reads mapped on reverse strand will be computed from the end of                                             the read (original position + read length - 1)

Parsers sam_parser module
-------------------------

   - parse_gem_3c:                           Parse gem 3c sam file using pysam tools.

   - `parse_sam <http://3dgenomes.github.io/TADbit/reference/reference_parser.html#pytadbit.parsers.sam_parser.parse_sam>`_: Parse sam/bam file using pysam tools.                                                                                          Keep a summary of the results into 2 tab-separated files that will contain 6                                             columns: read ID, Chromosome, position, strand (either 0 or 1), mapped                                             sequence lebgth, position of the closest upstream RE site, position of                                             the closest downstream RE site

Parsers tad_parser module
-------------------------

   - `parse_tads <http://3dgenomes.github.io/TADbit/reference/reference_parser.html#pytadbit.parsers.tad_parser.parse_tads>`_: Parse a tab separated value file that contains the list of TADs of a given                                             experiment. This file might have been generated whith the                                             print_result_R or with the R binding for tadbit

Tad_clustering tad_cmo module
-----------------------------

   - core_nw:                                Core of the fast Needleman-Wunsch algorithm that aligns matrices

   - core_nw_long:                           Core of the long Needleman-Wunsch algorithm that aligns matrices

   - `optimal_cmo <http://3dgenomes.github.io/TADbit/reference/reference_clustering.html#pytadbit.tad_clustering.tad_cmo.optimal_cmo>`_: Calculates the optimal contact map overlap between 2 matrices                                                                                          TODO: make the selection of number of eigen vectors automatic or relying on                                             the summed contribution (e.g. select the EVs that sum 80% of the info)

   - virgin_score:                           Fill a matrix with zeros, except first row and first column filled with     multiple values of penalty.

Tadbit module
-------------

   - `batch_tadbit <http://3dgenomes.github.io/TADbit/reference/reference_tadbit.html#pytadbit.tadbit.batch_tadbit>`_ [#second]_: Use tadbit on directories of data files.                                             All files in the specified directory will be considered data file. The                                             presence of non data files will cause the function to either crash or                                             produce aberrant results.                                                                                          Each file has to contain the data for a single unit/chromosome. The                                             files can be separated in sub-directories corresponding to single                                             experiments or any other organization. Data files that should be                                             considered replicates have to start with the same characters, until                                             the character sep. For instance, all replicates of the unit                                             'chr1' should start with 'chr1\_', using the default value of sep.                                                                                          The data files are read through read.delim. You can pass options                                             to read.delim through the list read_options. For instance                                             if the files have no header, use read_options=list(header=FALSE) and if                                             they also have row names, read_options=list(header=FALSE, row.names=1).                                                                                          Other arguments such as max_size, n_CPU and verbose are passed to                                             tadbit.                                                                                          NOTE: only used externally, not from Chromosome

   - `tadbit <http://3dgenomes.github.io/TADbit/reference/reference_tadbit.html#pytadbit.tadbit.tadbit>`_: The TADbit algorithm works on raw chromosome interaction count data.                                             The normalization is neither necessary nor recommended,                                             since the data is assumed to be discrete counts.                                                                                          TADbit is a breakpoint detection algorithm that returns the optimal                                             segmentation of the chromosome under BIC-penalized likelihood. The                                             model assumes that counts have a Poisson distribution and that the                                             expected value of the counts decreases like a power-law with the                                             linear distance on the chromosome. This expected value of the counts                                             at position (i,j) is corrected by the counts at diagonal positions                                             (i,i) and (j,j). This normalizes for different restriction enzyme                                             site densities and 'mappability' of the reads in case a bin contains                                             repeated regions.

Utils extraviews module
-----------------------

   - colorize:                               Colorize with ANSII colors a string for printing in shell. this acording to                                             a given number between 0 and 10

   - nicer:                                  writes resolution number for human beings.

   - add_subplot_axes:                       from https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib/35966183

   - augmented_dendrogram [#first]_:         

   - `chimera_view <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.extraviews.chimera_view>`_ [#first]_: Open a list of .cmm files with Chimera (http://www.cgl.ucsf.edu/chimera)                                             to view models.

   - `color_residues <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.extraviews.color_residues>`_: Function to color residues from blue to red.

   - `compare_models <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.extraviews.compare_models>`_: Plots the difference of contact maps of two group of structural models.

   - pcolormesh_45deg [#first]_:             Draw triangular matrix

   - `plot_2d_optimization_result <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.extraviews.plot_2d_optimization_result>`_ [#first]_: A grid of heatmaps representing the result of the optimization.                                             The maps will be divided in different pages depending on the 'scale' and 'kbending' values.                                             In each page there will be different maps depending the 'maxdist' values.                                             Each map has 'upfreq' values along the x-axes, and 'lowfreq' values along the y-axes.

   - `plot_3d_model <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.extraviews.plot_3d_model>`_ [#first]_: Given a 3 lists of coordinates (x, y, z) plots a three-dimentional model                                             using matplotlib

   - `plot_3d_optimization_result <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.extraviews.plot_3d_optimization_result>`_: Displays a three dimensional scatter plot representing the result of the                                             optimization.

   - plot_HiC_matrix [#first]_:              Plot HiC matrix with histogram of values inside color bar.

   - `tad_border_coloring <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.extraviews.tad_border_coloring>`_: Colors TAD borders from blue to red (bad to good score). TAD are displayed                                             in scale of grey, from light to dark grey (first to last particle in the                                             TAD)

   - `tad_coloring <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.extraviews.tad_coloring>`_: Colors TADs from blue to red (first to last TAD). TAD borders are displayed                                             in scale of grey, from light to dark grey (again first to last border)

Utils fastq_utils module
------------------------

   - `quality_plot <http://3dgenomes.github.io/TADbit/reference/reference_mapping.html#pytadbit.utils.fastq_utils.quality_plot>`_ [#first]_: Plots the sequencing quality of a given FASTQ file. If a restrinction enzyme                                             (RE) name is provided, can also represent the distribution of digested and                                             undigested RE sites and estimate an expected proportion of dangling-ends.                                                                                          Proportion of dangling-ends is inferred by counting the number of times a                                             dangling-end site, is found at the beginning of any of the reads (divided by                                             the number of reads).

Utils file_handling module
--------------------------

   - `magic_open <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.file_handling.magic_open>`_: To read uncompressed zip gzip bzip2 or tar.xx files

   - which:                                  stackoverflow: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

   - `get_free_space_mb <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.file_handling.get_free_space_mb>`_: Return folder/drive free space (in bytes)                                                                                          Based on stackoverflow answer: http://stackoverflow.com/questions/51658/cross-platform-space-remaining-on-volume-using-python

   - is_fastq:                               Check if a given file is in fastq format

   - wc:                                     Pythonic way to count lines

Utils hic_filtering module
--------------------------

   - `hic_filtering_for_modelling <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.hic_filtering.hic_filtering_for_modelling>`_ [#first]_: Call filtering function, to remove artifactual columns in a given Hi-C                                             matrix. This function will detect columns with very low interaction                                             counts; and columns with NaN values (in this case NaN will be replaced                                             by zero in the original Hi-C data matrix). Filtered out columns will be                                             stored in the dictionary Experiment._zeros.

   - filter_by_mean [#first]_:               fits the distribution of Hi-C interaction count by column in the matrix to                                             a polynomial. Then searches for the first possible

   - filter_by_zero_count:                   

   - filter_by_cis_percentage [#first]_:     Define artifactual columns with either too low or too high counts of                                             interactions by compraing their percentage of cis interactions                                             (inter-chromosomal).

Utils hmm module
----------------

   - best_path:                              Viterbi algorithm with backpointers

   - gaussian_prob:                          of x to follow the gaussian with given E                                             https://en.wikipedia.org/wiki/Normal_distribution

Utils normalize_hic module
--------------------------

   - iterative:                              Implementation of iterative correction Imakaev 2012

   - expected:                               Computes the expected values by averaging observed interactions at a given                                             distance in a given HiC matrix.

Utils tadmaths module
---------------------

   - zscore:                                 Calculates the log10, Z-score of a given list of values.

   - `calinski_harabasz <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.tadmaths.calinski_harabasz>`_: Implementation of the CH score [CalinskiHarabasz1974], that has shown to be                                             one the most accurate way to compare clustering methods                                             [MilliganCooper1985] [Tibshirani2001].                                                                                          The CH score is:

   - mad:                                    Median Absolute Deviation: a "Robust" version of standard deviation.                                             Indices variability of the sample.                                             https://en.wikipedia.org/wiki/Median_absolute_deviation

   - mean_none:                              Calculates the mean of a list of values without taking into account the None

   - newton_raphson:                         Newton-Raphson method as defined in:                                             http://www.maths.tcd.ie/~ryan/TeachingArchive/161/teaching/newton-raphson.c.html                                             used to search for the persistence length of a given model.

   - right_double_mad:                       Double Median Absolute Deviation: a 'Robust' version of standard deviation.                                             Indices variability of the sample.                                             http://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers

Interpolate class
+++++++++++++++++
                      Simple linear interpolation, to be used when the one from scipy is not
                      available.

Utils three_dim_stats module
----------------------------

   - `angle_between_3_points <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.angle_between_3_points>`_: Calculates the angle between 3 particles                                                                                          Given three particles A, B and C, the angle g (angle ACB, shown below):

   - build_mesh:                             Main function for the calculation of the accessibility of a model.

   - `calc_eqv_rmsd <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.calc_eqv_rmsd>`_: Calculates the RMSD, dRMSD, the number of equivalent positions and a score                                             combining these three measures. The measure are done between a group of                                             models in a one against all manner.

   - `dihedral <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.dihedral>`_: Calculates dihedral angle between 4 points in 3D (array with x,y,z)

   - fast_square_distance:                   Calculates the square distance between two coordinates.

   - find_angle_rotation_improve_x:          Finds the rotation angle needed to face the longest edge of the molecule

   - `generate_circle_points <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.generate_circle_points>`_: Returns list of 3d coordinates of points on a circle using the                                             Rodrigues rotation formula.                                                                                          see *Murray, G. (2013). Rotation About an Arbitrary Axis in 3 Dimensions*                                             for details

   - generate_sphere_points:                 Returns list of 3d coordinates of points on a sphere using the                                             Golden Section Spiral algorithm.

   - get_center_of_mass:                     get the center of mass of a given object with list of x, y, z coordinates

   - mass_center:                            Transforms coordinates according to the center of mass

   - mmp_score [#first]_:                    

   - rotate_among_y_axis:                    Rotate and object with a list of x, y, z coordinates among its center of                                             mass

   - `square_distance <http://3dgenomes.github.io/TADbit/reference/reference_utils.html#pytadbit.utils.three_dim_stats.square_distance>`_: Calculates the square distance between two particles.


.. [#first] functions generating plots

.. [#second] functions writing text files

