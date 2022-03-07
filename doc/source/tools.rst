TADbit tools
============

TADbit also provides a set of command line tools that are installed with
the library and that cover the main functionalities.

This tools are idependent but share a working directory where a local
database is created to store the input/outputs of each process or job,
and also some statistics.

Each of these tools has extensive help, so we will here review only
their general usage and function.

TADbit map
----------

.. code:: ipython3

    Map Hi-C reads and organize results in an output working directory
    
    usage: tadbit map [-h] [--skip_mapping] -w PATH --fastq PATH [--fastq2 PATH] --index PATH
                      [--genome PATH [PATH ...]] --read INT --renz STR [STR ...]
                      [--chr_name STR [STR ...]] [--tmp PATH] [--tmpdb PATH] [--noX] [--iterative]
                      [--fast_fragment] [--windows WINDOWS [WINDOWS ...]] [--species STR]
                      [--descr LIST [LIST ...]] [--skip] [--keep_tmp] [-C CPUS] [--mapper STR]
                      [--mapper_binary STR] [--mapper_param MAPPER_PARAM [MAPPER_PARAM ...]]
    
    optional arguments:
      -h, --help               show this help message and exit
    
    General options:
      --skip_mapping           generate a Hi-C specific quality plot from FASTQ and exits
      -w PATH, --workdir PATH  path to an output folder.
      --fastq PATH             path to a FASTQ files (can be compressed files)
      --fastq2 PATH            (beta) path to a FASTQ file of read 2 (can be compressed files).
                               Needed for fast_fragment
      --index PATH             paths to file(s) with indexed FASTA files of the reference genome.
      --genome PATH [PATH ...]
                               paths to file(s) with FASTA files of the reference genome. Needed
                               for fast_fragment mapping. If many, files will be concatenated.
                               I.e.: --genome chr_1.fa chr_2.fa In this last case, order is
                               important or the rest of the analysis. Note: it can also be the path
                               to a previously parsed genome in pickle format.
      --read INT               read number
      --renz STR [STR ...]     restriction enzyme name(s). Use "--renz CHECK" to search for most
                               probable and exit; and use "--renz NONE" to avoid using RE site
                               information.
      --chr_name STR [STR ...]
                               [fasta header] chromosome name(s). Used in the same order as data.
      --tmp PATH               path to a temporary directory (default next to "workdir" directory)
      --tmpdb PATH             if provided uses this directory to manipulate the database
      --noX                    no display server (X screen)
      --skip                   [DEBUG] in case already mapped.
      --keep_tmp               [DEBUG] keep temporary files.
    
    Mapping options:
      --iterative              default mapping strategy is fragment based use this flag for
                               iterative mapping
      --fast_fragment          (beta) use fast fragment mapping. Both fastq files are mapped using
                               fragment based mapping in GEM v3. The output file is an intersected
                               read file than can be used directly in tadbit filter (no tadbit
                               parse needed). Access to samtools is needed for fast_fragment to
                               work. --fastq2 and --genome needs to be specified and --read value
                               should be 0.
      --windows WINDOWS [WINDOWS ...]
                               defines windows to be used to trim the input FASTQ reads, for
                               example an iterative mapping can be defined as: "--windows 1:20 1:25
                               1:30 1:35 1:40 1:45 1:50". But this parameter can also be used for
                               fragment based mapping if for example pair-end reads are both in the
                               same FASTQ, for example: "--windows 1:50" (if the length of the
                               reads is 100). Note: that the numbers are both inclusive.
      -C CPUS, --cpus CPUS     [32] Maximum number of CPU cores available in the execution host. If
                               higher than 1, tasks with multi-threading capabilities will enabled
                               (if 0 all available) cores will be used
      --mapper STR             [gem] mapper used, options are gem, bowtie2 or hisat2
      --mapper_binary STR      [None] path to mapper binary
      --mapper_param MAPPER_PARAM [MAPPER_PARAM ...]
                               any parameter that could be passed to the GEM, BOWTIE2 or HISAT2
                               mapper. e.g. if we want to set the proportion of mismatches to 0.05
                               and the maximum indel length to 10, (in GEM v2 it would be: -e 0.05
                               --max-big-indel-length 10), here we could write: "--mapper_param
                               e:0.05 max-big-indel-length:10". For BOWTIE2, GEM3 and HISAT2 you
                               can also pass directly the parameters enclosed between quotes like:
                               --mapper_param "-e 0.05 --alignment-local-min-score 15" IMPORTANT:
                               some options are incompatible with 3C-derived experiments.
    
    Descriptive, optional arguments:
      --species STR            species name
      --descr LIST [LIST ...]  extra descriptive fields each filed separated by coma, and inside
                               each, name and value separated by column:
                               --descr=cell:lymphoblast,flowcell:C68AEACXX,index:24nf

TADbit parse
------------

.. code:: ipython3

    Parse mapped Hi-C reads and get the intersection
    
    usage: tadbit parse [-h] [-w PATH] [--type STR] [--read INT] [--mapped1 PATHs [PATHs ...]]
                        [--mapped2 PATHs [PATHs ...]] [--renz STR] [--filter_chrom FILTER_CHROM]
                        [--skip] [--compress_input] [--tmpdb PATH] [--genome PATH [PATH ...]]
                        [--jobids INT [INT ...]] [--noX]
    
    optional arguments:
      -h, --help               show this help message and exit
    
    General options:
      -w PATH, --workdir PATH  path to working directory (generated with the tool tadbit mapper)
      --type STR               [0map]file type to be parser, MAP (GEM-mapper), SAM or BAM
      --read INT               In case only one of the reads needs to be parsed
      --filter_chrom FILTER_CHROM
                               default: --filter_chrom
                               "^(chr)?[A-Za-z]?[0-9]{0,3}[XVI]{0,3}(?:ito)?[A-Z-a-z]?(_dna)?$",
                               regexp to consider only chromosome names passing
      --skip                   [DEBUG] in case already mapped.
      --compress_input         Compress input mapped files when parsing is done. This is done in
                               background, while next MAP file is processed, or while reads are
                               sorted.
      --tmpdb PATH             if provided uses this directory to manipulate the database
      --genome PATH [PATH ...]
                               paths to file(s) with FASTA files of the reference genome. If many,
                               files will be concatenated. I.e.: --genome chr_1.fa chr_2.fa In this
                               last case, order is important or the rest of the analysis. Note: it
                               can also be the path to a previously parsed genome in pickle format.
      --jobids INT [INT ...]   Use as input data generated by a job with a given jobid(s). Use
                               tadbit describe to find out which. In this case one jobid can be
                               passed per read.
      --noX                    no display server (X screen)
    
    Mapped outside TADbit options:
      --mapped1 PATHs [PATHs ...]
                               paths to mapped bam files (first read-end)
      --mapped2 PATHs [PATHs ...]
                               paths to mapped bam files (second read-end)
      --renz STR               restriction enzyme name

TADbit filter
-------------

.. code:: ipython3

    Filter parsed Hi-C reads and get valid pair of reads to work with
    
    usage: tadbit filter [-h] [--force] [--resume] [--apply INT [INT ...]] [-w PATH] [-C CPUS]
                         [--noX] [--over_represented NUM] [--strict_duplicates]
                         [--max_frag_size NUM] [--min_frag_size NUM] [--re_proximity NUM]
                         [--mad NUM] [--max_f NUM] [--median NUM] [--tmpdb PATH]
                         [--pathids INT [INT ...]] [--compress_input] [--format {short,mid,long}]
                         [--valid] [--clean] [--samtools PATH]
    
    optional arguments:
      -h, --help               show this help message and exit
    
    General options:
      --force                  overwrite previously run job
      --resume                 use filters of previously run job
      -w PATH, --workdir PATH  path to working directory (generated with the tool tadbit mapper)
      -C CPUS, --cpus CPUS     [32] Maximum number of CPU cores available in the execution host. If
                               higher than 1, tasks with multi-threading capabilities will enabled
                               (if 0 all available) cores will be used
      --noX                    no display server (X screen)
      --tmpdb PATH             if provided uses this directory to manipulate the database
      --pathids INT [INT ...]  Use as input data generated by a job under a given pathids. Use
                               tadbit describe to find out which. To filter an intersected file
                               produced with tadbit map --fast_fragment only one PATHid is needed
                               otherwise one per read is needed, first for read 1, second for read
                               2.
      --compress_input         Compress input mapped files when parsing is done. This is done in
                               background, while next MAP file is processed, or while reads are
                               sorted.
      --samtools PATH          path samtools binary
    
    Storage options:
      --format {short,mid,long}
                               [0mid] for compression into pseudo-BAM format. Short contains only
                               positions of reads mapped, mid everything but restriction sites.
      --valid                  stores only valid-pairs discards filtered out reads.
      --clean                  remove intermediate files. WARNING: together with format "short" or
                               valid options, this may results in losing data
    
    Filtering options:
      --apply INT [INT ...]    [[1, 2, 3, 4, 6, 7, 9, 10]] Use filters to define a set os valid
                               pair of reads e.g.: '--apply 1 2 3 4 6 7 8 9'. Where these
                               numberscorrespond to: 1: self-circle, 2: dangling-end, 3: error, 4:
                               extra dangling-end, 5: too close from RES, 6: too short, 7: too
                               large, 8: over-represented, 9: duplicated, 10: random breaks
      --over_represented NUM   [0.001%] percentage of restriction-enzyme (RE) genomic fragments
                               with more coverage to exclude (possible PCR artifact).
      --strict_duplicates      by default reads are considered duplicates if they coincide in
                               genomic coordinates and strand; with strict_duplicates enabled, we
                               also ask to consider read length (WARNING: this option is called
                               strict, but it is more permissive)
      --max_frag_size NUM      [100000] to exclude large genomic RE fragments (probably resulting
                               from gaps in the reference genome)
      --min_frag_size NUM      [50] to exclude small genomic RE fragments (smaller than sequenced
                               reads)
      --re_proximity NUM       [5] to exclude read-ends falling too close from RE site (pseudo-
                               dangling-ends)
      --mad NUM                MAD fragment length normally computed from observed distribution
      --max_f NUM              Maximum fragment length normally computed from observed distribution
      --median NUM             Median fragment length normally computed from observed distribution

TADbit normalize
----------------

.. code:: ipython3

    Normalize Hi-C data and generates array of biases
    
    usage: tadbit normalize [-h] -w PATH -r INT [--bam PATH] [-j INT] [--max_njobs INT]
                            [--tmpdb PATH] [-C CPUS] [--normalize_only] [--noX]
                            [--normalization STR] [--biases_path BIASES_PATH] [--mappability PATH]
                            [--fasta PATH] [--renz STR] [--factor NUM] [--prop_data FLOAT]
                            [--seed INT] [--min_count INT] [--cis_limit CIS_LIMIT]
                            [--trans_limit TRANS_LIMIT] [--ratio_limit RATIO_LIMIT]
                            [--cistrans_filter] [--filter_only]
                            [-B CHR:POS1-POS2 [CHR:POS1-POS2 ...]] [-F INT [INT ...]] [--valid]
    
    normalize Hi-C data and generates array of biases
    
    optional arguments:
      -h, --help               show this help message and exit
    
    Required options:
      -w PATH, --workdir PATH  path to working directory (generated with the tool tadbit mapper)
      -r INT, --resolution INT
                               resolution at which to output matrices
    
    General options:
      --bam PATH               path to a TADbit-generated BAM file with all reads (other wise the
                               tool will guess from the working directory database)
      -j INT, --jobid INT      Use as input data generated by a job with a given jobid. Use tadbit
                               describe to find out which.
      --max_njobs INT          [100] Define maximum number of jobs for reading BAM file (set to
                               higher numbers for large files and low RAM memory).
      --tmpdb PATH             if provided uses this directory to manipulate the database
      -C CPUS, --cpus CPUS     [32] Maximum number of CPU cores available in the execution host. If
                               higher than 1, tasks with multi-threading capabilities will enabled
                               (if 0 all available) cores will be used
      --normalize_only         skip calculation of Cis-percentage and decay
      --noX                    no display server (X screen)
    
    Bin filtering options:
      --min_count INT          [None] minimum number of reads mapped to a bin (recommended value
                               could be 2500). If set this option overrides the perc_zero
                               filtering... This option is slightly slower.
      --cis_limit CIS_LIMIT    Maximum distance in bins at which to consider an interaction cis for
                               the filtering. By default it is the number of bins corresponding to
                               1Mb
      --trans_limit TRANS_LIMIT
                               Maximum distance in bins at which to consider an interaction trans
                               for the filtering. By default it is five times the cis_limit (if
                               also default, it would correspond to the number of bins needed to
                               reach 5Mb).
      --ratio_limit RATIO_LIMIT
                               [1.0] Minimum cis/trans (as defined with cis_limit and trans_limit
                               parameters) to filter out bins.
      --cistrans_filter        filter using cis-trans ratio.
      --filter_only            skip normalization
      -B CHR:POS1-POS2 [CHR:POS1-POS2 ...], --badcols CHR:POS1-POS2 [CHR:POS1-POS2 ...]
                               extra regions to be added to bad-columns (ingenomic position). e.g.:
                               --badcols 1:150000000-160000000 2:1200000-1300000
    
    Read filtering options:
      -F INT [INT ...], --filter INT [INT ...]
                               [[1, 2, 3, 4, 6, 7, 9, 10]] Use filters to define a set os valid
                               pair of reads e.g.: '--apply 1 2 3 4 8 9 10'. Where these
                               numberscorrespond to: 1: self-circle, 2: dangling-end, 3: error, 4:
                               extra dangling-end, 5: too close from RES, 6: too short, 7: too
                               large, 8: over-represented, 9: duplicated, 10: random breaks, 11:
                               trans-chromosomic
      --valid                  input BAM file contains only valid pairs (already filtered).
    
    Normalization options:
      --normalization STR      [Vanilla] normalization(s) to apply. Order matters. Choices:
                               Vanilla, ICE, SQRT, oneD, custom
      --biases_path BIASES_PATH
                               biases file to compute decay. REQUIRED with "custom" normalization.
                               Format: single column with header
      --mappability PATH       Path to mappability bedGraph file, required for oneD normalization.
                               Mappability file can be generated with GEM (example from the genomic FASTA file hg38.fa):
                               
                                    gem-indexer -i hg38.fa -o hg38
                                    gem-mappability -I hg38.gem -l 50 -o hg38.50mer -T 8
                                    gem-2-wig -I hg38.gem -i hg38.50mer.mappability -o hg38.50mer
                                    wigToBigWig hg38.50mer.wig hg38.50mer.sizes hg38.50mer.bw
                                    bigWigToBedGraph hg38.50mer.bw  hg38.50mer.bedGraph
      --fasta PATH             Path to FASTA file with genome sequence, to compute GC content and
                               number of restriction sites per bin. Required for oneD normalization
      --renz STR               restriction enzyme name(s). Required for oneD normalization
      --factor NUM             [1] target mean value of a cell after normalization (can be used to
                               weight experiments before merging)
      --prop_data FLOAT        [1] Only for oneD normalization: proportion of data to be used in
                               fitting (for very large datasets). Number between 0 and 1.
      --seed INT               [1] Only for oneD normalization: seed number for the random picking
                               of data when using the "prop_data" parameter

TADbit bin
----------

.. code:: ipython3

    Bin Hi-C data into matrices
    
    usage: tadbit bin [-h] -w PATH [--noX] -r INT [--bam PATH] [-j INT] [--force] [-q]
                      [--tmpdb PATH] [--nchunks NCHUNKS] [-C CPUS] [--chr_name STR [STR ...]]
                      [--matrix] [--cooler] [--rownames] [--plot] [--force_plot] [--only_plot] [-i]
                      [--triangular] [--xtick_rotation XTICK_ROTATION] [--cmap CMAP]
                      [--bad_color BAD_COLOR] [--format FORMAT] [--zrange ZRANGE]
                      [--transform {log2,log,none}] [--figsize FIGSIZE] [--tad_def TAD_DEF] [-c]
                      [-c2] [--biases PATH] [--norm STR [STR ...]] [-F INT [INT ...]] [--only_txt]
    
    optional arguments:
      -h, --help               show this help message and exit
    
    Required options:
      -w PATH, --workdir PATH  path to working directory (generated with the tool tadbit mapper)
      -r INT, --resolution INT
                               resolution at which to output matrices
    
    General options:
      --noX                    no display server (X screen)
      --bam PATH               path to a TADbit-generated BAM file with all reads (other wise the
                               tool will guess from the working directory database)
      -j INT, --jobid INT      Use as input data generated by a job with a given jobid. Use tadbit
                               describe to find out which.
      --force                  overwrite previously run job
      -q, --quiet              remove all messages
      --tmpdb PATH             if provided uses this directory to manipulate the database
      --nchunks NCHUNKS        maximum number of chunks into which to cut the BAM
      -C CPUS, --cpus CPUS     [32] Maximum number of CPU cores available in the execution host. If
                               higher than 1, tasks with multi-threading capabilities will enabled
                               (if 0 all available) cores will be used
      --chr_name STR [STR ...]
                               [fasta header] chromosome name(s). Order of chromosomes in the
                               output matrices.
    
    Read filtering options:
      -F INT [INT ...], --filter INT [INT ...]
                               [[1, 2, 3, 4, 6, 7, 9, 10]] Use filters to define a set os valid
                               pair of reads e.g.: '--apply 1 2 3 4 8 9 10'. Where these
                               numberscorrespond to: 0: nothing, 1: self-circle, 2: dangling-end,
                               3: error, 4: extra dangling-end, 5: too close from RES, 6: too
                               short, 7: too large, 8: over-represented, 9: duplicated, 10: random
                               breaks, 11: trans-chromosomic
    
    Normalization options:
      --biases PATH            path to file with pre-calculated biases by columns
      --norm STR [STR ...]     [['raw']] normalization(s) to apply. Choices are: [norm, decay, raw,
                               raw&decay]
    
    Output options:
      --matrix                 Write text matrix in multiple columns (square). By defaults matrices
                               are written in BED-like format (also only way to get a raw matrix
                               with all values including the ones in masked columns).
      --cooler                 Write i,j,v matrix in cooler format instead of text.
      --rownames               To store row names in the output text matrix. WARNING: when non-
                               matrix, results in two extra columns
      --only_plot              [False] Skip writing matrix in text format.
      -i, --interactive        [False] Open matplotlib interactive plot (nothing will be saved).
      -c , --coord             Coordinate of the region to retrieve. By default all genome,
                               arguments can be either one chromosome name, or the coordinate in
                               the form: "-c chr3:110000000-120000000"
      -c2 , --coord2           Coordinate of a second region to retrieve the matrix in the
                               intersection with the first region.
      --only_txt               Save only text file for matrices, not images
    
    Plotting options:
      --plot                   Plot matrix in desired format.
      --force_plot             Force plotting even with demoniacally big matrices (more than
                               5000x5000, or 1500x1500with interactive option).
      --triangular             [False] represents only half matrix. Note that this also results in
                               truly vectorial images of matrix.
      --xtick_rotation XTICK_ROTATION
                               [-25] x-tick rotation
      --cmap CMAP              [viridis] Matplotlib color map to use.
      --bad_color BAD_COLOR    [white] Matplotlib color to use on bins filtered out (only used with
                               normalized matrices, not raw).
      --format FORMAT          [png] plot file format.
      --zrange ZRANGE          Range, in log2 scale of the color scale. i.e.: --zrange=-2,2
      --transform {log2,log,none}
                               [log2] can be any of [log2, log, none]
      --figsize FIGSIZE        Range, in log2 scale of the color scale. default for triangular
                               matrices: --figsize=16,10 and for square matrices: --figsize=16,14
      --tad_def TAD_DEF        jobid with the TAD segmentation, alternatively a tsv file with tad
                               definition, columns: # start end score density

TADbit segment
--------------

.. code:: ipython3

    
    Finds TAD or compartment segmentation in Hi-C data.
    
    usage: tadbit segment [-h] -w PATH [--tmpdb PATH] [--nosql] [--all_bins] [--mreads PATH]
                          [--biases PATH] -r INT [--norm_matrix PATH] [--raw_matrix PATH]
                          [-F INT [INT ...]] [--noX] [--rich_in_A PATH] [--fasta PATH] [--savecorr]
                          [--fix_corr_scale] [--format FORMAT] [--n_evs INT]
                          [--ev_index INT [INT ...]] [--only_compartments] [--only_tads] [-v]
                          [-j INT] [-c STR [STR ...]] [--max_tad_size INT] [-C CPUS] [--force]
    
    optional arguments:
      -h, --help               show this help message and exit
    
    General options:
      -w PATH, --workdir PATH  path to working directory (generated with the tool tadbit mapper)
      --tmpdb PATH             if provided uses this directory to manipulate the database
      --nosql                  do not load/store data from/in sqlite database
      --all_bins               skip the filtering of bins for detection of TADs
      --mreads PATH            path valid-pairs file (TADbit BAM format)
      --biases PATH            path to file with precalculated biases by columns
      -r INT, --resolution INT
                               resolution at which to output matrices
      --norm_matrix PATH       path to a matrix file with normalized read counts
      --raw_matrix PATH        path to a matrix file with raw read counts
      -F INT [INT ...], --filter INT [INT ...]
                               [[1, 2, 3, 4, 6, 7, 9, 10]] Use filters to define a set os valid
                               pair of reads e.g.: '--apply 1 2 3 4 8 9 10'. Where these
                               numberscorrespond to: 1: self-circle, 2: dangling-end, 3: error, 4:
                               extra dangling-end, 5: too close from RES, 6: too short, 7: too
                               large, 8: over-represented, 9: duplicated, 10: random breaks, 11:
                               trans-chromosomic
      --noX                    no display server (X screen)
      --only_compartments      search A/B compartments using first eigen vector of the correlation
                               matrix
      --only_tads              search TAD boundaries break-point detection algorithm
      -v, --verbose            print more messages
      -j INT, --jobid INT      Use as input data generated by a job with a given jobid. Use tadbit
                               describe to find out which.
      -c STR [STR ...], --chromosomes STR [STR ...]
                               Name of the chromosomes on which to search for TADs or compartments.
      -C CPUS, --cpu CPUS      [32] Maximum number of CPU cores available in the execution host. If
                               higher than 1, tasks with multi-threading capabilities will enabled
                               (if 0 all available) cores will be used
      --force                  overwrite previously run job
    
    Compartment calling options:
      --rich_in_A PATH         path to a BED or bedGraph file with list of protein coding gene or
                               other active epigenetic mark, to be used to label compartments
                               instead of using the relative interaction count.
      --fasta PATH             Path to fasta file with genome sequence, to compute GC content and
                               use it to label compartments
      --savecorr               Save correlation matrix used to predict compartments.
      --fix_corr_scale         Correlation matrix plot scaled between correlation 1 and -1 instead
                               of maximum observed values.
      --format FORMAT          [png] file format for figures
      --n_evs INT              [3] Number of eigenvectors to store. if "-1" all eigenvectors will
                               be calculated
      --ev_index INT [INT ...]
                               list of indexes of eigenvectors capturing compartments signal (one
                               index per chromosome, in the same order as chromosomes in fasta
                               file). Example picking the first eigenvector for all chromosomes but
                               for chromosome 3: '--ev_index 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                               1 1 1 1
    
    TAD calling options:
      --max_tad_size INT       an integer defining the maximum size of TAD. Default defines it as
                               the number of rows/columns

TADbit model
------------

.. code:: ipython3

    Generates 3D models given an input interaction matrix and a set of input parameters
    
    usage: tadbit model [-h] -w PATH [--input_matrix PATH] [--rand INT] [--nmodels INT]
                        [--nkeep INT] [-j INT] [--optimization_id INT] [--restart_id INT]
                        [--fig_format STR] [--noX] [--corr STR] [--species STRING]
                        [--assembly STRING] [--cell STRING] [--exp_type STRING] [--project STRING]
                        [--crm NAME] [--beg INT] [--end INT] [--matrix_beg INT] [-r INT]
                        [--perc_zero FLOAT] [--smooth_factor INT] [--optimize] [--model]
                        [--model_ptadbit] [--force] [--maxdist LIST [LIST ...]]
                        [--upfreq LIST [LIST ...]] [--lowfreq LIST [LIST ...]]
                        [--scale LIST [LIST ...]] [--dcutoff LIST [LIST ...]]
                        [--container LIST [LIST ...]] [--analyze] [-C CPUS] [--job_list]
                        [--nmodels_per_job INT] [--cpus_per_job INT] [--concurrent_jobs INT]
                        [--timeout_job INT] [--script_cmd STR] [--script_args STR]
                        [--script_template STR] [--tmpdb PATH] [--analyze_list INT [INT ...]]
                        [--not_write_cmm] [--not_write_xyz] [--not_write_json]
    
    Generates 3D models given an input interaction matrix and a set of input parameters
    
    optional arguments:
      -h, --help               show this help message and exit
    
    General options:
      -w PATH, --workdir PATH  path to working directory (generated with the tool TADbit mapper)
      --input_matrix PATH      In case input was not generated with the TADbit tools
      --rand INT               [1] random initial number. NOTE: when running single model at the
                               time, should be different for each run
      --nmodels INT            [5000] number of models to generate for modeling
      --nkeep INT              [1000] number of models to keep for modeling
      -j INT, --jobid INT      Use as input data generated by a job with a given jobid. Use tadbit
                               describe to find out which.
      --optimization_id INT    [None] ID of a pre-run optimization batch job
      --restart_id INT         [None] ID of a job to be restarted, for example after building the
                               models in a cluster
      --fig_format STR         file format and extension for figures and plots (can be any
                               supported by matplotlib, png, eps...)
      --noX                    no display server (X screen)
      --corr STR               correlation method to compare contact maps and original matrix
                               (options are speraman, pearson, kendall, logpearson, chi2, scc )
    
    Descriptive, optional arguments:
      --species STRING         species name, with no spaces, i.e.: homo_sapiens
      --assembly STRING        NCBI ID of the original assembly (i.e.: NCBI36 for human)
      --cell STRING            cell type name
      --exp_type STRING        experiment type name (i.e.: Hi-C)
      --project STRING         project name
    
    Modeling preparation:
      --crm NAME               chromosome name
      --beg INT                genomic coordinate from which to start modeling
      --end INT                genomic coordinate where to end modeling
      --matrix_beg INT         genomic coordinate of the first row/column of the input matrix. This
                               has to be specified if the input matrix is not the TADbit tools
                               generated abc format
      -r INT, --reso INT       resolution of the Hi-C experiment
      --perc_zero FLOAT
    
    Parameter optimization:
      --optimize               optimization run, store less info about models
      --model                  modelling run
      --model_ptadbit          modelling run using pTADbit
      --force                  use input parameters, and skip any precalculated optimization
      --maxdist LIST [LIST ...]
                               range of numbers for maxdist, i.e. 400:1000:100 -- or just a number
                               -- or a list of numbers
      --upfreq LIST [LIST ...]
                               range of numbers for upfreq, i.e. 0:1.2:0.3 -- or just a number --
                               or a list of numbers
      --lowfreq LIST [LIST ...]
                               range of numbers for lowfreq, i.e. -1.2:0:0.3 -- or just a number --
                               or a list of numbers
      --scale LIST [LIST ...]  [0.01] range of numbers to be test as optimal scale value, i.e.
                               0.005:0.01:0.001 -- Can also pass only one number -- or a list of
                               numbers
      --dcutoff LIST [LIST ...]
                               [2] range of numbers to be test as optimal distance cutoff parameter
                               (distance, in number of beads, from which to consider 2 beads as
                               being close), i.e. 1:1.5:0.5 -- Can also pass only one number -- or
                               a list of numbers
      --container LIST [LIST ...]
                               restrains particle to be within a given object. Can only be a
                               'cylinder', which is, in fact a cylinder of a given height to which
                               are added hemispherical ends. This cylinder is defined by a radius,
                               its height (with a height of 0 the cylinder becomes a sphere) and
                               the force applied to the restraint. E.g. for modeling E. coli genome
                               (2 micrometers length and 0.5 micrometer of width), these values
                               could be used: 'cylinder' 250 1500 50, and for a typical mammalian
                               nuclei (6 micrometers diameter): 'cylinder' 3000 0 50
      --analyze                analyze models.
    
    Analysis:
      --analyze_list INT [INT ...]
                               [2 3 4 5 6 7 8 9 10 11 12 13] list of numbers representing the
                               analysis to be done. Choose between: 0) do nothing 1) optimization
                               plot 2) correlation real/models 3) z-score plot 4) constraints 5)
                               objective function 6) centroid 7) consistency 8) density 9) contact
                               map 10) walking angle 11) persistence length 12) accessibility 13)
                               interaction
      --not_write_cmm          [False] do not generate cmm files for each model (Chimera input)
      --not_write_xyz          [False] do not generate xyz files for each model (3D coordinates)
      --not_write_json         [False] do not generate json file.
    
    Running jobs:
      --smooth_factor INT      Hi-C matrix smoothing value of the mean kernel for pTADbit. Useful
                               in case of using matrices with low sequencing depth
      -C CPUS, --cpu CPUS      [32] Maximum number of CPU cores available in the execution host. If
                               higher than 1, tasks with multi-threading capabilities will enabled
                               (if 0 all available) cores will be used
      --job_list               generate a list of commands stored in a file named joblist_HASH.q
                               (where HASH is replaced by a string specific to the parameters
                               used). note that dcutoff will never be split as it does not require
                               to re-run models.
      --nmodels_per_job INT    Number of models per distributed job.
      --cpus_per_job INT       Number of cpu nodes per distributed job.
      --concurrent_jobs INT    Number of concurrent jobs in distributed mode.
      --timeout_job INT        Time to wait for a concurrent jobs to finish before canceling it in
                               distributed mode.
      --script_cmd STR         Command to call the jobs in distributed mode.
      --script_args STR        Argumnets to script_cmd to call the jobs in distributed mode.
      --script_template STR    Template to generate a file that script_cmd will call for each job
                               in distributed mode. Each __file__ marker in the template will be
                               replacedby the job file __name__ with the name and __dir__ with the
                               folder.
      --tmpdb PATH             if provided uses this directory to manipulate the database

TADbit merge
------------

.. code:: ipython3

    Load two working directories with different Hi-C data samples and merges them into a new
    working directory generating some statistics.
    
    usage: tadbit merge [-h] [-w PATH] [-w1 PATH] [-w2 PATH] [--bam1 PATH] [--noX] [--bam2 PATH]
                        [-C CPUS] [-r INT] [--skip_comparison] [--skip_merge]
                        [--save STR [STR ...]] [--jobid1 INT] [--jobid2 INT] [--force] [--norm]
                        [--biases1 PATH] [--biases2 PATH] [--filter INT [INT ...]]
                        [--samtools PATH] [--tmpdb PATH]
    
    optional arguments:
      -h, --help               show this help message and exit
    
    General options:
      -w PATH, --workdir PATH  path to a new output folder
      -w1 PATH, --workdir1 PATH
                               path to working directory of the first HiC data sample to merge
      -w2 PATH, --workdir2 PATH
                               path to working directory of the second HiC data sample to merge
      --bam1 PATH              path to the first TADbit-generated BAM file with all reads (other
                               wise the tool will guess from the working directory database)
      --noX                    no display server (X screen)
      --bam2 PATH              path to the second TADbit-generated BAM file with all reads (other
                               wise the tool will guess from the working directory database)
      -C CPUS, --cpus CPUS     [32] Maximum number of CPU cores available in the execution host. If
                               higher than 1, tasks with multi-threading capabilities will enabled
                               (if 0 all available) cores will be used
      -r INT, --resolution INT
                               resolution at which to do the comparison, and generate the matrices.
      --skip_comparison        skip the comparison between replicates (faster). Comparisons are
                               performed at 3 levels 1- comparing first diagonals of each
                               experiment (and generating SCC score and standard deviation see
                               https://doi.org/10.1101/gr.220640.117) 2- Comparing the first
                               eigenvectors of input experiments 3- Generates reproducibility score
                               using function from https://doi.org/10.1093/bioinformatics/btx152
      --skip_merge             skip the merge of replicates (faster).
      --save STR [STR ...]     [genome] save genomic or chromosomic matrix.
      --jobid1 INT             Use as input data generated by a job with a given jobid. Use tadbit
                               describe to find out which.
      --jobid2 INT             Use as input data generated by a job with a given jobid. Use tadbit
                               describe to find out which.
      --force                  overwrite previously run job
      --norm                   compare normalized matrices
      --biases1 PATH           path to file with precalculated biases by columns
      --biases2 PATH           path to file with precalculated biases by columns
      --filter INT [INT ...]   [[1, 2, 3, 4, 6, 7, 9, 10]] Use filters to define a set os valid
                               pair of reads e.g.: '--apply 1 2 3 4 8 9 10'. Where these
                               numberscorrespond to: 1: self-circle, 2: dangling-end, 3: error, 4:
                               extra dangling-end, 5: too close from RES, 6: too short, 7: too
                               large, 8: over-represented, 9: duplicated, 10: random breaks, 11:
                               trans-chromosomic
      --samtools PATH          path samtools binary
      --tmpdb PATH             if provided uses this directory to manipulate the database

TADbit describe
---------------

.. code:: ipython3

    Describe jobs and results in a given working directory
    
    usage: tadbit describe [-h] [-w PATH] [--noX] [-t  [...]] [-T  [...]] [-j INT [INT ...]]
                           [-W STR [STR ...]] [-s STR [STR ...]] [--tmpdb PATH] [--tsv] [-o OUTPUT]
    
    optional arguments:
      -h, --help               show this help message and exit
    
    General options:
      -w PATH, --workdir PATH  path to working directory (generated with the tool tadbit map)
      --noX                    no display server (X screen)
      -t  [ ...], --tables  [ ...]
                               [['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
                               '13']] what tables to show, write either the sequence of names or
                               indexes, according to this list: 1: paths, 2: jobs, 3:
                               mapped_outputs, 4: mapped_inputs, 5: parsed_outputs, 6:
                               intersection_outputs, 7: filter_outputs, 8: normalize_outputs, 9:
                               merge_stats, 10: merge_outputs, 11: segment_outputs, 12: models, 13:
                               modeled_regions
      -T  [ ...], --skip_tables  [ ...]
                               [[]] what tables NOT to show, write either the sequence of names or
                               indexes, according to this list: 1: paths, 2: jobs, 3:
                               mapped_outputs, 4: mapped_inputs, 5: parsed_outputs, 6:
                               intersection_outputs, 7: filter_outputs, 8: normalize_outputs, 9:
                               merge_stats, 10: merge_outputs, 11: segment_outputs, 12: models, 13:
                               modeled_regions
      -j INT [INT ...], --jobids INT [INT ...]
                               Display only items matching these jobids.
      -W STR [STR ...], --where STR [STR ...]
                               Select rows. List pairs of keywords (column header) and values to
                               filter results. For example to get only results where "18" appears
                               in the column "Chromosome", the option should be set as: `-W
                               Chromosome,18`
      -s STR [STR ...], --select STR [STR ...]
                               Select columns. List the keyword (column headers) to be displayed.
                               E.g. to show only the colmns JobIds: `-s Jobids`
      --tmpdb PATH             if provided uses this directory to manipulate the database
      --tsv                    Print output in tab separated format
      -o OUTPUT, --output OUTPUT
                               Writes output in specified file.

TADbit clean
------------

.. code:: ipython3

    Delete jobs and results of a given list of jobids in a given directories
    
    usage: tadbit clean [-h] [-w PATH] [-j INT [INT ...]] [--delete] [--compress] [--noX]
                        [--change_workdir PATH] [--tmpdb PATH]
    
    optional arguments:
      -h, --help               show this help message and exit
      --change_workdir PATH    In case folder was moved, input the new path
    
    General options:
      -w PATH, --workdir PATH  path to working directory (generated with the tool tadbit mapper)
      -j INT [INT ...], --jobids INT [INT ...]
                               jobids of the files and entries to be removed
      --delete                 delete files, otherwise only DB entries.
      --compress               compress files and update paths accordingly
      --noX                    no display server (X screen)
      --tmpdb PATH             if provided uses this directory to manipulate the database

TADbit import
-------------

.. code:: ipython3

    Import Hi-C data to TADbit toy BAM
    
    usage: tadbit import [-h] -w PATH -r INT [--format {text,matrix,cooler}] -i STR [-c]
                         [--tmpdb PATH] [-C CPUS] [--samtools PATH]
    
    optional arguments:
      -h, --help               show this help message and exit
    
    Required options:
      -w PATH, --workdir PATH  path to working directory (generated with the tool tadbit mapper)
      -r INT, --resolution INT
                               resolution at which to output matrices
      --format {text,matrix,cooler}
                               [text] can be any of [text, matrix, cooler]
      -i STR, --input STR      path to input file
    
    General options:
      -c , --coord             Coordinate of the region to import. By default all genome, arguments
                               can be either one chromosome name, or the coordinate in the form:
                               "-c chr3:110000000-120000000"
      --tmpdb PATH             if provided uses this directory to manipulate the database
      -C CPUS, --cpus CPUS     [32] Maximum number of CPU cores available in the execution host. If
                               higher than 1, tasks with multi-threading capabilities will enabled
                               (if 0 all available) cores will be used
      --samtools PATH          path samtools binary

TADbit export
-------------

.. code:: ipython3

    Export Hi-C data to other formats
    
    usage: tadbit export [-h] -w PATH -r INT [--format {text,matrix,cooler,hic}] -o STR
                         [--bam PATH] [-j INT] [--force] [-q] [--tmpdb PATH] [--nchunks NCHUNKS]
                         [-C CPUS] [--chr_name STR [STR ...]] [--juicerjar PATH] [--rownames] [-c]
                         [-c2] [--biases PATH] [--norm] [-F INT [INT ...]]
    
    optional arguments:
      -h, --help               show this help message and exit
    
    Required options:
      -w PATH, --workdir PATH  path to working directory (generated with the tool tadbit mapper)
      -r INT, --resolution INT
                               resolution at which to output matrices
      --format {text,matrix,cooler,hic}
                               [text] can be any of [text, matrix, cooler, hic]
      -o STR, --output STR     path to output file
    
    General options:
      --bam PATH               path to a TADbit-generated BAM file with all reads (other wise the
                               tool will guess from the working directory database)
      -j INT, --jobid INT      Use as input data generated by a job with a given jobid. Use tadbit
                               describe to find out which.
      --force                  overwrite previously run job
      -q, --quiet              remove all messages
      --tmpdb PATH             if provided uses this directory to manipulate the database
      --nchunks NCHUNKS        maximum number of chunks into which to cut the BAM
      -C CPUS, --cpus CPUS     [32] Maximum number of CPU cores available in the execution host. If
                               higher than 1, tasks with multi-threading capabilities will enabled
                               (if 0 all available) cores will be used
      --chr_name STR [STR ...]
                               [fasta header] chromosome name(s). Order of chromosomes in the
                               output matrices.
      --juicerjar PATH         path to the juicer tools jar file needed to export matrices to hic
                               format (check https://github.com/aidenlab/juicer/wiki/Download).
                               Note that you also need java available in the path.
    
    Read filtering options:
      -F INT [INT ...], --filter INT [INT ...]
                               [[1, 2, 3, 4, 6, 7, 9, 10]] Use filters to define a set of valid
                               pair of reads e.g.: '--filter 1 2 3 4 8 9 10'. Where these numbers
                               correspond to: 0: nothing, 1: self-circle, 2: dangling-end, 3:
                               error, 4: extra dangling-end, 5: too close from RES, 6: too short,
                               7: too large, 8: over-represented, 9: duplicated, 10: random breaks,
                               11: trans-chromosomic
    
    Normalization options:
      --biases PATH            path to file with pre-calculated biases by columns
      --norm                   export normalized matrix
    
    Output options:
      --rownames               To store row names in the output text matrix. WARNING: when non-
                               matrix, results in two extra columns
      -c , --coord             Coordinate of the region to retrieve. By default all genome,
                               arguments can be either one chromosome name, or the coordinate in
                               the form: "-c chr3:110000000-120000000"
      -c2 , --coord2           Coordinate of a second region to retrieve the matrix in the
                               intersection with the first region.
