Retrieve HiC dataset from NCBI
==============================

We will use data from \ `(Stadhouders R, Vidal E, Serra F, Di Stefano B
et al. 2018) <#cite-ralph>`__, which comes from mouse cells where Hi-C
experiment where conducted in different states during highly-efficient
somatic cell reprogramming.

The data can be downloaded from:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53463

Once downloaded the files can be converted to the FASTQ format in order
for TADbit to read them.

The easiest way to download the data might be through the ``fastq-dump``
program from the SRA Toolkit
(http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software).

We download 100M reads for each of 4 replicates (2 replicates from B
cells and 2 from Pluripotent Stem Cells),and organize each in two files,
one per read-end (this step is long and can take **up to 6 hours**):

.. code:: bash

    %%bash
    
    mkdir -p FASTQs
    
    fastq-dump -A SRR5344921 -DQ '+' --defline-seq '@$ac.$si' -X 100000000 --split-files --outdir FASTQs/
    mv FASTQs/SRR5344921_1.fastq FASTQs/mouse_B_rep1_1.fastq
    mv FASTQs/SRR5344921_2.fastq FASTQs/mouse_B_rep1_2.fastq
    
    fastq-dump -A SRR5344925 -DQ '+' --defline-seq '@$ac.$si' -X 100000000 --split-files --outdir FASTQs/
    mv FASTQs/SRR5344925_1.fastq FASTQs/mouse_B_rep2_1.fastq
    mv FASTQs/SRR5344925_2.fastq FASTQs/mouse_B_rep2_2.fastq
    
    fastq-dump -A SRR5344969 -DQ '+' --defline-seq '@$ac.$si' -X 100000000 --split-files --outdir FASTQs
    mv FASTQs/SRR5344969_1.fastq FASTQs/mouse_PSC_rep1_1.fastq
    mv FASTQs/SRR5344969_2.fastq FASTQs/mouse_PSC_rep1_2.fastq
    
    fastq-dump -A SRR5344973 -DQ '+' --defline-seq '@$ac.$si' -X 100000000 --split-files --outdir FASTQs/
    mv FASTQs/SRR5344973_1.fastq FASTQs/mouse_PSC_rep2_1.fastq
    mv FASTQs/SRR5344973_2.fastq FASTQs/mouse_PSC_rep2_2.fastq


.. ansi-block::

    Read 100000000 spots for SRR5344921
    Written 100000000 spots for SRR5344921
    Read 100000000 spots for SRR5344925
    Written 100000000 spots for SRR5344925
    Read 100000000 spots for SRR5344969
    Written 100000000 spots for SRR5344969
    Read 100000000 spots for SRR5344973
    Written 100000000 spots for SRR5344973


Files are renamed for convenience.

*Note: the parameter used here for fastq-dump are for generating simple
FASTQ files, -DQ ‘+’ tells to put a single plus sign as header of the
quality input, ``--defline-seq ‘@$ac.$si’`` reduces the information in
the headers to the accession number and the read id, ``--split-files``
is to separate both read-ends in different files, finally
``-X 100000000`` is to download only the first 100 Million reads of each
replicate*

*Note: alternatively you can also directly download the FASTQ from
http://www.ebi.ac.uk/*

Compression
~~~~~~~~~~~

Each of these 8 files, contains 100M reads of 75 nucleotides each, and
occupies ~17 Gb (total 130 Gb).

Internally we use DSRC \ `(Roguski and Deorowicz,
2014) <#cite-roguski2014dsrc>`__ that allows better compression ration
and, more importantly, faster decompression:

.. code:: bash

    %%bash
    
    dsrc c -t8 FASTQs/mouse_B_rep1_1.fastq FASTQs/mouse_B_rep1_1.fastq.dsrc
    dsrc c -t8 FASTQs/mouse_B_rep1_2.fastq FASTQs/mouse_B_rep1_2.fastq.dsrc
    dsrc c -t8 FASTQs/mouse_B_rep2_1.fastq FASTQs/mouse_B_rep2_1.fastq.dsrc
    dsrc c -t8 FASTQs/mouse_B_rep2_2.fastq FASTQs/mouse_B_rep2_2.fastq.dsrc
    dsrc c -t8 FASTQs/mouse_PSC_rep1_1.fastq FASTQs/mouse_PSC_rep1_1.fastq.dsrc
    dsrc c -t8 FASTQs/mouse_PSC_rep1_2.fastq FASTQs/mouse_PSC_rep1_2.fastq.dsrc
    dsrc c -t8 FASTQs/mouse_PSC_rep2_1.fastq FASTQs/mouse_PSC_rep2_1.fastq.dsrc
    dsrc c -t8 FASTQs/mouse_PSC_rep2_2.fastq FASTQs/mouse_PSC_rep2_2.fastq.dsrc

After compression we reduce the total size to 27 Gb (**20% of the
original size, and dsrc ensures fast reading of the compressed data**)

*Note:* - *using gzip instead reduces size to ~38 Gb (occupies ~40% more
than dsrc compressed files)* - *using bzip2 instead reduces size to ~31
Gb (occupies ~15% more than dsrc compressed files)*

*Both are much slower to generate and read*

Cleanup
~~~~~~~

.. code:: bash

    %%bash
    
    rm -f FASTQs/mouse_B_rep1_1.fastq
    rm -f FASTQs/mouse_B_rep1_2.fastq
    rm -f FASTQs/mouse_B_rep2_1.fastq
    rm -f FASTQs/mouse_B_rep2_2.fastq
    rm -f FASTQs/mouse_PSC_rep1_1.fastq
    rm -f FASTQs/mouse_PSC_rep1_2.fastq
    rm -f FASTQs/mouse_PSC_rep2_1.fastq
    rm -f FASTQs/mouse_PSC_rep2_2.fastq

References
~~~~~~~~~~

[^](#ref-1) Stadhouders R, Vidal E, Serra F, Di Stefano B et al. 2018.
*Transcription factors orchestrate dynamic interplay between genome
topology and gene regulation during cell reprogramming*.

[^](#ref-4) Roguski, :raw-latex:`\Lukasz `and Deorowicz, Sebastian.
2014. *DSRC 2—Industry-oriented compression of FASTQ files*.
