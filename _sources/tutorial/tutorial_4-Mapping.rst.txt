
Iterative vs fragment-based mapping
===================================

Iterative mapping first proposed by \ `(Imakaev et al.,
2012) <#cite-Imakaev2012a>`__, allows to map usually a high number of
reads. However other methodologies, less "brute-force" can be used to
take into account the chimeric nature of the Hi-C reads.

A simple alternative is to allow split mapping.

Another way consists in *pre-truncating* \ `(Ay and Noble,
2015) <#cite-Ay2015a>`__ reads that contain a ligation site and map only
the longest part of the read \ `(Wingett et al.,
2015) <#cite-Wingett2015>`__.

Finally, an intermediate approach, *fragment-based*, consists in mapping
full length reads first, and than splitting unmapped reads at the
ligation sites \ `(Serra et al. 2017) <#cite-Serra2017>`__.

.. figure:: ../nbpictures/mapping.png
   :alt: mapping strategies

   mapping strategies

Advantages of iterative mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  It's the only solution when no restriction enzyme has been used (i.e.
   micro-C)
-  Can be faster when few windows (2 or 3) are used

Advantages of fragment-based mapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Generally faster
-  Safer: mapped reads are generally larger than 25-30 nm (the largest
   window used in iterative mapping). Less reads are mapped, but the
   difference is usually canceled or reversed when looking for
   "valid-pairs".

*Note:* We use **GEM2** \ `(Marco-Sola et al.
2012) <#cite-Marco-Sola2012>`__, performance are very similar to
Bowtie2, and in some cases slightly better.

*For now TADbit is only compatible with GEM2.*

Mapping
=======

.. code:: ipython2

    from pytadbit.mapping.full_mapper import full_mapping

The full mapping function can be used to perform either iterative or
fragment-based mapping, or a combination of both.

**It is important to note** that although the default mapping parameters
used by TADbit are relatively strict, **a non negligible proportion of
the reads will be mis-mapped**, and this applies **at each iteration of
the mapping.**

Iterative mapping
-----------------

Here an example of use as iterative mapping: (Estimated time 15h with 8
cores)

.. code:: ipython2

    cell = 'mouse_B'  # or mouse_PSC
    rep = 'rep1'  # or rep2

*Note: the execution of this notebook should be repeated for each of the
4 replicates*

.. code:: ipython2

    ! mkdir -p results/iterativ/$cell\_$rep
    ! mkdir -p results/iterativ/$cell\_$rep/01_mapping

.. code:: ipython2

    # for the first side of the reads
    full_mapping(gem_index_path='genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem',
                 out_map_dir='results/iterativ/{0}_{1}/01_mapping/mapped_{0}_{1}_r1/'.format(cell, rep),
                 fastq_path='FASTQs/%s_%s_1.fastq.dsrc' % (cell,rep),
                 frag_map=False,  clean=True, nthreads=8,
                 windows=((1,25),(1,35),(1,45),(1,55),(1,65),(1,75)),
                 temp_dir='results/iterativ/{0}_{1}/01_mapping/mapped_{0}_{1}_r1_tmp/'.format(cell, rep))


.. ansi-block::

    Preparing FASTQ file
      - conversion to MAP format
      - trimming reads 1-25
    Mapping reads in window 1-25...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_PSC_rep1/01_mapping/mapped_mouse_PSC_rep1_r1_tmp/mouse_PSC_rep1_1.fastq_khgGNE
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_PSC_rep1/01_mapping/mapped_mouse_PSC_rep1_r1_tmp/mouse_PSC_rep1_1.fastq_khgGNE -o /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_PSC_rep1/01_mapping/mapped_mouse_PSC_rep1_r1_tmp/mouse_PSC_rep1_1.fastq_khgGNE_full_1-25
    Parsing result...


And for the second side of the read-end:

.. code:: ipython2

    # for the second side of the reads
    full_mapping(gem_index_path='genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem',
                 out_map_dir='results/iterativ/{0}_{1}/01_mapping/mapped_{0}_{1}_r2/'.format(cell, rep),
                 fastq_path='FASTQs/%s_%s_2.fastq.dsrc' % (cell,rep),
                 frag_map=False,  clean=True, nthreads=8,
                 windows=((1,25),(1,35),(1,45),(1,55),(1,65),(1,75)), 
                 temp_dir='results/iterativ/{0}_{1}/01_mapping/mapped_{0}_{1}_r2_tmp/'.format(cell, rep))


.. ansi-block::

    Preparing FASTQ file
      - conversion to MAP format
      - trimming reads 1-25
    Mapping reads in window 1-25...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_rErDYK
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_rErDYK -o /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_rErDYK_full_1-25
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_rErDYK
       x removing map /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_rErDYK_full_1-25.map
    Preparing MAP file
      - trimming reads 1-35
       x removing original input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_rErDYK_filt_1-25.map
    Mapping reads in window 1-35...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_14i9tR
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_14i9tR -o /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_14i9tR_full_1-35
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_14i9tR
       x removing map /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_14i9tR_full_1-35.map
    Preparing MAP file
      - trimming reads 1-45
       x removing original input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_14i9tR_filt_1-35.map
    Mapping reads in window 1-45...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_xiEdeg
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_xiEdeg -o /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_xiEdeg_full_1-45
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_xiEdeg
       x removing map /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_xiEdeg_full_1-45.map
    Preparing MAP file
      - trimming reads 1-55
       x removing original input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_xiEdeg_filt_1-45.map
    Mapping reads in window 1-55...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_EYDFyI
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_EYDFyI -o /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_EYDFyI_full_1-55
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_EYDFyI
       x removing map /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_EYDFyI_full_1-55.map
    Preparing MAP file
      - trimming reads 1-65
       x removing original input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_EYDFyI_filt_1-55.map
    Mapping reads in window 1-65...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_pe6yQw
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_pe6yQw -o /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_pe6yQw_full_1-65
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_pe6yQw
       x removing map /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_pe6yQw_full_1-65.map
    Preparing MAP file
      - trimming reads 1-75
       x removing original input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_pe6yQw_filt_1-65.map
    Mapping reads in window 1-75...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_MkfcZ5
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_MkfcZ5 -o /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_MkfcZ5_full_1-75
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_MkfcZ5
       x removing map /scratch/workspace/MethodsMolBiol/Notebooks/results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_MkfcZ5_full_1-75.map




.. ansi-block::

    ['results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2/mouse_B_rep1_2.fastq_full_1-25.map',
     'results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2/mouse_B_rep1_2.fastq_full_1-35.map',
     'results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2/mouse_B_rep1_2.fastq_full_1-45.map',
     'results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2/mouse_B_rep1_2.fastq_full_1-55.map',
     'results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2/mouse_B_rep1_2.fastq_full_1-65.map',
     'results/iterativ/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2/mouse_B_rep1_2.fastq_full_1-75.map']



Fragment-based mapping
----------------------

The fragment-based mapping strategy works in 2 steps: 1. The read-ends
are mapped entirely, assuming that no ligation occurred in them. 2. For
unmapped read-ends, the function searches for a ligation site (e.g. in
the case of MboI this would correspond to ``GATCGATC`` and in the case
of HindIII to ``AAGCTAGCTT``). The read-end is split accordingly
replacing the ligation site by two RE sites:

::

                            read-end-part-one---AAGCTAGCTT----read-end-part-two
         will be split in:
                            read-end-part-one---AAGCTT
         and:
                                                    AAGCTT----read-end-part-two

*Note: **if no ligation site is found**, step two will be repeated using
digested RE site as split point (``AAGCT`` in the case of HindIII). This
is done in order to be protected against sequencing errors. When this
path is followed the digested RE site is removed, but not replaced. **If
digested RE sites are not found either, the read will be classified as
unmapped**.*

*Note: **both mapping strategies can be combined**, for example defining
the windows as previously (iterative mapping), but also give a RE name*
``r_enz=MboI`` *and setting* ``frag_map=True`` *like this if a read has
not been mapped in any window, TADbit will also try to apply the
fragment-based strategy.*

.. code:: ipython2

    ! mkdir -p results/fragment/$cell\_$rep
    ! mkdir -p results/fragment/$cell\_$rep/01_mapping

.. code:: ipython2

    # for the first side of the reads 
    full_mapping(gem_index_path='genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem',
                 out_map_dir='results/fragment/{0}_{1}/01_mapping/mapped_{0}_{1}_r1/'.format(cell, rep),
                 fastq_path='FASTQs/%s_%s_1.fastq.dsrc' % (cell, rep),
                 r_enz='MboI', frag_map=True, clean=True, nthreads=8, 
                 temp_dir='results/fragment/{0}_{1}/01_mapping/mapped_{0}_{1}_r1_tmp/'.format(cell, rep))


.. ansi-block::

    Preparing FASTQ file
      - conversion to MAP format
    Mapping reads in window 1-end...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_Q50vCN
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_Q50vCN -o /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_Q50vCN_full_1-end
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_Q50vCN
       x removing map /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_Q50vCN_full_1-end.map
      - splitting into restriction enzyme (RE) fragments using ligation sites
      - ligation sites are replaced by RE sites to match the reference genome
        * enzymes: MboI & MboI, ligation site: GATCGATC, RE site: GATC & GATC
    Preparing MAP file
       x removing pre-GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_Q50vCN_filt_1-end.map
    Mapping fragments of remaining reads...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_lYQC9s
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_lYQC9s -o /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_lYQC9s_frag_1-end
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_lYQC9s
       x removing failed to map /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_Q50vCN_fail.map
       x removing tmp mapped /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1_tmp/mouse_B_rep1_1.fastq_lYQC9s_frag_1-end.map




.. ansi-block::

    ['results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1/mouse_B_rep1_1.fastq_full_1-end.map',
     'results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r1/mouse_B_rep1_1.fastq_frag_1-end.map']



.. code:: ipython2

    # for the second side of the reads
    full_mapping(gem_index_path='genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem',
                 out_map_dir='results/fragment/{0}_{1}/01_mapping/mapped_{0}_{1}_r2/'.format(cell, rep),
                 fastq_path='FASTQs/%s_%s_2.fastq.dsrc' % (cell, rep),
                 r_enz='MboI', frag_map=True, clean=True, nthreads=8, 
                 temp_dir='results/fragment/{0}_{1}/01_mapping/mapped_{0}_{1}_r2_tmp/'.format(cell, rep))


.. ansi-block::

    Preparing FASTQ file
      - conversion to MAP format
    Mapping reads in window 1-end...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_aAwRcj
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_aAwRcj -o /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_aAwRcj_full_1-end
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_aAwRcj
       x removing map /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_aAwRcj_full_1-end.map
      - splitting into restriction enzyme (RE) fragments using ligation sites
      - ligation sites are replaced by RE sites to match the reference genome
        * enzymes: MboI & MboI, ligation site: GATCGATC, RE site: GATC & GATC
    Preparing MAP file
       x removing pre-GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_aAwRcj_filt_1-end.map
    Mapping fragments of remaining reads...
    TO GEM /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_V5kpMT
    /home/dcastillo/miniconda2/bin/gem-mapper -I /scratch/workspace/MethodsMolBiol/Notebooks/genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6_contigs.gem -q offset-33 -m 0.04 -s 0 --allow-incomplete-strata 0.00 --granularity 10000 --max-decoded-matches 1 --min-decoded-strata 0 --min-insert-size 0 --max-insert-size 0 --min-matched-bases 0.8 --gem-quality-threshold 26 --max-big-indel-length 15 --mismatch-alphabet ACGT -E 0.30 --max-extendable-matches 20 --max-extensions-per-match 1 -e 0.04 -T 8 -i /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_V5kpMT -o /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_V5kpMT_frag_1-end
    Parsing result...
       x removing GEM input /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_V5kpMT
       x removing failed to map /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_aAwRcj_fail.map
       x removing tmp mapped /scratch/workspace/MethodsMolBiol/Notebooks/results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2_tmp/mouse_B_rep1_2.fastq_V5kpMT_frag_1-end.map




.. ansi-block::

    ['results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2/mouse_B_rep1_2.fastq_full_1-end.map',
     'results/fragment/mouse_B_rep1/01_mapping/mapped_mouse_B_rep1_r2/mouse_B_rep1_2.fastq_frag_1-end.map']



.. raw:: html

   <!--bibtex
   @article{Imakaev2012a,
   abstract = {Extracting biologically meaningful information from chromosomal interactions obtained with genome-wide chromosome conformation capture (3C) analyses requires the elimination of systematic biases. We present a computational pipeline that integrates a strategy to map sequencing reads with a data-driven method for iterative correction of biases, yielding genome-wide maps of relative contact probabilities. We validate this ICE (iterative correction and eigenvector decomposition) technique on published data obtained by the high-throughput 3C method Hi-C, and we demonstrate that eigenvector decomposition of the obtained maps provides insights into local chromatin states, global patterns of chromosomal interactions, and the conserved organization of human and mouse chromosomes.},
   author = {Imakaev, Maxim V and Fudenberg, Geoffrey and McCord, Rachel Patton and Naumova, Natalia and Goloborodko, Anton and Lajoie, Bryan R and Dekker, Job and Mirny, Leonid A},
   doi = {10.1038/nmeth.2148},
   file = {:home/fransua/.local/share/data/Mendeley Ltd./Mendeley Desktop/Downloaded/Imakaev et al. - 2012 - Iterative correction of Hi-C data reveals hallmarks of chromosome organization.pdf:pdf},
   issn = {1548-7105},
   journal = {Nature methods},
   keywords = {Hi-C},
   mendeley-groups = {stats/Hi-C,Research articles},
   mendeley-tags = {Hi-C},
   month = {oct},
   number = {10},
   pages = {999--1003},
   pmid = {22941365},
   title = {{Iterative correction of Hi-C data reveals hallmarks of chromosome organization.}},
   url = {http://www.ncbi.nlm.nih.gov/pubmed/22941365},
   volume = {9},
   year = {2012}
   }
   @article{Ay2015a,
   author = {Ay, Ferhat and Noble, William Stafford},
   doi = {10.1186/s13059-015-0745-7},
   file = {:home/fransua/.local/share/data/Mendeley Ltd./Mendeley Desktop/Downloaded/Ay, Noble - 2015 - Analysis methods for studying the 3D architecture of the genome.pdf:pdf},
   issn = {1474-760X},
   journal = {Genome Biology},
   keywords = {Chromatin conformation capture,Genome architecture,Hi-C,Three-dimensional genome,Three-dimensional modeling,chromatin,conformation capture,genome architecture,in many other applications,ranging from genome assem-,review,three-dimensional genome,three-dimensional modeling},
   mendeley-groups = {Research articles},
   mendeley-tags = {Hi-C,review},
   number = {1},
   pages = {183},
   publisher = {Genome Biology},
   title = {{Analysis methods for studying the 3D architecture of the genome}},
   url = {http://genomebiology.com/2015/16/1/183},
   volume = {16},
   year = {2015}
   }
   @article{Wingett2015,
   abstract = {HiCUP is a pipeline for processing sequence data generated by Hi-C and Capture Hi-C (CHi-C) experiments, which are techniques used to investigate three-dimensional genomic organisation. The pipeline maps data to a specified reference genome and removes artefacts that would otherwise hinder subsequent analysis. HiCUP also produces an easy-to-interpret yet detailed quality control (QC) report that assists in refining experimental protocols for future studies. The software is freely available and has already been used for processing Hi-C and CHi-C data in several recently published peer-reviewed studies.},
   author = {Wingett, Steven and Ewels, Philip and Furlan-Magaril, Mayra and Nagano, Takashi and Schoenfelder, Stefan and Fraser, Peter and Andrews, Simon},
   doi = {10.12688/f1000research.7334.1},
   file = {:home/fransua/Downloads/f1000research-4-7903.pdf:pdf},
   issn = {2046-1402},
   journal = {F1000Research},
   mendeley-groups = {Computer programs/Hi-C/Hi-C processing},
   pages = {1310},
   pmid = {26835000},
   title = {{HiCUP: pipeline for mapping and processing Hi-C data.}},
   url = {http://f1000research.com/articles/4-1310/v1},
   volume = {4},
   year = {2015}
   }
   @article{Serra2016,
   abstract = {The sequence of a genome is insufficient to understand all genomic processes carried out in the cell nucleus. To achieve this, the knowledge of its three- dimensional architecture is necessary. Advances in genomic technologies and the development of new analytical methods, such as Chromosome Conformation Capture (3C) and its derivatives, now permit to investigate the spatial organization of genomes. However, inferring structures from raw contact data is a tedious process for shortage of available tools. Here we present TADbit, a computational framework to analyze and model the chromatin fiber in three dimensions. To illustrate the use of TADbit, we automatically modeled 50 genomic domains from the fly genome revealing differential structural features of the previously defined chromatin colors, establishing a link between the conformation of the genome and the local chromatin composition. More generally, TADbit allows to obtain three-dimensional models ready for visualization from 3C-based experiments and to characterize their relation to gene expression and epigenetic states. TADbit is open-source and available for download from http://www.3DGenomes.org.},
   author = {Serra, Fran{\c{c}}ois and Ba{\`{u}}, Davide and Filion, Guillaume and Marti-Renom, Marc A.},
   doi = {10.1101/036764},
   file = {:home/fransua/.local/share/data/Mendeley Ltd./Mendeley Desktop/Downloaded/Serra et al. - 2016 - Structural features of the fly chromatin colors revealed by automatic three-dimensional modeling.pdf:pdf},
   journal = {bioRxiv},
   keywords = {3d,Hi-C,capture,genome architecture,genome reconstruction chromosome conformation,modeling,optimization,resampling,restraint-based,three-dimensional genome reconstruction},
   mendeley-groups = {Research articles,projects/WIREs{\_}review/modeling{\_}tools,Computer programs/Hi-C/Hi-C modeling},
   mendeley-tags = {Hi-C,modeling,optimization,resampling},
   pages = {1--29},
   title = {{Structural features of the fly chromatin colors revealed by automatic three-dimensional modeling.}},
   url = {http://biorxiv.org/content/early/2016/01/15/036764},
   year = {2016}
   }
   @misc{Marco-Sola2012,
   abstract = {Because of ever-increasing throughput requirements of sequencing data, most existing short-read aligners have been designed to focus on speed at the expense of accuracy. The Genome Multitool (GEM) mapper can leverage string matching by filtration to search the alignment space more efficiently, simultaneously delivering precision (performing fully tunable exhaustive searches that return all existing matches, including gapped ones) and speed (being several times faster than comparable state-of-the-art tools).},
   author = {Marco-Sola, Santiago and Sammeth, Michael and Guig{\'{o}}, Roderic and Ribeca, Paolo},
   booktitle = {Nature Methods},
   doi = {10.1038/nmeth.2221},
   isbn = {1548-7105 (Electronic)$\backslash$r1548-7091 (Linking)},
   issn = {1548-7091},
   mendeley-groups = {Research articles},
   pmid = {23103880},
   title = {{The GEM mapper: fast, accurate and versatile alignment by filtration}},
   year = {2012}
   }

   -->

References
~~~~~~~~~~

[^](#ref-1) Imakaev, Maxim V and Fudenberg, Geoffrey and McCord, Rachel
Patton and Naumova, Natalia and Goloborodko, Anton and Lajoie, Bryan R
and Dekker, Job and Mirny, Leonid A. 2012. *Iterative correction of Hi-C
data reveals hallmarks of chromosome organization.*.
`URL <http://www.ncbi.nlm.nih.gov/pubmed/22941365>`__

[^](#ref-2) Ay, Ferhat and Noble, William Stafford. 2015. *Analysis
methods for studying the 3D architecture of the genome*.
`URL <http://genomebiology.com/2015/16/1/183>`__

[^](#ref-3) Wingett, Steven and Ewels, Philip and Furlan-Magaril, Mayra
and Nagano, Takashi and Schoenfelder, Stefan and Fraser, Peter and
Andrews, Simon. 2015. *HiCUP: pipeline for mapping and processing Hi-C
data.*. `URL <http://f1000research.com/articles/4-1310/v1>`__

[^](#ref-4) François Serra, Davide Baù, Mike Goodstadt, David Castillo,
Guillaume Filion and Marc A. Marti-Renom. 2017. *Automatic analysis and
3D-modelling of Hi-C data using TADbit reveals structural features of
the fly chromatin colors*.
`URL <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005665>`__

[^](#ref-5) Marco-Sola, Santiago and Sammeth, Michael and Guigó, Roderic
and Ribeca, Paolo. 2012. *The GEM mapper: fast, accurate and versatile
alignment by filtration*.
