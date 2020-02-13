Preparation of the reference genome
===================================

Usually NGS reads are mapped against a reference genome containing only
the assembled chromosomes, and not the remaining contigs. And this
methodology is perfectly valid. However in order to decrease the
probability of having mapping errors, adding all unassembled contigs may
help:

   For variant discovery, RNA-seq and ChIP-seq, it is recommended to use
   the entire primary assembly, including assembled chromosomes AND
   unlocalized/unplaced contigs, for the purpose of read mapping. Not
   including unlocalized and unplaced contigs potentially leads to more
   mapping errors.

   *from: http://lh3lh3.users.sourceforge.net/humanref.shtml*

We are thus going to download full chromosomes and unassembled contigs.
From these sequences we are then going to create two reference genomes:
- one “classic” reference genome with only assembled chromosomes, used
to compute statistics on the genome (GC content, number of restriction
sites or mappability) - one that would contain all chromosomes and
unassembled contigs, used exclusively for mapping.

*Mus musculus*\ ’s reference genome sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We search for the most recent reference genome corresponding to Mouse
(https://www.ncbi.nlm.nih.gov/genome?term=mus%20musculus).

From there we obtain these identifiers:

.. code:: ipython3

    species  = 'Mus_musculus'
    taxid    = '10090'
    assembly = 'GRCm38.p6'
    genbank  = 'GCF_000001635.26'

*The variables defined above can be modified for any other species,
resulting in new results for the following commands.*

Download from the NCBI
----------------------

List of chromosomes/contigs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    sumurl = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{0}/{1}/{2}/{3}/{4}_{5}/'
              '{4}_{5}_assembly_report.txt').format(genbank[:3], genbank[4:7], genbank[7:10],
                                                    genbank[10:13], genbank, assembly)

    crmurl = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
              '?db=nuccore&id=%s&rettype=fasta&retmode=text')

.. code:: ipython3

    print sumurl


.. ansi-block::

    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_assembly_report.txt


.. code:: ipython3

    ! wget -q $sumurl -O chromosome_list.txt

.. code:: ipython3

    ! head chromosome_list.txt


.. ansi-block::

    # Assembly name:  GRCm38.p6
    # Description:    Genome Reference Consortium Mouse Build 38 patch release 6 (GRCm38.p6)
    # Organism name:  Mus musculus (house mouse)
    # Infraspecific name:  strain=C57BL/6J
    # Taxid:          10090
    # BioProject:     PRJNA20689
    # Submitter:      Genome Reference Consortium
    # Date:           2017-9-15
    # Assembly type:  haploid-with-alt-loci
    # Release type:   patch


Sequences of each chromosome/contig
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    import os

.. code:: ipython3

    dirname = 'genome'
    ! mkdir -p {dirname}

For each contig/chromosome download the corresponding FASTA file from
NCBI

.. code:: ipython3

    contig = []
    for line in open('chromosome_list.txt'):
        if line.startswith('#'):
            continue
        seq_name, seq_role, assigned_molecule, _, genbank, _, refseq, _ = line.split(None, 7)
        if seq_role == 'assembled-molecule':
            name = 'chr%s.fasta' % assigned_molecule
        else:
            name = 'chr%s_%s.fasta' % (assigned_molecule, seq_name.replace('/', '-'))
        contig.append(name)

        outfile = os.path.join(dirname, name)
        if os.path.exists(outfile) and os.path.getsize(outfile) > 10:
            continue
        error_code = os.system('wget "%s" --no-check-certificate -O %s' % (crmurl % (genbank), outfile))
        if error_code:
            error_code = os.system('wget "%s" --no-check-certificate -O %s' % (crmurl % (refseq), outfile))
        if error_code:
            print genbank


.. ansi-block::

    CM001008.2


Concatenate all contigs/chromosomes into single files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def write_to_fasta(line):
        contig_file.write(line)

    def write_to_fastas(line):
        contig_file.write(line)
        simple_file.write(line)

.. code:: ipython3

    os.system('mkdir -p {}/{}-{}'.format(dirname, species, assembly))




.. ansi-block::

    0



.. code:: ipython3

    contig_file = open('{0}/{1}-{2}/{1}-{2}_contigs.fa'.format(dirname, species, assembly),'w')
    simple_file = open('{0}/{1}-{2}/{1}-{2}.fa'.format(dirname, species, assembly),'w')

    for molecule in contig:
        fh = open('{0}/{1}'.format(dirname, molecule))
        oline = '>%s\n' % (molecule.replace('.fasta', ''))
        _ = fh.next()
        # if molecule is an assembled chromosome we write to both files, otherwise only to the *_contigs one
        write = write_to_fasta if '_' in molecule else write_to_fastas
        for line in fh:
            write(oline)
            oline = line
        # last line usually empty...
        if line.strip():
            write(line)
    contig_file.close()
    simple_file.close()

Remove all the other files (with single chromosome/contig)

.. code:: ipython3

    ! rm -f {dirname}/*.fasta

Creation of an index file for GEM mapper
----------------------------------------

.. code:: ipython3

    ! gem-indexer -T 8 -i {dirname}/{species}-{assembly}/{species}-{assembly}_contigs.fa -o {dirname}/{species}-{assembly}/{species}-{assembly}_contigs


.. ansi-block::

    Welcome to GEM-indexer build 1.423 (beta) - (2013/04/01 01:02:13 GMT)
     (c) 2008-2013 Paolo Ribeca <paolo.ribeca@gmail.com>
     (c) 2010-2013 Santiago Marco Sola <santiagomsola@gmail.com>
     (c) 2010-2013 Leonor Frias Moya <leonor.frias@gmail.com>
    For the terms of use, run the program with the option --show-license.
    ************************************************************************
    * WARNING: this is a beta version, provided for testing purposes only; *
    *          check for updates at <http://gemlibrary.sourceforge.net>.   *
    ************************************************************************
    Creating sequence and location files... done.
    Computing DNA BWT (likely to take long)... done.
    Generating index (likely to take long)... done.
    Cleaning up... done.


The path to the index file will be:
**``{dirname}/{species}-{assembly}/{species}_contigs.gem``**

Compute mappability values needed for bias specific normalizations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case we can use the FASTA of the genome whithout contigs and
follow these step:

.. code:: ipython3

    ! gem-indexer -i {dirname}/{species}-{assembly}/{species}-{assembly}.fa \
       -o {dirname}/{species}-{assembly}/{species}-{assembly} -T 8

    ! gem-mappability -I {dirname}/{species}-{assembly}/{species}-{assembly}.gem -l 50 \
       -o {dirname}/{species}-{assembly}/{species}-{assembly}.50mer -T 8

    ! gem-2-wig -I {dirname}/{species}-{assembly}/{species}-{assembly}.gem \
       -i {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.mappability \
       -o {dirname}/{species}-{assembly}/{species}-{assembly}.50mer

    ! wigToBigWig {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.wig \
       {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.sizes \
       {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.bw

    ! bigWigToBedGraph {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.bw  \
       {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.bedGraph


.. ansi-block::

    Welcome to GEM-indexer build 1.423 (beta) - (2013/04/01 01:02:13 GMT)
     (c) 2008-2013 Paolo Ribeca <paolo.ribeca@gmail.com>
     (c) 2010-2013 Santiago Marco Sola <santiagomsola@gmail.com>
     (c) 2010-2013 Leonor Frias Moya <leonor.frias@gmail.com>
    For the terms of use, run the program with the option --show-license.
    ************************************************************************
    * WARNING: this is a beta version, provided for testing purposes only; *
    *          check for updates at <http://gemlibrary.sourceforge.net>.   *
    ************************************************************************
    Creating sequence and location files... done.
    Computing DNA BWT (likely to take long)... done.
    Generating index (likely to take long)... done.
    Cleaning up... done.
    Welcome to GEM-mappability build 1.315 (beta) - (2013/03/29 02:59:40 GMT)
     (c) 2008-2013 Paolo Ribeca <paolo.ribeca@gmail.com>
     (c) 2010-2013 Santiago Marco Sola <santiagomsola@gmail.com>
     (c) 2010-2013 Leonor Frias Moya <leonor.frias@gmail.com>
    For the terms of use, run the program with the option --show-license.
    ************************************************************************
    * WARNING: this is a beta version, provided for testing purposes only; *
    *          check for updates at <http://gemlibrary.sourceforge.net>.   *
    ************************************************************************
    Sat Jan 12 13:01:12 2019 -- Loading index (likely to take long)... done.
    Sat Jan 12 13:01:23 2019 -- Starting (2647560565 positions to go)...
    Sat Jan 12 20:34:17 2019 -- Writing frequencies to disk...
    Sat Jan 12 20:34:27 2019 -- ...done.
    Welcome to GEM-2-wig build 1.423 (beta) - (2013/04/01 01:02:13 GMT)
     (c) 2008-2013 Paolo Ribeca <paolo.ribeca@gmail.com>
     (c) 2010-2013 Santiago Marco Sola <santiagomsola@gmail.com>
     (c) 2010-2013 Leonor Frias Moya <leonor.frias@gmail.com>
    For the terms of use, run the program with the option --show-license.
    ************************************************************************
    * WARNING: this is a beta version, provided for testing purposes only; *
    *          check for updates at <http://gemlibrary.sourceforge.net>.   *
    ************************************************************************
    Sat Jan 12 20:34:28 2019 -- Loading index (likely to take long)... done.
    Sat Jan 12 20:34:41 2019 -- Inverting locations... done.
    genome/Mus_musculus-GRCm38.p6/Mus_musculus-GRCm38.p6.50mer.bw is not a bpt b-plus tree index file


Cleanup
~~~~~~~

.. code:: ipython3

    ! rm -f {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.mappability
    ! rm -f {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.wig
    ! rm -f {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.bw
    ! rm -f {dirname}/{species}-{assembly}/{species}-{assembly}.50mer.sizes
    ! rm -f {dirname}/{species}-{assembly}/*.log
