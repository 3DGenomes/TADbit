
Tadbit with Python
******************

.. contents::

Hi-C data format 
=================

Hi-C data is usually presented as a symmetric matrix in a tab separated file, usually like this:

::

  chrT_001	chrT_002	chrT_003	chrT_004	chrT_005	chrT_006
  chrT_001	629	164	88	105	10	35
  chrT_002	86	612	175	110	40	29
  chrT_003	159	216	437	105	43	73
  chrT_004	100	111	146	278	70	42
  chrT_005	16	36	26	79	243	39
  chrT_006	19	50	53	42	37	224


However, the number of extra columns or rows amy vary as no convention as been proposed yet. The tadbit library allows to load most of this matrices. Classically a Hi-C matrix is load as this:

::

  from pytadbit import Chromosome
  
  # initiate a chromosome object that will store all Hi-C data and analysis
  my_chrom = Chromosome(name='My fisrt chromsome', resolution=20000)

  # load Hi-C data
  test_chr.add_experiment("/some_path/hi-c_data.tsv", name='First Hi-C experiment')


In the case Tadbit is not able to parse the input file, user can create its own parser as pass it to Chromosome:

::

  
