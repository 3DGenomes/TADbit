
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

  from pytadbit import Slice
  
  # initiate a slice object that will store all Hi-C data and analysis
  my_chrom = Slice(name='My fisrt chromsome', resolution=20000)

  # load Hi-C data
  test_chr.add_experiment("/some_path/hi-c_data.tsv", name='First Hi-C experiment')

Strange data format
-------------------

In the case Tadbit is not able to parse the input file, user can create its own parser as pass it to Slice. For example one might be interested in using [Dixon2012]_ data that appears like this:

::

  chr21	0	20000	0	0	0	0	0	0	0	0
  chr21	20000	40000	0	0	0	0	0	0	0	0
  chr21	40000	60000	0	0	0	0	0	0	0	0
  chr21	60000	80000	0	0	0	0	0	0	0	0
  chr21	80000	100000	0	0	0	0	0	0	0	0
  chr21	100000	120000	0	0	0	0	0	0	0	0
  chr21	120000	140000	0	0	0	0	0	0	0	0
  chr21	140000	160000	0	0	0	0	0	0	0	0
  

In this case we could implement a simple parser like this one:

::

  def read_dixon_matrix(f_handle):
      """
      reads from file handler (already openned)
      """
      nums = []
      start = 3
      for line in f_handle:
          values = line.split()[start:]
          nums.append([int(v) for v in values])
      f_handle.close()
      return tuple([nums[j][i] for i in xrange(len(nums)) for j in xrange(len(nums))]), len(nums)
  
And call it as follow:

::
  
  test_chr.add_experiment("/some_path/hi-c_data.tsv", name='First Hi-C experiment', 
                           parser=read_dixon_matrix)


Run Tadbit core function
========================

Once loaded the location of topologically associating domains (TADs) can be estimated.

::

  test_chr.find_tads('First Hi-C experiment')


References
==========

.. [Dixon2012] Dixon, J. R., Selvaraj, S., Yue, F., Kim, A., Li, Y., Shen, Y., Hu, M., et al. (2012). Topological domains in mammalian genomes identified by analysis of chromatin interactions. Nature, 485(7398), 376–80. doi:10.1038/nature11082
