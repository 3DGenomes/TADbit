"""
02 Jul 2013

script that follows Tadbit tutorial presented in the documentation
"""

from pytadbit import Chromosome

# initiate a chromosome object that will store all Hi-C data and analysis
my_chrom = Chromosome(name='My fisrt chromsome')

# load Hi-C data
my_chrom.add_experiment('First Hi-C experiment', xp_handler="sample_data/HIC_k562_chr19_chr19_100000_obs.txt", resolution=100000)
my_chrom.add_experiment('Second Hi-C experiment', xp_handler="sample_data/HIC_gm06690_chr19_chr19_100000_obs.txt", resolution=100000)

# run core tadbit function to find TADs, on each experiment
my_chrom.find_tad('First Hi-C experiment')
my_chrom.find_tad('Second Hi-C experiment')

print my_chrom.experiments


my_chrom.align_experiments(names=["First Hi-C experiment", "Second Hi-C experiment"])

print my_chrom.alignment

ali = my_chrom.alignment[('First Hi-C experiment', 'Second Hi-C experiment')]


print ali.write_alignment(ftype='html')
