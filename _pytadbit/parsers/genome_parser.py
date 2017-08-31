"""
17 nov. 2014

convert a bunch of fasta files, or a single multi fasta file, into a dictionary
"""

from collections import OrderedDict
from os import path
import re

from pytadbit.utils.file_handling import magic_open

def parse_fasta(f_names, chr_names=None, chr_filter=None, chr_regexp=None,
                verbose=True, save_cache=True, reload_cache=False, only_length=False):
    """
    Parse a list of fasta files, or just one fasta.

    WARNING: The order is important

    :param f_names: list of pathes to files, or just a single path
    :param None chr_names: pass list of chromosome names, or just one. If None
       are passed, then chromosome names will be inferred from fasta headers
    :param None chr_filter: use only chromosome in the input list
    :param None chr_regexp: use only chromosome matching
    :param True save_cache: save a cached version of this file for faster
       loadings (~4 times faster)
    :param False reload_cache: reload cached genome
    :param False only_length: returns dictionary with length of genome,not sequence

    :returns: a sorted dictionary with chromosome names as keys, and sequences
       as values (sequence in upper case)
    """
    if isinstance(f_names, str):
        f_names = [f_names]

    if len(f_names) == 1:
        fname = f_names[0] + '_genome.TADbit'
    else:
        fname = path.join(path.commonprefix(f_names), 'genome.TADbit')
    if path.exists(fname) and not reload_cache:
        if verbose:
            print 'Loading cached genome'
        genome_seq = OrderedDict()
        for line in open(fname):
            if line.startswith('>'):
                c = line[1:].strip()
            else:
                if only_length:
                    genome_seq[c] = len(line.strip())
                else:
                    genome_seq[c] = line.strip()
        return genome_seq

    if isinstance(chr_names, str):
        chr_names = [chr_names]

    if chr_filter:
        bad_chrom = lambda x: not x in chr_filter
    else:
        bad_chrom = lambda x: False

    if chr_regexp:
        chr_regexp = re.compile(chr_regexp)
    else:
        chr_regexp = re.compile('.*')

    genome_seq = OrderedDict()
    if len(f_names) == 1:
        header = None
        seq = []
        for line in magic_open(f_names[0]):
            if line.startswith('>'):
                if header:
                    genome_seq[header] = ''.join(seq).upper()
                header = line[1:].split()[0]
                if bad_chrom(header) or not chr_regexp.match(header):
                    header = 'UNWANTED'
                elif not chr_names:
                    if verbose:
                        print 'Parsing %s' % (header)
                else:
                    header = chr_names.pop(0)
                    if verbose:
                        print 'Parsing %s as %s' % (line[1:].rstrip(), header)
                seq = []
                continue
            seq.append(line.rstrip())
        if only_length:
            genome_seq[header] = len(seq)
        else:
            genome_seq[header] = ''.join(seq).upper()
        if 'UNWANTED' in genome_seq:
            del(genome_seq['UNWANTED'])
    else:
        for fnam in f_names:
            fhandler = magic_open(fnam)
            try:
                while True:
                    if not chr_names:
                        header = fhandler.next()
                        if header.startswith('>'):
                            header = header[1:].split()[0]
                            if bad_chrom(header) or not chr_regexp.match(header):
                                header = 'UNWANTED'
                            genome_seq[header] = ''
                            break
                    else:
                        _ = fhandler.next()
                        header = chr_names.pop(0)
                        if bad_chrom(header):
                            header = 'UNWANTED'
                        genome_seq[header] = ''
                        break
            except StopIteration:
                raise Exception('No crocodiles found, is it fasta?')
            if only_length:
                genome_seq[header] = sum(len(l.rstrip()) for l in fhandler)
            else:
                genome_seq[header] = ''.join([l.rstrip() for l in fhandler]).upper()
        if 'UNWANTED' in genome_seq:
            del(genome_seq['UNWANTED'])
    if save_cache and not only_length:
        if verbose:
            print 'saving genome in cache'
        if len(f_names) == 1:
            fname = f_names[0] + '_genome.TADbit'
        else:
            fname = path.join(path.commonprefix(f_names), 'genome.TADbit')
        out = open(fname, 'w')
        for c in genome_seq:
            out.write('>%s\n%s\n' % (c, genome_seq[c]))
        out.close()
    return genome_seq
