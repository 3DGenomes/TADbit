"""
September 2, 2019.

"""

from __future__ import division
from __future__ import print_function
import os
import datetime
import json
import h5py
import numpy as np
from warnings       import warn
from math           import ceil
from collections    import OrderedDict
from time           import time

def printime(msg):
    print (msg +
           (' ' * (79 - len(msg.replace('\n', '')))) +
           '[' +
           str(datetime.datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
           ']')

def is_cooler(fname, resolution=None):
    """
    Check if file is a cooler and contains the wanted resolution

    :param f: an iterable (typically an open file).
    :param None resolution: matrix resolution.
    """

    try:
        with h5py.File(fname, "r") as f:
            resolution = resolution or list(f['resolutions'].keys())[0]
            if str(resolution) in f['resolutions']:
                return True
    except ValueError:
        warn('WARNING: cooler file exists but does not contain wanted resolution')
    except:
        pass
    return False

def parse_cooler(fname, resolution=None, normalized=False,
                 raw_values = False):
    """
    Read matrix stored in cooler

    :param f: an iterable (typically an open file).
    :param None resolution: matrix resolution.
    :param False normalized: whether to apply weights
    :param False raw_values: return separated raw and weights

    :returns: An iterator to be converted in dictionary, matrix size, raw_names
       as list of tuples (chr, pos), dictionary of masked bins, and boolean
       reporter of symetric transformation
    """

    with h5py.File(fname, "r") as f:

        resolution = resolution or list(f['resolutions'].keys())[0]
        root_grp = f['resolutions'][str(resolution)]

        chrom = root_grp["chroms"]["name"].value
        idregion = dict(list(zip(list(range(len(chrom))), [reg for reg in chrom])))

        try:
            chrom = [idregion[c].decode() for c in root_grp["bins"]["chrom"]]
        except (UnicodeDecodeError, AttributeError):
            chrom = [str(idregion[c]) for c in root_grp["bins"]["chrom"]]

        if raw_values:
            starts = root_grp["bins"]["start"]
            header = OrderedDict()
            for chromi, starti in zip(chrom, starts):
                header[chromi] = starti // resolution + 1
            size = 0
            for chromi in header:
                size += header[chromi] 
        else:
            starts = ['%d-%d'%(c+1,c+int(resolution)) for c in root_grp["bins"]["start"]]
            header = [(chromi, starti) for chromi, starti in zip(chrom, starts)]
            size = len(header)
        masked = {}
        if normalized and "weight" in root_grp["bins"]:
            weights = root_grp["bins"]["weight"].value
        else:
            weights = [1 for _ in range(size)]
        bin1_id = root_grp["pixels"]["bin1_id"].value
        bin2_id = root_grp["pixels"]["bin2_id"].value
        counti = root_grp["pixels"]["count"].value
        num = int if not normalized else float
        if raw_values:
            items = [(int(row) + int(col) * size, num(val))
                 for row, col, val in zip(bin1_id, bin2_id, counti)]
        else:
            items = [(int(row) + int(col) * size, num(val)*weights[row]*weights[col])
                     for row, col, val in zip(bin1_id, bin2_id, counti)]
    if raw_values:
        return items, weights, size, header
    else:
        return items, size, header, masked, False

class cooler_file(object):
    """
        Cooler file wrapper.
    """
    def __init__(self, outcool, resolution, sections, regions,
                 h5opts=None, verbose=False):
        """
        :param outcool: path to cool file to be created
        :param resolution: resolution of the matrix
        :param sections: dictionary with chromosomes and lengths
        :param regions: list of chromosomes present in the matrix, order matters
        :param None h5opts: options for the h5py outcool file
        :param True verbose: speak

        """

        self.root_grp = 'resolutions'
        self.resolution = resolution
        if not os.path.exists(outcool):
            with h5py.File(outcool, "w") as f:
                grp = f.create_group(self.root_grp)
                grp.create_group(str(self.resolution))
            if verbose:
                printime('Opening cooler %s'%outcool)
        else:
            with h5py.File(outcool, "r+") as f:
                try:
                    f[self.root_grp].create_group(str(self.resolution))
                except ValueError:
                    pass #existing cooler, maybe we want to modify or update
            if verbose:
                printime('Opening existing cooler %s'%outcool)

        self.outcool = outcool
        self.h5opts = _set_h5opts(h5opts)
        self.name = outcool
        self.sections = sections
        self.regions = regions
        self.nbins = sum([int(ceil(self.sections[reg]/self.resolution))
                          for reg in list(OrderedDict.fromkeys(self.regions))])
        self.sec_offset = 0
        for crm in self.sections:
            if crm == self.regions[0]:
                break
            self.sec_offset += int(ceil(self.sections[crm]/self.resolution))
        self.nnz = 0
        self.ncontacts = 0
        self.ichunk = 0
        self.buff = []
        self.nbuff = 0
        self.startj = 0
        self.startk = 0
        self.verbose = verbose
    
    def create_bins(self):
        """
        Write bins to cooler file.

        :param sections: dictionary with chromosomes and lengths

        """
        regions = list(OrderedDict.fromkeys(self.regions))
        bins = ((reg, p, p + self.resolution) for r, reg in enumerate(regions)
                for p in range(0,self.sections[reg],self.resolution))

        if self.verbose:
            printime('Writing regions')
        self.write_regions()

        if self.verbose:
            printime('Writing bins')
        self.write_bins(bins)

    def write_regions(self):
        """
        Write the regions table.

        """
        regions = list(OrderedDict.fromkeys(self.regions))
        with h5py.File(self.outcool, "r+") as f:
            root_grp = f[self.root_grp][str(self.resolution)]
            chr_names = np.array([reg for reg in regions], dtype=np.dtype("S"))
            grp = root_grp.create_group("chroms")
            grp.create_dataset("name", shape=(len(regions),), dtype=chr_names.dtype,
                               data=chr_names, **self.h5opts)
            grp.create_dataset("length", shape=(len(regions),), dtype=np.int32,
                               data=[self.sections[reg] for reg in regions], **self.h5opts)

    def write_bins(self, bins_it, chrom_as_enum=True):
        """
        Write the bins table.

        :param bins_it: iterator with list of tuples (chr, start, end)
            with the bins of the matrix

        """
        regions = list(OrderedDict.fromkeys(self.regions))
        bins = [bin_i for bin_i in bins_it]
        with h5py.File(self.outcool, "r+") as f:
            root_grp = f[self.root_grp][str(self.resolution)]
            grp = root_grp.create_group("bins")

            nregions = len(regions)
            idregion = dict(list(zip([reg for reg in regions], list(range(nregions)))))

            chrom_ids = [idregion[bin_i[0]] for bin_i in bins]
            if chrom_as_enum:
                chrom_dtype = h5py.special_dtype(enum=(np.int32, idregion))
            else:
                chrom_dtype = np.int32

            # Store bins
            try:
                chrom_dset = grp.create_dataset("chrom", shape=(self.nbins,), dtype=chrom_dtype,
                                                data=chrom_ids, **self.h5opts)
            except ValueError:
                # If too many scaffolds for HDF5 enum header,
                # try storing chrom IDs as raw int instead
                if chrom_as_enum:
                    chrom_as_enum = False
                    chrom_dtype = np.int32
                    chrom_dset = grp.create_dataset("chrom", shape=(self.nbins,), dtype=chrom_dtype,
                                                    data=chrom_ids, **self.h5opts)
                else:
                    raise
            if not chrom_as_enum:
                chrom_dset.attrs["enum_path"] = u"/chroms/name"
            bin_data = [bin_i[1] for bin_i in bins]
            grp.create_dataset("start", shape=(self.nbins,), dtype=np.int32,
                               data=bin_data, **self.h5opts)
            bin_data = [bin_i[2] for bin_i in bins]
            grp.create_dataset("end", shape=(self.nbins,), dtype=np.int32,
                               data=bin_data, **self.h5opts)

    def write_weights(self, weights_row, weights_col, start1=None,
                      end1=None, start2=None, end2=None):
        """
        Write the weights in the bins table.

        :param weights_row: list of row weights
        :param weights_col: list of column weights
        :param None start1: start bin of the first region of the matrix
        :param None start2: start bin of the second region of the matrix
        :param None end1: end bin of the first region of the matrix
        :param None end2: end bin of the second region of the matrix

        """
        startj = 0 if start1 is None else start1
        startk = 0 if start2 is None else start2
        endj = self.sections[self.regions[0]] if end1 is None else end1
        endk = self.sections[self.regions[-1]] if end2 is None else end2

        full_weights = np.zeros(shape=(self.nbins,), dtype=np.float64)
        full_weights[startj-self.sec_offset:endj-self.sec_offset] = weights_row
        full_weights[startk-self.sec_offset:endk-self.sec_offset] = weights_col
        with h5py.File(self.outcool, "r+") as f:
            root_grp = f[self.root_grp][str(self.resolution)]
            grp = root_grp["bins"]
            grp.create_dataset("weight", dtype=np.float64, data=full_weights, **self.h5opts)

    def prepare_matrix(self, start1=None, start2=None):
        """
        Prepare matrix datasets to be written as chunks.

        :param None start1: start bin of the first region of the matrix
        :param None start2: start bin of the second region of the matrix

        """

        self.startj = 0 if start1 is None else start1
        self.startk = 0 if start2 is None else start2
        max_size = self.nbins * (self.nbins - 1) // 2 + self.nbins

        if self.verbose:
            printime('Prepare matrix')

        with h5py.File(self.outcool, "r+") as f:
            root_grp = f[self.root_grp][str(self.resolution)]
            grp = root_grp.create_group("pixels")
            columns = ["bin1_id","bin2_id","count"]
            dsets_dtypes = [np.int64, np.int64, np.int32]
            for col, dset_dtype in zip(columns, dsets_dtypes):
                grp.create_dataset(col, shape=(max_size,), dtype=dset_dtype,
                                   maxshape=(None,), chunks=True, **self.h5opts)

    def write_iter(self, ichunk, j, k, v):
        """
        Write bin1, bin2, value to buffer. When the chunk number changes the buffer
        is written to the h5py file.

        :param ichunk: Chunk number
        :param j: row number
        :param k: column number
        :param v: interaction value

        """
        if self.ichunk != ichunk:
            self.buff.sort()
            with h5py.File(self.outcool, "r+") as f:
                root_grp = f[self.root_grp][str(self.resolution)]
                grp = root_grp["pixels"]
                dsets = ["bin1_id","bin2_id","count"]
                for i, dset in enumerate(dsets):
                    grp[dset].resize((self.nnz + self.nbuff,))
                    grp[dset][self.nnz : self.nnz + self.nbuff] = [int(ch[i]) for ch in self.buff]
            self.nnz += self.nbuff
            self.ncontacts += sum([int(ch[2]) for ch in self.buff])

            del self.buff[:]
            self.nbuff = 0
        vals = (j+(self.startj-self.sec_offset),k+(self.startk-self.sec_offset),v)
        self.buff.append(vals)
        self.nbuff += 1
        self.ichunk = ichunk

    def close(self):
        """
        Copy remaining buffer to file, index the pixelsand complete information
        """
        # copy remaining reads in buffer
        if self.nbuff > 0:
            self.buff.sort()
            with h5py.File(self.outcool, "r+") as f:
                root_grp = f[self.root_grp][str(self.resolution)]
                grp = root_grp["pixels"]
                dsets = ["bin1_id","bin2_id","count"]
                for i, dset in enumerate(dsets):
                    grp[dset].resize((self.nnz + self.nbuff,))
                    grp[dset][self.nnz : self.nnz + self.nbuff] = [int(ch[i]) for ch in self.buff]
            self.nnz += self.nbuff
            self.ncontacts += sum([int(ch[2]) for ch in self.buff])

            del self.buff[:]
            self.nbuff = 0
        self.ichunk = 0
        self.write_indexes()
        self.write_info()

    def write_indexes(self):
        """
        Write the indexes from existing bins and pixels.

        """
        if self.verbose:
            printime("Writing indexes")

        with h5py.File(self.outcool, "r+") as f:
            root_grp = f[self.root_grp][str(self.resolution)]
            grp = root_grp.create_group("indexes")

            chrom_offset = index_bins(root_grp["bins"], len(self.regions), self.nbins)
            bin1_offset = index_pixels(root_grp["pixels"], self.nbins, self.nnz)

            grp.create_dataset(
                "chrom_offset",
                shape=(len(chrom_offset),),
                dtype=np.int64,
                data=chrom_offset,
                **self.h5opts
            )
            grp.create_dataset(
                "bin1_offset",
                shape=(len(bin1_offset),),
                dtype=np.int64,
                data=bin1_offset,
                **self.h5opts
            )

    def write_info(self):
        """
        Write the file description and metadata attributes.

        """

        if self.verbose:
            printime("Writing info")
        info = {}
        info["bin-type"] = u"fixed"
        info["bin-size"] = self.resolution
        info["storage-mode"] = u"symmetric-upper"
        info["nchroms"] = len(self.regions)
        info["nbins"] = self.nbins
        info["sum"] = self.ncontacts
        info["nnz"] = self.nnz

        assert "nbins" in info
        assert "nnz" in info
        info.setdefault("genome-assembly", "unknown")
        info["metadata"] = json.dumps(info.get("metadata", {}))
        info["creation-date"] = datetime.datetime.now().isoformat()
        info["generated-by"] = "TADbit"
        info["format"] = u"HDF5::MCOOL"
        info["format-version"] = 2
        info["format-url"] = u"https://github.com/mirnylab/cooler"

        with h5py.File(self.outcool, "r+") as f:
            root_grp = f[self.root_grp][str(self.resolution)]
            root_grp.attrs.update(info)

### Kindly taken from cooler code, sharing for love
def _set_h5opts(h5opts):
    result = {}
    if h5opts is not None:
        result.update(h5opts)
    available_opts = {
        "chunks",
        "maxshape",
        "compression",
        "compression_opts",
        "scaleoffset",
        "shuffle",
        "fletcher32",
        "fillvalue",
        "track_times",
    }
    for key in list(result.keys()):
        if key not in available_opts:
            raise ValueError("Unknown storage option '{}'.".format(key))
    result.setdefault("compression", "gzip")
    if result["compression"] == "gzip" and "compression_opts" not in result:
        result["compression_opts"] = 6
    result.setdefault("shuffle", True)
    return result

def asarray_or_dataset(x):
    return x if isinstance(x, h5py.Dataset) else np.asarray(x)

def rlencode(array, chunksize=None):
    """
    Run length encoding.
    Based on http://stackoverflow.com/a/32681075, which is based on the rle
    function from R.

    Parameters
    ----------
    x : 1D array_like
        Input array to encode
    dropna: bool, optional
        Drop all runs of NaNs.

    Returns
    -------
    start positions, run lengths, run values

    """
    where = np.flatnonzero
    array = asarray_or_dataset(array)
    n = len(array)
    if n == 0:
        return (
            np.array([], dtype=int),
            np.array([], dtype=int),
            np.array([], dtype=array.dtype),
        )

    if chunksize is None:
        chunksize = n

    starts, values = [], []
    last_val = np.nan
    for i in range(0, n, chunksize):
        x = array[i : i + chunksize]
        locs = where(x[1:] != x[:-1]) + 1
        if x[0] != last_val:
            locs = np.r_[0, locs]
        starts.append(i + locs)
        values.append(x[locs])
        last_val = x[-1]
    starts = np.concatenate(starts)
    lengths = np.diff(np.r_[starts, n])
    values = np.concatenate(values)

    return starts, lengths, values

def index_pixels(grp, n_bins, nnz):
    bin1 = grp["bin1_id"]
    bin1_offset = np.zeros(n_bins + 1, dtype=np.int64)
    curr_val = 0
    for start, length, value in zip(*rlencode(bin1, 1000000)):
        bin1_offset[curr_val : value + 1] = start
        curr_val = value + 1
    bin1_offset[curr_val:] = nnz
    return bin1_offset

def index_bins(grp, n_chroms, n_bins):
    chrom_ids = grp["chrom"]
    chrom_offset = np.zeros(n_chroms + 1, dtype=np.int64)
    curr_val = 0
    for start, length, value in zip(*rlencode(chrom_ids)):
        chrom_offset[curr_val : value + 1] = start
        curr_val = value + 1
    chrom_offset[curr_val:] = n_bins
    return chrom_offset

### Kindly taken from cooler code, sharing for love