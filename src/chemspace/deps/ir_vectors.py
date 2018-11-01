import numpy as np

class IR_Spectra_Loader:

    def __init__(self, bin_size, wn_range, data):
        """
        :param bin_size: size for each bin
        :param wn_range: duple of minimum wavenumber and maximum wavenumber
        :param data: compressed vectors dictionary of ir data (duple of nonzero indices and values at those indices)
        """
        min_wn, max_wn = wn_range
        self.num_bins = int(np.ceil((max_wn - min_wn) / bin_size))
        self.wn_range = wn_range
        self.bin_size = bin_size
        self.data = data

    def rebuild_vectors(self):
        """
        Converts compressed vector data into full vector data

        :return: dictionary of vectors in full vector form
        """
        vectors = {}
        for molecule in self.data:
            indices, values = self.data[molecule][0], self.data[molecule][1]
            vectors[molecule] = np.zeros(self.num_bins)
            vectors[molecule].put(indices, values)
        return vectors
