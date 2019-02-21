import numpy as np

from chemspace.types import NMRVector, IRVector


class Fingerprints:
    smiles: str
    data: np.ndarray
    nmr: NMRVector
    ir: IRVector

    def __init__(self, nmr: NMRVector, ir: IRVector, **kwargs):
        assert nmr.smiles == ir.smiles
        self.nmr = nmr
        self.ir = ir
        self.smiles = nmr.smiles
        self.data = self._get_fingerprints(nmr, ir, **kwargs)

    def _get_fingerprints(self, nmr, ir, **kwargs):
        return self._get_fingerprints_by_width(nmr, ir, **kwargs)

    @staticmethod
    def _get_fingerprints_by_n(nmr, ir, nmr_bins_n=64, ir_bins_n=64):
        length = nmr_bins_n + ir_bins_n
        fingerprints = np.empty(length, dtype=float)
        nmr_per_bins = nmr.data.size // nmr_bins_n
        assert nmr_per_bins > 0
        for i in range(nmr_bins_n):
            fingerprints[i] = nmr.data[i * nmr_per_bins: (i + 1) * nmr_per_bins].mean()
        ir_per_bins = ir.data.size // ir_bins_n
        assert ir_per_bins > 0
        for i in range(ir_bins_n):
            fingerprints[nmr_bins_n + i] = ir.data[i * ir_per_bins: (i + 1) * ir_per_bins].mean()
        return fingerprints

    @staticmethod
    def _average_by_width(length, spectrum, width, fingerprints, start=0, scale=1.0):
        data = spectrum.data * scale
        for i in range(length):
            width_range = (np.array([0, 1]) + i) * width
            original_range = width_range / spectrum.bin_width
            index_range = [int(x) for x in (np.ceil(original_range[0]), np.floor(original_range[1]))]
            s = 0.0
            # w = 0.0
            if index_range[0] > original_range[0]:
                s += (index_range[0] - original_range[0]) * data[index_range[0] - 1]
                # w += index_range[0] - original_range[0]
            if index_range[1] < original_range[1] and index_range[1] < data.size:
                s += (original_range[1] - index_range[1]) * data[index_range[1]]
                # w += original_range[1] - index_range[1]
            if index_range[1] > data.size:
                index_range[1] = data.size
            s += data[index_range[0]:index_range[1]].sum()
            # w += (index_range[1] - index_range[0]) * spectrum.bin_width
            # s /= w
            s /= width
            fingerprints[start + i] = s

    @staticmethod
    def _get_fingerprints_by_width(nmr, ir, nmr_width=0.3, ir_width=111.7, ir_scale=0.1):
        nmr_length = int(np.ceil(nmr.bin_width * nmr.data.size / nmr_width))
        ir_length = int(np.ceil(ir.bin_width * ir.data.size / ir_width)) if ir_width <= 100000 else 0
        # maybe not include IR data.
        length = nmr_length + ir_length
        fingerprints = np.empty(length, dtype=float)
        Fingerprints._average_by_width(nmr_length, nmr, nmr_width, fingerprints)
        Fingerprints._average_by_width(ir_length, ir, ir_width, fingerprints, nmr_length, ir_scale)
        return fingerprints
