import numpy as np

from chemspace.types import NMRVector, IRVector


class Fingerprints:
    smiles: str
    data: np.ndarray

    def __init__(self, nmr: NMRVector, ir: IRVector, **kwargs):
        assert nmr.smiles == ir.smiles
        self.smiles = nmr.smiles
        self.data = self._get_fingerprints(nmr, ir, **kwargs)

    @staticmethod
    def _get_fingerprints_by_n(nmr, ir, nmr_bins_n=512, ir_bins_n=512):
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

    def _get_fingerprints(self, nmr, ir, **kwargs):
        return self._get_fingerprints_by_n(nmr, ir, **kwargs)
