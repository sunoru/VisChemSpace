import numpy as np

from chemspace.deps.ir_vectors import IR_Spectra_Loader
from chemspace.deps.nmr_vectors import Mol_NMR_Vector


class SpectrumVector:
    smiles: str
    data: np.ndarray
    bin_width: float

    def __init__(self, smiles, data, bin_width):
        self.smiles = smiles
        self.data = data
        self.bin_width = bin_width


class NMRVector(SpectrumVector):
    @classmethod
    def from_old(cls, old_data: Mol_NMR_Vector):
        return cls(old_data.smiles, np.array(old_data.nmr_vector), old_data.bin_width)


class IRVector(SpectrumVector):
    wave_number_range: (float, float)

    def __init__(self, smiles, data, bin_width, wave_number_range):
        super().__init__(smiles, data, bin_width)
        self.wave_number_range = wave_number_range

    @classmethod
    def from_old(cls, old_data: IR_Spectra_Loader, smiles: str):
        data = old_data.data[smiles]
        indices, values = data
        vector = np.zeros(old_data.num_bins)
        vector.put(indices, values)
        return cls(smiles, vector, old_data.bin_size, old_data.wn_range)
