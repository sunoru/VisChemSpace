import pickle
import sys
import chemspace.deps.nmr_vectors as nmr_vectors
from chemspace.deps.ir_vectors import IR_Spectra_Loader 


def load_file(filename):
    with open(filename, 'rb') as fi:
        return pickle.load(fi)


def load_data(nmr_file, ir_file):
    nmr_data = load_file(nmr_file)
    ir_data = load_file(ir_file)
    return nmr_data, ir_data


load_nmr = load_ir = load_file
sys.modules['vectorize'] = nmr_vectors
sys.modules['__main__'].IR_Spectra_Loader = IR_Spectra_Loader
