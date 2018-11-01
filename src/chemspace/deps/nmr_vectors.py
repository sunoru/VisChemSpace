class Mol_NMR_Vector():
	'''
	Class representing a NMR spectra vector
	'''

	def __init__(self, smiles, nmr_vector, bin_width):
		'''
		Constructor for a Mol_Point object
		@param smiles - the smiles string for the molecule
		@param nmr_vector - a numpy array representing the nmr spectra in vector form
		@param bin_width - the bin width for each bin in the vector
		'''

		self.smiles = smiles
		self.nmr_vector = nmr_vector
		self.bin_width = bin_width
