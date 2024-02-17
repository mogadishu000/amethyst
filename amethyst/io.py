# Logic for input processing
# https://www.rdkit.org/docs/source/rdkit.Chem.rdRGroupDecomposition.html#rdkit.Chem.rdRGroupDecomposition.RelabelMappedDummies

# Relabel dummy atoms bearing an R-group mapping (as atom map number, isotope or MDLRGroup label) 
# such that they will be displayed by the rendering code as R# rather than #*, :#, #:#, etc.

# Additional safe_mode parameter in case there are multiple types used for god know what reason

# For SMILES input use regex to find which type of dummy atoms was used
def rgroup_format():
    pass
