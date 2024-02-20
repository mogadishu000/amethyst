import re
from typing import Union, List
import logging

# Logic for input processing
# https://www.rdkit.org/docs/source/rdkit.Chem.rdRGroupDecomposition.html#rdkit.Chem.rdRGroupDecomposition.RelabelMappedDummies

# Relabel dummy atoms bearing an R-group mapping (as atom map number, isotope or MDLRGroup label) 
# such that they will be displayed by the rendering code as R# rather than #*, :#, #:#, etc.

# Additional safe_mode parameter in case there are multiple types used for god know what reason

# For SMILES input use regex to find which type of dummy atoms was used
# 1 - Atom map
atom_map = re.compile(r'(\[\*\:?\d+\])')
# 2 - Isotope
isotope = re.compile(r'(\[\d+\*\])')
# 3 - Mixed

def rgroup_format(smiles: str, safe_mode: bool = False) -> int:

    if safe_mode:
        if re.search(atom_map, smiles) is not None and re.search(isotope, smiles) is not None:
            logging.info("Mixed dummy atoms")
            return 3
    else:
        if re.search(atom_map, smiles) is not None:
            logging.info("Atom map dummy atoms")
            return 1
        elif re.search(isotope, smiles) is not None:
            logging.info("Isotope dummy atoms")
            return 2
        else:
            raise ValueError("No dummy atoms found")
        
def relabel_dummy_atoms(smiles: str, type) -> str:
    pass
