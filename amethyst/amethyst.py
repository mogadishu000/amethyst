from rdkit.Chem.rdchem import Atom, RWMol, Mol
from rdkit.Chem.rdmolops import SanitizeMol, AddHs, Cleanup, RemoveAllHs
from rdkit.Chem.AllChem import MolFromSmiles, ReplaceSubstructs
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from typing import List, Union, Dict

# Connecting atom is the first one
def placeholder_atom_sub(core_mol: Mol, placeholder_atom: str, r_groups: Union[List[str], List[Mol]]) -> List[Mol]:
    mols = []

    if type(r_groups[0]) == Mol:
        for i in r_groups:
            mod_mol = ReplaceSubstructs(core_mol, 
                                MolFromSmiles(placeholder_atom),
                                i,
                                replaceAll=True)
            SanitizeMol(mod_mol[0])
            mols.append(mod_mol[0])
    else:
        for i in r_groups:
            mod_mol = ReplaceSubstructs(core_mol, 
                                MolFromSmiles(placeholder_atom),
                                MolFromSmiles(i),
                                replaceAll=True)
            SanitizeMol(mod_mol[0])
            mols.append(mod_mol[0])        
    return mols

# Checks:
# - bond validation
# - insertion of Mol's at right indices
def general_sub(core_mol: Mol, subs: Dict[int, Union[str, Mol, Atom]]) -> List[Mol]:
    pass