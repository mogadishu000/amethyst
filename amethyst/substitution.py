from typing import Dict, List, Union

from loguru import logger
from rdkit.Chem.AllChem import MolFromSmiles, MolToSmiles, ReplaceSubstructs
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)
from rdkit.Chem.rdchem import Atom, Mol
from rdkit.Chem.rdmolops import SanitizeMol

logger.add(sink="substitution.log", level=10)


# Connecting atom is the first one
def placeholder_atom_sub(
    core_mol: Mol,
    placeholder_atom: str,
    r_groups: Union[List[str], List[Mol]],
    inner: bool = False,
) -> List[Mol]:
    mols = []
    core_dummy_idx = core_mol.GetSubstructMatch(MolFromSmiles(placeholder_atom))

    if type(r_groups[0]) == Mol:
        for i in r_groups:
            mod_mol = ReplaceSubstructs(
                core_mol,
                MolFromSmiles(placeholder_atom),
                i,
                replaceAll=True,
                replacementConnectionPoint=(0 if not inner else core_dummy_idx),
            )
            logger.debug(f"Generated SMILES: {MolToSmiles(mod_mol[0])}")
            SanitizeMol(mod_mol[0])
            mols.append(mod_mol[0])
    else:
        for i in r_groups:
            mod_mol = ReplaceSubstructs(
                core_mol,
                MolFromSmiles(placeholder_atom),
                MolFromSmiles(i),
                replaceAll=True,
                replacementConnectionPoint=(0 if not inner else core_dummy_idx),
            )
            logger.debug(f"Generated SMILES: {MolToSmiles(mod_mol[0])}")
            SanitizeMol(mod_mol[0])
            mols.append(mod_mol[0])
    return mols


# Checks:
# - actually use molzip with relabelMappedDummies
# - insertion of Mol's at right indices
# - if it works with rdRGroupDecomposition output
# -- https://www.rdkit.org/docs/source/rdkit.Chem.rdRGroupDecomposition.html#rdkit.Chem.rdRGroupDecomposition.RGroupDecompositionParameters
def general_sub(core_mol: Mol, subs: Dict[int, Union[str, Mol, Atom]]) -> List[Mol]:
    #
    pass
