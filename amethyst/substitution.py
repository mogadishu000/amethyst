import itertools
from typing import List, Union

from loguru import logger
from rdkit.Chem.AllChem import MolFromSmiles, MolToSmiles, ReplaceSubstructs
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import SanitizeMol, molzip

from amethyst.io import Substituents, parse_file_input

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
def general_sub(core_mol: Mol, subs: List[Substituents]) -> List[Mol]:
    r_groups: List[Mol] = [x.subs for x in subs]
    combinations: List[Mol] = [x for x in itertools.product(*r_groups)]
    output_mols: List[Mol] = []

    for sets in combinations:
        sub_mol = core_mol
        for r in sets:
            sub_mol = molzip(sub_mol, r)
        output_mols.append(sub_mol)

    return output_mols
