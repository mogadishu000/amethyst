from typing import List

import pytest
from rdkit.Chem.AllChem import MolFromSmiles, MolToInchi

from amethyst.substitution import placeholder_atom_sub

core_smi: str = "CCN"  # Ethylamine
placeholder_atom: str = "N"  # Nitrogen atom
r_groups = ["NCC", "NC(O)C"]
core_mol = MolFromSmiles(core_smi)


def test_placeholder_substitution():
    sub_mol = placeholder_atom_sub(core_mol, placeholder_atom, r_groups)
    assert MolToInchi(sub_mol[0]) == MolToInchi(MolFromSmiles("CCNCC"))
    assert MolToInchi(sub_mol[1]) == MolToInchi(MolFromSmiles("CCNC(O)C"))
