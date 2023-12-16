import pytest
from typing import List
from amethyst.utils import mol_or_atom
from rdkit.Chem.rdchem import Atom, Mol

def test_atom_recognition():
    assert(type(mol_or_atom('Br')) == Atom)

def test_molecule_recognition():
    assert(type(mol_or_atom('CC')) == Mol)