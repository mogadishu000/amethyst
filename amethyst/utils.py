from rdkit.Chem.rdchem import Atom, Mol
from rdkit.Chem.rdmolfiles import MolToSmiles
from rdkit.Chem.Draw import rdMolDraw2D
from typing import Tuple

# Checks wheter passed SMILES string is a molecule or an atom
def mol_or_atom():
    pass

# Returns an image of molecule with atom idx's
def depict_mol(mol: Mol, filename: str = '', size: Tuple(int) = (1000, 1000)):
    if filename is '':
        filename = MolToSmiles(mol)

    d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])

    d.drawOptions().addAtomIndices = True
    d.drawOptions().addAtomLabels = True
    d.DrawMolecule(mol)
    d.FinishDrawing()
    d.WriteDrawingText(f'png/{filename}.png')