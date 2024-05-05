# NOTE - checks for relative and absolute file paths
import pytest
from rdkit.Chem.AllChem import CanonSmiles, MolToSmiles

from amethyst.io import parse_file_input

final_smiles = [CanonSmiles(x) for x in ["*CC", "*CCC", "*CNCC", "*CC(O)C"]]


def test_file_parser_newline():
    newline = "tests\\newline_input.txt"
    newline_subs = parse_file_input(newline)
    assert len(newline_subs.subs) == 4
    smi = [CanonSmiles(MolToSmiles(x)) for x in newline_subs.subs]
    assert smi == final_smiles


def test_file_parser_delimiter():
    comma = "tests\\comma_input.txt"
    comma_subs = parse_file_input(comma, ",")
    assert len(comma_subs.subs) == 4
    smi_comma = [CanonSmiles(MolToSmiles(x)) for x in comma_subs.subs]
    assert smi_comma == final_smiles

    semicolon = "tests\\semicolon_input.txt"
    semicolon_subs = parse_file_input(semicolon, ";")
    assert len(semicolon_subs.subs) == 4
    smi_semicolon = [CanonSmiles(MolToSmiles(x)) for x in semicolon_subs.subs]
    assert smi_semicolon == final_smiles


def test_file_parser_r_num_in_filename():
    pass


def test_mol_parser():
    pass


def test_absolute_relative_paths():
    pass


def test_isotope_num():
    # NOTE - [2*]
    pass
