# NOTE - checks for relative and absolute file paths
import pytest
from rdkit.Chem.AllChem import CanonSmiles, MolFromSmiles, MolToSmiles

from amethyst.io import parse_file_input, parse_mol_input

final_smiles = [CanonSmiles(x) for x in ["*CC", "*CCC", "*CNCC", "*CC(O)C"]]


def test_file_parser_newline():
    newline = "tests\\newline_input.txt"
    newline_subs = parse_file_input(newline, r_num=1)
    assert len(newline_subs.subs) == 4
    smi = [CanonSmiles(MolToSmiles(x)) for x in newline_subs.subs]
    assert smi == final_smiles


def test_file_parser_delimiter():
    comma = "tests\\comma_input.txt"
    comma_subs = parse_file_input(comma, r_num=1, delimiter=",")
    assert len(comma_subs.subs) == 4
    smi_comma = [CanonSmiles(MolToSmiles(x)) for x in comma_subs.subs]
    assert smi_comma == final_smiles

    semicolon = "tests\\semicolon_input.txt"
    semicolon_subs = parse_file_input(semicolon, r_num=1, delimiter=";")
    assert len(semicolon_subs.subs) == 4
    smi_semicolon = [CanonSmiles(MolToSmiles(x)) for x in semicolon_subs.subs]
    assert smi_semicolon == final_smiles


def test_file_parser_r_num_in_filename():
    pass


mols_to_parse_single_str = [
    [CanonSmiles(x) for x in ["[*:2]CCCCCC", "CCC[*:2]C", "CCC[*:2]"]]
]
mols_to_parse_mult_str = [
    [CanonSmiles(x) for x in ["[*:2]CCCCCC", "CCC[*:1]C", "CCC[*:2]"]],
    [CanonSmiles(x) for x in ["[*:2]CCCCCC", "CCC[*:2]C", "CCC[*:2]"]],
]


def test_mol_parser_single_r_str():
    result = parse_mol_input(mols_to_parse_single_str)
    subs_1 = result[0]
    assert len(result) == 1
    assert subs_1.r_num == 1
    smi_result = [CanonSmiles(MolToSmiles(x)) for x in subs_1.subs]
    assert smi_result == [
        CanonSmiles(x) for x in ["[*:1]CCCCCC", "CCC[*:1]C", "CCC[*:1]"]
    ]


def test_mol_parser_multiple_r_str():
    results = parse_mol_input(mols_to_parse_mult_str)
    assert len(results) == 2
    subs_1 = results[0]
    assert subs_1.r_num == 1
    subs_2 = results[1]
    assert subs_2.r_num == 2
    smi_1_result = [CanonSmiles(MolToSmiles(x)) for x in subs_1.subs]
    smi_2_result = [CanonSmiles(MolToSmiles(x)) for x in subs_2.subs]
    assert smi_1_result == [
        CanonSmiles(x) for x in ["[*:1]CCCCCC", "CCC[*:1]C", "CCC[*:1]"]
    ]
    assert smi_2_result == [
        CanonSmiles(x) for x in ["[*:2]CCCCCC", "CCC[*:2]C", "CCC[*:2]"]
    ]


def test_mol_parser_single_r_mol():
    mols = [[MolFromSmiles(x) for x in mols_to_parse_single_str[0]]]
    result = parse_mol_input(mols)
    subs_1 = result[0]
    assert len(result) == 1
    assert subs_1.r_num == 1
    smi_result = [CanonSmiles(MolToSmiles(x)) for x in subs_1.subs]
    assert smi_result == [
        CanonSmiles(x) for x in ["[*:1]CCCCCC", "CCC[*:1]C", "CCC[*:1]"]
    ]


def test_mol_parser_multiple_r_mol():
    mols = []
    for i in mols_to_parse_mult_str:
        temp = []
        for j in i:
            temp.append(MolFromSmiles(j))
        mols.append(temp)
    result = parse_mol_input(mols)
    assert len(result) == 2
    subs_1 = result[0]
    assert subs_1.r_num == 1
    subs_2 = result[1]
    assert subs_2.r_num == 2
    smi_1_result = [CanonSmiles(MolToSmiles(x)) for x in subs_1.subs]
    smi_2_result = [CanonSmiles(MolToSmiles(x)) for x in subs_2.subs]
    assert smi_1_result == [
        CanonSmiles(x) for x in ["[*:1]CCCCCC", "CCC[*:1]C", "CCC[*:1]"]
    ]
    assert smi_2_result == [
        CanonSmiles(x) for x in ["[*:2]CCCCCC", "CCC[*:2]C", "CCC[*:2]"]
    ]


def test_absolute_relative_paths():
    pass


def test_isotope_num():
    # NOTE - [2*]
    pass
