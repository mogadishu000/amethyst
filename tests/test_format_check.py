import logging

import pytest

from amethyst.io import rgroup_format

with open("tests/smi_atom_map.txt", "r") as f:
    atom_map = f.read()
    atom_map = atom_map.split("\n")

with open("tests/smi_isotope.txt", "r") as f:
    isotope = f.read()
    isotope = isotope.split("\n")

with open("tests/smi_mixed_inline.txt", "r") as f:
    mixed_inline = f.read()
    mixed_inline = mixed_inline.split("\n")

with open("tests/smi_mixed_separate.txt", "r") as f:
    mixed_separate = f.read()
    mixed_separate = mixed_separate.split("\n")


def test_formats(caplog):
    caplog.set_level(logging.INFO)
    assert [rgroup_format(smi) for smi in atom_map] == [1] * len(atom_map)
    assert [rgroup_format(smi) for smi in isotope] == [2] * len(isotope)
    assert [rgroup_format(smi, safe_mode=True) for smi in mixed_inline] == [3, 3]
    assert [rgroup_format(smi) for smi in mixed_separate] == [2, 2, 1, 1]
    with pytest.raises(ValueError, match="No dummy atoms found"):
        rgroup_format("AAAAA")
