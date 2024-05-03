import pytest

from amethyst.amethyst import enumerate
from amethyst.io import parse_file_input

path = "tests\\newline_input.txt"
core = ""


def test_enumerate_file():
    enumerate(core, subs_path=path)
    pass


def test_enumerate_mol():
    pass
