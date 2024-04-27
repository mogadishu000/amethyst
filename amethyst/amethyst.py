from typing import Dict, List, Optional, Union

import click
from rdkit.Chem.AllChem import Mol, MolFromSmiles
from rdkit.Chem.rdRGroupDecomposition import RelabelMappedDummies, RGroupLabelling

from amethyst.io import (
    Substituents,
    parse_file_input,
    parse_mol_input,
)
from amethyst.substitution import general_sub, placeholder_atom_sub
from loguru import logger

logger.add("amethyst_main.log", level=10)


def enumerate(
    core: Union[str, Mol],
    subs_path: Optional[str] = None,
    delimiter: Optional[str] = None,
    r_num: Optional[int] = None,
    subs_mol: Optional[
        Union[List[List[Mol]], Dict[str, List[Mol]]]
    ] = None,  # what a fucking abomination
    enantiomers: bool = False,
):
    if subs_path is not None and subs_mol is not None:
        logger.error("Both sources for R-groups passed")
        raise ValueError("Only one source of R-groups must be passed")
    elif subs_path is not None:
        logger.debug("Filepath passed")
        r_groups: Substituents = parse_file_input(subs_path, delimiter, r_num)
    elif subs_mol is not None:
        logger.debug("Mol list passed")
        r_groups: Substituents = parse_mol_input()  # NYI
    else:
        logger.error("No R-groups were passed")
        raise ValueError("No R-groups passed!")

    if type(core) is not Mol:
        core = RelabelMappedDummies(
            MolFromSmiles(core), outputLabels=RGroupLabelling.AtomMap
        )
