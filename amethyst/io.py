import os
import re
from dataclasses import dataclass
from typing import List, Optional, Union

from loguru import logger
from rdkit.Chem.AllChem import Mol, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdRGroupDecomposition import RelabelMappedDummies, RGroupLabelling

from amethyst.utils import mols_to_str

logger.add(sink="io.log", level=10)


@dataclass
class Substituents:
    r_num: int
    subs: List[Mol]


def parse_file_input(path: str, delimiter: Optional[str], r_num: Optional[int] = None):
    if os.path.isfile(str):
        logger.debug("File path is good.")
        pass
    else:
        logger.error(f"{path} is not a file!")
        raise (FileNotFoundError(f"{path} is not a file!"))

    subs = Substituents()

    if r_num is not None:
        subs.r_num = r_num
    else:
        path_split = path.split("/").split("\\")
        m = re.match("[0-9]+", path_split[-1])
        subs.r_num = m.group(0)

    logger.debug(f"File give for R{r_num}.")

    subs_list: List[Mol] = []

    with open(path, "r") as file:
        if delimiter is None:
            for i in file:
                subs_list.append(MolFromSmiles(i))
                logger.debug(f"SMILES added: {i}")
        else:
            lines: List[str] = file.readlines()
            for i in lines:
                split_lines: List[str] = lines[i].split(
                    delimiter
                )  # idk if that works lol
                mols: List[Mol] = [MolFromSmiles(x) for x in split_lines]
                subs_list.append(*mols)
                logger.debug(f"SMILES added: {*split_lines,}")

    logger.debug(f"Final sub list: {mols_to_str(subs_list)}")

    subs.subs = list(
        map(
            lambda x: RelabelMappedDummies(x, outputLabels=RGroupLabelling.AtomMap),
            subs_list,
        )
    )

    return subs


def parse_string_input():
    pass
