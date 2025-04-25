import os
from typing import List, Tuple

import numpy as np


def read_bonds(path: str = "") -> List[Tuple[int, int]]:
    path = os.path.join(path, "BONDS")
    bonds = list()
    try:
        file_bonds = open(path, "r")
        title = file_bonds.readline().split()
        num_bonds = int(title[1])
        for _ in range(num_bonds):
            bonds.append(tuple(map(int, file_bonds.readline().split())))
        file_bonds.close()
    except:
        loging = open("Warning.log", "w")
        loging.write("File BONDS was not read")
        loging.close()
    return bonds


if __name__ == "__main__":
    bonds = read_bonds("data")
    print(max(bonds))
    print(max(max(bonds)))
