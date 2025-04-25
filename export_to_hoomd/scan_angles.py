import os
from typing import List, Tuple

import numpy as np


def read_angles(path: str = "") -> List[Tuple[int, int, int]]:
    path = os.path.join(path, "ANGLS")
    angls = list()
    try:
        file_angls = open(path, "r")
        title = file_angls.readline().split()
        num_angls = int(title[1])
        for _ in range(num_angls):
            angls.append(tuple(map(int, file_angls.readline().split())))
        file_angls.close()
    except:
        loging = open("Warning.log", "a")
        loging.write("File ANGLSS was not read")
        loging.close()
    return angls
