import os

import numpy as np


class ScanCoord:

    def __init__(self, path: str = "") -> None:
        path = os.path.join(path, "COORD")
        try:
            file = open(path, "r")
        except:
            raise FileNotFoundError(f"{str(path)} is not found")
        try:
            title = file.readline().split()
            self.num_atoms = int(title[1])
            self.box = (float(title[3]), float(title[4]), float(title[5]))
            self.x = np.zeros(self.num_atoms, dtype=float)
            self.y = np.zeros(self.num_atoms, dtype=float)
            self.z = np.zeros(self.num_atoms, dtype=float)
            self.btype = np.zeros(self.num_atoms, dtype=int)
            for i in range(self.num_atoms):
                record = file.readline().split()
                self.x[i] = float(record[1])
                self.y[i] = float(record[2])
                self.z[i] = float(record[3])
                self.btype[i] = int(record[4]) - 1
            file.close()
        except:
            file.close()
            raise FileNotFoundError(f"{str(path)} is not correct")


if __name__ == "__main__":
    coord = ScanCoord("data")
    print(coord.box, len(coord.x), coord.btype[-1])
