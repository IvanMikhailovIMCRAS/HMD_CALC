import os

import numpy as np

from .periodic_box import Box


class ReadTrack:
    def __init__(self, path=""):
        self.path = os.path.join(path, "TRACK")
        self.time_step = 0
        try:
            self.file_track = open(self.path, "r")
        except:
            raise FileNotFoundError(f"{str(self.path)} is not found")
        try:
            title = self.file_track.readline().split()
            self.num_atoms = int(title[1])
            self.box = Box(float(title[3]), float(title[4]), float(title[5]))
        except:
            self.file_track.close()
            raise FileNotFoundError(f"{str(self.path)} is not correct")
        self.x = np.zeros(self.num_atoms, dtype=float)
        self.y = np.zeros(self.num_atoms, dtype=float)
        self.z = np.zeros(self.num_atoms, dtype=float)
        self.btype = np.zeros(self.num_atoms, dtype=int)
        self.iterator = list(range(self.num_atoms))

    def one_step(self):
        try:
            title = self.file_track.readline().split()
            self.time_step = int(title[1])
            for i in self.iterator:
                record = self.file_track.readline().split()
                self.x[i] = float(record[1])
                self.y[i] = float(record[2])
                self.z[i] = float(record[3])
                self.btype[i] = int(record[4])
            return True
        except:
            self.file_track.close()
            return False


def read_bonds(path):
    path = os.path.join(path, "BONDS")
    bonds = list()
    try:
        file_bonds = open(path, "r")
        title = file_bonds.readline().split()
        num_bonds = int(title[1])
        for _ in range(num_bonds):
            bonds.append(list(map(int, file_bonds.readline().split())))
        file_bonds.close()
    except:
        loging = open("Warning.log", "w")
        loging.write("File BONDS was not read")
        loging.close()
    return bonds


if __name__ == "__main__":
    RT = ReadTrack("trajectory_sample")
    while RT.one_step():
        print(RT.time_step)
        x, y, z = [], [], []
        for i in range(RT.num_atoms):
            if RT.btype[i] == 1:
                x.append(RT.x[i])
                y.append(RT.y[i])
                z.append(RT.z[i])
        x, y, z = np.array(x), np.array(y), np.array(z)

    read_bonds("trajectory_sample")
