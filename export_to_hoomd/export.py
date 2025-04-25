from typing import Dict

import gsd
import gsd.hoomd
import numpy as np

from .scan_angles import read_angles
from .scan_bonds import read_bonds
from .scan_coord import ScanCoord


def gsd_compile(
    bead_types: Dict, path: str = "", name_output: str = "init.gsd"
) -> None:
    Coords = ScanCoord(path=path)
    coord = np.vstack([Coords.x, Coords.y, Coords.z]).T
    bonds = np.array(read_bonds(path=path)) - 1
    angls = np.array(read_angles(path=path)) - 1
    btype = Coords.btype
    snapshot = gsd.hoomd.Frame()
    snapshot.particles.N = len(coord)
    snapshot.configuration.box = [Coords.box[0], Coords.box[1], Coords.box[2], 0, 0, 0]
    snapshot.bonds.N = len(bonds)
    snapshot.angles.N = len(angls)

    print(f"particles: {snapshot.particles.N}")
    print(f"bonds: {snapshot.bonds.N}")
    print(f"angles: {snapshot.angles.N}")
    if set(btype.tolist()) != set(bead_types):
        raise ValueError("Error: Not all types ID are defined!")
    snapshot.particles.types = [bead_types[k] for k in sorted(list(bead_types))]
    print(f"bead types: {snapshot.particles.types}")

    snapshot.particles.typeid = Coords.btype
    snapshot.particles.position = coord
    snapshot.particles.mass = np.array([1.0] * snapshot.particles.N)
    if snapshot.bonds.N > 0:
        snapshot.bonds.group = bonds
    if snapshot.angles.N > 0:
        snapshot.angles.group = angls

    b_types = set()
    for b in bonds:
        if btype[b[0]] < btype[b[1]]:
            b_types.add(bead_types[btype[b[0]]] + bead_types[btype[b[1]]])
        else:
            b_types.add(bead_types[btype[b[1]]] + bead_types[btype[b[0]]])

    snapshot.bonds.types = sorted(list(b_types))
    print(f"bond types: {snapshot.bonds.types}")

    ang_types = set()
    for ang in angls:
        if btype[ang[0]] < btype[ang[2]]:
            ang_types.add(
                bead_types[btype[ang[0]]]
                + bead_types[btype[ang[1]]]
                + bead_types[btype[ang[2]]]
            )
        else:
            ang_types.add(
                bead_types[btype[ang[2]]]
                + bead_types[btype[ang[1]]]
                + bead_types[btype[ang[0]]]
            )

    snapshot.angles.types = sorted(list(ang_types))
    print(f"angle types: {snapshot.angles.types}")

    with gsd.hoomd.open(name=name_output, mode="w") as f:
        f.append(snapshot)
