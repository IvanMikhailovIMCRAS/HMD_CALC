from typing import Tuple

import numpy as np


class Box:
    """
    Sets up a 3D-box for simulation of a molecular system
    Its center has coordinates (0.0, 0.0, 0.0)
    Edge sizes: self.x, self.y, self.z
    """

    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z

    @property
    def volume(self) -> float:
        return self.x * self.y * self.z

    def periodic_correct(
        self, xb: float, yb: float, zb: float
    ) -> Tuple[float, float, float]:
        """
        Returns the coordinates to the box if they go out of it
        Args:
            xb (float), yb (float), zb (float): the beads coordinates
        Returns:
            Tuple[float, float, float]: changed (or old) coordinates
        """
        xb = Box.periodic(xb, self.x)
        yb = Box.periodic(yb, self.y)
        zb = Box.periodic(zb, self.z)
        return xb, yb, zb

    def check_in_box(self, xb: float, yb: float, zb: float) -> bool:
        """
        Checks if the bead with coordinates (xb, yb, zb) is in the box
        Args:
            xb (float), yb (float), zb (float): the beads coordinates
        Returns:
            bool: True if the bead is in the box, False otherwise
        """
        return all(
            [abs(xb) < 0.5 * self.x, abs(yb) < 0.5 * self.y, abs(zb) < 0.5 * self.z]
        )

    @staticmethod
    def periodic(coord: float, box: float) -> float:
        """
        Returns the coordinate to the box if it goes out of it
        Args:
            coord (float): coordinate
            box (float): box size along some axis
        Returns:
            float: changed (or old) coordinate
        """
        if abs(coord) > 1.5 * box:
            raise ValueError("Out of box")
        if abs(coord) > 0.5 * box:
            return coord - np.sign(coord) * box
        return coord


if __name__ == "__main__":
    box = Box(x=10.0, y=20.0, z=30.0)
    x, y, z = box.periodic_correct(2.0, 12.0, -16.0)
    print(x, y, z)
