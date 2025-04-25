import random

import numpy as np


class Box:
    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z

    def periodic(self, dx, dy, dz):
        if abs(dx) > self.x / 2:
            dx = dx - self.x * np.sign(dx)
        if abs(dy) > self.y / 2:
            dy = dy - self.y * np.sign(dy)
        if abs(dz) > self.z / 2:
            dz = dz - self.z * np.sign(dz)
        return dx, dy, dz


def force(r, f0=1.0, r0=1.0):
    if r < r0:
        return f0 * (1 - r / r0)
    else:
        return 0


def points_coord(N, box_x, box_y, box_z = 0.0, dt=0.1, max_shift=0.001, max_iteration=10000):
    if N < 1 or box_x < 0.0 or box_y < 0.0 or box_z < 0.0:
        raise ValueError("N or box_x or box_y or box_z are not correct")
    r0 = (box_x * box_y * box_z / N) ** (1/3)
    box = Box(box_x, box_y, box_z)

    # нач коорд всех точек
    x = (0.5 - np.random.random(N)) * box.x
    y = (0.5 - np.random.random(N)) * box.y
    z = (0.5 - np.random.random(N)) * box.z

    # red_or_blue = ['red', 'blue']
    # colors = random.choices(red_or_blue, k = N)

    # смотрим взаим точек
    max_shift = 1  # макс смещение за один шаг
    cycles = 0
    while max_shift > 0.001 and cycles < max_iteration:
        cycles += 1
        fx = np.zeros(N)
        fy = np.zeros(N)
        fz = np.zeros(N)
        for i in range(N - 1):
            for j in range(i + 1, N):
                dx = x[i] - x[j]
                dy = y[i] - y[j]
                dz = z[i] - z[j]
                dx, dy, dz = box.periodic(
                    dx, dy, dz
                )  # так как коробка замкнута как цилиндр выбираем меньшее (правильное) расстояние
                r = np.sqrt(dx**2 + dy**2 + dz**2)
                f = force(r, r0=r0) * dt
                fx[i] += f * dx / r
                fy[i] += f * dy / r
                fz[i] += f * dz / r
                fx[j] -= f * dx / r
                fy[j] -= f * dy / r
                fz[j] -= f * dz / r
        x = x + fx
        y = y + fy
        z = z + fz
        max_shift = np.max(fx)
        for i in range(N):
            x[i], y[i], z[i] = box.periodic(x[i], y[i], z[i])  # не даем вылетить за коробку
    # m = int(N*fraction)
    # best_indecies = []
    # best_distance = 0.0
    # for _ in range(N**2):
    #     indecies = random.choices(list(range(N)), k=m)
    #     sum_distance = 0.0
    #     for i in range(m-1):
    #         for j in range(i+1, m):
    #             ii = indecies[i]
    #             jj = indecies[j]
    #             dx = x[ii] - x[jj]
    #             dy = y[ii] - y[jj]
    #             dx, dy = box.periodic(
    #                 dx, dy
    #             )  # так как коробка замкнута как цилиндр выбираем меньшее (правильное) расстояние
    #             r = dx**2 + dy**2
    #             sum_distance += np.sqrt(r)
    #     if sum_distance > best_distance:
    #         best_distance = sum_distance
    #         best_indecies = indecies
    return x, y, z
