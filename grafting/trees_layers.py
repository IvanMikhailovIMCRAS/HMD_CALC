from trees import Dendron
from spread_points_3D import points_coord
import numpy as np

z0 = -5.0
z1 = -5.0
box_x, box_y, box_z = 15.0, 15.0, 15.0

N = 5 #lipid count in one layer; was 64
rho = 1 #density of particles; was 3
n_total = int(box_x*box_y*box_z*rho)
n_water = n_total - 2 * N * 11

x_head, y_head, z_head = points_coord(N, box_x, box_y)
x_head1, y_head1, z_head1 = points_coord(N, box_x, box_y) #second layer
D = Dendron(2,3)
x, y, z = [], [], []
types = []

for i in range(len(x_head)):
    D.create(x_head[i], y_head[i], z0, -2.0)
    x += list(D.coords[:,0])
    y += list(D.coords[:,1])
    z += list(D.coords[:,2])
    types += D.types
    if i == 0:
        bonds = D.bonds + 1
    else:
        bonds = np.vstack([bonds, D.bonds + 1 + i*len(D.coords)*2])
    D.create(x_head1[i], y_head1[i], z1, 2.0)
    x += list(D.coords[:,0])
    y += list(D.coords[:,1])
    z += list(D.coords[:,2])
    types += D.types
    bonds = np.vstack([bonds, D.bonds + 1 + len(D.coords) + i*len(D.coords)*2])

#xw, yw, zw = points_coord(15, box_x, box_y, box_z, max_iteration=100)

with open('dendron.ent', 'w') as f:
    n = 0
    for xt, yt, zt, t in zip(x, y, z, types):
        n += 1
        f.writelines(f'HETATM{n:5d}  {t}{n:12d}{xt:12.3f}{yt:8.3f}{zt:8.3f}\n')
    # for xw, yw, zw in zip(xw, yw, zw):
    #      n += 1
    #      f.writelines(f'HETATM{n:5d}  H{n:12d}{xw:12.3f}{yw:8.3f}{zw:8.3f}\n')
    for bond in bonds:
        f.writelines(f'CONECT{bond[0]:5d}{bond[1]:5d}\n')