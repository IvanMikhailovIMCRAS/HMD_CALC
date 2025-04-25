from scan_track import ReadTrack
from trees import Dendron
import numpy as np
import random
from scan_track import read_bonds


if __name__ == '__main__':
    num_grafting = 50 #num of polymers on ONE side
    track = ReadTrack('')
    dendron = Dendron(n=2, g=1)
    while track.one_step():
        pass
    z0 = np.mean(track.z[track.btype==1])
    num_lipid = len(track.z[track.btype==1])//3
    num_bot = []
    num_top = [] #список первых(нулевых) липдных бидов верхнего слоя
    for i in range(1, num_lipid*11+1, 11):
        if track.z[i]<z0:
            num_bot.append(i)
        else:
            num_top.append(i)
    list_top = random.sample(population=num_top, k=num_grafting)    
    list_bot = random.sample(population=num_bot, k=num_grafting) 
    x = track.x[:num_lipid*11] #в трэке липиды записываются первыми
    y = track.y[:num_lipid*11]
    z = track.z[:num_lipid*11]
    n1 = len(x) #количество липидных бидов
    dendron_bonds = []
    first = True
    dendron_angles = []
    for i in list_top:
        xt, yt, zt = dendron.create(x[i],y[i],z[i], z_direction=1.0)
        if first==True:
            dendron_bonds = [i, len(x)+1]
            dendron_bonds = np.vstack([dendron_bonds, dendron.bonds[1:]+len(x)])
            if len(dendron.angles) > 0:
                dendron_angles = [i, len(x)+1, len(x)+2]
                if len(dendron.angles) > 1:
                    dendron_angles = np.vstack([dendron_angles, np.asarray(dendron.angles[1:])+len(x)])
            first = False
        else:
            dendron_bonds = np.vstack([dendron_bonds, [i, len(x)+1]])
            dendron_bonds = np.vstack([dendron_bonds, dendron.bonds[1:]+len(x)])
            if len(dendron.angles) > 0:
                dendron_angles = np.vstack([dendron_angles,[i, len(x)+1, len(x)+2]])
                if len(dendron.angles) > 1:
                    dendron_angles = np.vstack([dendron_angles, np.asarray(dendron.angles[1:])+len(x)])
        x = np.hstack([x, xt])
        y = np.hstack([y, yt])
        z = np.hstack([z, zt])
    for i in list_bot:
        xt, yt, zt = dendron.create(x[i],y[i],z[i], z_direction=-1.0)
        dendron_bonds = np.vstack([dendron_bonds, [i, len(x)+1]])
        dendron_bonds = np.vstack([dendron_bonds, dendron.bonds[1:]+len(x)])
        if len(dendron.angles) > 0:
            dendron_angles = np.vstack([dendron_angles,[i, len(x)+1, len(x)+2]])
            if len(dendron.angles) > 1:
                dendron_angles = np.vstack([dendron_angles, np.asarray(dendron.angles[1:])+len(x)])
        x = np.hstack([x, xt])
        y = np.hstack([y, yt])
        z = np.hstack([z, zt])
    delta_n = len(x)-n1 #количество бидов дендронов
    x = np.hstack([x, track.x[num_lipid*11:track.num_atoms-delta_n]])
    y = np.hstack([y, track.y[num_lipid*11:track.num_atoms-delta_n]])
    z = np.hstack([z, track.z[num_lipid*11:track.num_atoms-delta_n]])
    track.btype[num_lipid*11:num_lipid*11+delta_n]=4

    all_bonds = np.vstack([read_bonds(''), dendron_bonds])
    with open('BONDS_exp', 'w') as f:
        f.writelines(f' num_bonds{len(all_bonds):13d}  num_atoms{len(x):13d}  BOX{track.box.x:14.7f}{track.box.y:17.7f}{track.box.z:17.7f}\n')
        for bond in all_bonds:
            f.writelines(f'{bond[0]:12d}{bond[1]:12d}\n')
   
    if len(dendron.angles) > 0:
        all_angles = []
        path = 'ANGLS'
        try:
            file_angles = open(path, "r")
            title = file_angles.readline().split()
            num_angles = int(title[1])
            for _ in range(num_angles):
                all_angles.append(list(map(int, file_angles.readline().split())))
            file_angles.close()
        except:
            loging = open("Warning.log", "w")
            loging.write("File BONDS was not read")
            loging.close()
        
        all_angles = np.vstack([all_angles, dendron_angles])    
        
        with open('ANGLS_exp', 'w') as f:
            f.writelines(f' num_angl{(len(all_angles)):13d}\n')
            for angle in all_angles:
                f.writelines(f'{angle[0]:12d}{angle[1]:12d}{angle[2]:12d}\n')

    with open('COORD_exp', 'w') as f:
        f.writelines(f'num_atoms{len(x):7d} box_size{track.box.x:15.10f}{track.box.y:17.10f}{track.box.z:17.10f}\n')
        n = 0
        for xt, yt, zt, t in zip(x, y, z, track.btype):
            n += 1
            xt, yt, zt = track.box.periodic_correct(xt,yt,zt)
            f.writelines(f'{n:12d}  {xt:13.10f}  {yt:13.10f}  {zt:13.10f} {t}\n')

    with open('FIXED_exp', 'w') as f:
        f.writelines(f'num_fixed{n1:7d}\n')
        for bead in range(1,n1+1):
            f.writelines(f'{bead:12d}\n')
           