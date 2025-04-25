from scan_track import ReadTrack
from trees import Dendron
import numpy as np
import random



if __name__ == '__main__':
    num_grafting = 100
    track = ReadTrack('')
    dendron = Dendron(n=2, g=0)
    while track.one_step():
        pass
    z0 = np.mean(track.z[track.btype==1])
    num_lipid = len(track.z[track.btype==1])//3
    num_bot = []
    num_top = []
    for i in range(0, num_lipid*11, 11):
        if track.z[i]<z0:
            num_bot.append(i)
        else:
            num_top.append(i)
    #print(len(num_bottom), len(num_top))
    list_top = random.sample(population=num_top, k=num_grafting)    
    list_bot = random.sample(population=num_bot, k=num_grafting) 
    x = track.x[:num_lipid*11]
    y = track.y[:num_lipid*11]
    z = track.z[:num_lipid*11]
    n1 = len(x)
    for i in list_top:
        xt, yt, zt = dendron.create(x[i],y[i],z[i], z_direction=1.0)
        x = np.hstack([x, xt])
        y = np.hstack([y, yt])
        z = np.hstack([z, zt])
    for i in list_bot:
        xt, yt, zt = dendron.create(x[i],y[i],z[i], z_direction=-1.0)
        x = np.hstack([x, xt])
        y = np.hstack([y, yt])
        z = np.hstack([z, zt])
    delta_n = len(x)-n1
    x = np.hstack([x, track.x[num_lipid*11:track.num_atoms+1-delta_n]])
    y = np.hstack([y, track.y[num_lipid*11:track.num_atoms+1-delta_n]])
    z = np.hstack([z, track.y[num_lipid*11:track.num_atoms+1-delta_n]])
    track.btype[num_lipid*11:num_lipid*11+delta_n+1]=4
    
    
#xt, yt, zt = dendron.create(x[i],y[i],z[i], z_direction=1.0)
#print(xt)    