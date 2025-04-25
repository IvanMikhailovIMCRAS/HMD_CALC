import numpy as np
import math

def random_vector(z_direction, length=1): #from lipid.py
    vector = 1.0 -2.0*np.random.random(3)
    vector[2] += z_direction
    r = np.sqrt(np.sum(vector**2))
    return vector/r*length


class Dendron():
    def __init__(self, n:int, g:int):
        self.n = n
        self.g = g
        num_bonds = (2**(g+1) - 1) * n
        connect = {i : [] for i in range(0,num_bonds+1)}
        if n < 1 or g < 0:
            raise ValueError("n or g are not correct")

        self.bonds = np.zeros(shape=(num_bonds, 2), dtype=int)
        if g == 0:
            for i in range(n):
                self.bonds[i] = [i,i+1]
                connect[i].append(i+1)
                connect[i+1].append(i)
        else:
            for i in range(n*2):
                self.bonds[i] = [i,i+1]
                connect[i].append(i+1)
                connect[i+1].append(i)
        i=n*2
        rep = 2
        while i < num_bonds:
            self.bonds[i] = [i-n*math.floor(rep/2), i+1]
            connect[i-n*math.floor(rep/2)].append(i+1)
            connect[i+1].append(i-n*math.floor(rep/2))
            i+=1
            j=1
            while j < n and i < num_bonds:
                self.bonds[i] = [i,i+1]
                connect[i].append(i+1)
                connect[i+1].append(i)
                i+=1
                j+=1
            rep+=1
        #self.bonds.astype(int)
        self.types = [4]*(1+num_bonds)
        self.coords = np.zeros(shape=(num_bonds+1,3), dtype = float)
        
        self.angles = []
        for i in range(num_bonds+1):
            c = connect[i]
            if len(c) == 2:
                self.angles.append([c[0], i, c[1]])
            if len(c) == 3:
                self.angles.append([c[0], i, c[1]])
                self.angles.append([c[0], i, c[2]])
                self.angles.append([c[1], i, c[2]])
    
    def create(self, x0, y0, z0, z_direction):
        self.coords[0,0] = x0
        self.coords[0,1] = y0
        self.coords[0,2] = z0
        for bond in self.bonds:
            self.coords[bond[1]] = self.coords[bond[0]] + random_vector(z_direction,length=0.4125 )
        return self.coords[1:,0], self.coords[1:,1], self.coords[1:,2]