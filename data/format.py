#! /usr/bin/env python

import numpy as np

fout = open('formatted/DY.txt','w')

with open('DY_v3.txt') as fin:
    for line in fin:
        sp = line.strip().split()
        if len(sp)!=11:
            continue
        q = float(sp[0])
        m = 105.
        x = 33.
        y = float(sp[9])
        z = float(sp[8])

        th = float(sp[5])
        thW = float(sp[6])
        thV = float(sp[7])
        p = float(sp[10])
        
        px = p*np.cos(th)
        py = p*np.cos(th)*np.tan(thV)
        pz = p*np.cos(th)*np.tan(thW)

        fout.write('{0:f}\t{1:f}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\n'.format(q,m,x,y,z,px,py,pz))

        #print p, np.sqrt(px**2+py**2+pz**2)

fout.close()
