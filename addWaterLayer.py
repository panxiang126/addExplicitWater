subName="sub.xsd"

#min(num, None) max(num, None) axis("x", "y", "z") cond("and", "or")
exCon=[1, 2, "z", "and"]

dirAxis=["x","y","z"]
if exCon[2] in dirAxis:
    exCon[2] = dirAxis.index(exCon[2])
else:
    raise Exception("Invalid axis name (optional: \"x\", \"y\", \"z\")", exCon[2])

import random
import numpy as np
from ase.io import read, write
from ase.data import atomic_numbers, atomic_names, covalent_radii


def getMinBondLeg(na1, na2, rTh=1.0):
    return rTh * (covalent_radii[atomic_numbers[na1]] +
                  covalent_radii[atomic_numbers[na2]])

def getInfo(struc, wrap=False):
    if wrap:
        #struc.center(vacuum=None, axis=(0, 1, 2), about=None)
        struc.wrap()
    cell = struc.cell.cellpar()
    x = [ eval(j)([ k.x for  k in struc]) for j in ["min","max"] ]
    y = [ eval(j)([ k.y for  k in struc]) for j in ["min","max"] ]
    z = [ eval(j)([ k.z for  k in struc]) for j in ["min","max"] ]
    fa = [ ((i[1]-i[0])/cell[n]) for n,i in enumerate([x,y,z]) ]
    fa_min = "x" if fa[0] < fa[1] else "y"
    fa_min = "z" if fa[0] > fa[2] and fa[1] > fa[2] else fa_min
    return( " cell: %s\n  x: %s\n  y: %s\n  z: %s\n" %(cell, x, y, z) )


struc = read(subName)
print("read struc:\n", getInfo(struc,wrap=True))

struc_tmp = struc[:0]
from ase import Atom
[ struc_tmp.append(Atom("H", [0.0, 0.0, 0.0])) for i in range(8) ]
fal=[0.001,0.999]
struc_tmp.set_scaled_positions(
    [[i,j,k] for i in fal for j in fal for k in fal]
)
pos = struc_tmp.get_positions()

max_x, max_y, max_z = [ pos[:,i].max()-pos[:,i].min() for i in range(3) ]
exTh = 1.5
th_x, th_y, th_z = max_x*exTh, max_y*exTh, max_z*exTh

print("Maximum possible range: th_x, th_y, th_z", th_x, th_y, th_z, "\n\n")

waters = read("iniWater.traj", index=":")
water  = random.choice(waters)
print("water:\n", getInfo(water))

superLa = [ int(i/water.cell.cellpar()[n])+1 for n,i in enumerate([th_x, th_y, th_z]) ]
print("supercell of water bulk:", superLa, "\n\n")
water = water * superLa

allElement = set(struc.get_chemical_symbols())
bondTh = dict([ [j+"-"+i, getMinBondLeg(i, j, rTh=1.5)] for i in allElement for j in ["H", "O"] ])
#allElement = set(struc.get_chemical_symbols()+["H", "O"])
#bondTh = dict([ [i+"-"+j, getMinBondLeg(i, j)] for i in allElement for j in allElement ])
# OO: 2.45 OH: 0.92 1.49 HH: 1.35
bondTh["O-H"],bondTh["H-O"],bondTh["H-H"],bondTh["O-O"] = 1.49*0.95,1.49*0.95,1.35*0.95,2.45*0.95
print("threshold for bond length: ", bondTh, "\n\n")

atom_O_index = [ [i.index] for i in water if i.symbol == "O" and
                                             i.x < th_x and
                                             i.y < th_y and
                                             i.z < th_z ]

print("Number of filtered oxygens:", len(atom_O_index) )

for index in range(len(atom_O_index)):
    allDis = water.get_distances(atom_O_index[index][0], range(len(water)), mic=True, vector=False)
    allDis_dic = dict(zip([i for i in range(len(allDis))], allDis))
    atom_O_index[index].append(sorted(allDis_dic.items(),key=lambda x:x[1])[1:3])


addList=[]
for index in range(len(atom_O_index[:])):
    struc_ori = struc.copy()
    h2o_index = atom_O_index[index]
    h2o = [water[h2o_index[0]],
           water[h2o_index[1][0][0]],
           water[h2o_index[1][1][0]]]
    break_flag = False
    for ni,i in enumerate(h2o):
        struc.append(i)
        struc.wrap()
        # jian ce jing zhi
        for j in struc[:-(1+ni)]:
            dl=struc.get_distance(-1, j.index, mic=True, vector=False)
            #print("   ",dl,end=" ")
            if dl < bondTh[struc[-1].symbol+"-"+j.symbol]:
                #print(f"too close", i, j, dl, bondTh[struc[-1].symbol+"-"+j.symbol])
                break_flag = True
                break
        if break_flag:
            break
    if break_flag:
        struc = struc_ori
    else:
        #print(f"add water:", h2o)
        addList.append(h2o)
    #print(len(struc))

print("Number of added water:", len(addList))

write("add_water_"+subName, struc)

