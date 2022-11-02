subName="sub.xsd"
#min(0<num<1, None) max(0<num<1, None) axis("x", "y", "z") cond("and", "or")
exCon=[0.9, 0.7, "z", "and"]

if exCon[2] not in ["x","y","z"]:
    raise Exception("Invalid axis name (optional: \"x\", \"y\", \"z\")", exCon[2])

axisNum=str(["x","y","z"].index(exCon[2]))
if exCon.count(None) == 0:
    ifcon = "struc[-1]." + exCon[2] + " < " + str(exCon[0]) + "*struc.cell.cellpar()[" + axisNum + "]" \
            + " " + exCon[3] + " " +                                                                   \
            "struc[-1]." + exCon[2] + " > " + str(exCon[1]) + "*struc.cell.cellpar()[" + axisNum + "]"
elif exCon.count(None) == 1:
    if exCon[1]  == None:
        ifcon = "struc[-1]." + exCon[2] + " < " + str(exCon[0]) + "*struc.cell.cellpar()[" + axisNum + "]"
    else:
        ifcon = "struc[-1]." + exCon[2] + " > " + str(exCon[1]) + "*struc.cell.cellpar()[" + axisNum + "]"
else:
    ifcon = "False"

import random
import numpy as np
from ase.io import read, write
from ase.data import atomic_numbers, atomic_names, covalent_radii
import time,sys


def printOut(po):
    with open("run_"+time.strftime("%Y%m%d%H%M%S")+"_log.txt", "a+") as fw:
        print(po, file=fw)
    #print(po)


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

import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
fig,ax = plt.subplots(1,3,dpi=300,figsize=(10, 7.65))


fig.text(0.1,0.5,"sub + water -> aimStr",fontsize=12, verticalalignment="center", horizontalalignment="center")
fig.tight_layout()
for i in range(3):
  ax[i].axis("off")


struc = read(subName)
print("read struc:\n", getInfo(struc,wrap=True))

plot_atoms(struc,ax[0],radii=0.5,rotation=("-90x,0y,0z"))


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

plot_atoms(water,ax[1],radii=0.3,rotation=("0x,0y,0z"))

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


print("Collecting water...")
O_num = len(atom_O_index[:])
for index in range(len(atom_O_index)):
    rp = index/O_num
    print("\r"+f"{rp*100:.0f}%: [ ","#"*int(rp*50)," "*(50-int(rp*50)) ,sep="", end="]")
    sys.stdout.flush()
    allDis = water.get_distances(atom_O_index[index][0], range(len(water)), mic=True, vector=False)
    allDis_dic = dict(zip([i for i in range(len(allDis))], allDis))
    atom_O_index[index].append(sorted(allDis_dic.items(),key=lambda x:x[1])[1:3])
print()


randomOffset = [random.random()*struc.cell.cellpar()[i] for i in  range(3)]
print("water offset:", randomOffset, "\n\n")

print("prohibited area:", ifcon.replace("struc[-1].","").replace("False","None").replace("struc.cell.cellpar()","cell"), "\n\n")



print("Adding water...")
addList=[]
O_num = len(atom_O_index[:])
for index in range(O_num):
    rp = index/O_num
    print("\r"+f"{rp*100:.0f}%: [ ","#"*int(rp*50)," "*(50-int(rp*50)) ,sep="", end="]")
    sys.stdout.flush()
    struc_ori = struc.copy()
    h2o_index = atom_O_index[index]
    h2o = [water[h2o_index[0]],
           water[h2o_index[1][0][0]],
           water[h2o_index[1][1][0]]]
    for i in h2o:
        i.position += randomOffset
    break_flag = False
    for ni,i in enumerate(h2o):
        struc.append(i)
        struc.wrap()
        if eval(ifcon):
            printOut(f"s{index}  exceeding the threshold: {struc[-1]}")
            break_flag = True
            break
        for j in struc[:-(1+ni)]:
            dl=struc.get_distance(-1, j.index, mic=True, vector=False)
            if dl < bondTh[struc[-1].symbol+"-"+j.symbol]:
                printOut(f"s{index}  too close: {i} {j} {dl}")
                break_flag = True
                break
        if break_flag:
            break
    if break_flag:
        struc = struc_ori
    else:
        printOut(f"s{index}  add water: {h2o}")
        addList.append(h2o)
print()

print("Number of added water:", len(addList))
plot_atoms(struc,ax[2],radii=0.3,rotation=("-90x,0y,0z"))

write("add_water_"+subName, struc)
fig.savefig("add_water_"+subName+".png")
