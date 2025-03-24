import sys
sys.path.append('d:/code/pywfn')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pywfn.base import Mol
from pywfn.reader import LogReader
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy.typing as npt
from typing import Any


### Figure 2
from pywfn.bondprop import order
root=f'paper/mols/Fig1.scan'
bonds=[(2,3),(3,10),(10,11),(10,12)]
orders=np.zeros(shape=(37,3,4))
orders[:,0,:]=np.loadtxt('paper/npys/piOrderSMO.txt')
for i in range(37):
    path=f'{root}/f{i+1:0>2}.log'
    print(path)
    mol=Mol(LogReader(f'{path}'))
    caler=order.Calculator(mol)
    # 轨道挑选法
    # print(f"{'pi bond order (SMO)':-^40}")
    # for b,(a1,a2) in enumerate(bonds): # spepd is too slowly
    #     result=caler.pi_smo([a1,a2])
    #     orders[i,0,b]=result
    #     print(a1,a2,result)
    
    
    for b,(a1,a2) in enumerate(bonds):
        val=orders[i,0,b]
        print(a1,a2,val)
    # pocv方法计算pi键键级
    print(f"{'pi bond order (POCV)':-^40}")
    result=caler.pi_pocv()
    index=0
    for a1,a2,val in result:
        a1,a2=int(a1),int(a2)
        if (a1,a2) not in bonds:continue
        print(a1,a2,val)
        orders[i,1,index]=val
        index+=1
    # mayer键键级
    print(f"{'mayer bond order':-^40}")
    result=caler.mayer()
    index=0
    for a1,a2,val in result:
        a1,a2=int(a1),int(a2)
        if (a1,a2) not in bonds:continue
        print(a1,a2,val)
        orders[i,2,index]=val
        index+=1


### Figure 3 pi electron population
print(f"{'Data of Figure 3':-^40}")
from pywfn.atomprop import charge
root='paper/mols/Fig2.piEle'
paths=[
    f'{root}/Ph2.log',
    f'{root}/CH2CHCHO.log',
    f'{root}/CO2.log',
    f'{root}/NO2.log',
    f'{root}/NNO.log',
    f'{root}/BF3.log',
]
mols=[Mol(LogReader(path)) for path in paths]

for m,mol in enumerate(mols):
    if m in [2,4]:
        mol.atom(1)._props['normal']=np.array([1.,0.,0.])
        mol.atom(2)._props['normal']=np.array([1.,0.,0.])
        mol.atom(3)._props['normal']=np.array([1.,0.,0.])
    caler=charge.Calculator(mol)
    elects=caler.piElectron('mulliken')
    print(f'{mol.reader.path}')
    for atm,dx,dy,dz,val in elects:
        if mol.atom(int(atm)).symbol=='H':continue
        print(f'{atm:>3.0f}{val:>10.4f}')
    print('total:',np.sum(elects[:,-1]))

### Figure 4
print(f"{'Data of Figure 4':-^40}")
from pywfn.bondprop import order
def get_dirs(angles):
    dirs:list[np.ndarray]=[]
    for a in angles:
        z=0
        y=np.sin(a)
        x=np.cos(a)
        dirs.append(np.array([x,y,z]))
    return dirs
angs=np.linspace(0,np.pi*2,45)
dirs=get_dirs(angs)
dirs=np.array(dirs)
data=np.zeros(shape=(2,len(dirs))) #键，方向，（键级1，键级2）
bonds=[
    [2,1],[1,5]
]
mol=Mol(reader=LogReader('paper/mols/bondOrderDirs/R-C3.log'))
caler=order.Calculator(mol)
for bond in bonds:
    orders=caler.dirMayer(bond,dirs)
    for a1,a2,dx,dy,dz,val in orders:
        print(f'{a1:>3.0f}{a2:>3.0f}{dx:>8.4f}{dy:>8.4f}{dz:>8.4f}{val:>10.4f}')

### Figure 5 Aromacity representation by standard deviation of pi Bond Orders
print(f"{'Data of Figure 5':-^40}")
from pywfn.bondprop import order
from pywfn.molprop import aromatic

root='paper/mols/aromacity'
paths=[
    f'{root}/C6H6.log',
    f'{root}/ph2.log',
    f'{root}/line6.log',
    f'{root}/rNH2.log',
    f'{root}/rNF2.log',
    f'{root}/C4H4.log',
    f'{root}/ph-ph.log',
    f'{root}/phR-phR.log',
    f'{root}/since.log',
]
mols=[Mol(LogReader(path)) for path in paths]
for m,mol in enumerate(mols):
    caler=order.Calculator(mol)
    arom=aromatic.Calculator(mol)
    print(f'{mol.reader.path}')
    print(arom.pisd())

### Figure 6
print(f"{'Data of Figure 6':-^40}")
from pywfn.atomprop import direction
from pywfn.atomprop import activity

root='paper/mols/neNHC'
paths=[
    f'{root}/NHC-wfn-.log',
    f'{root}/NHC-wfn0.log',
    f'{root}/NHC-wfn+.log',
    f'{root}/since4_s_wfn-.log',
    f'{root}/since4_s_wfn0.log',
    f'{root}/since4_s_wfn+.log',
]
atms=[1,1,1,13,13,13]
for i,path in enumerate(paths):
    mol=Mol(LogReader(path))
    valCaler=activity.Calculator(mol)
    dirCaler=direction.Calculator(mol)
    dirs=dirCaler.reactions(atms[i])
    vals=valCaler.freeValence(atms[i],dirs)
    print(path)
    for atm,x,y,z,val in vals:
        print(f'{atm:>3.0f}{x:>10.4f}{y:>10.4f}{z:>10.4f}{val:>10.4f}')

### Figure 7 Freevalece of different reaction sites
print(f"{'Data of Figure 7':-^40}")
from pywfn.atomprop import activity
from pywfn.atomprop import direction
path='paper/mols/freeValence/M4_wfn0.log'
mol=Mol(LogReader(path))
dirCaler=direction.Calculator(mol)
valCaler=activity.Calculator(mol)
atms=[3,25,26,27]
for atm in atms:
    dirs=dirCaler.reactions(atm)
    vals=valCaler.freeValence(atm,dirs)
    print(atm)
    for atm,x,y,z,val in vals:
        print(f'{atm:>3.0f}{x:>10.4f}{y:>10.4f}{z:>10.4f}{val:>10.4f}')