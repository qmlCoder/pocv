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
# print(f"{'Data of Figure 2':-^40}")
# piOrdersP=np.load('数据/piOrders-lianben.npy')
# piOrdersS=np.load('数据/SpiOrders-lianben.npy')
# myOrders =np.load('数据/myOrders-lianben.npy')
# np.loadtxt
# def dang(x):
#     if x<-180:
#         x+=360
#     if x>180:
#         x-=360
#     return x
# xs=np.linspace(-180,180,37)+90
# angs=[dang(float(x)) for x in xs]
# xticks=[f'{ang}°' for ang in angs]
# fig: Figure
# ax: Any
# fig, ax = plt.subplots(2,2,figsize=(10,6))  # type: ignore
# plt.rc('font',family='Times New Roman', size=13)
# plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.3)
# ax[0,0].axis('off')
# ax[0,1].plot(xs,piOrdersS[:,0],color='#e96d53',marker='.',linestyle='-',label='$P_{3,12}$')
# ax[0,1].plot(xs,piOrdersS[:,1],color='#2b9d8b',marker='.',linestyle='-',label='$P_{3,2}$')
# ax[1,0].plot(xs,piOrdersP[:,0],color='#e96d53',marker='.',linestyle='-',label="$P'_{3,12}$")
# ax[1,0].plot(xs,piOrdersP[:,1],color='#2b9d8b',marker='.',linestyle='-',label="$P'_{3,2}$")
# ax[1,1].plot(xs,myOrders[:,0]-1,color='#e96d53',marker='.',linestyle='-',label='$M_{3,12}-1$')
# ax[1,1].plot(xs,myOrders[:,1]-1,color='#2b9d8b',marker='.',linestyle='-',label='$M_{3,2}-1$')
# step=9
# def setFig(ax,label,title):
#     ax.set_xticks(xs[::step],labels=xticks[::step])
#     ax.legend(loc='upper left',prop={'size':10},bbox_to_anchor=(1.0, 1.0))
#     ax.set_xlabel('Dihedral Angle [°]')
# setFig(ax[0,1],'(b)','$\pi$ Bond Order cumputed by POCV method')
# setFig(ax[1,0],'(c)','Mayer Bond Order')
# setFig(ax[1,1],'(d)','Mayer Bond Order minus one')
# plt.subplots_adjust(hspace=0.4)
# fig.savefig('paper/Figure2.png',bbox_inches='tight',dpi=300)
# for each in piOrdersP:
#     print(each)
# for each in piOrdersS:
#     print(each)
# for each in myOrders:
#     print(each)


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