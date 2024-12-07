import sys
sys.path.append('d:/code/pywfn')
from pywfn.base import Mol
from pywfn.reader import LogReader
from pathlib import Path
import numpy as np


from pywfn.bondprop import order
root=f'paper/mols/Fig1.scan'
bonds=[(2,3),(3,10),(10,11),(10,12)]
orders=np.zeros(shape=(3,38,3))
i=0
for path in Path(root).iterdir():
    if not path.suffix=='.log':continue
    print(path)
    mol=Mol(LogReader(f'{path}'))
    caler=order.Calculator(mol)
    for b,bond in enumerate(bonds):
        result=caler.pi_smo(list(bond))
        orders[0,i,b]=result
    # 轨道挑选法
    # pocv方法计算pi键键级
    result=caler.pi_pocv()
    index=0
    for a1,a2,val in result:
        a1,a2=int(a1),int(a2)
        if (a1,a2) not in bonds:continue
        orders[1,i,index]=val
    # mayer键键级
    result=caler.mayer()
    index=0
    for a1,a2,val in result:
        a1,a2=int(a1),int(a2)
        if (a1,a2) not in bonds:continue
        orders[2,i,index]=val
    i+=1