import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
sys.path.append(rf"D:\code\pywfn")
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomprop import activity,direction
# from pywfn.bondprop import piDM
from pywfn.maths import points_rotate
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from pywfn.utils import printer
printer.ifShell=False

root="D:/code/POCV/paper/mols/Fig6.neNHC"
# root=r'D:\BaiduSyncdisk\Articles\HFV\gfile\NHC\NHC-wfn'
names=['NHC-wfn','since4_s_wfn']
atoms=[1,13]
results={}
for n,name in enumerate(names):
    mn=Mol(LogReader(f'{root}/{name}-.log'))
    m0=Mol(LogReader(f'{root}/{name}0.log'))
    mp=Mol(LogReader(f'{root}/{name}+.log'))
    # mols=[mn,m0,mp]
    
    for m,mol in enumerate([mn,m0,mp]):
        caler=activity.Calculator(mol)
        dirCaler=direction.Calculator(mol)
        dirs=dirCaler.reactions(atoms[n])
        result=caler.freeValence(atoms[n],dirs)
        results[f'{n},{m}']=result