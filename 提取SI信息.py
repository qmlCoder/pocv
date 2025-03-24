import sys;sys.path.append('d:/code/pywfn')
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.tools import extractSI
from pathlib import Path


root="D:/code/POCV/paper/mols/Fig3.piEle"
root="D:/code/POCV/paper/mols/Fig6.neNHC"
paths=[]
for path in Path(root).iterdir():
    if path.suffix!='.log':continue
    print(path)
    paths.append(f'{path}')

tool=extractSI.Tool(paths)
tool.save('Fig6_si.txt')