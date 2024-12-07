import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import sys

tree=ET.parse("D:\BaiduSyncdisk\Articles\HFV\cdxml\pi66.cdxml")
tree=ET.parse("test.cdxml")

root=tree.getroot()


frags=root.find('page').findall('fragment')

nodes={}
bonds={}

for frag in frags:
    ns=frag.findall('n')
    bs=frag.findall('b')
    for n in ns:
        x,y=n.attrib['p'].split(' ')
        id_=n.attrib['id']
        nodes[id_]=[float(x),float(y)]
    for b in bs:
        B=b.attrib['B']
        E=b.attrib['E']
        id_=b.attrib['id']
        bonds[id_]=[B,E]


print(nodes)
print(bonds)

for id_,(x,y) in nodes.items():
    plt.scatter(x,-y,color='#000000')
    # plt.text(x,y,id_)

for id_,(a1,a2) in bonds.items():
    x0,y0=nodes[a1]
    x1,y1=nodes[a2]
    y0*=-1
    y1*=-1
    plt.plot([x0,x1],[y0,y1])
    x=(x0+x1)/2
    y=(y0+y1)/2
    plt.text(int(x),int(y),id_)
# plt.show()

# sys.exit()

b2v={
    '32':0.66,
    '33':0.66,
    '34':0.66,
    '35':0.66,
    '36':0.66,
    '37':0.66,

    '78':0.88,
    '80':0.50,
    '82':0.81,
    '84':0.50,
    '86':0.88,

    '11':0.91,
    '13':0.40,
    '15':0.84,
    '17':0.40,
    '19':0.91,

    '50':0.94,
    '51':0.21,
    '52':0.19,
    '53':0.95,
    '54':0.19,
    '55':0.21,

    '64':1.00,
    '65':0.14,
    '66':1.00,
    '67':0.14,
    '72':0.14,
    '73':1.00,
    '74':0.14,

    '246':1.08,
    '247':0.11,
    '248':1.08,
    '249':0.11,
}
v2i={}
# 进行数值与颜色之间的转化
c0=np.array([255, 0, 0])
c1=np.array([0, 0, 255])
dc=c1-c0
v0=0.11
v1=1.08
dv=v1-v0
vs=[]
cs=[]
for b,v in b2v.items():
    if v in vs:continue
    vs.append(v)
    r=(v-v0)/dv
    c=(c0+r*dc)/255
    cs.append(c)
    v2i[v]=len(cs)
    

print('颜色数量',len(cs))
colors=root.find('colortable').clear()

for i,c in enumerate(cs):
    R,G,B=c
    ele=ET.Element('color')
    ele.set('r',f"{R:.4f}")
    ele.set('g',f"{G:.4f}")
    ele.set('b',f"{B:.4f}")
    ele.set('i',f"{i+3}")
    root.find('colortable').append(ele)
    # print(f'<color r="{R:.4f}" g="{G:.4f}" b="{B:.4f}"/>')

for frag in frags:
    bs=frag.findall('b')
    for b in bs:
        B=b.attrib['B']
        E=b.attrib['E']
        id_=b.attrib['id']
        v=b2v[id_]
        i=v2i[v]
        c=cs[i-1]
        b.set('color',f'{i+1}')
        # b.set('color','2')
        # print(id_,i,i+3,c,(c*255).round(0))
# tree.write('test.cdxml')

for id_,(a1,a2) in bonds.items():
    x0,y0=nodes[a1]
    x1,y1=nodes[a2]
    y0*=-1
    y1*=-1
    v=b2v[id_]
    i=v2i[v]
    c=cs[i-1]
    r,g,b=c
    plt.plot([x0,x1],[y0,y1],color=c,linewidth=4)
    x=(x0+x1)/2
    y=(y0+y1)/2
    plt.text(int(x),int(y),id_)
plt.show()