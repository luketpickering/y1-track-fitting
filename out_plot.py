#!/opt/local/bin/python2.7
import matplotlib.pyplot as plt
from sys import argv

#plt.xkcd()
print "starting."
color_ring = ['red','green','blue','purple']

def get_xy_from_m_c(ln):
    return [[0.0,80000.0],[ln[1],ln[0]*80000.0 + ln[1]]]

def construct_micron_grid():
    grid = []
    size = 8
    for i in range(size):
            grid.append([ ( y*10000.0 if not i%2 else (y*10000.0 + 5000.0) ) for y in range(size)])
    return grid
    
    
grid = construct_micron_grid()
data_file = None
try:
    data_file = open(argv[-1], 'r')
except:
    print "Invalid Data file"
    exit()
else:
    swx,swy = [], []
    for _i in range(len(grid)):
        for _y in grid[_i]:
            swx.append(_i*10000.0)
            swy.append(_y)
    plt.scatter(swx,swy, marker="*", color='red')
    
CA = []
eqs = []
for i, line in enumerate(data_file):
    raw = line.split()
    print i, len(raw), line
    if i < 8:
        CA.append([float(x) for x in raw])
    elif i < 12:
        eqs.append([float(x) for x in raw])
    else:
        raise IOError
print "trying to plot."
fig = plt.gcf()
for c in CA:
    print c
    c_art = plt.Circle((c[0],c[1]), c[2])
    fig.gca().add_artist(c_art)

for i, eq in enumerate(eqs):
    print eq
    xy = get_xy_from_m_c(eq)
    plt.plot(xy[0],xy[1],color=color_ring[i], linewidth=0.2)
plt.savefig('funky.pdf')