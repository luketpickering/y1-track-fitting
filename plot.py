#!/opt/local/bin/python2.7
import matplotlib.pyplot as plt
from math import sqrt, pow, asin, cos, sin, pi
from sys import argv

def drift_distance(phi, tdc):
    rtn = (phi*tdc)
    if rtn < 0:
        raise ValueError
    return rtn

def prop_track_distance(x0, y0, grad, kintercept):
    rtn = abs(y0 - grad*x0 + kintercept)/sqrt(1 + grad*grad)
    if rtn < 0:
        raise ValueError
    return rtn

#C = (R,X,Y)
def get_sinang_com_tan_ext(C1,C2):
    dx = C2[1] - C1[1]
    dy = C2[2] - C1[2]
    sinal = (C2[0] - C1[0])/sqrt(dx*dx + dy*dy)
    return sinal

def get_sinang_com_tan_int(C1,C2):
    dx = C2[1] - C1[1]
    dy = C2[2] - C1[2]
    sinal = (C2[0] + C1[0])/sqrt(dx*dx + dy*dy)
    if sinal < -1.0:
        raise ValueError
    return sinal

def rot_vector(V,theta):
    return ((V[0]*cos(theta) - V[1]*sin(theta)), V[0]*sin(theta) + V[1]*cos(theta))

def get_unit_vect_2P(P1,P2):
    dx = P2[0] - P1[0]
    dy = P2[1] - P1[1]
    len = sqrt(dx*dx + dy*dy)
    return (dx/len,dy/len)

def get_unit_vect_2C(C1,C2):
    dx = C2[1] - C1[1]
    dy = C2[2] - C1[2]
    len = sqrt(dx*dx + dy*dy)
    return (dx/len,dy/len)

def cotangents(C1, C2):
    #ext
    def cmp(x):
        return x[0]
    circs = sorted([C1,C2], key=cmp)
    
    dx = C2[1] - C1[1]
    dy = C2[2] - C1[2]
    len = sqrt(dx*dx + dy*dy)
    
    if len < (C1[0] + C2[0]):
        raise ValueError
    
    #print circs
    
    sin_ext = get_sinang_com_tan_ext(circs[0],circs[1])
    sin_int = get_sinang_com_tan_int(circs[0],circs[1])
    
    #print sin_ext, sin_int
    
    th_ext = asin(sin_ext)
    th_int = asin(sin_int)
    c2c_uv = get_unit_vect_2C(circs[0],circs[1])

    uvs = (rot_vector(c2c_uv,th_ext), rot_vector(c2c_uv, -1.0*th_ext),
            rot_vector(c2c_uv,th_int), rot_vector(c2c_uv, -1.0*th_int))
    return uvs

def vadd(v1,v2):
    return (v1[0] + v2[0], v1[1] + v2[1])
def vscale(v1,sf):
    return (v1[0] * sf, v1[1] * sf)

def get_pvects_cotangents(C1,C2):
    
    def cmp(x):
        return x[0]
    circs = sorted([C1,C2], key=cmp)
    
    sin_ext = get_sinang_com_tan_ext(circs[0],circs[1])
    sin_int = get_sinang_com_tan_int(circs[0],circs[1])

    #print sin_ext, sin_int
    
    th_ext = asin(sin_ext)
    th_int = asin(sin_int)
    c2c_uv = get_unit_vect_2C(circs[0],circs[1])

    uvs = (rot_vector(c2c_uv,th_ext), rot_vector(c2c_uv, -1.0*th_ext),
            rot_vector(c2c_uv,th_int), rot_vector(c2c_uv, -1.0*th_int))
    
    nty_rad = pi*90.0/180.0

    ah = rot_vector(c2c_uv, th_ext + nty_rad)
    A = vadd(vscale(ah,circs[0][0]), (circs[0][1],circs[0][2]))
    bh = rot_vector(c2c_uv, -1.0*(th_ext + nty_rad))
    B = vadd(vscale(bh,circs[0][0]), (circs[0][1],circs[0][2]))
    ch = rot_vector(c2c_uv, th_int - nty_rad )
    C = vadd(vscale(ch,circs[0][0]), (circs[0][1],circs[0][2]))
    dh = rot_vector(c2c_uv, nty_rad - th_int )
    D = vadd(vscale(dh,circs[0][0]), (circs[0][1],circs[0][2]))

    return (A,B,C,D)

def get_x_y_lists_P(P1,P2):
    return [ [P1[0],P2[0]], [P1[1],P2[1]] ]

def push_p_along_v_by_dx(P1, uv, dx):
    #dy = uv[1] * (dx/uv[0])
    #return vadd(P1, (dx,dy))
    numdx = dx/uv[0]
    return vadd(P1, vscale(uv,numdx))

def get_line_eq_uv_pt(uv,pt):
    m = uv[1]/uv[0]
    c = pt[1] - m* pt[0]
    return (m, c)
def get_xy_from_m_c(ln):
    return [[0.0,80000.0],[ln[1],ln[0]*80000.0 + ln[1]]]

def circ_c(C1):
    return (C1[1],C1[2])

def construct_micron_grid():
    grid = []
    size = 8
    for i in range(size):
            grid.append([ ( y*10000.0 if not i%2 else (y*10000.0 + 5000.0) ) for y in range(size)])
    return grid

def line_miss_all(circs, ln):
    for i,c in enumerate(circs):
        if i == 0 or i == 7:
            continue
        leng = prop_track_distance(c[1],c[2],ln[0], -1.0*ln[1])
        if leng < c[0]:
            return False
    return True

grid = construct_micron_grid()
data_file = None
try:
    data_file = open(argv[1], 'r')
except:
    print "Invalid Data file"
    exit()
else:
    swx,swy = [], []
    for _i in range(len(grid)):
        for _y in grid[_i]:
            swx.append(_i*10000.0)
            swy.append(_y)
    #plt.scatter(swx,swy, marker="*", color='red')

i_phi = 38.0

CA = []
for line in data_file:
    raw = line.split()
    CA.append([
    int(raw[2])*i_phi,
    int(raw[0])*10000.0,
    grid[int(raw[0])][int(raw[1])]
    ])
data_file.close()

def cmp(x):
        return x[0]
circs = sorted([CA[0],CA[7]], key=cmp)
print circs
color_ring = ['red','green','blue','purple']

n_phi = i_phi
one_line = 0
while (one_line != 1) and ( n_phi > 0.001):
    o_phi = n_phi
    n_phi -= 0.1
    for c in CA:
        c[0] = c[0]*n_phi/o_phi
    ct = cotangents(circs[0],circs[1])

    dx = circs[0][1] - circs[1][1]

    ABCD = get_pvects_cotangents(circs[0],circs[1])

    tan_eqs = []
    poss_lines = 0

    for i, tan_pt in enumerate(ABCD):
        eq = get_line_eq_uv_pt(ct[i],tan_pt)
        #wl = get_xy_from_m_c(eq)
        #plt.plot( wl[0],wl[1], color=color_ring[i])
        if line_miss_all(CA, eq):
            poss_lines ^= ( 1 << i )

    one_line = 0
    for i in range(len(ABCD)):
        #print i, ( 1 << i )
        if poss_lines & ( 1 << i ):
            one_line += 1
    #print "Possible lines = ", poss_lines & 1, poss_lines & 2,\
           #     poss_lines & 4, poss_lines & 8

print n_phi
for i, pt in enumerate(ABCD):
    if poss_lines & ( 1 << i ):
        ln_end = push_p_along_v_by_dx(pt,ct[i],-1.0*dx)
        plt.scatter([pt[0], ln_end[0]], [pt[1], ln_end[0]], color=color_ring[i])
        whole_ln = get_x_y_lists_P(ln_end, pt)
        plt.plot(whole_ln[0], whole_ln[1], color=color_ring[i])
print ABCD
l5 = get_x_y_lists_P(circ_c(circs[0]),circ_c(circs[1]))

fig = plt.gcf()
for c in CA:
    c_art = plt.Circle((c[1],c[2]), c[0])
    fig.gca().add_artist(c_art)

plt.plot(l5[0],l5[1], color='orange')
plt.show()


