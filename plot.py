#!/opt/local/bin/python2.7
import matplotlib.pyplot as plt
from math import sqrt, pow, asin, cos, sin, pi
from sys import argv

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
    
    print circs
    
    sin_ext = get_sinang_com_tan_ext(circs[0],circs[1])
    sin_int = get_sinang_com_tan_int(circs[0],circs[1])
    
    print sin_ext, sin_int
    
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
    ch = rot_vector(c2c_uv, nty_rad - th_int )
    C = vadd(vscale(ch,circs[0][0]), (circs[0][1],circs[0][2]))
    dh = rot_vector(c2c_uv, th_int - nty_rad )
    D = vadd(vscale(dh,circs[0][0]), (circs[0][1],circs[0][2]))

    return (A,B,C,D)

def get_x_y_lists_P(P1,P2):
    return [ [P1[0],P2[0]], [P1[1],P2[1]] ]

def push_p_along_v_by_dx(P1, uv, dx):
    dy = uv[1] * (dx/uv[0])
    return vadd(P1, (dx,dy))
def circ_c(C1):
    return (C1[1],C1[2])

phi = 26.0*0.0001
C1 = (142.0*phi,2.0, 3.0)
C2 = (111.0*phi,7.0,4.5)

ct = cotangents(C1,C2)

def cmp(x):
        return x[0]
circs = sorted([C1,C2], key=cmp)

dx = circs[0][1] - 1.0

ABCD = get_pvects_cotangents(C1,C2)

A_C2 = push_p_along_v_by_dx(ABCD[0],ct[0],-1.0*dx)
B_C2 = push_p_along_v_by_dx(ABCD[1],ct[1],-1.0*dx)
C_C2 = push_p_along_v_by_dx(ABCD[3],ct[2],-1.0*dx)
D_C2 = push_p_along_v_by_dx(ABCD[2],ct[3],-1.0*dx)


plt.scatter(ABCD[0][0],ABCD[0][1], color='red')
plt.scatter(ABCD[1][0],ABCD[1][1], color='green')
plt.scatter(ABCD[2][0],ABCD[2][1], color='blue')
plt.scatter(ABCD[3][0],ABCD[3][1], color='purple')

l1 = get_x_y_lists_P(A_C2, ABCD[0])
l2 = get_x_y_lists_P(B_C2, ABCD[1])
l3 = get_x_y_lists_P(C_C2, ABCD[3])
l4 = get_x_y_lists_P(D_C2, ABCD[2])
l5 =  get_x_y_lists_P(circ_c(circs[0]),circ_c(circs[1]))
circle1 = plt.Circle((circs[0][1],circs[0][2]), circs[0][0])
circle2 = plt.Circle((circs[1][1],circs[1][2]), circs[1][0])


fig = plt.gcf()
fig.gca().add_artist(circle1)
fig.gca().add_artist(circle2)
plt.scatter([0,8],[0,8])
plt.plot(l1[0],l1[1], color='red')
plt.plot(l2[0],l2[1], color='green')
plt.plot(l3[0],l3[1], color='blue')
plt.plot(l4[0],l4[1], color='purple')
plt.plot(l5[0],l5[1], color='orange')
plt.show()


