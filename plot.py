#!/opt/local/bin/python2.7
import matplotlib.pyplot as plt

f = open('out.txt', 'r')

sensex = [ x for x in range(0,8)]
sensey = [ y+0.5 for y in range(0,8)]

swx,swy = [], []
for _x in sensex:
    for _y in sensey:
        swx.append(_x)
        swy.append(_y)


xyz = [[],[],[]]
x, y, z = xyz[0], xyz[1], xyz[2]
for l in f:
    raw = l.split()
    print raw
    x.append(sensex[int(raw[0])])
    y.append(sensey[int(raw[1])])
    z.append(int(raw[2])*3.0)
f.close()

print xyz
plt.scatter(swx,swy, marker="*")
plt.scatter(x,y,z)
plt.show()