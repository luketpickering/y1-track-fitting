#!/opt/local/bin/python2.7
import matplotlib.pyplot as plt
from sys import argv

x = [[],[]]

data_file = None
try:
    data_file = open(argv[-5], 'r')
except:
    print "Invalid Data file: " ,argv[-5]
    exit()
    
lc = 0
for line in data_file:
    lc += 1
    raw = line.split()
    x[0].append(float(raw[0]))
    x[1].append(float(raw[1]))
    
plt.hist(x[int(argv[-4])],int(argv[-3]),(int(argv[-2]),int(argv[-1])))
plt.show()