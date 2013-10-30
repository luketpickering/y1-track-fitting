#!/opt/local/bin/python2.7
import matplotlib.pyplot as plt
from sys import argv

plt.xkcd()
x = []

data_file = None
try:
    data_file = open(argv[-1], 'r')
except:
    print "Invalid Data file: " ,argv[-1]
    exit()
    
lc = 0
for line in data_file:
    lc += 1
    x.append(float(line))
    
plt.hist(x,30, (30.0,100))
plt.show()