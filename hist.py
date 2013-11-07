#!/opt/local/bin/python2.7
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sys import argv
from settings import IAMSILLY

print "IAMSILLY = %s" % IAMSILLY

v_drift = []
angl = []

mu_vd = 52.82
sigma_vd = 0.07
bins_vd = 100
range_vd = (52.3,53.3)

mu_angl = 0.11
sigma_angl = 0.07
bins_angl = 56
range_angl = (-0.4,0.4)

data_file = None
try:
    data_file = open(argv[-1], 'r')
except:
    print "Invalid Data file: " ,argv[-1]
    exit()
    
lc = 0
for line in data_file:
    lc += 1
    raw = line.split()
    v_drift.append(float(raw[0]))
    angl.append(float(raw[1]))


if IAMSILLY:
    try:
        plt.xkcd()
    except:
        print "Couldn't Load the XKCD module for matplotlib."
    
fig = plt.figure()
ax = fig.add_subplot(2,1,1)

plt.tick_params(\
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
plt.tick_params(\
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        top='off')      # ticks along the bottom edge are off

a, xbins, pat = ax.hist(v_drift, bins_vd, range_vd, histtype='step')
max_c_vd = max(a)


if IAMSILLY:
    plt.annotate(
        'THE PEAK',
        xy=(mu_vd, max_c_vd), arrowprops=dict(arrowstyle='->'), xytext=(mu_vd+0.2, 45000 ))
    plt.annotate(
        "WHERE IT'S GONE DOWN\n          BY A BIT",
        xy=(mu_vd-sigma_vd,31000 ), arrowprops=dict(arrowstyle='->'), 
            xytext=(52.22, 40000 ))


y = mlab.normpdf(xbins, mu_vd, sigma_vd)*9500.0
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.set_xlabel('Drift Velocity/[km/s]')
ax.set_title(r'Drift Vel distribution: mean(%.2f), stdv(%.2f)' % (mu_vd,sigma_vd))
ax.plot(xbins, y, 'r--', color='red')

ax = fig.add_subplot(2,1,2)
a, xbins, pat = ax.hist(angl, bins_angl, range_angl, histtype='step',\
                    color='orange')
max_c_angl = max(a)
                    
if IAMSILLY:                    
    plt.annotate(
        "LOOKS PRETTY\n    NORMAL",
        xy=(mu_angl, max_c_angl), arrowprops=dict(arrowstyle='->'), xytext=(-0.3, 40000 ))

plt.tick_params(\
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
plt.tick_params(\
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        top='off')      # ticks along the bottom edge are off
y = mlab.normpdf(xbins, mu_angl, sigma_angl)*14000.0
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.set_xlabel('Incident Angle/[rads]')
ax.set_title(r'Incident angle distribution: mean(%.2f), stdv(%.2f)' % 
            (mu_angl,sigma_angl))
ax.plot(xbins, y, 'r--', color='green')

fig.subplots_adjust(hspace=0.6)

plt.savefig('silly-hist.pdf')