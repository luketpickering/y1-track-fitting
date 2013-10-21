#!/opt/local/bin/python2.7
import matplotlib.pyplot as plt
from math import sqrt, pow
from sys import argv
from random import gauss, uniform

step_p=0.5
step_g=0.05
step_c=50.0
phi_initial = 0.0

j_sig_p = 1.0
j_sig_g = 0.05
j_sig_k = 200.0

class param_blob:
    lhood_eval = 0
    samples = 0

    def __init__(self, phi, grad, interkept):
        self.phi = phi
        self.grad = grad
        self.interkept = interkept

    def sample(self, X,Y,T):
        self.lhood_eval = neg_llhood(X,Y,T,self)
        return self.lhood_eval

    def __lt__(self,other):
        return (self.lhood_eval < other.lhood_eval)
    def __le__(self,other):
        return (self.lhood_eval <= other.lhood_eval)
    def __gt__(self,other):
        return (self.lhood_eval > other.lhood_eval)
    def __le__(self,other):
        return (self.lhood_eval >= other.lhood_eval)
    def __repr__(self):
        return self.__str__()
    def __str__(self):
        return "\nlhood_eval: %s, Phi:%s, Grad:%s, K:%s" %(self.lhood_eval, self.phi, self.grad,self.interkept)

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

def find_grad(p1,p2):
    return (p2[1] - p1[1])/(p2[0] - p1[0])

def find_intercept(p1,grad):
    return grad*p1[0] - p1[1]

def init_guess(X,Y,T):
    
    t = T[:]
    t_min = min(t)
    t_min_ind = t.index(t_min)
    del t[t_min_ind]
    t_2nd_min = min(t)
    t_2nd_ind = t.index(t_2nd_min)
    
    min_inds = sorted([t_min_ind, t_2nd_ind])
    
    print min_inds
    
    grad = find_grad( (X[min_inds[0]],Y[min_inds[0]]),
                      (X[min_inds[1]],Y[min_inds[1]]))
    return ( grad,
             find_intercept((X[min_inds[0]],Y[min_inds[0]]), grad)
            )
            
def neg_llhood(X,Y,T,P):
    ss = 0
    l2p = 0
    dft = 0
    for i in range(len(X)):
        l2p = prop_track_distance(X[i],Y[i], P.grad, P.interkept)
        dft = drift_distance(P.phi, T[i])
        ss += (l2p - dft)*(l2p - dft)
    return ss

def MCMC_param_jump(P):
    nphi = gauss(P.phi,j_sig_p)
    while(nphi < 0):
        nphi = gauss(P.phi,j_sig_p)
    grad = gauss(P.grad,j_sig_g)
    interkept = gauss(P.interkept,j_sig_k)

    return param_blob(nphi,grad,interkept)

def MCMC_step(P, X,Y,T):

    prop_P = MCMC_param_jump(P)
    if P.lhood_eval == 0:
        print P, P.sample(X,Y,T)
    a = P.lhood_eval/prop_P.sample(X,T,Y)
    
    if( a >= 1.0):
        prop_P.samples += 1
        return prop_P
    else:
        uniform_pick = uniform(0.0,1.0)
        if ( uniform_pick < a ):
            prop_P.samples += 1
            return prop_P
        else:
            P.samples += 1
            return P


def construct_micron_grid():
    grid = []
    size = 8
    for i in range(size):
            grid.append([ ( y*10000.0 if not i%2 else (y*10000.0 + 5000.0) ) for y in range(size)])
    return grid

grid = construct_micron_grid()
data_file = None
try:
    data_file = open(argv[1], 'r')
except:
    print "Invalid Data file"
    exit()
xyz = [[],[],[]]
X, Y, T = xyz[0], xyz[1], xyz[2]

for line in data_file:
    raw = line.split()
    X.append(int(raw[0])*10000.0)
    Y.append(grid[int(raw[0])][int(raw[1])])
    T.append(int(raw[2])*1.0)
data_file.close()

initial_guess = init_guess(X,Y,T)

#ss_est = (neg_llhood(X,Y,T,initial_guess[0],initial_guess[1],
#                                    phi_initial),
#            initial_guess[0],
#            initial_guess[1],
#            phi_initial)

#phi, grad, kintercept = ss_est[3], ss_est[1], ss_est[2]
#print "Best initial fit: grad:%s, intercept:%s, phi:%s" % (grad, kintercept, phi)
"""
for i in range(100):
    print "I: ",i , "MIN_V:%s, MIN_g:%s, MIN_c:%s MIN_phi:%s" % ss_est
    for j in range(-50,51):
        for k in range(-50,51):
            
            phi = phi_initial + i*step_p
            grad = initial_guess[0] + j*step_g
            kintercept = initial_guess[1] + k*step_c
            
            new_est = neg_llhood(X,Y,T,grad,kintercept,phi)
            if new_est < ss_est[0]:
                ss_est = (new_est, grad, kintercept, phi)

#ss_est = (0,ig[0], ig[1] + 0.1, 1.0)
#igx, igy = [0,7], [ig[1], ig[1] + ig[0]*X[7]]
"""

steps = 0
samples = []
P = param_blob(phi_initial,initial_guess[0],initial_guess[1])
P.sample(X,Y,T)
print "Best initial fit: grad:%s, intercept:%s, phi:%s" % (P.grad, \
        P.interkept, P.phi)
while steps < 10000:
    np = MCMC_step(P,X,Y,T)
    if np.samples == 1:
        steps += 1
        P = np
        if steps == 1000:
            print "Finished burn in."
        if steps > 1000:
            samples.append(np)


#fitx, fity = [0,70000.0], \
#            [-1.0*ss_est[2], -1.0*ss_est[2] + ss_est[1]*70000.0]
#Dt = map(lambda x: float(x)*ss_est[3]/10.0, T)

min_P = min(samples)
#print samples

print "Min Found: ", min_P

fitx, fity = [0,70000.0], \
             [-1.0*min_P.interkept, \
                -1.0*min_P.interkept + min_P.grad*70000.0]

Dt = map(lambda x: float(x)*min_P.phi, T)


#print X,Y,T,ss_est

#for o in range(len(X)):
#    print "X: ", o, " l2p: ", prop_track_distance(X[o],Y[o],ss_est[1],ss_est[2]), " dft: ", drift_distance(ss_est[3],T[o])
#print "neg_llhood = ", neg_llhood(X,Y,T,ss_est[1],ss_est[2],ss_est[3])

swx,swy = [], []
for _i in range(len(grid)):
    for _y in grid[_i]:
        swx.append(_i*10000.0)
        swy.append(_y)

plt.scatter(X,Y,Dt, marker="o")
plt.scatter(swx,swy, marker="*", color='red')
plt.scatter(X,Y, color='green')
plt.plot(fitx,fity, color='red')
plt.savefig('out.pdf')


