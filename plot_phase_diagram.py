"""
This program reads the data in phase_diagram.dat
and draws the phase diagram of the alloy

Clovis Caetano - July/2015
"""
import matplotlib.pyplot as plt
import numpy as np

T = []
xb = []
xs = []
with open("phase_diagram.dat","r") as input_file:
    for line in input_file:
        data = line.split()
        T.append(float(data[0]))
        T.append(float(data[0]))
        xb.append(float(data[1]))
        xb.append(float(data[2]))
        xs.append(float(data[3]))
        xs.append(float(data[4]))

x_bin = [x for (x,y) in sorted(zip(xb,T))]
T_bin = [y for (x,y) in sorted(zip(xb,T))]
x_spin = [x for (x,y) in sorted(zip(xs,T))]
T_spin = [y for (x,y) in sorted(zip(xs,T))]

plt.plot(x_bin,T_bin,'b',lw=1)
plt.fill_between(x_bin,T_bin,color='b',alpha=0.3)
plt.plot(x_spin,T_spin,'r',lw=1)
plt.fill_between(x_spin,T_spin,color='r',alpha=0.3)

font = {'family': 'Times New Roman','size': 28}
plt.xlim([0.0,1.0])
plt.ylim([200.0,350.0])
plt.xticks(np.arange(0.0,1.1,0.1))
plt.yticks(np.arange(200.0,375.0,25.0))
plt.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('$x_{\mathregular{Br}}$',fontdict=font)
plt.ylabel('$T$ (K)',fontdict=font)
plt.savefig('phase_diagram.pdf',bbox_inches='tight')
plt.close()
