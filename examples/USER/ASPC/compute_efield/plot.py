#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

ea2db = 1./0.20819434

for damp in [0.76, 0.77, 0.78]:
    dat = np.loadtxt('dip.%4.2f.dat' %damp)
    dip = dat[:,4]*ea2db
    print 'average induced dipole moment for damping ', damp, ': ', np.average(dip), ' variance: ', np.var(dip), ' standard deviation: ', np.std(dip)
    plt.hist(dip, bins=50, label=str(damp), histtype='step', normed=True)

plt.legend()
plt.show()
