# Code for generating figure at the end of bin_edges_number.md

import numpy as np
import matplotlib.pyplot as plt

surveyArea_arr = np.logspace(-2, 1, int(1e4))
edge_size_arr = []
bin_number_arr = []

for surveyArea in surveyArea_arr:
    radius = np.sqrt(surveyArea / np.pi)  # degrees
    bin_edges_number = max(int(60 * 2 * radius) + 1, 3)
    bin_number = (bin_edges_number - 2)**2.
    edge_size = 2 * 3600 * (1.1 * radius - -1.1 * radius) / (bin_edges_number - 1)

    edge_size_arr.append(edge_size)
    bin_number_arr.append(bin_number)

fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.5))
ax[0].plot(surveyArea_arr, edge_size_arr)
ax[0].set_xlabel('surveyArea (sq. deg.)', fontsize=12)
ax[0].set_ylabel('Bigpatch Edge Size (arcseconds)', fontsize=12)
ax[0].set_xscale('log')
ax[0].axhline(132, label='132 arcseconds', c='C01')
ax[0].legend(fontsize=12)

ax[1].plot(surveyArea_arr, bin_number_arr)
ax[1].set_xlabel('surveyArea (sq. deg.)', fontsize=12)
ax[1].set_ylabel('Number of Bigpatches', fontsize=12)
ax[1].set_xscale('log')
ax[1].set_yscale('log')

fig.suptitle('Setting bin_edges_number = None', fontsize=12)
fig.tight_layout()
fig.subplots_adjust(top=.9)

fig.savefig('bin_edges_number.png')
