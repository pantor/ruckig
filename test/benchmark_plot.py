from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ticks = ('1 DoF', '3 DoF', '6 DoF', '7 DoF')
ruckig_average_mean = [1.034, 5.768, 13.29, 15.67]
ruckig_average_err = [0.113595, 0.2788, 0.386, 0.351]
ruckig_worst_mean = [40.8002, 73.2714, 106.79, 116.333]
ruckig_worst_err = [22.5748, 31.2448, 34.041, 35.1166]
reflexxes_average_mean = [5.421, 16.11, 32.80, 38.919]
reflexxes_average_err = [0.22, 0.465, 0.7919, 0.531]
reflexxes_worst_mean = [69.61, 127.043, 186.929, 218.70]
reflexxes_worst_err = [24.54, 51.463, 64.961, 76.97]

y_pos = np.arange(len(ticks))
height = 0.35

plt.figure(figsize=(9.0, 4.0), dpi=120)


# Average
plt.subplot(1, 2, 1)
plt.suptitle('Single-thread Benchmark on Intel i7-8700K CPU @ 3.70GHz (lower is better)', fontsize=10)

ax = plt.gca()

ax.grid(True, linestyle='--', color='grey', alpha=.25)
ax.barh(y_pos - height/2, ruckig_average_mean, height=height, xerr=ruckig_average_err, label='Ruckig')
ax.barh(y_pos + height/2, reflexxes_average_mean, height=height, xerr=reflexxes_average_err, label='Reflexxes Type IV')

ax.set_yticks(y_pos)
ax.set_yticklabels(ticks)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Mean Calculation Duration [µs]')
ax.legend()

# Worst
plt.subplot(1, 2, 2)
ax = plt.gca()

ax.grid(True, linestyle='--', color='grey', alpha=.25)
ax.barh(y_pos - height/2, ruckig_worst_mean, height=height, xerr=ruckig_worst_err, label='Ruckig')
ax.barh(y_pos + height/2, reflexxes_worst_mean, height=height, xerr=reflexxes_worst_err, label='Reflexxes Type IV')

ax.set_yticks(y_pos)
ax.set_yticklabels(ticks)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Worst Calculation Duration [µs]')
ax.legend()


plt.savefig(Path(__file__).parent.parent / 'build' / 'benchmark.png')
