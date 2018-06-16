import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('classic')
rc('font', family='serif')
rc('figure', facecolor='w')
import os
import math
from math import sqrt, pi, exp
import argparse

#####################################
# Usage:  store plots in the local directory in folder called 'output/'
#		  in terminal run:
# 			python plot.py dist
# 			python plot.py gif
#####################################

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Plot a file.')
	parser.add_argument("plot", action="store", type=str)
	# parser.add_argument("type", action="store", type=str)
	args = parser.parse_args()

	if 'dist' in args.plot:
		folder = 'output/'
		flist = sorted(os.listdir(folder))
		nfiles = len(flist)
		print('Plotting %s files...'%(nfiles))

		for i in range(nfiles):
			file = folder + flist[i]
			name = file.split('abspsi')[1].split('0000.txt')[0]
			vals = []
			with open(file) as f:
				for line in f:
					vals.append(float(line))

			xarray = np.arange(min(vals), max(vals), .1)
			yarray = np.exp(-xarray**2/20) / math.sqrt(math.pi)

			plt.hist(vals, facecolor='w', density=True, bins=50, histtype='stepfilled', label=r'$T=%s$'+name)
			# plt.plot(xarray, yarray, "--", color='r')
			plt.xlabel(r'$x$', fontsize=20)
			plt.ylabel(r'$|\psi|^{2}$', fontsize=20)
			plt.title('Distribution vs. Temperature', fontsize=15)
			plt.legend(loc='upper right')
			plt.xlim(-3.5,3.5)
			plt.ylim(0,.65)
			plt.savefig('psi_plots/dist'+str(i)+'.png')
			# plt.show()
			plt.close()


	if 'gif' in args.plot:
		l = len(os.listdir("psi_plots/")) - 1
		print("Converting {} images to gif...".format(l+1))

		os.system("convert -delay 0 $(for i in $(seq 0 1 %s); do echo psi_plots/dist${i}.png; done) \
			-loop 0  Psi2_vs_Temperature.gif"%(str(l)))
