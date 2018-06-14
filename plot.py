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


def linear_regression(x, y, yerr):

	ones = np.ones(len(x))
	Y = np.array(y)
	A = np.array([ones, x])
	C = np.diag(np.array(yerr)**2)
	Cinv = np.diag(1/np.array(yerr)**2)

	v1 = np.linalg.inv(np.dot(np.dot(A, Cinv), A.T))
	v2 = np.dot(np.dot(A, Cinv), Y)

	X = np.dot(v1, v2)

	return X


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Plot a file.')
	parser.add_argument("plot", action="store", type=str)
	# parser.add_argument("type", action="store", type=str)
	args = parser.parse_args()

	if 'dist' in args.plot:
		file = 'output/dist.dat'

		vals = []
		with open(file) as f:
			for line in f:
				vals.append(float(line))

		xarray = np.arange(min(vals), max(vals), .1)
		yarray = np.exp(-xarray**2/20) / math.sqrt(math.pi)

		plt.hist(vals, facecolor='w', density=True, bins=50, histtype='stepfilled')
		plt.plot(xarray, yarray, "--", color='r')
		# plt.legend(loc='upper right')
		plt.xlim(min(vals), max(vals))
		plt.savefig('output/dist.png')
		plt.show()
		plt.close()

	if 'walkers' in args.plot:
		file = 'output/dist.dat'

		vals = []
		with open(file) as f:
			for line in f:
				vals.append(float(line))
		steps = np.arange(0,len(vals))
		avg = round(np.mean(vals),3)
		std = round(np.std(vals)/math.sqrt(len(vals)),3)
		avg_array = np.array([avg for i in vals])

		plt.figure(figsize=[16,5])
		plt.step(steps, vals, color='k', alpha=.7, linewidth=1, label=r'$T=100$')
		plt.plot(steps, avg_array, color='r', label=r'$<E>=%s \pm %s$'%(str(avg), str(std)))
		plt.fill_between(steps, avg_array-std, avg_array+std, color='r', alpha=.1)
		plt.legend(loc='upper right')
		plt.xlim(steps[0], steps[-1])
		plt.xlabel('step', fontsize=15)
		plt.ylabel('Energy', fontsize=15)
		plt.title(' walker', fontsize=20)
		plt.savefig('output/walkers.png')
		plt.show()
		plt.close()

	if 'combined' in args.plot:
		file = 'output/dist.dat'

		vals = []
		with open(file) as f:
			for line in f:
				vals.append(float(line))
		steps = np.arange(0,len(vals))
		avg = round(np.mean(vals),3)
		std = round(np.std(vals)/math.sqrt(len(vals)),3)
		avg_array = np.array([avg for i in vals])

		# fig, (ax1, ax2) = plt.subplots(figsize=[18,6], sharey=True)
		plt.figure(figsize=[20,5])
		grid = plt.GridSpec(1, 10, wspace=0.2)

		ax1 = plt.subplot(grid[0, :8])
		ax2 = plt.subplot(grid[0, 8:], sharey=ax1)
		plt.setp(ax2.get_yticklabels(), visible=False)

		ax1.step(steps, vals, color='k', alpha=.7, linewidth=1)
		ax1.plot(steps, avg_array, color='r', label=r'$<x>=%s \pm %s$'%(str(avg), str(std)))
		ax1.fill_between(steps, avg_array-std, avg_array+std, color='r', alpha=.1)
		ax2.hist(vals, facecolor='w', density=True, bins=50, histtype='stepfilled', orientation="horizontal")
		ax2.axhline(avg, color='r')
		ax2.set_xlabel(r'$|\psi_{0}^2|$', fontsize=20)
		ax1.set_ylabel(r'$x$', fontsize=20)
		ax1.set_xlabel('lattice index', fontsize=15)
		ax1.set_xlim(min(steps), max(steps))
		ax1.set_ylim(-4,4)
		ax2.set_xlim(0,.6)
		ax1.legend(loc='upper left')
		plt.suptitle('MC Path Integral Demo (N=100)', fontsize=20)
		plt.savefig('output/combined.png')
		plt.show()
		plt.close()


	if 'animation' in args.plot:

		folder = os.listdir('output/xvec/')
		print('Plotting %s files...'%(len(folder)))
		for file in folder:
			vals = []
			with open('output/xvec/' + file) as f:
				for line in f:
					vals.append(float(line))
			i = file.split('x')[1].split('.csv')[0]

			steps = np.arange(0,len(vals))
			avg = round(np.mean(vals),3)
			std = round(np.std(vals)/math.sqrt(len(vals)),3)
			avg_array = np.array([avg for i in vals])

			plt.figure(figsize=[20,5])
			grid = plt.GridSpec(1, 10, wspace=0.2)

			ax1 = plt.subplot(grid[0, :8])
			ax2 = plt.subplot(grid[0, 8:], sharey=ax1)
			plt.setp(ax2.get_yticklabels(), visible=False)

			ax1.step(steps, vals, color='k', alpha=.7, linewidth=1, label='iteration='+str(i))
			ax1.plot(steps, avg_array, color='r', label=r'$<x>=%s \pm %s$'%(str(avg), str(std)))
			ax1.fill_between(steps, avg_array-std, avg_array+std, color='r', alpha=.1)
			ax1.set_ylabel(r'$x$', fontsize=20)
			ax1.set_xlabel('lattice index', fontsize=15)
			ax1.set_xlim(min(steps), max(steps))
			ax1.set_ylim(-4,4)

			ax2.hist(vals, facecolor='w', density=True, bins=50, histtype='stepfilled', orientation="horizontal")
			ax2.axhline(avg, color='r')
			ax2.set_xlabel(r'$|\psi_{0}^2|$', fontsize=20)
			ax2.set_xlim(0,.6)

			ax1.legend(loc='upper left')
			plt.suptitle('MC Path Integral Demo (N=100)', fontsize=20)
			plt.savefig('output/animation/dist'+str(i)+'.png')
			# plt.show()
			plt.close()

	if 'gif' in args.plot:
		l = len(os.listdir("output/animation/")) - 1
		print("Converting {} images to gif...".format(l+1))

		os.system("convert -delay 0 $(for i in $(seq 0 1 %s); do echo output/animation/dist${i}.png; done) \
			-loop 0  output/MCMC_animation.gif"%(str(l)))

	if 'expected' in args.plot:
		file1 = 'output/expected_energy.dat'
		file2 = 'output/expected_error.dat'

		vals, errs = [], []
		with open(file1) as f:
			for line in f:
				vals.append(float(line))
		with open(file2) as f:
			for line in f:
				errs.append(float(line))

		vals, errs = np.array(vals), np.array(errs)
		steps = np.arange(0,len(vals))
		avg = round(np.mean(vals),3)

		lin_reg = linear_regression(steps[1:], vals[1:], errs[1:])
		reg_y = lin_reg[1]*steps[1:] + lin_reg[0]

		plt.figure(figsize=[12,6])
		plt.plot(steps, vals, color='k', alpha=.7, label=r'$\langle E \rangle \pm error$')
		plt.fill_between(steps, vals-errs, vals+errs, color='k', alpha=.3)
		if lin_reg[0] > 0:
			plt.plot(steps[1:], reg_y, '--', color='r', label=r'$y = %s x + %s$'%(round(lin_reg[1],3), round(lin_reg[0],3)))
		else:
			plt.plot(steps[1:], reg_y, '--', color='r', label=r'$y = %s x %s$'%(round(lin_reg[1],3), round(lin_reg[0],3)))
		plt.legend(loc='upper left')
		plt.xlim(steps[0], steps[-1])
		plt.ylim(min(vals), max(vals))
		plt.xlabel('Tempertaure', fontsize=15)
		plt.ylabel(r'$<E>$', fontsize=15)
		plt.title(' expected energy', fontsize=20)
		plt.savefig('output/expected.png')
		plt.show()
		plt.close()


	if 'exp_sq' in args.plot:
		file1 = 'output/expected_energy.dat'

		vals = []
		with open(file1) as f:
			for line in f:
				vals.append(float(line))

		vals = np.array(vals)**2
		steps = np.arange(0,len(vals))
		avg = round(np.mean(vals),3)

		plt.figure(figsize=[12,6])
		plt.plot(steps, vals, color='k', alpha=.7)
		
		# plt.legend(loc='upper left')
		plt.xlim(steps[0], steps[-1])
		plt.ylim(min(vals), max(vals))
		plt.xlabel('Tempertaure', fontsize=15)
		plt.ylabel(r'$<E^2>$', fontsize=15)
		# plt.title(' expected energy ', fontsize=20)
		plt.savefig('output/expected_sq.png')
		plt.show()
		plt.close()