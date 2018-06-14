import os
import shutil

dirs = ['output/xvec', 'output/animation']

for path in dirs:
	print('Removing files from '+ path)
	shutil.rmtree(path)
	os.mkdir(path)