import os,glob,sys

sdir = '/home/adrian/Desktop/tornado1/meso'
os.chdir(sdir)

slist = glob.glob('*.nc')

for s in slist:

    print s

    os.system('ncl_convert2nc '+s)
