import os
import glob
import sys

files=glob.glob('dump_fatigue.peri*')
for i in range(len(files)):
    files[i]=files[i].split('dump_fatigue.peri')[1]

deleteFiles = [int(i) for i in files if int(i)>int(sys.argv[1])]
for i in deleteFiles:
    os.remove('dump_fatigue.peri'+str(i))

