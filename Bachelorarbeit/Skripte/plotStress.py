import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def f(r):
	return [4.85827/math.sqrt(i) for i in r]
	#return [1.12*0.007/1.25e-10/1e6*math.sqrt(math.pi*0.012)/math.sqrt(2*math.pi*i) for i in r]
	

file = open('stress.txt','r')
lines = file.readlines();
x = []
y = []
for line in lines:
	x.append(float(line.split(' ')[0])+0.0135)
	y.append(float(line.split(' ')[1])/1.25e-10/1e6)
plt.plot(x,y,'x',label=r'lammps')
other = f(x)
plt.plot(x,other,'x',label=r'theory')
plt.legend()
plt.title(r'Comparison of Stress Intensity Factor and LAMMPS stress')
plt.xlabel(r'Radius [$m$]')
plt.ylabel(r'Stress [$MPa$]')
plt.savefig('stressComparison.pdf',dpi=2000)
