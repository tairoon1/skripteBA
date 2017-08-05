from ovito.io import import_file, export_file
from ovito.modifiers import *
import numpy as np

# Define the custom modifier function:
def modify(frame, input, output):
    global N
    print(np.amax(input['c_3[2]'].array))
    maxIndex = np.argmax(input['c_3[2]'].array)
    print(input['lambda'].array[maxIndex],maxIndex,frame)
    if input['lambda'].array[maxIndex]==0:
        N = N + [frame]

N = []


node = import_file('dump_fatigue.peri*')
node.modifiers.append(SelectExpressionModifier(expression = 'abs(Position.Y)>0.001 || Position.Z>0.011 || Position.Z<0.009 || Position.X<-0.2399'))
node.modifiers.append(DeleteSelectedParticlesModifier())
node.modifiers.append(PythonScriptModifier(function = modify))

for frame in range(node.source.num_frames):
    node.compute(frame)
    tmp = '[1'
    for aaa in N:
        tmp = tmp+' '+str(aaa)
    tmp = tmp + ']'
    print(tmp)
