from ovito.io import import_file, export_file
from ovito.modifiers import *

# Define the custom modifier function:
def modify(frame, input, output):
    global i
    global N
    print(input['lambda'].array[i],i,frame)
    if input['lambda'].array[i]==0:
        N = N + [frame]
        i = i + 1

N = []
i = 0


node = import_file('dump_fatigue.peri*')
node.modifiers.append(SelectExpressionModifier(expression = 'abs(Position.Y)>0.001 || Position.Z>0.016 || Position.Z<0.011 || Position.X<-0.135'))
node.modifiers.append(DeleteSelectedParticlesModifier())
node.modifiers.append(PythonScriptModifier(function = modify))

for frame in range(node.source.num_frames):
    node.compute(frame)
    tmp = '[0'
    for aaa in N:
        tmp = tmp+' '+str(aaa)
    tmp = tmp + ']'
    print(tmp)