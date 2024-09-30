import numpy as np 
import ReadNode as rn 
import shutil

print("Calculate Displacement")
BladeName = "NREL5MW"

step_no = 1
Blade_no = 'B1-'

LeTeFile = "F:/Kuliah/Tesis/Program/pyAPDL/Output/"+BladeName+"-LeTeStatNodes.txt"
file = open(LeTeFile,'r')
lines=file.read()
lines=lines.split()

LeTeNodes = np.zeros((int(len(lines)/4),4))
for i in range ((int(len(lines)/4))):
    for j in range (4):
        LeTeNodes[i,j] = lines[i*4+j]
'''
target_name = BladeName+'-Displacement-'+str(step_no)+'.txt'
ori = r'F:/Kuliah/Tesis/Program/ANSYS_APDL/NREL5MW/Displacement-'+str(step_no)+'.txt'
target = r'F:/Kuliah/Tesis/Program/pyAPDL/Input/'+target_name
shutil.copyfile(ori, target)
'''
target_name = BladeName+'-Displacement-'+str(step_no)+'.txt'
dIn="F:/Kuliah/Tesis/Program/pyAPDL/Input/"+target_name
dOut="F:/Kuliah/Tesis/Program/pyAPDL/Output/"+BladeName+"-OutDisplacement.txt"
NodeDisplacement,LeTeDisplacement=rn.ReadDisplacement(dIn,dOut,LeTeNodes)

LeTeDeformed,DeformedAngle=rn.Deformed(LeTeNodes,LeTeDisplacement)
np.savetxt("Output/"+BladeName+"-"+Blade_no+"LeTeDisplacement-"+str(step_no)+".txt",LeTeDisplacement)
#np.savetxt("Output/"+BladeName+"-"+Blade_no+"LeTeDeformed-"+str(step_no)+".txt",LeTeDeformed)
#np.savetxt("Output/"+BladeName+"-"+Blade_no+"DeformedAngle-"+str(step_no)+".txt",DeformedAngle)


'''
bem_input='F:/Kuliah/Tesis/Program/pyBEMT-master - Copy/NREL5MW-'+Blade_no+'1.ini'
new_input='F:/Kuliah/Tesis/Program/pyBEMT-master - Copy/NREL5MW-'+Blade_no+str(step_no+1)+'.ini'
with open(bem_input,'r') as f:
    data=f.readlines()
    new_angle = ' '.join([str(item) for item in DeformedAngle])
    data[16]='pitch = '+new_angle

temp = ''.join([str(item) for item in data])
with open(new_input,'w') as f:
    f.write(temp)
'''