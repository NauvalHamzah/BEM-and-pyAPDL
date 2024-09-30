import numpy as np 
import ReadNode as rn 
import shutil
from ansys.mapdl.core import launch_mapdl

BladeName = "NREL5MW"

fileIn="F:/Kuliah/Tesis/Program/pyAPDL/Input/"+BladeName+"-InNode.txt"
fileOut="F:/Kuliah/Tesis/Program/pyAPDL/Output/"+BladeName+"-OutNode.txt"
NodeFloat,MaxNode=rn.NodeData(fileIn,fileOut)

LoadList="F:/Kuliah/Tesis/Program/pyAPDL/Input/"+BladeName+"-LoadList.txt"
StatList="F:/Kuliah/Tesis/Program/pyAPDL/Input/"+BladeName+"-StatList.txt"
r_load,twist_load,c_load,ac_load = rn.StatData(LoadList)
r_stat,twist_stat,c_stat,ac_stat = rn.StatData(StatList)

NodesPerElem="F:/Kuliah/Tesis/Program/pyAPDL/Output/"+BladeName+"-NodesPerElement.txt"
NodeFloatZ_load=rn.StatSort(r_load,NodesPerElem,NodeFloat)
NodeFloatZ_stat=rn.StatSort(r_stat,NodesPerElem,NodeFloat)

fLoadNodes="F:/Kuliah/Tesis/Program/pyAPDL/Output/"+BladeName+"-NodesAtLoad.txt"
fStatNodes="F:/Kuliah/Tesis/Program/pyAPDL/Output/"+BladeName+"-NodesAtStat.txt"
LoadLen,LoadNodes,LeTeLoadNodes,LoadChord,LoadAngle=rn.NodeAtStat(r_load,twist_load,fLoadNodes,NodeFloatZ_load)
StatLen,StatNodes,LeTeStatNodes,StatChord,StatAngle=rn.NodeAtStat(r_stat,twist_stat,fStatNodes,NodeFloatZ_stat)
np.savetxt("F:/Kuliah/Tesis/Program/pyAPDL/Output/"+BladeName+"-LeTeLoadNodes.txt",LeTeLoadNodes)
np.savetxt("F:/Kuliah/Tesis/Program/pyAPDL/Output/"+BladeName+"-LeTeStatNodes.txt",LeTeStatNodes)


step= 1
azimuth_inc = np.radians(-120)
azimuth = np.radians(0)
rpm = 11.443998288988665
blade_no = 'B1'
tilt = np.radians(-5)
time=np.zeros(step)
g_tilt = np.zeros((step,3))
g_azimuth = np.zeros((step,3))

#f = open('output/TransientCommand.txt', "w+")
#f.write('nsubst,1,1,1\noutres,erase\noutres,all,all\nkbc,0\n\n')
#f.close()

for i in range(step):
    LoadInput="F:/Kuliah/Tesis/Program/pyAPDL/Input/"+BladeName+"-Load-"+blade_no+'-'+str(i+1)+".txt"
    LoadOut="F:/Kuliah/Tesis/Program/pyAPDL/Output/"+BladeName+"-Load-"+blade_no+'-'+str(i+1)+".txt"
    rn.LoadPoint(LeTeLoadNodes,LoadInput,LoadOut,twist_load,c_load,ac_load,tilt,azimuth)    
    azimuth = rn.LoadCommand(LoadOut,step,rpm,tilt,azimuth,azimuth_inc,time,i,g_tilt,g_azimuth)

print('pyAPDL Azimuth = '+str(round(np.degrees(azimuth)-np.degrees(azimuth_inc))))
#f = open('output/TransientCommand.txt', "a") 
#f.write('lssolve,1,'+str(step)+',1')

ori_inp = r'F:/Kuliah/Tesis/Program/pyAPDL/Output/TransientCommand.txt'
target_inp = r'F:/Kuliah/Tesis/Program/ANSYS_APDL/NREL5MW/TransientCommand.txt'
ori_out = r'F:/Kuliah/Tesis/Program/pyAPDL/Output/OutCommand.txt'
target_out = r'F:/Kuliah/Tesis/Program/ANSYS_APDL/NREL5MW/OutCommand.txt'
shutil.copyfile(ori_inp, target_inp)
shutil.copyfile(ori_out, target_out)

print('Opening MAPDL')
'''
# Run MAPDL
mapdl = launch_mapdl(override=True, loglevel='ERROR')
#dbfile = 'F:/Kuliah/Tesis/Program/ANSYS_APDL/NREL5MW/master'
dbfile = 'F:/Kuliah/Tesis/Program/ANSYS_APDL/NREL5MW/Blade'
input_file = '/input,'+target_inp
output_file = '/input,'+target_out
rn.RunMAPDL(dbfile, input_file, output_file)

# Exit MAPDL
mapdl.exit()
'''
print("Calculation Complete")




