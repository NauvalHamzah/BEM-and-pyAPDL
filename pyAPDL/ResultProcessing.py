import numpy as np  
import shutil

BladeName = 'NREL5MW-'
Blade_no = 'B3-'
step = 13

all_df = np.zeros((step,3))

for i in range(step):
    df ='F:/Kuliah/Tesis/Program/pyBEMT-master - Copy/output/NREL5MW-'+Blade_no+str(i+1)+'-df.txt'
    with open(df,'r') as f:
        a=f.readlines()
        a[0]=a[0].split(' ')
    for j in range(3):
        a[0][1+j] = float(a[0][1+j])
        all_df[i,j] = a[0][1+j]    

df = 'F:/Kuliah/Tesis/Program/pyBEMT-master - Copy/output/NREL5MW-'+Blade_no+'All_df.txt'
np.savetxt(df, all_df)

file = 'F:/Kuliah/Tesis/Program/pyAPDL/Output/NREL5MW-'+Blade_no+'DeformedAngle-1.txt'
with open(file,'r') as f:
    a=f.readlines()

all_DeformedAngle = np.zeros((len(a),step))

for i in range (step):
    DeformedAngle = 'F:/Kuliah/Tesis/Program/pyAPDL/Output/NREL5MW-'+Blade_no+'DeformedAngle-'+str(i+1)+'.txt'
    with open(DeformedAngle,'r') as f:
        a=f.readlines()
        for j in range (len(a)):            
            all_DeformedAngle[j,i] = a[j]

DeformedAngle = 'F:/Kuliah/Tesis/Program/pyAPDL/Output/NREL5MW-'+Blade_no+'All_DeformedAngle.txt'
np.savetxt(DeformedAngle, all_DeformedAngle)




file = 'F:/Kuliah/Tesis/Program/pyAPDL/Output/NREL5MW-'+Blade_no+'LeTeDisplacement-1.txt'
with open(file,'r') as f:
    a=f.readlines()

Edgewise_LeTeDisplacement = np.zeros((len(a),step))
Flapwise_LeTeDisplacement = np.zeros((len(a),step))

for i in range (step):
    LeTeDisplacement = 'F:/Kuliah/Tesis/Program/pyAPDL/Output/NREL5MW-'+Blade_no+'LeTeDisplacement-'+str(i+1)+'.txt'
    with open(LeTeDisplacement,'r') as f:
        a=f.readlines()
        for j in range (len(a)):  
            a[j] = a[j].split(' ')          
            Edgewise_LeTeDisplacement[j,i] = a[j][1]
            Flapwise_LeTeDisplacement[j,i] = a[j][2]

Edge_file = 'F:/Kuliah/Tesis/Program/pyAPDL/Output/NREL5MW-'+Blade_no+'Edgewise_LeTeDisplacement.txt'
Flap_file = 'F:/Kuliah/Tesis/Program/pyAPDL/Output/NREL5MW-'+Blade_no+'Flapwise_LeTeDisplacement.txt'
np.savetxt(Edge_file, Edgewise_LeTeDisplacement)
np.savetxt(Flap_file, Flapwise_LeTeDisplacement)
