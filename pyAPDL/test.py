import numpy as np

'''StatList="Input/Stat_List_NREL5MW.txt"
f=open(StatList, "r")
stat=f.read()
stat=stat.split()
nstat=int(len(stat)/2)
stat_r=np.zeros(nstat)
stat_twist=np.zeros(nstat)


for i in range(nstat):
    stat_r[i]=float(stat[2*i])
    stat_twist[i]=float(stat[2*i+1])

print(stat_r)
a=5.0
a=int(a)
print(a)
'''

'''
#Coordinate Transformation
a=[1.5,2.86,5.6,8.33,11.75,15.85,19.95,24.05,28.15,32.25,36.35,40.45,44.55,48.65,52.75,56.17,58.9,61.63,63]
precone=np.radians(2.5)
azimuth=np.radians(-45)
tilt=np.radians(-5)
coor_init = np.zeros((19,3))
coor_precone = np.zeros((19,3))
coor_azimuth = np.zeros((19,3))
coor_tilt = np.zeros((19,3))

coor_init[:,2] = a

for i in range (19):
    coor_precone[i,0] = coor_init[i,0]
    coor_precone[i,1] = np.cos(precone)*(coor_init[i,1]-coor_init[0,1])-np.sin(precone)*(coor_init[i,2]-coor_init[0,2])+coor_init[0,1]
    coor_precone[i,2] = np.sin(precone)*(coor_init[i,1]-coor_init[0,1])+np.cos(precone)*(coor_init[i,2]-coor_init[0,2])+coor_init[0,2]

    coor_azimuth[i,0] = np.cos(azimuth)*(coor_precone[i,0])-np.sin(azimuth)*(coor_precone[i,2])
    coor_azimuth[i,1] = coor_precone[i,1]
    coor_azimuth[i,2] = np.sin(azimuth)*(coor_precone[i,0])+np.cos(azimuth)*(coor_precone[i,2])

    coor_tilt[i,0] = coor_azimuth[i,0]
    coor_tilt[i,1] = np.cos(tilt)*coor_azimuth[i,1]-np.sin(tilt)*coor_azimuth[i,2]
    coor_tilt[i,2] = np.sin(tilt)*coor_azimuth[i,1]+np.cos(tilt)*coor_azimuth[i,2]


#Gravity Acceleration Transformation

g_init = [0,5735.604773,0]
g_tilt = np.zeros(3)
g_azimuth = np.zeros(3)

g_tilt[0] = 0
g_tilt[1] = np.cos(-tilt)*g_init[1]
g_tilt[2] = np.sin(-tilt)*g_init[1]

g_azimuth[0] = -np.sin(-azimuth)*g_tilt[2]
g_azimuth[1] = g_tilt[1]
g_azimuth[2] = np.cos(-azimuth)*g_tilt[2]

print(g_azimuth)
#Wind Shear
h0 = 90
v0 = 10
shear = 0.15
h = np.zeros(19)
v = np.zeros(19)

h = 90+coor_tilt[:,2]
v = v0*(h/h0)**shear
'''


#Input Transient
step = 3
tilt = np.radians(-5)
azimuth_inc = np.radians(-45)
azimuth= np.radians(0)
rpm=12.1

omega=[0,1.267,0]
g_init = [0,0,9.8]
g_tilt = np.zeros((step,3))
g_azimuth = np.zeros((step,3))
time=np.zeros(step)

for i in range (step):
    g_tilt[i,0] = g_init[0]
    g_tilt[i,1] = np.cos(-tilt)*g_init[1]-np.sin(-tilt)*g_init[2]
    g_tilt[i,2] = np.sin(-tilt)*g_init[1]+np.cos(-tilt)*g_init[2]

    g_azimuth[i,0] = np.cos(-azimuth)*g_tilt[i,0]-np.sin(-azimuth)*g_tilt[i,2]
    g_azimuth[i,1] = g_tilt[i,1]
    g_azimuth[i,2] = np.sin(-azimuth)*g_tilt[i,0]+np.cos(-azimuth)*g_tilt[i,2]
    azimuth=azimuth+azimuth_inc

    time[i]=(60/rpm)*(-np.degrees(azimuth)/360)
g_azimuth=np.round(g_azimuth,3)
print(g_azimuth)
fileOut='output/TransientCommand.txt'
f = open(fileOut, "w+")
f.write('nsubst,1,1,1\noutres,erase\noutres,all,all\nkbc,0\n\n')
for i in range(step):
    f.write('time,'+str(time[i]))
    f.write('\nacel,'+str(g_azimuth[i,0])+','+str(g_azimuth[i,1])+','+str(g_azimuth[i,2]))
    f.write('\nomega,'+str(omega[0])+','+str(omega[1])+','+str(omega[2]))
    f.write('\nlswrite,'+str(i+1))
    f.write('\n\n')

f.write('lssolve,1,'+str(step)+',1')

