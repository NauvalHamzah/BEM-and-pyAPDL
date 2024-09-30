import numpy
import matplotlib.pyplot as plt
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
from pybemt.solver import Solver



file = open("pybemt/data.txt",'r')
lines = file.readlines() 

file1 = open("pybemt/data.txt",'r')
a = file1.read()
a = a.split()
nx = (int(len(a)/(len(lines)*2)))
ny = len(lines)
x = numpy.linspace(-10, 10, nx)
y = numpy.linspace(7, 12, ny)

pitch, rpm = numpy.meshgrid(x, y)

Power = numpy.zeros((ny,nx))
MaxPower = 0

for i in range (ny):
    n = i*2*nx+1
    for j in range (nx):
        Power[i,j] = float(a[n+j*2])/1000
        if MaxPower<Power[i,j]:
            MaxPower=Power[i,j]
            max_i = i
            max_j = j

pmax=numpy.amax(Power)
print(rpm[max_i,max_j])
print(pmax)

fig,ax = plt.subplots(1,1)
cp = ax.contourf(pitch, rpm, Power,50,cmap='jet')
fig.colorbar(cp,cmap='jet')
ax.set_title("Power")
ax.set_xlabel("Pitch")
ax.set_ylabel("RPM")
plt.show()


'''
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(pitch, rpm, Power, rstride=1, cstride=1, cmap='jet',linewidth=0, antialiased=False)
ax.set_xlim(-10, 10)
ax.set_ylim(7, 12)
ax.set_zlim(0,pmax)
ax.view_init(30, 225)
ax.set_xlabel('$Pitch$')
ax.set_ylabel('$RPM$')
pyplot.show()
'''
