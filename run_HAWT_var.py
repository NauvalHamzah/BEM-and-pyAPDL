import numpy as np
import matplotlib.pyplot as plt
from pybemt.solver import Solver
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

# Run sweep of tip speed ratio with the BEMT method
input_name = 'NREL5MW2'

var1 = 'tsr'
low1 = 8.5
high1 = 9.5
interval1 = 0.1
n1 = int((high1-low1)/(interval1)+1)
var1list = np.linspace(low1,high1,n1)
temp = [None]*n1

var2 = 'twist'
low2 = -3
high2 = -1
interval2 = 0.1
n2 = int((high2-low2)/(interval2)+1)

P = np.zeros((n2,n1))

fin = open(input_name+'.ini','rt')
data1 = fin.readlines()[0:4]

for i in range (4):
    a = data1[i].find(var1)
    if a != -1:
        for j in range (n1):
            temp[j] = var1+ " = " + str(var1list[j])

fin.close()

for j in range (n1):
    f = open(input_name+'.ini', 'rt')
    data = f.read()
    index1 = data.find(var1)
    index2 = data.find("\n",index1)
    old = data[index1:index2]
    data = data.replace(old,temp[j])
    f.close()

    f = open(input_name+'.ini','wt')
    f.write(data)
    f.close()

    s = Solver(input_name+".ini",input_name)
    
    df, section_df = s.run_sweep(var2, n2, low2, high2)   
    P[:,j] = df['P']/1000

for i in range (n2):
    for j in range (n1):
        P[i,j] = round(P[i,j],3)


maxP = np.amax(P)
maxP_Loc = np.where(P==maxP)
var1_max = maxP_Loc[1]
var2_max = maxP_Loc[0]


x = np.linspace(low1, high1, n1)
y = np.linspace(low2, high2, n2)

X, Y = np.meshgrid(x, y)

'''
fig,ax = plt.subplots(1,1)
cp = ax.contourf(X, Y, P,50,cmap='jet')
fig.colorbar(cp,cmap='jet')
ax.set_title("Power")
ax.set_xlabel(var1)
ax.set_ylabel(var2)
ax.grid(color='black')
ax.plot(X[var2_max,var1_max],Y[var2_max,var1_max],'o',color='blue')
ax.text(X[var2_max,var1_max],Y[var2_max,var1_max],'P = ' + str(maxP)+' kW\nTwist = ' +str(Y[var2_max,var1_max])+' deg\nTSR = '+str(X[var2_max,var1_max]))
plt.show()

'''
fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, P, rstride=1, cstride=1, cmap='jet',linewidth=0, antialiased=False)
ax.view_init(30, 225)
ax.set_xlabel(var1)
ax.set_ylabel(var2)
pyplot.show()


