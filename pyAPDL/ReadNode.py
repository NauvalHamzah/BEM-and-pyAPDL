import numpy as np 
import math
from ansys.mapdl.core import launch_mapdl

def NodeData(fileIn,fileOut):
    file = open(fileIn,'r')
    lines= file.readlines()[4:]

    row=len(lines)
    node=[]

    #Split Spacing
    for i in range(row):
        col=[]
        a=lines[i].split()
        for j in range(4):
            if (len(a))==4:
                b=a[j]
                col.append(b)
        if (len(a))==4:
            if col[0] != "NODE":
                node.append(col)
            elif col[0] == "NODE" and i==0:
                node.append(col)

    real_row=len(node)

    #Save Out file
    f = open(fileOut, "w+")
    for i in range(real_row):
        f.write(str(node[i][0])+ " " + str(node[i][1])+ " " + str(node[i][2])+ " " + str(node[i][3]))
        f.write("\n")

    #Node in float type
    Node_float = np.zeros((real_row-1,4))
    for i in range (1,real_row):
        for j in range (4):
            Node_float[i-1,j]=float(node[i][j])
    MaxNode=Node_float[len(Node_float)-1,0]
    return Node_float,MaxNode

def StatData(StatList):
    f=open(StatList, "r")
    stat=f.read()
    stat=stat.split()
    nstat=int(len(stat)/4)
    r=np.zeros(nstat)
    twist=np.zeros(nstat)
    c=np.zeros(nstat)
    ac=np.zeros(nstat)

    for i in range (nstat):
        r[i] = stat[i*4]
        twist[i] = stat[i*4+1]
        c[i] = stat[i*4+2]
        ac[i] = stat[i*4+3]

    return r,twist,c,ac

def StatSort(r,Nodes_Per_Element,Node_float):
    nstat=len(r)
    stat_r=np.zeros(nstat)
    for i in range(nstat):
        stat_r[i]=float(r[i])

    stat_r[nstat-1]=stat_r[nstat-1]+0.1
    sortZ=Node_float[Node_float[:,3].argsort()]
    

    #Save nodes coordinates per elemen
    f = open(Nodes_Per_Element, "w+")
    for i in range(nstat-1):
        f.write("Element number "+str(i+1))
        f.write("\n")
        for j in range (len(sortZ)):
            if sortZ[j,3]>=stat_r[i] and sortZ[j,3]<stat_r[i+1]:
                f.write("\t"+str(sortZ[j,0])+" "+str(sortZ[j,1])+" "+str(sortZ[j,2])+" "+str(sortZ[j,3])+" ")
                f.write("\n")
        f.write("\n\n")
    return sortZ

def NodeAtStat(r,twist,fStatNode,NodeFloatZ):
    nstat=len(r)
    stat_r=np.zeros(nstat)
    stat_twist=np.zeros(nstat)
    
    for i in range(nstat):
        stat_r[i]=float(r[i])
        stat_twist[i]=float(twist[i])

    
    NodeStat=[]

    #Save nodes coordinates at station
    f = open(fStatNode, "w+")
    for i in range(nstat):
        f.write("Station "+str(i+1))
        f.write("\n")
        for j in range (len(NodeFloatZ)):
            if NodeFloatZ[j,3]>(stat_r[i])-0.01 and NodeFloatZ[j,3]<stat_r[i]+0.01:
                NodeStat.append(NodeFloatZ[j])
                f.write("\t"+str(NodeFloatZ[j,0])+" "+str(NodeFloatZ[j,1])+" "+str(NodeFloatZ[j,2])+" "+str(NodeFloatZ[j,3])+" ")
                f.write("\n")
        f.write("\n\n")
    
    #Read nodes per station from output files
    StatNodeSize=np.zeros(nstat)
    f=open(fStatNode,"r")
    a=f.readlines()   
    for i in range(len(a)):a[i]=a[i].split()
    temp=[]
    for i in range(len(a)):
        if len(a[i])>0:
            temp.append(a[i])
    
    row=len(temp)
    lst=[]
 
    for i in range(row):
        char=temp[i][0]
        if char=="Station":
            lst.append(i)
    
    #Number of nodes per station
    StatLen=np.zeros(len(lst))

    for i in range (len(lst)-1):
        StatLen[i]=lst[i+1]-lst[i]-1
    StatLen[len(lst)-1]=len(temp)-1-lst[len(lst)-1]
    #print(StatLen)
    
    #Nodes per station
    NodeFloat3d=np.zeros((nstat,int(max(StatLen)),4))
    
    n_sum=0   
    for i in range(nstat):
        for j in range (int(StatLen[i])):
            index=int(n_sum+j)
            for k in range (4):
                NodeFloat3d[i,j,k]=NodeStat[index][k]
        n=StatLen[i]
        n_sum=n_sum+n
    
    
    #NodeFloat3d No Twist
    NodeFloat3d_NoTwist=np.zeros((nstat,int(max(StatLen)),4))
    for i in range(nstat):
        for j in range (int(StatLen[i])):
            for k in range (4):
                NodeFloat3d_NoTwist[i,j,k]=NodeFloat3d[i,j,k]

    
    for i in range(nstat):
        for j in range (int(StatLen[i])):
            NodeFloat3d_NoTwist[i,j,1]=NodeFloat3d[i,j,1]*np.cos(np.radians(stat_twist[i]))-NodeFloat3d[i,j,2]*np.sin(np.radians(stat_twist[i]))
            NodeFloat3d_NoTwist[i,j,2]=NodeFloat3d[i,j,1]*np.sin(np.radians(stat_twist[i]))+NodeFloat3d[i,j,2]*np.cos(np.radians(stat_twist[i]))

    #for i in range(nstat):
    #    np.savetxt(str(i)+".txt",NodeFloat3d_NoTwist[i,:,:])

    LeTeNodes=np.zeros((2*nstat,4))
    for i in range (nstat):
        a=np.zeros((int(StatLen[i]),4))
        b=np.zeros((int(StatLen[i]),4))
        for j in range (int(StatLen[i])):
            for k in range (4):
                a[j,k]=NodeFloat3d_NoTwist[i,j,k]
                b[j,k]=NodeFloat3d[i,j,k]

        Te=0
        Le=0
        TeIndex=0
        LeIndex=0

        for l in range (int(StatLen[i])):
            if a[l,1]>0:
                if Le<a[l,1]:
                    Le=a[l,1]
                    LeIndex=l                   
            
            elif a[l,1]<0:
                if Te>(a[l,1]):
                    Te=(a[l,1])
                    TeIndex=l 

        LeTeNodes[i,:]=b[TeIndex,:]
        LeTeNodes[nstat+i,:]=b[LeIndex,:]
        

    
    StatChord=np.zeros(nstat)
    StatAngle=np.zeros(nstat)
    for i in range (nstat):
        x=abs(LeTeNodes[i,1]-LeTeNodes[nstat+i,1])
        y=LeTeNodes[i,2]-LeTeNodes[nstat+i,2]
        StatChord[i]=np.sqrt(x**2+y**2)
        StatAngle[i]=np.degrees(np.arctan2(y,x))
                 
    return(StatLen,NodeStat,LeTeNodes,StatChord,StatAngle)

def LoadPoint(LeTeNodes,LoadInput,LoadOut,twist,c,ac,tilt,azimuth):
    file = open(LoadInput,'r')
    lines= file.readlines()[5:]
    row=len(lines)
    Fn=np.zeros(row)
    Ft=np.zeros(row)
    M=np.zeros(row)
    Xl=np.zeros(row)
    Yl=np.zeros(row)
    Yt=np.zeros(row)
    NodeNumber=np.array(range(2*row)).reshape((row,2))
    Load=np.zeros((row,6))
    nstat = int(len(LeTeNodes)/2)

    for i in range(row):
        a=lines[i].split()
        Fn[i]=a[1]   
        Ft[i]=a[2]
        M[i]=a[3]

        Yt[i] = (Fn[i]+(-M[i]+(np.sin(np.radians(twist[i+1]))*ac[i+1]*c[i+1]*Ft[i]))/(ac[i+1]*c[i+1]*np.cos(np.radians(twist[i+1]))))/(1+(1-ac[i+1])/ac[i+1])
        Yl[i] = Fn[i]-Yt[i]
        Xl[i] = Ft[i]

        NodeNumber[i,0] = int(LeTeNodes[nstat+i+1,0])
        NodeNumber[i,1] = int(LeTeNodes[i+1,0])
        
    YtTilt = np.zeros((len(Yl),3))
    YtAzim = np.zeros((len(Yt),3))

    YlTilt = np.zeros((len(Yl),3))
    YlAzim = np.zeros((len(Yt),3))

    YtTilt[:,1] = np.cos(-tilt)*Yt
    YtTilt[:,2] = np.sin(-tilt)*Yt

    YlTilt[:,1] = np.cos(-tilt)*Yl
    YlTilt[:,2] = np.sin(-tilt)*Yl

    YtAzim[:,0] = -np.sin(-azimuth)*YtTilt[:,2]
    YtAzim[:,1] = YtTilt[:,1]
    YtAzim[:,2] = np.cos(-azimuth)*YtTilt[:,2]

    YlAzim[:,0] = -np.sin(-azimuth)*YlTilt[:,2]
    YlAzim[:,1] = YlTilt[:,1]
    YlAzim[:,2] = np.cos(-azimuth)*YlTilt[:,2]

    LeLoad = np.zeros((len(Yl),3))
    TeLoad = np.zeros((len(Yt),3))

    LeLoad[:,0] = np.round((YlAzim[:,0] + Xl),3)
    LeLoad[:,1] = np.round(YlAzim[:,1],3)
    LeLoad[:,2] = np.round(YlAzim[:,2],3)

    TeLoad[:,0] = np.round(YtAzim[:,0],3)
    TeLoad[:,1] = np.round(YtAzim[:,1],3)
    TeLoad[:,2] = np.round(YtAzim[:,2],3)

    direction = ['Fx','Fy','Fz','Fx','Fy','Fz']
    Load[:,0:3] = LeLoad
    Load[:,3:6]  = TeLoad

    f = open(LoadOut,'w+')
    for i in range (6): 
        if i != 0:f.write('\n\n')   
        for j in range (row):   
            if i>2: a=0
            else: a=1
            f.write('F,'+str(NodeNumber[j,a])+','+str(direction[i])+','+str(Load[j,i])+'\n')
    f.close()
    return     

def ReadDisplacement(fileIn,fileOut,LeTeNodes):
    file = open(fileIn,'r')
    lines= file.readlines()[10:-7]
    row=len(lines)

    
    for i in range(row):
        lst=[]
        for pos,char in enumerate(lines[i]):
            if(char=="E"):
                lst.append(pos)
       
        for j in range(len(lst)):
            pos=lst[len(lst)-1-j]       
            if lines[i][pos+5]=="-":
                new=lines[i][0:pos+5]+" "+lines[i][pos+5:]
                lines[i]=new
    
    node_disp=[]

    for i in range(row):
        col=[]
        a=lines[i].split()
        for j in range(6):
            if (len(a))==6:
                b=a[j]
                col.append(b)
        if (len(a))==6:
            if col[0] != "NODE" and col[0] != "LOAD" and col[0] != "TIME=" and col[0] != "*****":
                node_disp.append(col)
            elif col[0] == "NODE" and i==0:
                node_disp.append(col)

    real_row=len(node_disp)

    f = open(fileOut, "w+")
    for i in range(real_row):
        f.write(str(node_disp[i][0])+ " " + str(node_disp[i][1])+ " " + str(node_disp[i][2])+ " " + str(node_disp[i][3])+ " " + str(node_disp[i][4])+ " " + str(node_disp[i][5]))
        f.write("\n")

    #Node in float type
    Disp_float = np.zeros((real_row-1,5))
    for i in range (1,real_row):
        for j in range (5):
            Disp_float[i-1,j]=float(node_disp[i][j])
    
    LeTeDisp=np.zeros((len(LeTeNodes),5))
    LeTeLoc=LeTeNodes[:,0]
    LeTeIndex=np.zeros(len(LeTeNodes))

    for i in range(len(LeTeLoc)):
        for j in range (len(Disp_float)):
            if abs(Disp_float[j,0]-LeTeLoc[i])<0.1:
                LeTeIndex[i]=j
                LeTeDisp[i]=Disp_float[j]

    return (Disp_float,LeTeDisp)

def Deformed(LeTeNodes,LeTeDisp):
    LeTeDeformed=LeTeNodes
    for i in range (len(LeTeNodes)):
        for j in range (1,4):
            LeTeDeformed[i,j]=LeTeNodes[i,j]+LeTeDisp[i,j]
    
    nstat=int(len(LeTeNodes)/2)
    StatAngle=np.zeros(nstat)
    for i in range (nstat):
        x=abs(LeTeDeformed[i,1]-LeTeDeformed[nstat+i,1])
        y=LeTeDeformed[i,2]-LeTeDeformed[nstat+i,2]
        StatAngle[i]=round(np.degrees(np.arctan2(y,x)),3)

    return(LeTeDeformed,StatAngle)
    
def LoadCommand(LoadOut,step,rpm,tilt,azimuth,azimuth_inc,time,i,g_tilt,g_azimuth):
    omega=[0,round((rpm*2*np.pi/60),3),0]
    g_init = [0,0,9.8]
    
    g_tilt[i,0] = g_init[0]
    g_tilt[i,1] = np.cos(-tilt)*g_init[1]-np.sin(-tilt)*g_init[2]
    g_tilt[i,2] = np.sin(-tilt)*g_init[1]+np.cos(-tilt)*g_init[2]

    g_azimuth[i,0] = round(np.cos(-azimuth)*g_tilt[i,0]-np.sin(-azimuth)*g_tilt[i,2],3)
    g_azimuth[i,1] = round(g_tilt[i,1],3)
    g_azimuth[i,2] = round(np.sin(-azimuth)*g_tilt[i,0]+np.cos(-azimuth)*g_tilt[i,2],3)
    

    if i == 0: time[i] = 1
    else : time[i]=round(1+(60/rpm)*(-np.degrees(azimuth)/360),3)
    azimuth=azimuth+azimuth_inc
    

    file = open(LoadOut,'r')
    lines= file.read()
    fileOut='output/TransientCommand.txt'
    f = open('output/TransientCommand.txt', "w+")
    f.write('/solu\nantype,4,new\nnsubst,1,1,1\noutres,erase\noutres,all,all\nkbc,0\n\n')
    
    f.write('time,'+str(time[i]))
    f.write('\nacel,'+str(g_azimuth[i,0])+','+str(g_azimuth[i,1])+','+str(g_azimuth[i,2]))
    f.write('\nomega,'+str(omega[0])+','+str(omega[1])+','+str(omega[2])+'\n')
    f.write(lines)
    #f.write('\nlswrite,'+str(i+1))
    f.write('\nlswrite,'+str(1))
    f.write('\n\n')
    f.write('lssolve,1,1,1\n')
    f.close()

    f=open('output/OutCommand.txt',"w+")
    f.write('/post1\n/output,F:/Kuliah/Tesis/Program/ANSYS_APDL/NREL5MW/Displacement-'+str(i+1)+'.txt\nPRNSOL,U,SUM\n/out')
    #f.write('lssolve,1,'+str(step)+',1')

    return azimuth

def RunMAPDL (dbfile,input_file,output_file):
    mapdl = launch_mapdl()
    # clear
    mapdl.finish()
    mapdl.clear()
    
    mapdl.filname('Blade')
    # resume dbfile
    resume_command = 'resume,'+dbfile+',db'
    mapdl.run(resume_command)
    mapdl.run(input_file)
    mapdl.run(output_file)
    mapdl.save('F:/Kuliah/Tesis/Program/ANSYS_APDL/NREL5MW/Blade','db')

    return
