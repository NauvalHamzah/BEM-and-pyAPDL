import numpy as np 

file = open("Displacement.txt",'r')
lines= file.readlines()[10:-3]
row=len(lines)

for i in range(row):
    lst=[]
    for pos,char in enumerate(lines[i]):
        if(char=="E"):
            lst.append(pos)
    for j in range(len(lst)):
        pos=lst[len(lst)-1-j]       
        if lines[i][pos+1]=="-" and lines[i][pos+5]=="-":
            new=lines[i][0:pos+5]+" "+lines[i][pos+5:]
            lines[i]=new

node_disp=[]

for i in range(row):
    col=[]
    a=lines[i].split()
    for j in range(5):
        if (len(a))==5:
            b=a[j]
            col.append(b)
    if (len(a))==5:
        if col[0] != "NODE" and col[0] != "LOAD" and col[0] != "TIME=" :
            node_disp.append(col)
        elif col[0] == "NODE" and i==0:
            node_disp.append(col)

real_row=len(node_disp)

f = open("Out_Displacement.txt", "w+")
for i in range(real_row):
    f.write(str(node_disp[i][0])+ " " + str(node_disp[i][1])+ " " + str(node_disp[i][2])+ " " + str(node_disp[i][3])+ " " + str(node_disp[i][4]))
    f.write("\n")
