import rhinoscriptsyntax as rs
import Rhino.Input as rip
import Rhino.Input.Custom as ric
#open point file and read
ptsList=[]
indxList=[]
pt=[]
path=rs.WorkingFolder()
filename=rip.RhinoGet.GetFileName(ric.GetFileNameMode.OpenTextFile, "*", "select file",None)

with open(filename,"r") as f:
    for _ in range(1): #skip first 1 lines
        next(f)
    for line in f:
        data=(line.split(" "))
        pt=[float(data[1]),float(data[2]),float(data[3])]
        ptsList.append(pt)
        indxList.append(data[4])
    f.close()
print(len(ptsList))

#draw point
for i in range(len(ptsList)):
    pnt=rs.CreatePoint(ptsList[i])
    ptID=rs.AddPoint(pnt) 
    rs.AddText(str(indxList[i]),pnt,height=0.2)
 

