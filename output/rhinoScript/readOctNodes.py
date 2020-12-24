import rhinoscriptsyntax as rs
import Rhino.Input as rip
import Rhino.Input.Custom as ric
import Rhino
import scriptcontext
import System.Guid


ptsList=[]
nodesList=[]

class node:
    def __init__(self, _indxList, _state):
        self.indxList = _indxList
        self.state = _state



def readPts():
    path=rs.WorkingFolder()
    filename=rip.RhinoGet.GetFileName(ric.GetFileNameMode.OpenTextFile, "*", "select file",None)

    with open(filename,"r") as f:
        for _ in range(11): #skip first 11 lines
            next(f)
        for line in f:
            data=(line.split(" "))
            if(len(data)==5):
                pt=[float(data[1]),float(data[2]),float(data[3])]
                ptsList.append(pt)
                #indxList.append(data[4])
        f.close()
    print(len(ptsList))

def readNodes():
    path=rs.WorkingFolder()
    filename=rip.RhinoGet.GetFileName(ric.GetFileNameMode.OpenTextFile, "*", "select file",None)

    with open(filename,"r") as f:
        for _ in range(1): #skip first 1 lines
            next(f)
        for line in f:
            dataTMP=(line.split(" "))
            data=dataTMP[2].split(",")
            if(len(data)==9):
                indxList=[int(data[0]),int(data[1]),int(data[2]),int(data[3]),
                               int(data[4]),int(data[5]),int(data[6]),int(data[7])]
                state= int(dataTMP[3])
                #node= node(indxList,state)
                nodesList.append(node(indxList,state))
                #indxList.append(data[4])
        f.close()
    print(len(nodesList))

def AddBrepBox(pt1, pt2):
    box = Rhino.Geometry.BoundingBox(pt1, pt2)
    brep=box.ToBrep()
    rc = Rhino.Commands.Result.Failure
    if( scriptcontext.doc.Objects.AddBrep(brep) != System.Guid.Empty ):
        rc = Rhino.Commands.Result.Success
        scriptcontext.doc.Views.Redraw()
    return rc

def addBox(node):
    box=[]
    for indx in node.indxList:
        box.append(ptsList[indx])
    boxBody=rs.AddBox(box)
    if (node.state==1):
        rs.ObjectColor(boxBody, [255,0,0])
    elif(node.state==3):
        rs.ObjectColor(boxBody, [0,255,0])
    elif(node.state==4):
        rs.ObjectColor(boxBody, [0,0,0])


#AddBrepBox(Rhino.Geometry.Point3d(0, 0, 0),Rhino.Geometry.Point3d(1, 1, 1))
#sp=rs.AddSphere([0,0,0],1)
#rs.ObjectColor(sp, [255,0,0])
#rs.ObjectLayer(objects, layername)
readPts()
readNodes()
for node in nodesList:
    addBox(node)
