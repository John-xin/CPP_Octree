import rhinoscriptsyntax as rs
import Rhino.Input as rip
import Rhino.Input.Custom as ric
#open point file and read
ptsList=[]
facesList=[]
ownerList=[]
nbrList=[]

indxList=[]
path=rs.WorkingFolder()

def readPtsFile(fileName):
    with open(fileName,"r") as f:
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

def readFacesFile(fileName):
    with open(fileName,"r") as f:
        for _ in range(11): #skip first 11 lines
            next(f)
        for line in f:
            data=(line.split(" "))
            if(len(data)>=7):
                nPts= int(data[0])
                face=[]
                for i in range(nPts):
                    face.append(int(data[i+2]))
                facesList.append(face)
        f.close()
    print(len(facesList))

def readOwnerFile(fileName):
    with open(fileName,"r") as f:
        lines = f.readlines()
        lines = lines[12:-1] #skip first 12 lines and remove last line
        for line in lines:
            data=(line.split(" "))
            if(len(data)==1):
                ownerList.append(int(data[0]))
        f.close()
    print(len(ownerList))
    
def readNbrFile(fileName):
    with open(fileName,"r") as f:
        lines = f.readlines()
        lines = lines[12:-1] #skip first 12 lines and remove last line
        for line in lines:
            data=(line.split(" "))
            if(data[0]!=")"):
                nbrList.append(int(data[0]))
        f.close()
    print(len(nbrList))



ptsFileName="../constant/polyMesh/points"
facesFileName="../constant/polyMesh/faces"
ownerFileName="../constant/polyMesh/owner"
nbrFileName="../constant/polyMesh/neighbour"
readPtsFile(ptsFileName)
readFacesFile(facesFileName)
readOwnerFile(ownerFileName)
readNbrFile(nbrFileName)


#draw cell

def drawCell(cellIndx):
    cellFaces=[]
    center=[]
    for i in range(len(ownerList)):
        if ownerList[i]==cellIndx:
            cellFaces.append(facesList[i])

    for i in range(len(nbrList)):
        if nbrList[i]==cellIndx:
            cellFaces.append(facesList[i])

    for face in cellFaces:
        points=[]
        for indx in face:
            points.append(ptsList[indx])
            rs.AddText(indx,ptsList[indx],0.2)
        points.append(points[0])
        rs.AddPolyline(points)

cellIndxList=[13,83 ]

for cellIndx in cellIndxList:
    drawCell(cellIndx)

 

