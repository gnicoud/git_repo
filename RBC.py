try: paraview.simple
except: from paraview.simple import *
import os
import os.path
import sys
from paraview import servermanager
from paraview import util
from paraview import vtk

#
# Useful functions !!!
#
def propagate(lisold,voisin,tagold,col):
#
# Propagate the front of colored nodes once
#
    num = len(lisold)
    tagnew = tagold[:]
    lisnew = lisold[:]
    for n in lisold:
        tmp = voisin[n]
        for m in tmp:
            if col[m]>0 and tagnew[m]==0:
                num = num+1
                lisnew.append(m)
                tagnew[m] = 1
    
    return lisnew, tagnew

def connexe(n0,voisin,col):
#
# Generate the list of nodes connected to n0 and colored
#
    coln0 = col[n0]
    nnode = len(col)
    lisold = []
    lisold.append(n0)
    tagold = [0 for n in col]
    tagold[n0] = 1
    var = 1

    while var != 0:
        [lisnew,tagnew] = propagate(lisold,voisin,tagold,col)
        lisold = lisnew[:]
        var = tagnew.count(1)-tagold.count(1)
        tagold = tagnew[:]

    return lisold
 
#
# End of useful functions
#

#INPUTFILE = "/home/nicoud/pyRBC/80GR/cell.sol000040.xmf"
INPUTFILE = sys.argv[1]
cell_sol000010_xmf = XDMFReader( FileName=INPUTFILE)

cell_sol000010_xmf.CellArrays = ['PROC_COLOR', 'SKEWNESS', 'TENSION', 'VOLUME_FLAG']
cell_sol000010_xmf.Sets = []
cell_sol000010_xmf.PointArrays = ['BND_FLAG', 'CURVATURE', 'CURVATURE_FORCE', 'ELASTIC_FORCE', 'EXTERNAL_FORCE', 'GAUSSIAN_CURVATURE', 'LAPL_CURV', 'MULTIPLICITY', 'SOLID_FORCE', 'SPONTANEOUS_CURVATURE', 'STRETCHING', 'TOTAL_FORCE', 'U', 'X_NODE_REF']
cell_sol000010_xmf.Grids = ['mesh']

# Restrict data to RBCs surface
surf = servermanager.filters.ExtractSurface(Input=cell_sol000010_xmf)
datasurf = servermanager.Fetch(surf)

# Store CURVATURE in curv arrays
datacurv1 = datasurf.GetPointData().GetArray("CURVATURE")
datacurv2 = datasurf.GetPointData().GetArray("GAUSSIAN_CURVATURE")

gauss_curv, curv = [], []
for i in range(datacurv1.GetNumberOfTuples()):
   #print(i,pcurv.GetTuple(i)[0])
   curv.append(datacurv1.GetTuple(i)[0])

for i in range(datacurv2.GetNumberOfTuples()):
   #print(i,pcurv.GetTuple(i)[0])
   gauss_curv.append(datacurv2.GetTuple(i)[0])

    
print('CURVATURE data read')    

# Put the points coordinates in x, y and z arrays
nnode = datasurf.GetNumberOfPoints()
x,y,z = [],[],[]
for i in range(nnode):
    coord = datasurf.GetPoint(i)
    xx, yy, zz = coord[:3]
    x.append(xx)
    y.append(yy)
    z.append(zz)

print('Coordinate data read')

# Put the node ids in n1, n2, n3 arrays
ntri = datasurf.GetNumberOfCells()
n1,n2,n3 = [], [], []
for i in range(ntri):
  n1.append(datasurf.GetCell(i).GetPointIds().GetId(0))
  n2.append(datasurf.GetCell(i).GetPointIds().GetId(1))
  n3.append(datasurf.GetCell(i).GetPointIds().GetId(2))

print('Cell2nodes data read')

#
# Generation des voisins des noeuds de surface
#
nbre_voisin = [0 for i in range(nnode)]
voisin = [[] for i in range(nnode)]

for tri in range(ntri):
    tn1 = n1[tri]
    tn2 = n2[tri]
    tn3 = n3[tri]
    
    istn2 = 0
    istn3 = 0
    for test in range(nbre_voisin[tn1]):
        if voisin[tn1][test] == tn2:
            istn2 = 1

        if voisin[tn1][test] == tn3:
            istn3 = 1 
  
    if istn2 == 0:
        nbre_voisin[tn1] = nbre_voisin[tn1] + 1
        voisin[tn1].append(tn2)
  
    if istn3 == 0:
        nbre_voisin[tn1] = nbre_voisin[tn1] + 1
        voisin[tn1].append(tn3)
  
    istn1 = 0
    istn3 = 0
    for test in range(nbre_voisin[tn2]):
        if voisin[tn2][test] == tn1:
            istn1 = 1

        if voisin[tn2][test] == tn3:
            istn3 = 1 
  
    if istn1 == 0:
        nbre_voisin[tn2] = nbre_voisin[tn2] + 1
        voisin[tn2].append(tn1)
  
    if istn3 == 0:
        nbre_voisin[tn2] = nbre_voisin[tn2] + 1
        voisin[tn2].append(tn3)

    istn1 = 0
    istn2 = 0
    for test in range(nbre_voisin[tn3]):
        if voisin[tn3][test] == tn1:
            istn1 = 1

        if voisin[tn3][test] == tn2:
            istn2 = 1 
  
    if istn1 == 0:
        nbre_voisin[tn3] = nbre_voisin[tn3] + 1
        voisin[tn3].append(tn1)
  
    if istn2 == 0:
        nbre_voisin[tn3] = nbre_voisin[tn3] + 1
        voisin[tn3].append(tn2)

print('Neighbor nodes computed')

# Compute the number of RBCs
flagRBC = [-1 for i in range(nnode)]
flag = -1

for n in range(nnode):
  if flagRBC[n]==-1:
    flag = flag + 1
    liste = connexe(n,voisin,[1 for i in range(nnode)])
    for nl in liste:
      flagRBC[nl] = flag                
        
number_of_RBC = flag+1
print('Number of RBCs detected: {}').format(number_of_RBC)

#for rbc in range(number_of_RBC):
#  print('RBC #',rbc,' : ',flagRBC.count(rbc),' nodes')
  
#print('number of non flaged nodes:',flagRBC.count(-1)) 



# Generation des taches a partir du signe de la courbure
col = [0 for n in range(nnode)]
for n in range(nnode):
    if (curv[n] < 0 and gauss_curv[n]>0):
        col[n] = 1

#
# Detection des regions coloriees connexes
#
number_of_region = [0 for rbc in range(number_of_RBC)]
number_of_region_filtered = [0 for rbc in range(number_of_RBC)]
size_of_region = [[] for rbc in range(number_of_RBC)]
tested = [0 for i in range(nnode)]
flag = 0

for n in range(nnode):
  if col[n]>0 and tested[n]==0:
    liste = connexe(n,voisin,col)
    for nl in liste:
      tested[nl] = 1
    number_of_region[flagRBC[n]] = number_of_region[flagRBC[n]] + 1
    size_of_region[flagRBC[n]].append(len(liste))
      
for rbc in range(number_of_RBC):
  toosmall = 5
  number_of_toosmall = 0
  for region in range(number_of_region[rbc]):  
    if size_of_region[rbc][region]<toosmall:
      number_of_toosmall = number_of_toosmall+1
  number_of_region_filtered[rbc] = number_of_region[rbc] - number_of_toosmall
  
  print('RBC #{}: {} nodes and {} concave zones containg {} nodes -- {} larger than {} nodes').format(rbc,flagRBC.count(rbc),number_of_region[rbc],size_of_region[rbc][:],number_of_region_filtered[rbc],toosmall)
  
print('==================')
print('STATISTICS')
print('  # of ellipsoids: {}').format(number_of_region_filtered.count(0))
print('  # of stomatocytes: {}').format(number_of_region_filtered.count(1))
print('  # of discocytes: {}').format(number_of_region_filtered.count(2))
print('  # of trilobes: {}').format(number_of_region_filtered.count(3))
print('  # of quadrilobes: {}').format(number_of_region_filtered.count(4))
total = number_of_region_filtered.count(0) + number_of_region_filtered.count(1) + number_of_region_filtered.count(2) + number_of_region_filtered.count(3) + number_of_region_filtered.count(4)
other = number_of_RBC - total
print('  # of others: {}').format(other)
  
  
