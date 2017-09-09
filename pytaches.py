
# -*- coding: cp1252 -*-
import numpy as np

#
# Definition des fonctions utiles
#

def propagate(lisold,voisin,tagold,col):
#
# Propage de 1 niveau le front des noeuds colories
#
    num = len(lisold)
    tagnew = tagold[:]
    lisnew = lisold[:]
    for n in lisold:
        tmp = voisin[n-1]
        for m in tmp:
            if col[m-1]==1 and tagnew[m-1]==0:
                num = num+1
                lisnew.append(m)
                tagnew[m-1] = 1
    
    return lisnew, tagnew

def connexe(n0,voisin,col):
#
# Cree la listes des noeuds colories et connectés à n0
#
    nnode = len(col)
    lisold = []
    lisold.append(n0+1)
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
# Fin de definition des fonctions utiles ----------------------------
#


#
# Lecture du maillage
#
f = open('RBC_002.msh', 'r')

for line in range(0,6):
    tmp=f.readline()
tmp = f.readline().split()
nnodeVol=int(tmp[3],16)
print('nnode de vol = ',nnodeVol)

x , y , z = [] , [] , []
for line in range(0,nnodeVol):
    tmp = f.readline().split()
    x.append(float(tmp[0]))
    y.append(float(tmp[1]))
    z.append(float(tmp[2]))

for line in range(0,4):
    tmp=f.readline()
tmp = f.readline().split()
ntri=int(tmp[2],16)
print('ntri = ',ntri)

n1 , n2 , n3 = [] , [] , []
for line in range(ntri):
    tmp = f.readline().split()
    n1.append(int(tmp[1],16))
    n3.append(int(tmp[2],16))
    n2.append(int(tmp[3],16))

f.close()

#
# Generation des voisins des noeuds de surface
#
nbre_voisin = [0 for i in range(nnodeVol)]
voisin = [[] for i in range(nnodeVol)]

for tri in range(ntri):
    tn1 = n1[tri]
    tn2 = n2[tri]
    tn3 = n3[tri]
    
    istn2 = 0
    istn3 = 0
    for test in range(nbre_voisin[tn1-1]):
        if voisin[tn1-1][test] == tn2:
            istn2 = 1

        if voisin[tn1-1][test] == tn3:
            istn3 = 1 
  
    if istn2 == 0:
        nbre_voisin[tn1-1] = nbre_voisin[tn1-1] + 1
        voisin[tn1-1].append(tn2)
  
    if istn3 == 0:
        nbre_voisin[tn1-1] = nbre_voisin[tn1-1] + 1
        voisin[tn1-1].append(tn3)
  
    istn1 = 0
    istn3 = 0
    for test in range(nbre_voisin[tn2-1]):
        if voisin[tn2-1][test] == tn1:
            istn1 = 1

        if voisin[tn2-1][test] == tn3:
            istn3 = 1 
  
    if istn1 == 0:
        nbre_voisin[tn2-1] = nbre_voisin[tn2-1] + 1
        voisin[tn2-1].append(tn1)
  
    if istn3 == 0:
        nbre_voisin[tn2-1] = nbre_voisin[tn2-1] + 1
        voisin[tn2-1].append(tn3)

    istn1 = 0
    istn2 = 0
    for test in range(nbre_voisin[tn3-1]):
        if voisin[tn3-1][test] == tn1:
            istn1 = 1

        if voisin[tn3-1][test] == tn2:
            istn2 = 1 
  
    if istn1 == 0:
        nbre_voisin[tn3-1] = nbre_voisin[tn3-1] + 1
        voisin[tn3-1].append(tn1)
  
    if istn2 == 0:
        nbre_voisin[tn3-1] = nbre_voisin[tn3-1] + 1
        voisin[tn3-1].append(tn2)

#
# Recuperer les nodes de bords seulement
#
nbr = [n for n in nbre_voisin if n!=0]
vois = [liste for liste in voisin if liste!=[]]
nnode = len(nbr)

#
# Generation des taches
#
col = [0 for n in range(nnode)]
Imax = y.index(max(y))
raycrit = np.sqrt(x[Imax]**2+z[Imax]**2)
raymax = max(x)**2
for n in range(nnode):
    dist = np.sqrt(x[n]**2+z[n]**2)
    if dist < raycrit:
        col[n] = 1
#    if x[n]>3.8:
#        col[n]=1

#
# Detection des regions coloriees connexes
#
number_of_region = []
for RBC in range(100):
    print('RBC = ',RBC)
    lenlist = []
    tested = [0 for i in range(nnode)]
    flag = 0

    for n in range(nnode):
        if col[n]==1 and tested[n]==0:
            flag = flag + 1
            liste = connexe(n,vois,col)
            lenlist.append(len(liste))
            for nl in liste:
                tested[nl-1] = flag                
        
number_of_region.append(flag)
print(lenlist)

