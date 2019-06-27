###########################################################
#                                                         #
#      __                                         __      #
#  ___( o)>  Author  : Léo GASPARD              <(o )___  #
#  \ <_. )   Twitter : @leo_gaspard              ( ._> /  #
#   `---'    Mail    : leo.gaspard@outlook.fr     `---'   #
#                                                         #
#                                                         #
###########################################################

import scipy.optimize as optimize
import numpy as np
import os
import operator
import sys
import datetime

###### GLOBAL VARIABLES #######
dr = 0.001
da = 0.01
DA = 1.5

rBath = 'x'
rPP = 'x'
center = 'x'
xAxis = 'x'
yAxis = 'x'
zAxis = 'x'
sym = 'x'
output_file = 'x'
pattern = []
npattern = []
atoms = 'x'
dist = 'x'
lattice = 'x'
a = 'x'
b = 'x'
c = 'x'
alpha = 'x'
beta = 'x'
gamma = 'x'
visu = 0
seefrag = 0
opti = 0
trsl = 'x'
notIn = 'x'
evj = 0
norep = 0
symop = []
generator = []
noinfrag = []

##############################
def big_cell(na,nb,nc):
    coords = []
    newCoords = []
    newNewCoords = []

    for i in generator:
        x = i[1]
        y = i[2]
        z = i[3]

        for j in symop:
            coords.append([eval(j[0]),eval(j[1]),eval(j[2]),i[0]])
    
    for i in coords:
        i[0] *= a
        i[1] *= b
        i[2] *= c
        if i[0] > a:
            i[0] -= a
        if i[1] > b:
            i[1] -= b
        if i[2] > c:
            i[2] -= c
        if i[0] < 0:
            i[0] += a
        if i[1] < 0:
            i[1] += b
        if i[2] < 0:
            i[2] += c
    for i in coords:
        if i not in newCoords:
            if i[0] < a and i[1] < b and i[2] < c:
                newCoords.append(i)
    
    for i in newCoords:
        newNewCoords.append(i)
        for j in range(1,na):
            newNewCoords.append([i[0]+a*j,i[1],i[2],i[3]])
            for k in range(1,nb):
                newNewCoords.append([i[0]+a*j,i[1]+b*k,i[2],i[3]])
                for l in range(1,nc):
                    newNewCoords.append([i[0]+a*j,i[1]+b*k,i[2]+c*l,i[3]])
            for k in range(1,nc):
                newNewCoords.append([i[0]+a*j,i[1],i[2]+c*k,i[3]])
        for j in range(1,nb):
            newNewCoords.append([i[0],i[1]+b*j,i[2],i[3]])
            for k in range(1,nc):
                newNewCoords.append([i[0],i[1]+b*j,i[2]+c*k,i[3]])
        for j in range(1,nc):
            newNewCoords.append([i[0],i[1],i[2]+c*j,i[3]])

    return newNewCoords

def printProgressBar (start, now, iteration, total, prefix = '', suffix = '', decimals = 3, length = 100, fill = '\u2588'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    
    dif = now-start
    
    print('\r%30s |%s| %7s%% %s Elapsed time : %s   ' % (prefix, bar, percent, suffix,str(dif)),end='\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

def read_input(inputFile):
    global rBath
    global rPP
    global center
    global xAxis
    global yAxis
    global zAxis
    global sym
    global output_file
    global pattern
    global npattern
    global atoms
    global dist
    global lattice
    global a
    global b
    global c
    global alpha
    global beta
    global gamma
    global visu
    global evj
    global seefrag
    global opti
    global trsl
    global notIn
    global norep
    global symop
    global generator
    f = open(inputFile,'r')
    line = 'x'
    
    while line != ['END_OF_INPUT']:
        line = f.readline()
        line = line.split()
        if line == []:
            continue
        elif line[0] == 'BATH':
            rBath = float(line[1])
            print("Bath radius : %f\n"%(rBath))
        elif line[0] == 'PSEUDO':
            rPP = float(line[1])
            print("First Shell radius : %f\n"%(rPP))
        elif line[0] == 'CENTER':
            center = []
            for i in range(1,len(line)):
                center.append(line[i])
        elif line[0] == 'X_AXIS':
            xAxis = []
            for i in range(1,len(line)):
                xAxis.append(line[i])
        elif line[0] == 'Y_AXIS':
            yAxis = []
            for i in range(1,len(line)):
                yAxis.append(line[i])
        elif line[0] == 'Z_AXIS':
            zAxis = []
            for i in range(1,len(line)):
                zAxis.append(line[i])
        elif line[0] == 'SYMMETRY':
            print("Will treat symmetry")
            sym = []
            for i in range(1,len(line)):
                sym.append(line[i])
        elif line[0] == 'OUTPUT':
            output_file = line[1]
        elif line[0] == 'PATTERN':
            pattern = []
            line = f.readline()
            while line.strip() != 'NRETTAP':
                pattern.append([])
                line = line.split()
                for i in range(len(line)):
                    if i%2 == 0:
                        pattern[-1].append(int(line[i]))
                    else:
                        pattern[-1].append(line[i])
                line = f.readline()
        elif line[0] == 'NPATTERN':
            for i in range(1,len(line)):
                npattern.append(int(line[i]))
        elif line[0] == 'LATTICE':
            a = float(f.readline().split()[1])
            b = float(f.readline().split()[1])
            c = float(f.readline().split()[1])
            alpha = float(f.readline().split()[1])
            beta = float(f.readline().split()[1])
            gamma = float(f.readline().split()[1])

            print("Lattice parameter : \na     = %f \nb     = %f \nc     = %f \nalpha = %f \nbeta  = %f \ngamma = %f \n"%(a,b,c,alpha,beta,gamma))
        elif line[0] == 'ATOMS':
            atoms = []
            for i in range(1,len(line)):
                if i%4 == 1:
                    atoms.append(line[i])
                elif i%4 == 2:
                    atoms.append(float(line[i]))
                elif i%4 == 3:
                    atoms.append(int(line[i]))
                elif i%4 == 0:
                    atoms.append(float(line[i]))
        elif line[0] == 'DIST':
            dist = []
            line = f.readline()
            while line.strip() != 'TSID':
                line = line.split()
                dist.append([line[0],line[1],float(line[2])])
                line = f.readline()
        elif line[0] == 'COLOR':
            visu = 1
        elif line[0] == 'NOCOLOR':
            visu = 2
        elif line[0] == 'TRANSLATE':
            trsl = [float(line[1]),float(line[2]),float(line[3])]
        elif line[0] == 'NOTINPP':
            notIn = []
            for i in range(1,len(line)):
                notIn.append(line[i])
        elif line[0] == 'OPTIMIZATION':
            opti = 1
        elif line[0] == 'SEEFRAG':
            seefrag = 1
        elif line[0] == 'EVJEN':
            evj = 1
        elif line[0] == 'NOREP':
            norep = 1
        elif line[0] == 'SYMOP':
            line = f.readline()
            while line.strip() != 'POMYS':
                symop.append(line.split(','))
                line = f.readline()
        elif line[0] == 'GENERATOR':
            line = f.readline()
            while line.strip() != 'ROTARENEG':
                gen = line.split()
                generator.append([gen[0],float(gen[1]),float(gen[2]),float(gen[3])])
                line = f.readline()
        elif line[0] == 'NOINFRAG':
            line = f.readline()
            while line.strip() != 'GARFNION':
                noin = line.split()
                noinfrag.append([noin[0],float(noin[1]),float(noin[2]),float(noin[3])])
                line = f.readline()

def parse(fileName):
    f = open(fileName,'r')
    line = 0
    dataTable = []

    f.readline()
    f.readline()

    while True:
        line = f.readline().strip().split()
        if line == [] :
            break
        line[1] = float(line[1])
        line[2] = float(line[2])
        line[3] = float(line[3])

        dataTable.append(line)

    return dataTable

def distance(u,v):  #Return the distance between two given points
    return np.sqrt((u[0]-v[0])**2+(u[1]-v[1])**2+(u[2]-v[2])**2)

def vect_product(u,v): #Return the vectorial product between two vectors
    return [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]

def dot_product(u,v): #Return the dot product between two vectors
    return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2])

def normalize(u): #Normalize a vector
    norm = np.sqrt(u[0]**2+u[1]**2+u[2]**2)
    if norm != 0:
        return [u[0]/norm, u[1]/norm, u[2]/norm]
    else:
        return u

def rot_matrix(oldAxis,newAxis): #Build a rotation matrix to change director vector
    newAxis = normalize(newAxis)
    vp = normalize(vect_product(oldAxis,newAxis))
    angle = np.arccos(dot_product(oldAxis,newAxis))

    rMat = [[np.cos(angle)+vp[0]*vp[0]*(1-np.cos(angle)),vp[0]*vp[1]*(1-np.cos(angle))-vp[2]*np.sin(angle),vp[0]*vp[2]*(1-np.cos(angle))+vp[1]*np.sin(angle)],[vp[0]*vp[1]*(1-np.cos(angle))+vp[2]*np.sin(angle),vp[1]*vp[1]*(1-np.cos(angle))+np.cos(angle),vp[1]*vp[2]*(1-np.cos(angle))-vp[0]*np.sin(angle)],[vp[0]*vp[2]*(1-np.cos(angle))-vp[1]*np.sin(angle),vp[1]*vp[2]*(1-np.cos(angle))+vp[0]*np.sin(angle),np.cos(angle)+vp[2]*vp[2]*(1-np.cos(angle))]]

    ##This is the formula for the rotation matrix to change axis

    return rMat

def rotation(coord, rMat): #Apply a rotation to coordinates
    newCoord = []

    for i in range(len(coord)):
        newX = coord[i][0]*rMat[0][0] + coord[i][1]*rMat[1][0] + coord[i][2]*rMat[2][0]
        newY = coord[i][0]*rMat[0][1] + coord[i][1]*rMat[1][1] + coord[i][2]*rMat[2][1]
        newZ = coord[i][0]*rMat[0][2] + coord[i][1]*rMat[1][2] + coord[i][2]*rMat[2][2]
        newCoord.append([newX,newY,newZ])

    return newCoord

def cut_bath(rBath, coords): #Select which atoms are in the bath with distance constraint (sphere)
    bath = []
    start = datetime.datetime.now()
    for i in range(len(coords)):
        now = datetime.datetime.now()
        printProgressBar(start,now,i+1,len(coords),prefix='Cutting the bath',length=50)
        if distance([0,0,0],[coords[i][0],coords[i][1],coords[i][2]]) <= rBath:
                bath.append([coords[i][0], coords[i][1], coords[i][2], coords[i][3], 'C'])
        else:
            printProgressBar(start,now,1,1,prefix='Cutting the bath',length=50)
            return bath 

def set_pp(rPP,coords, notIn): #Select which atoms are in the first shell of pseudopotential 

    pp = []
    start = datetime.datetime.now()
    for i in range(len(coords)):
        now = datetime.datetime.now()
        printProgressBar(start,now,i+1,len(coords),prefix='Finding the first shell',length=50)
        for j in range(len(coords)):
            if coords[i][4] != 'O' or dist_zero(coords[j]) > 2*(max([atoms[k] for k in range(3,len(atoms),4)])*sum(npattern)+rPP):
                break
            if coords[i][4] == 'O':
                if coords[j][4] == 'C' and coords[j][3] not in notIn:
                    if distance([coords[i][0],coords[i][1],coords[i][2]],[coords[j][0],coords[j][1],coords[j][2]]) <= rPP:
                        coords[j][4] == 'Cl'
                        pp.append(j)
    for i in pp:
        coords[i][4] = 'Cl'
    return coords

def find_frag(pattern, n, coords):                                                           #We mark the atoms in the bath corresponding to
                                                                                             #the fragment according to user input 
    inFrag = []
    start = datetime.datetime.now()
    inPattern = []
    for k in range(n):
        closest = [100,100,100]
        now = datetime.datetime.now()
        printProgressBar(start,now,k+1,n,prefix='Finding the fragment',length=50)
        for j in coords:
            if j[3] == pattern[1]:
                if distance(j,[0,0,0]) < distance([0,0,0],closest) and [j[0],j[1],j[2],distance(j,j), coords.index(j)] not in inFrag:
                    notin = 0
                    for noin in noinfrag:
                        if distance(j,[noin[1],noin[2],noin[3]]) < da:
                            notin = 1
                    if notin == 0:
                        closest = [j[0],j[1],j[2],distance(j,j), coords.index(j)]
        for i in range(len(pattern)//2):                                                      
            inPattern = [closest]
            for j in range(int(pattern[2*i])):
                inPattern.append([100,100,100,distance([100,100,100],closest)])
            for j in coords:
                inPattern = sorted(inPattern,key=operator.itemgetter(3))
                if j[3] == pattern[2*i+1]:
                    if distance(j,closest) <= inPattern[-1][3]:
                        inPattern[-1] = [j[0],j[1],j[2], distance(j,closest), coords.index(j)]
            for j in inPattern:
                inFrag.append(j)
    for j in inFrag:
        coords[j[4]][4] = 'O'

    return coords


def symmetry(coord,atoms,charges, operations): #Find symmetry elements in the coordinates
    newCoord = []
    total = len(coord)
    start = datetime.datetime.now()
    progress = 0
    coord = sorted(coord,key=dist_zero)

    while coord != []:
        toDel = []
        newCoord.append(coord[0])         #Add the atom to a new list
        name = atoms[0]
        a = newCoord[-1][:]               #label a = E
        b = [-a[0],a[1],a[2]]             #label b = yOz mirror plan
        c = [-a[0],-a[1],a[2]]            #label c =  C2 rotation around z
        d = [a[0],-a[1],a[2]]             #label d = xOz mirror plan
        e = [a[0],a[1],-a[2]]             #label e = xOy mirror plan
        f = [-a[0],a[1],-a[2]]            #label f = C2 rotation around y
        g = [a[0],-a[1],-a[2]]            #label g = C2 rotation around x
        h = [-a[0],-a[1],-a[2]]           #label h = i
        newCoord[-1].append(name+"a")
        newCoord[-1].append(charges[0])
        progress += 1
        del atoms[0]                      #Delete from old list
        del coord[0]
        del charges[0]
        
        now = datetime.datetime.now()
        printProgressBar(start,now,progress,total,prefix='Treating Symmetry',length=50,decimals=3)

        for t in coord:
            index = coord.index(t)
            if np.absolute(dist_zero(a) - dist_zero(t)) > da:
                break
            if name == atoms[index]:  #Check if it is the same atom
                if distance(t,a) == 0:
                    print("ERROR : Twice the same atom")
                if 'xOz' in operations:
                    if distance(t,d) <= da:
                        newCoord.append(d)
                        newCoord[-1].append(name+'d')
                        newCoord[-1].append(charges[index])
                        if index not in toDel:
                            progress += 1
                            toDel.append(index)
                    elif da < distance(t,d) and distance(t,d) < DA:
                        print("Error : This atom should not be there",distance(t,d),t,d,a,charges[index])
                        print("Are you sure about the xOz symmetry operation ?")
                        break
                if 'yOz' in operations:
                    if distance(t,b) <= da:
                        newCoord.append(b)
                        newCoord[-1].append(name+'b')
                        newCoord[-1].append(charges[index])
                        if index not in toDel:
                            progress += 1
                            toDel.append(index)
                    elif da < distance(t,b) and distance(t,b) < DA:
                        print("Error : This atom should not be there",distance(t,d),t,d,a,charges[index])
                        print("Are you sure about the yOz symmetry operation ?")
                        break
                if 'C2z' in operations:
                    if distance(t,c) <= da:
                        newCoord.append(c)
                        newCoord[-1].append(name+'c')
                        newCoord[-1].append(charges[index])
                        if index not in toDel:
                            progress += 1
                            toDel.append(index)
                    elif da < distance(t,c) and distance(t,c) < DA:
                        print("Error : This atom should not be there",distance(t,d),t,d,a,charges[index])
                        print("Are you sure about the C2z symmetry operation ?")
                        break
                if 'xOy' in operations:
                    if distance(t,e) <= da:
                        newCoord.append(e)
                        newCoord[-1].append(name+'e')
                        newCoord[-1].append(charges[index])
                        if index not in toDel:
                            progress += 1
                            toDel.append(index)
                    elif da < distance(t,e) and distance(t,e) < DA:
                        print("Error : This atom should not be there",distance(t,d),t,d,a,charges[index])
                        print("Are you sure about the xOy symmetry operation ?")
                        break
                if 'C2y' in operations:
                    if distance(t,f) <= da:
                        newCoord.append(f)
                        newCoord[-1].append(name+'f')
                        newCoord[-1].append(charges[index])
                        if index not in toDel:
                            progress += 1
                            toDel.append(index)
                    elif da < distance(t,f) and distance(t,f) < DA:
                        print("Error : This atom should not be there",distance(t,d),t,d,a,charges[index])
                        print("Are you sure about the C2y symmetry operation ?")
                        break
                if 'C2x' in operations:
                    if distance(t,g) <= da:
                        newCoord.append(g)
                        newCoord[-1].append(name+'g')
                        newCoord[-1].append(charges[index])
                        if index not in toDel:
                            progress += 1
                            toDel.append(index)
                    elif da < distance(t,g) and distance(t,g) < DA:
                        print("Error : This atom should not be there",distance(t,d),t,d,a,charges[index])
                        print("Are you sure about the C2x symmetry operation ?")
                        break
                if 'i' in operations:
                    if distance(t,h) <= da:
                        newCoord.append(h)
                        newCoord[-1].append(name+'h')
                        newCoord[-1].append(charges[index])
                        if index not in toDel:
                            progress += 1
                            toDel.append(index)
                    elif da < distance(t,h) and distance(t,h) < DA:
                        print("Error : This atom should not be there",distance(t,d),t,d,a,charges[index])
                        print("Are you sure about the i symmetry operation ?")
                        break
            now = datetime.datetime.now()
            printProgressBar(start,now,progress,total,prefix='Treating Symmetry',length=50,decimals=3)


        for m in range(len(toDel)):    #We delete the atoms seen in the simmetry
            del coord[toDel[m]-m]
            del atoms[toDel[m]-m]
            del charges[toDel[m]-m]
    return newCoord

def write_input(fragCoord,ppCoord,bathCoord,fileName, sym):
    g = open(fileName,'w')
    if sym != 'x':
        g.write('FRAGMENT\n')
        g.write('LABEL      X           Y            Z           CHARGE\n')
        for i in range(len(fragCoord)):
            if fragCoord[i][3][-1] == 'a':
                g.write('%8s      % 7.8f      % 7.8f       % 7.8f      % 8.5f\n' %(fragCoord[i][3][:-1]+str(i+1),fragCoord[i][0],fragCoord[i][1],fragCoord[i][2],fragCoord[i][4]))
        g.write("\n\n")
        g.write('PSEUDO \n')
        g.write('LABEL      X           Y            Z           CHARGE\n')
        for i in range(len(ppCoord)):
            if ppCoord[i][3][-1] == 'a':
                g.write('%8s      % 7.8f      % 7.8f       % 7.8f      % 8.5f\n' %(ppCoord[i][3][:-1]+str(i+1),ppCoord[i][0],ppCoord[i][1],ppCoord[i][2],ppCoord[i][4]))
        g.write("\n\n")
        g.write('CHARGES\n')
        g.write('LABEL      X           Y            Z           CHARGE\n')
        for i in range(len(bathCoord)):
           if bathCoord[i][3][-1] == 'a':
               g.write('%8s       % 7.8f      % 7.8f       % 7.8f       % 8.5f\n'%(bathCoord[i][3][:-1]+str(i+1),bathCoord[i][0],bathCoord[i][1],bathCoord[i][2],bathCoord[i][4]))
    if sym == 'x':
        g.write('FRAGMENT\n')
        g.write('LABEL      X           Y            Z           CHARGE\n')
        for i in range(len(fragCoord)):
            g.write('%8s      % 7.8f      % 7.8f       % 7.8f       % 8.5f\n' %(fragCoord[i][3]+str(i+1),fragCoord[i][0],fragCoord[i][1],fragCoord[i][2], fragCoord[i][4]))
        g.write('\n\n')
        g.write('PSEUDO\n')
        g.write('LABEL      X           Y            Z           CHARGE\n')
        for i in range(len(ppCoord)):
            g.write('%8s      % 7.8f      % 7.8f       % 7.8f       % 8.5f\n' %(ppCoord[i][3]+str(i+1),ppCoord[i][0],ppCoord[i][1],ppCoord[i][2],ppCoord[i][4]))
        g.write('\n\n')
        g.write('CHARGES\n')
        g.write('LABEL      X           Y            Z           CHARGE\n')
        for i in range(len(bathCoord)):
            g.write('%8s      % 7.8f      % 7.8f       % 7.8f       % 8.5f\n' %(bathCoord[i][3]+str(i+1),bathCoord[i][0],bathCoord[i][1],bathCoord[i][2], bathCoord[i][4]))
    g.close()

def translation(vec, coords): 
    for i in range(len(coords)):
        coords[i][0] -= vec[0]
        coords[i][1] -= vec[1]
        coords[i][2] -= vec[2]
    return coords

def get_charge(charge,numbers,const):  #Return the square of the charge of the atoms 
    result = 0
    for i in const:
        result += i[0]*i[1]

    for i in range(len(numbers)):
        result += numbers[i]*charge[i]

    return result**2

def optimization(coords):
    numbers = []
    const = []
    charge = []
    nneighbour = []

    for atom in coords:
        ini = 0
        if atom[5] == 'full':
            for i in const:
                if i[0] == atoms[atoms.index(atom[3])+1]:
                    i[1] += 1
                    ini = 1
            if ini == 0:
                const.append([atoms[atoms.index(atom[3])+1],1])
        else:
            ini = 0
            for i in range(len(charge)):
                if nneighbour[i] == atom[5] and charge[i] == atoms[atoms.index(atom[3])+1]:
                    numbers[i] += 1
                    ini = 1
            if ini == 0:
                nneighbour.append(atom[5])
                numbers.append(1)
                charge.append(atoms[atoms.index(atom[3])+1])

    results = optimize.minimize(get_charge,charge,args=(numbers,const))  #Scipy built in method thad uses gradient descent to find the local minima of a given function, here it works with the square of the total charge (so that the minima will be at 0)

    print('  CHARGE             OPTIMIZED      CHANGE (%)')
    newCharges = results.x
    for i in range(len(charge)):
            print('% 7.5f             % 7.5f       % 3.2f\n'%(charge[i],newCharges[i],100-newCharges[i]/charge[i]*100))

    for atom in coords:
        if atom[5] == 'full':
            atom[5] = atoms[atoms.index(atom[3])+1]
        else:
            for i in range(len(charge)):
                if atom[5] == nneighbour[i] and charge[i] == atoms[atoms.index(atom[3])+1]:
                    atom[5] = newCharges[i]
    return coords

def count_neighbours(coords):
    start = datetime.datetime.now()

    coords = sorted(coords,key=get_xyz)

    count = 0

    for i in coords:
        i.append(0)


    for i in range(len(coords)-1):
        for j in range(i+1,len(coords)):
            if np.absolute(coords[j][0] - coords[i][0]) > max([atoms[i] for i in range(3,len(atoms),4)]):
                count += len(coords)-j
                printProgressBar(start,now,count,(len(coords)*(len(coords)-1))/2,prefix='Counting neighbours',length=50)
                break
            elif np.absolute(coords[j][1] - coords[j][1]) > max([atoms[i] for i in range(3,len(atoms),4)]):
                count += len(coords)-j
                printProgressBar(start,now,count,(len(coords)*(len(coords)-1))/2,prefix='Counting neighbours',length=50)
                break
            elif np.absolute(coords[j][2] - coords[j][2]) > max([atoms[i] for i in range(3,len(atoms),4)]):
                count += len(coords)-j
                printProgressBar(start,now,count,(len(coords)*(len(coords)-1))/2,prefix='Counting neighbours',length=50)
                break
            elif coords[i][5] == atoms[atoms.index(coords[i][3])+2]:
                count += len(coords)-j
                printProgressBar(start,now,count,(len(coords)*(len(coords)-1))/2,prefix='Counting neighbours',length=50)
                break
            count += 1
            now = datetime.datetime.now()
            printProgressBar(start,now,count,(len(coords)*(len(coords)-1))/2,prefix='Counting neighbours',length=50)
            if distance(coords[i],coords[j]) < atoms[atoms.index(coords[i][3])+3] and coords[i][3] != coords[j][3]:
                coords[i][5] += 1 
                coords[j][5] += 1


    for i in coords:
        if i[5] == atoms[atoms.index(i[3])+2]:
            i[5] == 'full'

    return coords

def get_xyz(a):
    return (a[0],a[1],a[2])

def evjen(coords):
    for i in range(len(coords)):
        if coords[i][5] == 'full':
            coords[i][5] = atoms[atoms.index(coords[i][3])+1]
        else:
            coords[i][5] = atoms[atoms.index(coords[i][3])+1] * coords[i][5]/atoms[atoms.index(coords[i][3])+2]

    return coords

def dist_zero(atom):
    return np.sqrt(atom[0]**2+atom[1]**2+atom[2]**2)

def main():

    print("Input file is : %s"%(sys.argv[1]))

    read_input(sys.argv[1])


    nA = int(np.floor(2*rBath/a)+2)                                    #We chose the number of time we need to replicate
    nB = int(np.floor(2*rBath/(b*np.sin(np.radians(gamma))))+2)         #to be able to cut the bath 
    nC = int(np.floor(2*rBath/(c*np.sin(np.radians(beta))))+2)

    coords = big_cell(nA,nB,nC)



    coords = translation([nA*a/2,nB*b/2,nC*c/2],coords)                    #Putting the origin at the center of the cell
    

        
    labels = [i[3] for i in coords]

    if trsl != 'x':
        translation(trsl,coords)

    centers = []                                                                                            #input, and translating the coordinates
    for i in range(len(center)):
        centers.append([100,100,100])

    for i in centers:
        i.append(distance([0,0,0],i))

    for i in coords:
        centers = sorted(centers,key=operator.itemgetter(3))
        if i[3] in center:
            if distance(i,[0,0,0]) <= centers[-1][3]:
                centers[-1] = [i[0],i[1],i[2],distance(i,[0,0,0])]
    newOgn = np.mean(np.array(centers),axis=0)
    newOgn = [newOgn[0], newOgn[1], newOgn[2]]

    coords = translation(newOgn,coords)

    axis = ['x','y','z']

    for k in axis:
        if k == 'x':          #Loop to find the 3 axis 
            nAxis = xAxis
        if k == 'y':
            nAxis = yAxis
        if k == 'z':
            nAxis = zAxis
        
        nAxiss = []                                                                                               #according to user input, and      
                                                                                                                  #rotating the coordinates
        if nAxis[0] == 'x':
            continue 
        for i in range(len(nAxis)):
            nAxiss.append([100,100,100])

        for i in nAxiss:
            i.append(distance([0,0,0],i))

        for i in coords:
            nAxiss = sorted(nAxiss,key=operator.itemgetter(3))
            if i[3] in nAxis and dist_zero(i) > da:
                if distance(i,[0,0,0]) <= nAxiss[-1][3]:
                    nAxiss[-1] = [i[0],i[1],i[2],distance(i,[0,0,0])]
        newN = np.mean(np.array(nAxiss),axis=0)
        newN = [newN[0],newN[1],newN[2]]

        if k == 'x':
            oldN = [1,0,0]
        elif k == 'y':
            oldN = [0,1,0]
        elif k == 'z':
            oldN = [0,0,1]

        rMat = rot_matrix(oldN,newN)
        coords = rotation(coords,rMat)
        coords = [[coords[i][0],coords[i][1],coords[i][2],labels[i]] for i in range(len(coords))]
    
   
   #We now have one big cell oriented and centered as we want

    coords = sorted(coords,key=dist_zero)
    coords = cut_bath(rBath,coords)
    for i in range(len(pattern)):
        coords = find_frag(pattern[i], npattern[i],coords)
    coords = sorted(coords,key=operator.itemgetter(4))
    coords = set_pp(rPP,coords,notIn)
    coords = sorted(coords,key=dist_zero,reverse=True)

    
    if evj == 1:
        coords = count_neighbours(coords)
        coords = evjen(coords)
    elif opti == 1:
        coords = count_neighbours(coords)
        coords = optimization(coords)
    else:
        for i in coords:
            i.append(atoms[atoms.index(i[3])+1])



   # coords = sorted(coords,key=operator.itemgetter(5))
   # frag = sorted([i for i in coords if i[4] == 'O'],key=operator.itemgetter(3))
   # pp = sorted([i for i in coords if i[4] == 'Cl'],key=operator.itemgetter(3))
   # bath = sorted([i for i in coords if i[4] == 'C'],key=operator.itemgetter(3))

    frag = sorted([i for i in coords if i[4] == 'O'],key=dist_zero)
    pp = sorted([i for i in coords if i[4] == 'Cl'],key=dist_zero)
    bath = sorted([i for i in coords if i[4] == 'C'],key=dist_zero)



    if seefrag == 1:
        g = open('tmp.xyz','w')
        g.write('%i \n \n'%len(frag))
        for j in frag:
            g.write('%s   % 6.5f    % 6.5f    % 6.5f \n'%(j[3],j[0],j[1],j[2]))
        g.close()
        os.system('avogadro tmp.xyz')
        os.system('rm tmp.xyz')
    ch = 0
    for i in range(len(coords)):
        ch += coords[i][5]
    print("Total charge : % 8.5f"%ch)


    if sym != 'x':
        if norep == 0:
            rep = 0
            prog = 0
            start = datetime.datetime.now()
            for i in range(len(coords)-1):
                for j in range(i+1,len(coords)):
                    prog += 1
                    now = datetime.datetime.now()
                    printProgressBar(start,now,prog,((len(coords)-1)*len(coords))/2,prefix='Calculating nuclear repulsion',length=50)
                    rep += (coords[i][5]*coords[j][5])/distance(coords[i],coords[j])
            print("Nuclear repulsion before symmetry : %f"%rep)

        frag = symmetry([[i[0],i[1],i[2]] for i in frag],[i[3] for i in frag], [i[5] for i in frag], sym)
        pp = symmetry([[i[0],i[1],i[2]] for i in pp],[i[3] for i in pp], [i[5] for i in pp], sym)
        bath = symmetry([[i[0],i[1],i[2]] for i in bath],[i[3] for i in bath], [i[5] for i in bath], sym)

        coords = frag+pp+bath
        
        if norep == 0:
            start = datetime.datetime.now()
            rep = 0
            prog = 0
            for i in range(len(coords)-1):
                for j in range(i+1,len(coords)):
                    prog += 1
                    now = datetime.datetime.now()
                    printProgressBar(start,now,prog,((len(coords)-1)*len(coords))/2,prefix='Calculating nuclear repulsion',length=50)
                    rep += (coords[i][4]*coords[j][4])/distance(coords[i],coords[j])
            print("Nuclear repulsion after symmetry : %f"%rep)
    else:
        frag = [[i[0],i[1],i[2],i[3],i[5]] for i in frag]
        pp = [[i[0],i[1],i[2],i[3],i[5]] for i in pp]
        bath = [[i[0],i[1],i[2],i[3],i[5]] for i in bath]

    frag = sorted(frag,key=operator.itemgetter(3))
    pp = sorted(pp,key=operator.itemgetter(3))
    bath = sorted(bath,key=operator.itemgetter(3))

    write_input(frag,pp,bath,output_file,sym)

    if visu != 0:
        g = open('tmp.xyz','w')
        if visu == 1:
            g.write('%i \n\n'%(len(frag)+len(pp)+len(bath)))
            for i in frag:
                g.write('O    % 7.8f   % 7.6f   % 7.8f  \n'%(i[0],i[1],i[2]))
            for i in pp:
                g.write('Cl    % 7.8f   % 7.8f   % 7.8f  \n'%(i[0],i[1],i[2]))
            for i in bath:
                g.write('C    % 7.8f   % 7.8f   % 7.8f  \n'%(i[0],i[1],i[2]))
        elif visu == 2:
            coords = frag+pp+bath
            if sym != 'x':
                for i in coords:
                    l = ['a','b','c','d','e','f','g','h']
                    if i[3][-1] in l:
                        i[3] = i[3][:-1]
            g.write('%i \n \n'%len(coords))
            for i in coords:
                g.write('%s    % 7.8f   % 7.8f   % 7.8f  \n'%(i[3],i[0],i[1],i[2]))
        g.close()
        os.system('avogadro tmp.xyz')
        os.system('rm tmp.xyz')



main()
