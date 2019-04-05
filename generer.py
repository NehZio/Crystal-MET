import numpy as np
import os
import operator


###### CONSTANTS #######

a = 4.030
b = 4.030
c = 13.333

alpha = 90
beta = 90
gamma = 90

dr = 0.001
da = 0.01
DA = 1.5


chO = -2.00
chIr = 4.00
chSr = 2.00
chBa = 2.00
dIrO = 2.50
dSrO = 3.10
#######################

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

def distance(u,v):
    return np.sqrt((u[0]-v[0])**2+(u[1]-v[1])**2+(u[2]-v[2])**2)

def vect_product(u,v):
    return [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]

def dot_product(u,v):
    return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2])

def normalize(u):
    norm = np.sqrt(u[0]**2+u[1]**2+u[2]**2)
    return [u[0]/norm, u[1]/norm, u[2]/norm]

def rot_matrix(oldAxis,newAxis):
    newAxis = normalize(newAxis)
    vp = normalize(vect_product(oldAxis,newAxis))
    angle = np.arccos(dot_product(oldAxis,newAxis))

    rMat = [[np.cos(angle)+vp[0]*vp[0]*(1-np.cos(angle)),vp[0]*vp[1]*(1-np.cos(angle))-vp[2]*np.sin(angle),vp[0]*vp[2]*(1-np.cos(angle))+vp[1]*np.sin(angle)],[vp[0]*vp[1]*(1-np.cos(angle))+vp[2]*np.sin(angle),vp[1]*vp[1]*(1-np.cos(angle))+np.cos(angle),vp[1]*vp[2]*(1-np.cos(angle))-vp[0]*np.sin(angle)],[vp[0]*vp[2]*(1-np.cos(angle))-vp[1]*np.sin(angle),vp[1]*vp[2]*(1-np.cos(angle))+vp[0]*np.sin(angle),np.cos(angle)+vp[2]*vp[2]*(1-np.cos(angle))]]

    ##This is the formula for the rotation matrix to change axis

    return rMat

def rotation(coord, rMat):
    newCoord = []

    for i in range(len(coord)):
        newX = coord[i][0]*rMat[0][0] + coord[i][1]*rMat[1][0] + coord[i][2]*rMat[2][0]
        newY = coord[i][0]*rMat[0][1] + coord[i][1]*rMat[1][1] + coord[i][2]*rMat[2][1]
        newZ = coord[i][0]*rMat[0][2] + coord[i][1]*rMat[1][2] + coord[i][2]*rMat[2][2]
        newCoord.append([newX,newY,newZ])

    return newCoord

def cut_bath(rBath, coords):
    bath = []
    for i in range(len(coords)):
        if distance([0,0,0],[coords[i][0],coords[i][1],coords[i][2]]) <= rBath:
                bath.append([coords[i][0], coords[i][1], coords[i][2], coords[i][3], 'C'])

    return bath 

def set_pp(rPP,coords, notIn):

    pp = []

    for i in range(len(coords)):
        for j in range(len(coords)):
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
    for k in range(n):
        closest = [100,100,100]
        for j in coords:
            if j[3] == pattern[1]:
                if distance(j,[0,0,0]) < distance([0,0,0],closest) and [j[0],j[1],j[2],distance(j,j), coords.index(j)] not in inFrag:
                    closest = [j[0],j[1],j[2],distance(j,j), coords.index(j)]
        for i in range(1,len(pattern)/2):                                                      
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

def eivgen(coords):
    charges = []
    for i in coords:
        nIr = 0
        nSr = 0
        nO = 0
        for j in coords:
            if i[3] == 'Ir' and j[3] == 'O' and distance(i,j) < dIrO:
                nO += 1
            elif i[3] == 'Sr' and j[3] == 'O' and distance(i,j) < dSrO:
                nO += 1
        if i[3] == 'Ir':
            charges.append(chIr*(nO/6.00))
        if i[3] == 'Sr':
            charges.append(chSr*(nO/9.00))
        
        typeO = 0

        if i[3] == 'O':
            for j in coords:
                if j[3] == 'Ir' and distance(i,j) < dIrO:
                    nIr += 1
                    if distance(i,j) < 2.00:
                        typeO = 1
                elif j[3] == 'Sr' and distance(i,j) < dSrO:
                    nSr += 1
            if typeO == 0:
                charges.append(chO*((chSr*nSr+chIr*nIr)/(5*chSr+chIr)))
            else:
                charges.append(chO*((chSr*nSr+chIr*nIr)/(4*chSr+2*chIr)))


        

    coords = [[coords[i][0],coords[i][1],coords[i][2],coords[i][3],coords[i][4],charges[i]] for i in range(len(coords))]

    return coords

def symmetry(coord,atoms,charges, operations):
    newCoord = []

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
        del atoms[0]                      #Delete from old list
        del coord[0]
        del charges[0]

        for t in coord:
            index = coord.index(t)
            if name == atoms[index]:  #Check if it is the same atom
                if distance(t,a) == 0:
                    print("ERROR : Twice the same atom")
                if 'xOz' in operations:
                    if distance(t,d) <= da:
                        newCoord.append(d)
                        newCoord[-1].append(name+'d')
                        newCoord[-1].append(charges[index])
                        toDel.append(index)
                    elif da < distance(t,d) and distance(t,d) < DA:
                        print "Error : This atom should not be there",distance(t,d),t,d,a,charges[index]
                        print "Are you sure about the xOz symmetry operation ?"
                        break
                if 'yOz' in operations:
                    if distance(t,b) <= da:
                        newCoord.append(b)
                        newCoord[-1].append(name+'b')
                        newCoord[-1].append(charges[index])
                        toDel.append(index)
                    elif da < distance(t,b) and distance(t,b) < DA:
                        print "Error : This atom should not be there",distance(t,b),t,b,a,charges[index]
                        print "Are you sure about the yOz symmetry operation ?"
                        break
                if 'C2z' in operations:
                    if distance(t,c) <= da:
                        newCoord.append(c)
                        newCoord[-1].append(name+'c')
                        newCoord[-1].append(charges[index])
                        toDel.append(index)
                    elif da < distance(t,c) and distance(t,c) < DA:
                        print "Error : This atom should not be there",distance(t,c),t,c,a,charges[index]
                        print "Are you sure about the C2z axis ?"
                        break
                if 'xOy' in operations:
                    if distance(t,e) <= da:
                        newCoord.append(e)
                        newCoord[-1].append(name+'e')
                        newCoord[-1].append(charges[index])
                        toDel.append(index)
                    elif da < distance(t,e) and distance(t,e) < DA:
                        print "Error : This atom should not be there",distance(t,e),t,e,a,charges[index]
                        print "Are you sure about the xOy operation ?"
                        break
                if 'C2y' in operations:
                    if distance(t,f) <= da:
                        newCoord.append(f)
                        newCoord[-1].append(name+'f')
                        newCoord[-1].append(charges[index])
                        toDel.append(index)
                    elif da < distance(t,f) and distance(t,f) < DA:
                        print "Error : This atom should not be there",distance(t,f),t,f,a,charges[index]
                        print "Are you sure about the C2y axis ?"
                        break
                if 'C2x' in operations:
                    if distance(t,g) <= da:
                        newCoord.append(g)
                        newCoord[-1].append(name+'g')
                        newCoord[-1].append(charges[index])
                        toDel.append(index)
                    elif da < distance(t,g) and distance(t,g) < DA:
                        print "Error : This atom should not be there",distance(t,g),t,g,a,charges[index]
                        print "Are you sure about the C2x axis ?"
                        break
                if 'i' in operations:
                    if distance(t,h) <= da:
                        newCoord.append(h)
                        newCoord[-1].append(name+'h')
                        newCoord[-1].append(charges[index])
                        toDel.append(index)
                    elif da < distance(t,h) and distance(t,h) < DA:
                        print "Error : This atom should not be there",distance(t,h),t,h,a,charges[index]
                        print "Are you sure about the i operation ?"
                        break


        for m in range(len(toDel)):    #We delete the atoms seen in the simmetry
            del coord[toDel[m]-m]
            del atoms[toDel[m]-m]
            del charges[toDel[m]-m]

    return newCoord

def write_input(fragCoord,ppCoord,bathCoord,fileName, sym):
    g = open(fileName,'w')
    if sym == 'y':
        for i in range(len(fragCoord)):
            if fragCoord[i][3][-1] == 'a':
                g.write('%4s      % 5.3f      % 5.3f       % 5.3f     angstrom\n' %(fragCoord[i][3].replace('a','')+str(i+1),fragCoord[i][0],fragCoord[i][1],fragCoord[i][2]))
        g.write("\n\n")
        for i in range(len(ppCoord)):
            if ppCoord[i][3][-1] == 'a':
                g.write('%4s      % 5.3f      % 5.3f       % 5.3f     angstrom\n' %(ppCoord[i][3].replace('a','')+str(i+1),ppCoord[i][0],ppCoord[i][1],ppCoord[i][2]))
        g.write("\n\n")
        for i in range(len(bathCoord)):
           # if bathCoord[i][3][-1] == 'a':
           g.write('% 5.3f      % 5.3f     % 5.3f     % 5.3f    0.     0.    0.  \n'%(bathCoord[i][0],bathCoord[i][1],bathCoord[i][2],bathCoord[i][4]))
    if sym == 'n':
        for i in range(len(fragCoord)):
            g.write('%4s      % 5.3f      % 5.3f       % 5.3f     angstrom\n' %(fragCoord[i][3]+str(i+1),fragCoord[i][0],fragCoord[i][1],fragCoord[i][2]))
        g.write('\n\n')
        for i in range(len(ppCoord)):
            g.write('%4s      % 5.3f      % 5.3f       % 5.3f     angstrom\n' %(ppCoord[i][3]+str(i+1),ppCoord[i][0],ppCoord[i][1],ppCoord[i][2]))
        g.write('\n\n')
        for i in range(len(bathCoord)):
            g.write('% 5.3f      % 5.3f       % 5.3f    % 5.3f    0.    0.    0.  \n' %(bathCoord[i][0],bathCoord[i][1],bathCoord[i][2], bathCoord[i][4]))
    g.close()

def translation(vec, coords):
    for i in range(len(coords)):
        coords[i][0] -= vec[0]
        coords[i][1] -= vec[1]
        coords[i][2] -= vec[2]
    return coords

def main():

    fileName = raw_input("Name of the cif file : ")
    rBath = input("Chose the bath radius (in Angstrom) : ")
    rPP = input("Chose the pseudopotential radius (in Angstrom) : ")
    notIn = raw_input("What atoms should not be in the first shell ? ").split()

    nA = int(np.floor(2*rBath/a)+2)                                    #We chose the number of time we need to replicate
    nB = int(np.floor(2*rBath/(b*np.sin(np.radians(gamma))))+2)         #to be able to cut the bath 
    nC = int(np.floor(2*rBath/(c*np.sin(np.radians(beta))))+2)

    cmd = 'atomsk '+ fileName + ' ' +  '-duplicate ' + str(nA) + ' ' + str(nB) + ' ' + str(nC) + ' ' + fileName.replace('cif','xyz') + ' -v 0' #This is the command that calls the program that generates the big cell
    os.system(cmd)


    data = parse(fileName.replace('cif','xyz'))                               #Read the data from the xyz file

    coords = [[data[i][1],data[i][2],data[i][3]] for i in range(len(data))]
    labels = [data[i][0] for i in range(len(data))]

    coords = [[coords[i][0],coords[i][1],coords[i][2],labels[i]] for i in range(len(coords))]

    coords = translation([nA*a/2,nB*b/2,nC*c/2],coords)                    #Putting the origin at the center of the cell

        
    labels = [i[3] for i in coords]

    if raw_input("Do you want to make an initial translation ? (y/n) ") == "y":
        trsl = raw_input("Enter the coordinates of the translation ").split()
        trsl[0] = float(trsl[0])
        trsl[1] = float(trsl[1])
        trsl[2] = float(trsl[2])
        translation(trsl,coords)

    center = raw_input("Specify the center of the fragment (if between atoms, specify them all) ").split()  #Searching for the center according to user
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
        nAxis = raw_input("Where should the %s axis be headed ? (if between atoms, specify them all), x if none "%k).split() #Searching for the new orientation
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
            if i[3] in nAxis:
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
   #The rest of the code will cut what we want in this big cell

    pattern = raw_input("What pattern is in the fragment ? ").split()

    n = input("How many time does the pattern repeat in the fragment ? ")



    coords = cut_bath(rBath,coords)
    coords = find_frag(pattern, n ,coords)
    coords = set_pp(rPP,coords,notIn)


    coords = eivgen(coords)
    ch = 0
    for i in range(len(coords)):
        ch += coords[i][5]
    print("Total charge : % 8.5f"%ch)
    sym = 'x'

    os.system('rm '+fileName.replace('cif','xyz'))

    frag = sorted([i for i in coords if i[4] == 'O'],key=operator.itemgetter(3))
    pp = sorted([i for i in coords if i[4] == 'Cl'],key=operator.itemgetter(3))
    bath = sorted([i for i in coords if i[4] == 'C'],key=operator.itemgetter(3))

    if raw_input("Do you want to see the fragment ? (y/n) ") == 'y':
        g = open('tmp.xyz','w')
        g.write('%i \n \n'%len(frag))
        for j in frag:
            g.write('%s   % 6.2f    % 6.2f    % 6.2f \n'%(j[3],j[0],j[1],j[2]))
        g.close()
        os.system('avogadro tmp.xyz')
        os.system('rm tmp.xyz')

    while sym != 'y' and sym != 'n':
        sym = raw_input("Do you want to treat the symmetry ? (y/n) ")

    if sym == 'y':
        operation = raw_input("What are the symmetry operations you'd like to treat ? (C2(x,y,z), xOy xOz yOz, i) ").split()
        rep = 0
        for i in range(len(coords)-1):
            for j in range(i+1,len(coords)):
                rep += (coords[i][5]*coords[j][5])/distance(coords[i],coords[j])
        print("Nuclear repulsion before symmetry : %f"%rep)

        frag = symmetry([[i[0],i[1],i[2]] for i in frag],[i[3] for i in frag], [i[5] for i in frag], operation)
        pp = symmetry([[i[0],i[1],i[2]] for i in pp],[i[3] for i in pp], [i[5] for i in pp], operation)
        #bath = symmetry([[i[0],i[1],i[2]] for i in bath],[i[3] for i in bath], [i[5] for i in bath], operation)
        bath = [[i[0],i[1],i[2],i[3],i[5]] for i in bath]

        coords = frag+pp+bath
        
        rep = 0
        for i in range(len(coords)-1):
            for j in range(i+1,len(coords)):
                rep += (coords[i][4]*coords[j][4])/distance(coords[i],coords[j])

        print("Nuclear repulsion after symmetry : %f"%rep)

    if sym == 'n':
        frag = [[i[0],i[1],i[2],i[3],i[5]] for i in frag]
        pp = [[i[0],i[1],i[2],i[3],i[5]] for i in pp]
        bath = [[i[0],i[1],i[2],i[3],i[5]] for i in bath]
    fileName = raw_input("What name do you wish yo give to your file ? ")

    write_input(frag,pp,bath,fileName,sym)

    if raw_input("Do you want to visualize the bath ? (y/n) ") == 'y':
        g = open('tmp.xyz','w')
        if raw_input("Do you want it to be sorted with colors ? (y/n) ") == 'y':
            g.write('%i \n\n'%(len(frag)+len(pp)+len(bath)))
            for i in frag:
                g.write('O    % 5.3f   % 5.3f   % 5.3f  \n'%(i[0],i[1],i[2]))
            for i in pp:
                g.write('Cl    % 5.3f   % 5.3f   % 5.3f  \n'%(i[0],i[1],i[2]))
            for i in bath:
                g.write('C    % 5.3f   % 5.3f   % 5.3f  \n'%(i[0],i[1],i[2]))
        else:
            coords = frag+pp+bath
            if sym == 'y':
                for i in coords:
                    l = ['a','b','c','d','e','f','g','h']
                    if i[3][-1] in l:
                        i[3] = i[3][:-1]
            g.write('%i \n \n'%len(coords))
            for i in coords:
                g.write('%s    % 5.3f   % 5.3f   % 5.3f  \n'%(i[3],i[0],i[1],i[2]))
        g.close()
        os.system('avogadro tmp.xyz')
        os.system('rm tmp.xyz')


main()
