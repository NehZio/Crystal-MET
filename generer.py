import numpy as np
import os
import operator


###### CONSTANTS #######

a = 5.48463
b = 5.48463
c = 25.7977

dr = 0.001
da = 0.01
DA = 1.5


chO = -2.00
chIr = 4.00
chSr = 2.00
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

def set_pp(rPP,coords):

    pp = []

    for i in range(len(coords)):
        for j in range(len(coords)):
            if coords[i][4] == 'O':
                if coords[j][4] == 'C' and coords[j][3] != 'O':
                    if distance([coords[i][0],coords[i][1],coords[i][2]],[coords[j][0],coords[j][1],coords[j][2]]) <= rPP:
                        coords[j][4] == 'Cl'
                        pp.append(j)
    for i in pp:
        coords[i][4] = 'Cl'
    return coords

def find_frag(coords):
    for i in coords:
        if distance(i,[0,0,0]) < 1: 
            i[4] = 'O'
            for j in range(len(coords)):
                if distance(i,coords[j]) != 0 and distance(i,coords[j]) < 2.00:
                    coords[j][4] = 'O'
                    for k in range(len(coords)):
                        if distance(coords[j],coords[k]) != 0 and distance(coords[j],coords[k]) < 2.06:
                            coords[k][4] = 'O'


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
            if bathCoord[i][3][-1] == 'a':
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

    nA = int(np.floor(2*rBath/a)+2)                                    #We chose the number of time we need to replicate
    nB = int(np.floor(2*rBath/b)+2)                                    #to be able to cut the bath 
    nC = int(np.floor(2*rBath/c)+2)

    cmd = 'atomsk '+ fileName + ' ' +  '-duplicate ' + str(nA) + ' ' + str(nB) + ' ' + str(nC) + ' ' + fileName.replace('cif','xyz') + ' -v 0'
    os.system(cmd)


    data = parse(fileName.replace('cif','xyz'))                               #Read the data from the xyz file

    coords = [[data[i][1],data[i][2],data[i][3]] for i in range(len(data))]
    labels = [data[i][0] for i in range(len(data))]

    coords = [[coords[i][0],coords[i][1],coords[i][2],labels[i]] for i in range(len(coords))]

    coords = translation([nA*a/2,nB*b/2,nC*c/2],coords)                    #Putting the origin at the center of the cell

    ################ SPECIFIC FOR THE WANTED FRAG ###################
    
    labels = [i[3] for i in coords]

    nearestIridium = [100,100,100]
    for i in coords:
        if i[3] == 'Ir':
            if distance([0,0,0],[i[0],i[1],i[2]]) < distance([0,0,0],nearestIridium) and distance([0,0,0],[i[0],i[1],i[2]]) > 1:
                nearestIridium = i
    nearestIridium = [np.absolute(nearestIridium[0]), nearestIridium[1], nearestIridium[2]]

    rMat = rot_matrix([1,0,0],nearestIridium)
    coords = rotation(coords, rMat)
    print("Rotated to a new X axis")

    coords = translation([0.5*distance([0,0,0],nearestIridium),0,0],coords)

    print("Translation ok")


    nearestOxygen = [100,100,100]
    for i in range(len(labels)):
        if labels[i] == 'O':
            if distance([0,0,0],[coords[i][0],coords[i][1],coords[i][2]]) < distance([0,0,0],nearestOxygen) :
                    nearestOxygen = coords[i]

    rMat = rot_matrix([0,0,1],nearestOxygen)
    coords = rotation(coords,rMat)
    print("Rotated to a new Z axis")
    coords = [[coords[i][0],coords[i][1],coords[i][2],labels[i]] for i in range(len(coords))]

    #################################################################a
    
    
    coords = cut_bath(rBath,coords)
    coords = find_frag(coords)
    coords = set_pp(rPP,coords)


    coords = eivgen(coords)
    ch = 0
    for i in range(len(coords)):
        ch += coords[i][5]
    print("Total charge : % 8.5f"%ch)
    sym = 'x'

    os.system('rm '+fileName.replace('cif','xyz'))

    while sym != 'y' and sym != 'n':
        sym = raw_input("Do you want to treat the symmetry ? (y/n) ")

    
    frag = sorted([i for i in coords if i[4] == 'O'],key=operator.itemgetter(3))
    pp = sorted([i for i in coords if i[4] == 'Cl'],key=operator.itemgetter(3))
    bath = sorted([i for i in coords if i[4] == 'C'],key=operator.itemgetter(3))

    if sym == 'y':
        operation = raw_input("What are the symmetry operations you'd like to treat ? (C2(x,y,z), xOy xOz yOz, i) ").split()
        rep = 0
        for i in range(len(coords)-1):
            for j in range(i+1,len(coords)):
                rep += (coords[i][5]*coords[j][5])/distance(coords[i],coords[j])
        print("Nuclear repulsion before symmetry : %f"%rep)

        frag = symmetry([[i[0],i[1],i[2]] for i in frag],[i[3] for i in frag], [i[5] for i in frag], operation)
        pp = symmetry([[i[0],i[1],i[2]] for i in pp],[i[3] for i in pp], [i[5] for i in pp], operation)
        bath = symmetry([[i[0],i[1],i[2]] for i in bath],[i[3] for i in bath], [i[5] for i in bath], operation)

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
        coords = frag+pp+bath
        if sym == 'y':
            for i in coords:
                l = ['a','b','c','d','e','f','g','h']
                if i[3][-1] in l:
                    i[3] = i[3][:-1]
        g = open('tmp.xyz','w')
        g.write('%i \n \n'%len(coords))
        for i in coords:
            g.write('%s    % 5.3f   % 5.3f   % 5.3f  \n'%(i[3],i[0],i[1],i[2]))
        g.close()
        os.system('avogadro tmp.xyz')
        os.system('rm tmp.xyz')


main()
