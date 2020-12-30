'''
def mesh(self, division, nodeStartCount=1, elementStartCount = 1):
        """Subdivide nurbs into quadrangle finite elements.
        
        Arguments :
        division -- number of division (the same for each direction)
        nodeStartCount -- the number of the fist node to be generated
        elementStartCount -- the number of the first element to be generated
        
        Return nodes, connectivity :
        node -- array containing for each node : [ nodeNumber , nurbsNumber,[xiTilde, etaTilde],[x, y]]
        connectivity -- array containing for each element : [elNumber, node1, node2, ...]
        
        with:
        nodeNumber : number of the generated node
        nurbsNumber : number of the corresponding nurbs element
        xiTilde : xiTilde coordinate in the parent space element
        etaTilde : etaTilde coordinate in the parent element
        x : x coordinate in the physical space
        y : y coordinate in the physical space
        elNumber : number of the generated element
        node1, node2, ... : numbers of the node of the generated element
        """
        
        nodes = []
        connectivity = []
        nodeCount = nodeStartCount
        elementCount = elementStartCount
        
        
        if(dimension == 2):
            for j in range(0, division+1):
                for i in range(0, division+1):
                    parentCoords = [-1.+i*(2./division),-1.+j*(2./division)]
                    
                    coords = [0.,0.]
                    R = self.shape23D_only_R(ptGauss = parentCoords)
                    for k in range(0, self.numberOfNodes):
                        for l in range(0, self.dimension):
                            coords[l] = coords[l] + R[k]*self.nodesCoordinates[k][l]
                    thisNode = []
                    thisNode.append(nodeCount)
                    thisNode.append(self.number)
                    thisNode.append(parentCoords)
                    thisNode.append(coords)
                    nodes.append(thisNode)
                    
                    nodeCount = nodeCount+1
        elif(dimension == 3):
            for k in range(0, division+1):
                for j in range(0, division+1):
                    for i in range(0, division+1):
                        parentCoords = [-1.+i*(2./division), -1.+j*(2./division), -1.+k*(2./division)]
                        
                        coords = [0., 0., 0.]
                        R, dRdx, DetJac = self.shape23D(ptGauss = parentCoords)
                        for l in range(0, self.numberOfNodes):
                            for m in range(0, self.dimension):
                                coords[m] = coords[m] + R[l]*self.nodesCoordinates[l][m]
                        thisNode = []
                        thisNode.append(nodeCount)
                        thisNode.append(self.number)
                        thisNode.append(parentCoords)
                        thisNode.append(coords)
                        nodes.append(thisNode)
                        
                        nodeCount = nodeCount+1
                        
  
        counter = 0
        if(self.dimension == 2):
            for i in range(0, division):
                for j in range(0, division):
                    thisElement = []
                    thisElement.append(elementStartCount + counter)
                    thisElement.append(nodeStartCount+i*(division+1)+j)
                    thisElement.append(nodeStartCount+i*(division+1)+j+1)
                    thisElement.append(nodeStartCount+(i+1)*(division+1)+j+1)
                    thisElement.append(nodeStartCount+(i+1)*(division+1)+j)
                
                    connectivity.append(thisElement)
                    counter = counter + 1
            nz=1
            ny=1
            nx=2
        elif(self.dimension == 3):
            for i in range(0, nz):
                for j in range(0, ny):
                    for k in range(0, nx):
                        thisElement = []
                        thisElement.append(elementStartCount + counter)
                        thisElement.append(nodeStartCount+i*(nx+1)*(ny+1)+j*(nz+1)+k)
                        thisElement.append(nodeStartCount+i*(nx+1)*(ny+1)+j*(nz+1)+k+1)
                        thisElement.append(nodeStartCount+i*(nx+1)*(ny+1)+(j+1)*(nz+1)+k+1)
                        thisElement.append(nodeStartCount+i*(nx+1)*(ny+1)+(j+1)*(nz+1)+k)
                        thisElement.append(nodeStartCount+(i+1)*(nx+1)*(ny+1)+j*(nz+1)+k)
                        thisElement.append(nodeStartCount+(i+1)*(nx+1)*(ny+1)+j*(nz+1)+k+1)
                        thisElement.append(nodeStartCount+(i+1)*(nx+1)*(ny+1)+(j+1)*(nz+1)+k+1)
                        thisElement.append(nodeStartCount+(i+1)*(nx+1)*(ny+1)+(j+1)*(nz+1)+k)

                        connectivity.append(thisElement)
                        counter = counter+1

        return nodes, connectivity
'''
nodeStartCount=0
elementStartCount=0
counter=0
nz=1
ny=1
nx=3
connectivity=[]

for i in range(0, nz):
    for j in range(0, ny):
        for k in range(0, nx):
            thisElement = []
            thisElement.append(elementStartCount + counter)
            thisElement.append(nodeStartCount+i*(nx+1)*(ny+1)+j*(nz+1)+k)
            thisElement.append(nodeStartCount+i*(nx+1)*(ny+1)+j*(nz+1)+k+1)
            thisElement.append(nodeStartCount+i*(nx+1)*(ny+1)+(j+1)*(nz+1)+k+1)
            thisElement.append(nodeStartCount+i*(nx+1)*(ny+1)+(j+1)*(nz+1)+k)
            thisElement.append(nodeStartCount+(i+1)*(nx+1)*(ny+1)+j*(nz+1)+k)
            thisElement.append(nodeStartCount+(i+1)*(nx+1)*(ny+1)+j*(nz+1)+k+1)
            thisElement.append(nodeStartCount+(i+1)*(nx+1)*(ny+1)+(j+1)*(nz+1)+k+1)
            thisElement.append(nodeStartCount+(i+1)*(nx+1)*(ny+1)+(j+1)*(nz+1)+k)
            '''
            thisElement.append(nodeStartCount+i*(division+1)*(division+1)+j*(division+1)+k)
            thisElement.append(nodeStartCount+i*(division+1)*(division+1)+j*(division+1)+k+1)
            thisElement.append(nodeStartCount+i*(division+1)*(division+1)+(j+1)*(division+1)+k+1)
            thisElement.append(nodeStartCount+i*(division+1)*(division+1)+(j+1)*(division+1)+k)
            thisElement.append(nodeStartCount+(i+1)*(division+1)*(division+1)+j*(division+1)+k)
            thisElement.append(nodeStartCount+(i+1)*(division+1)*(division+1)+j*(division+1)+k+1)
            thisElement.append(nodeStartCount+(i+1)*(division+1)*(division+1)+(j+1)*(division+1)+k+1)
            thisElement.append(nodeStartCount+(i+1)*(division+1)*(division+1)+(j+1)*(division+1)+k)
            '''
            connectivity.append(thisElement)
            counter = counter+1                                                            

print(connectivity)
index=0
beam_coordinates=[]
'''
for k in range(0, nz+1):
    for j in range(0, ny+1):
        for i in range(0, nx+1):
            parentCoords.append([-1.+i*(2./nx), -1.+j*(2./ny), -1.+k*(2./nz),index])
            
            index+=1

print(parentCoords)
'''
length=1
height=1
width=1
for i in range(0, nz+1):
    for j in range(0, ny+1):
        for k in range(0, nx+1):
            beam_coordinates.append([(k)*(length/nx),j*(height/ny),i*(width/nz),1,index])
            index+=1

print(beam_coordinates)