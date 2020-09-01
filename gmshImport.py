import numpy

class gmshImport:
    
    def __init__(self,filename):
        self.filename = filename
        
        self.elems = []
        self.nodes = []
        
        self.v1 = []
        self.v2 = []
        self.v3 = []
        
    def quadraticTriangleGmshReader(self):
        
        f = open(self.filename)
        f.readline() # '$MeshFormat\n'
        f.readline() # '2.2 0 8\n'
        f.readline() # '$EndMeshFormat\n'
        f.readline() # '$Nodes\n'
        n_nodes = int(f.readline()) # '8\n'
        nodes_ = numpy.fromfile(f,count=n_nodes*4, sep=" ").reshape((n_nodes,4))
    
    
        f.readline() # '$EndNodes\n'
        f.readline() # '$Elements\n'
        n_elems = int(f.readline()) # '2\n'
        elems_ = numpy.fromfile(f,dtype=int,count=n_elems*11, sep=" ").reshape((n_elems,11))
    
        nodes = []
        for i in range(0, n_nodes):
            nodes.append([nodes_[i,1],
                          nodes_[i,2],
                          nodes_[i,3]])
        
        elems = []
        for i in range(0, n_elems):
            #elems.append([6,
            elems.append([elems_[i,5]-1, 
                          elems_[i,6]-1, 
                          elems_[i,7]-1,
                          elems_[i,8]-1,
                          elems_[i,9]-1,
                          elems_[i,10]-1])
        
        self.elems = numpy.array(elems)
        self.nodes = numpy.array(nodes)

    
    def createLocalCoordinates(self):
        vertices = self.nodes
        faces    = self.elems
        faces = faces[:,0:3]
        
        tris = vertices[faces]
            
        n = numpy.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
        
        gmshImport.normalize_v3(n)
        
        v20 = numpy.zeros( n.shape, dtype=n.dtype )
        v20[:,1] = 1
        
        dot23 = numpy.multiply(v20,n)
        dot23 = dot23.sum(axis=1)
        dot23 = numpy.reshape(dot23,(dot23.shape[0],1))
        
        ones_ = numpy.ones([1,3])
        
        dot23 = numpy.matmul(dot23,ones_)
        
        v2 = v20 - numpy.multiply(dot23,n)
        
        gmshImport.normalize_v3(v2)
        
        v1x = numpy.multiply(v2[:,1],n[:,2]) - numpy.multiply(v2[:,2],n[:,1])
        v1y = numpy.multiply(v2[:,2],n[:,0]) - numpy.multiply(v2[:,0],n[:,2])
        v1z = numpy.multiply(v2[:,0],n[:,1]) - numpy.multiply(v2[:,1],n[:,0])

        v1x = numpy.reshape(v1x,(v1x.shape[0],1))
        v1y = numpy.reshape(v1y,(v1y.shape[0],1))
        v1z = numpy.reshape(v1z,(v1z.shape[0],1))

        v1 = numpy.concatenate((v1x, v1y, v1z), axis=1)
    
        self.v1 = v1
        self.v2 = v2
        self.v3 = n
        
    def normalize_v3(arr):
        
        lens = numpy.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
        arr[:,0] /= lens
        arr[:,1] /= lens
        arr[:,2] /= lens                
        return arr