import numpy
from gmshImport import gmshImport
from FEMinterpolation import FEMinterpolation

# =============== IMPORT QUADRATIC MESH ===============


# import quadratic mesh
# file containing a mesh with quadratic triangles is loaded
mesh = gmshImport('quadraticMesh.msh') 
# elements and nodes are extracted from the file
mesh.quadraticTriangleGmshReader() 
# local coordinates are set on the triangles
mesh.createLocalCoordinates()


# import data to be interpolated - example file
mat = numpy.loadtxt('experimentalData.txt')
# dense point cloud, which describes the sample's shape
pointCloud = mat[:,0:3]
# "data" that was measured for every coordinate describing the point cloud
data = mat[:,3:4]


# =============== INTERPOLATE DATA & CREATE VTK FILE ===============


# Create new object containing both the customized mesh and the experimental data
finalMesh = FEMinterpolation(mesh.nodes, mesh.elems, mesh.v1, mesh.v2, pointCloud, data) 
# Interpolate the "data" located at the "point cloud" to the "nodes" of the mesh
finalMesh.interpolateData()
# Export the mesh with the "interpolated data" as a VTK file
finalMesh.VTKwriter('interpolatedData', 'meshWithExpData.vtk')