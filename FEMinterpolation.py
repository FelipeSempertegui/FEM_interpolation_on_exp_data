import numpy

import vtk

from vtk.util.numpy_support import numpy_to_vtk
    
class FEMinterpolation:
    
    def __init__(self, nodes, elems, e1, e2, PointCloud, data):
        
        self.vertices = nodes
        self.faces = elems
        self.v1 = e1
        self.v2 = e2
        
        self.nodes_exp = PointCloud
        self.data = data
        
    def interpolateData(self):
        
        faces = self.faces
        v1 = self.v1
        v2 = self.v2
        
        x = self.vertices[:,0:1]
        y = self.vertices[:,1:2]
        z = self.vertices[:,2:3]
        
        x1=x[faces[:,0]]
        y1=y[faces[:,0]]
        z1=z[faces[:,0]]
        
        x2=x[faces[:,1]]
        y2=y[faces[:,1]]
        z2=z[faces[:,1]]
        
        x3=x[faces[:,2]]
        y3=y[faces[:,2]]
        z3=z[faces[:,2]]
        
        xg= numpy.mean(x[faces[:,0:3]],axis=1)
        yg= numpy.mean(y[faces[:,0:3]],axis=1)
        zg= numpy.mean(z[faces[:,0:3]],axis=1)

        x_exp = self.nodes_exp[:,0:1]
        y_exp = self.nodes_exp[:,1:2]
        z_exp = self.nodes_exp[:,2:3]
        
        a = numpy.zeros((numpy.size(self.nodes_exp,axis=0),numpy.size(self.vertices,axis=0)))
        
        for k in range(0, len(faces)):
            alpha1 = (x1[k]-xg[k])*v1[k,0] + (y1[k]-yg[k])*v1[k,1] + (z1[k]-zg[k])*v1[k,2]
            beta1  = (x1[k]-xg[k])*v2[k,0] + (y1[k]-yg[k])*v2[k,1] + (z1[k]-zg[k])*v2[k,2]
                     
            alpha2 = (x2[k]-xg[k])*v1[k,0] + (y2[k]-yg[k])*v1[k,1] + (z2[k]-zg[k])*v1[k,2]
            beta2  = (x2[k]-xg[k])*v2[k,0] + (y2[k]-yg[k])*v2[k,1] + (z2[k]-zg[k])*v2[k,2]
    
            alpha3 = (x3[k]-xg[k])*v1[k,0] + (y3[k]-yg[k])*v1[k,1] + (z3[k]-zg[k])*v1[k,2]
            beta3  = (x3[k]-xg[k])*v2[k,0] + (y3[k]-yg[k])*v2[k,1] + (z3[k]-zg[k])*v2[k,2]
            
            alpha_exp_all = (x_exp-xg[k])*v1[k,0] + (y_exp-yg[k])*v1[k,1] + (z_exp-zg[k])*v1[k,2]
            beta_exp_all  = (x_exp-xg[k])*v2[k,0] + (y_exp-yg[k])*v2[k,1] + (z_exp-zg[k])*v2[k,2]    
            
            s = numpy.array([[1, alpha1, beta1],[1, alpha2, beta2],[1, alpha3, beta3]])
            det_s = abs(numpy.linalg.det(s))
            
            B1 = ( (alpha2*beta3 - alpha3*beta2) + (beta2-beta3)*alpha_exp_all + (alpha3-alpha2)*beta_exp_all )/det_s
            B2 = ( (alpha3*beta1 - alpha1*beta3) + (beta3-beta1)*alpha_exp_all + (alpha1-alpha3)*beta_exp_all )/det_s
            B3 = ( (alpha1*beta2 - alpha2*beta1) + (beta1-beta2)*alpha_exp_all + (alpha2-alpha1)*beta_exp_all )/det_s
            
            B1[B1 < 0] = 0
            B2[B2 < 0] = 0
            B3[B3 < 0] = 0
            
            B1[B1 > 1] = 0
            B2[B2 > 1] = 0
            B3[B3 > 1] = 0
            
            Bf_ = B1 + B2 + B3
            Bf_[Bf_ > 1.001] = 0
            Bf_[Bf_ < 0.999] = 0
            
            Bf = numpy.multiply(B1,B2,B3)
            Bf[Bf > 0] = 1
            
            test = numpy.multiply(Bf_,Bf)
            
            value_index = numpy.nonzero(test)[0]
            
            x_in_elem = x_exp[value_index[:]]
            y_in_elem = y_exp[value_index[:]]
            z_in_elem = z_exp[value_index[:]]
            
            alpha_exp = (x_in_elem-xg[k])*v1[k,0] + (y_in_elem-yg[k])*v1[k,1] + (z_in_elem-zg[k])*v1[k,2]
            beta_exp  = (x_in_elem-xg[k])*v2[k,0] + (y_in_elem-yg[k])*v2[k,1] + (z_in_elem-zg[k])*v2[k,2]    

            qsi1 = ( (alpha2*beta3 - alpha3*beta2) + (beta2-beta3)*alpha_exp + (alpha3-alpha2)*beta_exp )/det_s
            qsi2 = ( (alpha3*beta1 - alpha1*beta3) + (beta3-beta1)*alpha_exp + (alpha1-alpha3)*beta_exp )/det_s
            qsi3 = ( (alpha1*beta2 - alpha2*beta1) + (beta1-beta2)*alpha_exp + (alpha2-alpha1)*beta_exp )/det_s

            N1 = numpy.multiply(qsi1,2*qsi1-1)
            N2 = numpy.multiply(qsi2,2*qsi2-1)
            N3 = numpy.multiply(qsi3,2*qsi3-1)
            N4 = 4*numpy.multiply(qsi1,qsi2)
            N5 = 4*numpy.multiply(qsi2,qsi3)
            N6 = 4*numpy.multiply(qsi3,qsi1)
            
            Ns = numpy.concatenate((N1,N2,N3,N4,N5,N6),axis=1)
            
            for u in range(0, len(value_index)):
                a[value_index[u],faces[k,:]] = Ns[u,:]
            
        a_t = numpy.transpose(a)
        a_ = numpy.matmul(a_t,a)
        a_inv=numpy.linalg.inv(a_)
        lsam=numpy.matmul(a_inv,a_t)
        self.displMesh = numpy.matmul(lsam,self.data)
        

    def VTKwriter(self, dataName, fileName):
        points = vtk.vtkPoints()
        for i in range(0, len(self.vertices)):
            points.InsertPoint(i, self.vertices[i])
    
        ugrid = vtk.vtkUnstructuredGrid()
        ugrid.SetPoints(points)
        for i in range(0, len(self.faces)):
            ugrid.InsertNextCell(vtk.VTK_QUADRATIC_TRIANGLE, 6, self.faces[i])
    
        # Fields to be added to the mesh
        nodes_vtk = numpy_to_vtk(self.displMesh)
        nodes_vtk.SetName(dataName)
        ugrid.GetPointData().AddArray(nodes_vtk) 

        # write the vtk file
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(fileName)
        writer.SetInputData(ugrid)
        writer.Write()
        
