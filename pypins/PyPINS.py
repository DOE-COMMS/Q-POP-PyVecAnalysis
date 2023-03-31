"""
Python Programmatic Identification of Novel Structures
Author: Jacob Zorn
Institution: The Pennsylvania State University
Version: 0.5.0

Program Description:

"""

# ------------------------------------------------------------------------------
# Import Modules
import matplotlib.pyplot as plt
import sympy as sp
import numpy as np
import pandas as pd
import numericaloptimization as no
import time

from pynhhd import nHHD

# ------------------------------------------------------------------------------
# Analytical Vortex Model
"""
To test the various calculations, a simple analytical model is utilized to
generate a flow field.
"""

class analyticalField:

    def __init__(self, spacing=12, nx=30, ny=30, vortexFuncs=''):
        self.spacing = spacing
        self.nx = nx
        self.ny = ny

        x,y = sp.symbols('x y')

        if vortexFuncs == '':
            self.ux = 1.0 * sp.cos(np.pi/spacing * x) * sp.cos(np.pi/spacing * y)
            self.uy = 1.0 * sp.sin(np.pi/spacing * x) * sp.sin(np.pi/spacing * y)
        else:
            self.ux = sp.cos(x + 2 * y)
            self.uy = sp.sin(x - 2 * y)
            # self.ux = vortexFuncs[0]
            # self.uy = vortexFuncs[1]

        self.gradient = sp.Matrix([[self.ux.diff(x), self.ux.diff(y)],
                                    [self.uy.diff(x), self.uy.diff(y)]])

        self.gradient = sp.lambdify([x,y], self.gradient, 'numpy')

        self.ux = sp.lambdify([x,y],self.ux,'numpy')
        self.uy = sp.lambdify([x,y],self.uy,'numpy')

        x = np.arange(0,nx+1)
        y = np.arange(0,ny+1)
        x,y = np.meshgrid(x,y)

        self.meshX = x
        self.meshY = y
        self.vectorX = self.ux(self.meshX,self.meshY)
        self.vectorY = self.uy(self.meshX,self.meshY)
        self.grad = self.gradient(x,y)

    def curl(self):

        self.vectorCurl = self.grad[1,0] - self.grad[0,1]

    def divergence(self):

        self.vectorDivergence = self.grad[0,0] + self.grad[1,1]

    def kineticEnergy(self):

        self.vectorEnergy = np.sqrt(self.vectorX**2 + self.vectorY**2)

    def reynoldsDecomp(self):

        meanU, meanV = self.vectorX.mean(), self.vectorY.mean()
        self.reynoldsX = self.vectorX - meanU
        self.reynoldsY = self.vectorY - meanV

    def galileanDecomp(self, varValue=0.7):

        meanU, meanV = self.vectorX * varValue, self.vectorY * varValue
        self.galileanX  = self.vectorX - meanU
        self.galileanY = self.vectorY - meanV

    def windingNumber(self):

        print('Need to implement')

    def qCriterion(self):

        self.criterionQ = -self.grad[1,0] * self.grad[0,1] + self.grad[0,0] * self.grad[1,1]

    def deltaCriterion(self):

        q = self.grad[0,0] + self.grad[1,1] + self.grad[0,0]**2 + self.grad[1,1]**2

        R = np.zeros((self.nx+1, self.ny+1))
        for i in range(self.nx+1):
            for j in range(self.ny+1):
                R[i,j] = np.linalg.det(self.grad[:,:,i,j])

        self.criterionDelta = ((q**3)/27) + ((R**2)/4)

    def lambdaCriterion(self):

        self.criterionLambda = self.grad[0,1]**2 + self.grad[1,1] * self.grad[0,0]

    def helmholtzDecomposition(self):

        nhhd = nHHD(grid=self.vectorX.shape, spacings=(1.0,1.0))
        nhhd.decompose(np.stack([self.vectorX, self.vectorY],axis=2))

        self.rotationalField = nhhd.r
        self.harmonicField = nhhd.h
        self.divergenceField = nhhd.d

# ------------------------------------------------------------------------------

class discreteField:

    def __init__(self, field):

        self.field = field
        gradC = np.asarray(np.gradient(self.field, axis=range(1,len(self.field)+1)))
        self.grad = gradC
        self.assignVectors()

    def curl(self):

        if self.field.shape[0] == 2:
            #Two Dimensional Discrete Curl field
            self.vectorCurl = self.grad[1,1,:,:] - self.grad[0,0,:,:]
        elif self.field.shape[0] == 3:
            #Three dimensional Discrete Curl field
            self.vectorCurl = 'Non implemented'
        else:
            print('Unrecognized function size. Please address.')

    def divergence(self):

        if self.field.shape[0] == 2:
            #Two dimensional discrete divergence field
            self.vectorDivergence = self.grad[0,0,:,:] + self.grad[1,1,:,:]
        elif self.field.shape[0] == 3:
            #Three dimensional discrete divergence field
            self.vectorDivergence = 'Non implemented'
            # self.vectorDivergence = self.grad[0,0,:,:] + self.grad[1,1,:,:] + self.grad[2,2,:,:]
        else:
            print('Unrecognized function size. Please address.')

    def kineticEnergy(self):

        if self.field.shape[0] == 2:
            self.vectorEnergy = self.vectorX**2 + self.vectorY**2
            self.vectorEnergy = self.vectorEnergy**0.5
        elif self.field.shape[0] == 3:
            self.vectorEnergy = self.vectorX**2 + self.vectorY**2 + self.vectorZ**2
            self.vectorEnergy = self.vectorEnergy**0.5
        else:
            print('Unrecognized function size. Please address.')

    def assignVectors(self):

        if self.field.shape[0] == 2:
            self.vectorX = self.field[0,...]
            self.vectorY = self.field[1,...]
        elif self.field.shape[0] == 3:
            self.vectorX = self.field[0,...]
            self.vectorY = self.field[1,...]
            self.vectorZ = self.field[2,...]
        elif self.field.shape[0] == 4:
            self.vectorX = self.field[0,...]
            self.vectorY = self.field[1,...]
            self.vectorZ = self.field[2,...]
            self.vectorW = self.field[3,...]
        else:
            print('Unrecognized function size. Please address.')

    def reynoldsDecomp(self):

        if self.field.shape[0] == 2:
            meanU, meanV = self.vectorX.mean(), self.vectorY.mean()
            self.reynoldsX = self.vectorX - meanU
            self.reynoldsY = self.vectorY - meanV
        elif self.field.shape[0] == 3:
            meanU, meanV, meanW = self.vectorX.mean(), self.vectorY.mean(), self.vectorZ.mean()
            self.reynoldsX = self.vectorX - meanU
            self.reynoldsY = self.vectorY - meanV
            self.reynoldsZ = self.vectorZ - meanW
        else:
            print('Unrecognized function size. Please address.')

    def galileanDecomp(self, varValue=0.7):

        if self.field.shape[0] == 2:
            meanU, meanV = self.vectorX * varValue, self.vectorY * varValue
            self.galileanX = self.vectorX - meanU
            self.galileanY = self.vectorY - meanV
        elif self.field.shape[0] == 2:
            meanU, meanV, meanW = self.vectorX * varValue, self.vectorY * varValue, self.vectorZ * varValue
            self.galileanX = self.vectorX - meanU
            self.galileanY = self.vectorY - meanV
            self.galileanZ = self.vectorZ - meanW
        else:
            print('Unrecognized function size. Please address.')

    def qCriterion(self):

        self.criterionQ = -self.grad[1,0,...] * self.grad[0,1,...] + self.grad[0,0,...] * self.grad[1,1,...]

    def lambdaCriterion(self):

        self.criterionLambda = self.grad[0,1,...]**2 + self.grad[1,1,...] * self.grad[0,0,...]

    def deltaCriterion(self):

        if self.field.shape[0] == 2:
            q = self.grad[0,0,...] + self.grad[1,1,...] + self.grad[0,0,...]**2 + self.grad[1,1,...]**2
            R = np.zeros((self.field.shape[1], self.field.shape[2]))
            for i in range(self.field.shape[1]):
                for j in range(self.field.shape[2]):
                    R[i,j] = np.linalg.det(self.grad[:,:,i,j])

            self.criterionDelta = ((q**3)/27) + ((R**2)/4)

    def helmholtzDecomposition(self):


        nhhd = nHHD(grid=self.vectorX.shape, spacings=tuple([1.0]*self.field.shape[0]))
        if self.field.shape[0] == 2:
            nhhd.decompose(np.stack([self.vectorX, self.vectorY],axis=2))
            self.rotationalField = nhhd.r
            self.harmonicField = nhhd.h
            self.divergenceField = nhhd.d
        else:
            print('Unrecognized function size. Please address.')



# ------------------------------------------------------------------------------
class testing:

    def runTest():
        p1 = analyticalField(spacing=10)
        p1.curl()
        p1.kineticEnergy()
        p1.qCriterion()
        p1.deltaCriterion()
        p1.lambdaCriterion()
        p1.helmholtzDecomposition()
        p1.windingNumber()

        plt.quiver(p1.rotationalField[:,:,0], p1.rotationalField[:,:,1])
        plt.show()

        plt.quiver(p1.divergenceField[:,:,0], p1.divergenceField[:,:,1])
        plt.show()

        plt.pcolor(p1.criterionLambda)
        plt.quiver(p1.meshX, p1.meshY, p1.vectorX,p1.vectorY)
        plt.colorbar()
        plt.show()
