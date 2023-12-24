import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class ProperOrthogonalDecomposition():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFileName is None:
			self.Args.OutputFileName=self.Args.InputFolder+"/EigenValues.dat"			

	def Main(self):
                #Read all of the file name
		print ("--- Finding all of the .vtu files in: %s"%self.Args.InputFolder)
		InputFiles=sorted(glob(self.Args.InputFolder+"/*.vtu")) #volumetric files 
		print ("--- Number of Time Points: %d"%len(InputFiles))
                
		#Loop over all of the file names and store the values
		File1=ReadVTUFile(InputFiles[0])

                #Get the number of points
		NPoints=File1.GetNumberOfPoints()
		print ("--- Number of Mesh Points: %d"%NPoints)

		#Read all of the velocity data
		Vx=np.zeros(shape=(NPoints,len(InputFiles)))	
		Vy=np.zeros(shape=(NPoints,len(InputFiles)))	
		Vz=np.zeros(shape=(NPoints,len(InputFiles)))	
		
		#Loop over all of the files
		for i in range(len(InputFiles)): #Loop over all of the files
			print ("------ Looping over: %s"%InputFiles[i])
			VelocityFile_=ReadVTUFile(InputFiles[i]) #Read the velocity file
			VelocityData_=vtk_to_numpy(VelocityFile_.GetPointData().GetArray(self.Args.ArrayName))
			Vx[:,i]=VelocityData_[:,0]
			Vy[:,i]=VelocityData_[:,1]
			Vz[:,i]=VelocityData_[:,2]
	
		#Compute the autocorrelation matrix
		print ("--- Computing Auto-Correlation Matrix")
		AA=np.matmul(Vx.transpose(),Vx)
		AA+=np.matmul(Vy.transpose(),Vy)
		AA+=np.matmul(Vz.transpose(),Vz)
	
                # Compute eigenvalues
		print ("--- Computing POD Eigen Values and Eigen Vectors")
		eigval, eigvec = np.linalg.eigh(AA)

                # Sort eigenvalues and eigenvectors from high to low
		ind = eigval.argsort()
		indrev = ind[::-1]
		eigval = eigval[indrev]
		eigvec = eigvec[:,indrev]

		#Save the Eigen Values and Eigen Vectors
		print ("--- Storing the Normalized Eigen Values in %s"%self.Args.OutputFileName)
		eigval_normalized=eigval/np.sum(eigval)
		outfile=open(self.Args.OutputFileName,'w')
		for i in range(0,len(eigval)):
			outfile.write("%d %.08f %.08f\n"%(i+1,eigval_normalized[i], eigval[i]))
		outfile.close()




	

if __name__=="__main__":
        #Arguments 
	parser= argparse.ArgumentParser(description="This script will compute the eigen values and eigenvectors of velocity field using POD")   

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The input folder that contains all of the results file, taged as *.vtu")
                
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName", help="The file name of the to store the POD modes.")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName",default="velocity", help="The file name of the to store the POD modes.")


        #Put all the arguments together
	args=parser.parse_args()
                
        #Call your Class
	ProperOrthogonalDecomposition(args).Main()
