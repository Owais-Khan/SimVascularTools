import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import copy
class ExtractVelocityAwayFromWall():
	def __init__(self,Args):
		self.Args=Args
		os.system("mkdir %s"%self.Args.OutputFolder)        
	def Main(self):
		#Read the Files
		print ("--- Getting Velocity File Names...")
		InputFiles=sorted(glob(self.Args.InputFolder+"/*.vtu"))
		
		#Read the Surface file containing wall mesh	
		SurfaceFile=ReadVTPFile(self.Args.SurfaceFile)

		#Read First VTU File and Get Point Array
		MeshCoords_   =vtk_to_numpy(ReadVTUFile(InputFiles[0]).GetPoints().GetData())
		SurfaceCoords_=vtk_to_numpy(SurfaceFile.GetPoints().GetData())

		#File the Distance Tags
		CenterIDs=[]
		print ("---Tag the near wall surface")
		for i in range(len(MeshCoords_)):
			point_,id_,dist_=ClosestPoint(MeshCoords_[i],SurfaceCoords_)
			if dist_<self.Args.WallDistance: #Tag near wall nodes 
				CenterIDs.append(i)
		print ("Number of Near Wall Nodes Tagged: %d"%len(CenterIDs))
		print ("Total Number of Mesh Nodes: %s"%len(MeshCoords_))
	
		#Loop over all of the Centerlines File
		counter=0
		No_Sensors=0
		for FileName in InputFiles:
			#print ("--- Looping over: %s"%FileName)
			Data_=ReadVTUFile(FileName)
			velocity_=vtk_to_numpy(Data_.GetPointData().GetArray("Velocity"))

			#Make a deep copy of velocity
			velocity2_=np.zeros(shape=(len(velocity_),3))
			for i in range(len(velocity_)):
				velocity2_[i][0]=velocity_[i][0]
				velocity2_[i][1]=velocity_[i][1]
				velocity2_[i][2]=velocity_[i][2]
		
                       	#Make wall nodes zero 
			for value in CenterIDs:
				velocity2_[value][0]=0.0  
				velocity2_[value][1]=0.0 
				velocity2_[value][2]=0.0

			if counter==0:
				for i in range(len(velocity2_)):
					if np.linalg.norm(velocity2_[i])!=0.0: No_Sensors+=1
					

			velocityVTK_=numpy_to_vtk(velocity2_,deep=True)
			velocityVTK_.SetName("Velocity")
			Data_.GetPointData().AddArray(velocityVTK_)
			Data_.Modified()
			outfile_=self.Args.OutputFolder+"/"+FileName.split("/")[-1]
			#print ("---Writing Output File: %s"%outfile_)
			WriteVTUFile(outfile_,Data_)
			os.system("vmtkmeshscaling -ifile %s -ofile %s -scale 0.1"%(outfile_,outfile_))
			
			counter=+1

		print ("There are %d Sensor Points"%No_Sensors)	
			


if __name__=="__main__":
        #Arguments
        
	parser= argparse.ArgumentParser(description="This script will extract velocity along the centerline data.")

       
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="Folder containing the velocity data over a cardiac cycle.")
        
	parser.add_argument('-SurfaceFile', '--SurfaceFile', type=str, required=True, dest="SurfaceFile", help="Surface file of the wall mesh ")
	
	parser.add_argument('-WallDistance', '--WallDistance', type=float, required=False, default=2.5, dest="WallDistance", help="DistanceFromWall")

	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=True, dest="OutputFolder", help="The output folder to store the velocity and coordinate data.")

	
        #Put all the arguments together
	args=parser.parse_args()

    #Call your Class
	ExtractVelocityAwayFromWall(args).Main()

