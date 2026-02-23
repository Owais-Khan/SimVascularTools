import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class ExtractVelocityAlongCenterline():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			self.Args.OutputFolder=os.path.join(os.path.split(self.Args.InputFileName)[::-1]) 
	def Main(self):
		#Read the Files
		print ("--- Reading the PostProcessed Results: %s"%self.Args.InputFileName)
		InputFiles=ReadVTUFile(self.Args.InputFileName)
		
		print ("--- Reading Centerline File: %s"%self.Args.CenterlinesFile)	
		Centerlines=ReadVTPFile(self.Args.CenterlinesFile)
	
		#Create an array to store velocity data
		Velocity=[[] for i in range(len(InputFiles))]
		Pressure=[[] for i in range(len(InputFiles))]
		SPI=[[] for i in range(len(InputFiles))]
		PointIDs=[]
		XYZ=[]

		#Loop over all of the points in the centerline
		for i in range(Centerlines.GetNumberOfPoints()):
			print ("--- Looping over %d of %d Centerline Points"%(i,Centerlines.GetNumberOfPoints()))
			#Read the coordinate
			Point_=Centerlines.GetPoint(i)
			Radius_=Centerlines.GetPointData().GetArray("MaximumInscribedSphereRadius").GetValue(i)
				
				
			#Read the VTK or VTU files
			InputFile_=ReadVTUFile(InputFiles[j])
			#Get the IDs of Points within given radius
			if j==0: 
				IDs_,XYZ=self.GetPointsInRadius(InputFile_,Point_,Radius_,XYZ) #Get ids and XYZ
				for value in IDs_: PointIDs.append(value)#Keep track of point ids
				print (len(IDs_))

			for k in (len(IDs_)):
				Velocity[j].append(InputFile_.GetPointData().GetArray("Velocity").GetValue(IDs_[k]))
				Pressure[j].append(InputFile_.GetPointData().GetArray("Pressure").GetValue(IDs_[k]))
				SPI[j].append(InputFile_.GetPointData().GetArray("SPI").GetValue(IDs_[k]))
				ReynoldsNumber[j].append((Velocity[j]*Radius_*2)/self.Args.Viscosity)	
					
		#Now write out the data for all of the Timesteps
		print ("\n")
		print ("--- Writing the Output Data Along Centerlines")
		for i in range(len(InputFiles)):
			outfile_=self.Args.OutputFolder+"/HemodynamicsAlongCenterline.dat"%i
			print ("------ Writing: %s"%outfile_)
			outfile=open(outfile_,'w')
			outfile.write("X Y Z Distance ReynoldsNumber Velocity Pressure SPI\n")
			for j in range(len(XYZ)):
				outfile.write("%.05f %.05f %.05f %.05f %.05f %.05f\n"%(XYZ[j][0],XYZ[j][1],XYZ[j][2],ReynoldsNumber[j],Velocity[j],Pressure[j],SPI[j]))
			outfile.close()
					


	

	def GetPointsInRadius(self,Array,Point,Radius,XYZ):
		IDs_=[]
		for i in range (Array.GetNumberOfPoints()):
			Coord_=Array.GetPoint(i)
			if Coord_ in XYZ: continue
			else: 
				Dist_=((Coord_[0]-Point[0])**2+(Coord_[1]-Point[1])**2+(Coord_[2]-Point[2])**2)**0.5
				if Dist_<=Radius:
					IDs_.append(i)
					XYZ.append(Coord_)
		return IDs_,XYZ
			


if __name__=="__main__":
        #Arguments
        
	parser= argparse.ArgumentParser(description="This script will extract velocity along the centerline data.")

       
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="Folder containing the velocity data over a cardiac cycle.")
        
	parser.add_argument('-CenterlinesFile', '--CenterlinesFile', type=str, required=True, dest="CenterlinesFile", help="Centerline file of the geometry")

	
	parser.add_argument('-Viscosity', '--Viscosity', type=float, required=False, default=0.037735,  dest="Viscosity", help="The viscosity of blood. Default is 0.037735 cm/s2.")
	
	
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder", help="The output folder to store the velocity and coordinate data.")
	
	
        #Put all the arguments together
	args=parser.parse_args()

    #Call your Class
	ExtractVelocityAlongCenterline(args).Main()

