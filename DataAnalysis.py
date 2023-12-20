import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy as npvtk
import numpy as np
import vtk
import argparse
from utilities import *
import modred

class PlotAlongCenterline():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
                #Read the Centerline
		if self.Args.CenterlineFile is not None:
			Centerline=ReadVTPFile(self.Args.CenterlinesFile)
			CLCoords=np.array(Centerline.GetPoints(i) for i in range(Centerline.GetNumberOfPoints()))
		

		#Read the velocity data
		InputFiles=sorted(glob(self.Args.InputFolder+"/*.vtu"))

		#Get the number of points for source velocity
		Npts=ReadVTUFile(InputFiles[0]).GetNumberOfPoints() #Number of coordinate points
		Nt=len(InputFiles) #Number of time points
		
		#Get the number of points for the target velocity
		if self.Args.InputFolder2 is not None:
			InputFiles2=sorted(glob(self.Args.InputFolder2+"/*.vtu"))
			Npts_target=ReadVTUFile(InputFiles2).GetNumberOfPoints()
			Nt_target=len(InputFiles2)
			if Npts!=Npts_target: 
				print ("Mesh points for InputFolder and InputFolder2 are not equal")
				print ("Exiting...")
				exit(1)
			if Nt!=Nt_target:
				print ("The number of timepoints for InputFolder and InputFolder2 are not equal")
				print ("Exiting...")
				exit(1)
			

		#Create a Velocity Array
		Velocity=np.zeros(shape=(Npts,Nt))
		if self.Args.InputFolder2: 
			Velocity2=np.zeros(shape=(Npts_target,Nt_target))
			L2Norm=np.zeros(Nt)	

		#Assign the mesh point to the closest centerline point
		if self.Args.CenterlineFile is not None:
			CLMeshArray=[[] for i in range(len(CLCoords))]
			File1=ReadVTUFile(InputFiles)
			#First find the closest centerline point to each mesh point:
			for i in range(Npts):
				coord_=File1.GetPoint(i)
				ClosestPoint_,ID_=ClosestPoint(coord_,CLCoords)
				CLMeshArray[ID_].append(i)
		
			#Create a velocity array divided in centerline points
			VelocityCL=[np.zeros(shape=(len(CLMeshArray[i]),Nt)) for i in range(len(CLMeshArray))]
			L2NormCL=np.array(len(CLCoords))	

		#Compute global POD, L2Norm
		counter=0
		for i in range(len(InputFiles)):
			print ("------ Looping over: %s"%InputFiles[i])
			#Read the VTU Files for Velocity 1
			Vel_=vtk_to_numpy(ReadVTUFile(InputFiles[i]).GetPointData().GetArray(self.Args.ArrayName))
			VelMag_=np.power(Vel_[:,0]**2+Vel_[:,1]**2+Vel_[:,2]**2,0.5)
			Velocity[:,counter]=VelMag_
			
			#Compute the L2Norm at the current timestep
			if self.Args.InputFolder2:
				#Read the VTU Files for Velocity 2
				Vel2_=vtk_to_numpy(ReadVTUFile(InputFiles2[i]).GetPointData().GetArray(self.Args.ArrayName))
				VelMag2_=np.power(Vel2_[:,0]**2+Vel2_[:,1]**2+Vel2_[:,2]**2,0.5)
				Velocity2[:,counter]=VelMag2_

				#Computing the L2Norm
				L2Norm_=(np.sum(np.power(VelMag_[:,0]-VelMag2_[:,0],2)))**0.5
				L2Norm_/=(np.sum(np.power(VelMag2_,2)))**0.5
				L2Norm.append(L2Norm_)


			#Loop over the centerline points and compute quantities 
			#if (self.Args.CenterlineFile is not None):
				print ("--------- Centerline Analaysis ---------")
				#Loop over all of the centerline coordinates
				#for j in range(len(CLCoords)):
					#For each centerline coordinate, store velocity
				#	for k in range(len(CLMeshArray[j])):	
									
		
			
			counter+=1
	
		#Compute the POD Eigen-spectra
		print ("--- Computing POD Eigen-Spectra")
		POD_res=modred.compute_POD_arrays_snaps_method(Velocity,list(modred.range(60)))
		print(POD_res.modes)
		print(POD_res.eigvals)	


		print ("--- Writing out the results")
		#Write the L2 Norm if a ground-truth velocity file was provided.	
		if self.Args.InputFolder2 is not None: 
			print ("------ writing L2Norms: %s/L2norm_time.dat"%self.Args.OutputFolder)
			#Write L2Norm over time
			outfile=open(self.Args.OutFolder+"/L2norm_time.dat")
			outfile.write("Timestep Time L2Norm\n")	
			for i in range(len(L2Norm)): outfile.write("%.05f %.05f %.05f\n"%(i,float(i)/len(L2Norm),L2Norm[i]))
			outfile.close()


if __name__=="__main__":
         #Description
	parser = argparse.ArgumentParser(description="This script will compute quantities along a given centerline. A SimVascular velocity folder must be provided.")

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="The  input folder that contains the velocity files from SimVascular")
	
	parser.add_argument('-InputFolder2', '--InputFolder2', type=str, required=False, dest="InputFolder2",help="The  input folder that contains the ground-truth velocity for computing errors such as L2-norms etc.")
	
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False,default="velocity", dest="ArrayName",help="The array name to process. Default is velocity")
	
	parser.add_argument('-Period', '--Period', type=int, required=False,default=1, dest="Period",help="The duration of the cardiac cycle. Default=1 second.")
	
	parser.add_argument('-CenterlineFile', '--CenterlineFile', type=str, required=False, dest="CenterlineFile",help="The file containing the centerline.")

 	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="An output folder to store the centerlines and data")
	args=parser.parse_args()
	PlotAlongCenterline(args).Main()	
