#This script is written by Dr. Owais Khan on December 21, 2023
#This script does the following:
#	- Compute POD Eigen-spectra from SimVascular's results of simulated (coarse) velocity
#	- Compute POD Eigen-spectra from SimVascular's results of simulated (fine) velocity (optional)
# 	- Compute L2norm between coarse and fine velocity (optional)
#	- Compute L2norm along the centerline between coarse and fine velocity field (optional)
#	- Compute POD eigen-spectra along the centerline for coarse velocity (optional)
#	- Compute POD eigen-spectra along the centerline for fine velocity (optional)

########### NEEEDD TO FIX POD CL CODE ################
import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy 
from vtk.util.numpy_support import numpy_to_vtk 
import numpy as np
import vtk
import argparse
from utilities import *
import re
import modred

class DataAnalysis():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			self.Args.OutputFolder=self.Args.InputFolder+"/../DataAnalysis/"
			os.system("mkdir %s"%self.Args.OutputFolder)
	
		if self.Args.ClipPlane:
			self.Args.Plane=vtk.vtkPlane()
			self.Args.Plane.SetOrigin(self.Args.ClipPlaneOrigin)
			self.Args.Plane.SetNormal(self.Args.ClipPlaneNormal)
		

	def Main(self):
                #Read the Centerline
		if self.Args.CenterlineFile:
			print ("--- Reading Centerline File: %s"%self.Args.CenterlineFile)
			Centerline=ReadVTPFile(self.Args.CenterlineFile)
			CLCoords=np.array([Centerline.GetPoint(i) for i in range(Centerline.GetNumberOfPoints())])
			#Create a Distance Array as well
			CL_Length=[0]
			for i in range(1,len(CLCoords)):
				Distance_=CL_Length[i-1]+((CLCoords[i-1][0]-CLCoords[i][0])**2+(CLCoords[i-1][1]-CLCoords[i][1])**2+(CLCoords[i-1][2]-CLCoords[i][2])**2)**0.5
				CL_Length.append(Distance_)
			print ("------ The length of the geometry is is: %.05f"%CL_Length[-1])

		#Read the velocity data
		InputFiles=glob(self.Args.InputFolder+"/*.vtu")
		InputFiles=sorted(InputFiles, key=lambda s: int(re.search(r'\d+', s.split("/")[-1].split(".")[0]).group()))

		#Get the number of points for source velocity
		if self.Args.ClipPlane: Npts=ClipDataSet(ReadVTUFile(InputFiles[0]),self.Args.Plane).GetNumberOfPoints() #Number of coordinate points
		else: Npts=ReadVTUFile(InputFiles[0]).GetNumberOfPoints()
		Nt=len(InputFiles) #Number of time points
		
		#Get the number of points for the target velocity
		if self.Args.InputFolder2 is not None:
			InputFiles2=(glob(self.Args.InputFolder2+"/*.vtu"))
			InputFiles2=sorted(InputFiles2, key=lambda s: int(re.search(r'\d+', s.split("/")[-1]).group()))
			if self.Args.ClipPlane: Npts_target=ClipDataSet(ReadVTUFile(InputFiles2[0]),self.Args.Plane).GetNumberOfPoints()
			else: Npts_target=ReadVTUFile(InputFiles2[0]).GetNumberOfPoints()
			Nt_target=len(InputFiles2)
			if Npts!=Npts_target: 
				print ("Mesh points for InputFolder and InputFolder2 are not equal")
				print ("Exiting...")
				exit(1)
			elif Nt!=Nt_target:
				print ("The number of timepoints for InputFolder and InputFolder2 are not equal")
				print ("Exiting...")
				exit(1)
			else:
				print ("--- Number of time steps in InputFolder1 & InputFolder2:  %d"%Npts)
				print ("--- Number of Mesh Points in InputFolder1 & InputFolder2: %d"%Nt)
				print ("\n")

		#Create a Velocity Array
		Velocity=np.zeros(shape=(Npts,Nt))
		if self.Args.InputFolder2: 
			Velocity2=np.zeros(shape=(Npts_target,Nt_target))
			L2Norm=np.zeros(Nt)	

		#Assign the mesh point to the closest centerline point
		if self.Args.CenterlineFile:
			Mesh_ClosestCLId=[]
			File1=ReadVTUFile(InputFiles[0])
			#First find the closest centerline point to each mesh point:
			for i in range(Npts):
				coord_=File1.GetPoint(i)
				ClosestPoint_,ID_=ClosestPoint(coord_,CLCoords)
				Mesh_ClosestCLId.append(ID_)
		
			#Create a velocity array divided in centerline points
			VelocityCL=[[] for i in range(len(CLCoords))]
			if self.Args.InputFolder2:
				VelocityCL2=[[] for i in range(len(CLCoords))]
				L2NormCL=np.zeros(len(CLCoords))


		#WRITE A MESH FILE WITH CL TAGS ON FOR THE VOLUMETRIC MESH TO CHECK IF THE ASSIGNMENT IS CORRECT
		#Write the Centerline Sections to File.
		Mesh_ClosestCLId_VTK=numpy_to_vtk(np.array(Mesh_ClosestCLId))
		Mesh_ClosestCLId_VTK.SetName("CenterlineTags")
		#Remove extra arrays from the Velocity File
		NArrays_=File1.GetPointData().GetNumberOfArrays()
		for i in range(NArrays_): File1.GetPointData().RemoveArray(1) 
		File1.GetPointData().AddArray(Mesh_ClosestCLId_VTK)
		WriteVTUFile("%s/Velocity_CenterlineTags.vtu"%self.Args.OutputFolder,File1)
			

		#Compute global POD, L2Norm
		counter=0
		for i in range(len(InputFiles)):
			print ("------ Looping over: %s"%InputFiles[i])
			if self.Args.InputFolder2: print("--- Ground-truth File: %s"%InputFiles2[i])
			#Read the VTU Files for Velocity 1
			if self.Args.ClipPlane:
				Vel_=vtk_to_numpy(ClipDataSet(ReadVTUFile(InputFiles[i]),self.Args.Plane).GetPointData().GetArray(self.Args.ArrayName)) 
			else: Vel_=vtk_to_numpy(ReadVTUFile(InputFiles[i]).GetPointData().GetArray(self.Args.ArrayName))
			VelMag_=np.power(Vel_[:,0]**2+Vel_[:,1]**2+Vel_[:,2]**2,0.5)
			Velocity[:,counter]=VelMag_
			
			#Compute the L2Norm at the current timestep
			if self.Args.InputFolder2:
				#Read the VTU Files for Velocity 2
				if self.Args.ClipPlane: 
					Vel2_=vtk_to_numpy(ClipDataSet(ReadVTUFile(InputFiles2[i]),self.Args.Plane).GetPointData().GetArray(self.Args.ArrayName))
				else: Vel2_=vtk_to_numpy(ReadVTUFile(InputFiles2[i]).GetPointData().GetArray(self.Args.ArrayName))
				VelMag2_=np.power(Vel2_[:,0]**2+Vel2_[:,1]**2+Vel2_[:,2]**2,0.5)
				Velocity2[:,counter]=VelMag2_

				#Computing the L2Norm
				L2Norm_=(np.sum(np.power(VelMag_-VelMag2_,2)))**0.5
				L2Norm_/=(np.sum(np.power(VelMag2_,2)))**0.5
				L2Norm[i]=L2Norm_


			if self.Args.CenterlineFile:
				for j in range(0,Npts):
					CL_id_=Mesh_ClosestCLId[j]
					VelocityCL[CL_id_].append(VelMag_[j])
					if self.Args.InputFolder2: 
						VelocityCL2[CL_id_].append(VelMag2_[j])
				if self.Args.InputFolder2:
					for j in range(len(CLCoords)):
						L2NormCL_=(np.sum(np.power(np.array(VelocityCL[j])-np.array(VelocityCL2[j]),2)))**0.5
						L2NormCL_/=(np.sum(np.power(np.array(VelocityCL2[j]),2)))**0.5
						L2NormCL[j]+=L2NormCL_
						
			counter+=1

		#Get the time-averaged L2-norm along the centerline
		L2NormCL=L2NormCL/float(counter)
	
		#Compute the POD Eigen-spectra for InputFolder1
		print ("--- Computing POD Eigen-Spectra for InputFolder1")
		POD_res=modred.compute_POD_arrays_snaps_method(Velocity,list(modred.range(len(InputFiles))))
		EnergySpec=POD_res.eigvals/np.sum(POD_res.eigvals)

		print ("--- Writing POD sepectra for InputFolder1 in: %s"%self.Args.OutputFolder+"/POD_InFile1.dat")
		outfile=open(self.Args.OutputFolder+"/POD_InFile1.dat",'w')
		outfile.write("Mode(k) %Energy\n")
		for i in range(len(EnergySpec)):outfile.write(str(i)+" "+str(EnergySpec[i])+"\n")
		outfile.close()

		#Compute the POD Eigen-spectra for InputFolder2
		if self.Args.InputFolder2:
			print ("--- Computing POD Eigen-Spectra for InputFolder2")
			POD_res=modred.compute_POD_arrays_snaps_method(Velocity2,list(modred.range(len(InputFiles))))
			EnergySpec=POD_res.eigvals/np.sum(POD_res.eigvals)
			#Write the EigenSpectra
			print ("--- Writing POD spectra for InputFolder2 in: %s"%self.Args.OutputFolder+"/POD_InFile2.dat")
			outfile=open(self.Args.OutputFolder+"/POD_InFile2.dat",'w')
			outfile.write("Mode(k) %Energy\n")
			for i in range(len(EnergySpec)):outfile.write(str(i)+" "+str(EnergySpec[i])+"\n")
			outfile.close()

			#Compute the L2 Norms
			print ("------ Writing L2-Norm: %s/L2norm_time.dat"%self.Args.OutputFolder)
			#Write L2Norm over time
			outfile=open(self.Args.OutputFolder+"/L2norm_time.dat",'w')
			outfile.write("#Time-Averaged L2Norm: %.010f\n"%np.average(L2Norm))
			outfile.write("Timestep L2Norm\n")	
			for i in range(len(L2Norm)):outfile.write("%.010f %.010f\n"%(i,L2Norm[i]))
			outfile.close()

		if self.Args.CenterlineFile and self.Args.InputFolder2:
			print ("------ Writing time-averaged L2-Norm CL: %s/L2Norm_CL_timeaveraged.dat"%self.Args.OutputFolder)
			outfile=open(self.Args.OutputFolder+"/L2norm_CL_timeaveraged.dat",'w')
			outfile.write("Length L2Norm\n")
			for i in range(len(CLCoords)): outfile.write("%.010f %.010f\n"%(CL_Length[i],L2NormCL[i]))
			outfile.close()
				
		if self.Args.CenterlineFile:
			print ("------ Computing and Writing POD Eigen-Spectra (sum of > Mode2) Along the Centerline for InputFolder1")
			outfile=open(self.Args.OutputFolder+"/POD_CL_InFile1.dat",'w')
			outfile.write("Mode(k) Energy>%d\n"%self.Args.ModeCutOff)
			for i in range(len(CLCoords)):
				#POD Eigen-values after the cut-off mode 
				VelocityCL_=np.array(VelocityCL[i])
				VelocityCL_.resize(int(len(VelocityCL_)/Nt),Nt)
				POD_res_=modred.compute_POD_arrays_snaps_method(VelocityCL_,list(modred.range(Nt)))
				EnergySpec_=sum((POD_res.eigvals/np.sum(POD_res.eigvals))[2:])
				outfile.write("%.010f %.010f\n"%(CL_Length[i],EnergySpec_))
			outfile.close()
			
		if self.Args.CenterlineFile and self.Args.InputFolder2:
			print ("------ Computing and Writing POD Eigen-Spectra (sum of >Mode2) Along the Centerline for InputFolder2")
			outfile=open(self.Args.OutputFolder+"/POD_CL_InFile2.dat",'w')
			outfile.write("Mode(k) Energy>%d\n"%self.Args.ModeCutOff)
			for i in range(len(CLCoords)):
				#POD Eigen-values after the cut-off mode
				VelocityCL2_=np.array(VelocityCL2[i])
				VelocityCL2_.resize(int(len(VelocityCL2_)/Nt),Nt)
				POD_res_=modred.compute_POD_arrays_snaps_method(VelocityCL2_,list(modred.range(Nt)))
				EnergySpec_=sum((POD_res.eigvals/np.sum(POD_res.eigvals))[2:])
				outfile.write("%.010f %.010f\n"%(CL_Length[i],EnergySpec_))
			outfile.close()
	
					
				

if __name__=="__main__":
         #Description
	parser = argparse.ArgumentParser(description="This script will compute quantities along a given centerline. A SimVascular velocity folder must be provided.")

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="The  input folder that contains the velocity files from SimVascular")
	
	parser.add_argument('-InputFolder2', '--InputFolder2', type=str, required=False, dest="InputFolder2",help="The  input folder that contains the ground-truth velocity for computing errors such as L2-norms etc.")


	parser.add_argument('-ClipPlane', '--ClipPlane', type=int, required=False, dest="ClipPlane", help="Whether the dataset should be clipped before performing the analysis. Set to 1 to turn on.")
	
	parser.add_argument('-ClipPlaneOrigin', '--ClipPlaneOrigin', type=ListOfFloats, required=False, dest="ClipPlaneOrigin", default="0.06356052780154019,0.017335321237208618,0.5993404499350882",help="The origin of the plane to clip before analyzing the data.")
	
	parser.add_argument('-ClipPlaneNormal', '--ClipPlaneNormal', type=ListOfFloats, required=False, dest="ClipPlaneNormal", default="-0.09153368148760842,-0.03734511295045089,-0.9951014660284866",help="The normal of the plane to clip before analyzing the data")
	
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False,default="velocity", dest="ArrayName",help="The array name to process. Default is velocity")
	
	parser.add_argument('-Period', '--Period', type=int, required=False,default=1, dest="Period",help="The duration of the cardiac cycle. Default=1 second.")
	
	parser.add_argument('-CenterlineFile', '--CenterlineFile', type=str, required=False, dest="CenterlineFile",help="The file containing the centerline.")
	
	parser.add_argument('-ModeCutOff', '--ModeCutOff', type=int, required=False, default=2, dest="ModeCutOff",help="The POD modes after which we consider the flow to be turbulent.")

 	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="An output folder to store the centerlines and data")
	args=parser.parse_args()
	DataAnalysis(args).Main()	
