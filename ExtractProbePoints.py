import numpy as np
import vtk
from glob import glob
import argparse
from utilities import *
import os

class ExtractProbePoints():
	def __init__(self,args):
		self.Args=args
		if self.Args.OutputFolder is None:
			self.Args.OutputFolder=self.Args.InputFolder+"/ProbePointData"
			os.system("mkdir %s/ProbePointData"%self.Args.InputFolder)

	def Main(self):
		#Get the list of files
		FileNames=sorted(glob(self.Args.InputFolder+"/all_results.vtu_*.vtu"))

		#Read the Probe Points
		ProbePoints=[]
		Radius=[]
		ProbePointFile=open(self.Args.InputFile,'r')
		for Line in ProbePointFile:
			line=Line.split()
			ProbePoints.append([float(line[0]),float(line[1]),float(line[2])])
			Radius.append(float(line[1]))
		ProbePointFile.close() 
		
		#Get the Probe Points
		#Read the first velocity file
		VelocityFile0=ReadVTUFile(FileNames[0])
			
		#Get the Coordinates
		NPoints=VelocityFile0.GetNumberOfPoints()
		MeshPoints=np.array([VelocityFile0.GetPoint(i) for i in range(NPoints)])

		#Create a dictionary to store probe data
		ProbeData={}
		VelocityProbeData={}
		for i in range(len(ProbePoints)): #Loop over all of the probe points
			#For each probe point, find the closest nodes that fall within radius
			ProbeData[i]=[]
			ProbeData[i]={"Points":[],"PointIds":[],"Distance":[]}
			VelocityProbeData[i]=[]
				
			#Compute the distance from probe point
			Distance_=np.sum((MeshPoints - ProbePoints[i])**2, axis=1)
			SortArray=np.argsort(Distance_)
				
			#Find the closest point to the probe point within give radius
			for j in range(len(SortArray)):
				if Distance_[SortArray[j]]<Radius[i]:
					ProbeData[i]["Points"].append(MeshPoints[SortArray[j]])
					ProbeData[i]["PointIds"].append(SortArray[j])
					ProbeData[i]["Distance"].append(Distance_[SortArray[j]])
				else:
					break

			#Save the Dictionary to the output folder
		print ("Saving Probe Coords, Ids and Distances in %s\n"%(self.Args.OutputFolder))
		np.save("%s/ProbePointLocations.npy"%(self.Args.OutputFolder),ProbeData)

		#Create an array to store velocity data
		for i in range(len(ProbePoints)):
			VelocityProbeData[i]=np.zeros(shape=(len(ProbeData[i]["Points"]),len(FileNames)))

	
		#Loop over the velocity files and extract the data around each probe point
		counter=0
		for FileName in FileNames:
			print ("--- Looping over: %s"%FileName)
			#Read the data in the file
			Velocity_=ReadVTUFile(FileName)
			
			#Loop over all of the probe points
			for i in range(len(ProbePoints)):
				print ("------ Finished Probe Point %d"%i)

				#Loop over all the points with reach of the Radius,R
				for j in range(len(ProbeData[i]["Points"])):
					#Get Velocity Data from the vtu file
					Id_=ProbeData[i]["PointIds"][j]
					Vx_=Velocity_.GetPointData().GetArray("velocity").GetValue(Id_*3)
					Vy_=Velocity_.GetPointData().GetArray("velocity").GetValue(Id_*3+1)
					Vz_=Velocity_.GetPointData().GetArray("velocity").GetValue(Id_*3+2)
				
					#Compute the velocity magnitude
					Vmag=np.sqrt(Vx_**2+Vy_**2+Vz_**2)
					
					#Store the velocity data in Dictionary
					VelocityProbeData[i][j,counter]=Vmag	
			counter+=1

		#Store the Velocity Data
		print ("Saving Velocity data to %s\n"%(self.Args.OutputFolder))
		np.save("%s/VelocityData.npy"%(self.Args.OutputFolder),ProbeData)

		


if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will extrac pressure/velocity data from vtu files generated by SimVascular")

        #Provide a path to the Magnitude Images
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="This folder contains the velocity files in vtu format.")

	parser.add_argument('-InputFile', '--InputFile', type=str, required=True, dest="InputFile",help="The input file contains a column list of x y z R (coordinates and Radius) for which data is needed.")
	
	parser.add_argument('-OutputFolder', '--OutputFolder', type=int, required=False, dest="OutputFolder",help="The output folder to store the data")

	#Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	ExtractProbePoints(args).Main()

