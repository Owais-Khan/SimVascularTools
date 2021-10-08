import vtk
import numpy
import os
import argparse
from glob import glob

class PostProcess3DResults():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Get the Flow Rates at all of the outlets
		CapNames=self.ReadCapNames() #Read the CapNames
		CapData =[self.ReadVtpFile(CapName_) for CapName_ in CapNames]
		CapAreas=[self.ComputeSurfaceArea(CapData_) for CapData_ in CapData] #Compute the Areas of each cap

		#Create a Dictionary to Store the Velocity Data
		Velocity={}	
		Pressure={}		

		#Get All the results file
		for i in range(self.Args.StartTimestep,self.Args.StopTimestep,self.Args.Increment):
			
			#Read the filenames at the current timestep
			ResultVolumetricFileName_=glob("%s/*%.05d.vtu"%(self.Args.ResultsFolder,i))[0]
			ResultSurfaceFileName_   =glob("%s/*%.05d.vtp"%(self.Args.ResultsFolder,i))[0]
			if len(ResultVolumetricFileName_)==0 or len(ResultSurfaceFileName_)==0: 
				print ("Could Not Find the File at TimeStep: %.05d"%i)
		
			#Now Load the Volumetric File
			print ("------ Loading %s"%ResultVolumetricFileName_)
			VolumetricData_=self.ReadVtuFile(ResultVolumetricFileName_)	
	
			#Loop over all of the caps to compute the average velocity

	
	def ReadCapNames(self):
		#Grab all of the file names from the mesh-complete folder
		Caps_LCA  =sorted(glob(self.Args.InputMeshFolder+"/mesh-surfaces/lca_*.vtp"))	
		Caps_RCA  =sorted(glob(self.Args.InputMeshFolder+"/mesh-surfaces/rca_*.vtp"))	
		Caps_Aorta=sorted(glob(self.Args.InputMeshFolder+"/mesh-surfaces/aorta_*.vtp"))
		
		#Put all of these names in a singe array
		CapNames=[]+Caps_Aorta+Caps_LCA+Caps_RCA

		return CapNames	

	def ReadVtpFile(self,FileName):
		reader=vtk.vtkXMLPolyDataReader()
		reader.SetFileName(FileName)
		reader.Update()
		return reader.GetOutput()

	def ReadVtuFile(self,Filename):
		reader=vtk.vtkXMLUnstructuredGridReader()
		reader.SetFileName(Filename)
		reader.Update()
		return reader.GetOutput()
		
	def ComputeSurfaceArea(self,PolyData):
		MassProperties=vtk.vtkMassProperties()
		MassProperties.SetInputData(PolyData)
		MassProperties.Update()
		return MassProperties.GetSurfaceArea()			




if __name__=="__main__":
	#Arguments
	parser= argparse.ArgumentParser(description="This script will process the results from Simvascular")

	#Input filename for the mesh-complete folder
	parser.add_argument('-meshfolder', '--InputMeshFolder', type=str, required=True, dest="InputMeshFolder", help="The path to the mesh-complete folder that contains the volumetric and surface meshes") 

	#Input Folder where all of the results have been stored
	parser.add_argument('-resultsfolder', '--ResultsFolder', type=str, required=True, dest="ResultsFolder", help="The path to the folder that contains all all_data.vtp and all_data.vtu files postprocessed from SimVascular")

	#Start and End Argumenets		 
	parser.add_argument('-start', '--StartTimestep', type=int, required=True, dest="StartTimestep", help="The starting timestep to process the SimVascular files")
	
	#End Argument
	parser.add_argument('-stop', '--StopTimestep', type=int, required=True, dest="StopTimestep", help="The End timestep to process the SimVascular files")

	#Increment
	parser.add_argument('-incr', '--Increment', type=int, required=True, dest="Increment", help="The increment for the timestep files to process")



	#Put all the arguments together
	args=parser.parse_args()	

	#Call your Class
	PostProcess3DResults(args).Main()
