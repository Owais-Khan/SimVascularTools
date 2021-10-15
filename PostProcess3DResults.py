import vtk
import numpy
import os
import argparse
from glob import glob
import numpy as np
class PostProcess3DResults():
	def __init__(self,Args):
		self.Args=Args
		self.PressureScale=0.000750061683
	def Main(self):
		#Get the Flow Rates at all of the outlets
		CapNames=self.ReadCapNames() #Read the CapNames
		CapData =[self.ReadVtpFile(CapName_) for CapName_ in CapNames]
		CapAreas=[self.ComputeSurfaceArea(CapData_) for CapData_ in CapData] #Compute the Areas of each cap

		#Create a Dictionary to Store the Velocity Data
		Data_VelPres={}	
		for CapName_ in CapNames: Data_VelPres[CapName_]=[]
		

		#Get All the results file
		for i in range(self.Args.StartTimestep,self.Args.StopTimestep,self.Args.Increment):
	
			#Read the filenames at the current timestep
			#ResultVolumetricFileName_=glob("{0}.vtu_{1}.vtu".format(self.Args.ResultsFolder,i))[0]	
			ResultVolumetricFileName_=glob("%s/*%.05d.vtu"%(self.Args.ResultsFolder,i))[0]
			#ResultSurfaceFileName_   =glob("{0}.vtp_{1}.vtp".format(self.Args.ResultsFolder,i))[0]
			ResultSurfaceFileName_   =glob("%s/*%.05d.vtp"%(self.Args.ResultsFolder,i))[0]
			if len(ResultVolumetricFileName_)==0 or len(ResultSurfaceFileName_)==0: 
				print ("Could Not Find the File at TimeStep: %.05d"%i)
		
			#Now Load the Volumetric File
			print ("------ Loading %s"%ResultVolumetricFileName_)
			VolumetricData_=self.ReadVtuFile(ResultVolumetricFileName_)	
			
			
			#Loop over all of the caps to compute the average velocity
			for k in range(0,len(CapNames)):

				#Get the Surface Node IDs for the Given Cap
				Npts_=CapData[k].GetNumberOfPoints()

				#Loop over all of the nodes of the cap ids
				VelocityX_ =np.zeros(Npts_)
				VelocityY_ =np.zeros(Npts_)
				VelocityZ_ =np.zeros(Npts_)
				Pressure_ =np.zeros(Npts_)

				for j in range(Npts_):
					NodeId_=CapData[k].GetPointData().GetArray("GlobalNodeID").GetValue(j)

					VelocityX_[j]=VolumetricData_.GetPointData().GetArray("velocity").GetValue(NodeId_*3)
					VelocityY_[j]=VolumetricData_.GetPointData().GetArray("velocity").GetValue(NodeId_*3+1)
					VelocityZ_[j]=VolumetricData_.GetPointData().GetArray("velocity").GetValue(NodeId_*3+2)

					Pressure_[j] =VolumetricData_.GetPointData().GetArray("pressure").GetValue(NodeId_)


				#Get Velocity Magnitude
				VelocityMag_=np.power(np.power(VelocityX_,2)+np.power(VelocityY_,2)+np.power(VelocityZ_,2),0.5)
				#Get the Average
				VelocityX_avg_=np.average(VelocityX_)
				VelocityY_avg_=np.average(VelocityY_)
				VelocityZ_avg_=np.average(VelocityZ_)
				VelocityMag_avg_=np.average(VelocityMag_)
				Pressure_avg_ =np.average(Pressure_)
				
				#Get the Max Velocity
				VelocityX_max_=np.max(VelocityX_)
				VelocityY_max_=np.max(VelocityY_)
				VelocityZ_max_=np.max(VelocityZ_)
				VelocityMag_max_=np.max(VelocityMag_)
			
				#Store this inside the global Velocity and Pressure Array
				Data_VelPres[CapNames[k]].append([VelocityX_avg_,VelocityY_avg_,VelocityZ_avg_,VelocityMag_avg_,VelocityX_max_,VelocityY_max_,VelocityZ_max_,VelocityMag_max_,Pressure_avg_,VelocityMag_avg_*CapAreas[k]])


		#Write the Output Tecplot File
		self.WriteTecplotOutlets(Data_VelPres)

	
	def WriteTecplotOutlets(self,Data):
		outfile=open(self.Args.OutFolder+"/OutletData.dat",'w')
		#In the file header, add all of the time-averaged quantities
		outfile.write("#Time Averaged Quantities\n")
		outfile.write("VelocityMag, VelocityPeak, FlowRate, Pressure\n")
		for key in Data:
			print (key,Data[key])
			#print (np.array(Data[key])[:])
			exit(1)
			
			outfile.write("#%.05f %.05f %.05f %.0f5\n"%(np.average(Data[key][:][3]),np.average(Data[key][:][7]),Data[key][:][9],Data[key][:][8]*self.PressureScale))
		
		outfile.write('TITLE="Outlet FLow Rates and Velocities"\n')
		outfile.write('VARIABLES="Time", "VelocityMag","VelocityMagPeak","FlowRate","Pressure"\n')

		#Define the period
		Period=60/self.Args.HeartBeat


		#Loop over all of the cap names
		counter=0
		for key in Data:
			#Define the Time array
			if counter==0: Time=np.linspace(0,Period,Data[key])

			outfile.write('Zone T="%s", I=%d, F=POINT\n'%(key,len(Data[key])))
			for i in range(Data[key]):
				outfile.write("%.06f %.06f %.06f %.06f %.06f\n"%(Time[i],Data[key][i][3],Data[key][i][7],Data[key][i][9],Data[key][i][8]*self.PressureScale))	

			counter+=1
		outfile.close()

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
	
	#Define the output folder
	parser.add_argument('-outfolder', '--OutFolder', type=str, required=True, dest="OutFolder", help="The path to the folder to store the post-processed files")

	#Start and End Argumenets		 
	parser.add_argument('-start', '--StartTimestep', type=int, required=True, dest="StartTimestep", help="The starting timestep to process the SimVascular files")
	
	#End Argument
	parser.add_argument('-stop', '--StopTimestep', type=int, required=True, dest="StopTimestep", help="The End timestep to process the SimVascular files")

	#Increment
	parser.add_argument('-incr', '--Increment', type=int, required=True, dest="Increment", help="The increment for the timestep files to process")

	parser.add_argument('-HeartBeat', '--HeartBeat', type=int, required=True, dest="HeartBeat", help="The Heart Beat for the patient to calculate the period of the cycle")

	#Put all the arguments together
	args=parser.parse_args()	

	#Call your Class
	PostProcess3DResults(args).Main()
