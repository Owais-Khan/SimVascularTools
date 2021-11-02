import vtk
import numpy
import os
import argparse
from glob import glob
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk

class PostProcess3DResults():
	def __init__(self,Args):
		self.Args=Args
		self.PressureScale=0.000750061683

		#Define Input Variables if not available
		#Input folder
		if self.Args.InputMeshFolder is None:
			self.Args.InputMeshFolder="./mesh-complete/"
		#Input and Output timesteps
		if self.Args.ResultsFolder is None:
			self.Args.ResultsFolder=glob("*-procs_case")[0]
		#Define Start and End Time step and Increment
		if self.Args.StartTimestep is None:
			self.Args.StartTimestep=int(sorted(glob(self.Args.ResultsFolder+"/all_results*.vtp"))[0].split("/")[-1][-9:-4])
		if self.Args.StopTimestep is None:
			self.Args.StopTimestep=int(sorted(glob(self.Args.ResultsFolder+"/all_results*.vtp"))[-1].split("/")[-1][-9:-4])
		if self.Args.Increment is None:
			self.Args.Increment=int(sorted(glob(self.Args.ResultsFolder+"/all_results*.vtp"))[1].split("/")[-1][-9:-4])-self.Args.StartTimestep
		
		#If outfolder is not defined
		if self.Args.OutFolder is None:
			os.system("mkdir ./PostProcessedData")
			self.Args.OutFolder="./PostProcessedData"

		#Define Heart Beat
		if self.Args.HeartBeat is None:
			self.Args.HeartBeat=60
		
		#Define tag for manually cut planes
		if (self.Args.NonCapPlanes is True) and (self.Args.NonCapPlanesFolder is None):
			self.Args.NonCapPlanesFolder="./ClippedPlanes"


	def Main(self):
		#Compute the averaged quantities over the Surface
		#Also compute velocity/pressure over outlet slices
		self.ComputeTemporalStatistics("Surface")

		#Compute the Average Quantities over the Volume
		#Also compute Velocity/Pressure over non-outlet slices in ./ClippedPlanes
		self.ComputeTemporalStatistics("Volume")
	
		#Now We want to add slices througout the model
		


	def ParticleTracer(self):
		#Create a folder for velocity vector
		os.system("%s/PartcileTracer"%self.ResultsFolder)
		os.system("%s/ParticleTracer/VelocityData"%self.ResultsFolder)
	
		N_cycles=6
	
		#Create a Symbolic Link to Velocity Data for 
                for i in range(self.Args.StartTimestep,self.Args.StopTimestep,self.Args.Increment):
                        VolumeFileName_ =glob("%s/*%.05d.vtu"%(self.Args.ResultsFolder,i))[0]
                        if len(VolumeFileName_)==0: print ("Could Not Find the File at TimeStep: %.05d"%i)
				
				

		os.system("vmtkmeshmergetimesteps -directory ./ -firststep 4000 -laststep 7000 -intervalstep 10 -pattern all_results.vtu_%5s.vtu -ofile mesh_timesteps.vtu -velocityvector 1 -vector velocity")


		os.system("vmtkparticletracer -ifile mesh_timesteps_temp.vtu -sfile source_append.vtp -ofile traces.vtp -maximumnumberofsteps 1000")
	
		os.system("vmtkrenderer -background 1 1 1 --pipe vmtksurfaceviewer -ifile all_results.vtp_04000.vtp  -opacity 0.03 --pipe vmtksurfaceviewer -ifile wall_svg_to_lcx_pre.tec -opacity 0.1 --pipe vmtksurfaceviewer -ifile wall_svg_to_ramus.tec -opacity 0.12 --pipe vmtkpathlineanimator -ifile traces.vtp -timestep 0.02 -legend 0 -maxtime 6 -pointsize 12 -colormap blackbody -screenshot 1 -imagesdirectory ./animations/")

	

	def ComputeTemporalStatistics(self,Tag):
		#Read the VTP and Vtu FileNames
		SurfacePointArray={}	
	
                #Get the Cap Names to compute outlet quantities
		if Tag=="Surface": CapNames=self.ReadCapNames() #Read the CapNames
		else:
			CapNames=glob(self.Args.NonCapPlanesFolder+"/*.vtp")
			for CapName_ in CapNames:
				os.system("vmtksurfacetomesh -ifile %s -ofile %s"%(CapName_,CapName_.replace(".vtp",".vtu")))

		
		CapAreas={}
		for CapName_ in CapNames:
			CapAreas[CapName_]=self.ComputeSurfaceArea(self.ReadVtpFile(CapName_)) 

                #Create a Dictionary to Store the Velocity Data
		CapVelAverage={}
		CapVelMax={}
		CapPres={}
		for CapName_ in CapNames: 
			CapVelAverage[CapName_]=[]
			CapVelMax[CapName_]=[]
			CapPres[CapName_]=[]
	
		#Loop over all the Files
		print ("--- Processing Surface Data")
		counter=0
		for i in range(self.Args.StartTimestep,self.Args.StopTimestep+self.Args.Increment,self.Args.Increment):
			if Tag=="Surface": SurfaceFileName_   =glob("%s/*%.05d.vtp"%(self.Args.ResultsFolder,i))[0]
			else: SurfaceFileName_   =glob("%s/*%.05d.vtu"%(self.Args.ResultsFolder,i))[0]
			if len(SurfaceFileName_)==0:
				print ("Could Not Find the File at TimeStep: %.05d"%i)
			print ("--- ---Looping over: %s"%SurfaceFileName_)
			
			
			if Tag=="Surface":SurfaceData_=self.ReadVtpFile(SurfaceFileName_)
			else: SurfaceData_=self.ReadVtuFile(SurfaceFileName_)
			

			#################### Compute the Values at Outlets #########################
			if i==self.Args.StartTimestep: os.system("mkdir %s/CapSurface/"%self.Args.OutFolder)
			for CapName_ in CapNames:
				CapOutFileName_="%s/CapSurfaceData/%s_%05d.vtp"%(self.Args.OutFolder,CapName_.split("/")[-1].replace(".vtp",""),i)
				if Tag=="Surface": 
					os.system("vmtksurfaceprojection -rfile %s -ifile %s -ofile %s"%(SurfaceFileName_,CapName_,CapOutFileName_))	
					CapSurface_=self.ReadVtpFile(CapOutFileName_)
				else:
					os.system("vmtkmeshprojection -rfile %s -ifile %s -ofile %s"%(SurfaceFileName_,CapName_.replace(".vtp",".vtu"),CapOutFileName_.replace(".vtp",".vtu")))
					CapSurface_=self.ReadVtuFile(CapOutFileName_.replace(".vtp",".vtu"))
				#Read Cap Data and Store Values in Array
				N_arrays=CapSurface_.GetPointData().GetNumberOfArrays()
				VelMag_=np.linalg.norm(CapSurface_.GetPointData().GetArray("velocity"),axis=1)
				CapVelMax[CapName_].append(np.max(VelMag_))
				CapVelAverage[CapName_].append(np.average(VelMag_))
				CapPres[CapName_].append(np.average(CapSurface_.GetPointData().GetArray("pressure")))
					

			################### Computing the Average File ###########################
			N_arrays=SurfaceData_.GetPointData().GetNumberOfArrays()
			for j in range(0,N_arrays):
				#Initialize the Arrays at the 0th timestep
				ArrayName_=SurfaceData_.GetPointData().GetArrayName(j)
				if i==self.Args.StartTimestep:
					SurfacePointArray[ArrayName_]=vtk_to_numpy(SurfaceData_.GetPointData().GetArray(ArrayName_))
					if ArrayName_.find("Velocity")>=0 or ArrayName_.find("velocity")>=0 or ArrayName_.find("WSS")>=0 or ArrayName_.find("wss")>=0:
						SurfacePointArray["%s_mag"%ArrayName_]=np.linalg.norm(SurfaceData_.GetPointData().GetArray(ArrayName_),axis=1)
				else:
					SurfacePointArray[ArrayName_]+=vtk_to_numpy(SurfaceData_.GetPointData().GetArray(ArrayName_))
					if ArrayName_.find("Velocity")>=0 or ArrayName_.find("velocity")>=0 or ArrayName_.find("WSS")>=0 or ArrayName_.find("wss")>=0:
						SurfacePointArray["%s_mag"%ArrayName_]+=np.linalg.norm(SurfaceData_.GetPointData().GetArray(ArrayName_),axis=1)


			counter+=1	


		###################### Store the Average File #########################
		SurfaceData_=self.UpdatePolyData(SurfaceData_,SurfacePointArray,factor=1./counter)	
		if Tag=="Surface":
			self.WriteVtpFile(self.Args.OutFolder+"/SurfaceData_TimeAveragedResults.vtp",SurfaceData_)
		else:
			self.WriteVtuFile(self.Args.OutFolder+"/VolumeData_TimeAveragedResults.vtu",SurfaceData_)
	

		###################### Write the Tecplot File #########################
		self.WriteTecplotCaps(CapVelAverage,CapVelMax,CapPres,CapAreas,Tag)

	def UpdatePolyData(self,SurfaceData,SurfacePointArrays,factor=1):
		#Remove all of the point arrays
		N_arrays=SurfaceData.GetPointData().GetNumberOfArrays()
		for i in range(N_arrays):
			SurfaceData.GetPointData().RemoveArray(SurfaceData.GetPointData().GetArrayName(i))
		#Append New Array
		for key in SurfacePointArrays:
			NewArray_=numpy_to_vtk(SurfacePointArrays[key]*factor)
			NewArray_.SetName(key)
			SurfaceData.GetPointData().AddArray(NewArray_)
		return SurfaceData
				
	
	def WriteTecplotCaps(self,DataVelAverage,DataVelMax,DataPres,CapAreas,Tag):
		if Tag=="Surface": outfile=open(self.Args.OutFolder+"/CapData.dat",'w')
		else: outfile=open(self.Args.OutFolder+"/NonCapData.dat",'w')
		#In the file header, add all of the time-averaged quantities
		outfile.write("#Temporal Averages\n")
		outfile.write("#CapName, Velocity, VelocityMaximum, CapArea, FlowRate, Pressure, SystolicPressure, DiastolicPressure \n")
		for key in DataVelAverage:
			CapName_=key.split("/")[-1].replace(".vtp","")
			outfile.write("#%s %.05f %.05f %.05f %.05f %.05f %.05f %.05f\n"%(CapName_,np.average(DataVelAverage[key]),np.average(DataVelMax[key]),CapAreas[key],np.average(DataVelAverage[key])*CapAreas[key],np.average(DataPres[key])*self.PressureScale,np.max(DataPres[key])*self.PressureScale,np.min(DataPres[key])*self.PressureScale))	
			
		outfile.write('TITLE="Outlet FLow Rates and Velocities"\n')
		outfile.write('VARIABLES="Time", "Velocity","VelocityMaximum","FlowRate","Pressure"\n')

		#Define the period
		Period=60/self.Args.HeartBeat

		#Loop over all of the cap names
		counter=0
		for key in DataVelAverage:
			#Define the Time array
			if counter==0: Time=np.linspace(0,Period,len(DataVelAverage[key]))

			outfile.write('Zone T="%s", I=%d, F=POINT\n'%(key,len(DataVelAverage[key])))
			for i in range(len(DataVelAverage[key])):
				outfile.write("%.06f %.06f %.06f %.06f %.06f %.06f %.06f\n"%(Time[i],DataVelAverage[key][i],DataVelMax[key][i],DataVelAverage[key][i]*CapAreas[key],DataPres[key][i]*self.PressureScale))
			

			counter+=1
		outfile.close()

	def ReadCapNames(self):
		#Grab all of the file names from the mesh-complete folder
		Caps_LCA  =sorted(glob(self.Args.InputMeshFolder+"/mesh-surfaces/lca_*.vtp"))	
		Caps_RCA  =sorted(glob(self.Args.InputMeshFolder+"/mesh-surfaces/rca_*.vtp"))	
		Caps_Aorta=sorted(glob(self.Args.InputMeshFolder+"/mesh-surfaces/aorta_*.vtp"))
		Caps_inflow=sorted(glob(self.Args.InputMeshFolder+"/mesh-surfaces/inflow.vtp"))
		
		#Put all of these names in a singe array
		CapNames=Caps_inflow+Caps_Aorta+Caps_LCA+Caps_RCA

		return CapNames	

	def ReadVtpFile(self,Filename):
		reader=vtk.vtkXMLPolyDataReader()
		reader.SetFileName(Filename)
		reader.Update()
		return reader.GetOutput()

	def ReadVtuFile(self,Filename):
		reader=vtk.vtkXMLUnstructuredGridReader()
		reader.SetFileName(Filename)
		reader.Update()
		return reader.GetOutput()

	def WriteVtpFile(self,Filename,Data):
		writer=vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(Filename)
		writer.SetInputData(Data)
		writer.Update()	

	def WriteVtuFile(self,Filename,Data):
		writer=vtk.vtkXMLUnstructuredGridWriter()
		writer.SetFileName(Filename)
		writer.SetInputData(Data)
		writer.Update()
		
	def ComputeSurfaceArea(self,PolyData):
		MassProperties=vtk.vtkMassProperties()
		MassProperties.SetInputData(PolyData)
		MassProperties.Update()
		return MassProperties.GetSurfaceArea()			
	

	def ComputeCentroid(self,PolyData):
		num_points = PolyData.GetNumberOfPoints()
		x_list = []
		y_list = []
		z_list = []
		for i in range(num_points):
			x_list.append(float(PolyData.GetPoints().GetPoint(i)[0]))
			y_list.append(float(PolyData.GetPoints().GetPoint(i)[1]))
			z_list.append(float(PolyData.GetPoints().GetPoint(i)[2]))
		return [np.average(x_list), np.average(y_list), np.average(z_list)]
	

if __name__=="__main__":
	#Arguments
	parser= argparse.ArgumentParser(description="This script will process the results from Simvascular")

	#Input filename for the mesh-complete folder
	parser.add_argument('-meshfolder', '--InputMeshFolder', type=str, required=False, dest="InputMeshFolder", help="The path to the mesh-complete folder that contains the volumetric and surface meshes") 

	#Input Folder where all of the results have been stored
	parser.add_argument('-resultsfolder', '--ResultsFolder', type=str, required=False, dest="ResultsFolder", help="The path to the folder that contains all all_data.vtp and all_data.vtu files postprocessed from SimVascular")
	
	#Define the output folder
	parser.add_argument('-outfolder', '--OutFolder', type=str, required=False, dest="OutFolder", help="The path to the folder to store the post-processed files")

	#Start and End Argumenets		 
	parser.add_argument('-start', '--StartTimestep', type=int, required=False, dest="StartTimestep", help="The starting timestep to process the SimVascular files")
	
	#End Argument
	parser.add_argument('-stop', '--StopTimestep', type=int, required=False, dest="StopTimestep", help="The End timestep to process the SimVascular files")

	#Increment
	parser.add_argument('-incr', '--Increment', type=int, required=False, dest="Increment", help="The increment for the timestep files to process")

	parser.add_argument('-HeartBeat', '--HeartBeat', type=int, required=False, dest="HeartBeat", help="The Heart Beat for the patient to calculate the period of the cycle")
	
	#Define the argument for non-cap planes 
	parser.add_argument('-NonCapPlanes', '--NonCapPlanes', type=bool, required=False, dest="NonCapPlanes", default=True, help="Tag to incidate whether there are cut-planes to processed that are not outlets")
	parser.add_argument('-NonCapPlanesFolder', '--NonCapPlanesFolder', type=str, required=False, dest="NonCapPlanesFolder", help="Folder that contains the slices for the planes that are not outlets but inside the domain. If not provided, default is to lool for 'ClippedPlanes' folder. ")

	#Put all the arguments together
	args=parser.parse_args()	

	#Call your Class
	PostProcess3DResults(args).Main()
